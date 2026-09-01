/*
 AD-aware dense LU solver for cppDE: LAPACK backend.

 Uses dgetrf/dgetrs from R's bundled LAPACK for the base case.
 Nested AD types (dual<dual<T,N>,N>) are handled by recursive IFT peeling.

 BLAS-3 optimization: IFT derivative propagation uses batched
 dgetrs (nrhs = n_derivs) instead of n_derivs separate solves.
 The matvec phase fuses extraction and subtraction in a single
 pass over W_stored for cache-optimal access.

 Copyright (C) 2026 Simon Beyer
 */

#ifndef CPPDE_AD_LU_HPP
#define CPPDE_AD_LU_HPP

#include <type_traits>
#include <algorithm>
#include <vector>
#include <cstring>
#include <cstdlib>

#include <cppde/cppde_types.hpp>
#include <cppde/cppde_ad_traits.hpp>
#include <cppde/cppde_dual2nd.hpp>

// LAPACK declarations via R's headers (handles ILP64 automatically)
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>


namespace cppde {
namespace ad_lu {

// ============================================================================
//  AD type traits / bulk helpers: pulled in from cppde_ad_traits.hpp.
//  Re-exported here so existing consumers using `ad_lu::is_ad`,
//  `ad_lu::scalar_value`, etc. compile unchanged.
// ============================================================================

using cppde::ad_traits::is_ad;
using cppde::ad_traits::inner_type;
using cppde::ad_traits::inner_type_t;
using cppde::ad_traits::scalar_type;
using cppde::ad_traits::scalar_type_t;
using cppde::ad_traits::scalar_value;
using cppde::ad_traits::extract_values;
using cppde::ad_traits::extract_derivs;
using cppde::ad_traits::max_deriv_size;
using cppde::ad_traits::any_deriv;
using cppde::ad_traits::bulk_extract_derivs;
using cppde::ad_traits::bulk_inject_results;

// ============================================================================
//  dense_lu_solver<double>: Base case: LAPACK dgetrf + dgetrs
// ============================================================================

template<class T, class Enable = void>
class dense_lu_solver;

template<class Scalar>
class dense_lu_solver<Scalar, std::enable_if_t<!is_ad<Scalar>::value>>
{
public:
  void factorize(const dense_matrix<Scalar>& W)
  {
    const int n = W.rows();
    m_n = n;
    m_lu_data = W.data;
    m_ipiv.resize(n);

    int info = 0;
    F77_CALL(dgetrf)(&m_n, &m_n, m_lu_data.data(), &m_n, m_ipiv.data(), &info);
  }

  // Factorize by MOVING the data out of the source matrix (zero-copy).
  // The source matrix is left in a valid but unspecified state.
  // Use when the caller no longer needs the original W (e.g. m_W_temp
  // which gets rebuilt from scratch on every factorize_W call).
  void factorize_move(dense_matrix<Scalar>& W)
  {
    const int n = W.rows();
    m_n = n;
    m_lu_data.swap(W.data);   // O(1) pointer swap, no memcpy
    m_ipiv.resize(n);

    int info = 0;
    F77_CALL(dgetrf)(&m_n, &m_n, m_lu_data.data(), &m_n, m_ipiv.data(), &info);
  }

  // Solve in-place: b ← W⁻¹ b
  void solve(std::vector<Scalar>& b) const
  {
    char trans = 'N';
    int nrhs = 1;
    int info = 0;
    F77_CALL(dgetrs)(&trans, &m_n, &nrhs,
             const_cast<double*>(m_lu_data.data()), &m_n,
             const_cast<int*>(m_ipiv.data()),
             b.data(), &m_n, &info
                       FCONE);
  }

  // Batched solve: B ← W⁻¹ B  (B is column-major n × nrhs)
  //
  // Uses BLAS-3 internally (dtrsm) via dgetrs with nrhs > 1.
  // This is the key optimization: cache-blocked triangular solves
  // instead of nrhs separate memory-bound triangular back-subs.
  void solve_batch(std::vector<Scalar>& B, int nrhs) const
  {
    if (nrhs <= 0) return;
    char trans = 'N';
    int info = 0;
    F77_CALL(dgetrs)(&trans, &m_n, &nrhs,
             const_cast<double*>(m_lu_data.data()), &m_n,
             const_cast<int*>(m_ipiv.data()),
             B.data(), &m_n, &info
                       FCONE);
  }

  // Solve with separate output
  std::vector<Scalar> solve_copy(const std::vector<Scalar>& b) const
  {
    std::vector<Scalar> x = b;
    solve(x);
    return x;
  }

  // Scalar-only solve (identity for base case: same as solve)
  void solve_scalar(std::vector<double>& b) const
  { solve(b); }

private:
  mutable int m_n = 0;
  std::vector<double> m_lu_data;   // LU-factorized copy of W
  std::vector<int>    m_ipiv;      // pivot indices
};

    // ============================================================================
    //  dense_lu_solver<AD_T>: recursive AD case (IFT peeling)
    //
    //    value:       W_val * x_val = b_val
    //    direction j: W_val * dx_j  = db_j - dW_j * x_val, reusing the factors
    //
    //  The scalar matrix, the stored AD matrix and the solve buffers persist
    //  across calls, and the derivative layer goes out as one batched solve.
    // ============================================================================

// Generic AD specialization for cppde::dual<Inner,N>.
// The body uses only the standard accessor surface (.x(), .d(j),
// .size(), .diff(), operator[], .depend()).
template<class AD_T>
class dense_lu_solver<AD_T, std::enable_if_t<is_ad<AD_T>::value>>
{
  using F = AD_T;
  using Inner = inner_type_t<AD_T>;
  static constexpr unsigned N = ad_traits::detail::ad_static_size<AD_T>::value;

public:

  void factorize(const dense_matrix<F>& W)
  {
    const int n = W.rows();
    m_n = n;
    const int nn = n * n;

    // --- Extract scalar values into persistent m_W_val ---
    // First call: allocate.  Subsequent: only update values.
    if (m_W_val.rows() != n)
      m_W_val.resize(n, n);
    if constexpr (cppde::ad_traits::is_dual2nd<F>::value) {
      // dual2nd: synthesise the value-layer dual<S, N> from the inline
      // gradient (outer.tan_[k].x()), bypassing the redundant val_tan_block.
      for (int k = 0; k < nn; ++k)
        m_W_val.data[k] = first_order_view(W.data[k]);
    } else {
      for (int k = 0; k < nn; ++k)
        m_W_val.data[k] = const_cast<F&>(W.data[k]).x();
    }

    if constexpr (!is_ad<Inner>::value) {
    // ============================================================
    //  Inner = double: the IFT matvec needs the derivative part of W as a
    //  column-major block. W does not change between factorize() and solve(),
    //  so it is extracted once here and m_W_stored is not needed at all.
    // ============================================================

    // Determine n_derivs. The cached width has to reflect actual dependence,
    // not the compile-time width: the error-weight vector is sized from the
    // seeded directions and indexed per active direction.
      unsigned nd;
      if constexpr (N > 0) {
        nd = any_deriv(W) ? N : 0u;
      } else {
        nd = max_deriv_size(W);
      }
      m_n_derivs_cached = nd;

      if (nd > 0) {
        int m = n * static_cast<int>(nd);
        m_dW_block.resize(static_cast<size_t>(m) * n);

        // AoS to SoA deinterleave: extract all derivative components
        for (int col = 0; col < n; ++col) {
          double* col_dst = m_dW_block.data() + static_cast<size_t>(col) * m;
          for (int row = 0; row < n; ++row) {
            auto& w = const_cast<F&>(W(row, col));
            unsigned wsz = w.size();
            for (unsigned j = 0; j < nd; ++j)
              col_dst[j * n + row] = (j < wsz) ?
                static_cast<double>(w.d(j)) : 0.0;
          }
        }
      }

    } else {
      // ============================================================
      //  Inner = dual<...> (nested AD): store full AD matrix for the
      //  generic element-wise IFT path.
      // ============================================================
      if (m_W_stored.rows() != n) {
        m_W_stored = W;
      } else {
        std::copy(W.data.begin(), W.data.end(), m_W_stored.data.begin());
      }
    }

    // Recursive: inner solver factorizes the scalar matrix.
    m_inner.factorize(m_W_val);
  }

  void solve(std::vector<F>& b) const
  {
    const int n = m_n;

    // 1. Extract and solve value part (reuse buffer)
    m_b_val.resize(n);
    if constexpr (cppde::ad_traits::is_dual2nd<F>::value) {
      for (int i = 0; i < n; ++i)
        m_b_val[i] = first_order_view(b[i]);
    } else {
      for (int i = 0; i < n; ++i)
        m_b_val[i] = const_cast<F&>(b[i]).x();
    }
    m_inner.solve(m_b_val);

    // 2. Determine derivative directions. A static width N is claimed only when
    // some input is seeded: injecting results activates tangents on b, and the
    // active directions must match the error-weight vector of this state.
    unsigned n_derivs;
    if constexpr (!is_ad<Inner>::value) {
      // Inner = double: n_derivs was determined at factorize time.
      // For RHS-only derivatives (W has none), check b too.
      if constexpr (N > 0) {
        n_derivs = (m_n_derivs_cached > 0 || max_deriv_size(b) > 0) ? N : 0u;
      } else {
        n_derivs = std::max(m_n_derivs_cached, max_deriv_size(b));
      }
    } else {
      // Inner = dual<...>: need to scan both b and W_stored
      if constexpr (N > 0) {
        n_derivs = (max_deriv_size(b) > 0 || max_deriv_size(m_W_stored) > 0)
                     ? N : 0u;
      } else {
        n_derivs = max_deriv_size(b);
        n_derivs = std::max(n_derivs, max_deriv_size(m_W_stored));
      }
    }
    if (n_derivs == 0) {
      for (int i = 0; i < n; ++i)
        b[i].x() = m_b_val[i];
      return;
    }

    // 3. Bulk-extract ALL derivative RHS into column-major n × n_derivs
    m_rhs_all.resize(static_cast<size_t>(n) * n_derivs);
    for (int i = 0; i < n; ++i) {
      auto& bi = const_cast<F&>(b[i]);
      unsigned sz = bi.size();
      for (unsigned j = 0; j < n_derivs; ++j)
        m_rhs_all[j * n + i] = (j < sz) ? bi.d(j) : Inner(0);
    }

    // 4. IFT matvec: rhs_all -= dW · b_val
    ift_matvec(n, n_derivs);

    // 5. Batched solve: solve ALL n_derivs RHS in one call
    m_inner.solve_batch(m_rhs_all, static_cast<int>(n_derivs));

    // 6. Bulk-inject: write values and derivatives in one pass
    bulk_inject_results(b, m_b_val, m_rhs_all, n, n_derivs);
  }

    // Batched solve over nrhs column-major right-hand sides: one value-layer
    // solve, one IFT correction as a dgemm or a fused pass, one derivative-layer
    // solve. Saves nrhs - 1 inner dispatches and amortises the buffers.
  void solve_batch(std::vector<F>& B_flat, int nrhs) const
  {
    if (nrhs <= 0) return;
    const int n = m_n;
    const std::size_t total = static_cast<std::size_t>(n) * nrhs;

    // 1. Bulk-extract value layer for all nrhs columns into a flat buffer.
    m_b_val_batch.resize(total);
    if constexpr (cppde::ad_traits::is_dual2nd<F>::value) {
      for (int col = 0; col < nrhs; ++col) {
        for (int i = 0; i < n; ++i)
          m_b_val_batch[col * n + i] =
            first_order_view(B_flat[col * n + i]);
      }
    } else {
      for (int col = 0; col < nrhs; ++col) {
        for (int i = 0; i < n; ++i)
          m_b_val_batch[col * n + i] =
            const_cast<F&>(B_flat[col * n + i]).x();
      }
    }

    // 2. Batched value-layer solve.
    m_inner.solve_batch(m_b_val_batch, nrhs);

    // 3. Determine n_derivs: scan all RHS columns plus W (when nested AD).
    // Static width is claimed only when some right-hand side or W is
    // seeded; see dense_lu_solver::solve().
    auto any_seeded_rhs = [&]() {
      for (int col = 0; col < nrhs; ++col)
        for (int i = 0; i < n; ++i)
          if (const_cast<F&>(B_flat[col * n + i]).size() > 0) return true;
      return false;
    };
    unsigned n_derivs;
    if constexpr (!is_ad<Inner>::value) {
      if constexpr (N > 0) {
        n_derivs = (m_n_derivs_cached > 0 || any_seeded_rhs()) ? N : 0u;
      } else {
        n_derivs = m_n_derivs_cached;
        for (int col = 0; col < nrhs; ++col) {
          for (int i = 0; i < n; ++i) {
            unsigned sz = const_cast<F&>(B_flat[col * n + i]).size();
            if (sz > n_derivs) n_derivs = sz;
          }
        }
      }
    } else {
      if constexpr (N > 0) {
        n_derivs = (any_seeded_rhs() || max_deriv_size(m_W_stored) > 0) ? N : 0u;
      } else {
        n_derivs = max_deriv_size(m_W_stored);
        for (int col = 0; col < nrhs; ++col) {
          for (int i = 0; i < n; ++i) {
            unsigned sz = const_cast<F&>(B_flat[col * n + i]).size();
            if (sz > n_derivs) n_derivs = sz;
          }
        }
      }
    }

    if (n_derivs == 0) {
      // No derivatives anywhere: just write back the value layer.
      for (int col = 0; col < nrhs; ++col)
        for (int i = 0; i < n; ++i)
          B_flat[col * n + i].x() = m_b_val_batch[col * n + i];
      return;
    }

    // 4. Bulk-extract all derivative RHS into a flat buffer.
    //    Layout: column-major (n * n_derivs) x nrhs. Within each column,
    //    the n_derivs blocks of length n are laid out by tangent index j
    //    (consistent with the existing m_rhs_all layout for solve()).
    const std::size_t per_col = static_cast<std::size_t>(n) * n_derivs;
    m_rhs_all_batch.resize(per_col * nrhs);
    for (int col = 0; col < nrhs; ++col) {
      for (int i = 0; i < n; ++i) {
        auto& bi = const_cast<F&>(B_flat[col * n + i]);
        unsigned sz = bi.size();
        for (unsigned j = 0; j < n_derivs; ++j)
          m_rhs_all_batch[col * per_col + j * n + i] =
            (j < sz) ? bi.d(j) : Inner(0);
      }
    }

    // 5. IFT matvec (batched): rhs_all_batch -= dW · b_val_batch.
    ift_matvec_batch(n, n_derivs, nrhs);

    // 6. Batched derivative-layer solve: n_derivs * nrhs columns at once.
    m_inner.solve_batch(m_rhs_all_batch,
                        static_cast<int>(n_derivs * nrhs));

    // 7. Bulk-inject values + derivatives back into B_flat.
    constexpr unsigned StaticN = N;
    for (int col = 0; col < nrhs; ++col) {
      for (int i = 0; i < n; ++i) {
        auto& bi = B_flat[col * n + i];
        bi.x() = m_b_val_batch[col * n + i];
        if constexpr (StaticN > 0) {
          if (!bi.depend()) bi.diff(0);
        } else {
          if (!bi.depend()) bi.diff(0, n_derivs);
        }
        for (unsigned j = 0; j < n_derivs; ++j)
          bi[j] = m_rhs_all_batch[col * per_col + j * n + i];
      }
    }
  }

  // Scalar-only solve: use the value-level factorization on plain doubles.
  // No IFT, no derivative propagation: just the base-case LAPACK solve.
  void solve_scalar(std::vector<double>& b) const
  {
    m_inner.solve_scalar(b);
  }

private:

    // ================================================================
    //  IFT matvec: rhs_all -= dW * b_val. For Inner = double this is the dgemv
    //  on the block extracted in factorize(); for a nested inner type it is one
    //  element-wise pass over m_W_stored.
    // ================================================================

  void ift_matvec(int n, unsigned n_derivs) const
  {
    if constexpr (!is_ad<Inner>::value) {
    // =============================================================
    //  BLAS-2 path: m_dW_block is (n * n_derivs_cached) x n, column-major. A
    //  right-hand side carrying more directions than W has dW = 0 for those, so
    //  dgemv covers the cached range and the rest stays untouched.
    // =============================================================

      unsigned nd_W = m_n_derivs_cached;
      if (nd_W == 0) return;  // W has no derivatives, no IFT correction

      int m = n * static_cast<int>(nd_W);
      double alpha = -1.0, beta = 1.0;
      int inc = 1;
      char trans = 'N';
      F77_CALL(dgemv)(&trans, &m, &n, &alpha,
               const_cast<double*>(m_dW_block.data()), &m,
               const_cast<double*>(m_b_val.data()), &inc,
               &beta, m_rhs_all.data(), &inc  FCONE);

    } else {
      // =============================================================
      //  Generic path: nested AD types (Inner = dual<...>)
      //
      //  Single pass over W_stored: each dual<Inner,N> element is
      //  touched exactly once.
      // =============================================================

      for (int col = 0; col < n; ++col) {
        const Inner b_col = m_b_val[col];
        for (int row = 0; row < n; ++row) {
          auto& w_entry = const_cast<F&>(m_W_stored(row, col));
          unsigned wsz = w_entry.size();
          for (unsigned j = 0; j < n_derivs; ++j) {
            Inner dw = (j < wsz) ? w_entry.d(j) : Inner(0);
            m_rhs_all[j * n + row] -= dw * b_col;
          }
        }
      }
    }
  }

    // Batched IFT matvec: rhs_all[col, j, row] -= sum_k dW[row, k, j] *
    // b_val[col, k]. One dgemm for Inner = double, otherwise an outer loop over
    // the batch columns with the single-RHS pass inside.
  void ift_matvec_batch(int n, unsigned n_derivs, int nrhs) const
  {
    if constexpr (!is_ad<Inner>::value) {
      unsigned nd_W = m_n_derivs_cached;
      if (nd_W == 0) return;

      int m = n * static_cast<int>(nd_W);
      double alpha = -1.0, beta = 1.0;
      char trans = 'N';
      // dgemm: C[m × nrhs] = beta * C + alpha * A[m × n] * B[n × nrhs].
      // m_rhs_all_batch is laid out column-major (n*n_derivs) × nrhs. When
      // n_derivs > nd_W the unused tail of each column stays untouched
      // (correct: dW = 0 in those tangent directions).
      const std::size_t per_col = static_cast<std::size_t>(n) * n_derivs;
      // Pass column stride for C as per_col so dgemm steps over the full
      // n_derivs blocks even though it only writes the first nd_W of them.
      int ldc = static_cast<int>(per_col);
      F77_CALL(dgemm)(&trans, &trans, &m, &nrhs, &n, &alpha,
               const_cast<double*>(m_dW_block.data()), &m,
               m_b_val_batch.data(), &n,
               &beta, m_rhs_all_batch.data(), &ldc
                       FCONE FCONE);

    } else {
      const std::size_t per_col = static_cast<std::size_t>(n) * n_derivs;
      for (int col_b = 0; col_b < nrhs; ++col_b) {
        Inner* rhs_col = m_rhs_all_batch.data() + col_b * per_col;
        const Inner* b_val_col = m_b_val_batch.data() + col_b * n;
        for (int col = 0; col < n; ++col) {
          const Inner b_col = b_val_col[col];
          for (int row = 0; row < n; ++row) {
            auto& w_entry = const_cast<F&>(m_W_stored(row, col));
            unsigned wsz = w_entry.size();
            for (unsigned j = 0; j < n_derivs; ++j) {
              Inner dw = (j < wsz) ? w_entry.d(j) : Inner(0);
              rhs_col[j * n + row] -= dw * b_col;
            }
          }
        }
      }
    }
  }

  int m_n = 0;
  dense_matrix<F>     m_W_stored;      // Full AD matrix (nested AD path only)
  dense_matrix<Inner> m_W_val;         // Persistent scalar extraction buffer
  dense_lu_solver<Inner> m_inner;

  // Pre-extracted derivative block for BLAS-2 IFT (Inner = double only).
  // Built once in factorize(), reused across all solve() calls.
  // Layout: (n*n_derivs_cached) × n, column-major.
  std::vector<double> m_dW_block;
  unsigned m_n_derivs_cached = 0;     // n_derivs at last factorize

  // Persistent solve buffers: avoid per-call heap allocations.
  mutable std::vector<Inner>  m_b_val;          // n scalars: value part (single solve)
  mutable std::vector<Inner>  m_rhs_all;        // n × n_derivs: derivs (single solve)
  mutable std::vector<F>      m_col_buf;        // n entries: legacy column buffer
  mutable std::vector<Inner>  m_b_val_batch;    // n × nrhs: value part (batched)
  mutable std::vector<Inner>  m_rhs_all_batch;  // n × n_derivs × nrhs (batched)
};

} // namespace ad_lu
} // namespace cppde

#endif // CPPDE_AD_LU_HPP
