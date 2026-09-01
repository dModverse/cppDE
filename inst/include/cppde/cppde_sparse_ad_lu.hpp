/*
 AD-aware sparse LU solver for cppDE: raw CSC + KLU.

 Backend: KLU (KLU) or Eigen::SparseLU fallback.
 Operates on csc_matrix<T> with raw Ap/Ai/Ax arrays.

 Copyright (C) 2026 Simon Beyer
 */

#ifndef CPPDE_SPARSE_AD_LU_HPP
#define CPPDE_SPARSE_AD_LU_HPP

#include <type_traits>
#include <algorithm>
#include <vector>

#include <cppde/cppde_types.hpp>
#include <cppde/cppde_ad_lu.hpp>   // for is_ad, extract_values, etc.

#ifdef KLU
#include <cppde/cppde_klu_solver.hpp>
#endif

namespace cppde {
namespace ad_lu {

// ============================================================================
//  extract_csc_values: csc_matrix<dual<Inner>> to csc_matrix<Inner>
//
//  Copies the CSC structure (Ap, Ai) and extracts value part of Ax.
// ============================================================================

// Generic over any AD type with the standard accessor surface
// (cppde::dual, cppde::dual2nd).
template<class AD,
         std::enable_if_t<is_ad<AD>::value, int> = 0>
inline csc_matrix<inner_type_t<AD>>
extract_csc_values(const csc_matrix<AD>& W)
{
  using Inner = inner_type_t<AD>;
  csc_matrix<Inner> out;
  out.n   = W.n;
  out.nnz = W.nnz;
  out.Ap  = W.Ap;
  out.Ai  = W.Ai;
  out.Ax.resize(W.nnz);
  out.pattern_built = W.pattern_built;

  for (int k = 0; k < W.nnz; ++k)
    out.Ax[k] = const_cast<AD&>(W.Ax[k]).x();

  return out;
}

// ============================================================================
//  sparse_lu_solver<T>: Base case: T is a non-AD scalar
//
//  KLU: klu_lu_solver (preferred)
//  otherwise: compile error (KLU is now required for sparse)
// ============================================================================

template<class T, class Enable = void>
class sparse_lu_solver;

template<class Scalar>
class sparse_lu_solver<Scalar, std::enable_if_t<!is_ad<Scalar>::value>>
{
public:

#if defined(KLU)

  sparse_lu_solver() = default;
  sparse_lu_solver(sparse_lu_solver&& o) noexcept
    : m_solver(std::move(o.m_solver)) {}

  void factorize(const csc_matrix<Scalar>& W)
  {
    m_solver.factorize(W.n, W.Ap.data(), W.Ai.data(), W.Ax.data());
  }

  void analyze_pattern(const csc_matrix<Scalar>& W)
  {
    m_solver.analyze_pattern(W.n, W.Ap.data(), W.Ai.data());
  }

  void factorize_numeric(const csc_matrix<Scalar>& W)
  {
    m_solver.factorize(W.n, W.Ap.data(), W.Ai.data(), W.Ax.data());
  }

  void solve(std::vector<Scalar>& b) const
  { m_solver.solve(b); }

  // Batched solve: B ← W⁻¹ B (column-major n × nrhs)
  void solve_batch(std::vector<Scalar>& B, int nrhs) const
  { m_solver.solve_batch(B.data(), nrhs); }

  std::vector<Scalar> solve_copy(const std::vector<Scalar>& b) const
  {
    std::vector<Scalar> x = b;
    m_solver.solve(x);
    return x;
  }

  // Scalar-only solve (identity for base case: same as solve)
  void solve_scalar(std::vector<double>& b) const
  { m_solver.solve(b); }

  void reset_pattern() { m_solver.reset_pattern(); }

private:
  klu_lu_solver m_solver;

#else
  // Without KLU, sparse is not supported: static_assert at compile time
  static_assert(sizeof(Scalar) == 0,
                "cppDE sparse LU requires KLU (KLU). "
                "Rebuild the package or use dense mode.");
#endif
};

  // ============================================================================
  //  sparse_lu_solver<AD_T>: recursive AD case (IFT peeling)
  //
  //  As in the dense case: extract the value part, factorise, propagate the
  //  derivatives through the IFT. The CSC structure is shared between levels and
  //  copied once, and the solve buffers persist across calls.
  // ============================================================================

// Generic AD specialization for cppde::dual<Inner,N>.
template<class AD_T>
class sparse_lu_solver<AD_T, std::enable_if_t<is_ad<AD_T>::value>>
{
  using F = AD_T;
  using Inner = inner_type_t<AD_T>;
  static constexpr unsigned N = ad_traits::detail::ad_static_size<AD_T>::value;

  // Re-adopt W's CSC structure into the persistent buffers whenever it is not
  // the one already held.  Keying this on nnz alone would reuse a stale Ap/Ai
  // for a different pattern with the same entry count -- and would leave KLU
  // holding a symbolic factorization of the old graph.
  bool adopt_structure(const csc_matrix<F>& W)
  {
    if (m_W_val.n == W.n && m_W_val.nnz == W.nnz &&
        m_W_val.Ap == W.Ap && m_W_val.Ai == W.Ai)
      return false;

    m_W_val.n   = W.n;
    m_W_val.nnz = W.nnz;
    m_W_val.Ap  = W.Ap;
    m_W_val.Ai  = W.Ai;
    m_W_val.Ax.resize(W.nnz);
    m_W_val.pattern_built = W.pattern_built;
    m_Ap_cached = W.Ap;
    m_Ai_cached = W.Ai;
    m_inner.reset_pattern();
    return true;
  }

public:

  void factorize(const csc_matrix<F>& W)
  {
    m_n = W.n;
    const int nnz = W.nnz;

    // --- Extract scalar values into persistent m_W_val ---
    // First call (and any pattern change): copy structure (Ap, Ai).
    // Subsequent calls on the same pattern: only update Ax.
    adopt_structure(W);

    if constexpr (!is_ad<Inner>::value) {
  // ============================================================
  //  Inner = double: values and derivative components come out of W.Ax in one
  //  pass. The derivative block m_dW_ax is nnz by n_derivs, so the IFT matvec
  //  reads each CSC entry's derivatives contiguously. m_W_stored is not needed.
  // ============================================================

      // Determine n_derivs
      // Static width is claimed only when some entry of W is seeded; see
      // the dense factorize() in cppde_ad_lu.hpp.
      unsigned nd;
      if constexpr (N > 0) {
        nd = 0u;
        for (int k = 0; k < nnz; ++k) {
          if (const_cast<F&>(W.Ax[k]).size() > 0) { nd = N; break; }
        }
      } else {
        nd = 0;
        for (int k = 0; k < nnz; ++k) {
          unsigned sz = const_cast<F&>(W.Ax[k]).size();
          if (sz > nd) nd = sz;
        }
      }
      m_n_derivs_cached = nd;

      if (nd > 0)
        m_dW_ax.resize(static_cast<size_t>(nnz) * nd);

      // Fused extraction: values + derivatives in one pass
      for (int k = 0; k < nnz; ++k) {
        auto& wk = const_cast<F&>(W.Ax[k]);
        m_W_val.Ax[k] = static_cast<double>(wk.x());
        if (nd > 0) {
          unsigned wsz = wk.size();
          double* dst = m_dW_ax.data() + static_cast<size_t>(k) * nd;
          for (unsigned j = 0; j < nd; ++j)
            dst[j] = (j < wsz) ? static_cast<double>(wk.d(j)) : 0.0;
        }
      }

      // Ap/Ai for the IFT matvec are cached by adopt_structure() above.
      m_nnz_cached = nnz;

    } else {
      // ============================================================
      //  Inner = dual<...> (nested AD): store full AD matrix for the
      //  generic element-wise IFT path.
      // ============================================================
      if (m_W_stored.nnz != nnz) {
        m_W_stored = W;
      } else {
        std::copy(W.Ax.begin(), W.Ax.end(), m_W_stored.Ax.begin());
      }
      if constexpr (cppde::ad_traits::is_dual2nd<F>::value) {
        for (int k = 0; k < nnz; ++k)
          m_W_val.Ax[k] = first_order_view(W.Ax[k]);
      } else {
        for (int k = 0; k < nnz; ++k)
          m_W_val.Ax[k] = const_cast<F&>(W.Ax[k]).x();
      }
    }

    m_inner.factorize(m_W_val);
  }

  void analyze_pattern(const csc_matrix<F>& W)
  {
    // Structural analysis only: reuse m_W_val buffer
    const int nnz = W.nnz;
    m_n = W.n;
    adopt_structure(W);
    if constexpr (cppde::ad_traits::is_dual2nd<F>::value) {
      for (int k = 0; k < nnz; ++k)
        m_W_val.Ax[k] = first_order_view(W.Ax[k]);
    } else {
      for (int k = 0; k < nnz; ++k)
        m_W_val.Ax[k] = const_cast<F&>(W.Ax[k]).x();
    }
    m_inner.analyze_pattern(m_W_val);
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

    // 2. Determine derivative directions
    // Static width is claimed only when b or W is seeded; see
    // dense_lu_solver::solve() in cppde_ad_lu.hpp.
    unsigned n_derivs;
    if constexpr (!is_ad<Inner>::value) {
      if constexpr (N > 0) {
        n_derivs = (m_n_derivs_cached > 0 || max_deriv_size(b) > 0) ? N : 0u;
      } else {
        n_derivs = std::max(m_n_derivs_cached, max_deriv_size(b));
      }
    } else {
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

    // 4. IFT sparse matvec: rhs_all -= dW · b_val
    ift_sparse_matvec(n, n_derivs);

    // 5. Batched solve
    m_inner.solve_batch(m_rhs_all, static_cast<int>(n_derivs));

    // 6. Bulk-inject results
    bulk_inject_results(b, m_b_val, m_rhs_all, n, n_derivs);
  }

  // Batched solve for nrhs RHS vectors (column-major). Fully BLAS-3 / KLU-batched
  // at the recursion level: amortises buffer extraction and IFT work across the
  // batch and feeds the inner solver one bulk solve_batch call per layer.
  void solve_batch(std::vector<F>& B_flat, int nrhs) const
  {
    if (nrhs <= 0) return;
    const int n = m_n;
    const std::size_t total = static_cast<std::size_t>(n) * nrhs;

    // 1. Bulk-extract value layer.
    m_b_val_batch.resize(total);
    if constexpr (cppde::ad_traits::is_dual2nd<F>::value) {
      for (int col = 0; col < nrhs; ++col)
        for (int i = 0; i < n; ++i)
          m_b_val_batch[col * n + i] =
            first_order_view(B_flat[col * n + i]);
    } else {
      for (int col = 0; col < nrhs; ++col)
        for (int i = 0; i < n; ++i)
          m_b_val_batch[col * n + i] =
            const_cast<F&>(B_flat[col * n + i]).x();
    }

    // 2. Batched value-layer solve.
    m_inner.solve_batch(m_b_val_batch, nrhs);

    // 3. Determine n_derivs.
    // Static width is claimed only when b or W is seeded; see
    // dense_lu_solver::solve() in cppde_ad_lu.hpp.
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
        for (int col = 0; col < nrhs; ++col)
          for (int i = 0; i < n; ++i) {
            unsigned sz = const_cast<F&>(B_flat[col * n + i]).size();
            if (sz > n_derivs) n_derivs = sz;
          }
      }
    } else {
      if constexpr (N > 0) {
        n_derivs = (any_seeded_rhs() || max_deriv_size(m_W_stored) > 0) ? N : 0u;
      } else {
        n_derivs = max_deriv_size(m_W_stored);
        for (int col = 0; col < nrhs; ++col)
          for (int i = 0; i < n; ++i) {
            unsigned sz = const_cast<F&>(B_flat[col * n + i]).size();
            if (sz > n_derivs) n_derivs = sz;
          }
      }
    }

    if (n_derivs == 0) {
      for (int col = 0; col < nrhs; ++col)
        for (int i = 0; i < n; ++i)
          B_flat[col * n + i].x() = m_b_val_batch[col * n + i];
      return;
    }

    // 4. Bulk-extract derivative RHS.
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

    // 5. IFT sparse matvec batched.
    ift_sparse_matvec_batch(n, n_derivs, nrhs);

    // 6. Batched derivative-layer solve.
    m_inner.solve_batch(m_rhs_all_batch,
                        static_cast<int>(n_derivs * nrhs));

    // 7. Bulk-inject results.
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

  void reset_pattern() { m_inner.reset_pattern(); }

  // Scalar-only solve: forward to inner solver, no IFT
  void solve_scalar(std::vector<double>& b) const
  {
    m_inner.solve_scalar(b);
  }

private:

    // ================================================================
    //  IFT sparse matvec: rhs_all -= dW * b_val.
    //
    //  For Inner = double the derivatives were pre-extracted in factorize(),
    //  laid out as m_dW_ax[p * nd + j] for the p-th CSC entry, so this is plain
    //  double arithmetic. A nested inner type takes one element-wise pass.
    // ================================================================

    // Batched IFT sparse matvec: an outer loop over the batch columns, each
    // reading its own slice of the batch buffers, with the same inner CSC pass
    // as the single right-hand side.
  void ift_sparse_matvec_batch(int n, unsigned n_derivs, int nrhs) const
  {
    const std::size_t per_col = static_cast<std::size_t>(n) * n_derivs;
    if constexpr (!is_ad<Inner>::value) {
      unsigned nd_W = m_n_derivs_cached;
      if (nd_W == 0) return;

      const int* Ap = m_Ap_cached.data();
      const int* Ai = m_Ai_cached.data();

      for (int col_b = 0; col_b < nrhs; ++col_b) {
        double* rhs_col = m_rhs_all_batch.data() + col_b * per_col;
        const double* b_val_col = m_b_val_batch.data() + col_b * n;
        for (int col = 0; col < n; ++col) {
          const double b_col = b_val_col[col];
          for (int p = Ap[col]; p < Ap[col + 1]; ++p) {
            int row = Ai[p];
            const double* dw = m_dW_ax.data() + static_cast<size_t>(p) * nd_W;
            for (unsigned j = 0; j < nd_W; ++j)
              rhs_col[j * n + row] -= dw[j] * b_col;
          }
        }
      }
    } else {
      for (int col_b = 0; col_b < nrhs; ++col_b) {
        Inner* rhs_col = m_rhs_all_batch.data() + col_b * per_col;
        const Inner* b_val_col = m_b_val_batch.data() + col_b * n;
        for (int col = 0; col < n; ++col) {
          const Inner b_col = b_val_col[col];
          for (int p = m_W_stored.Ap[col]; p < m_W_stored.Ap[col + 1]; ++p) {
            int row = m_W_stored.Ai[p];
            auto& w_entry = const_cast<F&>(m_W_stored.Ax[p]);
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

  void ift_sparse_matvec(int n, unsigned n_derivs) const
  {
    if constexpr (!is_ad<Inner>::value) {
      // =============================================================
      //  Pre-extracted path: iterate CSC with flat double derivatives
      //
      //  m_dW_ax was built in factorize().  If n_derivs > nd_W
      //  (RHS has more derivs than W), the extra directions have
      //  dW = 0, so no correction is needed for them.
      // =============================================================

      unsigned nd_W = m_n_derivs_cached;
      if (nd_W == 0) return;

      const int* Ap = m_Ap_cached.data();
      const int* Ai = m_Ai_cached.data();

      for (int col = 0; col < n; ++col) {
        const double b_col = m_b_val[col];
        for (int p = Ap[col]; p < Ap[col + 1]; ++p) {
          int row = Ai[p];
          const double* dw = m_dW_ax.data() + static_cast<size_t>(p) * nd_W;
          for (unsigned j = 0; j < nd_W; ++j)
            m_rhs_all[j * n + row] -= dw[j] * b_col;
        }
      }

    } else {
      // =============================================================
      //  Generic path: nested AD types (Inner = dual<...>)
      //
      //  Single pass over W_stored: each dual<Inner,N> element is
      //  touched exactly once.
      // =============================================================

      for (int col = 0; col < n; ++col) {
        const Inner b_col = m_b_val[col];
        for (int p = m_W_stored.Ap[col]; p < m_W_stored.Ap[col + 1]; ++p) {
          int row = m_W_stored.Ai[p];
          auto& w_entry = const_cast<F&>(m_W_stored.Ax[p]);
          unsigned wsz = w_entry.size();
          for (unsigned j = 0; j < n_derivs; ++j) {
            Inner dw = (j < wsz) ? w_entry.d(j) : Inner(0);
            m_rhs_all[j * n + row] -= dw * b_col;
          }
        }
      }
    }
  }

  int m_n = 0;
  csc_matrix<F>             m_W_stored;  // Full AD matrix (nested AD path only)
  csc_matrix<Inner> m_W_val;             // Persistent scalar extraction buffer
  sparse_lu_solver<Inner> m_inner;

  // Pre-extracted derivative block for IFT (Inner = double only).
  // Built once in factorize(), reused across all solve() calls.
  // Layout: m_dW_ax[p * n_derivs + j] = derivative j of CSC entry p.
  std::vector<double> m_dW_ax;
  unsigned m_n_derivs_cached = 0;
  int m_nnz_cached = 0;
  std::vector<int> m_Ap_cached;       // CSC column pointers (for IFT matvec)
  std::vector<int> m_Ai_cached;       // CSC row indices (for IFT matvec)

  // Persistent solve buffers.
  mutable std::vector<Inner> m_b_val;          // n scalars: value part (single)
  mutable std::vector<Inner> m_rhs_all;        // n × n_derivs (single)
  mutable std::vector<F>     m_col_buf;        // legacy column buffer
  mutable std::vector<Inner> m_b_val_batch;    // n × nrhs: value part (batched)
  mutable std::vector<Inner> m_rhs_all_batch;  // n × n_derivs × nrhs (batched)
};

// ============================================================================
//  Sparse Jacobian matrix-vector product: y += W * x
// ============================================================================

template<class T>
void sparse_jac_matvec(const csc_matrix<T>& W,
                       const std::vector<T>& x,
                       std::vector<T>& y)
{
  csc_matvec_add(W, x, y);
}

} // namespace ad_lu
} // namespace cppde

#endif // CPPDE_SPARSE_AD_LU_HPP
