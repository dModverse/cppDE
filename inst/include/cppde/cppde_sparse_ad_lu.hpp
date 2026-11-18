/*
 AD-aware sparse LU solver for cppDE: raw CSC + KLU.

 Backend: KLU; without it the sparse path does not compile.
 Operates on csc_matrix<T> with raw Ap/Ai/Ax arrays.

 Copyright (C) 2026 Simon Beyer
 */

#ifndef CPPDE_SPARSE_AD_LU_HPP
#define CPPDE_SPARSE_AD_LU_HPP

#include <type_traits>
#include <algorithm>
#include <vector>

#include <cppde/cppde_types.hpp>
#include <cppde/cppde_ad_lu.hpp>   // for is_ad, inner_type_t, etc.

#ifdef KLU
#include <cppde/cppde_klu_solver.hpp>
#endif

namespace cppde {
namespace ad_lu {

// ============================================================================
//  sparse_lu_solver<T>: Base case: T is a non-AD scalar
//
//  KLU through klu_lu_solver; a compile error without KLU.
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

  void solve(std::vector<Scalar>& b) const
  { m_solver.solve(b); }

  void solve_transposed(std::vector<Scalar>& b) const
  { m_solver.solve_transposed(b); }

  // Batched solve: B ← W⁻¹ B (column-major n × nrhs)
  void solve_batch(std::vector<Scalar>& B, int nrhs) const
  { m_solver.solve_batch(B.data(), nrhs); }

  void reset_pattern() { m_solver.reset_pattern(); }

  // See klu_lu_solver::refactor_solve_ratio.
  double refactor_solve_ratio() { return m_solver.refactor_solve_ratio(); }

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
  //  As in the dense case: extract the value part, factorize, propagate the
  //  derivatives through the IFT. The CSC structure is shared between levels and
  //  copied once, and the solve buffers persist across calls.
  // ============================================================================

// Generic AD specialization for cppde::dual<Inner>.
template<class AD_T>
class sparse_lu_solver<AD_T, std::enable_if_t<is_ad<AD_T>::value>>
{
  using F = AD_T;
  using Inner = inner_type_t<AD_T>;

  // Re-adopt W's CSC structure into the persistent buffers whenever it is not
  // the one already held.  Keying this on nnz alone would reuse a stale Ap/Ai
  // for a different pattern with the same entry count, and would leave KLU
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

  // A frozen matrix is factorized on the values alone: a solve passes the
  // tangents of the right-hand side through it and never the tangents of W.
  void set_frozen(bool f) {
    m_frozen = f;
    if constexpr (is_ad<Inner>::value) m_inner.set_frozen(f);
  }

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

      // Determine n_derivs from the seeded entries of W.
      unsigned nd = 0;
      for (int k = 0; k < nnz && !m_frozen; ++k) {
        unsigned sz = const_cast<F&>(W.Ax[k]).size();
        if (sz > nd) nd = sz;
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
      if (m_frozen) {
        m_W_stored = csc_matrix<F>();
      } else if (m_W_stored.nnz != nnz) {
        m_W_stored = W;
      } else {
        std::copy(W.Ax.begin(), W.Ax.end(), m_W_stored.Ax.begin());
      }
      if constexpr (cppde::ad_traits::is_dual2nd<F>::value) {
        for (int k = 0; k < nnz; ++k)
          assign_first_order(m_W_val.Ax[k], W.Ax[k]);
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
        assign_first_order(m_W_val.Ax[k], W.Ax[k]);
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
        assign_first_order(m_b_val[i], b[i]);
    } else {
      for (int i = 0; i < n; ++i)
        m_b_val[i] = const_cast<F&>(b[i]).x();
    }
    m_inner.solve(m_b_val);

    // 2. Determine derivative directions from the seeded b and W.
    unsigned n_derivs;
    if constexpr (!is_ad<Inner>::value) {
      n_derivs = std::max(m_n_derivs_cached, max_deriv_size(b));
    } else {
      n_derivs = max_deriv_size(b);
      n_derivs = std::max(n_derivs, max_deriv_size(m_W_stored));
    }

    if (n_derivs == 0) {
      for (int i = 0; i < n; ++i)
        cppde::ad_traits::set_value(b[i], m_b_val[i]);
      return;
    }

    // 3. Bulk-extract ALL derivative RHS into column-major n × n_derivs
    m_rhs_all.resize(static_cast<size_t>(n) * n_derivs);
    for (int i = 0; i < n; ++i) {
      auto& bi = const_cast<F&>(b[i]);
      unsigned sz = bi.size();
      for (unsigned j = 0; j < n_derivs; ++j)
        if (j < sz) m_rhs_all[j * n + i] = bi.d(j);
        else m_rhs_all[j * n + i] = 0.0;
    }

    // 4. IFT sparse matvec: rhs_all -= dW · b_val
    ift_sparse_matvec(n, n_derivs);

    // 5. Batched solve
    m_inner.solve_batch(m_rhs_all, static_cast<int>(n_derivs));

    // 6. Bulk-inject results
    bulk_inject_results(b, m_b_val, m_rhs_all, n, n_derivs);
  }

  // W^-T b on the same factorization, the peeling transposed with the matrix.
  // See dense_lu_solver::solve_transposed in cppde_ad_lu.hpp.
  void solve_transposed(std::vector<F>& b) const
  {
    const int n = m_n;

    // 1. Extract and solve value part (reuse buffer)
    m_b_val.resize(n);
    if constexpr (cppde::ad_traits::is_dual2nd<F>::value) {
      for (int i = 0; i < n; ++i)
        assign_first_order(m_b_val[i], b[i]);
    } else {
      for (int i = 0; i < n; ++i)
        m_b_val[i] = const_cast<F&>(b[i]).x();
    }
    m_inner.solve_transposed(m_b_val);

    // 2. Determine derivative directions from the seeded b and W.
    unsigned n_derivs;
    if constexpr (!is_ad<Inner>::value) {
      n_derivs = std::max(m_n_derivs_cached, max_deriv_size(b));
    } else {
      n_derivs = max_deriv_size(b);
      n_derivs = std::max(n_derivs, max_deriv_size(m_W_stored));
    }

    if (n_derivs == 0) {
      for (int i = 0; i < n; ++i)
        cppde::ad_traits::set_value(b[i], m_b_val[i]);
      return;
    }

    // 3. Bulk-extract ALL derivative RHS into column-major n × n_derivs
    m_rhs_all.resize(static_cast<size_t>(n) * n_derivs);
    for (int i = 0; i < n; ++i) {
      auto& bi = const_cast<F&>(b[i]);
      unsigned sz = bi.size();
      for (unsigned j = 0; j < n_derivs; ++j)
        if (j < sz) m_rhs_all[j * n + i] = bi.d(j);
        else m_rhs_all[j * n + i] = 0.0;
    }

    // 4. IFT sparse matvec: rhs_all -= dW · b_val
    ift_sparse_matvec_transposed(n, n_derivs);

    // 5. One transposed solve per direction.
    {
      m_col.resize(static_cast<std::size_t>(n));
      for (unsigned j = 0; j < n_derivs; ++j) {
        for (int i = 0; i < n; ++i) m_col[i] = m_rhs_all[j * n + i];
        m_inner.solve_transposed(m_col);
        for (int i = 0; i < n; ++i) m_rhs_all[j * n + i] = m_col[i];
      }
    }

    // 6. Bulk-inject results
    bulk_inject_results(b, m_b_val, m_rhs_all, n, n_derivs);
  }

  // Batched solve for nrhs RHS vectors (column-major). Fully BLAS-3 / KLU-batched
  // at the recursion level: amortizes buffer extraction and IFT work across the
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
          assign_first_order(m_b_val_batch[col * n + i], B_flat[col * n + i]);
    } else {
      for (int col = 0; col < nrhs; ++col)
        for (int i = 0; i < n; ++i)
          m_b_val_batch[col * n + i] =
            const_cast<F&>(B_flat[col * n + i]).x();
    }

    // 2. Batched value-layer solve.
    m_inner.solve_batch(m_b_val_batch, nrhs);

    // 3. Determine n_derivs from the seeded right-hand sides and W.
    unsigned n_derivs;
    if constexpr (!is_ad<Inner>::value) {
      n_derivs = m_n_derivs_cached;
      for (int col = 0; col < nrhs; ++col)
        for (int i = 0; i < n; ++i) {
          unsigned sz = const_cast<F&>(B_flat[col * n + i]).size();
          if (sz > n_derivs) n_derivs = sz;
        }
    } else {
      n_derivs = max_deriv_size(m_W_stored);
      for (int col = 0; col < nrhs; ++col)
        for (int i = 0; i < n; ++i) {
          unsigned sz = const_cast<F&>(B_flat[col * n + i]).size();
          if (sz > n_derivs) n_derivs = sz;
        }
    }

    if (n_derivs == 0) {
      for (int col = 0; col < nrhs; ++col)
        for (int i = 0; i < n; ++i)
          cppde::ad_traits::set_value(B_flat[col * n + i], m_b_val_batch[col * n + i]);
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
          if (j < sz) m_rhs_all_batch[col * per_col + j * n + i] = bi.d(j);
          else m_rhs_all_batch[col * per_col + j * n + i] = 0.0;
      }
    }

    // 5. IFT sparse matvec batched.
    ift_sparse_matvec_batch(n, n_derivs, nrhs);

    // 6. Batched derivative-layer solve.
    m_inner.solve_batch(m_rhs_all_batch,
                        static_cast<int>(n_derivs * nrhs));

    // 7. Bulk-inject results.
    for (int col = 0; col < nrhs; ++col) {
      for (int i = 0; i < n; ++i) {
        auto& bi = B_flat[col * n + i];
        cppde::ad_traits::set_value(bi, m_b_val_batch[col * n + i]);
        if (!bi.depend()) bi.diff(0, n_derivs);
        for (unsigned j = 0; j < n_derivs; ++j)
          bi[j] = m_rhs_all_batch[col * per_col + j * n + i];
      }
    }
  }

  void reset_pattern() { m_inner.reset_pattern(); }

private:

    // ================================================================
    //  IFT sparse matvec: rhs_all -= dW * b_val.
    //
    //  For Inner = double the derivatives were pre-extracted in factorize(),
    //  laid out as m_dW_ax[p * nd + j] for the p-th CSC entry, so this is plain
    //  double arithmetic. A nested inner type takes one element-wise pass.
    // ================================================================

    // Batched IFT sparse matvec: an outer loop over the batch columns, each
    // reading a separate slice of the batch buffers, with the same inner CSC pass
    // as the single right-hand side.
  void ift_sparse_matvec_batch(int n, unsigned n_derivs, int nrhs) const
  {
    if (m_frozen) return;
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
          const Inner& b_col = b_val_col[col];
          for (int p = m_W_stored.Ap[col]; p < m_W_stored.Ap[col + 1]; ++p) {
            int row = m_W_stored.Ai[p];
            auto& w_entry = const_cast<F&>(m_W_stored.Ax[p]);
            unsigned wsz = w_entry.size();
            for (unsigned j = 0; j < n_derivs && j < wsz; ++j)
              rhs_col[j * n + row] -= w_entry.d(j) * b_col;
          }
        }
      }
    }
  }

  // The same correction with dW transposed. One CSC pass either way: the
  // transpose swaps which index reads b_val and which one is written, so the
  // pattern is walked once and never built a second time.
  void ift_sparse_matvec_transposed(int n, unsigned n_derivs) const
  {
    if (m_frozen) return;
    if constexpr (!is_ad<Inner>::value) {
      const unsigned nd_W = m_n_derivs_cached;
      if (nd_W == 0) return;
      const int* Ap = m_Ap_cached.data();
      const int* Ai = m_Ai_cached.data();
      for (int col = 0; col < n; ++col) {
        for (int p = Ap[col]; p < Ap[col + 1]; ++p) {
          const int row = Ai[p];
          const double b_row = m_b_val[row];
          const double* dw = m_dW_ax.data() + static_cast<size_t>(p) * nd_W;
          for (unsigned j = 0; j < nd_W && j < n_derivs; ++j)
            m_rhs_all[j * n + col] -= dw[j] * b_row;
        }
      }
    } else {
      for (int col = 0; col < n; ++col) {
        for (int p = m_W_stored.Ap[col]; p < m_W_stored.Ap[col + 1]; ++p) {
          const int row = m_W_stored.Ai[p];
          const Inner& b_row = m_b_val[row];
          auto& w_entry = const_cast<F&>(m_W_stored.Ax[p]);
          unsigned wsz = w_entry.size();
          for (unsigned j = 0; j < n_derivs && j < wsz; ++j)
            m_rhs_all[j * n + col] -= w_entry.d(j) * b_row;
        }
      }
    }
  }

  void ift_sparse_matvec(int n, unsigned n_derivs) const
  {
    if (m_frozen) return;
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
      //  Single pass over W_stored: each dual<Inner> element is
      //  touched exactly once.
      // =============================================================

      for (int col = 0; col < n; ++col) {
        const Inner& b_col = m_b_val[col];
        for (int p = m_W_stored.Ap[col]; p < m_W_stored.Ap[col + 1]; ++p) {
          int row = m_W_stored.Ai[p];
          auto& w_entry = const_cast<F&>(m_W_stored.Ax[p]);
          unsigned wsz = w_entry.size();
          for (unsigned j = 0; j < n_derivs && j < wsz; ++j)
            m_rhs_all[j * n + row] -= w_entry.d(j) * b_col;
        }
      }
    }
  }

  int m_n = 0;
  bool                      m_frozen = false;
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
  mutable std::vector<Inner> m_col;            // one direction of a transposed solve
  mutable std::vector<Inner> m_b_val_batch;    // n × nrhs: value part (batched)
  mutable std::vector<Inner> m_rhs_all_batch;  // n × n_derivs × nrhs (batched)
};

} // namespace ad_lu
} // namespace cppde

#endif // CPPDE_SPARSE_AD_LU_HPP
