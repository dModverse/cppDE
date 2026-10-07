/**
 * @file cppde_extrapolated_flow.hpp
 * @brief The flow between two grid points and its adjoint, by extrapolated
 *        linearly implicit Euler.
 *
 * The checked sweep of a multistep model takes the flow over each interval of
 * the forward grid from the state stored at its start. Here that flow is the
 * extrapolation of linearly implicit Euler chains,
 *
 *   (I/h_j - J) D_k = f(x_k, t_k),   x_{k+1} = x_k + D_k,   h_j = H / n_j,
 *
 * on one Jacobian J at the interval's start for every chain and every
 * substep, as SEULEX takes it [Hairer, Wanner II, IV.9]. A chain is a
 * W-method: it is consistent for any J, so its error has an expansion in h_j
 * and the Aitken-Neville tableau over n_j = 1, 2, 3, ... raises the order by
 * one per column. A chain costs one factorisation, n_j right-hand sides and
 * n_j solves; the factorisations of all chains share the one Jacobian.
 *
 * The extrapolated state is a fixed linear combination of the chains' end
 * states, so its adjoint is the same combination of the chains' adjoints, and
 * the tableau is built on those directly: the cotangent of the start state
 * and the interval's share of the gradient. Two neighbouring entries of the
 * tableau give the error estimate the local test of CVODES' backward problem
 * reads. An interval the last column does not pass is halved.
 *
 * Copyright (C) 2026 Simon Beyer
 */

#ifndef CPPDE_EXTRAPOLATED_FLOW_HPP
#define CPPDE_EXTRAPOLATED_FLOW_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include <cppde/cppde_lu.hpp>
#include <cppde/cppde_adjoint_step.hpp>

namespace cppde {
namespace adjoint {

template<class LU>
class extrapolated_flow {
public:
  using T = double;
  static constexpr int max_columns = 8;
  static constexpr int max_depth = 10;

  /// The local error test of CVODES' backward problem: lambda to rtol, atol,
  /// the share of the gradient to rtol, gradtol when gradtol > 0. `max_pieces`
  /// bounds the pieces one interval is halved into.
  void set_refine(double rtol, double atol, double gradtol = 0.0,
                  int max_pieces = default_max_sub) {
    m_rtol = rtol; m_atol = atol; m_gradtol = gradtol; m_max_pieces = max_pieces;
  }

  void begin_flow(std::size_t n_phi) { m_wp.assign(n_phi, 0.0); m_ref_fail = 0; }
  const std::vector<T>& wp() const { return m_wp; }
  /// The pieces since begin_flow() that met neither test at the finest halving.
  std::size_t refine_failures() const { return m_ref_fail; }

  /// The interval [t0, t0 + dt] from `x0`, `dt` of either sign: `w_end` the
  /// cotangent on its end state, `obs` the seeds observed inside it or at
  /// either end; `w_in` receives the cotangent on `x0`, wp() the share of the
  /// gradient. Returns the substeps taken in the accepted chains.
  template<class System, class AdjTerms>
  int flow_interval(System& sys, const AdjTerms& adj, std::size_t n,
                    std::size_t n_phi, const std::vector<T>& x0, double t0, double dt,
                    const std::vector<obs_ref<T>>& obs,
                    const std::vector<T>& w_end, std::vector<T>& w_in)
  {
    m_n = n; m_nphi = n_phi;
    const double t1 = t0 + dt;
    const double eps = 1e-12 * std::max(1.0, std::abs(t1));
    // Distance from t0 in the direction of integration.
    const double len = std::abs(dt);
    auto along = [&](double t) { return dt > 0 ? t - t0 : t0 - t; };
    // The interval is cut at the observations inside it: each cut is a point
    // a seed lands on.
    m_cut.assign(1, t0);
    for (const auto& o : obs) {
      const double s = along(o.t);
      if (s > eps && s < len - eps) m_cut.push_back(o.t);
    }
    std::sort(m_cut.begin() + 1, m_cut.end(),
              [&](double a, double b) { return along(a) < along(b); });
    m_cut.erase(std::unique(m_cut.begin(), m_cut.end(),
                            [&](double a, double b) { return std::abs(a - b) <= eps; }),
                m_cut.end());
    m_cut.push_back(t1);
    const std::size_t K = m_cut.size() - 1;
    auto add_seeds = [&](std::size_t p, std::vector<T>& lam) {
      for (const auto& o : obs) {
        const double s = along(o.t);
        std::size_t q;
        if (s <= eps) q = 0;
        else if (s >= len - eps) q = K;
        else q = static_cast<std::size_t>(
          std::lower_bound(m_cut.begin(), m_cut.end(), s - eps,
                           [&](double c, double v) { return along(c) < v; }) - m_cut.begin());
        if (q == p) for (std::size_t i = 0; i < n; ++i) lam[i] += o.w[i];
      }
    };
    // The start states of the pieces, forward under the state's own test.
    m_starts.resize(K);
    m_starts[0] = x0;
    for (std::size_t p = 1; p < K; ++p)
      forward(sys, m_starts[p - 1], m_cut[p - 1], m_cut[p], m_starts[p], 0);
    std::vector<T> lam = w_end;
    add_seeds(K, lam);
    int sub = 0;
    for (std::size_t p = K; p-- > 0;) {
      std::vector<T> lam_in;
      sub += backward(sys, adj, m_starts[p], m_cut[p], m_cut[p + 1], lam, lam_in, 0);
      lam.swap(lam_in);
      add_seeds(p, lam);
    }
    w_in = lam;
    return sub;
  }

private:
  // The Jacobian at (x, t), kept for the factorisations of every chain.
  template<class System>
  void linearise(System& sys, const std::vector<T>& x, double t) {
    m_lu.resize(x);
    m_lu.call_jacobian(sys.second, x, t);
    m_lu.cache_jacobian(m_n);
    m_fresh = true;
  }
  // The chains of one Jacobian come in shrinking steps, so after the first
  // only the diagonal grows; backwards in time it shrinks and is factorised anew.
  void factorise(double h) {
    if (m_fresh || h < 0.0) m_lu.factorize_W(m_n, 1.0 / h);
    else m_lu.refactorize_W_larger_diagonal(m_n, 1.0 / h);
    m_fresh = false;
  }

  // One chain of nj steps from (x, t) on the current factorisation, its
  // states and increments kept for the adjoint.
  template<class System>
  void chain(System& sys, const std::vector<T>& x, double t, double h, int nj) {
    if (m_X.size() < static_cast<std::size_t>(nj) + 1) m_X.resize(nj + 1);
    if (m_D.size() < static_cast<std::size_t>(nj)) m_D.resize(nj);
    m_X[0] = x;
    for (int k = 0; k < nj; ++k) {
      m_D[k].resize(m_n);
      sys.first(m_X[k], m_D[k], t + k * h);
      m_lu.solve(m_D[k]);
      m_X[k + 1].resize(m_n);
      for (std::size_t i = 0; i < m_n; ++i) m_X[k + 1][i] = m_X[k][i] + m_D[k][i];
    }
  }

  // The chain's adjoint from `lam_end` on its end state: into `a`, the cotangent
  // on the start state, then the share of the gradient. The Jacobian's own
  // dependence on state and theta is one pullback of sum_k mu_k' (dJ) D_k.
  template<class AdjTerms>
  void chain_adjoint(const AdjTerms& adj, const std::vector<T>& x, double t, double h,
                     int nj, const std::vector<T>& lam_end, std::vector<T>& a) {
    const std::size_t n = m_n;
    // `share` is indexed as the rows of the gradient, states first.
    a.assign(n + m_nphi, 0.0);
    T* share = a.data() + n;
    m_lam = lam_end;
    m_wJ.assign(n, 0.0);
    if constexpr (has_jvp_bundle<AdjTerms>::value)
      m_M.assign(AdjTerms::jvp_bundle_size, 0.0);
    for (int k = nj; k-- > 0;) {
      m_mu = m_lam;
      m_lu.solve_transposed(m_mu);
      const double tk = t + k * h;
      rhs_pullback(adj, m_X[k], m_mu, tk, 1.0, m_jv, share);
      for (std::size_t i = 0; i < n; ++i) m_lam[i] += m_jv[i];
      if constexpr (has_jvp_bundle<AdjTerms>::value) {
        const int* r = AdjTerms::jvp_bundle_rows();
        const int* c = AdjTerms::jvp_bundle_cols();
        const std::vector<T>& D = m_D[k];
        for (int q = 0; q < AdjTerms::jvp_bundle_size; ++q) m_M[q] += m_mu[r[q]] * D[c[q]];
      } else {
        adj.jvp_x_t_vec(x, m_D[k], m_mu, t, m_jv);
        for (std::size_t i = 0; i < n; ++i) m_wJ[i] += m_jv[i];
        adj.jvp_p_t_vec_axpy(x, m_D[k], m_mu, t, 1.0, share);
      }
    }
    if constexpr (has_jvp_bundle<AdjTerms>::value) {
      bundle_pullback(adj, x, m_M, t, m_jv, share);
      for (std::size_t i = 0; i < n; ++i) m_wJ[i] += m_jv[i];
    }
    for (std::size_t i = 0; i < n; ++i) a[i] = m_lam[i] + m_wJ[i];
  }

  // Aitken-Neville: row j of the tableau from its first entry `a` and row
  // j - 1, entry k of order k + 1.
  void extrapolate(int j, std::vector<T>& a) {
    if (m_tab.size() < static_cast<std::size_t>(max_columns)) m_tab.resize(max_columns);
    auto& row = m_tab[j];
    if (row.size() < static_cast<std::size_t>(j) + 1) row.resize(j + 1);
    row[0].swap(a);
    const std::size_t sz = row[0].size();
    for (int k = 1; k <= j; ++k) {
      const double inv_r = 1.0 / (static_cast<double>(seq(j)) / seq(j - k) - 1.0);
      row[k].resize(sz);
      const T* CPPDE_RESTRICT up = row[k - 1].data();
      const T* CPPDE_RESTRICT lo = m_tab[j - 1][k - 1].data();
      T* CPPDE_RESTRICT out = row[k].data();
      for (std::size_t i = 0; i < sz; ++i) out[i] = up[i] + (up[i] - lo[i]) * inv_r;
    }
  }
  static int seq(int j) { return j + 1; }

  // The weighted root-mean-square norm of the difference of two adjoints:
  // lambda against rtol |lambda| + atol, and under gradtol the share against
  // rtol |running gradient| + gradtol; the larger of the two.
  double adjoint_error(const std::vector<T>& a, const std::vector<T>& b) const {
    const std::size_t n = m_n;
    double sl = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
      const double w = m_rtol * std::abs(a[i]) + m_atol;
      if (w > 0.0) { const double r = (a[i] - b[i]) / w; sl += r * r; }
    }
    double e = n > 0 ? std::sqrt(sl / static_cast<double>(n)) : 0.0;
    if (m_gradtol > 0.0 && m_nphi > n) {
      double sq = 0.0;
      for (std::size_t i = n; i < m_nphi; ++i) {
        const double q = m_wp[i] + a[n + i];
        const double r = (a[n + i] - b[n + i]) / (m_rtol * std::abs(q) + m_gradtol);
        sq += r * r;
      }
      e = std::max(e, std::sqrt(sq / static_cast<double>(m_nphi - n)));
    }
    return e;
  }
  double state_error(const std::vector<T>& a, const std::vector<T>& b) const {
    double s = 0.0;
    for (std::size_t i = 0; i < m_n; ++i) {
      const double r = (a[i] - b[i]) / (m_rtol * std::abs(a[i]) + m_atol);
      s += r * r;
    }
    return m_n > 0 ? std::sqrt(s / static_cast<double>(m_n)) : 0.0;
  }

  // The adjoint of the flow over [a, b] from `x`: the cotangent `lam_out` on
  // the end state to `lam_in` on `x`, the share added to wp(). Returns the
  // substeps of the chains it accepted.
  template<class System, class AdjTerms>
  int backward(System& sys, const AdjTerms& adj, const std::vector<T>& x, double a,
               double b, const std::vector<T>& lam_out, std::vector<T>& lam_in, int depth)
  {
    const double H = b - a;
    if (!(std::abs(H) > 0.0)) { lam_in = lam_out; return 0; }
    linearise(sys, x, a);
    const bool can_split = depth < max_depth && (2 << depth) <= m_max_pieces;
    int used = 0, j = 0;
    bool ok = false;
    double e_prev = 0.0;
    for (; j < max_columns; ++j) {
      const int nj = seq(j);
      const double h = H / nj;
      factorise(h);
      chain(sys, x, a, h, nj);
      chain_adjoint(adj, x, a, h, nj, lam_out, m_a);
      used += nj;
      extrapolate(j, m_a);
      if (j < 1) continue;
      // The difference estimates the error of T[j][j-1]; T[j][j] is returned,
      // its error is the rate of the last two differences applied one column on.
      const double e = adjoint_error(m_tab[j][j], m_tab[j][j - 1]);
      const double est = (j >= 2 && e_prev > 0.0) ? e * std::min(1.0, e / e_prev) : e;
      if (est <= 1.0) { ok = true; break; }
      e_prev = e;
    }
    if (!ok && can_split) {
      // Halved: the right half first, from the midpoint the state's own test
      // reaches, then the left half on the cotangent it hands down.
      const double mid = a + 0.5 * H;
      std::vector<T> xm, lam_mid;
      forward(sys, x, a, mid, xm, depth + 1);
      int s = backward(sys, adj, xm, mid, b, lam_out, lam_mid, depth + 1);
      s += backward(sys, adj, x, a, mid, lam_mid, lam_in, depth + 1);
      return s;
    }
    if (!ok) { ++m_ref_fail; j = std::min(j, max_columns - 1); }
    const auto& acc = m_tab[j][j];
    lam_in.assign(acc.begin(), acc.begin() + static_cast<std::ptrdiff_t>(m_n));
    for (std::size_t i = m_n; i < m_nphi; ++i) m_wp[i] += acc[m_n + i];
    return used;
  }

  // The flow's end state over [a, b] from `x`, under the state's own test.
  template<class System>
  void forward(System& sys, const std::vector<T>& x, double a, double b,
               std::vector<T>& out, int depth)
  {
    const double H = b - a;
    if (!(std::abs(H) > 0.0)) { out = x; return; }
    linearise(sys, x, a);
    const bool can_split = depth < max_depth && (2 << depth) <= m_max_pieces;
    bool ok = false;
    int j = 0;
    for (; j < max_columns; ++j) {
      const int nj = seq(j);
      factorise(H / nj);
      chain(sys, x, a, H / nj, nj);
      m_a = m_X[nj];
      extrapolate(j, m_a);
      if (j >= 1 && state_error(m_tab[j][j], m_tab[j][j - 1]) <= 1.0) { ok = true; break; }
    }
    if (!ok && can_split) {
      const double mid = a + 0.5 * H;
      std::vector<T> xm;
      forward(sys, x, a, mid, xm, depth + 1);
      forward(sys, xm, mid, b, out, depth + 1);
      return;
    }
    if (!ok) { ++m_ref_fail; j = max_columns - 1; }
    out = m_tab[j][j];
  }

  LU m_lu;
  bool m_fresh = false;
  double m_rtol = 0.0, m_atol = 0.0, m_gradtol = 0.0;
  int m_max_pieces = default_max_sub;
  std::size_t m_n = 0, m_nphi = 0, m_ref_fail = 0;
  std::vector<T> m_wp, m_lam, m_mu, m_jv, m_wJ, m_M, m_a;
  std::vector<double> m_cut;
  std::vector<std::vector<T>> m_X, m_D, m_starts;
  std::vector<std::vector<std::vector<T>>> m_tab;
};

}  // namespace adjoint
}  // namespace cppde

#endif  // CPPDE_EXTRAPOLATED_FLOW_HPP
