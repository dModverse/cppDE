/**
 * @file cppde_extrapolated_flow.hpp
 * @brief The flow between two grid points and the adjoint, by extrapolated
 *        linearly implicit Euler.
 *
 * The checked sweep of a multistep model takes the flow over each interval of
 * the forward grid from the state stored at the interval start. Here that flow is the
 * extrapolation of linearly implicit Euler chains,
 *
 *   (I/h_j - J) D_k = f(x_k, t_k),   x_{k+1} = x_k + D_k,   h_j = H / n_j,
 *
 * on one matrix J for every chain and every substep, as SEULEX takes it
 * [Hairer, Wanner II, IV.9]. A chain is a W-method: it is consistent for any
 * J, so the error has an expansion in h_j and the Aitken-Neville tableau over
 * n_j = 1, 2, 3, ... raises the order by one per column. A chain costs one
 * factorization, n_j right-hand sides and n_j solves.
 *
 * J is the Jacobian at the start of the interval or of a later interval of the
 * same length, at most max_jacobian_age intervals back, so the factorization
 * of each chain serves a run of equal intervals. J is a constant of the map:
 * the extrapolated state is a fixed linear combination of the chains' end
 * states, so the adjoint is the same combination of the chains' adjoints, and
 * the tableau is built on those directly: the cotangent of the start state
 * and the interval's share of the gradient. Two neighboring entries of the
 * tableau give the error estimate the local test of CVODES' backward problem
 * reads. An interval the last column does not pass is halved.
 *
 * Several seed columns share the chains and the factorizations: the tableau
 * holds the adjoint of every column, and an interval passes when every column
 * passes.
 *
 * The tableau also extrapolates the chains' end states, so the accepted column
 * gives the map the adjoint differentiates. Its defects against the computed
 * trajectory, weighted by the cotangents, are the first-order change of the
 * seeded functional, correction().
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
  /// Intervals one Jacobian serves before the next is taken.
  static constexpr int max_jacobian_age = 10;

  /// The local error test of CVODES' backward problem: lambda to rtol, atol,
  /// the share of the gradient to rtol, gradtol when gradtol > 0. `max_pieces`
  /// bounds the pieces one interval is halved into.
  void set_refine(double rtol, double atol, double gradtol = 0.0,
                  int max_pieces = default_max_sub) {
    m_rtol = rtol; m_atol = atol; m_gradtol = gradtol; m_max_pieces = max_pieces;
  }

  /// Starts a sweep of `n_seed` seed columns.
  void begin_flow(std::size_t n_phi, std::size_t n_seed = 1) {
    m_S = n_seed;
    m_wp.assign(n_seed, std::vector<T>(n_phi, 0.0));
    m_corr.assign(n_seed, 0.0);
    m_ref_fail = 0; m_jac_valid = false;
  }
  const std::vector<T>& wp(std::size_t s = 0) const { return m_wp[s]; }
  /// The value correction of seed column s since begin_flow(): the cotangents
  /// times the defects of the accepted maps against the computed trajectory.
  double correction(std::size_t s = 0) const { return m_corr[s]; }
  /// The pieces since begin_flow() that met neither test at the finest
  /// halving, counted per seed column.
  std::size_t refine_failures() const { return m_ref_fail; }

  /// The interval [t0, t0 + dt] from `x0`, `dt` of either sign, for every seed
  /// column s: `obs[s]` the seeds observed inside it or at either end, `w[s]`
  /// the cotangent on the end state, replaced by the one on `x0`; wp(s)
  /// receives the share of the gradient. `xh_end`, the end state the run
  /// computed, may be nullptr when every `w[s]` is zero. Returns the substeps
  /// taken in the accepted chains.
  template<class System, class AdjTerms>
  int flow_interval(System& sys, const AdjTerms& adj, std::size_t n,
                    std::size_t n_phi, const std::vector<T>& x0, double t0, double dt,
                    const std::vector<const std::vector<obs_ref<T>>*>& obs,
                    const std::vector<std::vector<T>*>& w,
                    const double* xh_end = nullptr)
  {
    m_n = n; m_nphi = n_phi;
    const std::size_t S = m_S;
    const double t1 = t0 + dt;
    const double eps = 1e-12 * std::max(1.0, std::abs(t1));
    // Distance from t0 in the direction of integration.
    const double len = std::abs(dt);
    auto along = [&](double t) { return dt > 0 ? t - t0 : t0 - t; };
    // The interval is cut at the observations inside it: each cut is a point
    // a seed lands on.
    m_cut.assign(1, t0);
    for (std::size_t c = 0; c < S; ++c)
      for (const auto& o : *obs[c]) {
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
    // The cut an observation lands on.
    auto cut_of = [&](double t) -> std::size_t {
      const double s = along(t);
      if (s <= eps) return 0;
      if (s >= len - eps) return K;
      return static_cast<std::size_t>(
        std::lower_bound(m_cut.begin(), m_cut.end(), s - eps,
                         [&](double c, double v) { return along(c) < v; }) - m_cut.begin());
    };
    auto add_seeds = [&](std::size_t p, std::vector<std::vector<T>>& lam) {
      for (std::size_t c = 0; c < S; ++c)
        for (const auto& o : *obs[c])
          if (cut_of(o.t) == p)
            for (std::size_t i = 0; i < n; ++i) lam[c][i] += o.w[i];
    };
    // The start states of the pieces, forward under the state's own test.
    m_starts.resize(K);
    m_starts[0] = x0;
    for (std::size_t p = 1; p < K; ++p)
      forward(sys, m_starts[p - 1], m_cut[p - 1], m_cut[p], m_starts[p], 0);
    // The state at the interval end the pieces are compared with: the run's
    // own, or, without it, the value observed there.
    const double* s_end = xh_end;
    for (std::size_t c = 0; c < S; ++c)
      for (const auto& o : *obs[c])
        if (!s_end && o.xh && cut_of(o.t) == K) s_end = o.xh;
    auto state_at = [&](std::size_t q) -> const T* {
      return q == K ? s_end : m_starts[q].data();
    };
    // An observation the run read at cut q, against the value the pieces put there.
    for (std::size_t c = 0; c < S; ++c)
      for (const auto& o : *obs[c]) {
        const T* sq = state_at(cut_of(o.t));
        if (!o.xh || !sq) continue;
        for (std::size_t i = 0; i < n; ++i) m_corr[c] += o.w[i] * (sq[i] - o.xh[i]);
      }
    std::vector<std::vector<T>> lam(S), lam_in(S);
    for (std::size_t c = 0; c < S; ++c) lam[c] = *w[c];
    add_seeds(K, lam);
    int sub = 0;
    std::vector<T> xe;
    for (std::size_t p = K; p-- > 0;) {
      sub += backward(sys, adj, m_starts[p], m_cut[p], m_cut[p + 1], lam, lam_in, 0, xe);
      // The piece's defect against the state the next piece starts from.
      if (const T* s1 = state_at(p + 1))
        for (std::size_t c = 0; c < S; ++c)
          for (std::size_t i = 0; i < n; ++i) m_corr[c] += lam[c][i] * (xe[i] - s1[i]);
      lam.swap(lam_in);
      add_seeds(p, lam);
    }
    for (std::size_t c = 0; c < S; ++c) w[c]->swap(lam[c]);
    return sub;
  }

  /// The same for a sweep of one seed column, `w_in` the cotangent on `x0`.
  template<class System, class AdjTerms>
  int flow_interval(System& sys, const AdjTerms& adj, std::size_t n,
                    std::size_t n_phi, const std::vector<T>& x0, double t0, double dt,
                    const std::vector<obs_ref<T>>& obs,
                    const std::vector<T>& w_end, std::vector<T>& w_in,
                    const double* xh_end = nullptr)
  {
    w_in = w_end;
    return flow_interval(sys, adj, n, n_phi, x0, t0, dt, {&obs}, {&w_in}, xh_end);
  }

private:
  // The Jacobian for an interval of length H from (x, t): the one in use while
  // H stays and it is younger than max_jacobian_age, else, or on `fresh`, a new
  // one at (x, t). Returns whether it is new.
  bool linearize(const std::vector<T>& x, double t, double H, bool fresh) {
    if (!fresh && m_jac_valid && H == m_jac_H && m_jac_age < max_jacobian_age) {
      ++m_jac_age;
      return false;
    }
    m_jac_x = x; m_jac_t = t; m_jac_H = H; m_jac_age = 1;
    m_jac_valid = true;
    ++m_jac_gen;
    return true;
  }
  // The factorization of I/h - J for chain j, kept while J and h stay.
  template<class System>
  LU& factorize(System& sys, int j, double h) {
    LU& lu = m_lus[static_cast<std::size_t>(j)];
    if (m_col_gen[j] == m_jac_gen && m_col_h[j] == h) return lu;
    lu.resize(m_jac_x);
    lu.call_jacobian(sys.second, m_jac_x, m_jac_t);
    lu.factorize_W(m_n, 1.0 / h);
    m_col_gen[j] = m_jac_gen; m_col_h[j] = h;
    return lu;
  }

  // One chain of nj steps from (x, t) on the factorization `lu`, with the
  // states kept for the adjoint.
  template<class System>
  void chain(System& sys, LU& lu, const std::vector<T>& x, double t, double h, int nj) {
    if (m_X.size() < static_cast<std::size_t>(nj) + 1) m_X.resize(nj + 1);
    if (m_D.size() < static_cast<std::size_t>(nj)) m_D.resize(nj);
    m_X[0] = x;
    for (int k = 0; k < nj; ++k) {
      m_D[k].resize(m_n);
      sys.first(m_X[k], m_D[k], t + k * h);
      lu.solve(m_D[k]);
      m_X[k + 1].resize(m_n);
      for (std::size_t i = 0; i < m_n; ++i) m_X[k + 1][i] = m_X[k][i] + m_D[k][i];
    }
  }

  // The chain's adjoint from `lam_end` on the chain end state, on the
  // factorization `lu`: into the zeroed block `a`, the cotangent on the start
  // state, then the share of the gradient.
  template<class AdjTerms>
  void chain_adjoint(const AdjTerms& adj, LU& lu, double t, double h, int nj,
                     const std::vector<T>& lam_end, T* a) {
    const std::size_t n = m_n;
    // `share` is indexed as the rows of the gradient, states first.
    T* share = a + n;
    m_lam = lam_end;
    for (int k = nj; k-- > 0;) {
      m_mu = m_lam;
      lu.solve_transposed(m_mu);
      rhs_pullback(adj, m_X[k], m_mu, t + k * h, 1.0, m_jv, share);
      for (std::size_t i = 0; i < n; ++i) m_lam[i] += m_jv[i];
    }
    std::copy(m_lam.begin(), m_lam.end(), a);
  }

  // Aitken-Neville: row j of the tableau from the first entry `a` and row
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

  // The weighted root-mean-square norm of the difference of two adjoints of
  // seed column c: lambda against rtol |lambda| + atol, and under gradtol the
  // share against rtol |running gradient| + gradtol; the larger of the two.
  double adjoint_error(const std::vector<T>& ta, const std::vector<T>& tb, std::size_t c) const {
    const std::size_t n = m_n;
    const T* a = ta.data() + c * (n + m_nphi);
    const T* b = tb.data() + c * (n + m_nphi);
    const std::vector<T>& wp = m_wp[c];
    double sl = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
      const double w = m_rtol * std::abs(a[i]) + m_atol;
      if (w > 0.0) { const double r = (a[i] - b[i]) / w; sl += r * r; }
    }
    double e = n > 0 ? std::sqrt(sl / static_cast<double>(n)) : 0.0;
    if (m_gradtol > 0.0 && m_nphi > n) {
      double sq = 0.0;
      for (std::size_t i = n; i < m_nphi; ++i) {
        const double q = wp[i] + a[n + i];
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

  // The adjoint of the flow over [a, b] from `x` for every seed column c: the
  // cotangent `lam_out[c]` on the end state to `lam_in[c]` on `x`, the share
  // added to wp(c), and in `xe` the end state of the accepted map, the defects
  // of inner halvings added to the correction. Returns the substeps of the
  // chains it accepted.
  template<class System, class AdjTerms>
  int backward(System& sys, const AdjTerms& adj, const std::vector<T>& x, double a,
               double b, const std::vector<std::vector<T>>& lam_out,
               std::vector<std::vector<T>>& lam_in, int depth, std::vector<T>& xe)
  {
    const std::size_t S = m_S, n = m_n, blk = m_n + m_nphi;
    const double H = b - a;
    if (!(std::abs(H) > 0.0)) { lam_in = lam_out; xe = x; return 0; }
    const bool can_split = depth < max_depth && (2 << depth) <= m_max_pieces;
    int used = 0, j = 0;
    bool ok = false;
    m_e_prev.resize(S);
    m_est.resize(S);
    auto columns = [&]() {
      std::fill(m_e_prev.begin(), m_e_prev.end(), 0.0);
      ok = false;
      used = 0;
      for (j = 0; j < max_columns; ++j) {
        const int nj = seq(j);
        const double h = H / nj;
        LU& lu = factorize(sys, j, h);
        chain(sys, lu, x, a, h, nj);
        // The adjoints of every column, then the chain's end state, extrapolated with them.
        m_a.assign(S * blk + n, 0.0);
        for (std::size_t c = 0; c < S; ++c)
          chain_adjoint(adj, lu, a, h, nj, lam_out[c], m_a.data() + c * blk);
        std::copy(m_X[nj].begin(), m_X[nj].end(), m_a.begin() + static_cast<std::ptrdiff_t>(S * blk));
        used += nj;
        extrapolate(j, m_a);
        if (j < 1) continue;
        // The difference estimates the error of T[j][j-1]; T[j][j] is returned,
        // the error of T[j][j] is the rate of the last two differences applied one column on.
        ok = true;
        for (std::size_t c = 0; c < S; ++c) {
          const double e = adjoint_error(m_tab[j][j], m_tab[j][j - 1], c);
          m_est[c] = (j >= 2 && m_e_prev[c] > 0.0) ? e * std::min(1.0, e / m_e_prev[c]) : e;
          m_e_prev[c] = e;
          ok = ok && m_est[c] <= 1.0;
        }
        if (ok) return;
      }
    };
    const bool kept = !linearize(x, a, H, false);
    columns();
    // A Jacobian kept from another interval is replaced before the interval is halved.
    if (!ok && kept) { linearize(x, a, H, true); columns(); }
    if (!ok && can_split) {
      // Halved: the right half first, from the midpoint the state's own test
      // reaches, then the left half on the cotangent it hands down.
      const double mid = a + 0.5 * H;
      std::vector<T> xm, xe_left;
      std::vector<std::vector<T>> lam_mid(S);
      forward(sys, x, a, mid, xm, depth + 1);
      int sub = backward(sys, adj, xm, mid, b, lam_out, lam_mid, depth + 1, xe);
      sub += backward(sys, adj, x, a, mid, lam_mid, lam_in, depth + 1, xe_left);
      // The right half was linearized at xm, the left half ends elsewhere.
      for (std::size_t c = 0; c < S; ++c)
        for (std::size_t i = 0; i < n; ++i) m_corr[c] += lam_mid[c][i] * (xe_left[i] - xm[i]);
      return sub;
    }
    if (!ok) {
      for (std::size_t c = 0; c < S; ++c) m_ref_fail += m_est[c] > 1.0;
      j = std::min(j, max_columns - 1);
    }
    const auto& acc = m_tab[j][j];
    lam_in.resize(S);
    for (std::size_t c = 0; c < S; ++c) {
      const auto b0 = acc.begin() + static_cast<std::ptrdiff_t>(c * blk);
      lam_in[c].assign(b0, b0 + static_cast<std::ptrdiff_t>(n));
      for (std::size_t i = n; i < m_nphi; ++i) m_wp[c][i] += acc[c * blk + n + i];
    }
    const auto e0 = acc.begin() + static_cast<std::ptrdiff_t>(S * blk);
    xe.assign(e0, e0 + static_cast<std::ptrdiff_t>(n));
    return used;
  }

  // The flow's end state over [a, b] from `x`, under the state's own test.
  template<class System>
  void forward(System& sys, const std::vector<T>& x, double a, double b,
               std::vector<T>& out, int depth)
  {
    const double H = b - a;
    if (!(std::abs(H) > 0.0)) { out = x; return; }
    const bool can_split = depth < max_depth && (2 << depth) <= m_max_pieces;
    bool ok = false;
    int j = 0;
    auto columns = [&]() {
      ok = false;
      for (j = 0; j < max_columns; ++j) {
        const int nj = seq(j);
        LU& lu = factorize(sys, j, H / nj);
        chain(sys, lu, x, a, H / nj, nj);
        m_a = m_X[nj];
        extrapolate(j, m_a);
        if (j >= 1 && state_error(m_tab[j][j], m_tab[j][j - 1]) <= 1.0) { ok = true; return; }
      }
    };
    const bool kept = !linearize(x, a, H, false);
    columns();
    if (!ok && kept) { linearize(x, a, H, true); columns(); }
    if (!ok && can_split) {
      const double mid = a + 0.5 * H;
      std::vector<T> xm;
      forward(sys, x, a, mid, xm, depth + 1);
      forward(sys, xm, mid, b, out, depth + 1);
      return;
    }
    if (!ok) { m_ref_fail += m_S; j = max_columns - 1; }
    out = m_tab[j][j];
  }

  // One factorization per column, each with the Jacobian generation and the h
  // it was built for.
  LU m_lus[max_columns];
  int m_col_gen[max_columns] = {};
  double m_col_h[max_columns] = {};
  std::vector<T> m_jac_x;
  double m_jac_t = 0.0, m_jac_H = 0.0;
  int m_jac_gen = 0, m_jac_age = 0;
  bool m_jac_valid = false;
  double m_rtol = 0.0, m_atol = 0.0, m_gradtol = 0.0;
  int m_max_pieces = default_max_sub;
  std::size_t m_S = 1, m_n = 0, m_nphi = 0, m_ref_fail = 0;
  std::vector<double> m_corr, m_e_prev, m_est;
  std::vector<std::vector<T>> m_wp;
  std::vector<T> m_lam, m_mu, m_jv, m_a;
  std::vector<double> m_cut;
  std::vector<std::vector<T>> m_X, m_D, m_starts;
  std::vector<std::vector<std::vector<T>>> m_tab;
};

}  // namespace adjoint
}  // namespace cppde

#endif  // CPPDE_EXTRAPOLATED_FLOW_HPP
