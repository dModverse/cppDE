/**
 * @file cppde_sliding_adjoint.hpp
 * @brief The model's adjoint terms where a switch slides (cppde_sliding.hpp).
 *
 * While the switches S slide, f = sum_c w_c(a) f_c with the rates
 * R_k(a) = sum_c w_c(a) r'_{k,c} at zero, so grad a = -R_a^{-1} R_x and the
 * pullback of lambda is
 *
 *   J' lambda = sum_c w_c J_c' lambda - sum_k u_k sum_c w_c grad r'_{k,c},
 *   u = R_a^{-T} q,  q_j = df/da_j . lambda,
 *
 * and the same in the parameters. grad r'_{k,c} is the generated gradient of
 * root k's rate with the slides held at combination c. The terms of
 * the second order the Rosenbrock and refined sweeps read are central
 * differences of these pullbacks.
 *
 * Copyright (C) 2026 Simon Beyer
 */

#ifndef CPPDE_SLIDING_ADJOINT_HPP
#define CPPDE_SLIDING_ADJOINT_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include <cppde/cppde_jump_adjoint.hpp>
#include <cppde/cppde_sliding.hpp>

namespace cppde {
namespace adjoint {

template<class Adj, class EvAdj, class Slider>
struct sliding_adjoint {
  using core_type = typename Slider::core_type;
  using T = typename core_type::V;
  using Vec = std::vector<T>;

  const Adj& adj;
  const EvAdj& eadj;
  core_type* c;

  // The slides S at (x, t): their weights a, the rates r'_k under every
  // combination c, f of each combination, and R_a, the derivative of the
  // rates' mix in a.
  struct mix {
    std::vector<int> S;
    std::vector<T> a, Ra;
    std::vector<std::vector<T>> r;
    std::vector<Vec> fm;
  };

  bool prepare(const Vec& x, const T& t, mix& m) const {
    m.S = cppde::detail::sliding_switches();
    if (m.S.empty()) return false;
    if (c->inexact)
      throw std::runtime_error(
          "a switch slides along its surface, which the reverse mode of a "
          "sparse model does not support");
    if (!c->solve(m.S, x, t, m.a, &m.r)) core_type::no_mix();
    m.fm = c->fm;
    const std::size_t n = m.S.size();
    m.Ra.assign(n * n, T(0.0));
    for (std::size_t k = 0; k < n; ++k)
      for (std::size_t cc = 0; cc < m.fm.size(); ++cc)
        for (std::size_t j = 0; j < n; ++j)
          m.Ra[k * n + j] = m.Ra[k * n + j] +
                            m.r[k][cc] * cppde::detail::combo_weight(m.a, cc, static_cast<int>(j));
    return true;
  }

  // `f` evaluated with the slides held at combination cc.
  template<class F>
  static void held(const std::vector<int>& S, std::size_t cc, F&& f) {
    signed char* md = switch_mode_sink().m;
    for (std::size_t q = 0; q < S.size(); ++q) md[S[q]] = ((cc >> q) & 1u) ? 1 : 0;
    f();
    for (int k : S) md[k] = cppde::sliding_mode;
  }

  static T dot(const Vec& u, const Vec& v) {
    T s(0.0);
    for (std::size_t i = 0; i < u.size(); ++i) s += u[i] * v[i];
    return s;
  }

  // u = R_a^{-T} q, q_j = df/da_j . lambda: the multipliers of grad a.
  static std::vector<T> multipliers(const mix& m, const Vec& lam) {
    const std::size_t n = m.S.size();
    std::vector<T> q(n, T(0.0)), At(n * n);
    for (std::size_t j = 0; j < n; ++j)
      for (std::size_t cc = 0; cc < m.fm.size(); ++cc)
        q[j] = q[j] + cppde::detail::combo_weight(m.a, cc, static_cast<int>(j)) * dot(m.fm[cc], lam);
    for (std::size_t k = 0; k < n; ++k)
      for (std::size_t j = 0; j < n; ++j) At[j * n + k] = m.Ra[k * n + j];
    cppde::detail::small_solve(At, q, [](const auto& v) { return cppde::detail::scalar_value(v); });
    return q;
  }

  // J' lambda = sum_c w_c J_c' lambda - sum_k u_k sum_c w_c grad r'_{k,c}.
  void jac_t_vec(const Vec& x, const Vec& lam, const T& t, Vec& out) const {
    mix m;
    if (!prepare(x, t, m)) {
      adj.jac_t_vec(x, lam, t, out);
      return;
    }
    const std::vector<T> u = multipliers(m, lam);
    out.assign(x.size(), T(0.0));
    // Armed, so the generated callees write their tangents into these.
    Vec jc, gk;
    for (std::size_t cc = 0; cc < m.fm.size(); ++cc) {
      const T w = cppde::detail::combo_weight(m.a, cc);
      held(m.S, cc, [&] {
        zero_armed(jc, x.size());
        adj.jac_t_vec(x, lam, t, jc);
        for (std::size_t i = 0; i < x.size(); ++i) out[i] = out[i] + w * jc[i];
        for (std::size_t k = 0; k < m.S.size(); ++k) {
          zero_armed(gk, x.size());
          eadj.root_gdot_dx(c->root_of[m.S[k]], x, t, gk);
          for (std::size_t i = 0; i < x.size(); ++i) out[i] = out[i] - u[k] * w * gk[i];
        }
      });
    }
  }

  void dfdp_t_vec_axpy(const Vec& x, const Vec& lam, const T& t, const T& sc, T* out) const {
    mix m;
    if (!prepare(x, t, m)) {
      adj.dfdp_t_vec_axpy(x, lam, t, sc, out);
      return;
    }
    const std::vector<T> u = multipliers(m, lam);
    for (std::size_t cc = 0; cc < m.fm.size(); ++cc) {
      const T w = cppde::detail::combo_weight(m.a, cc);
      held(m.S, cc, [&] {
        adj.dfdp_t_vec_axpy(x, lam, t, sc * w, out);
        for (std::size_t k = 0; k < m.S.size(); ++k)
          eadj.root_gdot_dp_axpy(c->root_of[m.S[k]], x, t, T(0.0) - sc * u[k] * w, out);
      });
    }
  }

  void vjp_t_vec(const Vec& x, const Vec& lam, const T& t, const T& sc, Vec& out,
                 T* outp) const {
    if (cppde::detail::sliding_switch() < 0) {
      adj.vjp_t_vec(x, lam, t, sc, out, outp);
      return;
    }
    jac_t_vec(x, lam, t, out);
    dfdp_t_vec_axpy(x, lam, t, sc, outp);
  }

  // df/dt . lambda; da/dt by a central difference of the weights at fixed x.
  T dfdt_dot(const Vec& x, const Vec& lam, const T& t) const {
    mix m;
    if (!prepare(x, t, m)) return adj.dfdt_dot(x, lam, t);
    T out(0.0);
    for (std::size_t cc = 0; cc < m.fm.size(); ++cc) {
      const T w = cppde::detail::combo_weight(m.a, cc);
      held(m.S, cc, [&] { out = out + w * adj.dfdt_dot(x, lam, t); });
    }
    const double h = 1e-6 * std::max(1.0, std::abs(cppde::detail::scalar_value(t)));
    std::vector<T> up, dn;
    if (!c->solve(m.S, x, t + T(h), up) || !c->solve(m.S, x, t - T(h), dn))
      core_type::no_mix();
    for (std::size_t j = 0; j < m.S.size(); ++j) {
      T q(0.0);
      for (std::size_t cc = 0; cc < m.fm.size(); ++cc)
        q = q + cppde::detail::combo_weight(m.a, cc, static_cast<int>(j)) * dot(m.fm[cc], lam);
      out = out + q * (up[j] - dn[j]) / T(2.0 * h);
    }
    return out;
  }

  // The terms of the second order where a switch slides: central differences
  // of the exact pullbacks above along v, or in t.
  static double step(const Vec& x, const Vec& v) {
    double n = 0.0, s = 1.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
      n = std::max(n, std::abs(cppde::detail::scalar_value(v[i])));
      s = std::max(s, std::abs(cppde::detail::scalar_value(x[i])));
    }
    return n == 0.0 ? 0.0 : 1e-6 * s / n;
  }

  static Vec shifted(const Vec& x, const Vec& v, double h) {
    Vec y(x);
    for (std::size_t i = 0; i < x.size(); ++i) y[i] = x[i] + T(h) * v[i];
    return y;
  }

  static double tstep(const T& t) {
    return 1e-6 * std::max(1.0, std::abs(cppde::detail::scalar_value(t)));
  }

  void jvp_x_t_vec(const Vec& x, const Vec& v, const Vec& lam, const T& t, Vec& out) const {
    if (cppde::detail::sliding_switch() < 0) {
      adj.jvp_x_t_vec(x, v, lam, t, out);
      return;
    }
    const double h = step(x, v);
    out.assign(x.size(), T(0.0));
    if (h == 0.0) return;
    Vec up, dn;
    jac_t_vec(shifted(x, v, h), lam, t, up);
    jac_t_vec(shifted(x, v, -h), lam, t, dn);
    for (std::size_t i = 0; i < x.size(); ++i) out[i] = (up[i] - dn[i]) / T(2.0 * h);
  }

  void jvp_p_t_vec_axpy(const Vec& x, const Vec& v, const Vec& lam, const T& t, const T& sc,
                        T* out) const {
    if (cppde::detail::sliding_switch() < 0) {
      adj.jvp_p_t_vec_axpy(x, v, lam, t, sc, out);
      return;
    }
    const double h = step(x, v);
    if (h == 0.0) return;
    dfdp_t_vec_axpy(shifted(x, v, h), lam, t, sc / T(2.0 * h), out);
    dfdp_t_vec_axpy(shifted(x, v, -h), lam, t, T(0.0) - sc / T(2.0 * h), out);
  }

  void dfdt_x_t_vec(const Vec& x, const Vec& lam, const T& t, Vec& out) const {
    if (cppde::detail::sliding_switch() < 0) {
      adj.dfdt_x_t_vec(x, lam, t, out);
      return;
    }
    const double h = tstep(t);
    Vec up, dn;
    jac_t_vec(x, lam, t + T(h), up);
    jac_t_vec(x, lam, t - T(h), dn);
    out.assign(x.size(), T(0.0));
    for (std::size_t i = 0; i < x.size(); ++i) out[i] = (up[i] - dn[i]) / T(2.0 * h);
  }

  void dfdt_p_t_vec_axpy(const Vec& x, const Vec& lam, const T& t, const T& sc, T* out) const {
    if (cppde::detail::sliding_switch() < 0) {
      adj.dfdt_p_t_vec_axpy(x, lam, t, sc, out);
      return;
    }
    const double h = tstep(t);
    dfdp_t_vec_axpy(x, lam, t + T(h), sc / T(2.0 * h), out);
    dfdp_t_vec_axpy(x, lam, t - T(h), T(0.0) - sc / T(2.0 * h), out);
  }

  // The fused matrix forms are not offered, so the sweeps take the vector forms.
};

// The event terms of a model where switches slide. A root read while others
// slide has the rate of the mix, sum_c w_c r'_{i,c}, whose gradient adds
// sum_j dR_i/da_j grad a_j; an exit root, a_k - 1 or -a_k, has the parameter
// derivative of a_k. Its rate's gradients enter only scaled by a quantity zero
// at an exit, where f is continuous, so they are left zero.
template<class EvAdj, class Slider>
struct sliding_event_adjoint {
  using core_type = typename Slider::core_type;
  using T = typename core_type::V;
  using Vec = std::vector<T>;

  const EvAdj& eadj;
  core_type* c;

  int exit_of(int ev) const { return c->roots[ev].slide_exit; }

  // The slides, their weights, rates and R_a^{-1} at (x, t); false without any.
  bool prepare(const Vec& x, const T& t, std::vector<int>& S, std::vector<T>& a,
               std::vector<std::vector<T>>& r, std::vector<T>& Rinv) const {
    S = cppde::detail::sliding_switches();
    if (S.empty()) return false;
    if (!c->solve(S, x, t, a, &r)) core_type::no_mix();
    const std::size_t n = S.size(), M = std::size_t(1) << n;
    std::vector<T> Ra(n * n, T(0.0));
    for (std::size_t k = 0; k < n; ++k)
      for (std::size_t cc = 0; cc < M; ++cc)
        for (std::size_t j = 0; j < n; ++j)
          Ra[k * n + j] = Ra[k * n + j] +
                          r[k][cc] * cppde::detail::combo_weight(a, cc, static_cast<int>(j));
    Rinv.assign(n * n, T(0.0));
    for (std::size_t j = 0; j < n; ++j) {
      std::vector<T> e(n, T(0.0));
      e[j] = T(1.0);
      cppde::detail::small_solve(Ra, e, [](const auto& v) { return cppde::detail::scalar_value(v); });
      for (std::size_t i = 0; i < n; ++i) Rinv[i * n + j] = e[i];
    }
    return true;
  }

  template<class F>
  static void held(const std::vector<int>& S, std::size_t cc, F&& f) {
    signed char* md = switch_mode_sink().m;
    for (std::size_t q = 0; q < S.size(); ++q) md[S[q]] = ((cc >> q) & 1u) ? 1 : 0;
    f();
    for (int k : S) md[k] = cppde::sliding_mode;
  }

  void root_dg_dp_axpy(int ev, const Vec& x, const T& t, const T& sc, T* out) const {
    const int k = exit_of(ev);
    if (k < 0) {
      eadj.root_dg_dp_axpy(ev, x, t, sc, out);
      return;
    }
    std::vector<int> S;
    std::vector<T> a, Rinv;
    std::vector<std::vector<T>> r;
    if (!prepare(x, t, S, a, r, Rinv)) return;
    const std::size_t n = S.size(), M = std::size_t(1) << n;
    std::size_t kk = 0;
    while (kk < n && S[kk] != k) ++kk;
    if (kk == n) return;
    // g = a_k - 1 or -a_k, da/dp = -R_a^{-1} R_p.
    const T sign = c->roots[ev].exit_on ? T(1.0) : T(-1.0);
    for (std::size_t j = 0; j < n; ++j)
      for (std::size_t cc = 0; cc < M; ++cc) {
        const T scale = T(0.0) - sc * sign * Rinv[kk * n + j] * cppde::detail::combo_weight(a, cc);
        held(S, cc, [&] { eadj.root_gdot_dp_axpy(c->root_of[S[j]], x, t, scale, out); });
      }
  }

  void root_gdot_dx(int ev, const Vec& x, const T& t, Vec& out) const {
    std::vector<int> S;
    std::vector<T> a, Rinv;
    std::vector<std::vector<T>> r;
    if (exit_of(ev) >= 0) {
      zero_armed(out, x.size());
      return;
    }
    if (!prepare(x, t, S, a, r, Rinv)) {
      eadj.root_gdot_dx(ev, x, t, out);
      return;
    }
    const std::size_t n = S.size(), M = std::size_t(1) << n;
    zero_armed(out, x.size());
    // d rate_i / d a_j, with rate_i = sum_c w_c r'_{i,c}.
    std::vector<T> ri(M), dRi(n, T(0.0));
    double terms = 0.0;
    for (std::size_t cc = 0; cc < M; ++cc)
      ri[cc] = cppde::detail::root_rate(c->roots[ev], x, t, c->fm[cc], terms);
    for (std::size_t j = 0; j < n; ++j)
      for (std::size_t cc = 0; cc < M; ++cc)
        dRi[j] = dRi[j] + ri[cc] * cppde::detail::combo_weight(a, cc, static_cast<int>(j));
    // grad = sum_c w_c grad r'_{i,c} - sum_{j,l} dRi_j Rinv_{jl} sum_c w_c grad r'_{S_l,c}.
    std::vector<T> v(n, T(0.0));
    for (std::size_t l = 0; l < n; ++l)
      for (std::size_t j = 0; j < n; ++j) v[l] = v[l] + dRi[j] * Rinv[j * n + l];
    Vec g;
    for (std::size_t cc = 0; cc < M; ++cc) {
      const T w = cppde::detail::combo_weight(a, cc);
      held(S, cc, [&] {
        zero_armed(g, x.size());
        eadj.root_gdot_dx(ev, x, t, g);
        for (std::size_t i = 0; i < x.size(); ++i) out[i] = out[i] + w * g[i];
        for (std::size_t l = 0; l < n; ++l) {
          zero_armed(g, x.size());
          eadj.root_gdot_dx(c->root_of[S[l]], x, t, g);
          for (std::size_t i = 0; i < x.size(); ++i) out[i] = out[i] - v[l] * w * g[i];
        }
      });
    }
  }

  void root_gdot_dp_axpy(int ev, const Vec& x, const T& t, const T& sc, T* out) const {
    std::vector<int> S;
    std::vector<T> a, Rinv;
    std::vector<std::vector<T>> r;
    if (exit_of(ev) >= 0) return;
    if (!prepare(x, t, S, a, r, Rinv)) {
      eadj.root_gdot_dp_axpy(ev, x, t, sc, out);
      return;
    }
    const std::size_t n = S.size(), M = std::size_t(1) << n;
    std::vector<T> ri(M), dRi(n, T(0.0)), v(n, T(0.0));
    double terms = 0.0;
    for (std::size_t cc = 0; cc < M; ++cc)
      ri[cc] = cppde::detail::root_rate(c->roots[ev], x, t, c->fm[cc], terms);
    for (std::size_t j = 0; j < n; ++j)
      for (std::size_t cc = 0; cc < M; ++cc)
        dRi[j] = dRi[j] + ri[cc] * cppde::detail::combo_weight(a, cc, static_cast<int>(j));
    for (std::size_t l = 0; l < n; ++l)
      for (std::size_t j = 0; j < n; ++j) v[l] = v[l] + dRi[j] * Rinv[j * n + l];
    for (std::size_t cc = 0; cc < M; ++cc) {
      const T w = cppde::detail::combo_weight(a, cc);
      held(S, cc, [&] {
        eadj.root_gdot_dp_axpy(ev, x, t, sc * w, out);
        for (std::size_t l = 0; l < n; ++l)
          eadj.root_gdot_dp_axpy(c->root_of[S[l]], x, t, T(0.0) - sc * v[l] * w, out);
      });
    }
  }

  template<class... A> void fixed_dh_dx(A&&... a) const { eadj.fixed_dh_dx(std::forward<A>(a)...); }
  template<class... A> void fixed_dh_dp_axpy(A&&... a) const { eadj.fixed_dh_dp_axpy(std::forward<A>(a)...); }
  template<class... A> void fixed_dtime_dp_axpy(A&&... a) const { eadj.fixed_dtime_dp_axpy(std::forward<A>(a)...); }
  template<class... A> void fixed_dh_dt_axpy(A&&... a) const { eadj.fixed_dh_dt_axpy(std::forward<A>(a)...); }
  template<class... A> void root_dh_dx(A&&... a) const { eadj.root_dh_dx(std::forward<A>(a)...); }
  template<class... A> void root_dh_dp_axpy(A&&... a) const { eadj.root_dh_dp_axpy(std::forward<A>(a)...); }
  template<class... A> void root_dh_dt_axpy(A&&... a) const { eadj.root_dh_dt_axpy(std::forward<A>(a)...); }
};

template<class EvAdj, class Slider>
sliding_event_adjoint<EvAdj, Slider>
make_sliding_event_adjoint(const EvAdj& eadj, const Slider& slider) {
  return {eadj, slider.first.c.get()};
}

template<class Adj, class EvAdj, class Slider>
sliding_adjoint<Adj, EvAdj, Slider>
make_sliding_adjoint(const Adj& adj, const EvAdj& eadj, const Slider& slider) {
  return {adj, eadj, slider.first.c.get()};
}

}  // namespace adjoint
}  // namespace cppde

#endif  // CPPDE_SLIDING_ADJOINT_HPP
