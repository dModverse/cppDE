/**
 * @file cppde_jump_adjoint.hpp
 * @brief The adjoint across the jumps between two steps: root and fixed events,
 * switches, restarts after a stall, and the observations that sit on them.
 *
 * Copyright (C) 2026 Simon Beyer
 */

#ifndef CPPDE_JUMP_ADJOINT_HPP
#define CPPDE_JUMP_ADJOINT_HPP

#include <cmath>
#include <cstddef>
#include <vector>

#include <cppde/cppde_ad_traits.hpp>
#include <cppde/cppde_dual_math.hpp>
#include <cppde/cppde_events.hpp>
#include <cppde/cppde_saltation.hpp>
#include <cppde/cppde_switch_modes.hpp>

namespace cppde {
namespace adjoint {

// ---------------------------------------------------------------------------
//  n zeros with their tangent storage bound, so a callee that opens its own
//  dual_arena::scope writes into them in place. Sized to exactly what the
//  callee writes: a later growth copy-constructs the new elements unarmed.
// ---------------------------------------------------------------------------
template<class T>
inline void zero_armed(std::vector<T>& v, std::size_t n) {
  v.assign(n, T(0.0));
  for (std::size_t i = 0; i < n; ++i) cppde::ad_traits::arm_tangents(v[i]);
}

/// What a trajectory needs to take a jump apart: the right-hand side, the
/// model's derivatives of the event expressions, and the events themselves.
/// A model without any passes `no_jumps`, and the boundary compiles out.
struct no_jumps { static constexpr bool active = false; };

template<class System, class EvAdj, class FixedEvents, class RootEvents>
struct jump_terms {
  static constexpr bool active = true;
  System& sys;
  const EvAdj& eadj;
  const FixedEvents& fixed;
  const RootEvents& root;
};

template<class System, class EvAdj, class FixedEvents, class RootEvents>
jump_terms<System, EvAdj, FixedEvents, RootEvents>
make_jump_terms(System& s, const EvAdj& e, const FixedEvents& f,
                const RootEvents& r) { return {s, e, f, r}; }

/// The restart, transposed. initialize() builds the whole Nordsieck history out
/// of one state, zn[0] = x and zn[1] = h f(x, t) with the rest zero, so the
/// history's cotangent collapses onto that state and nothing is left above it.
/// The trajectory start and every event boundary apply the same map.
template<class AdjTerms, class T>
void collapse_restart(const std::vector<T>& w_carry, std::size_t n,
                      const std::vector<T>& x0, double t0, double h0,
                      const AdjTerms& adj, T* w_state, T* w_theta,
                      std::vector<T>& mu, std::vector<T>& jv)
{
  for (std::size_t i = 0; i < n; ++i)
    w_state[i] = (i < w_carry.size()) ? w_carry[i] : T(0.0);
  if (w_carry.size() < 2 * n) return;
  mu.assign(w_carry.begin() + static_cast<std::ptrdiff_t>(n),
            w_carry.begin() + static_cast<std::ptrdiff_t>(2 * n));
  adj.jac_t_vec(x0, mu, t0, jv);
  for (std::size_t i = 0; i < n; ++i) w_state[i] += h0 * jv[i];
  adj.dfdp_t_vec_axpy(x0, mu, t0, h0, w_theta);
}

// ---------------------------------------------------------------------------
//  The adjoint of one jump. In value the Heun sandwich collapses to the
//  saltation dx = R'(dx_before + f_before s) - f_after s, with s the event
//  time's differential. See vignette("Methods"), "The adjoint of an event".
// ---------------------------------------------------------------------------
template<class T = double>
struct jump_workspace {
  std::vector<T> fb, fa, g, wy;
  std::vector<std::vector<T> > path;
  // The reset chain on one sandwich's surface, rebuilt per event.
  std::vector<std::vector<T> > micro;
};

/// One reset, transposed. `w_z` is the cotangent on the state it wrote, `w_y`
/// the one on the state it read, which starts as a copy of `w_z`.
///     Replace   z[k] = h(y, t)
///     Add       z[k] = y[k] + h(y, t)
///     Multiply  z[k] = y[k] * h(y, t)
/// `w_time` takes what the reset contributes to the event time's cotangent.
template<class GradX, class GradP, class GradT, class T>
void reset_transpose(int k, cppde::detail::EventMethod method, const T& h,
                     int idx,
                     const std::vector<T>& y, const T& t, std::size_t n,
                     const T* w_z, T* w_y, T* w_theta,
                     GradX&& dh_dx, GradP&& dh_dp_axpy, std::vector<T>& g,
                     GradT&& dh_dt_axpy, T& w_time)
{
  if (k < 0) return;
  const T wk = w_z[k];
  T c = T(1.0);
  using cppde::detail::EventMethod;
  switch (method) {
    case EventMethod::Replace:  w_y[k] -= wk; break;
    case EventMethod::Add:      break;
    case EventMethod::Multiply: w_y[k] += wk * (h - T(1.0)); c = y[k]; break;
  }
  // A zero cotangent contributes nothing, and asking is not free of the
  // scalar type: under tangents the value alone does not say so.
  zero_armed(g, n);
  dh_dx(idx, y, t, g);
  for (std::size_t i = 0; i < n; ++i) w_y[i] += wk * c * g[i];
  dh_dp_axpy(idx, y, t, wk * c, w_theta);
  // A reset that reads the clock moves with the event time, so dh/dt belongs
  // on the scalar the shift is scaled by.
  dh_dt_axpy(idx, y, t, wk * c, &w_time);
}

/// A batch of root events, which share one surface and one dt*. f reads the
/// modes `mb` on the near side of the surface and `ma` on the far side.
template<class System, class RootEvents, class EvAdj, class AdjTerms, class T>
void apply_root_jump_adjoint(const std::vector<T>& x_before,
                             const std::vector<T>& x_after,
                             double t, const RootEvents& root_events,
                             const std::vector<cppde::detail::TriggeredEvent>& triggered,
                             System& sys, const EvAdj& eadj,
                             const AdjTerms& adj, std::size_t n,
                             const T* w_out, T* w_in, T* w_theta,
                             jump_workspace<T>& ws,
                             const std::vector<signed char>& mb = {},
                             const std::vector<signed char>& ma = {})
{
  using cppde::detail::use_switch_modes;
  use_switch_modes(mb);
  zero_armed(ws.fb, n);
  sys.first(x_before, ws.fb, t);

  // Which event dt* came from, picked the way the forward run picks it: the
  // first triggered, non-terminal one whose gradients are there.
  std::size_t src = triggered.size();
  T g_dot = T(0.0);
  for (std::size_t j = 0; j < triggered.size(); ++j) {
    const auto& evt = root_events[triggered[j].index];
    if (evt.terminal) continue;
    if (evt.dg_dx && evt.dg_dt) {
      zero_armed(ws.g, n);
      evt.dg_dx(x_before, t, ws.g);
      T gd = evt.dg_dt(x_before, t);
      for (std::size_t i = 0; i < n; ++i) gd += ws.g[i] * ws.fb[i];
      if (std::abs(ad_traits::scalar_value(gd)) >= 1e-15) { src = j; g_dot = gd; }
      break;
    }
  }
  const bool shifted = (src < triggered.size());

  const T t_ad = T(t);

  // The shift, s = -(g - value(g)) / g_dot: zero in value, dt*/dv in the
  // tangent. It reads off the stored state because the localisation is
  // tangent-blind, so that state holds dx/dv at a fixed time and g contracts
  // it to the event time's. Zero without gradients, which collapses everything
  // below to the bare reset transpose.
  T s = T(0.0);
  if (shifted) {
    const T gv = root_events[triggered[src].index].func(x_before, t_ad);
    s = (T(ad_traits::scalar_value(gv)) - gv) / g_dot;
  }
  const T t_s = t_ad + s;

  // Forward: x_e = x_b + f_b s, x_* = x_b + (f_b + f_e) s / 2, x_a = R(x_*),
  // x_k = x_a - f_a s = the stored state, x_out = x_a - (f_a + f_k) s / 2.
  // Transposed below; see vignette("Methods"), "Transposing the root sandwich".
  use_switch_modes(ma);
  zero_armed(ws.fa, n);
  sys.first(x_after, ws.fa, t);
  use_switch_modes(mb);
  std::vector<T> xe, xst, xa, xk, fe, fa, fk;
  std::vector<T> wfa, wfk, wfb, wfe, wxb, wv, wa;
  for (auto* v : {&xe, &xst, &xa, &xk, &fe, &fa, &fk,
                  &wfa, &wfk, &wfb, &wfe, &wxb, &wv, &wa})
    zero_armed(*v, n);

  // The forward sandwich, rebuilt as saltation_root_analytical_batch builds it.
  // The transpose below reverses that straight-line program assignment by
  // assignment: pruning on the value of s would drop an s^2 factor, which
  // leaves the gradient and stays in the Hessian.
  for (std::size_t i = 0; i < n; ++i) xe[i] = x_before[i] + ws.fb[i] * s;
  sys.first(xe, fe, t_s);
  for (std::size_t i = 0; i < n; ++i)
    xst[i] = x_before[i] + T(0.5) * (ws.fb[i] + fe[i]) * s;
  for (std::size_t i = 0; i < n; ++i) xa[i] = xst[i];
  for (std::size_t j = 0; j < triggered.size(); ++j) {
    const auto& evt = root_events[triggered[j].index];
    if (!evt.terminal)
      cppde::detail::apply_event_action(xa, xa, t_s, evt);
  }
  use_switch_modes(ma);
  sys.first(xa, fa, t_s);
  for (std::size_t i = 0; i < n; ++i) xk[i] = xa[i] - fa[i] * s;
  sys.first(xk, fk, t_ad);

  //  x_out = x_a - (f_a + f_k) s / 2
  T w_s = T(0.0);
  for (std::size_t i = 0; i < n; ++i) {
    wa[i]  = w_out[i];
    wfa[i] = (T(0.0) - T(0.5)) * s * w_out[i];
    wfk[i] = wfa[i];
    w_s   -= T(0.5) * (fa[i] + fk[i]) * w_out[i];
  }

  //  f_k = f(x_k, t)
  adj.jac_t_vec(xk, wfk, t, wv);
  adj.dfdp_t_vec_axpy(xk, wfk, t, T(1.0), w_theta);

  //  x_k = x_a - f_a s
  for (std::size_t i = 0; i < n; ++i) {
    wa[i]  += wv[i];
    wfa[i] -= s * wv[i];
    w_s    -= fa[i] * wv[i];
  }

  //  f_a = f(x_a, t_s), so t_s takes df/dt
  adj.jac_t_vec(xa, wfa, t, wv);
  for (std::size_t i = 0; i < n; ++i) wa[i] += wv[i];
  adj.dfdp_t_vec_axpy(xa, wfa, t, T(1.0), w_theta);
  w_s += adj.dfdt_dot(xa, wfa, t_s);

  //  x_a = R(x_*), every reset reading the surface state the forward reset
  for (std::size_t i = 0; i < n; ++i) w_in[i] = wa[i];
  for (std::size_t j = triggered.size(); j-- > 0;) {
    const std::size_t idx = triggered[j].index;
    const auto& evt = root_events[idx];
    if (evt.terminal) continue;
    const T h = (evt.state_index >= 0 && evt.value_func)
                   ? evt.value_func(xst, t_s) : T(0.0);
    for (std::size_t i = 0; i < n; ++i) wa[i] = w_in[i];
    reset_transpose(
        evt.state_index, evt.method, h, static_cast<int>(idx), xst, t_s, n,
        wa.data(), w_in, w_theta,
        [&](int e, const std::vector<T>& y, const T& tt, std::vector<T>& o)
          { eadj.root_dh_dx(e, y, tt, o); },
        [&](int e, const std::vector<T>& y, const T& tt, const T& sc, T* o)
          { eadj.root_dh_dp_axpy(e, y, tt, sc, o); },
        ws.g,
        [&](int e, const std::vector<T>& y, const T& tt, const T& sc, T* o)
          { eadj.root_dh_dt_axpy(e, y, tt, sc, o); },
        w_s);
  }

  //  x_* = x_b + (f_b + f_e) s / 2
  use_switch_modes(mb);
  for (std::size_t i = 0; i < n; ++i) {
    wxb[i] = w_in[i];
    wfb[i] = T(0.5) * s * w_in[i];
    wfe[i] = wfb[i];
    w_s   += T(0.5) * (ws.fb[i] + fe[i]) * w_in[i];
  }

  //  f_e = f(x_e, t_s), the same pair the other way
  adj.jac_t_vec(xe, wfe, t, wv);
  adj.dfdp_t_vec_axpy(xe, wfe, t, T(1.0), w_theta);
  w_s += adj.dfdt_dot(xe, wfe, t_s);

  //  x_e = x_b + f_b s
  for (std::size_t i = 0; i < n; ++i) {
    wxb[i] += wv[i];
    wfb[i] += s * wv[i];
    w_s    += ws.fb[i] * wv[i];
  }
  for (std::size_t i = 0; i < n; ++i) w_in[i] = wxb[i];

  if (!shifted) return;

  // Where the shift comes from: s = -(grad g . dx + dg/dp . dp) / g_dot, and
  // the forward run corrects it, so ds/ds_lin is 1 + 2 c2 s, one in value and
  // a tangent in the second order.
  const std::size_t idx = triggered[src].index;
  const double c2 = cppde::detail::root_ift_curvature(
      x_before, t_ad, sys, ws.fb, root_events[idx],
      ad_traits::scalar_value(g_dot));
  const T c = (T(0.0) - w_s / g_dot) * (T(1.0) + T(2.0 * c2) * s);
  zero_armed(ws.g, n);
  root_events[idx].dg_dx(x_before, t, ws.g);
  for (std::size_t i = 0; i < n; ++i) w_in[i] += c * ws.g[i];
  eadj.root_dg_dp_axpy(static_cast<int>(idx), x_before, t, c, w_theta);

  // g_dot reads the state too, so the quotient has a second half:
  //
  //   ds/dx = -grad g / g_dot - (s / g_dot) grad g_dot,
  //
  // which contains g, zero in value: absent from the gradient and full
  // strength in the Hessian. The model differentiates g_dot whole, because
  // J' grad g alone holds only for a g linear in x and blind to the clock.
  const T cs = c * s;
  zero_armed(wv, n);
  eadj.root_gdot_dx(static_cast<int>(idx), x_before, t, wv);
  for (std::size_t i = 0; i < n; ++i) w_in[i] += cs * wv[i];
  eadj.root_gdot_dp_axpy(static_cast<int>(idx), x_before, t, cs, w_theta);

  //  f_b = f(x_b, t). grad g_dot already includes the route through f_b that
  //  the shift takes, so only the two Heun legs feed this one.
  adj.jac_t_vec(x_before, wfb, t, wv);
  for (std::size_t i = 0; i < n; ++i) w_in[i] += wv[i];
  adj.dfdp_t_vec_axpy(x_before, wfb, t, T(1.0), w_theta);
}

/// A restart after a stall, transposed. No event fired; the state followed
/// f(t_before, x_before) from t_before to t, past a jump of f in t.
template<class AdjTerms, class T>
void apply_crossing_adjoint(const std::vector<T>& x_before, double t_before,
                            double t, const AdjTerms& adj, std::size_t n,
                            const T* w_after, T* w_before, T* w_theta,
                            jump_workspace<T>& ws)
{
  const double d = t - t_before;
  zero_armed(ws.wy, n);
  for (std::size_t i = 0; i < n; ++i) ws.wy[i] = w_after[i];
  zero_armed(ws.g, n);
  adj.jac_t_vec(x_before, ws.wy, t_before, ws.g);
  for (std::size_t i = 0; i < n; ++i) w_before[i] = w_after[i] + d * ws.g[i];
  adj.dfdp_t_vec_axpy(x_before, ws.wy, t_before, d, w_theta);
}

/// The fixed events at one time, each its own sandwich, applied in order. The
/// last of them applies the root resets a jump switched on.
template<class System, class FixedEvents, class RootEvents, class EvAdj,
         class AdjTerms, class T>
void apply_fixed_jump_adjoint(const std::vector<T>& x_before,
                              double t, const FixedEvents& fixed_events,
                              const RootEvents& root_events,
                              const std::vector<std::size_t>& switched,
                              System& sys, const EvAdj& eadj,
                              const AdjTerms& adj, std::size_t n,
                              const T* w_out, T* w_in, T* w_theta,
                              jump_workspace<T>& ws,
                              const std::vector<signed char>& mb = {},
                              const std::vector<signed char>& ma = {})
{
  using cppde::detail::use_switch_modes;
  // Which of them fire here, and which is last: the same test the engine makes.
  std::vector<int> fired;
  for (std::size_t j = 0; j < fixed_events.size(); ++j)
    if (std::abs(ad_traits::scalar_value(fixed_events[j].time) - t) < 1e-14)
      fired.push_back(static_cast<int>(j));

  // The path the forward took, which the store does not keep: one whole
  // sandwich per event, the switched resets on the surface of the last. A bare
  // reset agrees in value and misses dx_after/dt* = f_before - f_after.
  const std::size_t nf = fired.size();
  std::vector<T> tmp;
  zero_armed(tmp, n);
  // The modes change on the surface of the last of them, as in the forward run.
  auto at_surface = [&](std::vector<T>& xs, const T& te) {
    for (std::size_t q = 0; q < switched.size(); ++q) {
      tmp = xs;
      cppde::detail::apply_event_action(xs, tmp, te, root_events[switched[q]]);
    }
    use_switch_modes(ma);
  };
  ws.path.assign(nf + 1, std::vector<T>());
  ws.path[0] = x_before;
  for (std::size_t j = 0; j < nf; ++j) {
    use_switch_modes(mb);
    zero_armed(ws.path[j + 1], n);
    ws.path[j + 1] = ws.path[j];
    const auto& evt = fixed_events[fired[j]];
    if constexpr (std::is_arithmetic<T>::value) {
      cppde::detail::apply_event_action_fixed(ws.path[j + 1], ws.path[j], evt);
      if (j + 1 == nf) at_surface(ws.path[j + 1], evt.time);
    } else if (j + 1 == nf) {
      cppde::detail::saltation_fixed_analytical(ws.path[j + 1], ws.path[j], sys,
                                                evt, at_surface);
    } else {
      cppde::detail::saltation_fixed_analytical(ws.path[j + 1], ws.path[j], sys,
                                                evt);
    }
  }

  // Backwards through the same chain: the root path's sandwich with the event
  // time's residual tau in place of the shift. Products of two factors that
  // vanish in value are dropped.
  std::vector<T> w(w_out, w_out + n);
  std::vector<T> f1, f2, g1, g2, xe, xk, jv;
  zero_armed(f1, n); zero_armed(f2, n);
  zero_armed(g1, n); zero_armed(g2, n);
  zero_armed(xe, n); zero_armed(xk, n); zero_armed(jv, n);
  ws.wy.assign(n, T(0.0));
  zero_armed(ws.wy, n);

  for (std::size_t j = nf; j-- > 0;) {
    const auto& evt = fixed_events[fired[j]];
    const std::vector<T>& y = ws.path[j];
    const T tau = evt.time - T(ad_traits::scalar_value(evt.time));

    // The half that reaches the surface, replayed so the resets are
    // transposed where the forward applied them.
    use_switch_modes(mb);
    sys.first(y, f1, t);
    for (std::size_t i = 0; i < n; ++i) xe[i] = y[i] + f1[i] * tau;
    sys.first(xe, f2, evt.time);
    const std::size_t n_res = 1 + ((j + 1 == nf) ? switched.size() : 0);
    ws.micro.assign(n_res + 1, std::vector<T>());
    zero_armed(ws.micro[0], n);
    for (std::size_t i = 0; i < n; ++i)
      ws.micro[0][i] = y[i] + T(0.5) * (f1[i] + f2[i]) * tau;
    zero_armed(ws.micro[1], n);
    ws.micro[1] = ws.micro[0];
    cppde::detail::apply_event_action_fixed(ws.micro[1], ws.micro[0], evt);
    for (std::size_t q = 1; q < n_res; ++q) {
      zero_armed(ws.micro[q + 1], n);
      ws.micro[q + 1] = ws.micro[q];
      cppde::detail::apply_event_action(ws.micro[q + 1], ws.micro[q], evt.time,
                                        root_events[switched[q - 1]]);
    }
    const std::vector<T>& xa = ws.micro[n_res];

    // Leaving the surface: x_out = x_a - (g1 + g2) tau / 2, g1 at the event
    // time and g2 at the grid time.
    if (j + 1 == nf) use_switch_modes(ma);
    sys.first(xa, g1, evt.time);
    for (std::size_t i = 0; i < n; ++i) xk[i] = xa[i] - g1[i] * tau;
    sys.first(xk, g2, t);
    for (std::size_t i = 0; i < n; ++i) ws.wy[i] = w[i];
    adj.jac_t_vec(xa, ws.wy, t, jv);
    adj.dfdp_t_vec_axpy(xa, ws.wy, t, T(0.0) - tau, w_theta);
    T w_s = T(0.0);
    for (std::size_t i = 0; i < n; ++i)
      w_s += T(0.5) * (tau * (g1[i] * jv[i]) - (g1[i] + g2[i]) * w[i]);
    // df/dt of g1 at the event time, cotangent -tau w / 2.
    for (std::size_t i = 0; i < n; ++i) ws.wy[i] = (T(0.0) - T(0.5) * tau) * w[i];
    w_s += adj.dfdt_dot(xa, ws.wy, evt.time);
    for (std::size_t i = 0; i < n; ++i) w[i] -= tau * jv[i];

    // The resets on that surface, switched ones first because they came last.
    for (std::size_t q = n_res; q-- > 1;) {
      const auto& re = root_events[switched[q - 1]];
      const std::vector<T>& yy = ws.micro[q];
      const T h = (re.state_index >= 0 && re.value_func)
                     ? re.value_func(yy, evt.time) : T(0.0);
      ws.wy = w;
      reset_transpose(
          re.state_index, re.method, h, static_cast<int>(switched[q - 1]), yy,
          evt.time, n, w.data(), ws.wy.data(), w_theta,
          [&](int e, const std::vector<T>& yv, const T& tt, std::vector<T>& o)
            { eadj.root_dh_dx(e, yv, tt, o); },
          [&](int e, const std::vector<T>& yv, const T& tt, const T& sc, T* o)
            { eadj.root_dh_dp_axpy(e, yv, tt, sc, o); },
          ws.g,
          [&](int e, const std::vector<T>& yv, const T& tt, const T& sc, T* o)
            { eadj.root_dh_dt_axpy(e, yv, tt, sc, o); },
          w_s);
      w.swap(ws.wy);
    }
    {
      const std::vector<T>& yy = ws.micro[0];
      const T h = (evt.state_index >= 0 && evt.value_func)
                     ? evt.value_func(yy, evt.time) : T(0.0);
      ws.wy = w;
      reset_transpose(
          evt.state_index, evt.method, h, fired[j], yy, evt.time, n,
          w.data(), ws.wy.data(), w_theta,
          [&](int e, const std::vector<T>& yv, const T& tt, std::vector<T>& o)
            { eadj.fixed_dh_dx(e, yv, tt, o); },
          [&](int e, const std::vector<T>& yv, const T& tt, const T& sc, T* o)
            { eadj.fixed_dh_dp_axpy(e, yv, tt, sc, o); },
          ws.g,
          [&](int e, const std::vector<T>& yv, const T& tt, const T& sc, T* o)
            { eadj.fixed_dh_dt_axpy(e, yv, tt, sc, o); },
          w_s);
      w.swap(ws.wy);
    }

    // Entering it: x_* = x_b + (f1 + f2) tau / 2, the same pair the other way.
    use_switch_modes(mb);
    for (std::size_t i = 0; i < n; ++i) ws.wy[i] = w[i];
    adj.jac_t_vec(y, ws.wy, t, jv);
    adj.dfdp_t_vec_axpy(xe, ws.wy, t, tau, w_theta);
    for (std::size_t i = 0; i < n; ++i)
      w_s += T(0.5) * ((f1[i] + f2[i]) * w[i] + tau * (f1[i] * jv[i]));
    // df/dt of f2 at the event time, cotangent tau w / 2.
    for (std::size_t i = 0; i < n; ++i) ws.wy[i] = T(0.5) * tau * w[i];
    w_s += adj.dfdt_dot(xe, ws.wy, evt.time);
    for (std::size_t i = 0; i < n; ++i) w[i] += tau * jv[i];

    // kappa scales the event time's derivative, which reads the parameters.
    eadj.fixed_dtime_dp_axpy(fired[j], w_s, w_theta);
  }

  for (std::size_t i = 0; i < n; ++i) w_in[i] = w[i];
}

/// One seed inside a step: the time and the row of seeds.
template<class T>
struct obs_ref { double t; const T* w; };

/// The most substeps a checked step is swept in, the most pieces for the
/// extrapolated flow.
constexpr int default_max_sub = 1024;
/// From this many substeps on, a halving that does not shrink the difference
/// ends the refinement.
constexpr int stall_substeps = 32;

/// The cotangent on the state entering a jump that was read off the dense
/// output of the step below it, waiting for that step's sweep.
template<class T>
struct jump_handover {
  std::vector<T> w;
  double t = 0.0;
  bool armed = false;

  /// Hands `w_before` on to the step below, or onto the initial state `w_x`
  /// when the jump is at the run's own start.
  template<class Event>
  void take(const Event& e, const std::vector<T>& w_before, std::vector<T>& w_x) {
    if (e.after_step > 0) {
      w = w_before; t = e.t_before; armed = true;
    } else {
      for (std::size_t i = 0; i < w_x.size(); ++i) w_x[i] += w_before[i];
    }
  }
};

/// The seeds step `k` holds into `out`: the handover a jump left, then the
/// observations inside it. One a jump produced is a value and is seeded at the
/// boundary instead. Returns whether anything observes inside the step.
template<class Store, class T>
bool step_seeds(const Store& store, std::size_t k, const T* seeds, std::size_t n,
                std::size_t& next_obs, jump_handover<T>& h, std::vector<obs_ref<T>>& out)
{
  std::size_t obs_lo = next_obs;
  while (obs_lo > 0 && store.obs(obs_lo - 1).step == k + 1) --obs_lo;
  const bool observed = obs_lo < next_obs || h.armed;
  out.clear();
  if (h.armed) { out.push_back({h.t, h.w.data()}); h.armed = false; }
  for (std::size_t o = obs_lo; o < next_obs; ++o) {
    if (store.obs(o).event < store.n_events()) continue;
    out.push_back({store.obs(o).t, seeds + o * n});
  }
  next_obs = obs_lo;
  return observed;
}

/// The observations before the first step, which are the initial state itself,
/// except the value a jump at the start produced, which that jump seeded.
template<class Store, class T>
void seed_initial_state(const Store& store, std::size_t next_obs, const T* seeds,
                        std::size_t n, std::vector<T>& w_x)
{
  while (next_obs > 0) {
    --next_obs;
    if (store.obs(next_obs).event < store.n_events()) continue;
    const T* w = seeds + next_obs * n;
    for (std::size_t i = 0; i < n; ++i) w_x[i] += w[i];
  }
}

/// Switches triggered together, and nothing else, each took their own sandwich
/// in the forward run, in the order triggered; their transposes run backwards.
/// A switch moves no state, so every sandwich starts from the stored state.
template<class System, class RootEvents, class EvAdj, class AdjTerms, class T>
void apply_root_jumps_adjoint(const std::vector<T>& x_before,
                              const std::vector<T>& x_after,
                              double t, const RootEvents& root_events,
                              const std::vector<cppde::detail::TriggeredEvent>& triggered,
                              System& sys, const EvAdj& eadj,
                              const AdjTerms& adj, std::size_t n,
                              const T* w_out, T* w_in, T* w_theta,
                              jump_workspace<T>& ws,
                              const std::vector<signed char>& mb,
                              const std::vector<signed char>& ma)
{
  std::vector<std::size_t> sw;
  bool other = false;
  for (std::size_t j = 0; j < triggered.size(); ++j) {
    const auto& evt = root_events[triggered[j].index];
    if (evt.terminal) continue;
    if (evt.mode >= 0 && evt.dg_dx && evt.dg_dt) sw.push_back(j);
    else other = true;
  }
  if (other || sw.size() < 2 || mb.empty()) {
    apply_root_jump_adjoint(x_before, x_after, t, root_events, triggered, sys, eadj,
                            adj, n, w_out, w_in, w_theta, ws, mb, ma);
    return;
  }
  std::vector<std::vector<signed char>> modes(sw.size() + 1, mb);
  for (std::size_t k = 0; k < sw.size(); ++k) {
    modes[k + 1] = modes[k];
    const int m = root_events[triggered[sw[k]].index].mode;
    modes[k + 1][m] = ma[m];
  }
  std::vector<T> wa(w_out, w_out + n), wb(n);
  for (std::size_t k = sw.size(); k-- > 0;) {
    const std::vector<cppde::detail::TriggeredEvent> one{triggered[sw[k]]};
    apply_root_jump_adjoint(x_before, x_after, t, root_events, one, sys, eadj, adj, n,
                            wa.data(), wb.data(), w_theta, ws, modes[k], modes[k + 1]);
    wa.swap(wb);
  }
  for (std::size_t i = 0; i < n; ++i) w_in[i] = wa[i];
}

/// The adjoint of the intervention `ei`: `w_after`, the cotangent on the state
/// it produced, gets the observations of that value and goes through the jump
/// into `w_before` on the state entering it, the share of theta into `wp`.
template<class Store, class AdjTerms, class Jumps, class T>
void event_adjoint(const Store& store, std::size_t ei, const T* seeds, std::size_t n,
                   const AdjTerms& adj, const Jumps& jumps, std::vector<T>& w_after,
                   std::vector<T>& w_before, T* wp, jump_workspace<T>& ws)
{
  const auto& e = store.event(ei);
  for (std::size_t o = 0; o < store.n_obs(); ++o)
    if (store.obs(o).event == ei)
      for (std::size_t i = 0; i < n; ++i) w_after[i] += seeds[o * n + i];
  w_before.assign(n, T(0.0));
  if (e.crossing)
    apply_crossing_adjoint(e.x_before, e.t_before, e.t, adj, n,
                           w_after.data(), w_before.data(), wp, ws);
  else if (e.root)
    apply_root_jumps_adjoint(e.x_before, e.x_after, e.t, jumps.root,
                            e.triggered, jumps.sys, jumps.eadj, adj, n,
                            w_after.data(), w_before.data(), wp, ws,
                            e.modes_before, e.modes_after);
  else
    apply_fixed_jump_adjoint(e.x_before, e.t, jumps.fixed, jumps.root,
                             e.switched, jumps.sys, jumps.eadj, adj, n,
                             w_after.data(), w_before.data(), wp, ws,
                             e.modes_before, e.modes_after);
}

// ---------------------------------------------------------------------------
//  A run that took no step, which is a run over a single time. Every
//  observation there is the initial state itself, or the value a jump at that
//  time produced, which reaches the initial state and the parameters through
//  that jump's adjoint. The flow in between has length zero and is the
//  identity.
// ---------------------------------------------------------------------------
template<class Store, class AdjTerms, class Jumps, class T>
void sweep_without_steps(const Store& store, std::size_t n, const T* seeds,
                         const AdjTerms& adj, const Jumps& jumps,
                         std::vector<T>& w_x, std::vector<T>& w_theta,
                         jump_workspace<T>& ws)
{
  seed_initial_state(store, store.n_obs(), seeds, n, w_x);
  if constexpr (Jumps::active) {
    std::vector<T> w_after, w_before;
    for (std::size_t ei = 0; ei < store.n_events(); ++ei) {
      w_after.assign(n, T(0.0));
      event_adjoint(store, ei, seeds, n, adj, jumps, w_after, w_before,
                    w_theta.data(), ws);
      for (std::size_t i = 0; i < n; ++i) w_x[i] += w_before[i];
    }
  } else {
    (void)adj; (void)jumps; (void)w_theta; (void)ws;
  }
}

}  // namespace adjoint
}  // namespace cppde

#endif  // CPPDE_JUMP_ADJOINT_HPP
