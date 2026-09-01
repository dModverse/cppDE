/*
 Sensitivity transport across an event.

 The reset is sandwiched between a forward and a backward Heun shift, so that
 the parameter dependence of the firing time reaches the sensitivities. See
 vignette("Methods"), section "Sensitivity transport across a discontinuity".

 Copyright (C) 2026 Simon Beyer
 */

#ifndef CPPDE_SALTATION_HPP
#define CPPDE_SALTATION_HPP

#include <vector>
#include <functional>
#include <cmath>
#include <limits>
#include <type_traits>
#include <cppde/cppde_events.hpp>
#include <cppde/cppde_ad_traits.hpp>

namespace cppde {
namespace detail {

// ============================================================================
// Analytical saltation correction for root events (AD path)
//
// dt* comes from the implicit function theorem as a dual quotient, the state
// is carried across the surface on Heun shifts so that the second-order
// components survive, and events triggered together share one roundtrip.
// See vignette("Methods"), section "Root-triggered events".
// ============================================================================

// --- Helper: compute dt* for a root event (IFT + 2nd-order correction) ---
template<class state_type, class time_type, class System>
inline typename state_type::value_type compute_dt_star(
    const state_type& x_before,
    const time_type& t_event,
    System& sys,
    const state_type& f_before,
    const RootEvent<state_type, time_type>& evt)
{
  using value_type = typename state_type::value_type;
  const size_t n = x_before.size();

  // Analytical g_dot = sum_i(dg/dx_i · f_i) + dg/dt
  state_type grad_g(n);
  evt.dg_dx(x_before, t_event, grad_g);
  value_type g_dot = evt.dg_dt(x_before, t_event);
  for (size_t i = 0; i < n; ++i)
    g_dot += grad_g[i] * f_before[i];

  double g_dot_s = scalar_value(g_dot);
  if (std::abs(g_dot_s) < 1e-15) {
    return value_type(0.0);
  }

  // dt* from IFT (dual quotient rule)
  value_type g_val = evt.func(x_before, t_event);
  g_val = g_val - value_type(scalar_value(g_val));
  value_type dt_star = -g_val / g_dot;

  // Second-order IFT correction
  {
    double G_tt;
    if (evt.g_dot_dot) {
      G_tt = evt.g_dot_dot(x_before, t_event);
    } else {
      constexpr double eps = 1e-8;
      state_type x_fwd(n);
      for (size_t i = 0; i < n; ++i)
        x_fwd[i] = x_before[i] + f_before[i] * value_type(eps);

      state_type f_fwd(n);
      sys.first(x_fwd, f_fwd, t_event + time_type(eps));

      state_type grad_g_fwd(n);
      evt.dg_dx(x_fwd, t_event + time_type(eps), grad_g_fwd);
      value_type g_dot_fwd = evt.dg_dt(x_fwd, t_event + time_type(eps));
      for (size_t i = 0; i < n; ++i)
        g_dot_fwd += grad_g_fwd[i] * f_fwd[i];

      G_tt = (scalar_value(g_dot_fwd) - g_dot_s) / eps;
    }

    double corr_coeff = -0.5 * G_tt / g_dot_s;
    dt_star = dt_star + value_type(corr_coeff) * dt_star * dt_star;
  }

  return dt_star;
}

// --- Batch saltation: one Heun roundtrip, N event actions in the middle ---
//
// Events triggered together share one dt* and one surface, so the roundtrip is
// paid once: four right-hand-side evaluations whatever their number. dt* is
// taken from the first of them that has usable gradients.

template<class state_type, class time_type, class System>
inline void saltation_root_analytical_batch(
    state_type& x,
    const state_type& x_before,
    const time_type& t_event,
    System& sys,
    const std::vector<RootEvent<state_type, time_type>>& root_events,
    const std::vector<TriggeredEvent>& triggered)
{
  using value_type = typename state_type::value_type;
  const size_t n = x.size();
  if (n == 0) return;

  // --- 1. RHS before event (full AD) ---
  state_type f_before(n);
  sys.first(x_before, f_before, t_event);

  // --- 2. Compute dt* from the first event with valid gradients ---
  value_type dt_star(0.0);
  bool have_dt_star = false;
  for (const auto& te : triggered) {
    const auto& evt = root_events[te.index];
    if (evt.terminal) continue;
    if (evt.dg_dx && evt.dg_dt) {
      state_type grad_g(n);
      evt.dg_dx(x_before, t_event, grad_g);
      value_type g_dot = evt.dg_dt(x_before, t_event);
      for (size_t i = 0; i < n; ++i)
        g_dot += grad_g[i] * f_before[i];

      if (std::abs(scalar_value(g_dot)) >= 1e-15) {
        dt_star = compute_dt_star(x_before, t_event, sys, f_before, evt);
        have_dt_star = true;
      }
      break;
    }
  }

  if (!have_dt_star) {
    for (const auto& te : triggered) {
      if (!root_events[te.index].terminal) {
        apply_event_action(x, x_before, t_event, root_events[te.index]);
      }
    }
    return;
  }

  // --- 3. Forward Heun shift to event surface ---
  //     f2 is evaluated at t_event + dt_star, not at t_event. The scalar time
  //     is the same either way, but the AD components of the shift are what
  //     carry the curvature term when the right-hand side depends on time.
  state_type x_euler(n);
  for (size_t i = 0; i < n; ++i)
    x_euler[i] = x_before[i] + f_before[i] * dt_star;

  time_type t_star = t_event + dt_star;
  state_type f_euler(n);
  sys.first(x_euler, f_euler, t_star);

  value_type half(0.5);
  state_type x_star(n);
  for (size_t i = 0; i < n; ++i)
    x_star[i] = x_before[i] + half * (f_before[i] + f_euler[i]) * dt_star;

  // --- 4. Apply ALL event actions at the event surface ---
  state_type x_after(n);
  for (size_t i = 0; i < n; ++i) x_after[i] = x_star[i];

  for (const auto& te : triggered) {
    const auto& evt = root_events[te.index];
    if (evt.terminal) continue;
    const int k = evt.state_index;
    if (k >= 0) {
      value_type h = evt.value_func(x_star, t_event);
      switch (evt.method) {
      case EventMethod::Replace:  x_after[k] = h; break;
      case EventMethod::Add:      x_after[k] = x_star[k] + h; break;
      case EventMethod::Multiply: x_after[k] = x_star[k] * h; break;
      }
    }
  }

  // --- 5. Backward Heun shift to grid time ---
  //     x_after lives at t* = t_event + dt_star.  We transport backward
  //     by -dt_star.  f_after at departure (t_star), f_back at arrival (t_event).
  state_type f_after(n);
  sys.first(x_after, f_after, t_star);

  state_type x_back(n);
  for (size_t i = 0; i < n; ++i)
    x_back[i] = x_after[i] - f_after[i] * dt_star;

  state_type f_back(n);
  sys.first(x_back, f_back, t_event);

  for (size_t i = 0; i < n; ++i)
    x[i] = x_after[i] - half * (f_after[i] + f_back[i]) * dt_star;
}

// ============================================================================
// Analytical saltation correction for fixed-time events (AD path)
//
// The residual dt_corr = t_event - scalar(t_event) has scalar part zero and
// carries dt_event/dp in its AD components. The sandwich around it is the one
// the root path uses.
// ============================================================================

// Resets that ride along on an event surface. The default adds none.
struct no_extra_reset {
  template<class state_type, class time_type>
  void operator()(state_type&, const time_type&) const {}
};

template<class state_type, class System, class AtSurface = no_extra_reset>
inline void saltation_fixed_analytical(
    state_type& x,
    const state_type& x_before,
    System& sys,
    const FixedEvent<state_type, typename state_type::value_type>& evt,
    const AtSurface& at_surface = AtSurface())
{
  using value_type = typename state_type::value_type;
  const size_t n = x.size();
  if (n == 0) return;

  value_type dt_corr = evt.time - value_type(scalar_value(evt.time));
  value_type half(0.5);
  value_type t_grid = value_type(scalar_value(evt.time));

  // --- 1. Forward Heun shift to true event time ---
  //     x_before lives at t_grid.  f1 at departure (t_grid),
  //     f2 at arrival (evt.time).
  state_type f1(n);
  sys.first(x_before, f1, t_grid);

  state_type x_euler(n);
  for (size_t i = 0; i < n; ++i)
    x_euler[i] = x_before[i] + f1[i] * dt_corr;

  state_type f2(n);
  sys.first(x_euler, f2, evt.time);

  state_type x_star(n);
  for (size_t i = 0; i < n; ++i)
    x_star[i] = x_before[i] + half * (f1[i] + f2[i]) * dt_corr;

  // --- 2. Apply event action ---
  state_type x_after(n);
  for (size_t i = 0; i < n; ++i) x_after[i] = x_star[i];
  const int k = evt.state_index;
  if (k >= 0) {
    value_type h = evt.value_func(x_star, evt.time);
    switch (evt.method) {
    case EventMethod::Replace:  x_after[k] = h; break;
    case EventMethod::Add:      x_after[k] = x_star[k] + h; break;
    case EventMethod::Multiply: x_after[k] = x_star[k] * h; break;
    }
  }
  // Resets that the jump switches on belong on the same surface, so the whole
  // discontinuity is carried across by one forward and one backward shift.
  at_surface(x_after, evt.time);

  // --- 3. Backward Heun shift to grid time ---
  //     x_after lives at evt.time.  g1 at departure (evt.time),
  //     g2 at arrival (t_grid).
  state_type g1(n);
  sys.first(x_after, g1, evt.time);

  state_type x_back(n);
  for (size_t i = 0; i < n; ++i)
    x_back[i] = x_after[i] - g1[i] * dt_corr;

  state_type g2(n);
  sys.first(x_back, g2, t_grid);

  for (size_t i = 0; i < n; ++i)
    x[i] = x_after[i] - half * (g1[i] + g2[i]) * dt_corr;
}

// ============================================================================
// Apply fixed events at a given time: SFINAE dispatch
//
// double path:  plain event action (no saltation needed)
// AD path:      analytical saltation correction
// ============================================================================

template<class state_type, class Time, class System, class AtSurface = no_extra_reset>
bool apply_fixed_events_at_time(
    state_type& x,
    const Time& t,
    const std::vector<FixedEvent<state_type, typename state_type::value_type>>& evs,
    System& sys,
    const AtSurface& at_surface = AtSurface())
{
  using value_type = typename state_type::value_type;
  const double tt = scalar_value(t);

  // Several events can share a time. The extra resets go on the surface of the
  // last of them, which is the state the jump ends on.
  size_t last = evs.size();
  for (size_t j = 0; j < evs.size(); ++j)
    if (std::abs(scalar_value(evs[j].time) - tt) < 1e-14) last = j;
  if (last == evs.size()) return false;

  for (size_t j = 0; j < evs.size(); ++j) {
    if (std::abs(scalar_value(evs[j].time) - tt) >= 1e-14) continue;
    const auto& e = evs[j];
    if constexpr (std::is_arithmetic_v<value_type>) {
      apply_event_action_fixed(x, x, e);
      if (j == last) at_surface(x, e.time);
    } else {
      state_type x_before = x;
      if (j == last) saltation_fixed_analytical(x, x_before, sys, e, at_surface);
      else           saltation_fixed_analytical(x, x_before, sys, e);
    }
  }
  return true;
}

// ============================================================================
// Merge user times with event times
// ============================================================================

template<class Time, class It, class state_type, class V>
std::vector<Time> merge_user_and_event_times(
    It ubegin, It uend,
    const std::vector<FixedEvent<state_type, V>>& fix)
{
  std::vector<Time> out(ubegin, uend);
  for (const auto& e : fix) out.push_back(e.time);

  std::sort(out.begin(), out.end(),
            [](const auto& a, const auto& b) {
              return scalar_value(a) < scalar_value(b);
            });

  out.erase(std::unique(out.begin(), out.end(),
                        [](const auto& a, const auto& b) {
                          return std::abs(scalar_value(a) - scalar_value(b)) < 1e-14;
                        }),
                        out.end());
  return out;
}

} // namespace detail
} // namespace cppde

#endif // CPPDE_SALTATION_HPP
