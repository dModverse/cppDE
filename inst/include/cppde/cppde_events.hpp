/*
 Event specifications and the reset maps applied without transport.

 A fixed event fires at a time the caller supplies, a root event when its
 condition crosses zero. Both carry one reset on one state.

 Copyright (C) 2026 Simon Beyer
 */

#ifndef CPPDE_EVENTS_HPP
#define CPPDE_EVENTS_HPP

#include <vector>
#include <functional>
#include <cmath>
#include <cppde/cppde_ad_traits.hpp>
#include <cppde/cppde_utils.hpp>

namespace cppde {
namespace detail {

// ============================================================================
// Event application methods
// ============================================================================

enum class EventMethod {
  Replace,
  Add,
  Multiply
};

// ============================================================================
// Fixed-time event specification
// ============================================================================

template<class state_type, class value_type>
struct FixedEvent {
  value_type time;
  int state_index;
  std::function<value_type(const state_type&, const value_type&)> value_func;
  EventMethod method;
};

// ============================================================================
// Root-finding event specification
//
// The partials of the root function are provided by the codegen. An AD model
// needs them for the saltation correction, a pure-double model leaves them
// unset.
// ============================================================================

template<class state_type, class time_type>
struct RootEvent {
  using value_t = typename state_type::value_type;

  // Root condition g(x, t) = 0
  std::function<value_t(const state_type&, const time_type&)> func;

  // Affected state index (-1 for none / terminal-only)
  int state_index;

  // Event action value h(x, t)
  std::function<value_t(const state_type&, const time_type&)> value_func;

  EventMethod method;
  bool terminal = false;
  int direction = 0;

  // Partials of g provided by the codegen: dg_dx writes the gradient with
  // respect to the state, dg_dt returns the explicit time derivative. Both
  // carry the model's value_type, so an AD model gets AD partials.
  std::function<void(const state_type&, const time_type&, state_type&)> dg_dx;
  std::function<value_t(const state_type&, const time_type&)> dg_dt;

  // G_tt, the second total time derivative of g along the trajectory. Scalar,
  // because it only scales the second-order correction of dt*. A null callback
  // falls back to a finite difference of g_dot.
  std::function<double(const state_type&, const time_type&)> g_dot_dot;
};

// ============================================================================
// Event tracking
// ============================================================================

struct TriggeredEvent {
  size_t index;
  double last_val;
  double curr_val;
};

// g is sampled at step ends and at output times, so a root can land exactly on
// a sample. A strict sign change misses it, and the stored zero then blocks the
// next interval too. A zero at the lower end was already reported there.
inline bool root_crossed(double g_lo, double g_hi) {
  if (g_lo == 0.0) return false;
  return g_hi == 0.0 || g_lo * g_hi < 0.0;
}

// Called only for an interval root_crossed() accepted, so the sign of g before
// the crossing determines the direction; g at the upper end may be zero.
inline bool direction_matches(double last_val, int direction) {
  if (direction == 0) return true;
  if (direction < 0) return last_val > 0.0;
  return last_val < 0.0;
}

// ============================================================================
// Steady-state termination helper (threshold check, no root-finding)
// ============================================================================

template<class System, class State, class Time>
std::function<bool(const State&, const Time&)> make_steady_state_termination(System& sys, double tol) {
  return [&sys, tol](const State& x, const Time& t) -> bool {
    State dxdt(x.size());
    sys(x, dxdt, t);
    return cppde::max_abs_all_levels_vec(dxdt) < tol;
  };
}

// ============================================================================
// Plain event action (no saltation): works for any value_type
// ============================================================================

template<class state_type, class time_type>
inline void apply_event_action(
    state_type& x,
    const state_type& x_ref,
    const time_type& t,
    const RootEvent<state_type, time_type>& evt)
{
  const int k = evt.state_index;
  if (k >= 0) {
    auto h = evt.value_func(x_ref, t);
    switch (evt.method) {
    case EventMethod::Replace:  x[k] = h; break;
    case EventMethod::Add:      x[k] = x_ref[k] + h; break;
    case EventMethod::Multiply: x[k] = x_ref[k] * h; break;
    }
  }
}

template<class state_type, class value_type>
inline void apply_event_action_fixed(
    state_type& x,
    const state_type& x_ref,
    const FixedEvent<state_type, value_type>& evt)
{
  const int k = evt.state_index;
  if (k >= 0) {
    auto h = evt.value_func(x_ref, evt.time);
    switch (evt.method) {
    case EventMethod::Replace:  x[k] = h; break;
    case EventMethod::Add:      x[k] = x_ref[k] + h; break;
    case EventMethod::Multiply: x[k] = x_ref[k] * h; break;
    }
  }
}

} // namespace detail
} // namespace cppde

#endif // CPPDE_EVENTS_HPP
