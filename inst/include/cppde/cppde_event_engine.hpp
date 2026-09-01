/*
 The integration loops that carry events: one over a controlled stepper, one
 over dense output.

 The loop structure is derived from Boost.Odeint's integrate_times by Karsten
 Ahnert, Mario Mulansky, and Christoph Koke (2011-2015), distributed under the
 Boost Software License, Version 1.0.

 Modified work (root finding, event application, restarts, AD support):
 Copyright (C) 2026 Simon Beyer
 */

#ifndef CPPDE_EVENT_ENGINE_HPP
#define CPPDE_EVENT_ENGINE_HPP

#include <vector>
#include <functional>
#include <cmath>
#include <limits>
#include <algorithm>
#include <type_traits>
#include <memory>
#include <cppde/cppde_step_checker.hpp>
#include <cppde/cppde_stepper_traits.hpp>
#include <cppde/cppde_profiler.hpp>
#include <cppde/cppde_utils.hpp>
#include <cppde/cppde_blas_threads.hpp>
#include <cppde/cppde_ad_traits.hpp>
#include <cppde/cppde_stepper_diagnostics.hpp>
#include <cppde/cppde_events.hpp>
#include <cppde/cppde_saltation.hpp>

namespace cppde {
namespace detail {

// ============================================================================
// Unified stepper reset
//
// In this order: a multistep method discards its history and restarts at order
// one, a dense-output stepper reinitialises both buffers, a controlled stepper
// resets its PI controller, anything else does nothing.
// ============================================================================

template<class S, class State, class Time, class = void>
struct has_reinitialize_at_event : std::false_type {};

template<class S, class State, class Time>
struct has_reinitialize_at_event<S, State, Time,
                                std::void_t<decltype(std::declval<S&>().reinitialize_at_event(
                                    std::declval<State&>(), std::declval<Time>(), std::declval<Time&>()))>
> : std::true_type {};

template<class S, class State, class Time, class = void>
struct has_restart_from_order1 : std::false_type {};

template<class S, class State, class Time>
struct has_restart_from_order1<S, State, Time,
                               std::void_t<decltype(std::declval<S&>().restart_from_order1(
                                   std::declval<State&>(), std::declval<Time>(), std::declval<Time&>()))>
> : std::true_type {};

template<class S, class Time, class = void>
struct has_reset_after_event : std::false_type {};

template<class S, class Time>
struct has_reset_after_event<S, Time,
                            std::void_t<decltype(std::declval<S&>().reset_after_event(
                                std::declval<Time>()))>
> : std::true_type {};

template<class S, class State, class Time>
inline void reset_stepper_unified(S& st, State& x, Time t, Time& dt) {
 if constexpr (::cppde::needs_restart_after_event_v<S>) {
   if constexpr (has_restart_from_order1<S, State, Time>::value) {
     st.restart_from_order1(x, t, dt);
   } else if constexpr (has_reinitialize_at_event<S, State, Time>::value) {
     st.reinitialize_at_event(x, t, dt);
   }
 } else if constexpr (has_reset_after_event<S, Time>::value) {
   st.reset_after_event(dt);
 } else if constexpr (has_reinitialize_at_event<S, State, Time>::value) {
   st.reinitialize_at_event(x, t, dt);
 } else {
   (void)st; (void)x; (void)t; (void)dt;
 }
}

// ============================================================================
// no_dt_estimator: sentinel type: "use dt as-is (no re-estimation)"
// ============================================================================

struct no_dt_estimator {};

// ============================================================================
// EventEngine
// ============================================================================

template<class Stepper, class System, class State, class Time,
        class DtEstimator = no_dt_estimator>
class EventEngine {
public:
 using state_type = State;
 using time_type  = Time;
 using value_type = typename State::value_type;

 using TerminationFunc = std::function<bool(const State&, const Time&)>;

 EventEngine(Stepper& st, System& sys,
             const std::vector<FixedEvent<State, value_type>>& fixed,
             const std::vector<RootEvent<State, Time>>& root,
             DtEstimator dt_est = DtEstimator())
   : m_st(st), m_sys(sys), m_fixed(fixed), m_root(root),
     m_dt_estimator(std::move(dt_est)) {}

 void set_termination(TerminationFunc f) { m_termination = std::move(f); }

 cppde::profiler& get_profiler() const {
   if constexpr (has_controlled_stepper_method<Stepper>::value) {
     return m_st.controlled_stepper().m_prof;
   } else {
     return m_st.m_prof;
   }
 }

private:
 void check_root_triggers(
     std::vector<TriggeredEvent>& triggered,
     const std::vector<double>& last_val,
     const std::vector<double>& curr_val,
     const std::vector<size_t>& fired,
     size_t max_trigger)
 {
   triggered.clear();
   for (size_t i = 0; i < m_root.size(); ++i) {
     if (fired[i] < max_trigger &&
         !std::isnan(last_val[i]) &&
         root_crossed(last_val[i], curr_val[i]) &&
         direction_matches(last_val[i], m_root[i].direction))
       triggered.push_back({i, last_val[i], curr_val[i]});
   }
 }

 // --------------------------------------------------------------------------
 // Apply root events with SFINAE-dispatched saltation correction
 // --------------------------------------------------------------------------

 bool apply_root_events(
     State& x_root, const State& x_before, const Time& t_root,
     const std::vector<TriggeredEvent>& triggered,
     std::vector<size_t>& fired)
 {
   bool has_terminal = false;
   for (const auto& te : triggered)
     if (m_root[te.index].terminal) { has_terminal = true; break; }

   if constexpr (std::is_arithmetic_v<value_type>) {
     for (const auto& te : triggered) {
       size_t i = te.index;
       if (m_root[i].terminal) { fired[i]++; continue; }
       apply_event_action(x_root, x_before, t_root, m_root[i]);
       fired[i]++;
     }
   } else {
     bool has_gradients = false;
     bool has_non_terminal = false;
     for (const auto& te : triggered) {
       if (m_root[te.index].terminal) continue;
       has_non_terminal = true;
       if (m_root[te.index].dg_dx && m_root[te.index].dg_dt) {
         has_gradients = true;
         break;
       }
     }

     if (has_non_terminal && has_gradients) {
       saltation_root_analytical_batch(x_root, x_before, t_root, m_sys,
                                       m_root, triggered);
     } else if (has_non_terminal) {
       for (const auto& te : triggered) {
         if (!m_root[te.index].terminal) {
           apply_event_action(x_root, x_before, t_root, m_root[te.index]);
         }
       }
     }

     for (const auto& te : triggered) {
       fired[te.index]++;
     }
   }
   return has_terminal;
 }

 // An event puts the state on a different trajectory, so the step size is
 // re-estimated from the post-event state the way the solve estimated its
 // first. A multistep method also loses its history at this point.
 void recalibrate_dt(State& x, Time t, Time& dt) {
   if constexpr (!std::is_same_v<DtEstimator, no_dt_estimator>) {
     dt = m_dt_estimator(x, t);
   } else if constexpr (has_controlled_stepper_method<Stepper>::value) {
     dt = Time(odeint_utils::cppde_hin<value_type>(
         m_sys.first, x, t, m_t_final,
         m_st.controlled_stepper().atol(), m_st.controlled_stepper().rtol()));
   } else if constexpr (has_tolerances<Stepper>::value) {
     dt = Time(odeint_utils::cppde_hin<value_type>(
         m_sys.first, x, t, m_t_final, m_st.atol(), m_st.rtol()));
   }
 }

 void init_stepper_after_event(State& x, Time t, Time& dt) {
   recalibrate_dt(x, t, dt);
   m_st.initialize(x, t, dt);
   if constexpr (::cppde::needs_restart_after_event_v<Stepper>) {
     State f0(x.size());
     m_sys.first(x, f0, t);
     m_st.controlled_stepper().stepper().initialize(x, t, f0, dt);
     m_st.controlled_stepper().reset_after_event(dt);
   }
 }

 template<class Checker>
 void reinit_after_event(
     State& x, const Time& t_event, Time& dt,
     Time& t_start, Time& t_end, State& x_at_start,
     std::vector<double>& last_val,
     std::vector<State>& last_state, std::vector<Time>& last_time,
     const std::vector<TriggeredEvent>& triggered,
     size_t& steps, Checker& checker)
 {
   init_stepper_after_event(x, t_event, dt);
   // The restarted step is interpolated at t_start below.
   m_st.set_dense_demand(t_event, true, true);
   m_st.do_step(m_sys);
   ++steps; checker(); checker.reset();
   checker.set_last_order(get_stepper_order(m_st));

   t_start = m_st.previous_time();
   t_end = m_st.current_time();
   dt = m_st.current_time_step();
   checker.set_last_dt(scalar_value(t_end) - scalar_value(t_start));

   m_st.calc_state(t_start, x_at_start);
   for (size_t j = 0; j < m_root.size(); ++j) {
     last_val[j] = scalar_value(m_root[j].func(x_at_start, t_start));
     last_state[j] = x_at_start;
     last_time[j] = t_start;
   }
   // The step restarts on the surface of the events that just fired. Their
   // round-off residual there carries a sign, and keeping it would let the
   // crossing that was just handled be detected a second time.
   for (const auto& te : triggered) last_val[te.index] = 0.0;
 }

 void eval_root_funcs(std::vector<double>& cv, const State& x, const Time& t) {
   for (size_t i = 0; i < m_root.size(); ++i)
     cv[i] = scalar_value(m_root[i].func(x, t));
 }

 // The reset alone, without the transport. Used to read the state a jump ends
 // on, which decides which root conditions it switches.
 bool apply_fixed_events_plain(State& x, const Time& t) {
   const double tt = scalar_value(t);
   bool any = false;
   for (const auto& e : m_fixed)
     if (std::abs(scalar_value(e.time) - tt) < 1e-14) {
       apply_event_action_fixed(x, x, e);
       any = true;
     }
   return any;
 }

 // Fixed events at t, plus the root events their jump switches on. A jump is
 // not a crossing, so every condition is read on both sides of it and fires on
 // a change of sign, one per sweep. Terminal conditions are left out.
 bool apply_fixed_events(State& x, const Time& t,
                         std::vector<size_t>& fired, size_t max_trigger)
 {
   if (m_root.empty())
     return apply_fixed_events_at_time(x, t, m_fixed, m_sys);

   // Called at every output time, so nothing is read or copied until a reset
   // actually sits here.
   const double tt = scalar_value(t);
   bool here = false;
   for (const auto& e : m_fixed)
     if (std::abs(scalar_value(e.time) - tt) < 1e-14) { here = true; break; }
   if (!here) return false;

   std::vector<double> before(m_root.size());
   eval_root_funcs(before, x, t);

   // Which conditions switch is a property of the scalar states, and the
   // transported jump has the same scalar part as the plain one.
   State probe = x;
   apply_fixed_events_plain(probe, t);

   std::vector<size_t> order;
   std::vector<double> after(m_root.size());
   std::vector<size_t> count = fired;
   for (size_t sweep = 0; sweep <= m_root.size(); ++sweep) {
     eval_root_funcs(after, probe, t);
     size_t pick = m_root.size();
     for (size_t i = 0; i < m_root.size(); ++i)
       if (!m_root[i].terminal && count[i] < max_trigger &&
           !std::isnan(before[i]) && root_crossed(before[i], after[i]) &&
           direction_matches(before[i], m_root[i].direction)) { pick = i; break; }
     before = after;
     if (pick == m_root.size()) break;
     order.push_back(pick);
     count[pick]++;
     apply_event_action(probe, probe, t, m_root[pick]);
   }

   auto at_surface = [&](State& xs, const auto& ts) {
     const Time t_surface(ts);
     for (size_t i : order) apply_event_action(xs, xs, t_surface, m_root[i]);
   };
   apply_fixed_events_at_time(x, t, m_fixed, m_sys, at_surface);
   for (size_t i : order) fired[i]++;
   return true;
 }

public:
 // ========================================================================
 // Controlled stepper
 // ========================================================================
 template<class Obs, class Checker>
 size_t process_controlled(
     state_type& x, const std::vector<Time>& times, Time dt,
     Obs& obs_raw, Checker& checker, double root_tol, size_t max_trigger)
 {
   auto& prof = get_profiler();
   auto obs = [&](const auto& x_, const auto& t_) {
     auto _tp = prof.timer(cppde::prof_cat::observer);
     obs_raw(x_, t_);
   };

   size_t steps = 0;
   auto it = times.begin();
   const auto end = times.end();
   Time t = *it;
   m_t_final = scalar_value(times.back());

   std::vector<size_t> fired(m_root.size(), 0);
   if (apply_fixed_events(x, t, fired, max_trigger)) {
     recalibrate_dt(x, t, dt);
     reset_stepper_unified(m_st, x, t, dt);
   }

   obs(x, t); ++it;
   if (it == end) return 0;

   std::vector<double> last_val(m_root.size(), std::numeric_limits<double>::quiet_NaN());
   std::vector<State>  last_state(m_root.size(), x);
   std::vector<Time>   last_time(m_root.size(), t);
   std::vector<TriggeredEvent> triggered;
   triggered.reserve(m_root.size());
   std::vector<double> curr_val(m_root.size());

   while (it != end) {
     Time t_target = *it;
     double t_target_s = scalar_value(t_target);

     while (scalar_value(t) < t_target_s - 1e-14) {
       double rem = t_target_s - scalar_value(t);
       Time dt_step = (scalar_value(dt) < rem) ? dt : Time(rem);
       auto result = m_st.try_step(m_sys, x, t, dt_step);

       if (result == success) {
         ++steps; checker(); checker.reset(); dt = dt_step;
         checker.set_last_order(get_stepper_order(m_st));
         checker.set_last_dt(scalar_value(dt_step));

         if (m_termination && m_termination(x, t)) {
           obs(x, t); return steps;
         }

         eval_root_funcs(curr_val, x, t);
         check_root_triggers(triggered, last_val, curr_val, fired, max_trigger);

         if (!triggered.empty()) {
           localize_root_controlled(
             triggered[0].index,
             last_state[triggered[0].index], last_time[triggered[0].index],
                                                      x, t, triggered[0].last_val, triggered[0].curr_val,
                                                      root_tol, checker);

           State x_before = x;
           Time t_before = t - Time(1e-15);
           obs(x_before, t_before);

           State x_after = x;
           if (apply_root_events(x_after, x, t, triggered, fired)) {
             obs(x_after, t); x = x_after; return steps;
           }
           obs(x_after, t); x = x_after;
           recalibrate_dt(x, t, dt);
           reset_stepper_unified(m_st, x, t, dt);

           for (size_t j = 0; j < m_root.size(); ++j) {
             last_val[j] = std::numeric_limits<double>::quiet_NaN();
             last_state[j] = x; last_time[j] = t;
           }
         } else {
           for (size_t i = 0; i < m_root.size(); ++i) {
             last_val[i] = curr_val[i]; last_state[i] = x; last_time[i] = t;
           }
         }
       } else {
         checker(); dt = dt_step;
       }
     }

     t = t_target;
     if (apply_fixed_events(x, t, fired, max_trigger)) {
       obs(x, t);
       recalibrate_dt(x, t, dt);
       reset_stepper_unified(m_st, x, t, dt);
       for (size_t j = 0; j < m_root.size(); ++j) {
         last_val[j] = std::numeric_limits<double>::quiet_NaN();
         last_state[j] = x; last_time[j] = t; fired[j] = 0;
       }
     } else { obs(x, t); }
     ++it;
   }
   return steps;
 }

 // ========================================================================
 // Dense output stepper
 // ========================================================================
 template<class Obs, class Checker>
 size_t process_dense(
     state_type& x, const std::vector<Time>& times, Time dt,
     Obs& obs_raw, Checker& checker, double root_tol, size_t max_trigger)
 {
   auto& prof = get_profiler();
   auto obs = [&](const auto& x_, const auto& t_) {
     auto _tp = prof.timer(cppde::prof_cat::observer);
     obs_raw(x_, t_);
   };

   size_t steps = 0;
   auto it = times.begin(); auto end = times.end();
   m_t_final = scalar_value(times.back());

   std::vector<size_t> fired(m_root.size(), 0);
   if (apply_fixed_events(x, *it, fired, max_trigger)) {
     recalibrate_dt(x, *it, dt);
     m_st.initialize(x, *it, dt);
   }
   obs(x, *it); ++it;
   if (it == end) return 0;

   const bool fwd = scalar_value(dt) >= 0.0;
   const bool dense_always = !m_root.empty() || static_cast<bool>(m_termination);

   m_st.initialize(x, times.front(), dt);
   m_st.set_dense_demand(*it, dense_always, fwd);
   m_st.do_step(m_sys); ++steps; checker(); checker.reset();
   checker.set_last_order(get_stepper_order(m_st));

   Time t_start = m_st.previous_time();
   Time t_end = m_st.current_time();
   dt = m_st.current_time_step();
   checker.set_last_dt(scalar_value(t_end) - scalar_value(t_start));

   // x_at_start is read by root localisation and the termination
   // predicate only.
   const bool track_bracket = dense_always;

   State x_at_start(x.size());
   if (track_bracket) m_st.calc_state(t_start, x_at_start);

   if (m_termination) {
     State x_check(x.size());
     m_st.calc_state(t_end, x_check);
     if (m_termination(x_check, t_end)) {
       obs(x_check, t_end); x = x_check; return steps;
     }
   }

   std::vector<double> last_val(m_root.size());
   std::vector<State>  last_state(m_root.size(), x_at_start);
   std::vector<Time>   last_time(m_root.size(), t_start);
   for (size_t i = 0; i < m_root.size(); ++i)
     last_val[i] = scalar_value(m_root[i].func(x_at_start, t_start));

   std::vector<TriggeredEvent> triggered;
   triggered.reserve(m_root.size());
   std::vector<double> curr_val(m_root.size());

   while (it != end) {
     while (!less_eq_with_sign(*it, t_end, dt)) {
       if (track_bracket) m_st.calc_state(t_end, x);
       eval_root_funcs(curr_val, x, t_end);
       check_root_triggers(triggered, last_val, curr_val, fired, max_trigger);

       if (!triggered.empty()) {
         State x_root = x_at_start; Time t_root = t_start;
         localize_root_dense(triggered[0].index, x_root, t_root, t_end,
                             triggered[0].last_val, triggered[0].curr_val, root_tol);
         State x_before = x_root;
         Time t_before = t_root - Time(1e-15);
         obs(x_before, t_before);
         if (apply_root_events(x_root, x_before, t_root, triggered, fired)) {
           obs(x_root, t_root); x = x_root; return steps;
         }
         obs(x_root, t_root); x = x_root;
         reinit_after_event(x, t_root, dt, t_start, t_end, x_at_start,
                            last_val, last_state, last_time, triggered,
                            steps, checker);
         continue;
       }

       for (size_t i = 0; i < m_root.size(); ++i) {
         last_val[i] = curr_val[i]; last_state[i] = x; last_time[i] = t_end;
       }

       if (m_termination) {
         m_st.calc_state(t_end, x);
         if (m_termination(x, t_end)) {
           obs(x, t_end); return steps;
         }
       }

       m_st.set_dense_demand(*it, dense_always, fwd);
       m_st.do_step(m_sys); ++steps; checker(); checker.reset();
       checker.set_last_order(get_stepper_order(m_st));
       t_start = m_st.previous_time(); t_end = m_st.current_time();
       dt = m_st.current_time_step();
       checker.set_last_dt(scalar_value(t_end) - scalar_value(t_start));
       if (track_bracket) m_st.calc_state(t_start, x_at_start);
       for (size_t j = 0; j < m_root.size(); ++j) {
         last_val[j] = scalar_value(m_root[j].func(x_at_start, t_start));
         last_state[j] = x_at_start; last_time[j] = t_start;
       }
     }

     while (it != end && less_eq_with_sign(*it, t_end, dt)) {
       Time t_eval = *it;
       if (scalar_value(t_eval) < scalar_value(t_start)) {
         m_st.calc_state(t_start, x); obs(x, t_eval); ++it; continue;
       }
       Time t_eval_s = Time(scalar_value(t_eval));
       m_st.calc_state(t_eval_s, x);
       eval_root_funcs(curr_val, x, t_eval);
       check_root_triggers(triggered, last_val, curr_val, fired, max_trigger);

       if (!triggered.empty()) {
         State x_root = last_state[triggered[0].index];
         Time t_root = last_time[triggered[0].index];
         localize_root_dense(triggered[0].index, x_root, t_root, t_eval,
                             triggered[0].last_val, triggered[0].curr_val, root_tol);
         State x_before = x_root;
         Time t_before = t_root - Time(1e-15);
         if (std::abs(scalar_value(t_eval) - scalar_value(t_before)) >= 1e-14)
           obs(x_before, t_before);
         if (apply_root_events(x_root, x_before, t_root, triggered, fired)) {
           if (std::abs(scalar_value(t_eval) - scalar_value(t_root)) >= 1e-14)
             obs(x_root, t_root);
           x = x_root; return steps;
         }
         if (std::abs(scalar_value(t_eval) - scalar_value(t_root)) >= 1e-14)
           obs(x_root, t_root);
         x = x_root;
         reinit_after_event(x, t_root, dt, t_start, t_end, x_at_start,
                            last_val, last_state, last_time, triggered,
                            steps, checker);
         break;
       }

       bool fef = apply_fixed_events(x, t_eval, fired, max_trigger);
       obs(x, t_eval); ++it;

       if (fef) {
         init_stepper_after_event(x, t_eval_s, dt);
         // The restarted step is interpolated at t_start below.
         m_st.set_dense_demand(t_eval_s, true, fwd);
         m_st.do_step(m_sys); ++steps; checker(); checker.reset();
         checker.set_last_order(get_stepper_order(m_st));
         t_start = m_st.previous_time(); t_end = m_st.current_time();
         dt = m_st.current_time_step();
         checker.set_last_dt(scalar_value(t_end) - scalar_value(t_start));
         State x_ns = x; m_st.calc_state(t_start, x_ns);
         for (size_t j = 0; j < m_root.size(); ++j) {
           last_val[j] = scalar_value(m_root[j].func(x_ns, t_start));
           last_state[j] = x_ns; last_time[j] = t_start; fired[j] = 0;
         }
         break;
       }
       for (size_t i = 0; i < m_root.size(); ++i) {
         last_val[i] = curr_val[i]; last_state[i] = x; last_time[i] = t_eval;
       }
     }
   }
   return steps;
 }

private:
 template<class Checker>
 void localize_root_controlled(
     size_t idx, State& x_lo, Time& t_lo, State& x_hi, Time& t_hi,
     double g_lo, double g_hi, double tol, Checker& checker)
 {
   for (int iter = 0; iter < 50; ++iter) {
     double dti = scalar_value(t_hi) - scalar_value(t_lo);
     if (dti < tol) break;
     double alpha = std::max(0.1, std::min(0.9, -g_lo / (g_hi - g_lo)));
     Time t_mid = t_lo + Time(alpha * dti);
     State x_mid = x_lo; Time t_tmp = t_lo; Time dt_tmp = t_mid - t_lo;
     while (scalar_value(t_tmp) < scalar_value(t_mid) - 1e-15) {
       double r = scalar_value(t_mid) - scalar_value(t_tmp);
       Time sd = (scalar_value(dt_tmp) < r) ? dt_tmp : Time(r);
       if (m_st.try_step(m_sys, x_mid, t_tmp, sd) == success) {
         checker(); checker.reset();
       }
       dt_tmp = sd;
     }
     double g_mid = scalar_value(m_root[idx].func(x_mid, t_mid));
     // An exact zero is the root itself. Moving the lower end onto it and
     // carrying on would step past it and lose the bracket.
     if (g_mid == 0.0) { x_lo = x_mid; t_lo = t_mid; break; }
     if (g_lo * g_mid < 0.0) { x_hi = x_mid; t_hi = t_mid; g_hi = g_mid; }
     else { x_lo = x_mid; t_lo = t_mid; g_lo = g_mid; }
   }
   x_hi = x_lo; t_hi = t_lo;
 }

 void localize_root_dense(
     size_t idx, State& x_root, Time& t_root, Time t_hi,
     double g_lo, double g_hi, double tol)
 {
   Time t_lo = t_root;
   if (g_hi == 0.0) t_lo = t_hi;
   for (int iter = 0; iter < 50; ++iter) {
     double dti = scalar_value(t_hi) - scalar_value(t_lo);
     if (dti < tol) break;
     double alpha = std::max(0.1, std::min(0.9, -g_lo / (g_hi - g_lo)));
     Time t_mid = t_lo + Time(alpha * dti);
     State x_mid(x_root.size()); m_st.calc_state(t_mid, x_mid);
     double g_mid = scalar_value(m_root[idx].func(x_mid, t_mid));
     // An exact zero is the root itself. Moving the lower end onto it and
     // carrying on would step past it and lose the bracket.
     if (g_mid == 0.0) { t_lo = t_mid; break; }
     if (g_lo * g_mid < 0.0) { t_hi = t_mid; g_hi = g_mid; }
     else { t_lo = t_mid; g_lo = g_mid; x_root = x_mid; t_root = t_mid; }
   }
   m_st.calc_state(t_lo, x_root); t_root = t_lo;
 }

 Stepper& m_st; System& m_sys;
 const std::vector<FixedEvent<State, value_type>>& m_fixed;
 const std::vector<RootEvent<State, Time>>& m_root;
 DtEstimator m_dt_estimator;
 TerminationFunc m_termination;
 // Integration endpoint, set at the start of process_{controlled,dense}.
 // Used as the upper-bound hint for cppde_hin when re-estimating the
 // initial step size after an event restart on multistep methods.
 double m_t_final = 0.0;
};

} // namespace detail
} // namespace cppde

#endif // CPPDE_EVENT_ENGINE_HPP
