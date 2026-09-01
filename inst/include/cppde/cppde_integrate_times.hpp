/*
 Event-aware integration at specified times.

 The public entry points. The event types, the transport across a
 discontinuity and the two integration loops live in the headers included
 below; what is exported from cppde:: is listed at the end of this file.

 The integrate_times signature follows Boost.Odeint by Karsten Ahnert, Mario
 Mulansky, and Christoph Koke (2011-2015), distributed under the Boost
 Software License, Version 1.0.

 Modified work: Copyright (C) 2026 Simon Beyer
 */

#ifndef CPPDE_INTEGRATE_TIMES_HPP
#define CPPDE_INTEGRATE_TIMES_HPP

#include <vector>

#include <cppde/cppde_events.hpp>
#include <cppde/cppde_saltation.hpp>
#include <cppde/cppde_event_engine.hpp>
#include <cppde/cppde_step_checker.hpp>
#include <cppde/cppde_stepper_traits.hpp>

namespace cppde {
namespace detail {

// ============================================================================
// Public API
// ============================================================================

template<class Stepper, class System, class State,
        class TimeIterator, class Time, class Observer,
        class DtEstimator = no_dt_estimator>
size_t integrate_times(
   Stepper& stepper, System system, State& x,
   TimeIterator t_begin, TimeIterator t_end, Time dt, Observer obs,
   const std::vector<FixedEvent<State, typename State::value_type>>& fixed,
   const std::vector<RootEvent<State, Time>>& root,
   StepChecker& checker, double root_tol = 1e-8, size_t max_trigger_root = 1,
   DtEstimator dt_est = DtEstimator(),
   std::function<bool(const State&, const Time&)> termination = nullptr,
   cppde::controlled_stepper_tag = cppde::controlled_stepper_tag())
{
 // Pin BLAS to one thread for this solve, see cppde_blas_threads.hpp. The
 // guard sits here because every stepper reaches one of these two overloads,
 // while the dense LU is skipped by sparse and explicit models.
 cppde::detail::single_thread_blas_scope _cppde_blas_guard;
 auto times = merge_user_and_event_times<Time>(t_begin, t_end, fixed);
 EventEngine<Stepper, System, State, Time, DtEstimator> eng(stepper, system, fixed, root, std::move(dt_est));
 if (termination) eng.set_termination(std::move(termination));
 try {
   size_t steps = eng.process_controlled(x, times, dt, obs, checker, root_tol, max_trigger_root);
   transfer_stepper_diagnostics(stepper, checker);
   return steps;
 } catch (...) {
   transfer_stepper_diagnostics(stepper, checker);
   throw;
 }
}

template<class Stepper, class System, class State,
        class TimeIterator, class Time, class Observer,
        class DtEstimator = no_dt_estimator>
size_t integrate_times_dense(
   Stepper& stepper, System system, State& x,
   TimeIterator t_begin, TimeIterator t_end, Time dt, Observer obs,
   const std::vector<FixedEvent<State, typename State::value_type>>& fixed,
   const std::vector<RootEvent<State, Time>>& root,
   StepChecker& checker, double root_tol = 1e-8, size_t max_trigger_root = 1,
   DtEstimator dt_est = DtEstimator(),
   std::function<bool(const State&, const Time&)> termination = nullptr,
   cppde::dense_output_stepper_tag = cppde::dense_output_stepper_tag())
{
 // Pin BLAS to one thread for this solve, see cppde_blas_threads.hpp. The
 // guard sits here because every stepper reaches one of these two overloads,
 // while the dense LU is skipped by sparse and explicit models.
 cppde::detail::single_thread_blas_scope _cppde_blas_guard;
 auto times = merge_user_and_event_times<Time>(t_begin, t_end, fixed);
 EventEngine<Stepper, System, State, Time, DtEstimator> eng(stepper, system, fixed, root, std::move(dt_est));
 if (termination) eng.set_termination(std::move(termination));
 try {
   size_t steps = eng.process_dense(x, times, dt, obs, checker, root_tol, max_trigger_root);
   transfer_stepper_diagnostics(stepper, checker);
   return steps;
 } catch (...) {
   transfer_stepper_diagnostics(stepper, checker);
   throw;
 }
}

} // namespace detail

using detail::FixedEvent;
using detail::RootEvent;
using detail::EventMethod;
using detail::integrate_times;
using detail::integrate_times_dense;
using detail::make_steady_state_termination;

} // namespace cppde

#endif // CPPDE_INTEGRATE_TIMES_HPP
