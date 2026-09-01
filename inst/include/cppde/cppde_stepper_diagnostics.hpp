/*
 Stepper introspection: the current method order and the counters a solve
 reports back. Every stepper exposes the same names, so the traits below pick
 out what a given one actually provides.

 Copyright (C) 2026 Simon Beyer
 */

#ifndef CPPDE_STEPPER_DIAGNOSTICS_HPP
#define CPPDE_STEPPER_DIAGNOSTICS_HPP

#include <type_traits>

namespace cppde {
namespace detail {

// ============================================================================
// Diagnostics helpers
//
// Steppers, controllers and dense-output wrappers expose the same names:
// n_accepted(), n_rejected(), n_fevals(), n_jevals(), current_method_order().
// The traits below pick out the ones a given stepper actually has.
// ============================================================================

// --- Extract current method order (used per accepted step) ---

template<class T, class = void>
struct has_method_order : std::false_type {};

template<class T>
struct has_method_order<T, std::void_t<
  decltype(std::declval<const T&>().current_method_order())>>
  : std::true_type {};

template<class Stepper>
inline auto get_stepper_order(const Stepper& st)
  -> std::enable_if_t<has_method_order<Stepper>::value, int>
  {
    return st.current_method_order();
  }

// Fallback for unknown steppers
template<class Stepper>
inline auto get_stepper_order(const Stepper&)
  -> std::enable_if_t<!has_method_order<Stepper>::value, int>
  {
    return 0;
  }

// --- Transfer exact counters from stepper into StepChecker (called once) ---

template<class T, class = void>
struct has_diagnostics_counters : std::false_type {};

template<class T>
struct has_diagnostics_counters<T, std::void_t<
  decltype(std::declval<const T&>().n_accepted()),
  decltype(std::declval<const T&>().n_rejected()),
  decltype(std::declval<const T&>().n_fevals()),
  decltype(std::declval<const T&>().n_jevals())>>
    : std::true_type {};

// --- Profiler report dispatch (SFINAE) ---
// Must be declared before transfer_stepper_diagnostics which calls it.
template<class T, class = void>
struct has_controlled_stepper_method : std::false_type {};

template<class T>
struct has_controlled_stepper_method<T, std::void_t<
 decltype(std::declval<const T&>().controlled_stepper())>>
 : std::true_type {};

template<class T, class = void>
struct has_tolerances : std::false_type {};

template<class T>
struct has_tolerances<T, std::void_t<
 decltype(std::declval<const T&>().atol()),
 decltype(std::declval<const T&>().rtol())>>
 : std::true_type {};

template<class T, class = void>
struct has_report_profiler_method : std::false_type {};

template<class T>
struct has_report_profiler_method<T, std::void_t<
 decltype(std::declval<const T&>().report_profiler())>>
 : std::true_type {};

template<class Stepper>
inline void report_profiler_if_available(const Stepper& st) {
 if constexpr (has_controlled_stepper_method<Stepper>::value) {
   auto& ctrl = st.controlled_stepper();
   if constexpr (has_report_profiler_method<std::decay_t<decltype(ctrl)>>::value) {
     ctrl.report_profiler();
   }
 } else if constexpr (has_report_profiler_method<Stepper>::value) {
   st.report_profiler();
 }
}

template<class T, class = void>
struct has_n_setups_method : std::false_type {};

template<class T>
struct has_n_setups_method<T, std::void_t<
  decltype(std::declval<const T&>().n_setups())>>
  : std::true_type {};

template<class Stepper, class Checker>
inline auto transfer_stepper_diagnostics(const Stepper& st, Checker& checker)
  -> std::enable_if_t<has_diagnostics_counters<Stepper>::value>
  {
    checker.add_accepted(st.n_accepted());
    checker.add_rejected(st.n_rejected());
    checker.add_fevals(st.n_fevals());
    checker.add_jevals(st.n_jevals());
    if constexpr (has_n_setups_method<Stepper>::value) {
      checker.add_setups(st.n_setups());
    }
    checker.set_last_order(get_stepper_order(st));
    report_profiler_if_available(st);
  }

template<class Stepper, class Checker>
inline auto transfer_stepper_diagnostics(const Stepper&, Checker&)
  -> std::enable_if_t<!has_diagnostics_counters<Stepper>::value>
  {
  }

} // namespace detail
} // namespace cppde

#endif // CPPDE_STEPPER_DIAGNOSTICS_HPP
