/*
 Stepper traits for cppDE: compile-time dispatch between
 single-step (Rosenbrock) and multi-step methods.

 Copyright (C) 2026 Simon Beyer

 This header defines stepper_traits<Stepper>, a type-trait struct that
 the EventEngine, integrate_times, and integrate_times_dense use to
 adapt their behaviour to different stepper families:

 - is_multistep:              true for multistep steppers (NDF/BDF/Adams
                              and variants: all instantiations of
                              cppde::multistepper)
 - needs_restart_after_event: true if the stepper must discard its
                              history and restart at order 1 after
                              a state discontinuity (events)
 - min_order / max_order:     order range of the method

 The default (unspecialised) traits match Rosenbrock4 and other
 single-step methods.  The multistepper specialisation is matched via
 a nested tag typedef: see below.

 Distributed under the Boost Software License, Version 1.0.
 */

#ifndef CPPDE_STEPPER_TRAITS_HPP
#define CPPDE_STEPPER_TRAITS_HPP

#include <type_traits>

namespace cppde {

// ============================================================================
//  Primary template: defaults for single-step methods (Rosenbrock, etc.)
// ============================================================================

template<class Stepper, class = void>
struct stepper_traits {
  // Single-step method: no history to manage
  static constexpr bool is_multistep = false;

  // No restart needed after events: single-step methods are self-starting
  static constexpr bool needs_restart_after_event = false;

  // Fixed order (not applicable for single-step, but provided for uniformity)
  static constexpr int min_order = 0;
  static constexpr int max_order = 0;
};
// ============================================================================
//  Trait extraction helpers (for use in if-constexpr / enable_if)
// ============================================================================

// True if the stepper is a multi-step method (NDF, BDF, Adams, etc.)
template<class Stepper>
inline constexpr bool is_multistep_v = stepper_traits<Stepper>::is_multistep;

// True if the stepper needs a full restart (history discard + order 1)
// after a state discontinuity (event)
template<class Stepper>
inline constexpr bool needs_restart_after_event_v =
stepper_traits<Stepper>::needs_restart_after_event;
// ============================================================================
//  Multistepper detection via tag
//
//  cppde::multistepper declares a nested is_multistepper_tag typedef, detected
//  here by SFINAE, which covers both bdf and adams.
// ============================================================================

// SFINAE detector for the multistepper tag
template<class T, class = void>
struct has_multistepper_tag : std::false_type {};

template<class T>
struct has_multistepper_tag<T, std::void_t<typename T::is_multistepper_tag>>
: std::true_type {};

// Specialisation for any type carrying the multistepper tag
template<class Stepper>
struct stepper_traits<Stepper,
                     std::enable_if_t<has_multistepper_tag<Stepper>::value>>
                     {
                       static constexpr bool is_multistep = true;
                       static constexpr bool needs_restart_after_event = true;
                       static constexpr int  min_order = 1;
                       static constexpr int  max_order = 5;
                     };
// ============================================================================
//  Trait propagation through wrapper layers
//
//  A wrapper exposing a stepper_type typedef inherits that stepper's traits
//  unless it specialises stepper_traits itself, so a caller can query the
//  outermost wrapper without knowing the nesting depth.
// ============================================================================

// SFINAE detector for nested stepper_type
template<class T, class = void>
struct has_inner_stepper_type : std::false_type {};

template<class T>
struct has_inner_stepper_type<T, std::void_t<typename T::stepper_type>>
: std::true_type {};

// Detector for multistepper tag on inner stepper (through wrapper chain)
template<class T, class = void>
struct inner_has_multistepper_tag : std::false_type {};

template<class T>
struct inner_has_multistepper_tag<T,
                         std::enable_if_t<
                           has_inner_stepper_type<T>::value &&
                           !has_multistepper_tag<T>::value &&
                           !std::is_same<T, typename T::stepper_type>::value>>  // guard against self-referential stepper_type
                           : has_multistepper_tag<typename T::stepper_type> {};

// Recursive version, chasing through several wrapper layers. The recursion
// stops when T::stepper_type is T itself, which some Boost steppers declare and
// which would otherwise not terminate.
template<class T, class = void>
struct deep_has_multistepper_tag : has_multistepper_tag<T> {};

template<class T>
struct deep_has_multistepper_tag<T,
                        std::enable_if_t<
                          has_inner_stepper_type<T>::value &&
                          !has_multistepper_tag<T>::value &&
                          !std::is_same<T, typename T::stepper_type>::value>>
                          : deep_has_multistepper_tag<typename T::stepper_type> {};

// Propagated traits for wrapper types (controller, dense output)
// that wrap a multistepper somewhere in their chain
template<class Wrapper>
struct stepper_traits<Wrapper,
                      std::enable_if_t<
                        !has_multistepper_tag<Wrapper>::value &&
                        deep_has_multistepper_tag<Wrapper>::value>>
                        {
                          static constexpr bool is_multistep = true;
                          static constexpr bool needs_restart_after_event = true;
                          static constexpr int  min_order = 1;
                          static constexpr int  max_order = 5;
                        };

} // namespace cppde

#endif // CPPDE_STEPPER_TRAITS_HPP
