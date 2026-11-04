/*
 Stepper traits for cppDE: compile-time dispatch between
 single-step (Rosenbrock, Tsit5) and multi-step methods.

 Copyright (C) 2026 Simon Beyer

 This header defines stepper_traits<Stepper>, a type-trait struct that
 the EventEngine uses to adapt its behaviour to the stepper family:

 - needs_restart_after_event: true if the stepper must discard its
                              history and restart at order 1 after
                              a state discontinuity (events)

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
  // No restart needed after events: single-step methods are self-starting
  static constexpr bool needs_restart_after_event = false;
};
// ============================================================================
//  Trait extraction helpers (for use in if-constexpr / enable_if)
// ============================================================================

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

// Specialisation for any type with the multistepper tag
template<class Stepper>
struct stepper_traits<Stepper,
                     std::enable_if_t<has_multistepper_tag<Stepper>::value>>
                     {
                       static constexpr bool needs_restart_after_event = true;
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

// Recursive version, chasing through several wrapper layers. The recursion
// stops when T::stepper_type is T itself, which would otherwise not terminate.
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
                          static constexpr bool needs_restart_after_event = true;
                        };

} // namespace cppde

#endif // CPPDE_STEPPER_TRAITS_HPP
