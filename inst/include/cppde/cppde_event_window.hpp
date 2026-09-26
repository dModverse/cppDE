/*
 Which fixed-time events a solve applies.

 A solve runs from the first to the last time of its grid, zero included when
 the model injects it. A fixed event takes effect at the first time and inside
 the grid, and never at or beyond the last time. Two solves over [t0, t1] and
 [t1, t2], the second started from the state the first ends on, then apply
 every event exactly once and reproduce one solve over [t0, t2]; an event at
 t1 belongs to the second. An event outside the window adds no output row and
 moves neither end of the integration. A grid of a single time applies the
 events at that time, as the start of a longer grid would.

 Plain double and free of R, so the solvers, both code generators and the
 batch preallocation read the same rule.

 Copyright (C) 2026 Simon Beyer
 */

#ifndef CPPDE_EVENT_WINDOW_HPP
#define CPPDE_EVENT_WINDOW_HPP

#include <cmath>

namespace cppde {
namespace detail {

// Two times closer than this are the same time, the tolerance the event
// engine matches an event time against an output time with.
constexpr double event_time_tol = 1e-14;

inline bool same_event_time(double a, double b) {
  return std::abs(a - b) < event_time_tol;
}

// Whether a fixed event at t_e fires on a grid running from t_first to t_last.
inline bool fixed_event_in_window(double t_e, double t_first, double t_last) {
  if (same_event_time(t_e, t_first)) return true;
  return t_e > t_first && t_e < t_last && !same_event_time(t_e, t_last);
}

} // namespace detail
} // namespace cppde

#endif // CPPDE_EVENT_WINDOW_HPP
