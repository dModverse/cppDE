/*
 Switching times of the right-hand side.

 A condition of the model on time and parameters alone, such as the bounds of
 piecewise(1, time > ts && time <= t2, 0), switches f at times the code
 generator solves for before the solve. The integration loops stop in front of
 each of them inside the grid and cross it as a jump of f in t, see
 cppde_event_engine.hpp; nothing is observed there.

 The generated solve hands its times over per thread, since a batch runs
 several solves at once. Plain double and free of R.

 Copyright (C) 2026 Simon Beyer
 */

#ifndef CPPDE_SWITCH_TIMES_HPP
#define CPPDE_SWITCH_TIMES_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

namespace cppde {

struct switch_times {
  const double* t = nullptr;
  std::size_t   n = 0;
};

// Where the integration loops look for the switching times of the running
// solve. Trivially destructible, see cppde_tls.hpp.
inline switch_times& switch_time_sink() {
  thread_local switch_times s;
  return s;
}

struct switch_time_scope {
  switch_times prev;
  switch_time_scope(const double* t, std::size_t n) : prev(switch_time_sink()) {
    switch_time_sink() = switch_times{t, n};
  }
  ~switch_time_scope() { switch_time_sink() = prev; }
  switch_time_scope(const switch_time_scope&) = delete;
  switch_time_scope& operator=(const switch_time_scope&) = delete;
};

namespace detail {

// How far in front of a switching time a step ends, and how far past it the
// crossing reaches: a few doubles of the grid's largest time, which covers the
// rounding of a solved switching time and where a shortened step lands.
inline double switch_margin(double t_first, double t_last) {
  return 16.0 * std::numeric_limits<double>::epsilon() *
         std::max(std::abs(t_first), std::abs(t_last));
}

// The switching times strictly inside the grid, ordered in the direction of
// integration, of which those closer than the margin to the one before or to
// either end are dropped. Non-finite times are dropped too.
inline std::vector<double> switch_times_in_window(const switch_times& s,
                                                  double t_first, double t_last) {
  std::vector<double> out;
  const double dir = t_last >= t_first ? 1.0 : -1.0;
  const double m = switch_margin(t_first, t_last);
  for (std::size_t i = 0; i < s.n; ++i) {
    const double ts = s.t[i];
    if (std::isfinite(ts) && dir * (ts - t_first) > 2.0 * m &&
        dir * (t_last - ts) > 2.0 * m)
      out.push_back(ts);
  }
  std::sort(out.begin(), out.end(),
            [dir](double a, double b) { return dir * a < dir * b; });
  std::vector<double> kept;
  for (double ts : out)
    if (kept.empty() || dir * (ts - kept.back()) > m) kept.push_back(ts);
  return kept;
}

} // namespace detail
} // namespace cppde

#endif // CPPDE_SWITCH_TIMES_HPP
