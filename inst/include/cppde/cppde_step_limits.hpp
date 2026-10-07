/*
 An upper bound on the step size as a function of time, piecewise constant; a
 step never exceeds the smallest bound of the intervals it overlaps. Plain
 double: the controller is not differentiated. Internal: only the tests set it,
 through attr(times, "hmax"), to force a step sequence.

 Copyright (C) 2026 Simon Beyer
 */

#ifndef CPPDE_STEP_LIMITS_HPP
#define CPPDE_STEP_LIMITS_HPP

#include <algorithm>
#include <cstddef>
#include <limits>
#include <vector>

namespace cppde {

// Interval k starts at t[k] and holds the bound h[k], up to t[k + 1]; the last
// one reaches to the end of the run, times before t[0] take h[0].
struct step_limits {
  std::vector<double> t, h;

  bool empty() const { return t.empty(); }

  // The largest step from `t0` in direction `dir` that stays within the bound
  // of every interval it reaches into.
  double cap(double t0, double dir) const
  {
    if (t.empty()) return std::numeric_limits<double>::infinity();
    const std::size_t n = t.size();
    auto at = [&](double s) {
      const auto it = std::upper_bound(t.begin(), t.end(), s);
      return it == t.begin() ? std::size_t(0)
                             : static_cast<std::size_t>(it - t.begin()) - 1;
    };
    std::size_t k = at(t0);
    double c = h[k];
    // Later intervals the step would reach into, forward runs only.
    if (dir > 0)
      for (std::size_t j = k + 1; j < n && t[j] < t0 + c; ++j) c = std::min(c, h[j]);
    return c;
  }
};

// Where the dense loop looks for the bounds, per thread and per solve. Null
// means no bound. A pointer, so the sink costs nothing unset.
inline const step_limits*& step_limit_sink() {
  thread_local const step_limits* p = nullptr;
  return p;
}

struct step_limit_scope {
  const step_limits* prev;
  explicit step_limit_scope(const step_limits& s) : prev(step_limit_sink()) {
    step_limit_sink() = s.empty() ? nullptr : &s;
  }
  ~step_limit_scope() { step_limit_sink() = prev; }
  step_limit_scope(const step_limit_scope&) = delete;
  step_limit_scope& operator=(const step_limit_scope&) = delete;
};

} // namespace cppde

#endif  // CPPDE_STEP_LIMITS_HPP
