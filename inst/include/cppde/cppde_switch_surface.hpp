/**
 * @file cppde_switch_surface.hpp
 * @brief Where a state switch stands on its surface, read the same way by the
 * native event engine and the generated CVODE source.
 *
 * Copyright (C) 2026 Simon Beyer
 */

#ifndef CPPDE_SWITCH_SURFACE_HPP
#define CPPDE_SWITCH_SURFACE_HPP

#include <cmath>
#include <cstddef>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>

namespace cppde {
namespace detail {

// A value is zero where it vanishes against the rounding of the terms it is
// summed from, `terms` being the sum of their sizes.
inline bool negligible(double v, double terms) {
  return std::abs(v) <= 64.0 * std::numeric_limits<double>::epsilon() * terms;
}

// The sign of r', zero where it is negligible.
inline int rate_sign(double rdot, double terms) {
  if (negligible(rdot, terms)) return 0;
  return rdot > 0.0 ? 1 : -1;
}

enum class surface_side { keep, on, off, slide };

// The side a switch on its surface takes, from the signs of r' with its mode
// off and on: the one both branches move to, a slide where both move onto the
// surface, and the current mode where they part or one stays on it.
inline surface_side surface_choice(int off, int on) {
  if (on < 0 && off > 0) return surface_side::slide;
  if (on > 0 && off > 0) return surface_side::on;
  if (on < 0 && off < 0) return surface_side::off;
  return surface_side::keep;
}

// The weight of combination c of the slides in a, the product over slide i of
// a_i where bit i of c is set and 1 - a_i where it is not; with `d` >= 0 its
// derivative in a_d.
template<class V>
V combo_weight(const std::vector<V>& a, std::size_t c, int d = -1) {
  V w(1.0);
  for (std::size_t i = 0; i < a.size(); ++i) {
    const bool on = (c >> i) & 1u;
    if (static_cast<int>(i) == d) w = w * V(on ? 1.0 : -1.0);
    else w = w * (on ? a[i] : V(1.0) - a[i]);
  }
  return w;
}

// x = A^{-1} b for a small dense A, row-major, by elimination with pivoting on
// the values `val` reads.
template<class V, class Val>
void small_solve(std::vector<V> A, std::vector<V>& b, Val val) {
  const std::size_t n = b.size();
  for (std::size_t k = 0; k < n; ++k) {
    std::size_t p = k;
    for (std::size_t i = k + 1; i < n; ++i)
      if (std::abs(val(A[i * n + k])) > std::abs(val(A[p * n + k]))) p = i;
    if (p != k) {
      for (std::size_t j = 0; j < n; ++j) std::swap(A[k * n + j], A[p * n + j]);
      std::swap(b[k], b[p]);
    }
    for (std::size_t i = k + 1; i < n; ++i) {
      const V l = A[i * n + k] / A[k * n + k];
      for (std::size_t j = k; j < n; ++j) A[i * n + j] = A[i * n + j] - l * A[k * n + j];
      b[i] = b[i] - l * b[k];
    }
  }
  for (std::size_t k = n; k-- > 0;) {
    V s = b[k];
    for (std::size_t j = k + 1; j < n; ++j) s = s - A[k * n + j] * b[j];
    b[k] = s / A[k * n + k];
  }
}

inline void small_solve(const std::vector<double>& A, std::vector<double>& b) {
  small_solve(A, b, [](double v) { return v; });
}

}  // namespace detail
}  // namespace cppde

#endif  // CPPDE_SWITCH_SURFACE_HPP
