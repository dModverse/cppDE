/*
 lambda as a step-size weight: the goal-oriented half of the error norm.

 A grid chosen from the state error alone need not be good enough for the
 adjoint. lambda solves backwards, and where it grows a step's contribution to
 the gradient error is large even where its own state error is small: the order
 carries over to the gradient, the constant does not.

 lambda at step k is not available while step k runs, because it depends on
 every step after it, so the forward pass cannot compute this term itself. It
 is available across *runs*, which is what this header carries: a sweep leaves
 lambda on its own grid, and the next run reads it back as a weight.

 Two properties make that sound rather than merely cheap.

   err = max(err_state, err_lambda)

 is the shape the dual path already uses. A max can only shrink a step, so the
 grid stays at least as fine as atol/rtol demand and is finer only where lambda
 asks. The trajectory therefore keeps its own accuracy, and a wrong weight costs
 time and never accuracy: too large means needlessly small steps, too small
 means the state term governs as before. That is what makes an *estimated*
 lambda, from a previous run at a nearby theta, admissible at all.

 Plain double throughout, and deliberately: this is a controller input, and the
 controller is not differentiated. See dev/adjoint-plan.md, stage 9.

 Copyright (C) 2026 Simon Beyer
 */

#ifndef CPPDE_ERR_WEIGHTS_HPP
#define CPPDE_ERR_WEIGHTS_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace cppde {

// ============================================================================
//  err_weights
//
//  lambda sampled at times t[0..n), n_x components each, row-major:
//  lam[i * n_x + j] is component j at t[i]. Read back by linear interpolation.
//
//  Linear and not PCHIP, for three reasons. A forcing is AD-aware because a
//  forcing belongs on the tape, and this must not; lambda jumps where the
//  objective seeds it and where an event resets the state, so an interpolant
//  spanning those points would smear the jump; and shape preservation buys
//  nothing for a weight, whose overshoots cost only time under the max.
//
//  `breaks` names sample indices the interpolant must not span. A query inside
//  a broken interval takes the nearer endpoint rather than a blend of two
//  values that belong to different sides of a jump.
// ============================================================================

class err_weights {
public:
  void clear() {
    m_t.clear(); m_lam.clear(); m_breaks.clear(); m_n_x = 0;
  }

  // Samples must be ascending in t. n_x is the state count, and lam.size() has
  // to be t.size() * n_x or the weights are dropped rather than half-read.
  void set(std::vector<double> t, std::vector<double> lam, std::size_t n_x,
           std::vector<std::size_t> breaks = {})
  {
    if (n_x == 0 || t.empty() || lam.size() != t.size() * n_x) { clear(); return; }
    m_t = std::move(t);
    m_lam = std::move(lam);
    m_breaks = std::move(breaks);
    m_n_x = n_x;
    std::sort(m_breaks.begin(), m_breaks.end());
  }

  bool empty() const { return m_n_x == 0 || m_t.empty(); }
  std::size_t n_states() const { return m_n_x; }
  std::size_t n_samples() const { return m_t.size(); }

  // The objective's own tolerance. Without it the two halves of the norm are
  // not comparable, and a large lambda would make the weighted term govern
  // everywhere, which is only atol/rtol set tighter, by a route that hides
  // what it did.
  double gradtol() const { return m_gradtol; }
  void gradtol(double g) { if (g > 0.0) m_gradtol = g; }

  // Weights below this fraction of the largest are lifted to it. lambda is the
  // *linearised* influence of a state, so a state it weights near zero can
  // still drift nonlinearly; under the max this is caution and not necessity.
  double floor() const { return m_floor; }
  void floor(double f) { if (f >= 0.0 && f < 1.0) m_floor = f; }

  // lambda at time t, clamped to the sample range at either end.
  void at(double t, std::vector<double>& out) const {
    out.assign(m_n_x, 0.0);
    if (empty()) return;

    const std::size_t n = m_t.size();
    if (n == 1 || t <= m_t.front()) { copy_row(0, out); return; }
    if (t >= m_t.back())            { copy_row(n - 1, out); return; }

    // First sample strictly greater than t, so the query sits in [hi-1, hi).
    const std::size_t hi =
        static_cast<std::size_t>(
            std::upper_bound(m_t.begin(), m_t.end(), t) - m_t.begin());
    const std::size_t lo = hi - 1;

    if (spans_break(lo, hi)) {
      // Two values from opposite sides of a jump do not average into anything;
      // the nearer sample is the honest answer.
      copy_row((t - m_t[lo] <= m_t[hi] - t) ? lo : hi, out);
      return;
    }

    const double h = m_t[hi] - m_t[lo];
    const double w = (h > 0.0) ? (t - m_t[lo]) / h : 0.0;
    for (std::size_t j = 0; j < m_n_x; ++j)
      out[j] = (1.0 - w) * m_lam[lo * m_n_x + j] + w * m_lam[hi * m_n_x + j];
    apply_floor(out);
  }

private:
  void copy_row(std::size_t i, std::vector<double>& out) const {
    for (std::size_t j = 0; j < m_n_x; ++j) out[j] = m_lam[i * m_n_x + j];
    apply_floor(out);
  }

  void apply_floor(std::vector<double>& out) const {
    if (m_floor <= 0.0) return;
    double mx = 0.0;
    for (double v : out) mx = std::max(mx, std::abs(v));
    if (mx <= 0.0) return;
    const double lo = m_floor * mx;
    for (double& v : out)
      if (std::abs(v) < lo) v = (v < 0.0) ? -lo : lo;
  }

  // Whether a break index sits strictly inside (lo, hi]: a jump at hi belongs
  // to the interval that ends there.
  bool spans_break(std::size_t lo, std::size_t hi) const {
    if (m_breaks.empty()) return false;
    auto it = std::lower_bound(m_breaks.begin(), m_breaks.end(), lo + 1);
    return it != m_breaks.end() && *it <= hi;
  }

  std::vector<double>      m_t, m_lam;
  std::vector<std::size_t> m_breaks;
  std::size_t              m_n_x = 0;
  double                   m_gradtol = 1e-6;
  double                   m_floor   = 0.0;
};

// Where the controllers look for the weights. Per-thread and per-solve, the way
// the step trace's sink is, and for the same reason: a batch runs several
// solves at once and they do not share a grid. Null means "no weighting", which
// is the shipped state.
//
// Pointer, not thread_local object: see cppde_tls.hpp.
inline const err_weights*& err_weight_sink() {
  thread_local const err_weights* p = nullptr;
  return p;
}

struct err_weight_scope {
  const err_weights* prev;
  explicit err_weight_scope(const err_weights& w) : prev(err_weight_sink()) {
    err_weight_sink() = &w;
  }
  ~err_weight_scope() { err_weight_sink() = prev; }
  err_weight_scope(const err_weight_scope&) = delete;
  err_weight_scope& operator=(const err_weight_scope&) = delete;
};

namespace detail {

// The weighted half of the error norm: |lambda(t)' e| / gradtol, the step's
// contribution to the error in the objective measured against its own
// tolerance. Zero when no weights are set, so the caller's max is unchanged.
//
// `get` scalarises, because an AD run weights the value layer: the tangent
// columns have their own term in the norm already.
template<class V, class Get>
inline double weighted_error(const std::vector<V>& xerr, double t, Get get) {
  const err_weights* w = err_weight_sink();
  if (w == nullptr || w->empty()) return 0.0;
  static thread_local std::vector<double>* buf = nullptr;
  if (buf == nullptr) buf = new std::vector<double>();  // leaked, see cppde_tls.hpp
  w->at(t, *buf);
  const std::size_t n = std::min(xerr.size(), buf->size());
  double s = 0.0;
  for (std::size_t i = 0; i < n; ++i)
    s += (*buf)[i] * static_cast<double>(get(xerr[i]));
  return std::abs(s) / w->gradtol();
}

}  // namespace detail
}  // namespace cppde

#endif  // CPPDE_ERR_WEIGHTS_HPP
