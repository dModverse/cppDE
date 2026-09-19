/*
 lambda as a step-size weight: the goal-oriented half of the error norm.

 A sweep leaves lambda on its own grid, and the next run reads it back as a
 weight in err = max(err_state, err_lambda). The max only shrinks steps, so a
 wrong weight costs time and never accuracy. Plain double: the controller is
 not differentiated. See vignette("Methods"), "Goal-oriented step-size control".

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
//  `breaks` names sample indices where lambda jumps (seeds, events), which the
//  interpolant must not span: a query inside a broken interval takes the
//  nearer endpoint.
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

  // The objective's own tolerance, which puts the weighted half of the norm on
  // the same scale as the state half.
  double gradtol() const { return m_gradtol; }
  void gradtol(double g) { if (g > 0.0) m_gradtol = g; }

  // Weights below this fraction of the largest are lifted to it, so a state
  // the linearised lambda weights near zero still counts.
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
      // Two values from opposite sides of a jump are not averaged; the nearer
      // sample is used.
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

// Where the controllers look for the weights, per thread and per solve because
// a batch runs several solves at once. Null means "no weighting". A pointer,
// not a thread_local object: see cppde_tls.hpp.
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

// The weighted half of the error norm, |lambda(t)' e| / gradtol, and zero
// without weights. `get` scalarises: an AD run weights the value layer only,
// the tangent columns having their own term in the norm.
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
