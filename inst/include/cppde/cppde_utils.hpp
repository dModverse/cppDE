/*
 Utility functions for ODE integration and sensitivity calculation with automatic differentiation

 This header provides:
 - Safe scalar extraction from nested AD types (cppde::dual2nd<T,N>)
 - Weighted sup-norms for step-size control
 - estimate_initial_dt: unified initial step-size estimator used by all cppDE solvers.
   Combines an HNW-style phase-1 rough h with caller-supplied ÿ (analytic via J·f + dfdt,
   or finite-difference when no Jacobian is available) and CVODE-style hlb/hub bounds.

 Copyright (C) 2026 Simon Beyer
 */

#ifndef CPPDE_UTILS_HPP
#define CPPDE_UTILS_HPP

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
#include <cppde/cppde_types.hpp>
#include <cppde/cppde_ad_traits.hpp>

namespace odeint_utils {

// =========================================================================================
//  Scalar extraction: pulled in from cppde_ad_traits.hpp
// =========================================================================================

using cppde::ad_traits::scalar_value;

// =========================================================================================
//  Weighted infinity norm (includes AD-derivative components in the weighted max)
// =========================================================================================

inline double weighted_sup_norm(
    const std::vector<double>& v,
    const std::vector<double>& x0,
    double atol,
    double rtol)
{
  double nrm = 0.0;
  for (std::size_t i = 0; i < v.size(); ++i) {
    double scale = atol + rtol * std::abs(x0[i]);
    nrm = std::max(nrm, std::abs(v[i]) / scale);
  }
  return nrm;
}

// Generic AD overload (works for cppde::dual, cppde::dual2nd).
template<class AD,
         std::enable_if_t<cppde::ad_traits::is_ad<AD>::value, int> = 0>
inline double weighted_sup_norm(
    const std::vector<AD>& v,
    const std::vector<AD>& x0,
    double atol,
    double rtol)
{
  double nrm = 0.0;
  for (std::size_t i = 0; i < v.size(); ++i) {
    auto& vi = const_cast<AD&>(v[i]);
    auto& yi = const_cast<AD&>(x0[i]);

    double vi_val = std::abs(scalar_value(vi));
    double yi_val = std::abs(scalar_value(yi));
    nrm = std::max(nrm, vi_val / (atol + rtol * yi_val));

    unsigned nd = vi.size();
    for (unsigned j = 0; j < nd; ++j) {
      double vd = std::abs(scalar_value(vi.d(j)));
      double yd = std::abs(scalar_value(yi.d(j)));
      nrm = std::max(nrm, vd / (atol + rtol * yd));
    }
  }
  return nrm;
}

// =========================================================================================
//  Weighted RMS 2-norm, matching the CVODES N_VWrmsNorm:
//  sqrt(mean((v_i*ewt_i)^2)) with ewt_i = 1/(atol + rtol*|x0_i|). Derivative
//  components of an AD value count as further entries with their own weights.
// =========================================================================================

inline double weighted_rms_norm(
    const std::vector<double>& v,
    const std::vector<double>& x0,
    double atol,
    double rtol)
{
  double sum = 0.0;
  for (std::size_t i = 0; i < v.size(); ++i) {
    double scale = atol + rtol * std::abs(x0[i]);
    double r = v[i] / scale;
    sum += r * r;
  }
  std::size_t n = v.size();
  return (n > 0) ? std::sqrt(sum / static_cast<double>(n)) : 0.0;
}

template<class AD,
         std::enable_if_t<cppde::ad_traits::is_ad<AD>::value, int> = 0>
inline double weighted_rms_norm(
    const std::vector<AD>& v,
    const std::vector<AD>& x0,
    double atol,
    double rtol)
{
  double sum = 0.0;
  std::size_t count = 0;
  for (std::size_t i = 0; i < v.size(); ++i) {
    auto& vi = const_cast<AD&>(v[i]);
    auto& yi = const_cast<AD&>(x0[i]);

    double vi_val = scalar_value(vi);
    double yi_val = std::abs(scalar_value(yi));
    double r = vi_val / (atol + rtol * yi_val);
    sum += r * r;
    ++count;

    unsigned nd = vi.size();
    for (unsigned j = 0; j < nd; ++j) {
      double vd = scalar_value(vi.d(j));
      double yd = std::abs(scalar_value(yi.d(j)));
      double rd = vd / (atol + rtol * yd);
      sum += rd * rd;
      ++count;
    }
  }
  return (count > 0) ? std::sqrt(sum / static_cast<double>(count)) : 0.0;
}

// =========================================================================================
//  cvUpperBoundH0 max-ratio: max_i |f_i| / (HUB·|y_i| + atol + rtol·|y_i|)
//
//  Matches cvUpperBoundH0 in sundials/src/cvodes/cvodes.c.  AD-aware: derivative
//  slices contribute additional (y, f) pairs, each with their own denominator.
// =========================================================================================

inline double cvhub_max_ratio(
    const std::vector<double>& f0,
    const std::vector<double>& x0,
    double hub_factor, double atol, double rtol)
{
  double m = 0.0;
  for (std::size_t i = 0; i < x0.size(); ++i) {
    double yi = std::abs(x0[i]);
    double fi = std::abs(f0[i]);
    double denom = hub_factor * yi + atol + rtol * yi;
    if (denom > 0.0) m = std::max(m, fi / denom);
  }
  return m;
}

template<class AD,
         std::enable_if_t<cppde::ad_traits::is_ad<AD>::value, int> = 0>
inline double cvhub_max_ratio(
    const std::vector<AD>& f0,
    const std::vector<AD>& x0,
    double hub_factor, double atol, double rtol)
{
  double m = 0.0;
  for (std::size_t i = 0; i < x0.size(); ++i) {
    auto& fi_ = const_cast<AD&>(f0[i]);
    auto& yi_ = const_cast<AD&>(x0[i]);
    double yi = std::abs(scalar_value(yi_));
    double fi = std::abs(scalar_value(fi_));
    double denom = hub_factor * yi + atol + rtol * yi;
    if (denom > 0.0) m = std::max(m, fi / denom);

    unsigned nd = fi_.size();
    for (unsigned j = 0; j < nd; ++j) {
      double yd = std::abs(scalar_value(yi_.d(j)));
      double fd = std::abs(scalar_value(fi_.d(j)));
      double dn = hub_factor * yd + atol + rtol * yd;
      if (dn > 0.0) m = std::max(m, fd / dn);
    }
  }
  return m;
}

// =========================================================================================
//  estimate_initial_dt: the first step, for every cppDE solver.
//
//  A rough scale from ||y0||/||f0||, a curvature from the caller's second
//  derivative, the HNW h1 formula, then an Euler-change cap above and the
//  rounding of the time variable below. See vignette("Methods"), "The first step".
// =========================================================================================

constexpr double HLB_FACTOR = 100.0;
constexpr double HUB_FACTOR = 0.1;

template<class Value, class System, class YddFn>
inline double estimate_initial_dt(
    System      system,
    YddFn       compute_ydd,
    std::vector<Value>& x0,
    Value       t0,
    double      t_final,
    double      atol,
    double      rtol,
    int         order)
{
  const std::size_t n = x0.size();
  const double t0_s = scalar_value(t0);

  // --- 1. f0 = f(t0, x0) ---
  std::vector<Value> f0(n);
  system(x0, f0, t0);

  // --- 2. phase-1 rough h0 for FD trial step (HNW §II.4) ---
  double d0 = weighted_sup_norm(x0, x0, atol, rtol);
  double d1 = weighted_sup_norm(f0, x0, atol, rtol);

  double h0;
  if (d0 < 1e-5 || d1 < 1e-5) {
    h0 = 1e-6;
  } else {
    h0 = 0.01 * d0 / d1;
  }

  // --- 3. ÿ estimate (analytic or FD: caller's choice) ---
  std::vector<Value> ydd(n);
  compute_ydd(x0, t0, f0, h0, ydd);

  double d2 = weighted_sup_norm(ydd, x0, atol, rtol);

  // --- 4. HNW phase-2 h1 formula ---
  //
  // ||y''|| stands in for ||y^(p+1)||, which is exact only at p = 1, so callers
  // pass order = 1 whatever the method. The multistep methods start at q = 1
  // anyway, and a single-step method ramps up within a few steps.
  double h1;
  double max_d = std::max(d1, d2);
  if (max_d <= 1e-15) {
    h1 = std::max(1e-6, h0 * 1e-3);
  } else {
    h1 = std::pow(0.01 / max_d, 1.0 / (order + 1));
  }

  // --- 5. Hub cap on the relative Euler change ---
  //
  // That scale is max(|y_i|, |f_i|*h0), not |y_i|: a component starting at zero
  // with a non-zero rate would otherwise collapse the denominator onto atol and
  // cap h far below what the integrator can take. State only, no AD components.
  double hub_inv = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    double yi = std::abs(scalar_value(x0[i]));
    double fi = std::abs(scalar_value(f0[i]));
    double yi_ref = std::max(yi, fi * h0);
    double denom = HUB_FACTOR * yi_ref + atol + rtol * yi;
    if (denom > 0.0) {
      hub_inv = std::max(hub_inv, fi / denom);
    }
  }
  double hub = (hub_inv > 0.0) ? 1.0 / hub_inv
                               : std::numeric_limits<double>::infinity();

  // tdist cap: no step bigger than HUB·|t_final - t0|
  double t_abs = std::abs(t0_s);
  double tdist = std::abs(t_final - t0_s);
  if (tdist > 0.0) {
    hub = std::min(hub, HUB_FACTOR * tdist);
  }

  // --- 6. hlb floor ---
  //
  // h is not clamped up to the CVODES seed HLB*eps*max(|t0|,|t_final|). That
  // quantity only seeds cvHin's refinement, whose result can be far smaller,
  // and on an AD-extended system it has to be. Only t0 + h > t0 is enforced.

  // --- 7. Combine ---
  double h = std::min({100.0 * h0, h1, hub});

  // Minimal floor: h must advance t in double precision.
  double h_tick = std::numeric_limits<double>::epsilon() * std::max(t_abs, 1.0);
  if (h < h_tick) h = h_tick;

  if (!std::isfinite(h) || h <= 0.0) {
    h = std::max(h_tick, atol);
  }
  return h;
}

// ---------------------------------------------------------------------------
//  Default FD ÿ-computer: used by explicit methods that have no Jacobian.
//
//  Computes  ÿ ≈ [f(t0 + h0, x0 + h0·f0) - f0] / h0  via one extra f-call.
//  This is the classical Hairer-Nørsett-Wanner phase-2 estimate.
// ---------------------------------------------------------------------------

template<class System>
struct fd_ydd {
  System sys;
  explicit fd_ydd(System s) : sys(s) {}

  template<class Value>
  void operator()(
      const std::vector<Value>& x0,
      Value       t0,
      const std::vector<Value>& f0,
      double      h0,
      std::vector<Value>& ydd)
  {
    const std::size_t n = x0.size();
    std::vector<Value> y1(n), f1(n);
    for (std::size_t i = 0; i < n; ++i) y1[i] = x0[i] + Value(h0) * f0[i];
    sys(y1, f1, Value(scalar_value(t0) + h0));
    for (std::size_t i = 0; i < n; ++i) ydd[i] = (f1[i] - f0[i]) / Value(h0);
  }
};

template<class System>
inline fd_ydd<System> make_fd_ydd(System sys) { return fd_ydd<System>(sys); }

// =========================================================================================
//  cppde_hin: port of CVODES cvHin (sundials/src/cvodes/cvodes.c).
//
//  Used by the multistep methods so that the first step cppDE takes is the one
//  CVODES would take on the same problem. The refinement and its two bounds are
//  in vignette("Methods"), "The first step". The WRMS norm sweeps AD components.
// =========================================================================================

template<class Value, class System>
inline double cppde_hin(
    System      system,
    std::vector<Value>& x0,
    Value       t0,
    double      t_final,
    double      atol,
    double      rtol)
{
  const std::size_t n   = x0.size();
  const double t0_s     = scalar_value(t0);
  const double eps      = std::numeric_limits<double>::epsilon();

  constexpr double HLB_F    = 100.0;
  constexpr double HUB_F    = 0.1;
  constexpr double H_BIAS   = 0.5;
  constexpr int    MAX_ITER = 4;

  // --- 1. signed / tdist / tround / too-close guard ---
  double tdiff  = t_final - t0_s;
  int    sign   = (tdiff >= 0.0) ? 1 : -1;
  double tdist  = std::abs(tdiff);
  double tround = eps * std::max(std::abs(t0_s), std::abs(t_final));

  if (tdist < 2.0 * tround) {
    // TOO_CLOSE: return a tick large enough to advance t in double.
    double h = eps * std::max(std::abs(t0_s), 1.0);
    return sign * (h > 0.0 ? h : std::numeric_limits<double>::min());
  }

  // --- 2. hlb, f0, hub ---
  double hlb = HLB_F * tround;

  std::vector<Value> f0(n);
  system(x0, f0, t0);

  double hub_inv = cvhub_max_ratio(f0, x0, HUB_F, atol, rtol);
  double hub     = HUB_F * tdist;
  if (hub * hub_inv > 1.0) hub = 1.0 / hub_inv;

  // --- 3. Short-interval / huge-rate shortcut ---
  if (hub < hlb) {
    return sign * hub;
  }

  // --- 4. Geometric-mean seed ---
  double hg   = std::sqrt(hlb * hub);
  double hnew = hg;

  // --- 5. Iterative refinement ---
  std::vector<Value> y1(n), f1(n), ydd(n);
  for (int iter = 1; iter <= MAX_ITER; ++iter) {
    double hgs = hg * sign;
    for (std::size_t i = 0; i < n; ++i) y1[i] = x0[i] + Value(hgs) * f0[i];
    system(y1, f1, Value(t0_s + hgs));
    for (std::size_t i = 0; i < n; ++i) ydd[i] = (f1[i] - f0[i]) / Value(hgs);

    double yddnrm = weighted_rms_norm(ydd, x0, atol, rtol);

    hnew = (yddnrm * hub * hub > 2.0)
         ? std::sqrt(2.0 / yddnrm)
         : std::sqrt(hg * hub);

    if (iter == MAX_ITER) break;
    double hrat = hnew / hg;
    if (hrat > 0.5 && hrat < 2.0) break;
    if (iter >= 2 && hnew > hg) { hnew = hg; break; }
    hg = hnew;
  }

  // --- 6. Bias + bounds + sign ---
  double h0 = H_BIAS * hnew;
  if (h0 < hlb) h0 = hlb;
  if (h0 > hub) h0 = hub;
  return sign * h0;
}

} // namespace odeint_utils

#endif // CPPDE_UTILS_HPP
