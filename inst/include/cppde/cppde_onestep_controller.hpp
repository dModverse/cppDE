/*
 One-step controller with Gustafsson–Söderlind PI step-size control.

 Generic PI controller for any single-step stepper satisfying the
 do_step(sys, x, t, xout, dt, xerr [, hint]) interface.  Works with
 both double and AD types.  Error norm includes derivative components
 for AD types, equivalent to solving the sensitivity-augmented system.

 The error_order of the embedded pair is read from Stepper::error_order
 and drives the PI gain defaults.

 Originally derived from Boost.Odeint's rosenbrock4_controller
 by Karsten Ahnert, Mario Mulansky, and Christoph Koke (2011-2012),
 distributed under the Boost Software License, Version 1.0.

 Modified work:
 Copyright (C) 2026 Simon Beyer

 Modifications:
 - Generalized to any single-step stepper (rosenbrock4, tsit5, ...)
 - Error order derived from Stepper::error_order
 - SFINAE-based Jacobian hint dispatch
 - Added reset_after_event() for event-driven integration
 - Implemented Gustafsson–Söderlind PI control algorithm
 - Unified double and AD handling in a single class

 PI controller based on:
 Gustafsson, K., Lundh, M. & Söderlind, G. (1988).
 "A PI stepsize control for the numerical solution of ordinary differential equations".
 BIT 28, 270–287. https://doi.org/10.1007/BF01934091
 */

#ifndef CPPDE_ONESTEP_CONTROLLER_HPP
#define CPPDE_ONESTEP_CONTROLLER_HPP

#include <cmath>
#include <algorithm>
#include <cstdio>
#include <type_traits>
#include <vector>
#include <cppde/cppde_types.hpp>
#include <cppde/cppde_ad_lu.hpp>   // for ad_lu::is_ad
#include <cppde/cppde_ad_traits.hpp>
#include <cppde/cppde_dual_slab.hpp>
#include <cppde/cppde_profiler.hpp>
#include <cppde/cppde_err_weights.hpp>
#include <cppde/cppde_step_trace.hpp>

namespace cppde {

// ============================================================================
//  Scalar value extraction: pulled in from cppde_ad_traits.hpp
// ============================================================================

namespace controller_detail {
using cppde::ad_traits::scalar_value;
} // namespace controller_detail

// ============================================================================
//  SFINAE: detect has_valid_jacobian() on the stepper (implicit methods)
// ============================================================================

namespace onestep_detail {

template<class S, class = void>
struct has_jacobian_api : std::false_type {};

template<class S>
struct has_jacobian_api<S, std::void_t<decltype(std::declval<const S&>().has_valid_jacobian())>>
: std::true_type {};

// ----------------------------------------------------------------------------
//  The error norm and the control law, written on the value type.
//
//  Templated and spelled with the cppde:: math names so one statement of each
//  formula serves every value type the forward run instantiates. The reverse
//  replay does not call them: the grid is read off the checkpoints and the
//  controller stays off the tape.
//
//  wrms_state takes an accessor because an AD number has to be scalarised down
//  to its value before it enters the norm.
// ----------------------------------------------------------------------------

template<class V, class Get>
inline auto wrms_state(const std::vector<V>& x, const std::vector<V>& xold,
                       const std::vector<V>& xerr,
                       double atol, double rtol, Get get)
  -> decltype(get(xerr[0]))
{
  using S = decltype(get(xerr[0]));
  const size_t n = xerr.size();
  S sumsq = S(0.0);
  for (size_t i = 0; i < n; ++i) {
    S scale = atol + rtol * cppde::max(cppde::abs(get(xold[i])),
                                       cppde::abs(get(x[i])));
    S r = cppde::abs(get(xerr[i])) / scale;
    sumsq = sumsq + r * r;
  }
  return cppde::sqrt(sumsq / static_cast<double>(n));
}

// Step-size factor after an accepted step. pi_form is false on the first step
// and after a rejection, where err_old carries nothing usable; cap_at_one holds
// the step back after a rejection.
template<class V>
inline V accept_factor(const V& err, const V& err_old, bool pi_form, bool cap_at_one,
                       double order, double alpha, double beta, double safety,
                       double min_factor, double max_factor)
{
  V factor = pi_form
    ? safety * cppde::pow(err_old / err, beta) * cppde::pow(V(1.0) / err, alpha)
    : safety * cppde::pow(V(1.0) / err, 1.0 / (order + 1.0));
  if (cap_at_one) factor = cppde::min(factor, 1.0);
  return cppde::min(cppde::max(factor, min_factor), max_factor);
}

// Step-size factor after a rejected attempt. Never grows the step.
template<class V>
inline V reject_factor(const V& err, double order, double safety, double min_factor)
{
  V factor = safety * cppde::pow(V(1.0) / err, 1.0 / (order + 1.0));
  return cppde::min(cppde::max(factor, min_factor), 0.9);
}

} // namespace onestep_detail

// ============================================================================
//  onestep_controller<Stepper>
//
//  Generic PI controller for single-step ODE methods.
//  Error norm includes derivative components for AD types.
// ============================================================================

template<class Stepper>
class onestep_controller
{
public:

  using stepper_type       = Stepper;
  using state_type         = typename stepper_type::state_type;
  using value_type         = typename stepper_type::value_type;
  using time_type          = typename stepper_type::value_type;   // external: value_type for API compat
  using scalar_time_type   = typename stepper_type::time_type;     // internal: scalar for stepper
  using wrapped_state_type = typename stepper_type::wrapped_state_type;
  using resizer_type       = typename stepper_type::resizer_type;

  using stepper_category = controlled_stepper_tag;
  using controller_type  = onestep_controller<Stepper>;

  // Error order from the embedded pair (drives PI gain defaults)
  static constexpr double order = static_cast<double>(stepper_type::error_order);

  // Default proportional gain: alpha = 0.7 / (order + 1)
  static constexpr double default_alpha = 0.7 / (order + 1.0);

  // Default integral gain: beta = 0.4 / (order + 1)
  static constexpr double default_beta  = 0.4 / (order + 1.0);

  // Default safety factor for step-size selection
  static constexpr double default_safety = 0.9;

  // Default maximum step increase factor
  static constexpr double default_max_factor = 5.0;

  // Default minimum step decrease factor
  static constexpr double default_min_factor = 0.2;

  // ====================================================================
  //  Constructors
  // ====================================================================

  onestep_controller(
    double atol = 1.0e-6,
    double rtol = 1.0e-6,
    double alpha = default_alpha,
    double beta = default_beta,
    double safety = default_safety,
    double max_factor = default_max_factor,
    double min_factor = default_min_factor)
    : m_stepper()
    , m_atol(atol), m_rtol(rtol)
    , m_max_dt(0)
    , m_alpha(alpha), m_beta(beta)
    , m_safety(safety)
    , m_max_factor(max_factor), m_min_factor(min_factor)
    , m_first_step(true), m_err_old(1.0), m_dt_old(1.0)
    , m_last_rejected(false)
    , m_n_accepted(0), m_n_rejected(0)
  {}

  onestep_controller(
    double atol, double rtol, double max_dt,
    double alpha, double beta,
    double safety = default_safety,
    double max_factor = default_max_factor,
    double min_factor = default_min_factor)
    : m_stepper()
    , m_atol(atol), m_rtol(rtol)
    , m_max_dt(max_dt)
    , m_alpha(alpha), m_beta(beta)
    , m_safety(safety)
    , m_max_factor(max_factor), m_min_factor(min_factor)
    , m_first_step(true), m_err_old(1.0), m_dt_old(1.0)
    , m_last_rejected(false)
    , m_n_accepted(0), m_n_rejected(0)
  {}

  // Move-only, as in the multistep family: copying the slab members would build
  // a fresh tangent block while the duals still point at the original. A vector
  // move preserves data(), so the embedded pointers stay valid.
  onestep_controller(const onestep_controller&)            = delete;
  onestep_controller& operator=(const onestep_controller&) = delete;
  onestep_controller(onestep_controller&&)                 = default;
  onestep_controller& operator=(onestep_controller&&)      = default;

  // ====================================================================
  //  Jacobian hint: SFINAE dispatch
  //
  //  For implicit steppers (rosenbrock4) that expose has_valid_jacobian(),
  //  reuse the Jacobian on rejection.  For explicit steppers (tsit5),
  //  always return recompute_all (the hint is ignored by the stepper).
  // ====================================================================

  jacobian_hint compute_hint() const
  {
    if constexpr (onestep_detail::has_jacobian_api<stepper_type>::value) {
      if (m_first_step || !m_stepper.has_valid_jacobian())
        return jacobian_hint::recompute_all;
      if (m_last_rejected)
        return jacobian_hint::reuse_jacobian;
      return jacobian_hint::recompute_all;
    } else {
      return jacobian_hint::recompute_all;
    }
  }

  // ====================================================================
  //  Error norm: AD-aware WRMS, the maximum over the state and each sensitivity
  //  vector, the CVODES cvSensUpdateNorm convention. For a non-AD type the
  //  derivative loop compiles out.
  // ====================================================================

  // t_end is where the step lands, which is where the lambda weight is read.
  // Defaulted so a caller that sets no weights need not know about them.
  double error(const state_type& x, const state_type& xold, const state_type& xerr,
               double t_end = 0.0)
  {
    auto _tp = m_prof.timer(cppde::prof_cat::error_norm);
    using controller_detail::scalar_value;

    const size_t n = x.size();
    if (n == 0) return 0.0;

    // The value half. Stage 9 adds a lambda-weighted term beside the
    // sensitivity one below, under the same max.
    double max_norm = onestep_detail::wrms_state(
        x, xold, xerr, m_atol, m_rtol,
        [](const value_type& v) { return scalar_value(v); });

    if constexpr (ad_lu::is_ad<value_type>::value) {
      const unsigned nd = m_sens_err_con
        ? const_cast<value_type&>(xerr[0]).size() : 0u;
      std::vector<double> sens_sumsq(nd, 0.0);
      for (size_t i = 0; i < n; ++i) {
        auto& xerr_ad = const_cast<value_type&>(xerr[i]);
        auto& xold_ad = const_cast<value_type&>(xold[i]);
        auto& xnew_ad = const_cast<value_type&>(x[i]);
        for (unsigned j = 0; j < nd; ++j) {
          double err_d = std::abs(scalar_value(xerr_ad.d(j)));
          double xold_d = std::abs(scalar_value(xold_ad.d(j)));
          double xnew_d = std::abs(scalar_value(xnew_ad.d(j)));
          double scale_d = m_atol + m_rtol * std::max(xold_d, xnew_d);
          double rd = err_d / scale_d;
          sens_sumsq[j] += rd * rd;
        }
      }
      for (unsigned j = 0; j < nd; ++j) {
        double sens_norm = std::sqrt(sens_sumsq[j] / n);
        if (sens_norm > max_norm) max_norm = sens_norm;
      }
    }

    // The goal-oriented term, stage 9: |lambda(t)' e| / gradtol, this step's
    // share of the error in the objective. Zero unless a sweep left weights,
    // and under the same max, so it can only make a step smaller.
    const double lam_norm = ::cppde::detail::weighted_error(
        xerr, t_end, [](const value_type& v) { return scalar_value(v); });
    if (lam_norm > max_norm) max_norm = lam_norm;
    return max_norm;
  }

  // ====================================================================
  //  try_step: in-place
  // ====================================================================

  template<class System>
  controlled_step_result
  try_step(System& sys, state_type& x, time_type& t, time_type& dt)
  {
    m_xnew_resizer.adjust_size(x, [this](const state_type& s){ return this->resize_m_xnew(s); });
    controlled_step_result res = try_step(sys, x, t, m_xnew.m_v, dt);
    if (res == success) {
      std::swap(x, m_xnew.m_v);
    }
    return res;
  }

  // ====================================================================
  //  try_step: separate input/output: main entry point
  // ====================================================================

  template<class System>
  controlled_step_result
  try_step(System& sys, const state_type& x, time_type& t, state_type& xout, time_type& dt)
  {
    if (m_max_dt != 0.0) {
      double dt_val = controller_detail::scalar_value(dt);
      if (std::abs(dt_val) > m_max_dt) {
        dt = time_type(std::copysign(m_max_dt, dt_val));
        return fail;
      }
    }

    m_xerr_resizer.adjust_size(x, [this](const state_type& s){ return this->resize_m_xerr(s); });

    jacobian_hint hint = compute_hint();

    m_stepper.do_step(sys, x, t, xout, dt, m_xerr.m_v, hint);
    double err = error(xout, x, m_xerr.m_v,
                       controller_detail::scalar_value(t) +
                       controller_detail::scalar_value(dt));

    // Prevent division by zero
    err = std::max(err, 1e-15);

#ifdef CPPDE_STEP_TRACE
    {
  // The row uses the schema of the multistepper trace, so both step kinds land
  // in one buffer. Fields that do not apply to a single-step method get neutral
  // defaults, tq2 = 1 so that dsm equals err.
      ::cppde::ndf_detail::TraceBuffer* _tbp = ::cppde::ndf_detail::trace_sink();
      if (_tbp) {
      auto& tb = *_tbp;
      const double dt_d = controller_detail::scalar_value(dt);
      const double t_d  = controller_detail::scalar_value(t);
      tb.nst.push_back(m_n_accepted + m_n_rejected);
      tb.t.push_back(t_d + dt_d);
      tb.h.push_back(dt_d);
      tb.q.push_back(static_cast<int>(stepper_type::stepper_order));
      tb.dsm.push_back(err);
      tb.acnrm.push_back(err);
      tb.acnrm_state.push_back(err);
      tb.tq2.push_back(1.0);
      tb.gamma.push_back(0.0);
      tb.gamrat.push_back(1.0);
      tb.newton_conv.push_back(1);
      tb.mode.emplace_back("ONESTEP");
      tb.nfe.push_back(0);
      tb.njev.push_back(0);
      tb.nsetups.push_back(0);
      tb.setup_reason.emplace_back("");
      tb.pece_iters.push_back(0);
      tb.pece_diverged.push_back(0);
      }
    }
#endif

    return update_stepsize(err, t, dt);
  }

  // ====================================================================
  //  try_step: mutable input (for dense output compatibility)
  // ====================================================================

  template<class System>
  controlled_step_result
  try_step(System& sys, state_type& x_in, time_type& t, state_type& x_out, time_type& dt)
  {
    return try_step(sys, static_cast<const state_type&>(x_in), t, x_out, dt);
  }

  // ====================================================================
  //  Event reset
  // ====================================================================

  void reset_after_event(time_type /*dt_before*/)
  {
    m_first_step = true;
    m_err_old = 1.0;
    m_dt_old = 1.0;
    m_last_rejected = false;
    m_stepper.invalidate_lu();
  }

  // ====================================================================
  //  Accessors
  // ====================================================================

  stepper_type&       stepper()       { return m_stepper; }
  const stepper_type& stepper() const { return m_stepper; }

  // Slab-prime forwarder. It has to run before the move into the dense-output
  // wrapper, so that the wrapper inherits a primed slab. The controller's own
  // buffers may still be empty here; the resize helpers prime again once they
  // grow.
  void prepare_sensitivities(unsigned n_sens)
  {
    m_n_sens = n_sens;
    m_stepper.prepare_sensitivities(n_sens);
    if constexpr (detail::is_dynamic_dual<value_type>::value) {
      if (n_sens == 0) return;
      if (!m_xerr.m_v.empty())
        m_xerr_slab.prime(m_xerr.m_v,
                          static_cast<unsigned>(m_xerr.m_v.size()), n_sens);
      if (!m_xnew.m_v.empty())
        m_xnew_slab.prime(m_xnew.m_v,
                          static_cast<unsigned>(m_xnew.m_v.size()), n_sens);
    }
  }

  template<class StateType>
  void adjust_size(const StateType& x)
  {
    resize_m_xerr(x);
    resize_m_xnew(x);
  }

  double last_error() const { return m_err_old; }
  bool first_step() const { return m_first_step; }
  bool last_rejected() const { return m_last_rejected; }
  double dt_old() const { return m_dt_old; }

  // --- Diagnostics ---
  int n_accepted() const { return m_n_accepted; }
  int n_rejected() const { return m_n_rejected; }
  int n_fevals()   const { return m_stepper.n_fevals(); }
  int n_jevals()   const { return m_stepper.n_jevals(); }
  int n_setups()   const { return m_stepper.n_setups(); }
  static constexpr int current_method_order() { return static_cast<int>(stepper_type::stepper_order); }
  void reset_counters() { m_n_accepted = 0; m_n_rejected = 0; m_stepper.reset_counters(); }

  // Tolerance accessors
  double atol() const { return m_atol; }
  double rtol() const { return m_rtol; }
  void set_tolerances(double atol, double rtol) { m_atol = atol; m_rtol = rtol; }

  // With the tangents out of the error test the step sequence reads value
  // arithmetic only, so it is the one a scalar run takes whatever the tangent
  // count is. The goal-oriented term below is unaffected; it is a double.
  void set_sens_err_con(bool v) { m_sens_err_con = v; }
  bool sens_err_con() const { return m_sens_err_con; }

  // PI controller parameter accessors
  double alpha() const { return m_alpha; }
  double beta() const { return m_beta; }
  double safety() const { return m_safety; }
  double max_factor() const { return m_max_factor; }
  double min_factor() const { return m_min_factor; }
  void set_pi_gains(double alpha, double beta) { m_alpha = alpha; m_beta = beta; }
  void set_pi_parameters(double alpha, double beta, double safety,
                         double max_factor, double min_factor)
  {
    m_alpha = alpha; m_beta = beta; m_safety = safety;
    m_max_factor = max_factor; m_min_factor = min_factor;
  }

private:

  // ====================================================================
  //  PI step-size control
  // ====================================================================

  controlled_step_result update_stepsize(double err, time_type& t, time_type& dt)
  {
    if (err <= 1.0) {
      // === Step accepted ===
      const double factor = onestep_detail::accept_factor(
          err, m_err_old, !(m_first_step || m_last_rejected), m_last_rejected,
          order, m_alpha, m_beta, m_safety, m_min_factor, m_max_factor);

      // Update controller state
      m_dt_old = controller_detail::scalar_value(dt);
      m_err_old = std::max(0.01, err);
      m_first_step = false;
      m_last_rejected = false;

      // Advance time and update step size
      t += dt;
      dt *= factor;

      // Apply max_dt limit
      if (m_max_dt != 0.0) {
        double dt_val = controller_detail::scalar_value(dt);
        if (std::abs(dt_val) > m_max_dt) {
          dt = time_type(std::copysign(m_max_dt, dt_val));
        }
      }

      ++m_n_accepted;
      return success;
    }
    else {
      // === Step rejected ===
      ++m_n_rejected;
      const double factor = onestep_detail::reject_factor(
          err, order, m_safety, m_min_factor);

      m_last_rejected = true;
      dt *= factor;

      return fail;
    }
  }

  // ====================================================================
  //  Resize helpers
  // ====================================================================

  template<class StateIn>
  bool resize_m_xerr(const StateIn& x)
  {
    bool resized = adjust_size_by_resizeability(m_xerr, x);
    if (resized && m_n_sens != 0)
      m_xerr_slab.prime(m_xerr.m_v,
                        static_cast<unsigned>(m_xerr.m_v.size()), m_n_sens);
    return resized;
  }

  template<class StateIn>
  bool resize_m_xnew(const StateIn& x)
  {
    bool resized = adjust_size_by_resizeability(m_xnew, x);
    if (resized && m_n_sens != 0)
      m_xnew_slab.prime(m_xnew.m_v,
                        static_cast<unsigned>(m_xnew.m_v.size()), m_n_sens);
    return resized;
  }

  // ====================================================================
  //  Members
  // ====================================================================

  stepper_type       m_stepper;
  resizer_type       m_xerr_resizer;
  resizer_type       m_xnew_resizer;
  wrapped_state_type m_xerr;
  wrapped_state_type m_xnew;

  // SoA tangent slabs for the dynamic-dual heap path (empty stubs otherwise).
  detail::tangent_slab<value_type> m_xerr_slab;
  detail::tangent_slab<value_type> m_xnew_slab;
  unsigned m_n_sens = 0;

  double m_atol, m_rtol;
  bool m_sens_err_con = true;   // tangents count in the error test
  double m_max_dt;

  // PI controller parameters
  double m_alpha, m_beta;
  double m_safety;
  double m_max_factor, m_min_factor;

  // Controller state
  bool   m_first_step;
  double m_err_old, m_dt_old;
  bool   m_last_rejected;

  // Diagnostics
  int m_n_accepted;
  int m_n_rejected;

public:
  mutable cppde::profiler m_prof;
  void finalize_profiler() const { m_prof.merge(m_stepper.m_prof); }
  void report_profiler(const char* label = nullptr) const {
    finalize_profiler();
    m_prof.report(label ? label : "cppDE one-step");
  }
};

} // namespace cppde

#endif // CPPDE_ONESTEP_CONTROLLER_HPP
