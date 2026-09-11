/*
 Rosenbrock4 stepper

 The stepper architecture is derived from Boost.Odeint's rosenbrock4
 by Karsten Ahnert, Mario Mulansky, and Christoph Koke (2011–2013),
 distributed under the Boost Software License, Version 1.0.  See
 inst/COPYRIGHTS for the full license text and original copyright
 notice.

 Copyright (C) 2026 Simon Beyer
 */

#ifndef CPPDE_ROSENBROCK4_HPP
#define CPPDE_ROSENBROCK4_HPP

#include <cstddef>
#include <type_traits>

#include <cppde/cppde_types.hpp>
#include <cppde/cppde_odeint_compat.hpp>
#include <cppde/cppde_lu.hpp>
#include <cppde/cppde_ad_lu.hpp>   // for ad_lu::scalar_type_t
#include <cppde/cppde_dual_slab.hpp>
#include <cppde/cppde_stage_matrix.hpp>
#include <cppde/cppde_profiler.hpp>

namespace cppde {

// ============================================================================
//  Rosenbrock4 coefficients (Hairer-Wanner)
// ============================================================================

template<class Value>
struct default_rosenbrock_coefficients
{
  typedef Value value_type;
  typedef unsigned short order_type;

  default_rosenbrock_coefficients( void )
    : gamma ( static_cast<value_type>( 0.25 ) ) ,
      d1 ( static_cast<value_type>( 0.25 ) ) ,
      d2 ( static_cast<value_type>( -0.1043 ) ) ,
      d3 ( static_cast<value_type>( 0.1035 ) ) ,
      d4 ( static_cast<value_type>( 0.3620000000000023e-01 ) ) ,
      c2 ( static_cast<value_type>( 0.386 ) ) ,
      c3 ( static_cast<value_type>( 0.21 ) ) ,
      c4 ( static_cast<value_type>( 0.63 ) ) ,
      c21 ( static_cast<value_type>( -0.5668800000000000e+01 ) ) ,
      a21 ( static_cast<value_type>( 0.1544000000000000e+01 ) ) ,
      c31 ( static_cast<value_type>( -0.2430093356833875e+01 ) ) ,
      c32 ( static_cast<value_type>( -0.2063599157091915e+00 ) ) ,
      a31 ( static_cast<value_type>( 0.9466785280815826e+00 ) ) ,
      a32 ( static_cast<value_type>( 0.2557011698983284e+00 ) ) ,
      c41 ( static_cast<value_type>( -0.1073529058151375e+00 ) ) ,
      c42 ( static_cast<value_type>( -0.9594562251023355e+01 ) ) ,
      c43 ( static_cast<value_type>( -0.2047028614809616e+02 ) ) ,
      a41 ( static_cast<value_type>( 0.3314825187068521e+01 ) ) ,
      a42 ( static_cast<value_type>( 0.2896124015972201e+01 ) ) ,
      a43 ( static_cast<value_type>( 0.9986419139977817e+00 ) ) ,
      c51 ( static_cast<value_type>( 0.7496443313967647e+01 ) ) ,
      c52 ( static_cast<value_type>( -0.1024680431464352e+02 ) ) ,
      c53 ( static_cast<value_type>( -0.3399990352819905e+02 ) ) ,
      c54 ( static_cast<value_type>(  0.1170890893206160e+02 ) ) ,
      a51 ( static_cast<value_type>( 0.1221224509226641e+01 ) ) ,
      a52 ( static_cast<value_type>( 0.6019134481288629e+01 ) ) ,
      a53 ( static_cast<value_type>( 0.1253708332932087e+02 ) ) ,
      a54 ( static_cast<value_type>( -0.6878860361058950e+00 ) ) ,
      c61 ( static_cast<value_type>( 0.8083246795921522e+01 ) ) ,
      c62 ( static_cast<value_type>( -0.7981132988064893e+01 ) ) ,
      c63 ( static_cast<value_type>( -0.3152159432874371e+02 ) ) ,
      c64 ( static_cast<value_type>( 0.1631930543123136e+02 ) ) ,
      c65 ( static_cast<value_type>( -0.6058818238834054e+01 ) ) ,
      d21 ( static_cast<value_type>( 0.1012623508344586e+02 ) ) ,
      d22 ( static_cast<value_type>( -0.7487995877610167e+01 ) ) ,
      d23 ( static_cast<value_type>( -0.3480091861555747e+02 ) ) ,
      d24 ( static_cast<value_type>( -0.7992771707568823e+01 ) ) ,
      d25 ( static_cast<value_type>( 0.1025137723295662e+01 ) ) ,
      d31 ( static_cast<value_type>( -0.6762803392801253e+00 ) ) ,
      d32 ( static_cast<value_type>( 0.6087714651680015e+01 ) ) ,
      d33 ( static_cast<value_type>( 0.1643084320892478e+02 ) ) ,
      d34 ( static_cast<value_type>( 0.2476722511418386e+02 ) ) ,
      d35 ( static_cast<value_type>( -0.6594389125716872e+01 ) )
  {}

  const value_type gamma;
  const value_type d1 , d2 , d3 , d4;
  const value_type c2 , c3 , c4;
  const value_type c21;
  const value_type a21;
  const value_type c31 , c32;
  const value_type a31 , a32;
  const value_type c41 , c42 , c43;
  const value_type a41 , a42 , a43;
  const value_type c51 , c52 , c53 , c54;
  const value_type a51 , a52 , a53 , a54;
  const value_type c61 , c62 , c63 , c64 , c65;
  const value_type d21 , d22 , d23 , d24 , d25;
  const value_type d31 , d32 , d33 , d34 , d35;

  static const order_type stepper_order = 4;
  static const order_type error_order = 3;
};

// ============================================================================
//  Tags (re-exported from cppde_lu.hpp for external use)
// ============================================================================

// dense_lu_tag and sparse_lu_tag are defined in cppde_lu.hpp

// ============================================================================
//  jacobian_hint
// ============================================================================

enum class jacobian_hint : unsigned char {
  recompute_all   = 0,
    reuse_jacobian  = 1,
    reuse_lu        = 2
};

// ============================================================================
//  rosenbrock4<Value, JacobianPattern, Coefficients, Resizer>
// ============================================================================

template<
  class Value,
  class JacobianPattern = dense_lu_tag,
  class Coefficients = default_rosenbrock_coefficients<Value>,
  class Resizer = initially_resizer
>
class rosenbrock4
{
public:

  typedef Value                              value_type;
  typedef std::vector<value_type>            state_type;
  typedef state_type                         deriv_type;
  typedef ad_lu::scalar_type_t<value_type>   time_type;
  typedef dense_matrix<value_type>           matrix_type;
  typedef csc_matrix<value_type>             sparse_matrix_type;
  typedef Resizer                            resizer_type;
  typedef Coefficients                       rosenbrock_coefficients;
  typedef stepper_tag                        stepper_category;
  typedef unsigned short                     order_type;

  typedef state_wrapper<state_type>          wrapped_state_type;
  typedef state_wrapper<deriv_type>          wrapped_deriv_type;
  typedef state_wrapper<matrix_type>         wrapped_matrix_type;

  typedef rosenbrock4<Value, JacobianPattern, Coefficients, Resizer> stepper_type;

  // Same method on another scalar type. The coefficients follow Value2, so a
  // hand-supplied set does not survive the rebind.
  template<class Value2> using rebind_value =
    rosenbrock4<Value2, JacobianPattern,
                default_rosenbrock_coefficients<Value2>, Resizer>;

  static constexpr bool is_sparse = is_sparse_tag<JacobianPattern>::value;
  static constexpr bool is_dense = !is_sparse;

  using lu_type = lu_W<Value, is_sparse>;

  const static order_type stepper_order = rosenbrock_coefficients::stepper_order;
  const static order_type error_order = rosenbrock_coefficients::error_order;

  rosenbrock4( void )
    : m_lu() ,
      m_resizer() , m_x_err_resizer() ,
      m_dxdt() ,
      m_g1() , m_g2() , m_g3() , m_g4() , m_g5() ,
      m_cont3() , m_cont4() , m_xtmp() , m_x_err() ,
      m_coef() ,
      m_n_fevals( 0 ) , m_n_jevals( 0 ) , m_n_lu_setups( 0 )
  { }

  // Same move/non-copy semantics as the multistepper: copying the slab
  // would build a fresh tangent block while the dual elements still
  // point at the original: UB. Move-assignment preserves data() of
  // both the slab storage and the std::vector<dual> elements.
  rosenbrock4(const rosenbrock4&)            = delete;
  rosenbrock4& operator=(const rosenbrock4&) = delete;
  rosenbrock4(rosenbrock4&&)                 = default;
  rosenbrock4& operator=(rosenbrock4&&)      = default;

  // ====================================================================
  //  prepare_sensitivities
  //
  //  Stores n_sens and primes the slabs this stepper owns. The controller primes
  //  its own xerr. No-op unless the value type is a heap-backed dual.
  // ====================================================================

  void prepare_sensitivities(unsigned n_sens)
  {
    m_n_sens = n_sens;
    if constexpr (detail::is_dynamic_dual<value_type>::value) {
      if (n_sens == 0) return;
      auto prime = [n_sens](auto& wrapped, auto& slab) {
        if (!wrapped.m_v.empty())
          slab.prime(wrapped.m_v,
                     static_cast<unsigned>(wrapped.m_v.size()), n_sens);
      };
      prime(m_dxdt,  m_dxdt_slab);

      // The five Rosenbrock stages g1..g5 live in one contiguous tangent
      // block via m_G. Per-stage slab views are accessed as m_G.slab(j)
      // for the existing vec_*_with_slab call sites.
      const std::size_t n_g = m_g1.m_v.size();
      bool gs_ready = (n_g > 0)
                   && (m_g2.m_v.size() == n_g) && (m_g3.m_v.size() == n_g)
                   && (m_g4.m_v.size() == n_g) && (m_g5.m_v.size() == n_g);
      if (gs_ready) {
        std::array<std::vector<value_type>*, 5> facades{
          &m_g1.m_v, &m_g2.m_v, &m_g3.m_v, &m_g4.m_v, &m_g5.m_v
        };
        m_G.prime(facades, static_cast<unsigned>(n_g), n_sens);
      }

      prime(m_cont3, m_cont3_slab);
      prime(m_cont4, m_cont4_slab);
      prime(m_xtmp,  m_xtmp_slab);
      prime(m_x_err, m_x_err_slab);
    }
  }

  unsigned n_sens() const noexcept { return m_n_sens; }

  order_type order() const { return stepper_order; }

  int n_fevals() const { return m_n_fevals; }
  int n_jevals() const { return m_n_jevals; }
  int n_setups() const { return m_n_lu_setups; }
  void reset_counters() { m_n_fevals = 0; m_n_jevals = 0; m_n_lu_setups = 0; }

  // ====================================================================
  //  stages: the six Rosenbrock stages, with the linear solve left open.
  //
  //  The solve is a parameter so the stage arithmetic is stated once. dfdt is
  //  what the Jacobian evaluation filled; the caller owns it.
  // ====================================================================

  template<class DerivFunc, class Solve>
  void stages(DerivFunc& deriv_func, const state_type& x,
              time_type t_s, time_type dt_s,
              state_type& xout, state_type& xerr,
              state_type& dfdt, Solve& solve)
  {
    const size_t n = x.size();

    // --- Stage 1 ---
    vec_copy_with_slab(m_g1.m_v, m_G.slab(0), m_dxdt.m_v, m_dxdt_slab);
    vec_axpy_with_slab(m_g1.m_v, m_G.slab(0),
                       dt_s * ad_lu::scalar_value(m_coef.d1),
                       dfdt, m_dfdt_unslabbed);
    { auto _tp = m_prof.timer(prof_cat::lu_solve);
      solve(m_g1.m_v); }

    // --- Stage 2 ---
    vec_copy_with_slab(m_xtmp.m_v, m_xtmp_slab, x, m_x_in_unslabbed);
    vec_axpy_with_slab(m_xtmp.m_v, m_xtmp_slab,
                       ad_lu::scalar_value(m_coef.a21),
                       m_g1.m_v, m_G.slab(0));
    { auto _tp = m_prof.timer(prof_cat::f_eval);
      deriv_func(m_xtmp.m_v, m_g2.m_v, value_type(t_s + m_coef.c2 * dt_s)); }
    ++m_n_fevals;
    vec_axpy_with_slab(m_g2.m_v, m_G.slab(1),
                       dt_s * ad_lu::scalar_value(m_coef.d2),
                       dfdt, m_dfdt_unslabbed);
    vec_axpy_with_slab(m_g2.m_v, m_G.slab(1),
                       ad_lu::scalar_value(m_coef.c21) / dt_s,
                       m_g1.m_v, m_G.slab(0));
    { auto _tp = m_prof.timer(prof_cat::lu_solve);
      solve(m_g2.m_v); }

    // --- Stage 3 ---
    vec_copy_with_slab(m_xtmp.m_v, m_xtmp_slab, x, m_x_in_unslabbed);
    vec_axpy_with_slab(m_xtmp.m_v, m_xtmp_slab,
                       ad_lu::scalar_value(m_coef.a31),
                       m_g1.m_v, m_G.slab(0));
    vec_axpy_with_slab(m_xtmp.m_v, m_xtmp_slab,
                       ad_lu::scalar_value(m_coef.a32),
                       m_g2.m_v, m_G.slab(1));
    { auto _tp = m_prof.timer(prof_cat::f_eval);
      deriv_func(m_xtmp.m_v, m_g3.m_v, value_type(t_s + m_coef.c3 * dt_s)); }
    ++m_n_fevals;
    vec_axpy_with_slab(m_g3.m_v, m_G.slab(2),
                       dt_s * ad_lu::scalar_value(m_coef.d3),
                       dfdt, m_dfdt_unslabbed);
    vec_axpy_with_slab(m_g3.m_v, m_G.slab(2),
                       ad_lu::scalar_value(m_coef.c31) / dt_s,
                       m_g1.m_v, m_G.slab(0));
    vec_axpy_with_slab(m_g3.m_v, m_G.slab(2),
                       ad_lu::scalar_value(m_coef.c32) / dt_s,
                       m_g2.m_v, m_G.slab(1));
    { auto _tp = m_prof.timer(prof_cat::lu_solve);
      solve(m_g3.m_v); }

    // --- Stage 4 ---
    vec_copy_with_slab(m_xtmp.m_v, m_xtmp_slab, x, m_x_in_unslabbed);
    vec_axpy_with_slab(m_xtmp.m_v, m_xtmp_slab,
                       ad_lu::scalar_value(m_coef.a41),
                       m_g1.m_v, m_G.slab(0));
    vec_axpy_with_slab(m_xtmp.m_v, m_xtmp_slab,
                       ad_lu::scalar_value(m_coef.a42),
                       m_g2.m_v, m_G.slab(1));
    vec_axpy_with_slab(m_xtmp.m_v, m_xtmp_slab,
                       ad_lu::scalar_value(m_coef.a43),
                       m_g3.m_v, m_G.slab(2));
    { auto _tp = m_prof.timer(prof_cat::f_eval);
      deriv_func(m_xtmp.m_v, m_g4.m_v, value_type(t_s + m_coef.c4 * dt_s)); }
    ++m_n_fevals;
    vec_axpy_with_slab(m_g4.m_v, m_G.slab(3),
                       dt_s * ad_lu::scalar_value(m_coef.d4),
                       dfdt, m_dfdt_unslabbed);
    vec_axpy_with_slab(m_g4.m_v, m_G.slab(3),
                       ad_lu::scalar_value(m_coef.c41) / dt_s,
                       m_g1.m_v, m_G.slab(0));
    vec_axpy_with_slab(m_g4.m_v, m_G.slab(3),
                       ad_lu::scalar_value(m_coef.c42) / dt_s,
                       m_g2.m_v, m_G.slab(1));
    vec_axpy_with_slab(m_g4.m_v, m_G.slab(3),
                       ad_lu::scalar_value(m_coef.c43) / dt_s,
                       m_g3.m_v, m_G.slab(2));
    { auto _tp = m_prof.timer(prof_cat::lu_solve);
      solve(m_g4.m_v); }

    // --- Stage 5 ---
    vec_copy_with_slab(m_xtmp.m_v, m_xtmp_slab, x, m_x_in_unslabbed);
    vec_axpy_with_slab(m_xtmp.m_v, m_xtmp_slab,
                       ad_lu::scalar_value(m_coef.a51),
                       m_g1.m_v, m_G.slab(0));
    vec_axpy_with_slab(m_xtmp.m_v, m_xtmp_slab,
                       ad_lu::scalar_value(m_coef.a52),
                       m_g2.m_v, m_G.slab(1));
    vec_axpy_with_slab(m_xtmp.m_v, m_xtmp_slab,
                       ad_lu::scalar_value(m_coef.a53),
                       m_g3.m_v, m_G.slab(2));
    vec_axpy_with_slab(m_xtmp.m_v, m_xtmp_slab,
                       ad_lu::scalar_value(m_coef.a54),
                       m_g4.m_v, m_G.slab(3));
    { auto _tp = m_prof.timer(prof_cat::f_eval);
      deriv_func(m_xtmp.m_v, m_g5.m_v, value_type(t_s + dt_s)); }
    ++m_n_fevals;
    vec_axpy_with_slab(m_g5.m_v, m_G.slab(4),
                       ad_lu::scalar_value(m_coef.c51) / dt_s,
                       m_g1.m_v, m_G.slab(0));
    vec_axpy_with_slab(m_g5.m_v, m_G.slab(4),
                       ad_lu::scalar_value(m_coef.c52) / dt_s,
                       m_g2.m_v, m_G.slab(1));
    vec_axpy_with_slab(m_g5.m_v, m_G.slab(4),
                       ad_lu::scalar_value(m_coef.c53) / dt_s,
                       m_g3.m_v, m_G.slab(2));
    vec_axpy_with_slab(m_g5.m_v, m_G.slab(4),
                       ad_lu::scalar_value(m_coef.c54) / dt_s,
                       m_g4.m_v, m_G.slab(3));
    { auto _tp = m_prof.timer(prof_cat::lu_solve);
      solve(m_g5.m_v); }

    // --- Error estimate (stage 6) ---
    // Uses the Hairer-Wanner 6-stage formulation: an additional
    // f-evaluation and W⁻¹ solve to produce the embedded error.
    // The error is added to the solution (local extrapolation).
    vec_axpy_with_slab(m_xtmp.m_v, m_xtmp_slab,
                       1.0,
                       m_g5.m_v, m_G.slab(4));
    { auto _tp = m_prof.timer(prof_cat::f_eval);
      deriv_func(m_xtmp.m_v, xerr, value_type(t_s + dt_s)); }
    ++m_n_fevals;
  // xerr belongs to the controller and is not slab-bound. The alphas go in as
  // plain doubles, the coefficients having no tangents anyway, which keeps the
  // axpy free of arena allocations.
    vec_axpy(xerr, ad_lu::scalar_value(m_coef.c61) / dt_s, m_g1.m_v);
    vec_axpy(xerr, ad_lu::scalar_value(m_coef.c62) / dt_s, m_g2.m_v);
    vec_axpy(xerr, ad_lu::scalar_value(m_coef.c63) / dt_s, m_g3.m_v);
    vec_axpy(xerr, ad_lu::scalar_value(m_coef.c64) / dt_s, m_g4.m_v);
    vec_axpy(xerr, ad_lu::scalar_value(m_coef.c65) / dt_s, m_g5.m_v);
    { auto _tp = m_prof.timer(prof_cat::lu_solve);
      solve(xerr); }

    // --- Solution ---
    for (size_t i = 0; i < n; ++i)
      xout[i] = m_xtmp.m_v[i] + xerr[i];
  }

  // ====================================================================
  //  do_step (with error output)
  // ====================================================================

  template<class Sys, class TimeArg>
  void do_step(
      Sys& system,
      const state_type& x, TimeArg t,
      state_type& xout, TimeArg dt, state_type& xerr,
      jacobian_hint hint = jacobian_hint::recompute_all)
  {
    auto& deriv_func = system.first;
    auto& jacobi_func = system.second;

    const size_t n = x.size();

    // Extract scalar time values (strip AD derivatives from incoming time args)
    const time_type t_s = static_cast<time_type>(ad_lu::scalar_value(t));
    const time_type dt_s = static_cast<time_type>(ad_lu::scalar_value(dt));

    m_resizer.adjust_size(x, [this](auto&& arg) {
      return this->resize_impl(std::forward<decltype(arg)>(arg));
    });

    // --- Initial derivative ---
    { auto _tp = m_prof.timer(prof_cat::f_eval);
      deriv_func(x, m_dxdt.m_v, value_type(t_s)); }
    ++m_n_fevals;

    const value_type inv_gamma_dt =
    static_cast<value_type>(1) / (m_coef.gamma * dt_s);

    // --- Jacobian / LU ---
    if (hint == jacobian_hint::reuse_lu && m_lu.has_valid_lu())
    {
      // Reuse everything
    }
    else if (hint == jacobian_hint::reuse_jacobian && m_lu.has_valid_jacobian())
    {
      { auto _tp = m_prof.timer(prof_cat::lu_factor);
        m_lu.refactorize_W_from_cache(n, inv_gamma_dt); }
      ++m_n_lu_setups;
      m_lu.set_lu_valid(dt_s);
    }
    else
    {
      { auto _tp = m_prof.timer(prof_cat::jac_eval);
        m_lu.call_jacobian(jacobi_func, x, t); }
      ++m_n_jevals;
      m_lu.cache_jacobian(n);
      { auto _tp = m_prof.timer(prof_cat::lu_factor);
        m_lu.factorize_W(n, inv_gamma_dt); }
      ++m_n_lu_setups;
      m_lu.set_jacobian_valid();
      m_lu.set_lu_valid(dt_s);
    }
    auto lu_solve = [this](state_type& v) { m_lu.solve(v); };
    stages(deriv_func, x, t_s, dt_s, xout, xerr, m_lu.dfdt_mut(), lu_solve);
  }

  // ====================================================================
  //  Convenience overloads
  // ====================================================================

  template<class Sys, class TimeArg>
  void do_step(Sys& system, state_type& x, TimeArg t,
               TimeArg dt, state_type& xerr,
               jacobian_hint hint = jacobian_hint::recompute_all)
  { do_step(system, x, t, x, dt, xerr, hint); }

  template<class Sys, class TimeArg>
  void do_step(Sys& system, const state_type& x, TimeArg t,
               state_type& xout, TimeArg dt,
               jacobian_hint hint = jacobian_hint::recompute_all)
  {
    m_x_err_resizer.adjust_size(x, [this](auto&& arg) {
      return this->resize_x_err<state_type>(std::forward<decltype(arg)>(arg));
    });
    do_step(system, x, t, xout, dt, m_x_err.m_v, hint);
  }

  template<class Sys, class TimeArg>
  void do_step(Sys& system, state_type& x, TimeArg t, TimeArg dt,
               jacobian_hint hint = jacobian_hint::recompute_all)
  {
    m_x_err_resizer.adjust_size(x, [this](auto&& arg) {
      return this->resize_x_err<state_type>(std::forward<decltype(arg)>(arg));
    });
    do_step(system, x, t, dt, m_x_err.m_v, hint);
  }

  // ====================================================================
  //  Dense output
  // ====================================================================

  void prepare_dense_output()
  {
    auto _tp = m_prof.timer(prof_cat::dense_snapshot);
    // cont3 = d21*g1 + d22*g2 + d23*g3 + d24*g4 + d25*g5
    vec_zero_with_slab(m_cont3.m_v, m_cont3_slab);
    vec_axpy_with_slab(m_cont3.m_v, m_cont3_slab,
                       ad_lu::scalar_value(m_coef.d21), m_g1.m_v, m_G.slab(0));
    vec_axpy_with_slab(m_cont3.m_v, m_cont3_slab,
                       ad_lu::scalar_value(m_coef.d22), m_g2.m_v, m_G.slab(1));
    vec_axpy_with_slab(m_cont3.m_v, m_cont3_slab,
                       ad_lu::scalar_value(m_coef.d23), m_g3.m_v, m_G.slab(2));
    vec_axpy_with_slab(m_cont3.m_v, m_cont3_slab,
                       ad_lu::scalar_value(m_coef.d24), m_g4.m_v, m_G.slab(3));
    vec_axpy_with_slab(m_cont3.m_v, m_cont3_slab,
                       ad_lu::scalar_value(m_coef.d25), m_g5.m_v, m_G.slab(4));
    // cont4 = d31*g1 + d32*g2 + d33*g3 + d34*g4 + d35*g5
    vec_zero_with_slab(m_cont4.m_v, m_cont4_slab);
    vec_axpy_with_slab(m_cont4.m_v, m_cont4_slab,
                       ad_lu::scalar_value(m_coef.d31), m_g1.m_v, m_G.slab(0));
    vec_axpy_with_slab(m_cont4.m_v, m_cont4_slab,
                       ad_lu::scalar_value(m_coef.d32), m_g2.m_v, m_G.slab(1));
    vec_axpy_with_slab(m_cont4.m_v, m_cont4_slab,
                       ad_lu::scalar_value(m_coef.d33), m_g3.m_v, m_G.slab(2));
    vec_axpy_with_slab(m_cont4.m_v, m_cont4_slab,
                       ad_lu::scalar_value(m_coef.d34), m_g4.m_v, m_G.slab(3));
    vec_axpy_with_slab(m_cont4.m_v, m_cont4_slab,
                       ad_lu::scalar_value(m_coef.d35), m_g5.m_v, m_G.slab(4));
  }

  // ====================================================================
  //  What a written adjoint reads.
  //
  //  The stage recursion belongs to the method and is stated once, in stages().
  //  The adjoint transposes it, and reads the coefficients and the stage
  //  vectors from here rather than carrying a second copy of either.
  // ====================================================================

  /// Stage vectors g1..g5; the sixth solve is the embedded error, which
  /// do_step hands back separately.
  static constexpr int n_stages_used = 5;

  /// The explicit combination weights a(i, j), i > j, both one-based.
  static double stage_a(int i, int j) {
    const rosenbrock_coefficients c;
    switch (i * 10 + j) {
      case 21: return ad_lu::scalar_value(c.a21);
      case 31: return ad_lu::scalar_value(c.a31);
      case 32: return ad_lu::scalar_value(c.a32);
      case 41: return ad_lu::scalar_value(c.a41);
      case 42: return ad_lu::scalar_value(c.a42);
      case 43: return ad_lu::scalar_value(c.a43);
      case 51: return ad_lu::scalar_value(c.a51);
      case 52: return ad_lu::scalar_value(c.a52);
      case 53: return ad_lu::scalar_value(c.a53);
      case 54: return ad_lu::scalar_value(c.a54);
      default: return 0.0;
    }
  }

  /// The stage-coupling weights c(i, j), i > j, both one-based. Row six
  /// belongs to the error estimate. The recursion divides them by the step.
  static double stage_c(int i, int j) {
    const rosenbrock_coefficients c;
    switch (i * 10 + j) {
      case 21: return ad_lu::scalar_value(c.c21);
      case 31: return ad_lu::scalar_value(c.c31);
      case 32: return ad_lu::scalar_value(c.c32);
      case 41: return ad_lu::scalar_value(c.c41);
      case 42: return ad_lu::scalar_value(c.c42);
      case 43: return ad_lu::scalar_value(c.c43);
      case 51: return ad_lu::scalar_value(c.c51);
      case 52: return ad_lu::scalar_value(c.c52);
      case 53: return ad_lu::scalar_value(c.c53);
      case 54: return ad_lu::scalar_value(c.c54);
      case 61: return ad_lu::scalar_value(c.c61);
      case 62: return ad_lu::scalar_value(c.c62);
      case 63: return ad_lu::scalar_value(c.c63);
      case 64: return ad_lu::scalar_value(c.c64);
      case 65: return ad_lu::scalar_value(c.c65);
      default: return 0.0;
    }
  }

  /// The df/dt weights d_i, one-based; stages five and six carry none.
  static double stage_d(int i) {
    const rosenbrock_coefficients c;
    switch (i) {
      case 1: return ad_lu::scalar_value(c.d1);
      case 2: return ad_lu::scalar_value(c.d2);
      case 3: return ad_lu::scalar_value(c.d3);
      case 4: return ad_lu::scalar_value(c.d4);
      default: return 0.0;
    }
  }

  /// The nodes, one-based; the first is zero and the last two are one.
  static double stage_node(int i) {
    const rosenbrock_coefficients c;
    switch (i) {
      case 2: return ad_lu::scalar_value(c.c2);
      case 3: return ad_lu::scalar_value(c.c3);
      case 4: return ad_lu::scalar_value(c.c4);
      case 5: case 6: return 1.0;
      default: return 0.0;
    }
  }

  /// The two continuous-extension vectors as combinations of the stages:
  /// which = 3 gives cont3's weights, which = 4 gives cont4's, j one-based.
  static double dense_stage_weight(int which, int j) {
    const rosenbrock_coefficients c;
    if (which == 3) switch (j) {
      case 1: return ad_lu::scalar_value(c.d21);
      case 2: return ad_lu::scalar_value(c.d22);
      case 3: return ad_lu::scalar_value(c.d23);
      case 4: return ad_lu::scalar_value(c.d24);
      case 5: return ad_lu::scalar_value(c.d25);
      default: return 0.0;
    }
    switch (j) {
      case 1: return ad_lu::scalar_value(c.d31);
      case 2: return ad_lu::scalar_value(c.d32);
      case 3: return ad_lu::scalar_value(c.d33);
      case 4: return ad_lu::scalar_value(c.d34);
      case 5: return ad_lu::scalar_value(c.d35);
      default: return 0.0;
    }
  }

  /// 1 / (gamma h), the diagonal W is built with.
  static double inv_gamma_dt_of(double dt_s) {
    const rosenbrock_coefficients c;
    return 1.0 / (ad_lu::scalar_value(c.gamma) * dt_s);
  }

  /// Stage vector i, one-based, valid after do_step.
  const state_type& stage_g(int i) const {
    switch (i) {
      case 1: return m_g1.m_v; case 2: return m_g2.m_v; case 3: return m_g3.m_v;
      case 4: return m_g4.m_v; default: return m_g5.m_v;
    }
  }

  /// W^-T b on the factorisation the step itself used. The stage solves are
  /// direct, so the adjoint has to transpose that matrix and no other.
  void stage_solve_transposed(state_type& b) { m_lu.solve_transposed(b); }

  /// f(x, t) at the step start, valid after do_step.
  const state_type& stage_f0() const { return m_dxdt.m_v; }

  /// df/dt at the step start, as the Jacobian evaluation filled it.
  const state_type& stage_dfdt() const { return m_lu.dfdt(); }

  /// The continuous extension's four weights at theta, over
  /// (x_old, x_new, cont3, cont4). One statement of the basis, two readers:
  /// calc_state contracts it forward, the written adjoint transposes it.
  template<class TimeArg>
  static void dense_weights(const TimeArg& s, TimeArg* w) {
    const TimeArg s1 = TimeArg(1.0) - s;
    w[0] = s1;
    w[1] = s;
    w[2] = s * s1;
    w[3] = s * s1 * s;
  }

  template<class TimeArg>
  void calc_state(TimeArg t, state_type& x,
                  const state_type& x_old, TimeArg t_old,
                  const state_type& x_new, TimeArg t_new)
  {
    auto _tp = m_prof.timer(prof_cat::dense_interp);
    const size_t n = m_g1.m_v.size();
    TimeArg dt = t_new - t_old;
    TimeArg s  = (t - t_old) / dt;
    TimeArg s1 = 1.0 - s;
    for (size_t i = 0; i < n; ++i)
      x[i] = x_old[i] * s1 + s * (x_new[i] + s1 * (m_cont3.m_v[i] + s * m_cont4.m_v[i]));
  }

  template<class StateType>
  void adjust_size(const StateType& x)
  { resize_impl(x); resize_x_err(x); }

  // LU access (for controller)
  void invalidate_lu() { m_lu.invalidate(); }
  bool has_valid_jacobian() const { return m_lu.has_valid_jacobian(); }
  bool has_valid_lu() const { return m_lu.has_valid_lu(); }
  double last_factorized_dt() const {
    using ad_lu::scalar_value;
    return static_cast<double>(scalar_value(m_lu.last_factorized_dt()));
  }

  lu_type& lu() { return m_lu; }
  const lu_type& lu() const { return m_lu; }

protected:

  template<class StateIn>
  bool resize_impl(const StateIn& x)
  {
    bool resized = false;
    resized |= adjust_size_by_resizeability(m_dxdt,    x);
    resized |= adjust_size_by_resizeability(m_xtmp,    x);
    resized |= adjust_size_by_resizeability(m_g1,      x);
    resized |= adjust_size_by_resizeability(m_g2,      x);
    resized |= adjust_size_by_resizeability(m_g3,      x);
    resized |= adjust_size_by_resizeability(m_g4,      x);
    resized |= adjust_size_by_resizeability(m_g5,      x);
    resized |= adjust_size_by_resizeability(m_cont3,   x);
    resized |= adjust_size_by_resizeability(m_cont4,   x);
    resized |= m_lu.resize(x);
    if (resized) m_lu.invalidate();
    if (resized && m_n_sens != 0) prepare_sensitivities(m_n_sens);
    return resized;
  }

  template<class StateIn>
  bool resize_x_err(const StateIn& x)
  {
    bool resized = adjust_size_by_resizeability(m_x_err, x);
    if (resized && m_n_sens != 0) prepare_sensitivities(m_n_sens);
    return resized;
  }

private:

  lu_type m_lu;

  resizer_type m_resizer;
  resizer_type m_x_err_resizer;

  wrapped_deriv_type   m_dxdt;
  wrapped_state_type   m_g1, m_g2, m_g3, m_g4, m_g5;
  wrapped_state_type   m_cont3, m_cont4;
  wrapped_state_type   m_xtmp;
  wrapped_state_type   m_x_err;

  // SoA tangent storage for the dynamic-dual heap path. Empty stubs
  // for non-dual value_type. The 5 Rosenbrock stages g1..g5 share one
  // contiguous tangent block via m_G; per-stage slab views are
  // m_G.slab(j). Other vectors keep individual slabs.
  detail::tangent_slab<value_type> m_dxdt_slab;
  detail::stage_matrix<value_type, 5> m_G;
  detail::tangent_slab<value_type> m_cont3_slab, m_cont4_slab;
  detail::tangent_slab<value_type> m_xtmp_slab;
  detail::tangent_slab<value_type> m_x_err_slab;
  // Placeholder slabs for the vectors living outside this stepper: the
  // controller's xerr, lu_W's m_dfdt and the input x. Permanently unprimed, so
  // the slab-aware helpers fall through to their per-element path. Mutable only
  // to satisfy the helper signatures.
  mutable detail::tangent_slab<value_type> m_dfdt_unslabbed;
  mutable detail::tangent_slab<value_type> m_x_in_unslabbed;
  unsigned m_n_sens = 0;

  const rosenbrock_coefficients m_coef;

  int m_n_fevals;
  int m_n_jevals;
  int m_n_lu_setups;

public:
  mutable cppde::profiler m_prof;
};

} // namespace cppde

#endif // CPPDE_ROSENBROCK4_HPP
