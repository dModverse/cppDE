/*
 Main header for cppDE: ODE integration and sensitivity calculation
 using in-tree forward-mode dual numbers for automatic differentiation.

 This is the ODE surface, not the whole library. The batch entry point, the
 chain-rule kernels and the return codes are included by the generated sources
 that need them, so they are not reachable from here.

 The stepper architecture (Rosenbrock4, NDF/BDF) is derived from Boost.Odeint
 by Karsten Ahnert, Mario Mulansky, and Christoph Koke (2011–2015),
 distributed under the Boost Software License, Version 1.0.
 Substantially rewritten: LAPACK/KLU linear algebra, AD-aware LU
 decomposition, event handling, NDF/BDF support, PI step-size control.

 Copyright (C) 2026 Simon Beyer
 */

#ifndef CPPDE_HPP
#define CPPDE_HPP

// ============================================================================
//  Standard library
// ============================================================================
#include <algorithm>
#include <vector>
#include <cmath>
#include <functional>

// ============================================================================
//  cppde::dual / cppde::dual2nd (forward-mode AD backend)
// ============================================================================
#include <cppde/cppde_dual_math.hpp>
#include <cppde/cppde_dual_expr.hpp>
#include <cppde/cppde_dual2nd.hpp>
#include <cppde/cppde_dual2nd_math.hpp>
#include <cppde/cppde_dual2nd_expr.hpp>

// ============================================================================
//  cppDE types and infrastructure
// ============================================================================
#include <cppde/cppde_types.hpp>
#include <cppde/cppde_odeint_compat.hpp>

// ============================================================================
//  cppDE AD-aware LU solvers (IFT for dense and sparse)
// ============================================================================
#include <cppde/cppde_ad_lu.hpp>
#include <cppde/cppde_sparse_ad_lu.hpp>

// ============================================================================
//  cppDE unified LU iteration matrix solver
// ============================================================================
#include <cppde/cppde_lu.hpp>

// ============================================================================
//  cppDE Newton solver
// ============================================================================
#include <cppde/cppde_newton.hpp>

// ============================================================================
//  cppDE utilities
// ============================================================================
#include <cppde/cppde_utils.hpp>
#include <cppde/cppde_pchip_forcing.hpp>

// ============================================================================
//  cppDE stepper traits (multi-step vs single-step dispatch)
// ============================================================================
#include <cppde/cppde_stepper_traits.hpp>

// ============================================================================
//  cppDE single-step methods (Rosenbrock4, Tsit5)
//
//  Unified onestep_controller and onestep_dense_output work with any
//  single-step stepper.  The old rosenbrock4_controller / _dense_output
//  headers are thin wrappers that include these.
// ============================================================================
#include <cppde/cppde_rosenbrock4.hpp>
#include <cppde/cppde_tsit5.hpp>
#include <cppde/cppde_onestep_controller.hpp>
#include <cppde/cppde_onestep_dense_output.hpp>
#include <cppde/cppde_integrate_times.hpp>
#include <cppde/cppde_step_checker.hpp>

// ============================================================================
//  cppDE multistep family (BDF / Adams)
//
//  One multistepper class with a method selector. The default coefficients are
//  the Klopfenstein-Shampine NDF family; classical BDF and Adams-Moulton are
//  instantiations of the same class.
// ============================================================================
#include <cppde/cppde_multistepper.hpp>
#include <cppde/cppde_multistepper_controller.hpp>
#include <cppde/cppde_multistepper_dense_output.hpp>

#endif // CPPDE_HPP
