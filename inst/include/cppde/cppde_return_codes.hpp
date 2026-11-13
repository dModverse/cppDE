/*
 cppDE return-code constants.

 Values are kept numerically identical to the SUNDIALS CVODE flag scheme
 (see sundials/cvode.h) so that solveODE()'s `return_code` can be
 interpreted uniformly across the CVODE backend and the native cppDE
 backend.  The CVODE backend simply passes the raw `CVode()` flag through;
 the native backend maps the exceptions it throws to the closest CVODE flag.

 Copyright (C) 2026 Simon Beyer
 */

#ifndef CPPDE_RETURN_CODES_HPP
#define CPPDE_RETURN_CODES_HPP

namespace cppde {

// --- Success ---
constexpr int RC_SUCCESS            =   0;  // CV_SUCCESS

// --- Error codes (CVODE negative codes) ---
constexpr int RC_TOO_MUCH_WORK      =  -1;  // CV_TOO_MUCH_WORK   (maxsteps)
constexpr int RC_CONV_FAILURE       =  -4;  // CV_CONV_FAILURE    (maxprogress / repeated conv failures)
constexpr int RC_LINIT_FAIL         =  -5;  // CV_LINIT_FAIL
constexpr int RC_LSETUP_FAIL        =  -6;  // CV_LSETUP_FAIL     (KLU factor failure in native backend)
constexpr int RC_ILL_INPUT          = -22;  // CV_ILL_INPUT
constexpr int RC_NO_MALLOC          = -23;  // CV_NO_MALLOC
constexpr int RC_UNRECOGNIZED_ERR   = -99;  // CV_UNRECOGNIZED_ERR (fallback for C++ exceptions)

} // namespace cppde

#endif // CPPDE_RETURN_CODES_HPP
