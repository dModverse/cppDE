/*
 R-visible half of the BLAS fork guard.

 The guard is installed from R_init_cppDE() because a pthread_atfork() handler
 cannot be registered from R and has to be in place before anything forks.  The
 handler and the vendor probe live in inst/include/cppde/cppde_blas_threads.hpp.

 Copyright (C) 2026 Simon Beyer
 */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include <cppde/cppde_blas_threads.hpp>

namespace {
bool guard_installed = false;
}

// list(api, threads, guard) for forkGuard().  NA where nothing resolved, so the
// R side can tell "pinned" from "no lever found" and say so.
extern "C" SEXP cppde_blas_info(void) {
  const char* nm = cppde::detail::blas_name();
  const int   nt = cppde::detail::blas_threads();

  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  SEXP nms = PROTECT(Rf_allocVector(STRSXP, 3));
  SET_STRING_ELT(nms, 0, Rf_mkChar("api"));
  SET_STRING_ELT(nms, 1, Rf_mkChar("threads"));
  SET_STRING_ELT(nms, 2, Rf_mkChar("guard"));
  Rf_setAttrib(out, R_NamesSymbol, nms);

  SET_VECTOR_ELT(out, 0, Rf_mkString(nm != nullptr ? nm : "NA"));
  if (nm == nullptr) SET_STRING_ELT(VECTOR_ELT(out, 0), 0, NA_STRING);
  SET_VECTOR_ELT(out, 1, Rf_ScalarInteger(nt > 0 ? nt : NA_INTEGER));
  SET_VECTOR_ELT(out, 2, Rf_ScalarLogical(guard_installed));

  UNPROTECT(2);
  return out;
}

extern "C" SEXP cppde_ptc(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" int cppde_ptc_v1(void (*)(double*, double*, double*, double*, double*, double*, int*, int*,
                                     int*, int*, int*),
                            int, int, int, double*, const double*, int, const int*, const int*, int,
                            const double*, const double*, const double*, int, int*, char*, int);
extern "C" SEXP cppde_ptc_entry(void);
extern "C" SEXP cppde_forcing_values(SEXP, SEXP, SEXP);
extern "C" SEXP cppde_eval_prep(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP cppde_eval_wrap(SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP cppde_eval_many(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef callMethods[] = {
  {"cppde_blas_info", (DL_FUNC) &cppde_blas_info, 0},
  {"cppde_ptc",       (DL_FUNC) &cppde_ptc,       9},
  {"cppde_ptc_entry", (DL_FUNC) &cppde_ptc_entry, 0},
  {"cppde_forcing_values", (DL_FUNC) &cppde_forcing_values, 3},
  {"cppde_eval_prep",      (DL_FUNC) &cppde_eval_prep,      6},
  {"cppde_eval_wrap",      (DL_FUNC) &cppde_eval_wrap,      4},
  {"cppde_eval_many",      (DL_FUNC) &cppde_eval_many,      4},
  {NULL, NULL, 0}
};

extern "C" void R_init_cppDE(DllInfo* dll) {
  R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_RegisterCCallable("cppDE", "cppde_ptc_v1", (DL_FUNC) &cppde_ptc_v1);
  guard_installed = cppde::detail::install_blas_fork_guard();
}
