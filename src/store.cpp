/*
 Early release of the checkpoints of a reverse solve, which R would free only at a garbage
 collection: the model's finalizer sits in the protected slot of the store and runs at once.

 Copyright (C) 2026 Simon Beyer
 */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern "C" SEXP cppde_release_store(SEXP p) {
  if (TYPEOF(p) != EXTPTRSXP || R_ExternalPtrAddr(p) == nullptr ||
      R_ExternalPtrTag(p) != Rf_install("cppde_rev_state"))
    return Rf_ScalarLogical(FALSE);
  SEXP fn = R_ExternalPtrProtected(p);
  if (TYPEOF(fn) != EXTPTRSXP) return Rf_ScalarLogical(FALSE);
  DL_FUNC f = R_ExternalPtrAddrFn(fn);
  if (f == nullptr) return Rf_ScalarLogical(FALSE);
  reinterpret_cast<void (*)(SEXP)>(f)(p);
  return Rf_ScalarLogical(TRUE);
}
