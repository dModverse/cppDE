// Values of PCHIP forcings at given times, the interpolant the generated
// models evaluate.

#include <R.h>
#include <Rinternals.h>

#include <cstdio>
#include <stdexcept>
#include <vector>

#include <cppde/cppde_pchip_forcing.hpp>

// times: numeric; ftimes, fvalues: lists of numeric, one entry per forcing.
// Returns an n_times x n_forcings matrix.
extern "C" SEXP cppde_forcing_values(SEXP times, SEXP ftimes, SEXP fvalues) {
  const int nt = Rf_length(times);
  const int nf = Rf_length(ftimes);
  const double* t = REAL(times);
  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, nt, nf));
  double* o = REAL(out);
  char msg[256] = "";
  for (int f = 0; f < nf && !msg[0]; ++f) {
    SEXP ti = VECTOR_ELT(ftimes, f), vi = VECTOR_ELT(fvalues, f);
    std::vector<double> ft(REAL(ti), REAL(ti) + Rf_length(ti));
    std::vector<double> fv(REAL(vi), REAL(vi) + Rf_length(vi));
    try {
      cppde::PchipForcing<double> F(ft, fv);
      for (int i = 0; i < nt; ++i) o[i + static_cast<R_xlen_t>(nt) * f] = F(t[i]);
    } catch (const std::exception& e) {
      std::snprintf(msg, sizeof(msg), "%s", e.what());
    }
  }
  UNPROTECT(1);
  if (msg[0]) Rf_error("%s", msg);
  return out;
}
