// Argument sets of the AD batch entries of cppFUN() objects aligned to the function's internal
// order, and their results named, with no R call per set. Sets without seeds, with attached inputs,
// or with inputs other than double matrices and named vectors take the R path of .evaluateBatch_impl().

// no length(), error() ... macros: they break libc++ headers included below
#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

#include <algorithm>
#include <cstring>
#include <limits>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif
#ifndef _WIN32
#include <pthread.h>
#endif

namespace {

SEXP elt(SEXP list, const char* name) {
  SEXP nm = Rf_getAttrib(list, R_NamesSymbol);
  if (Rf_isNull(nm)) return R_NilValue;
  for (R_xlen_t i = 0; i < Rf_xlength(list); ++i)
    if (std::strcmp(CHAR(STRING_ELT(nm, i)), name) == 0) return VECTOR_ELT(list, i);
  return R_NilValue;
}

bool same(SEXP a, SEXP b) { return a == b || std::strcmp(CHAR(a), CHAR(b)) == 0; }

// The position of every entry of `want` in `have`, -1 where absent.
std::vector<int> positions(SEXP want, SEXP have) {
  const int nw = Rf_isNull(want) ? 0 : Rf_length(want);
  const int nh = Rf_isNull(have) ? 0 : Rf_length(have);
  std::vector<int> out(nw, -1);
  for (int i = 0; i < nw; ++i) {
    SEXP w = STRING_ELT(want, i);
    if (i < nh && same(w, STRING_ELT(have, i))) { out[i] = i; continue; }
    for (int j = 0; j < nh; ++j)
      if (same(w, STRING_ELT(have, j))) { out[i] = j; break; }
  }
  return out;
}

SEXP dimnamesAt(SEXP x, int k) {
  SEXP dn = Rf_getAttrib(x, R_DimNamesSymbol);
  if (Rf_isNull(dn) || Rf_length(dn) <= k) return R_NilValue;
  return VECTOR_ELT(dn, k);
}

int dimAt(SEXP x, int k) {
  SEXP d = Rf_getAttrib(x, R_DimSymbol);
  return (Rf_isNull(d) || Rf_length(d) <= k) ? -1 : INTEGER(d)[k];
}

bool isDouble(SEXP x, int ndim) {
  return !Rf_isNull(x) && TYPEOF(x) == REALSXP && dimAt(x, ndim - 1) >= 0 &&
         Rf_length(Rf_getAttrib(x, R_DimSymbol)) == ndim;
}

// libgomp's thread pool does not survive fork(), so a forked child evaluates serially, as the
// batch entries of the models do.
bool childOfFork = false;
[[maybe_unused]] const bool forkGuard = [] {
#ifndef _WIN32
  ::pthread_atfork(nullptr, nullptr, +[] { childOfFork = true; });
#endif
  return true;
}();

int threadsFor(int K, int requested) {
#ifdef _OPENMP
  // the requested team even for fewer sets, as batch_threads() of the solves
  if (childOfFork || omp_in_parallel() || K <= 1) return 1;
  return std::max(1, std::min(requested, omp_get_num_procs()));
#else
  (void) K; (void) requested;
  return 1;
#endif
}

}  // namespace

// innames, parameters, fixed: of the function; sets: the argument sets of evaluateBatch(); n_out
// the output count. With a list of parameters, the four are lists and a vector, one entry per set, for
// sets of different functions. Returns list(sets, theta), the sets as the AD batch entry reads
// them and the direction names per set, or NULL where a set needs the R path.
extern "C" SEXP cppde_eval_prep(SEXP inn, SEXP par, SEXP fix, SEXP sets, SEXP deriv2SEXP, SEXP nOutSEXP) {
  const bool deriv2 = Rf_asLogical(deriv2SEXP) == TRUE;
  const bool perSet = TYPEOF(par) == VECSXP;
  const int K = Rf_length(sets);
  const SEXP derivSym = Rf_install("deriv"), deriv2Sym = Rf_install("deriv2");

  // every object goes into `out` as soon as it is allocated, so `out` alone is protected
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
  SEXP outSets = Rf_allocVector(VECSXP, K);
  SET_VECTOR_ELT(out, 0, outSets);
  SEXP thetas = Rf_allocVector(VECSXP, K);
  SET_VECTOR_ELT(out, 1, thetas);
  auto none = [&]() { UNPROTECT(1); return R_NilValue; };
  auto valid = [](const std::vector<int>& v) { for (int x : v) if (x < 0) return false; return true; };

  for (int k = 0; k < K; ++k) {
    SEXP innames = perSet ? VECTOR_ELT(inn, k) : inn;
    SEXP parameters = perSet ? VECTOR_ELT(par, k) : par;
    SEXP stFixed = perSet ? VECTOR_ELT(fix, k) : fix;
    const int n_vars = Rf_isNull(innames) ? 0 : Rf_length(innames);
    const int n_params = Rf_isNull(parameters) ? 0 : Rf_length(parameters);
    const int n_out = perSet ? INTEGER(nOutSEXP)[k] : Rf_asInteger(nOutSEXP);
    SEXP a = VECTOR_ELT(sets, k);
    if (Rf_asLogical(elt(a, "attach.input")) == TRUE) return none();
    SEXP vars = elt(a, "vars"), params = elt(a, "params");
    SEXP tX = elt(a, "tangentX"), tP = elt(a, "tangentP");
    if (Rf_isNull(tX) && !Rf_isNull(vars)) tX = Rf_getAttrib(vars, derivSym);
    if (Rf_isNull(tP) && !Rf_isNull(params)) tP = Rf_getAttrib(params, derivSym);
    SEXP hX = R_NilValue, hP = R_NilValue;
    if (deriv2) {
      hX = elt(a, "hessianX"); hP = elt(a, "hessianP");
      if (Rf_isNull(hX) && !Rf_isNull(vars)) hX = Rf_getAttrib(vars, deriv2Sym);
      if (Rf_isNull(hP) && !Rf_isNull(params)) hP = Rf_getAttrib(params, deriv2Sym);
    }
    if (Rf_isNull(tX) && Rf_isNull(tP) && Rf_isNull(hX) && Rf_isNull(hP)) return none();
    if ((!Rf_isNull(tX) && !isDouble(tX, 3)) || (!Rf_isNull(tP) && !isDouble(tP, 2)) ||
        (!Rf_isNull(hX) && !isDouble(hX, 4)) || (!Rf_isNull(hP) && !isDouble(hP, 3)))
      return none();

    // the inputs in the order of the function
    int n_obs = 1;
    std::vector<int> vi;
    if (n_vars > 0) {
      if (!isDouble(vars, 2) || Rf_isNull(dimnamesAt(vars, 1))) return none();
      n_obs = dimAt(vars, 0);
      vi = positions(innames, dimnamesAt(vars, 1));
      if (!valid(vi)) return none();
    } else if (!Rf_isNull(vars)) {
      if (!isDouble(vars, 2)) return none();
      n_obs = dimAt(vars, 0) > 0 ? dimAt(vars, 0) : 1;
    }
    std::vector<int> pi;
    if (n_params > 0) {
      if (Rf_isNull(params) || TYPEOF(params) != REALSXP) return none();
      pi = positions(parameters, Rf_getAttrib(params, R_NamesSymbol));
      if (!valid(pi)) return none();
    }

    // the directions: the first seed that names them, the others name the same in any order
    SEXP cands[4] = {Rf_isNull(tX) ? R_NilValue : dimnamesAt(tX, 2), Rf_isNull(tP) ? R_NilValue : dimnamesAt(tP, 1),
                     Rf_isNull(hX) ? R_NilValue : dimnamesAt(hX, 2), Rf_isNull(hP) ? R_NilValue : dimnamesAt(hP, 1)};
    SEXP theta = R_NilValue;
    for (SEXP c : cands) if (!Rf_isNull(c) && Rf_length(c) > 0) { theta = c; break; }
    const int nt = Rf_isNull(theta) ? 0 : Rf_length(theta);
    for (SEXP c : cands)
      if (!Rf_isNull(c) && Rf_length(c) > 0 && c != theta && (Rf_length(c) != nt || !valid(positions(theta, c))))
        return none();
    // a seed with directions has them named
    if (nt > 0 && ((!Rf_isNull(tX) && !valid(positions(theta, dimnamesAt(tX, 2)))) ||
                   (!Rf_isNull(tP) && !valid(positions(theta, dimnamesAt(tP, 1)))) ||
                   (!Rf_isNull(hX) && (!valid(positions(theta, dimnamesAt(hX, 2))) ||
                                       !valid(positions(theta, dimnamesAt(hX, 3))))) ||
                   (!Rf_isNull(hP) && (!valid(positions(theta, dimnamesAt(hP, 1))) ||
                                       !valid(positions(theta, dimnamesAt(hP, 2)))))))
      return none();
    if ((!Rf_isNull(tX) && dimAt(tX, 0) != n_obs) || (!Rf_isNull(hX) && dimAt(hX, 0) != n_obs)) return none();

    // fixed parameters seed zero
    std::vector<char> fixedRow(n_params, 0);
    SEXP fx = elt(a, "fixed");
    for (SEXP f : {stFixed, fx}) {
      if (Rf_isNull(f)) continue;
      for (int p : positions(f, parameters)) if (p >= 0) fixedRow[p] = 1;
    }

    SEXP s = Rf_allocVector(VECSXP, deriv2 ? 11 : 9);
    SET_VECTOR_ELT(outSets, k, s);
    SET_VECTOR_ELT(thetas, k, theta);
    SEXP M = Rf_allocVector(REALSXP, (R_xlen_t) n_obs * n_vars);
    SET_VECTOR_ELT(s, 0, M);
    for (int j = 0; j < n_vars; ++j)
      std::memcpy(REAL(M) + (size_t) n_obs * j, REAL(vars) + (size_t) n_obs * vi[j], sizeof(double) * n_obs);
    SEXP P = Rf_allocVector(REALSXP, n_params);
    SET_VECTOR_ELT(s, 1, P);
    for (int j = 0; j < n_params; ++j) REAL(P)[j] = REAL(params)[pi[j]];

    if (nt > 0 && n_vars > 0) {
      SEXP dX = Rf_allocVector(REALSXP, (R_xlen_t) n_obs * n_vars * nt);
      SET_VECTOR_ELT(s, 2, dX);
      std::fill(REAL(dX), REAL(dX) + Rf_xlength(dX), 0.0);
      if (!Rf_isNull(tX)) {
        const int ns = dimAt(tX, 1);
        std::vector<int> si = positions(innames, dimnamesAt(tX, 1));
        std::vector<int> ti = positions(theta, dimnamesAt(tX, 2));
        for (int t = 0; t < nt; ++t)
          for (int j = 0; j < n_vars; ++j)
            if (si[j] >= 0)
              std::memcpy(REAL(dX) + (size_t) n_obs * (j + (size_t) n_vars * t),
                          REAL(tX) + (size_t) n_obs * (si[j] + (size_t) ns * ti[t]), sizeof(double) * n_obs);
      }
    }
    if (nt > 0 && n_params > 0) {
      SEXP dP = Rf_allocVector(REALSXP, (R_xlen_t) n_params * nt);
      SET_VECTOR_ELT(s, 3, dP);
      std::fill(REAL(dP), REAL(dP) + Rf_xlength(dP), 0.0);
      if (!Rf_isNull(tP)) {
        const int nr = dimAt(tP, 0);
        std::vector<int> ri = positions(parameters, dimnamesAt(tP, 0));
        std::vector<int> ti = positions(theta, dimnamesAt(tP, 1));
        for (int t = 0; t < nt; ++t)
          for (int j = 0; j < n_params; ++j)
            if (ri[j] >= 0 && !fixedRow[j]) REAL(dP)[j + (size_t) n_params * t] = REAL(tP)[ri[j] + (size_t) nr * ti[t]];
      }
    }
    int at = 4;
    if (deriv2) {
      if (!Rf_isNull(hX) && nt > 0 && n_vars > 0) {
        SEXP d = Rf_allocVector(REALSXP, (R_xlen_t) n_obs * n_vars * nt * nt);
        SET_VECTOR_ELT(s, 4, d);
        std::fill(REAL(d), REAL(d) + Rf_xlength(d), 0.0);
        const int ns = dimAt(hX, 1), n3 = dimAt(hX, 2);
        std::vector<int> si = positions(innames, dimnamesAt(hX, 1));
        std::vector<int> t3 = positions(theta, dimnamesAt(hX, 2));
        std::vector<int> t4 = positions(theta, dimnamesAt(hX, 3));
        for (int u = 0; u < nt; ++u)
          for (int t = 0; t < nt; ++t)
            for (int j = 0; j < n_vars; ++j)
              if (si[j] >= 0)
                std::memcpy(REAL(d) + (size_t) n_obs * (j + (size_t) n_vars * (t + (size_t) nt * u)),
                            REAL(hX) + (size_t) n_obs * (si[j] + (size_t) ns * (t3[t] + (size_t) n3 * t4[u])),
                            sizeof(double) * n_obs);
      }
      if (!Rf_isNull(hP) && nt > 0 && n_params > 0) {
        SEXP d = Rf_allocVector(REALSXP, (R_xlen_t) n_params * nt * nt);
        SET_VECTOR_ELT(s, 5, d);
        std::fill(REAL(d), REAL(d) + Rf_xlength(d), 0.0);
        const int nr = dimAt(hP, 0), n2 = dimAt(hP, 1);
        std::vector<int> ri = positions(parameters, dimnamesAt(hP, 0));
        std::vector<int> t2 = positions(theta, dimnamesAt(hP, 1));
        std::vector<int> t3 = positions(theta, dimnamesAt(hP, 2));
        for (int u = 0; u < nt; ++u)
          for (int t = 0; t < nt; ++t)
            for (int j = 0; j < n_params; ++j)
              if (ri[j] >= 0 && !fixedRow[j])
                REAL(d)[j + (size_t) n_params * (t + (size_t) nt * u)] =
                    REAL(hP)[ri[j] + (size_t) nr * (t2[t] + (size_t) n2 * t3[u])];
      }
      at = 6;
    }
    const int sizes[5] = {n_obs, n_vars, n_params, n_out, nt};
    for (int q = 0; q < 5; ++q) SET_VECTOR_ELT(s, at + q, Rf_ScalarInteger(sizes[q]));
  }

  SEXP nm = Rf_allocVector(STRSXP, 2);
  Rf_setAttrib(out, R_NamesSymbol, nm);
  SET_STRING_ELT(nm, 0, Rf_mkChar("sets"));
  SET_STRING_ELT(nm, 1, Rf_mkChar("theta"));
  UNPROTECT(1);
  return out;
}

// raw: the results of the AD batch entry; outnames of the function, or a list of them per set;
// thetas: the direction names per set from cppde_eval_prep(). Returns per set
// list(y, tangent[, hessian]) with dimnames.
extern "C" SEXP cppde_eval_wrap(SEXP raw, SEXP outnamesIn, SEXP thetas, SEXP deriv2SEXP) {
  const bool deriv2 = Rf_asLogical(deriv2SEXP) == TRUE;
  const bool perSet = TYPEOF(outnamesIn) == VECSXP;
  const int K = Rf_length(raw);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, K));
  SEXP dnShared = PROTECT(Rf_allocVector(VECSXP, 2));
  if (!perSet) SET_VECTOR_ELT(dnShared, 1, outnamesIn);
  SEXP nm = PROTECT(Rf_allocVector(STRSXP, deriv2 ? 3 : 2));
  SET_STRING_ELT(nm, 0, Rf_mkChar("y"));
  SET_STRING_ELT(nm, 1, Rf_mkChar("tangent"));
  if (deriv2) SET_STRING_ELT(nm, 2, Rf_mkChar("hessian"));

  for (int k = 0; k < K; ++k) {
    SEXP outnames = perSet ? VECTOR_ELT(outnamesIn, k) : outnamesIn;
    const int n_out = Rf_length(outnames);
    SEXP r = VECTOR_ELT(raw, k);
    SEXP theta = VECTOR_ELT(thetas, k);
    const int nt = Rf_isNull(theta) ? 0 : Rf_length(theta);
    SEXP y = VECTOR_ELT(r, 0);
    const int n_obs = dimAt(y, 0);
    SEXP res = PROTECT(Rf_allocVector(VECSXP, deriv2 ? 3 : 2));
    SEXP dnY = dnShared;
    if (perSet) {
      dnY = Rf_allocVector(VECSXP, 2);
      SET_VECTOR_ELT(res, 0, dnY);
      SET_VECTOR_ELT(dnY, 1, outnames);
    }
    Rf_setAttrib(y, R_DimNamesSymbol, dnY);
    SET_VECTOR_ELT(res, 0, y);
    for (int order = 1; order <= (deriv2 ? 2 : 1); ++order) {
      SEXP d = VECTOR_ELT(r, order);
      const int nd = 2 + order;
      SEXP dn = PROTECT(Rf_allocVector(VECSXP, nd));
      SET_VECTOR_ELT(dn, 1, outnames);
      if (nt == 0) {
        // the entry writes one direction where there is none; the result has zero extent
        SEXP dims = PROTECT(Rf_allocVector(INTSXP, nd));
        INTEGER(dims)[0] = n_obs; INTEGER(dims)[1] = n_out;
        for (int q = 2; q < nd; ++q) INTEGER(dims)[q] = 0;
        d = PROTECT(Rf_allocArray(REALSXP, dims));
        Rf_setAttrib(d, R_DimNamesSymbol, dn);
        SET_VECTOR_ELT(res, order, d);
        UNPROTECT(3);
        continue;
      }
      for (int q = 2; q < nd; ++q) SET_VECTOR_ELT(dn, q, theta);
      Rf_setAttrib(d, R_DimNamesSymbol, dn);
      SET_VECTOR_ELT(res, order, d);
      UNPROTECT(1);
    }
    Rf_setAttrib(res, R_NamesSymbol, nm);
    SET_VECTOR_ELT(out, k, res);
    UNPROTECT(1);
  }
  UNPROTECT(3);
  return out;
}

using EvalAd = void (*)(double*, double*, double*, double*, double*, double*, int*, int*, int*, int*, int*);
using EvalAd2 = void (*)(double*, double*, double*, double*, double*, double*, int*, int*, double*, double*,
                         double*, int*, int*, int*, int*, int*);

// fns: the addresses of the raw AD entries (`<model>_eval_ad`, `_eval_ad2` with deriv2), one per
// set; sets: from cppde_eval_prep(). Each set by the function given for it, all in one parallel region;
// the results as the AD batch entries return them, NaN where a function throws.
extern "C" SEXP cppde_eval_many(SEXP fns, SEXP sets, SEXP deriv2SEXP, SEXP nthreadsSEXP) {
  const bool deriv2 = Rf_asLogical(deriv2SEXP) == TRUE;
  const int K = Rf_length(sets);
  if (Rf_length(fns) != K) Rf_error("cppde_eval_many: one function per set");
  const int at = deriv2 ? 6 : 4;

  std::vector<DL_FUNC> F(K);
  std::vector<double*> X(K), P(K), DX(K), DP(K), DX2(K), DP2(K), Y(K), DY(K), D2Y(K);
  std::vector<int> NOBS(K), NVAR(K), NPAR(K), NOUT(K), NTH(K), H2X(K), H2P(K);
  // an address that dyn.unload() has cleared returns NULL, for the caller to resolve it anew
  for (int k = 0; k < K; ++k) {
    SEXP f = VECTOR_ELT(fns, k);
    if (TYPEOF(f) != EXTPTRSXP || (F[k] = R_ExternalPtrAddrFn(f)) == nullptr) return R_NilValue;
  }
  SEXP out = PROTECT(Rf_allocVector(VECSXP, K));
  for (int k = 0; k < K; ++k) {
    SEXP s = VECTOR_ELT(sets, k);
    auto ptr = [&](int i) { SEXP v = VECTOR_ELT(s, i); return Rf_isNull(v) ? nullptr : REAL(v); };
    X[k] = ptr(0); P[k] = ptr(1); DX[k] = ptr(2); DP[k] = ptr(3);
    DX2[k] = deriv2 ? ptr(4) : nullptr; DP2[k] = deriv2 ? ptr(5) : nullptr;
    H2X[k] = DX2[k] != nullptr; H2P[k] = DP2[k] != nullptr;
    NOBS[k] = INTEGER(VECTOR_ELT(s, at))[0];
    NVAR[k] = INTEGER(VECTOR_ELT(s, at + 1))[0];
    NPAR[k] = INTEGER(VECTOR_ELT(s, at + 2))[0];
    NOUT[k] = INTEGER(VECTOR_ELT(s, at + 3))[0];
    NTH[k] = INTEGER(VECTOR_ELT(s, at + 4))[0];
    const int nth = NTH[k] > 0 ? NTH[k] : 1;

    SEXP res = Rf_allocVector(VECSXP, deriv2 ? 3 : 2);
    SET_VECTOR_ELT(out, k, res);
    SEXP y = Rf_allocMatrix(REALSXP, NOBS[k], NOUT[k]);
    SET_VECTOR_ELT(res, 0, y);
    Y[k] = REAL(y);
    SEXP d1 = PROTECT(Rf_allocVector(INTSXP, 3));
    INTEGER(d1)[0] = NOBS[k]; INTEGER(d1)[1] = NOUT[k]; INTEGER(d1)[2] = nth;
    SEXP dy = Rf_allocArray(REALSXP, d1);
    SET_VECTOR_ELT(res, 1, dy);
    DY[k] = REAL(dy);
    UNPROTECT(1);
    D2Y[k] = nullptr;
    if (deriv2) {
      SEXP d2 = PROTECT(Rf_allocVector(INTSXP, 4));
      INTEGER(d2)[0] = NOBS[k]; INTEGER(d2)[1] = NOUT[k]; INTEGER(d2)[2] = nth; INTEGER(d2)[3] = nth;
      SEXP d2y = Rf_allocArray(REALSXP, d2);
      SET_VECTOR_ELT(res, 2, d2y);
      D2Y[k] = REAL(d2y);
      UNPROTECT(1);
    }
  }

  const int nt = threadsFor(K, Rf_asInteger(nthreadsSEXP));
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(nt) if (nt > 1)
#endif
  for (int k = 0; k < K; ++k) {
    try {
      if (deriv2)
        reinterpret_cast<EvalAd2>(F[k])(X[k], P[k], DX[k], DP[k], DX2[k], DP2[k], &H2X[k], &H2P[k], Y[k], DY[k],
                                        D2Y[k], &NOBS[k], &NVAR[k], &NPAR[k], &NOUT[k], &NTH[k]);
      else
        reinterpret_cast<EvalAd>(F[k])(X[k], P[k], DX[k], DP[k], Y[k], DY[k], &NOBS[k], &NVAR[k], &NPAR[k],
                                       &NOUT[k], &NTH[k]);
    } catch (...) {
      const size_t ny = (size_t) NOBS[k] * NOUT[k], nth = NTH[k] > 0 ? NTH[k] : 1;
      const double nan = std::numeric_limits<double>::quiet_NaN();
      std::fill(Y[k], Y[k] + ny, nan);
      std::fill(DY[k], DY[k] + ny * nth, nan);
      if (D2Y[k]) std::fill(D2Y[k], D2Y[k] + ny * nth * nth, nan);
    }
  }
  UNPROTECT(1);
  return out;
}
