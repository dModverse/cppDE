/*
 R marshalling for the generated solve entry points.

 The ONLY R-API-touching header under cppde/.  Never include it from
 cppde.hpp: the numerical core stays free of R so it can run off-thread.

 Layering:
   solve_impl(args, res)          pure C++, no R API, noexcept
   solve_<model>(SEXP x14)        thin wrapper, one condition
   solve_<model>_batch(conds, nt) reads all conditions, runs solve_impl over
                                  them (OpenMP), then builds the result list

 Phase rule: R API only before and after the parallel region.  R's allocator
 is not thread-safe, so a worker must never allocate an R object.  Plain
 C++ new/malloc from workers is fine.

 Copyright (C) 2026 Simon Beyer
 */

#ifndef CPPDE_R_BATCH_HPP
#define CPPDE_R_BATCH_HPP

#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

#include <algorithm>
#include <cstring>
#include <string>
#include <vector>

#ifndef _WIN32
#include <pthread.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

#include <cppde/cppde_blas_threads.hpp>
#include <cppde/cppde_return_codes.hpp>
#include <cppde/cppde_step_trace.hpp>

namespace cppde {
namespace rbatch {

// Per-condition inputs. Raw pointers into R memory, read serially in phase A
// and valid for the caller's .Call frame.
struct solve_args {
  const double* times    = nullptr;  int n_times = 0;
  const double* params   = nullptr;
  const double* sens1ini = nullptr;  int n_sens1 = 0;  int n_sens1_cols = 0;
  const double* sens2ini = nullptr;  int n_sens2 = 0;
  const int*    fixed    = nullptr;  int n_fixed = 0;
  double abstol = 1e-6, reltol = 1e-6, hini = 0.0, root_tol = 1e-6;
  int maxprogress = 50, maxsteps = 1000000, maxroot = 1;
  std::vector<const double*> ftimes, fvalues;
  std::vector<int> flen;
};

// Per-condition outputs. double only: AD tangents live in the thread-local
// arena and die with solve_impl's scope, so results must be flattened there.
struct solve_result {
  std::vector<double> time, variable, sens1, sens2;
  int n_out = 0, n_sens = 0;
  int return_code = RC_SUCCESS;
  std::string message;
  int accepted = 0, rejected = 0, fevals = 0, jevals = 0, setups = 0, last_order = 0;
  double last_dt = 0.0, t_reached = 0.0;
  ndf_detail::TraceBuffer trace;

  // Where the flatten writes.  Defaults into the vectors above; a caller
  // already on the R thread can redirect it straight into R memory through
  // `acquire`, which saves one full copy of the result.
  double *t_out = nullptr, *v_out = nullptr, *s1_out = nullptr, *s2_out = nullptr;
  bool (*acquire)(void*, solve_result&) = nullptr;
  void* acquire_ctx = nullptr;
  bool used_sink = false;   // true when the results live in caller-owned memory

  int fail(int rc, std::string msg) {
    return_code = rc;
    message = std::move(msg);
    return rc;
  }

  // Called by solve_impl once n_out and n_sens are known, before the flatten.
  void prepare_out(int n_out_, int n_sens_, int n_variables,
                   bool deriv, bool deriv2) {
    n_out = n_out_;
    n_sens = n_sens_;
    if (acquire && acquire(acquire_ctx, *this)) { used_sink = true; return; }

    time.resize((size_t)n_out);
    variable.resize((size_t)n_out * n_variables);
    if (deriv)  sens1.resize((size_t)n_out * n_variables * n_sens);
    if (deriv2) sens2.resize((size_t)n_out * n_variables * n_sens * n_sens);
    t_out  = time.data();
    v_out  = variable.data();
    s1_out = sens1.empty() ? nullptr : sens1.data();
    s2_out = sens2.empty() ? nullptr : sens2.data();
  }

};

  // ---------------------------------------------------------------------------
  //  Fork guard
  //
  //  libgomp's thread pool does not survive fork(): a child reusing it with more
  //  than one thread deadlocks. dMod forks in mstrust, so the child is marked at
  //  fork time and runs serially. The handler is installed at DSO load.
  // ---------------------------------------------------------------------------
namespace detail {

inline bool& child_of_fork() {
  static bool flag = false;
  return flag;
}

inline bool install_fork_guard() {
#ifndef _WIN32
  ::pthread_atfork(nullptr, nullptr, +[] { child_of_fork() = true; });
#endif
  return true;
}

inline const bool fork_guard_installed = install_fork_guard();

}  // namespace detail

  // The thread count actually used for a batch of K conditions. It is 1 inside a
  // forked child, and 1 inside an existing OpenMP region: a caller parallelising
  // over the wider axis owns the threads, and a second level oversubscribes.
inline int batch_threads(int K, int requested) {
#ifdef _OPENMP
  if (detail::child_of_fork()) return 1;
  if (omp_in_parallel()) return 1;
  int nt = (requested > 0) ? requested : omp_get_max_threads();
  if (nt > K) nt = K;
  return (nt < 1) ? 1 : nt;
#else
  (void)requested;
  (void)K;
  return 1;
#endif
}

  // libgomp is not instrumented, so ThreadSanitizer does not see the implicit
  // barrier at the end of a parallel region and reports every value written
  // inside and read after it as a race. These annotations supply the
  // happens-before edge the barrier already provides and compile to nothing
  // without -fsanitize=thread.
  // GCC before 14 does not define __has_feature, and #if does not short-circuit
  // at token level, so the call is routed through a macro that always exists.
#if defined(__has_feature)
#  define CPPDE_HAS_FEATURE(x) __has_feature(x)
#else
#  define CPPDE_HAS_FEATURE(x) 0
#endif

#if defined(__SANITIZE_THREAD__) || CPPDE_HAS_FEATURE(thread_sanitizer)
extern "C" void AnnotateHappensBefore(const char*, int, void*);
extern "C" void AnnotateHappensAfter(const char*, int, void*);
#  define CPPDE_TSAN_HB(a) AnnotateHappensBefore(__FILE__, __LINE__, (a))
#  define CPPDE_TSAN_HA(a) AnnotateHappensAfter(__FILE__, __LINE__, (a))
#else
#  define CPPDE_TSAN_HB(a) ((void)0)
#  define CPPDE_TSAN_HA(a) ((void)0)
#endif

// Report the thread count the batch actually used. `batch_threads` bails out
// silently inside a forked child and inside an enclosing OpenMP region, so a
// caller cannot tell from the requested value alone.
inline void set_threads_attr(SEXP out, int K, int requested) {
  SEXP v = PROTECT(Rf_ScalarInteger(batch_threads(K, requested)));
  Rf_setAttrib(out, Rf_install("threads"), v);
  UNPROTECT(1);
}

// Run fn(k) for k in [0, K). fn must be noexcept and must not touch the R API.
template <class Fn>
inline void run_batch(int K, int nthreads, Fn&& fn) {
  // One scope for the whole region: the BLAS thread count is process-global,
  // so the workers' own guards stand down and this serial-phase one restores
  // the caller's setting once, after the join.
  cppde::detail::single_thread_blas_scope blas_guard;
  const int nt = batch_threads(K, nthreads);
  char sync_in = 0, sync_out = 0;
  CPPDE_TSAN_HB(&sync_in);    // phase A is visible to every worker
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(nt) if (nt > 1)
#endif
  for (int k = 0; k < K; ++k) {
    CPPDE_TSAN_HA(&sync_in);
    fn(k);
    CPPDE_TSAN_HB(&sync_out);
  }
  CPPDE_TSAN_HA(&sync_out);   // every result is visible to phase C
  (void)nt; (void)sync_in; (void)sync_out;
}

// ---------------------------------------------------------------------------
//  Phase A: SEXP -> solve_args
// ---------------------------------------------------------------------------
inline solve_args read_solve_args(SEXP timesSEXP, SEXP paramsSEXP,
                                  SEXP sens1iniSEXP, SEXP sens2iniSEXP,
                                  SEXP fixedSEXP, SEXP abstolSEXP,
                                  SEXP reltolSEXP, SEXP maxprogressSEXP,
                                  SEXP maxstepsSEXP, SEXP hiniSEXP,
                                  SEXP root_tolSEXP, SEXP maxrootSEXP,
                                  SEXP forcingTimesSEXP, SEXP forcingValuesSEXP) {
  solve_args a;
  a.times   = REAL(timesSEXP);
  a.n_times = Rf_length(timesSEXP);
  a.params  = REAL(paramsSEXP);

  if (!Rf_isNull(sens1iniSEXP)) {
    a.sens1ini     = REAL(sens1iniSEXP);
    a.n_sens1      = Rf_length(sens1iniSEXP);
    a.n_sens1_cols = static_cast<int>(Rf_ncols(sens1iniSEXP));
  }
  if (!Rf_isNull(sens2iniSEXP)) {
    a.sens2ini = REAL(sens2iniSEXP);
    a.n_sens2  = Rf_length(sens2iniSEXP);
  }
  if (!Rf_isNull(fixedSEXP)) {
    a.fixed   = INTEGER(fixedSEXP);
    a.n_fixed = Rf_length(fixedSEXP);
  }

  a.abstol      = REAL(abstolSEXP)[0];
  a.reltol      = REAL(reltolSEXP)[0];
  a.hini        = REAL(hiniSEXP)[0];
  a.root_tol    = REAL(root_tolSEXP)[0];
  a.maxprogress = INTEGER(maxprogressSEXP)[0];
  a.maxsteps    = INTEGER(maxstepsSEXP)[0];
  a.maxroot     = INTEGER(maxrootSEXP)[0];

  const int nf = Rf_length(forcingTimesSEXP);
  a.ftimes.resize(nf);
  a.fvalues.resize(nf);
  a.flen.resize(nf);
  for (int fi = 0; fi < nf; ++fi) {
    SEXP ti = VECTOR_ELT(forcingTimesSEXP, fi);
    SEXP vi = VECTOR_ELT(forcingValuesSEXP, fi);
    a.ftimes[fi]  = REAL(ti);
    a.fvalues[fi] = REAL(vi);
    a.flen[fi]    = Rf_length(ti);
  }
  return a;
}

// Read one element of a batch condition list, falling back to R_NilValue.
inline SEXP cond_elt(SEXP cond, int i) {
  return (i < Rf_length(cond)) ? VECTOR_ELT(cond, i) : R_NilValue;
}

inline solve_args read_cond_args(SEXP cond) {
  return read_solve_args(cond_elt(cond, 0),  cond_elt(cond, 1),
                         cond_elt(cond, 2),  cond_elt(cond, 3),
                         cond_elt(cond, 4),  cond_elt(cond, 5),
                         cond_elt(cond, 6),  cond_elt(cond, 7),
                         cond_elt(cond, 8),  cond_elt(cond, 9),
                         cond_elt(cond, 10), cond_elt(cond, 11),
                         cond_elt(cond, 12), cond_elt(cond, 13));
}

// ---------------------------------------------------------------------------
//  Phase C: solve_result -> SEXP
// ---------------------------------------------------------------------------
inline SEXP build_diagnostics(const solve_result& r) {
  static const char* nms[10] = {"return_code", "message",   "accepted", "rejected",
                                "fevals",      "jevals",    "setups",   "last_dt",
                                "last_order",  "t_reached"};
  SEXP diag  = PROTECT(Rf_allocVector(VECSXP, 10));
  SEXP dn    = PROTECT(Rf_allocVector(STRSXP, 10));
  for (int i = 0; i < 10; ++i) SET_STRING_ELT(dn, i, Rf_mkChar(nms[i]));
  Rf_setAttrib(diag, R_NamesSymbol, dn);

  SET_VECTOR_ELT(diag, 0, Rf_ScalarInteger(r.return_code));
  SEXP msg = PROTECT(Rf_allocVector(STRSXP, 1));
  SET_STRING_ELT(msg, 0, Rf_mkChar(r.message.empty() ? "Integration was successful."
                                                     : r.message.c_str()));
  SET_VECTOR_ELT(diag, 1, msg);
  UNPROTECT(1);
  SET_VECTOR_ELT(diag, 2, Rf_ScalarInteger(r.accepted));
  SET_VECTOR_ELT(diag, 3, Rf_ScalarInteger(r.rejected));
  SET_VECTOR_ELT(diag, 4, Rf_ScalarInteger(r.fevals));
  SET_VECTOR_ELT(diag, 5, Rf_ScalarInteger(r.jevals));
  SET_VECTOR_ELT(diag, 6, Rf_ScalarInteger(r.setups));
  SET_VECTOR_ELT(diag, 7, Rf_ScalarReal(r.last_dt));
  SET_VECTOR_ELT(diag, 8, Rf_ScalarInteger(r.last_order));
  SET_VECTOR_ELT(diag, 9, Rf_ScalarReal(r.t_reached));
  UNPROTECT(2);
  return diag;
}

inline SEXP build_trace(const ndf_detail::TraceBuffer& tb) {
  const int n = static_cast<int>(tb.size());
  // Not tracing: 18 zero-length vectors plus a names vector, per condition, is
  // pure allocation in the serial phase. R drops the element anyway.
  if (n == 0) return R_NilValue;
  static const char* cn[18] = {
      "nst",   "t",      "h",       "q",        "dsm",          "acnrm",
      "acnrm_state", "tq2", "gamma", "gamrat",  "newton_conv",  "mode",
      "nfe",   "njev",   "nsetups", "setup_reason", "pece_iters", "pece_diverged"};
  SEXP lst = PROTECT(Rf_allocVector(VECSXP, 18));
  SEXP nm  = PROTECT(Rf_allocVector(STRSXP, 18));
  for (int i = 0; i < 18; ++i) SET_STRING_ELT(nm, i, Rf_mkChar(cn[i]));
  Rf_setAttrib(lst, R_NamesSymbol, nm);
  UNPROTECT(1);

  auto put_int = [&](int slot, const std::vector<int>& v) {
    SEXP s = PROTECT(Rf_allocVector(INTSXP, n));
    std::copy(v.begin(), v.end(), INTEGER(s));
    SET_VECTOR_ELT(lst, slot, s);
    UNPROTECT(1);
  };
  auto put_dbl = [&](int slot, const std::vector<double>& v) {
    SEXP s = PROTECT(Rf_allocVector(REALSXP, n));
    std::memcpy(REAL(s), v.data(), sizeof(double) * (size_t)n);
    SET_VECTOR_ELT(lst, slot, s);
    UNPROTECT(1);
  };
  auto put_str = [&](int slot, const std::vector<std::string>& v) {
    SEXP s = PROTECT(Rf_allocVector(STRSXP, n));
    for (int i = 0; i < n; ++i) SET_STRING_ELT(s, i, Rf_mkChar(v[i].c_str()));
    SET_VECTOR_ELT(lst, slot, s);
    UNPROTECT(1);
  };

  put_int(0, tb.nst);          put_dbl(1, tb.t);            put_dbl(2, tb.h);
  put_int(3, tb.q);            put_dbl(4, tb.dsm);          put_dbl(5, tb.acnrm);
  put_dbl(6, tb.acnrm_state);  put_dbl(7, tb.tq2);          put_dbl(8, tb.gamma);
  put_dbl(9, tb.gamrat);       put_int(10, tb.newton_conv); put_str(11, tb.mode);
  put_int(12, tb.nfe);         put_int(13, tb.njev);        put_int(14, tb.nsetups);
  put_str(15, tb.setup_reason);put_int(16, tb.pece_iters);  put_int(17, tb.pece_diverged);

  UNPROTECT(1);
  return lst;
}

// Build list(time, variable, [sens1], [sens2], diagnostics, trace).
inline SEXP build_result_sexp(const solve_result& r, int n_variables,
                              bool deriv, bool deriv2) {
  const int n_out  = r.n_out;
  const int n_sens = r.n_sens;
  const int n_el   = 3 + (deriv ? 1 : 0) + (deriv2 ? 1 : 0) + 1;

  SEXP ans   = PROTECT(Rf_allocVector(VECSXP, n_el));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, n_el));
  int slot = 0;
  SET_STRING_ELT(names, slot++, Rf_mkChar("time"));
  SET_STRING_ELT(names, slot++, Rf_mkChar("variable"));
  if (deriv)  SET_STRING_ELT(names, slot++, Rf_mkChar("sens1"));
  if (deriv2) SET_STRING_ELT(names, slot++, Rf_mkChar("sens2"));
  SET_STRING_ELT(names, slot++, Rf_mkChar("diagnostics"));
  SET_STRING_ELT(names, slot++, Rf_mkChar("trace"));
  Rf_setAttrib(ans, R_NamesSymbol, names);
  UNPROTECT(1);

  slot = 0;
  SEXP tv = PROTECT(Rf_allocVector(REALSXP, n_out));
  std::memcpy(REAL(tv), r.time.data(), sizeof(double) * (size_t)n_out);
  SET_VECTOR_ELT(ans, slot++, tv);
  UNPROTECT(1);

  SEXP vm = PROTECT(Rf_allocMatrix(REALSXP, n_out, n_variables));
  std::memcpy(REAL(vm), r.variable.data(), sizeof(double) * r.variable.size());
  SET_VECTOR_ELT(ans, slot++, vm);
  UNPROTECT(1);

  if (deriv) {
    SEXP d = PROTECT(Rf_allocVector(INTSXP, 3));
    INTEGER(d)[0] = n_out; INTEGER(d)[1] = n_variables; INTEGER(d)[2] = n_sens;
    SEXP a = PROTECT(Rf_allocArray(REALSXP, d));
    std::memcpy(REAL(a), r.sens1.data(), sizeof(double) * r.sens1.size());
    SET_VECTOR_ELT(ans, slot++, a);
    UNPROTECT(2);
  }
  if (deriv2) {
    SEXP d = PROTECT(Rf_allocVector(INTSXP, 4));
    INTEGER(d)[0] = n_out;  INTEGER(d)[1] = n_variables;
    INTEGER(d)[2] = n_sens; INTEGER(d)[3] = n_sens;
    SEXP a = PROTECT(Rf_allocArray(REALSXP, d));
    std::memcpy(REAL(a), r.sens2.data(), sizeof(double) * r.sens2.size());
    SET_VECTOR_ELT(ans, slot++, a);
    UNPROTECT(2);
  }

  SEXP diag = PROTECT(build_diagnostics(r));
  SET_VECTOR_ELT(ans, slot++, diag);
  UNPROTECT(1);

  SEXP tr = PROTECT(build_trace(r.trace));
  SET_VECTOR_ELT(ans, slot++, tr);
  UNPROTECT(1);

  UNPROTECT(1);
  return ans;
}

  // ---------------------------------------------------------------------------
  //  Single-condition path
  //
  //  n_out and n_sens are known before the flatten runs, so the R objects are
  //  allocated there and solve_impl gathers straight into them, which keeps the
  //  single solve free of the staging copy the batch path pays.
  // ---------------------------------------------------------------------------
  // Dimnames are attached while the array is still unaliased. Doing it in R on
  // the returned list would duplicate every array, the element being referenced
  // by both the list and the replacement call.
inline void set_dimnames(SEXP x, SEXP var_nm, SEXP sens_nm, int n_sens_dims) {
  if (Rf_isNull(var_nm)) return;
  const int nd = 2 + n_sens_dims;
  SEXP d = PROTECT(Rf_allocVector(VECSXP, nd));
  SET_VECTOR_ELT(d, 0, R_NilValue);
  SET_VECTOR_ELT(d, 1, var_nm);
  for (int i = 0; i < n_sens_dims; ++i) SET_VECTOR_ELT(d, 2 + i, sens_nm);
  if (n_sens_dims > 0) {          // the serial path names these; match it
    SEXP nm = PROTECT(Rf_allocVector(STRSXP, nd));
    SET_STRING_ELT(nm, 0, Rf_mkChar("time"));
    SET_STRING_ELT(nm, 1, Rf_mkChar("variable"));
    if (n_sens_dims == 1) {
      SET_STRING_ELT(nm, 2, Rf_mkChar("sens"));
    } else {
      SET_STRING_ELT(nm, 2, Rf_mkChar("sens1"));
      SET_STRING_ELT(nm, 3, Rf_mkChar("sens2"));
    }
    Rf_setAttrib(d, R_NamesSymbol, nm);
    UNPROTECT(1);
  }
  Rf_setAttrib(x, R_DimNamesSymbol, d);
  UNPROTECT(1);
}

// Splits list(variables, sens_names) into its two parts; both R_NilValue when
// the caller passed nothing.
inline void split_dimnames(SEXP dn, SEXP* var_nm, SEXP* sens_nm) {
  *var_nm = R_NilValue; *sens_nm = R_NilValue;
  if (!Rf_isNull(dn) && Rf_length(dn) >= 2) {
    *var_nm  = VECTOR_ELT(dn, 0);
    *sens_nm = VECTOR_ELT(dn, 1);
  }
}

// `dn` is either list(variables, sens) shared by every condition, or a list of
// K such pairs. Conditions rarely agree on their sensitivity labels -- a
// reparametrised chain gives each one its own -- and setting them in R
// afterwards duplicates the whole sensitivity array.
inline SEXP dimnames_for(SEXP dn, int k, int K) {
  if (Rf_isNull(dn)) return R_NilValue;
  if (Rf_length(dn) == K && K != 2) return VECTOR_ELT(dn, k);
  // length 2 is ambiguous when K == 2: a pair of pairs is per-condition.
  if (Rf_length(dn) == 2 && K == 2) {
    SEXP first = VECTOR_ELT(dn, 0);
    if (!Rf_isNull(first) && TYPEOF(first) == VECSXP) return VECTOR_ELT(dn, k);
  }
  return dn;
}

struct r_out_ctx {
  int n_variables = 0;
  bool deriv = false, deriv2 = false;
  SEXP var_nm = R_NilValue, sens_nm = R_NilValue;
  SEXP time = R_NilValue, variable = R_NilValue;
  SEXP sens1 = R_NilValue, sens2 = R_NilValue;
  int nprot = 0;
};

inline bool r_out_acquire(void* vctx, solve_result& r) {
  r_out_ctx* c = static_cast<r_out_ctx*>(vctx);
  const int n_out = r.n_out, n_sens = r.n_sens;

  c->time     = PROTECT(Rf_allocVector(REALSXP, n_out));                   ++c->nprot;
  c->variable = PROTECT(Rf_allocMatrix(REALSXP, n_out, c->n_variables));   ++c->nprot;
  set_dimnames(c->variable, c->var_nm, R_NilValue, 0);
  r.t_out = REAL(c->time);
  r.v_out = REAL(c->variable);

  if (c->deriv) {
    SEXP d = PROTECT(Rf_allocVector(INTSXP, 3));                           ++c->nprot;
    INTEGER(d)[0] = n_out; INTEGER(d)[1] = c->n_variables; INTEGER(d)[2] = n_sens;
    c->sens1 = PROTECT(Rf_allocArray(REALSXP, d));                         ++c->nprot;
    set_dimnames(c->sens1, c->var_nm, c->sens_nm, 1);
    r.s1_out = REAL(c->sens1);
  }
  if (c->deriv2) {
    SEXP d = PROTECT(Rf_allocVector(INTSXP, 4));                           ++c->nprot;
    INTEGER(d)[0] = n_out;  INTEGER(d)[1] = c->n_variables;
    INTEGER(d)[2] = n_sens; INTEGER(d)[3] = n_sens;
    c->sens2 = PROTECT(Rf_allocArray(REALSXP, d));                         ++c->nprot;
    set_dimnames(c->sens2, c->var_nm, c->sens_nm, 2);
    r.s2_out = REAL(c->sens2);
  }
  return true;
}

inline SEXP solve_one(const solve_args& a, int n_variables, bool deriv, bool deriv2,
                      int (*impl)(const solve_args&, solve_result&),
                      SEXP dn = R_NilValue) {
  r_out_ctx ctx;
  ctx.n_variables = n_variables;
  ctx.deriv = deriv;
  ctx.deriv2 = deriv2;
  split_dimnames(dn, &ctx.var_nm, &ctx.sens_nm);

  solve_result r;
  r.acquire     = &r_out_acquire;
  r.acquire_ctx = &ctx;
  impl(a, r);

  // Failing before the flatten leaves ctx empty; hand back zero-row results
  // so the R side always sees the same shape.
  if (ctx.time == R_NilValue) {
    ctx.time     = PROTECT(Rf_allocVector(REALSXP, 0));                    ++ctx.nprot;
    ctx.variable = PROTECT(Rf_allocMatrix(REALSXP, 0, n_variables));       ++ctx.nprot;
    if (deriv) {
      SEXP d = PROTECT(Rf_allocVector(INTSXP, 3));                         ++ctx.nprot;
      INTEGER(d)[0] = 0; INTEGER(d)[1] = n_variables; INTEGER(d)[2] = r.n_sens;
      ctx.sens1 = PROTECT(Rf_allocArray(REALSXP, d));                      ++ctx.nprot;
    }
    if (deriv2) {
      SEXP d = PROTECT(Rf_allocVector(INTSXP, 4));                         ++ctx.nprot;
      INTEGER(d)[0] = 0;         INTEGER(d)[1] = n_variables;
      INTEGER(d)[2] = r.n_sens;  INTEGER(d)[3] = r.n_sens;
      ctx.sens2 = PROTECT(Rf_allocArray(REALSXP, d));                      ++ctx.nprot;
    }
  }

  const int n_el = 3 + (deriv ? 1 : 0) + (deriv2 ? 1 : 0) + 1;
  SEXP ans   = PROTECT(Rf_allocVector(VECSXP, n_el));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, n_el));
  int slot = 0;
  SET_STRING_ELT(names, slot++, Rf_mkChar("time"));
  SET_STRING_ELT(names, slot++, Rf_mkChar("variable"));
  if (deriv)  SET_STRING_ELT(names, slot++, Rf_mkChar("sens1"));
  if (deriv2) SET_STRING_ELT(names, slot++, Rf_mkChar("sens2"));
  SET_STRING_ELT(names, slot++, Rf_mkChar("diagnostics"));
  SET_STRING_ELT(names, slot++, Rf_mkChar("trace"));
  Rf_setAttrib(ans, R_NamesSymbol, names);

  slot = 0;
  SET_VECTOR_ELT(ans, slot++, ctx.time);
  SET_VECTOR_ELT(ans, slot++, ctx.variable);
  if (deriv)  SET_VECTOR_ELT(ans, slot++, ctx.sens1);
  if (deriv2) SET_VECTOR_ELT(ans, slot++, ctx.sens2);

  SEXP diag = PROTECT(build_diagnostics(r));
  SET_VECTOR_ELT(ans, slot++, diag);
  SEXP tr = PROTECT(build_trace(r.trace));
  SET_VECTOR_ELT(ans, slot++, tr);

  UNPROTECT(4 + ctx.nprot);   // tr, diag, names, ans, then everything ctx took
  return ans;
}

  // ---------------------------------------------------------------------------
  //  Pre-allocated batch output
  //
  //  With the output grid fixed by `times`, n_out is known before the solve, so
  //  the R objects are allocated on the main thread and the workers gather into
  //  them. Allocation stays single-threaded, R's allocator is not thread-safe;
  //  only the writes are parallel, which also moves the first-touch page faults
  //  into the parallel region.
  // ---------------------------------------------------------------------------
struct pre_ctx {
  int n_out = 0, n_sens = 0;
  double *t = nullptr, *v = nullptr, *s1 = nullptr, *s2 = nullptr;
};

// Falls back to solve_result's own vectors if the prediction was wrong, so a
// mispredicted n_out degrades in performance rather than writing out of bounds.
inline bool pre_acquire(void* vctx, solve_result& r) {
  pre_ctx* c = static_cast<pre_ctx*>(vctx);
  if (r.n_out != c->n_out || r.n_sens != c->n_sens) return false;
  r.t_out = c->t; r.v_out = c->v; r.s1_out = c->s1; r.s2_out = c->s2;
  return true;
}

// Number of output points `times` will produce after the generated code's
// zero-injection, sort and unique. A fixed event whose time is not already a
// requested time adds a row of its own, so those times count too; the caller
// passes them when it can evaluate them up front (no root event, no forcing).
inline int processed_time_count(const double* t, int n, bool include_zero,
                                const double* ev = nullptr, int n_ev = 0) {
  std::vector<double> v(t, t + n);
  if (include_zero && std::find(v.begin(), v.end(), 0.0) == v.end()) v.push_back(0.0);
  for (int i = 0; i < n_ev; ++i) v.push_back(ev[i]);
  std::sort(v.begin(), v.end());
  v.erase(std::unique(v.begin(), v.end()), v.end());
  return static_cast<int>(v.size());
}

// Mirrors the generated code's n_sens: explicit Phi' width, else the
// compile-time count minus the distinct runtime-fixed indices.
inline int sens_width(const solve_args& a, int n_sens_total) {
  if (a.sens1ini != nullptr) return a.n_sens1_cols;
  if (n_sens_total <= 0) return 0;
  std::vector<char> fixed(n_sens_total, 0);
  for (int i = 0; i < a.n_fixed; ++i)
    if (a.fixed[i] >= 0 && a.fixed[i] < n_sens_total) fixed[a.fixed[i]] = 1;
  int nf = 0;
  for (int i = 0; i < n_sens_total; ++i) nf += fixed[i];
  return n_sens_total - nf;
}

  // Builds each condition's result skeleton into `out`, already protected by the
  // caller, and returns the buffers the workers write into. `dn` is R_NilValue
  // or list(variables, sens_names); setting the dimnames here rather than in R
  // avoids duplicating every array, which a refcount above one would force.
  // `ev_times` is empty or one vector of fixed-event times per condition.
inline std::vector<pre_ctx> prealloc_batch(SEXP out, const std::vector<solve_args>& A,
                                           int n_variables, int n_sens_total,
                                           bool deriv, bool deriv2, bool include_zero,
                                           SEXP dn = R_NilValue,
                                           const std::vector<std::vector<double> >* ev_times = nullptr) {
  SEXP var_nm = R_NilValue, sens_nm = R_NilValue;
  split_dimnames(dn, &var_nm, &sens_nm);
  const int K = static_cast<int>(A.size());
  std::vector<pre_ctx> ctx(K);
  const int n_el = 3 + (deriv ? 1 : 0) + (deriv2 ? 1 : 0) + 1;

  for (int k = 0; k < K; ++k) {
    SEXP dn_k = dimnames_for(dn, k, K);
    split_dimnames(dn_k, &var_nm, &sens_nm);
    const bool has_ev = ev_times != nullptr && (int) ev_times->size() > k &&
                        !(*ev_times)[k].empty();
    const int n_out  = processed_time_count(
        A[k].times, A[k].n_times, include_zero,
        has_ev ? (*ev_times)[k].data() : nullptr,
        has_ev ? (int) (*ev_times)[k].size() : 0);
    const int n_sens = deriv ? sens_width(A[k], n_sens_total) : 0;
    ctx[k].n_out = n_out;
    ctx[k].n_sens = n_sens;

    SEXP ans   = PROTECT(Rf_allocVector(VECSXP, n_el));
    SEXP names = PROTECT(Rf_allocVector(STRSXP, n_el));
    int slot = 0;
    SET_STRING_ELT(names, slot++, Rf_mkChar("time"));
    SET_STRING_ELT(names, slot++, Rf_mkChar("variable"));
    if (deriv)  SET_STRING_ELT(names, slot++, Rf_mkChar("sens1"));
    if (deriv2) SET_STRING_ELT(names, slot++, Rf_mkChar("sens2"));
    SET_STRING_ELT(names, slot++, Rf_mkChar("diagnostics"));
    SET_STRING_ELT(names, slot++, Rf_mkChar("trace"));
    Rf_setAttrib(ans, R_NamesSymbol, names);
    UNPROTECT(1);  // names

    slot = 0;
    SEXP tv = PROTECT(Rf_allocVector(REALSXP, n_out));
    ctx[k].t = REAL(tv);
    SET_VECTOR_ELT(ans, slot++, tv);
    UNPROTECT(1);

    SEXP vm = PROTECT(Rf_allocMatrix(REALSXP, n_out, n_variables));
    ctx[k].v = REAL(vm);
    set_dimnames(vm, var_nm, R_NilValue, 0);
    SET_VECTOR_ELT(ans, slot++, vm);
    UNPROTECT(1);

    if (deriv) {
      SEXP d = PROTECT(Rf_allocVector(INTSXP, 3));
      INTEGER(d)[0] = n_out; INTEGER(d)[1] = n_variables; INTEGER(d)[2] = n_sens;
      SEXP a1 = PROTECT(Rf_allocArray(REALSXP, d));
      ctx[k].s1 = REAL(a1);
      set_dimnames(a1, var_nm, sens_nm, 1);
      SET_VECTOR_ELT(ans, slot++, a1);
      UNPROTECT(2);
    }
    if (deriv2) {
      SEXP d = PROTECT(Rf_allocVector(INTSXP, 4));
      INTEGER(d)[0] = n_out;  INTEGER(d)[1] = n_variables;
      INTEGER(d)[2] = n_sens; INTEGER(d)[3] = n_sens;
      SEXP a2 = PROTECT(Rf_allocArray(REALSXP, d));
      ctx[k].s2 = REAL(a2);
      set_dimnames(a2, var_nm, sens_nm, 2);
      SET_VECTOR_ELT(ans, slot++, a2);
      UNPROTECT(2);
    }

    SET_VECTOR_ELT(out, k, ans);
    UNPROTECT(1);  // ans, now kept alive by out
  }
  return ctx;
}

// Phase C for the pre-allocated path: only diagnostics and trace are left.
inline void finish_prealloc(SEXP out, int k, const solve_result& r,
                            bool deriv, bool deriv2) {
  SEXP ans = VECTOR_ELT(out, k);
  const int base = 2 + (deriv ? 1 : 0) + (deriv2 ? 1 : 0);
  SEXP diag = PROTECT(build_diagnostics(r));
  SET_VECTOR_ELT(ans, base, diag);
  UNPROTECT(1);
  SEXP tr = PROTECT(build_trace(r.trace));
  SET_VECTOR_ELT(ans, base + 1, tr);
  UNPROTECT(1);
}

}  // namespace rbatch
}  // namespace cppde

#endif  // CPPDE_R_BATCH_HPP
