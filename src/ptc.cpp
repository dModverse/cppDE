/*
 Pseudo-transient continuation for the steady state of a cppFUN model.

 Solves G(x) = [f(x); C x - T] = 0 for the variables in `act`, f the model
 outputs `rows`, one per solved variable and in the same order. Each step solves

   (M / dt - J) d = G      (M the identity on the rows of f, zero on C)

 in the least-squares sense after scaling rows by their residual scale and
 columns by |x| + atol. Small dt gives implicit Euler steps along dx/dt = f,
 large dt Newton steps (Kelley & Keyes 1998; Coffey, Kelley & Keyes 2003). dt is
 controlled by the local error dt/2 (f(x+) - f(x)), filtered by (M - dt J)^-1,
 relative to flowTol |x| + atol; a step that halves the residual is accepted
 regardless. With flow = false the step follows the Newton flow,
 d = -J^+ G dt / (1 + dt). Positive variables are updated as x exp(d / x),
 exponent capped at +-zmax, and projected onto C x = T by x exp(C' lambda).
 Converged: every row within rtol of its turnover sum_k |dG_i/dx_k| |x_k| plus
 atol times its largest |dG_i/dx_k|, or, for dt >= dtNewton, a step within
 rtol |x| + atol.

 Copyright (C) 2026 Simon Beyer
 */

// no length(), error() ... macros: they break libc++ headers included below
#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Lapack.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>

#ifndef FCONE
#define FCONE
#endif

namespace {

typedef void (*eval_ad_fn)(double* x, double* p, double* dX, double* dP, double* y,
                           double* dy, int* n_obs, int* n_vars, int* n_params,
                           int* n_out, int* n_theta);

// The model with every variable outside `act` held at its value.
struct System {
  eval_ad_fn fn;
  int nv, np, no, na, nf, k;
  std::vector<double> xfull, p, dX, dP, y, dy;
  std::vector<int> act, rows;
  std::vector<double> C, T;  // k x na, column-major

  // G (nf + k) and J ((nf + k) x na, column-major) at the active values xa.
  bool eval(const std::vector<double>& xa, std::vector<double>& G, std::vector<double>& J) {
    for (int t = 0; t < na; ++t) xfull[act[t]] = xa[t];
    int one = 1, nth = na;
    fn(xfull.data(), p.data(), dX.data(), dP.data(), y.data(), dy.data(),
       &one, &nv, &np, &no, &nth);
    const int m = nf + k;
    for (int r = 0; r < nf; ++r) {
      G[r] = y[rows[r]];
      for (int t = 0; t < na; ++t) J[r + m * t] = dy[rows[r] + no * t];
    }
    for (int i = 0; i < k; ++i) {
      double s = -T[i];
      for (int t = 0; t < na; ++t) {
        s += C[i + k * t] * xa[t];
        J[nf + i + m * t] = C[i + k * t];
      }
      G[nf + i] = s;
    }
    for (double v : G) if (!std::isfinite(v)) return false;
    for (double v : J) if (!std::isfinite(v)) return false;
    return true;
  }
};

// QR of a tall matrix scaled by rows rs and columns cs; solves A d = b in the
// least-squares sense for several right-hand sides. Rank test per column as in
// R's qr(tol = 1e-14).
struct ScaledQR {
  int m, n;
  std::vector<double> A, tau, work, rs, cs;
  bool ok = false;

  bool factor(const std::vector<double>& A0, int m_, int n_,
              const std::vector<double>& rs_, const std::vector<double>& cs_) {
    m = m_; n = n_; rs = rs_; cs = cs_;
    A.assign(A0.begin(), A0.end());
    std::vector<double> norm(n, 0.0);
    for (int j = 0; j < n; ++j)
      for (int i = 0; i < m; ++i) {
        A[i + m * j] *= cs[j] / rs[i];
        norm[j] += A[i + m * j] * A[i + m * j];
      }
    tau.assign(n, 0.0);
    int lwork = -1, info = 0;
    double wq = 0;
    F77_CALL(dgeqrf)(&m, &n, A.data(), &m, tau.data(), &wq, &lwork, &info);
    lwork = std::max(1, static_cast<int>(wq));
    work.assign(lwork, 0.0);
    F77_CALL(dgeqrf)(&m, &n, A.data(), &m, tau.data(), work.data(), &lwork, &info);
    ok = (info == 0);
    for (int j = 0; ok && j < n; ++j)
      if (!(std::fabs(A[j + m * j]) > 1e-14 * std::sqrt(norm[j]))) ok = false;
    return ok;
  }

  bool solve(const std::vector<double>& b, std::vector<double>& d) {
    if (!ok) return false;
    std::vector<double> c(m);
    for (int i = 0; i < m; ++i) c[i] = b[i] / rs[i];
    int one = 1, info = 0, lwork = -1;
    double wq = 0;
    F77_CALL(dormqr)("L", "T", &m, &one, &n, A.data(), &m, tau.data(), c.data(), &m,
                     &wq, &lwork, &info FCONE FCONE);
    lwork = std::max(1, static_cast<int>(wq));
    std::vector<double> w(lwork);
    F77_CALL(dormqr)("L", "T", &m, &one, &n, A.data(), &m, tau.data(), c.data(), &m,
                     w.data(), &lwork, &info FCONE FCONE);
    if (info != 0) return false;
    F77_CALL(dtrtrs)("U", "N", "N", &n, &one, A.data(), &m, c.data(), &m, &info
                     FCONE FCONE FCONE);
    if (info != 0) return false;
    d.assign(n, 0.0);
    for (int j = 0; j < n; ++j) {
      d[j] = cs[j] * c[j];
      if (!std::isfinite(d[j])) return false;
    }
    return true;
  }
};

// x exp(C' lambda) with C x = T, lambda by damped Newton.
void project(std::vector<double>& x, const std::vector<double>& C, const std::vector<double>& T,
             int k, int n) {
  if (k == 0) return;
  std::vector<double> lam(k, 0.0), e(n), r(k), H(k * k), step(k), rn(k);
  auto residual = [&](const std::vector<double>& l, std::vector<double>& out) {
    for (int t = 0; t < n; ++t) {
      double s = 0;
      for (int i = 0; i < k; ++i) s += C[i + k * t] * l[i];
      e[t] = x[t] * std::exp(s);
    }
    for (int i = 0; i < k; ++i) {
      double s = -T[i];
      for (int t = 0; t < n; ++t) s += C[i + k * t] * e[t];
      out[i] = s;
    }
  };
  std::vector<double> tol(k);
  for (int i = 0; i < k; ++i) {
    double s = 0;
    for (int t = 0; t < n; ++t) s += std::fabs(C[i + k * t]) * x[t];
    tol[i] = 1e-15 * std::max(std::fabs(T[i]), s);
  }
  residual(lam, r);
  for (int it = 0; it < 100; ++it) {
    bool done = true;
    for (int i = 0; i < k; ++i) if (std::fabs(r[i]) > tol[i]) done = false;
    if (done) break;
    residual(lam, r);  // refreshes e at lam
    for (int a = 0; a < k; ++a)
      for (int b = 0; b < k; ++b) {
        double s = 0;
        for (int t = 0; t < n; ++t) s += C[a + k * t] * e[t] * C[b + k * t];
        H[a + k * b] = s;
      }
    step = r;
    std::vector<int> ipiv(k);
    int one = 1, info = 0;
    F77_CALL(dgesv)(&k, &one, H.data(), &k, ipiv.data(), step.data(), &k, &info);
    if (info != 0) break;
    double r2 = 0;
    for (double v : r) r2 += v * v;
    double alpha = 1;
    std::vector<double> lt(k);
    for (;;) {
      for (int i = 0; i < k; ++i) lt[i] = lam[i] - alpha * step[i];
      residual(lt, rn);
      double n2 = 0;
      bool fin = true;
      for (double v : rn) { n2 += v * v; if (!std::isfinite(v)) fin = false; }
      if ((fin && n2 < r2) || alpha < 1e-8) break;
      alpha /= 2;
    }
    lam = lt; r = rn;
  }
  for (int t = 0; t < n; ++t) {
    double s = 0;
    for (int i = 0; i < k; ++i) s += C[i + k * t] * lam[i];
    x[t] = x[t] * std::exp(s);
  }
}

double getControl(SEXP ctrl, const char* name) {
  SEXP nms = Rf_getAttrib(ctrl, R_NamesSymbol);
  for (int i = 0; i < Rf_length(ctrl); ++i)
    if (std::string(CHAR(STRING_ELT(nms, i))) == name)
      return Rf_asReal(VECTOR_ELT(ctrl, i));
  Rf_error("ptc: control '%s' missing", name);
  return 0;
}

}  // namespace

// fnS: external pointer to <model>_eval_ad; xS all variables; pS parameters;
// actS, rowsS 0-based indices; CS k x na; TS k; ctrlS named list.
extern "C" SEXP cppde_ptc(SEXP fnS, SEXP nOutS, SEXP xS, SEXP pS, SEXP actS, SEXP rowsS,
                          SEXP CS, SEXP TS, SEXP ctrlS) {
  System sys;
  sys.fn = reinterpret_cast<eval_ad_fn>(R_ExternalPtrAddrFn(fnS));
  if (sys.fn == nullptr) Rf_error("ptc: null function pointer");
  sys.nv = Rf_length(xS); sys.np = Rf_length(pS); sys.no = Rf_asInteger(nOutS);
  sys.na = Rf_length(actS); sys.nf = Rf_length(rowsS); sys.k = Rf_length(TS);
  if (sys.nf != sys.na) Rf_error("ptc: one output row per solved variable");
  sys.xfull.assign(REAL(xS), REAL(xS) + sys.nv);
  sys.p.assign(REAL(pS), REAL(pS) + sys.np);
  sys.act.assign(INTEGER(actS), INTEGER(actS) + sys.na);
  sys.rows.assign(INTEGER(rowsS), INTEGER(rowsS) + sys.nf);
  sys.C.assign(REAL(CS), REAL(CS) + sys.k * sys.na);
  sys.T.assign(REAL(TS), REAL(TS) + sys.k);
  sys.dX.assign(static_cast<size_t>(sys.nv) * sys.na, 0.0);
  for (int t = 0; t < sys.na; ++t) sys.dX[sys.act[t] + sys.nv * t] = 1.0;
  sys.dP.assign(static_cast<size_t>(std::max(sys.np, 1)) * sys.na, 0.0);
  sys.y.assign(sys.no, 0.0);
  sys.dy.assign(static_cast<size_t>(sys.no) * sys.na, 0.0);

  const double rtol = getControl(ctrlS, "rtol"), atol = getControl(ctrlS, "atol");
  const double flowTol = getControl(ctrlS, "flowTol");
  const int maxit = static_cast<int>(getControl(ctrlS, "maxit"));
  const double dtInit = getControl(ctrlS, "dtInit");
  const bool flow = getControl(ctrlS, "flow") != 0, positive = getControl(ctrlS, "positive") != 0;
  const double zmax = getControl(ctrlS, "zmax");

  const int n = sys.na, nf = sys.nf, m = nf + sys.k;
  std::vector<double> x(n), G(m), J(static_cast<size_t>(m) * n), Gn(m), Jn(J.size());
  for (int t = 0; t < n; ++t) x[t] = sys.xfull[sys.act[t]];

  std::string reason;
  bool ok = false;
  int it = 0;
  auto scaleOf = [&](const std::vector<double>& xv, const std::vector<double>& Jv,
                     std::vector<double>& sc) {
    sc.assign(m, 0.0);
    for (int r = 0; r < m; ++r) {
      double tv = 0;
      for (int t = 0; t < n; ++t) tv += std::fabs(Jv[r + m * t]) * std::fabs(xv[t]);
      double own = 1.0;
      if (r < nf) {
        own = 0;
        for (int t = 0; t < n; ++t) own = std::max(own, std::fabs(Jv[r + m * t]));
      }
      sc[r] = rtol * tv + atol * std::max(own, std::numeric_limits<double>::epsilon());
    }
  };
  auto rms = [&](const std::vector<double>& g, const std::vector<double>& sc) {
    double s = 0;
    for (int r = 0; r < m; ++r) s += (g[r] / sc[r]) * (g[r] / sc[r]);
    return std::sqrt(s / m);
  };
  auto maxsc = [&](const std::vector<double>& g, const std::vector<double>& sc) {
    double s = 0;
    for (int r = 0; r < m; ++r) s = std::max(s, std::fabs(g[r]) / sc[r]);
    return s;
  };

  std::vector<double> sc, scn;
  if (!sys.eval(x, G, J)) {
    reason = "non-finite residual at the start";
  } else {
    double rate = 1.0;
    if (flow) {
      rate = std::numeric_limits<double>::epsilon();
      for (int t = 0; t < nf; ++t) rate = std::max(rate, std::fabs(J[t + m * t]));
    }
    const double dt0 = dtInit / rate, dtNewton = 1e8 / rate;
    double dt = dt0;
    scaleOf(x, J, sc);
    double rn = rms(G, sc);
    ScaledQR qr;
    std::vector<double> A(J.size()), d, cs(n), xn(n), raw(m), lte;
    for (it = 0; it < maxit; ++it) {
      if (maxsc(G, sc) <= 1) { ok = true; reason = "converged"; break; }
      for (size_t q = 0; q < A.size(); ++q) A[q] = -J[q];
      if (flow) for (int t = 0; t < nf; ++t) A[t + m * t] += 1.0 / dt;
      for (int t = 0; t < n; ++t) cs[t] = std::fabs(x[t]) + atol;
      if (!qr.factor(A, m, n, sc, cs) || !qr.solve(G, d)) { dt /= 10; continue; }
      if (!flow) for (double& v : d) v *= dt / (1 + dt);
      bool small = dt >= dtNewton;
      for (int t = 0; small && t < n; ++t)
        small = std::fabs(d[t]) <= rtol * std::fabs(x[t]) + atol;
      if (positive) {
        for (int t = 0; t < n; ++t)
          xn[t] = x[t] * std::exp(std::min(zmax, std::max(-zmax, d[t] / x[t])));
        project(xn, sys.C, sys.T, sys.k, n);
      } else {
        for (int t = 0; t < n; ++t) xn[t] = x[t] + d[t];
      }
      if (!sys.eval(xn, Gn, Jn)) { dt /= 10; continue; }
      scaleOf(xn, Jn, scn);
      double rnn = rms(Gn, scn);
      bool progress = rnn <= 0.5 * rn;
      double grow;
      if (flow) {
        for (int r = 0; r < m; ++r) raw[r] = r < nf ? 0.5 * (Gn[r] - G[r]) : 0.0;
        if (!qr.solve(raw, lte)) {
          lte.assign(n, 0.0);
          for (int t = 0; t < n; ++t) lte[t] = dt * raw[t];
        }
        double err = 0;
        for (int t = 0; t < n; ++t) {
          double v = lte[t] / (flowTol * std::fabs(xn[t]) + atol);
          err += v * v;
        }
        err = std::sqrt(err / n);
        double fac = err > 0 ? 0.9 / std::sqrt(err) : 10.0;
        if (err > 1 && !progress) { dt *= std::max(0.1, fac); continue; }
        grow = progress ? std::max(2.0, std::min(1e3, rn / std::max(rnn, 1e-300)))
                        : std::min(10.0, std::max(0.2, fac));
      } else {
        if (rnn > 2 * rn && dt > 10 * dt0) { dt /= 4; continue; }
        grow = std::max(2.0, std::min(1e3, rn / std::max(rnn, 1e-300)));
      }
      x = xn; G = Gn; J = Jn; sc = scn;
      if (small) { ok = true; reason = "converged"; ++it; break; }
      dt = std::min(1e20, dt * grow);
      rn = rnn;
    }
    if (!ok) {
      if (maxsc(G, sc) <= 1) { ok = true; reason = "converged"; }
      else {
        char buf[160];
        std::snprintf(buf, sizeof buf, "no convergence in %d iterations (scaled residual %.2e)",
                      maxit, maxsc(G, sc));
        reason = buf;
      }
    }
  }

  SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
  SEXP nms = PROTECT(Rf_allocVector(STRSXP, 4));
  SEXP xo  = PROTECT(Rf_allocVector(REALSXP, n));
  for (int t = 0; t < n; ++t) REAL(xo)[t] = x[t];
  SET_VECTOR_ELT(out, 0, xo);
  SET_VECTOR_ELT(out, 1, Rf_ScalarLogical(ok));
  SET_VECTOR_ELT(out, 2, Rf_ScalarInteger(it));
  SET_VECTOR_ELT(out, 3, Rf_mkString(reason.c_str()));
  SET_STRING_ELT(nms, 0, Rf_mkChar("x"));
  SET_STRING_ELT(nms, 1, Rf_mkChar("converged"));
  SET_STRING_ELT(nms, 2, Rf_mkChar("iterations"));
  SET_STRING_ELT(nms, 3, Rf_mkChar("message"));
  Rf_setAttrib(out, R_NamesSymbol, nms);
  UNPROTECT(3);
  return out;
}
