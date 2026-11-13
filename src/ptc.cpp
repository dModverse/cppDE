/*
 Pseudo-transient continuation for the steady state of a cppFUN model.

 Solves G(x) = [f(x); C x - T] = 0 for the variables in `act`, f the model
 outputs `rows`, one per solved variable and in the same order. Each step solves

   (M / dt - J) d = G      (M the identity on the rows of f, zero on C)

 in the least-squares sense after scaling rows by their residual scale and
 columns by |x| + atol. Small dt gives implicit Euler steps along dx/dt = f,
 large dt Newton steps (Kelley & Keyes 1998; Coffey, Kelley & Keyes 2003). dt is
 controlled by the local error dt/2 (f(x+) - f(x)), filtered by (M - dt J)^-1,
 relative to flowTol |x| + atol; once dt reaches the Newton range, a step that
 halves the residual is accepted regardless. With flow = false the step follows
 the Newton flow, d = -J^+ G dt / (1 + dt). Positive variables are updated as
 x exp(d / x), exponent capped at +-zmax, and projected onto C x = T by
 x exp(C' lambda).

 Converged: every row within rtol of the row turnover sum_k |dG_i/dx_k| |x_k| plus
 atol times the largest |dG_i/dx_k| of the row, or, for dt >= dtNewton, a step within
 rtol |x| + atol. Newton steps from there bring every variable to within
 rtol (|x| + atol) of the root, as long as they contract.

 With flow and stability, the root is checked on ker C: where the leading
 eigenvalue lambda of J has a positive real part, the search leaves along the
 eigenvector, first to the side the iterates came from while their steps kept
 dt lambda < 1, with dt <= 1 / (2 lambda) until it is away from the root.

 Copyright (C) 2026 Simon Beyer
 */

// no length(), error() ... macros: they break libc++ headers included below
#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Lapack.h>
#include <R_ext/Rdynload.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
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

const double kInf = std::numeric_limits<double>::infinity();

struct Controls {
  double rtol = 1e-10, atol = 1e-14, flowTol = 1, dtInit = 1e-2, zmax = 5;
  int maxit = 400;
  bool flow = true, positive = true, stability = true;
  // a root whose leading eigenvalue on ker C has at most this real part counts as stable
  double stabilityTol = 0;
};

// The model with every variable outside `act` held at the value in `xfull`.
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

// The leading eigenvalue of the rows of f in J on ker C: the real and imaginary part, the
// largest modulus, and the real part of the eigenvector. Without ker C, re = -inf.
struct Lead {
  double re = -kInf, im = 0, scale = 0;
  std::vector<double> v;
  bool ok = true;
};

Lead leading(const std::vector<double>& J, int m, int n, const std::vector<double>& C, int k) {
  Lead L;
  // an orthonormal basis N of ker C: the right singular vectors past the rank
  std::vector<double> N;
  int mm = n;
  if (k > 0) {
    std::vector<double> A(C), s(std::min(k, n)), vt((std::size_t) n * n), u(1);
    int lwork = -1, info = 0, one = 1;
    double wq = 0;
    F77_CALL(dgesvd)("N", "A", &k, &n, A.data(), &k, s.data(), u.data(), &one, vt.data(), &n,
                     &wq, &lwork, &info FCONE FCONE);
    lwork = std::max(1, static_cast<int>(wq));
    std::vector<double> w(lwork);
    F77_CALL(dgesvd)("N", "A", &k, &n, A.data(), &k, s.data(), u.data(), &one, vt.data(), &n,
                     w.data(), &lwork, &info FCONE FCONE);
    if (info != 0) { L.ok = false; return L; }
    const double tol = std::max(k, n) * (s.empty() ? 0.0 : s[0]) * std::numeric_limits<double>::epsilon();
    int r = 0;
    for (double sv : s) if (sv > tol) ++r;
    mm = n - r;
    N.assign((std::size_t) n * mm, 0.0);
    for (int c = 0; c < mm; ++c)
      for (int t = 0; t < n; ++t) N[t + (std::size_t) n * c] = vt[(r + c) + (std::size_t) n * t];
  } else {
    N.assign((std::size_t) n * n, 0.0);
    for (int t = 0; t < n; ++t) N[t + (std::size_t) n * t] = 1.0;
  }
  if (mm == 0) return L;
  // B = N' Jf N, Jf the first n rows of J
  std::vector<double> JN((std::size_t) n * mm, 0.0), B((std::size_t) mm * mm, 0.0);
  for (int c = 0; c < mm; ++c)
    for (int t = 0; t < n; ++t) {
      const double nt = N[t + (std::size_t) n * c];
      if (nt == 0.0) continue;
      for (int r = 0; r < n; ++r) JN[r + (std::size_t) n * c] += J[r + (std::size_t) m * t] * nt;
    }
  for (int c = 0; c < mm; ++c)
    for (int a = 0; a < mm; ++a) {
      double s = 0;
      for (int r = 0; r < n; ++r) s += N[r + (std::size_t) n * a] * JN[r + (std::size_t) n * c];
      B[a + (std::size_t) mm * c] = s;
    }
  std::vector<double> wr(mm), wi(mm), vr((std::size_t) mm * mm), vl(1);
  int lwork = -1, info = 0, one = 1;
  double wq = 0;
  F77_CALL(dgeev)("N", "V", &mm, B.data(), &mm, wr.data(), wi.data(), vl.data(), &one, vr.data(), &mm,
                  &wq, &lwork, &info FCONE FCONE);
  lwork = std::max(1, static_cast<int>(wq));
  std::vector<double> w(lwork);
  F77_CALL(dgeev)("N", "V", &mm, B.data(), &mm, wr.data(), wi.data(), vl.data(), &one, vr.data(), &mm,
                  w.data(), &lwork, &info FCONE FCONE);
  if (info != 0) { L.ok = false; return L; }
  int j = 0;
  for (int i = 0; i < mm; ++i) {
    if (wr[i] > wr[j]) j = i;
    L.scale = std::max(L.scale, std::hypot(wr[i], wi[i]));
  }
  L.re = wr[j];
  L.im = std::fabs(wi[j]);
  // the real part of a complex pair's eigenvector is the first column of the pair
  const int col = (wi[j] < 0 && j > 0) ? j - 1 : j;
  L.v.assign(n, 0.0);
  for (int c = 0; c < mm; ++c) {
    const double z = vr[c + (std::size_t) mm * col];
    for (int t = 0; t < n; ++t) L.v[t] += N[t + (std::size_t) n * c] * z;
  }
  return L;
}

// A run of the solver on one system.
struct Ptc {
  System& sys;
  Controls c;
  int n, nf, m;
  int it = 0;
  std::vector<double> G, J;
  // the accepted iterates of the last descent and the dt of the step that reached each
  std::vector<std::vector<double>> hist;
  std::vector<double> histDt;

  Ptc(System& s, const Controls& cc) : sys(s), c(cc), n(s.na), nf(s.nf), m(s.nf + s.k) {
    G.assign(m, 0.0);
    J.assign((std::size_t) m * n, 0.0);
  }

  void scaleOf(const std::vector<double>& xv, const std::vector<double>& Jv, std::vector<double>& sc) const {
    sc.assign(m, 0.0);
    for (int r = 0; r < m; ++r) {
      double tv = 0;
      for (int t = 0; t < n; ++t) tv += std::fabs(Jv[r + m * t]) * std::fabs(xv[t]);
      double own = 1.0;
      if (r < nf) {
        own = 0;
        for (int t = 0; t < n; ++t) own = std::max(own, std::fabs(Jv[r + m * t]));
      }
      sc[r] = c.rtol * tv + c.atol * std::max(own, std::numeric_limits<double>::epsilon());
    }
  }
  double rms(const std::vector<double>& g, const std::vector<double>& sc) const {
    double s = 0;
    for (int r = 0; r < m; ++r) s += (g[r] / sc[r]) * (g[r] / sc[r]);
    return std::sqrt(s / m);
  }
  double maxsc(const std::vector<double>& g, const std::vector<double>& sc) const {
    double s = 0;
    for (int r = 0; r < m; ++r) s = std::max(s, std::fabs(g[r]) / sc[r]);
    return s;
  }
  void update(const std::vector<double>& x, const std::vector<double>& d, std::vector<double>& xn) const {
    if (c.positive) {
      for (int t = 0; t < n; ++t)
        xn[t] = x[t] * std::exp(std::min(c.zmax, std::max(-c.zmax, d[t] / x[t])));
      project(xn, sys.C, sys.T, sys.k, n);
    } else {
      for (int t = 0; t < n; ++t) xn[t] = x[t] + d[t];
    }
  }

  // PTC from x, which becomes the last iterate, G and J at it. With `xc`, dt stays at most
  // dtCap while the iterate is within `capDist` of `xc` along `vc`.
  bool descend(std::vector<double>& x, std::string& reason, double dtCap = kInf,
               const std::vector<double>* xc = nullptr, const std::vector<double>* vc = nullptr,
               double capDist = 0) {
    hist.clear(); histDt.clear();
    if (!sys.eval(x, G, J)) { reason = "non-finite residual at the start"; return false; }
    hist.push_back(x); histDt.push_back(0);
    double rate = 1.0;
    if (c.flow) {
      rate = std::numeric_limits<double>::epsilon();
      for (int t = 0; t < nf; ++t) rate = std::max(rate, std::fabs(J[t + m * t]));
    }
    const double dt0 = c.dtInit / rate, dtNewton = 1e8 / rate;
    double dt = std::min(dt0, dtCap);
    std::vector<double> sc, scn, Gn(m), Jn(J.size());
    scaleOf(x, J, sc);
    double rn = rms(G, sc);
    ScaledQR qr;
    std::vector<double> A(J.size()), d, cs(n), xn(n), raw(m), lte;
    auto capped = [&](const std::vector<double>& xv) {
      if (xc == nullptr) return false;
      double s = 0;
      for (int t = 0; t < n; ++t) s += (*vc)[t] * (xv[t] - (*xc)[t]);
      return std::fabs(s) < capDist;
    };
    const int stop = it + c.maxit;
    for (; it < stop; ++it) {
      if (maxsc(G, sc) <= 1) { reason = "converged"; return true; }
      if (capped(x)) dt = std::min(dt, dtCap);
      for (std::size_t q = 0; q < A.size(); ++q) A[q] = -J[q];
      if (c.flow) for (int t = 0; t < nf; ++t) A[t + m * t] += 1.0 / dt;
      for (int t = 0; t < n; ++t) cs[t] = std::fabs(x[t]) + c.atol;
      if (!qr.factor(A, m, n, sc, cs) || !qr.solve(G, d)) { dt /= 10; continue; }
      if (!c.flow) for (double& v : d) v *= dt / (1 + dt);
      bool small = dt >= dtNewton;
      for (int t = 0; small && t < n; ++t)
        small = std::fabs(d[t]) <= c.rtol * std::fabs(x[t]) + c.atol;
      update(x, d, xn);
      if (!sys.eval(xn, Gn, Jn)) { dt /= 10; continue; }
      scaleOf(xn, Jn, scn);
      double rnn = rms(Gn, scn);
      bool progress = rnn <= 0.5 * rn;
      double grow;
      if (c.flow) {
        for (int r = 0; r < m; ++r) raw[r] = r < nf ? 0.5 * (Gn[r] - G[r]) : 0.0;
        if (!qr.solve(raw, lte)) {
          lte.assign(n, 0.0);
          for (int t = 0; t < n; ++t) lte[t] = dt * raw[t];
        }
        double err = 0;
        for (int t = 0; t < n; ++t) {
          double v = lte[t] / (c.flowTol * std::fabs(xn[t]) + c.atol);
          err += v * v;
        }
        err = std::sqrt(err / n);
        double fac = err > 0 ? 0.9 / std::sqrt(err) : 10.0;
        // in the flow range the error decides; a halved residual alone may leave the basin
        if (err > 1 && !(progress && dt >= dtNewton)) { dt *= std::max(0.1, fac); continue; }
        grow = progress ? std::max(2.0, std::min(1e3, rn / std::max(rnn, 1e-300)))
                        : std::min(10.0, std::max(0.2, fac));
      } else {
        if (rnn > 2 * rn && dt > 10 * dt0) { dt /= 4; continue; }
        grow = std::max(2.0, std::min(1e3, rn / std::max(rnn, 1e-300)));
      }
      x = xn; G = Gn; J = Jn; sc = scn;
      hist.push_back(x); histDt.push_back(dt);
      if (small) { ++it; reason = "converged"; return true; }
      dt = std::min(1e20, dt * grow);
      if (capped(x)) dt = std::min(dt, dtCap);
      rn = rnn;
    }
    if (maxsc(G, sc) <= 1) { reason = "converged"; return true; }
    char buf[160];
    std::snprintf(buf, sizeof buf, "no convergence in %d iterations (scaled residual %.2e)",
                  c.maxit, maxsc(G, sc));
    reason = buf;
    return false;
  }

  // Newton steps from a converged x while they contract, until every variable is within
  // rtol (|x| + atol) of the root; each step keeps the residual converged.
  void polish(std::vector<double>& x) {
    std::vector<double> sc, scn, Gn(m), Jn(J.size()), A(J.size()), d, cs(n), xn(n);
    ScaledQR qr;
    double prev = kInf;
    for (int q = 0; q < 20; ++q) {
      scaleOf(x, J, sc);
      for (std::size_t i = 0; i < A.size(); ++i) A[i] = -J[i];
      for (int t = 0; t < n; ++t) cs[t] = std::fabs(x[t]) + c.atol;
      if (!qr.factor(A, m, n, sc, cs) || !qr.solve(G, d)) return;
      double s = 0;
      for (int t = 0; t < n; ++t) s = std::max(s, std::fabs(d[t]) / (c.rtol * cs[t]));
      if (s <= 1 || s > 0.5 * prev) return;
      update(x, d, xn);
      if (!sys.eval(xn, Gn, Jn)) return;
      scaleOf(xn, Jn, scn);
      if (maxsc(Gn, scn) > 1) return;
      x = xn; G = Gn; J = Jn;
      prev = s;
      ++it;
    }
  }

  // The side of x* along v the iterates of the last descent came from, while their steps
  // kept dt lambda < 1.
  double side(const std::vector<double>& xs, const std::vector<double>& v, double lambda) const {
    std::size_t K = 0;
    while (K + 1 < hist.size() && histDt[K + 1] * lambda < 1) ++K;
    double s = 0;
    for (int t = 0; t < n; ++t) s += v[t] * (hist[K][t] - xs[t]);
    return s < 0 ? -1.0 : 1.0;
  }

  // A converged root of the flow that is stable, leaving an unstable one along the leading
  // direction. False with the reason where there is none.
  bool stable(std::vector<double>& x, std::string& reason) {
    Lead L = leading(J, m, n, sys.C, sys.k);
    if (!L.ok) { reason = "eigenvalues failed at the root"; return false; }
    if (L.re <= std::max(1e-8 * L.scale, c.stabilityTol)) return true;
    const std::vector<double> xs = x;
    const double first = side(xs, L.v, L.re);
    double vmax = 0;
    for (double vt : L.v) vmax = std::max(vmax, std::fabs(vt));
    if (vmax == 0) { reason = "unstable root without a direction to leave it"; return false; }
    double eta = kInf;
    for (int t = 0; t < n; ++t) {
      if (L.v[t] == 0) continue;
      const double room = c.positive ? xs[t] : std::max(std::fabs(xs[t]), 1.0);
      if (room > 0) eta = std::min(eta, 0.1 * room / std::fabs(L.v[t]));
    }
    if (!std::isfinite(eta)) eta = 0.1 / vmax;
    double vv = 0;
    for (double vt : L.v) vv += vt * vt;
    char buf[200];
    std::snprintf(buf, sizeof buf, "unstable root, leading eigenvalue %.2e%s, and no stable root along the eigenvector",
                  L.re, L.im > 0 ? " (oscillatory)" : "");
    for (double s : {first, -first}) {
      std::vector<double> x2(n);
      for (int t = 0; t < n; ++t) x2[t] = xs[t] + s * eta * L.v[t];
      if (c.positive) project(x2, sys.C, sys.T, sys.k, n);
      std::string why;
      if (!descend(x2, why, 0.5 / L.re, &xs, &L.v, 10 * eta * vv)) continue;
      polish(x2);
      Lead L2 = leading(J, m, n, sys.C, sys.k);
      if (!L2.ok || L2.re > std::max(1e-8 * L2.scale, c.stabilityTol)) continue;
      x = x2;
      reason = "converged, leaving an unstable root";
      return true;
    }
    reason = buf;
    return false;
  }

  bool solve(std::vector<double>& x, std::string& reason) {
    if (!descend(x, reason)) return false;
    polish(x);
    if (c.flow && c.stability) return stable(x, reason);
    return true;
  }
};

Controls controlsOf(const double* v, int nv) {
  Controls c;
  double* f[] = {&c.rtol, &c.atol, &c.flowTol, &c.dtInit, &c.zmax};
  for (int i = 0; i < 5 && i < nv; ++i) *f[i] = v[i];
  if (nv > 5) c.maxit = static_cast<int>(v[5]);
  if (nv > 6) c.flow = v[6] != 0;
  if (nv > 7) c.positive = v[7] != 0;
  if (nv > 8) c.stability = v[8] != 0;
  if (nv > 9) c.stabilityTol = v[9];
  return c;
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

// The solver for other packages, R_GetCCallable("cppDE", "cppde_ptc_v1"); no R API, safe
// in threads. `fn` is a <model>_eval_ad entry, `x` all nv variables, the solved ones `act`
// (0-based) replaced by the result, `rows` the outputs paired with them; C is k x na
// (column-major) with C x = T. `ctrl` holds rtol, atol, flowTol, dtInit, zmax, maxit, flow,
// positive, stability, stabilityTol; entries past `nctrl` keep their defaults. Returns 1
// when converged, with the iterations and a message of at most `msglen` bytes.
extern "C" int cppde_ptc_v1(eval_ad_fn fn, int nv, int np, int no, double* x, const double* p,
                            int na, const int* act, const int* rows, int k, const double* C,
                            const double* T, const double* ctrl, int nctrl, int* iterations,
                            char* message, int msglen) {
  System sys;
  sys.fn = fn;
  sys.nv = nv; sys.np = np; sys.no = no; sys.na = na; sys.nf = na; sys.k = k;
  sys.xfull.assign(x, x + nv);
  sys.p.assign(p, p + np);
  sys.act.assign(act, act + na);
  sys.rows.assign(rows, rows + na);
  sys.C.assign(C, C + (std::size_t) k * na);
  sys.T.assign(T, T + k);
  sys.dX.assign(static_cast<size_t>(nv) * na, 0.0);
  for (int t = 0; t < na; ++t) sys.dX[act[t] + (std::size_t) nv * t] = 1.0;
  sys.dP.assign(static_cast<size_t>(std::max(np, 1)) * na, 0.0);
  sys.y.assign(no, 0.0);
  sys.dy.assign(static_cast<size_t>(no) * na, 0.0);

  Ptc run(sys, controlsOf(ctrl, nctrl));
  std::vector<double> xa(na);
  for (int t = 0; t < na; ++t) xa[t] = x[act[t]];
  std::string reason;
  const bool ok = run.solve(xa, reason);
  for (int t = 0; t < na; ++t) x[act[t]] = xa[t];
  if (iterations != nullptr) *iterations = run.it;
  if (message != nullptr && msglen > 0) {
    std::strncpy(message, reason.c_str(), msglen - 1);
    message[msglen - 1] = '\0';
  }
  return ok ? 1 : 0;
}

// cppde_ptc_v1() as an external pointer, for packages that reach cppDE from R.
extern "C" SEXP cppde_ptc_entry(void) {
  return R_MakeExternalPtrFn((DL_FUNC) &cppde_ptc_v1, R_NilValue, R_NilValue);
}

// fnS: external pointer to <model>_eval_ad; xS all variables; pS parameters;
// actS, rowsS 0-based indices; CS k x na; TS k; ctrlS named list.
extern "C" SEXP cppde_ptc(SEXP fnS, SEXP nOutS, SEXP xS, SEXP pS, SEXP actS, SEXP rowsS,
                          SEXP CS, SEXP TS, SEXP ctrlS) {
  eval_ad_fn fn = reinterpret_cast<eval_ad_fn>(R_ExternalPtrAddrFn(fnS));
  if (fn == nullptr) Rf_error("ptc: null function pointer");
  const int nv = Rf_length(xS), np = Rf_length(pS), no = Rf_asInteger(nOutS);
  const int na = Rf_length(actS), nf = Rf_length(rowsS), k = Rf_length(TS);
  if (nf != na) Rf_error("ptc: one output row per solved variable");
  const char* names[] = {"rtol", "atol", "flowTol", "dtInit", "zmax", "maxit", "flow", "positive",
                         "stability", "stabilityTol"};
  double ctrl[10];
  for (int i = 0; i < 10; ++i) ctrl[i] = getControl(ctrlS, names[i]);
  std::vector<double> x(REAL(xS), REAL(xS) + nv);
  int iterations = 0;
  char message[256];
  const int ok = cppde_ptc_v1(fn, nv, np, no, x.data(), REAL(pS), na, INTEGER(actS), INTEGER(rowsS), k,
                              REAL(CS), REAL(TS), ctrl, 10, &iterations, message, sizeof message);

  SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
  SEXP nms = PROTECT(Rf_allocVector(STRSXP, 4));
  SEXP xo  = PROTECT(Rf_allocVector(REALSXP, na));
  const int* act = INTEGER(actS);
  for (int t = 0; t < na; ++t) REAL(xo)[t] = x[act[t]];
  SET_VECTOR_ELT(out, 0, xo);
  SET_VECTOR_ELT(out, 1, Rf_ScalarLogical(ok));
  SET_VECTOR_ELT(out, 2, Rf_ScalarInteger(iterations));
  SET_VECTOR_ELT(out, 3, Rf_mkString(message));
  SET_STRING_ELT(nms, 0, Rf_mkChar("x"));
  SET_STRING_ELT(nms, 1, Rf_mkChar("converged"));
  SET_STRING_ELT(nms, 2, Rf_mkChar("iterations"));
  SET_STRING_ELT(nms, 3, Rf_mkChar("message"));
  Rf_setAttrib(out, R_NamesSymbol, nms);
  UNPROTECT(3);
  return out;
}
