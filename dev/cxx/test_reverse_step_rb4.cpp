// One Rosenbrock step backwards: stage 3b of dev/adjoint-plan.md.
//
// The oracle is the forward mode on the same step. Forward carries the full
// sensitivity S of (x, theta) -> x_out, reverse carries w' S for one w. Both
// differentiate the same discrete step, so they agree to rounding.
//
// Six linear solves against one shared W = -J + I/(gamma*dt), and the reverse
// pass never solves them the way the forward does. Each becomes the equation
// W g = rhs with g recovered from the same factorisation the sweep uses
// transposed, so the checkpoint holds nothing but the step start.
//
// Covered: the unit seeds and a mixed one, over one, two and four chained steps,
// against a forward reference driven on the same step sizes.
//
// Every number is printed at %.17g, so the output is the assertion as well.
//
// Build and run:  dev/cxx/run.sh --reverse-step-rb4
//
// Copyright (C) 2026 Simon Beyer

#include <cstdio>
#include <cmath>
#include <string>
#include <utility>
#include <vector>

#include <cppde/cppde.hpp>
#include <cppde/cppde_adjoint_step.hpp>

using cppde::codual;
using cppde::dual;

static int g_failures = 0;

static void check(bool ok, const std::string& what) {
  if (!ok) { std::printf("FAIL  %s\n", what.c_str()); ++g_failures; }
}

static void close(double a, double b, const std::string& what, double tol = 1e-11) {
  const double scale = std::fabs(a) > 1.0 ? std::fabs(a) : 1.0;
  char buf[160];
  std::snprintf(buf, sizeof buf, "  forward %.17g  reverse %.17g  rel %.2e",
                a, b, std::fabs(a - b) / scale);
  check(std::fabs(a - b) <= tol * scale, what + buf);
}

// ---------------------------------------------------------------------------
//  A small nonlinear, non-autonomous model, positivity-preserving so the run
//  stays on one smooth branch, with a hand-written Jacobian beside it.
// ---------------------------------------------------------------------------

static constexpr std::size_t NX = 3;
static constexpr std::size_t NP = 3;

template<class V>
struct model {
  std::vector<V> p;
  void operator()(const std::vector<V>& x, std::vector<V>& d, const V& t) const {
    d[0] = -p[0] * x[0] + p[1] * x[1] * x[2];
    d[1] =  p[0] * x[0] - p[1] * x[1] * x[2] - p[2] * x[1] * x[1];
    d[2] =  p[2] * x[1] * x[1] * cppde::cos(t);
  }
};

// Minus the Jacobian, which is what the iteration matrix is built from:
// factorize_W only adds the diagonal, W = -J + inv_gamma_dt * I.
template<class V>
struct jacobian {
  std::vector<V> p;
  void operator()(const std::vector<V>& x, cppde::dense_matrix<V>& J, const V& t,
                  std::vector<V>& dfdt) const {
    J(0,0) =  p[0];   J(0,1) = -p[1] * x[2];  J(0,2) = -p[1] * x[1];
    J(1,0) = -p[0];   J(1,1) =  p[1] * x[2] + V(2.0) * p[2] * x[1];
                      J(1,2) =  p[1] * x[1];
    J(2,0) =  V(0.0); J(2,1) = -V(2.0) * p[2] * x[1] * cppde::cos(t);
                      J(2,2) =  V(0.0);
    dfdt[0] = V(0.0);
    dfdt[1] = V(0.0);
    dfdt[2] = -p[2] * x[1] * x[1] * cppde::sin(t);
  }
};

template<class V>
static std::pair<model<V>, jacobian<V>> make_system(const std::vector<V>& p) {
  return std::make_pair(model<V>{p}, jacobian<V>{p});
}

static const double X0[NX] = {1.0, 0.35, 0.15};
static const double P [NP] = {0.9, 1.4, 0.6};

static const double T0 = 0.1;
static const double DT = 0.05;

static constexpr unsigned ND = NX + NP;
using D = dual<double, ND>;

using StepD = cppde::rosenbrock4<double>;
using cp_type = cppde::reverse::step_checkpoint<StepD, double>;

// ---------------------------------------------------------------------------
//  Forward reference: n steps under dual, seeded with the identity, so the
//  tangents are the sensitivity matrix.
// ---------------------------------------------------------------------------

static void forward_steps(unsigned n_steps, std::vector<double>& S,
                          std::vector<double>& x_end)
{
  std::vector<D> p(NP);
  for (std::size_t j = 0; j < NP; ++j) { p[j] = D(P[j]); p[j].diff(NX + j); }
  auto sys = make_system<D>(p);

  std::vector<D> x(NX), xout(NX), xerr(NX);
  for (std::size_t i = 0; i < NX; ++i) { x[i] = D(X0[i]); x[i].diff(i); }

  cppde::rosenbrock4<D> st;
  st.prepare_sensitivities(ND);
  for (unsigned s = 0; s < n_steps; ++s) {
    st.do_step(sys, x, D(T0 + s * DT), xout, D(DT), xerr);
    x = xout;
  }

  S.assign(NX * ND, 0.0);
  x_end.assign(NX, 0.0);
  for (std::size_t i = 0; i < NX; ++i) {
    x_end[i] = x[i].x();
    for (unsigned j = 0; j < ND; ++j) S[i * ND + j] = x[i][j];
  }
}

// ---------------------------------------------------------------------------
//  The forward value run, one checkpoint per step.
// ---------------------------------------------------------------------------

static void checkpoints(unsigned n_steps, std::vector<cp_type>& cps,
                        std::vector<double>& x_end)
{
  std::vector<double> p(P, P + NP);
  auto sys = make_system<double>(p);

  std::vector<double> x(X0, X0 + NX), xout(NX), xerr(NX);
  cppde::rosenbrock4<double> st;

  cps.assign(n_steps, cp_type());
  for (unsigned s = 0; s < n_steps; ++s) {
    const double t = T0 + s * DT;
    cps[s].capture(st, x, t, DT);
    st.do_step(sys, x, t, xout, DT, xerr);
    x = xout;
  }
  x_end = x;
}

// What the code generator emits beside the model for a reverse rb4 build.
// Written out here for the same model, so the written step adjoint is checkable
// without a generated source. A Rosenbrock stage solves against a matrix built
// from the Jacobian, so two of these four are second derivatives.
struct adjoint_terms {
  std::vector<double> p;

  void jac_t_vec(const std::vector<double>& x, const std::vector<double>& lam,
                 const double& t, std::vector<double>& out) const {
    out.assign(NX, 0.0);
    out[0] = (-p[0]) * lam[0] + (p[0]) * lam[1];
    out[1] = (p[1] * x[2]) * lam[0]
           + (-p[1] * x[2] - 2.0 * p[2] * x[1]) * lam[1]
           + (2.0 * p[2] * x[1] * std::cos(t)) * lam[2];
    out[2] = (p[1] * x[1]) * lam[0] + (-p[1] * x[1]) * lam[1];
  }

  void dfdp_t_vec_axpy(const std::vector<double>& x,
                       const std::vector<double>& lam,
                       const double& t, const double& sc,
                       double* out) const {
    out[NX + 0] += sc * ((-x[0]) * lam[0] + (x[0]) * lam[1]);
    out[NX + 1] += sc * ((x[1] * x[2]) * lam[0] + (-x[1] * x[2]) * lam[1]);
    out[NX + 2] += sc * ((-x[1] * x[1]) * lam[1]
                         + (x[1] * x[1] * std::cos(t)) * lam[2]);
  }

  void jvp_x_t_vec(const std::vector<double>& x, const std::vector<double>& v,
                   const std::vector<double>& lam, const double& t,
                   std::vector<double>& out) const {
    out.assign(NX, 0.0);
    out[0] = 0.0;
    out[1] = (p[1] * v[2]) * lam[0]
           + (-2.0 * p[2] * v[1] - p[1] * v[2]) * lam[1]
           + (2.0 * p[2] * std::cos(t) * v[1]) * lam[2];
    out[2] = (p[1] * v[1]) * lam[0] + (-p[1] * v[1]) * lam[1];
  }

  void jvp_p_t_vec_axpy(const std::vector<double>& x,
                        const std::vector<double>& v,
                        const std::vector<double>& lam, const double& t,
                        const double& sc, double* out) const {
    const double q = x[2] * v[1] + x[1] * v[2];
    out[NX + 0] += sc * ((-v[0]) * lam[0] + (v[0]) * lam[1]);
    out[NX + 1] += sc * ((q) * lam[0] + (-q) * lam[1]);
    out[NX + 2] += sc * ((-2.0 * x[1] * v[1]) * lam[1]
                         + (2.0 * x[1] * std::cos(t) * v[1]) * lam[2]);
  }

  void dfdt_x_t_vec(const std::vector<double>& x, const std::vector<double>& lam,
                    const double& t, std::vector<double>& out) const {
    out.assign(NX, 0.0);
    out[1] = (-2.0 * p[2] * x[1] * std::sin(t)) * lam[2];
  }

  void dfdt_p_t_vec_axpy(const std::vector<double>& x,
                         const std::vector<double>& lam,
                         const double& t, const double& sc,
                         double* out) const {
    out[NX + 2] += sc * ((-x[1] * x[1] * std::sin(t)) * lam[2]);
  }
};

// ---------------------------------------------------------------------------
//  The same steps, written rather than recorded.
// ---------------------------------------------------------------------------

static void closed_steps(const std::vector<cp_type>& cps,
                         std::vector<double>& wx, std::vector<double>& wp)
{
  std::vector<double> pv(P, P + NP);
  auto sysd = make_system<double>(pv);
  adjoint_terms adj{pv};
  cppde::adjoint::rosenbrock_workspace ws;
  StepD st;

  std::vector<double> wphi(NX + NP, 0.0), w_in(NX, 0.0);
  for (std::size_t s = cps.size(); s-- > 0;) {
    w_in.assign(NX, 0.0);
    cppde::adjoint::apply_rosenbrock_adjoint(
        st, sysd, cps[s].x.data(), cps[s].t, cps[s].dt, NX, NX + NP,
        wx.data(), adj, w_in.data(), wphi.data(), ws);
    wx = w_in;
  }
  wp.assign(NP, 0.0);
  for (std::size_t j = 0; j < NP; ++j) wp[j] = wphi[NX + j];
}

static void reverse_steps(const std::vector<cp_type>& cps,
                          std::vector<double>& wx, std::vector<double>& wp,
                          std::vector<double>& replayed_end)
{
  using C = codual<double>;
  std::vector<double> pv(P, P + NP);
  auto jac_d = jacobian<double>{pv};

  wp.assign(NP, 0.0);
  for (std::size_t s = cps.size(); s-- > 0;) {
    cppde::reverse::step_recorder<StepD, double> rec;
    rec.begin();

    std::vector<C> p(NP);
    for (std::size_t j = 0; j < NP; ++j) p[j] = C(P[j]);
    rec.independent(p);
    auto sys = make_system<C>(p);

    rec.load(cps[s], cps[s].dt);

    // W at the step start, which is where rosenbrock4 evaluates it.
    cppde::reverse::equation_solver<jacobian<double>, double> solver(jac_d);
    solver.prepare(cps[s].x, cps[s].t,
                   rec.stepper().replay_inv_gamma_dt(cps[s].dt));

    rec.attempt_staged(sys, rec.dt_in(), solver);

    if (s + 1 == cps.size()) {
      replayed_end.assign(NX, 0.0);
      for (std::size_t i = 0; i < NX; ++i) replayed_end[i] = rec.xout()[i].x();
    }

    rec.seed(wx);
    rec.sweep();
    rec.accumulate(p, wp);
    wx = rec.wx();
  }
}

// ---------------------------------------------------------------------------

static void compare(const char* name, unsigned n_steps, const double* w)
{
  std::vector<double> S, fwd_end;
  forward_steps(n_steps, S, fwd_end);

  std::vector<cp_type> cps;
  std::vector<double>  cp_end;
  checkpoints(n_steps, cps, cp_end);

  std::vector<double> wx(w, w + NX), wp, replayed_end;
  reverse_steps(cps, wx, wp, replayed_end);

  std::vector<double> cwx(w, w + NX), cwp;
  closed_steps(cps, cwx, cwp);

  std::printf("%-20s", name);
  for (std::size_t i = 0; i < NX; ++i) std::printf(" %.17g", wx[i]);
  std::printf("  |");
  for (std::size_t j = 0; j < NP; ++j) std::printf(" %.17g", wp[j]);
  std::printf("\n");

  for (unsigned j = 0; j < ND; ++j) {
    double wS = 0.0;
    for (std::size_t i = 0; i < NX; ++i) wS += w[i] * S[i * ND + j];
    const double got = (j < NX) ? wx[j] : wp[j - NX];
    const std::string tag = (j < NX) ? "  dx" + std::to_string(j)
                                     : "  dp" + std::to_string(j - NX);
    close(wS, got, std::string(name) + tag);
    const double gotc = (j < NX) ? cwx[j] : cwp[j - NX];
    close(wS, gotc, std::string(name) + " written" + tag);
  }

  // The stage values come back out of the factorisation rather than a
  // checkpoint, so the replayed step end is the check that they came back right.
  for (std::size_t i = 0; i < NX; ++i) {
    close(fwd_end[i], cp_end[i], std::string(name) + " value run x" + std::to_string(i));
    close(cp_end[i], replayed_end[i], std::string(name) + " replay x" + std::to_string(i));
  }
}

int main() {
  std::printf("%-20s %-20s %s\n", "case", "w'dx/dx0", "| w'dx/dtheta");

  static const double e0[NX] = {1.0, 0.0, 0.0};
  static const double e1[NX] = {0.0, 1.0, 0.0};
  static const double e2[NX] = {0.0, 0.0, 1.0};
  static const double wm[NX] = {0.6, -1.3, 2.2};

  compare("1 step e0", 1, e0);
  compare("1 step e1", 1, e1);
  compare("1 step e2", 1, e2);
  compare("1 step mixed", 1, wm);
  compare("2 steps mixed", 2, wm);
  compare("4 steps mixed", 4, wm);

  static const double z[NX] = {0.0, 0.0, 0.0};
  compare("1 step zero", 1, z);

  std::printf(g_failures == 0 ? "\nOK\n" : "\n%d FAILURES\n", g_failures);
  return g_failures == 0 ? 0 : 1;
}
