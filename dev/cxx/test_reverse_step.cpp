// One integration step backwards on tsit5: stage 3a of dev/adjoint-plan.md.
//
// The oracle is the forward mode on the same step. Forward carries the full
// sensitivity S of (x, theta, t, h) -> x_out, reverse carries w' S for one w.
// Both differentiate the same discrete step, so they agree to rounding and not
// to a solver tolerance.
//
// The step size is an independent, not a constant: the adjoint runs through the
// step-size control, so dy/dh has to come out of the sweep. The forward
// reference only carries h symbolically under CPPDE_SYMBOLIC_STEPSIZE, which is
// why this file defines it; a shipped build never does.
//
// Covered: w' S against S' w for the unit seeds and for a mixed one; the frozen
// path against the same run's x and theta block, with zero time cotangents; the
// replayed step end against the forward run's; chained steps, where h is shared
// and its cotangent sums over them while t chains through t = t0 + s*h.
//
// The step end agrees to rounding, not bit for bit: a vector of double takes the
// BLAS overload of vec_axpy for every stage and that kernel fuses the multiply
// and the add, which a tape of recorded operations cannot.
//
// Not covered here: the control law itself, dt_{k+1} = Ctrl(err_k). That is
// stage 4, and seed_err() is where it attaches.
//
// Every number is printed at %.17g, so the output is the assertion as well.
//
// Build and run:  dev/cxx/run.sh --reverse-step
//
// Copyright (C) 2026 Simon Beyer

// The forward reference must keep h symbolic to be an oracle for dy/dh.
#define CPPDE_SYMBOLIC_STEPSIZE 1

#include <cstdio>
#include <cmath>
#include <string>
#include <utility>
#include <vector>

#include <cppde/cppde.hpp>
#include <cppde/cppde_adjoint_step.hpp>

using cppde::dual;

static int g_failures = 0;

static void check(bool ok, const std::string& what) {
  if (!ok) { std::printf("FAIL  %s\n", what.c_str()); ++g_failures; }
}

static void close(double a, double b, const std::string& what, double tol = 1e-14) {
  const double scale = std::fabs(a) > 1.0 ? std::fabs(a) : 1.0;
  check(std::fabs(a - b) <= tol * scale,
        what + "  forward " + std::to_string(a) + "  reverse " + std::to_string(b));
}

// ---------------------------------------------------------------------------
//  A small nonlinear, non-autonomous model, so no partial is accidentally
//  constant and the cotangent of t is not trivially zero. Parameters are held
//  as the value type, the way the emitted ode_system holds them.
// ---------------------------------------------------------------------------

static constexpr std::size_t NX = 3;
static constexpr std::size_t NP = 4;

template<class V>
struct model {
  std::vector<V> params;

  explicit model(const std::vector<V>& p) : params(p) {}

  void operator()(const std::vector<V>& x, std::vector<V>& dxdt, const V& t) const {
    const V& p0 = params[0];
    const V& p1 = params[1];
    const V& p2 = params[2];
    const V& p3 = params[3];
    const V  s  = cppde::sqrt(x[1]);
    dxdt[0] = -p0 * x[0] * x[1] + p1 * x[2];
    dxdt[1] =  p0 * x[0] * x[1] - p2 * s;
    dxdt[2] =  p2 * s - p1 * x[2] + p3 * cppde::exp(-x[0]) * cppde::cos(t);
  }
};

// The stepper takes a (deriv, jacobian) pair; tsit5 never reads the second.
template<class V>
static std::pair<model<V>, int> make_system(const std::vector<V>& p) {
  return std::make_pair(model<V>(p), 0);
}

// What the code generator emits beside the model for a reverse build. Written
// out here for the same model, so the written step adjoint is checkable
// without a generated source.
struct adjoint_terms {
  std::vector<double> p;

  void jac_t_vec(const std::vector<double>& x, const std::vector<double>& lam,
                 const double& t, std::vector<double>& out) const {
    out.assign(NX, 0.0);
    const double s   = std::sqrt(x[1]);
    const double ds  = 0.5 / s;                 // d sqrt(x1) / d x1
    const double ex  = std::exp(-x[0]) * std::cos(t);
    out[0] = (-p[0] * x[1]) * lam[0] + (p[0] * x[1]) * lam[1]
           + (-p[3] * ex) * lam[2];
    out[1] = (-p[0] * x[0]) * lam[0] + (p[0] * x[0] - p[2] * ds) * lam[1]
           + (p[2] * ds) * lam[2];
    out[2] = (p[1]) * lam[0] + (-p[1]) * lam[2];
  }

  void dfdp_t_vec_axpy(const std::vector<double>& x,
                       const std::vector<double>& lam,
                       const double& t, const double& sc,
                       double* out) const {
    const double s  = std::sqrt(x[1]);
    const double ex = std::exp(-x[0]) * std::cos(t);
    out[NX + 0] += sc*((-x[0] * x[1]) * lam[0] + (x[0] * x[1]) * lam[1]);
    out[NX + 1] += sc*((x[2]) * lam[0] + (-x[2]) * lam[2]);
    out[NX + 2] += sc*((-s) * lam[1] + (s) * lam[2]);
    out[NX + 3] += sc*((ex) * lam[2]);
  }
};

static const double X0[NX] = {1.4, 0.9, 0.3};
static const double P [NP] = {0.7, 0.35, 1.1, 0.25};

static const double T0 = 0.2;
static const double DT = 0.13;

// Direction layout: the states, the parameters, then the trajectory start and
// the step size.
static constexpr unsigned ND = NX + NP;

using D = dual<double, ND>;

// ---------------------------------------------------------------------------
//  Forward reference: n steps under dual, seeded with the identity, so the
//  tangents are the sensitivity matrix. seed_time selects the frozen reference:
//  with t and h unseeded their tangents stay zero and the x/theta block is
//  unchanged, which is what the frozen reverse path has to meet.
// ---------------------------------------------------------------------------

// S is [NX x ND] row-major, x_end is the trajectory end.
static void forward_steps(unsigned n_steps,
                          std::vector<double>& S, std::vector<double>& x_end)
{
  std::vector<D> p(NP);
  for (std::size_t j = 0; j < NP; ++j) { p[j] = D(P[j]); p[j].diff(NX + j); }
  auto sys = make_system<D>(p);

  D t0(T0), h(DT);

  std::vector<D> x(NX), xout(NX), xerr(NX);
  for (std::size_t i = 0; i < NX; ++i) { x[i] = D(X0[i]); x[i].diff(i); }

  cppde::tsit5<D> st;
  // The controller does this; driving do_step directly must too, or the FSAL
  // branch copies k7's value into k1 and leaves its tangents stale.
  st.prepare_sensitivities(ND);
  for (unsigned s = 0; s < n_steps; ++s) {
    const D t = t0 + static_cast<double>(s) * h;
    st.do_step(sys, x, t, xout, h, xerr);
    st.prepare_dense_output();   // acceptance: k7 becomes reusable as the next k1
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
//  Reverse: the adjoint of each step from its checkpoint, newest first, handing
//  the step-start cotangent to the previous step and summing the shared ones.
// ---------------------------------------------------------------------------

using cp_type  = cppde::reverse::step_checkpoint<cppde::tsit5<double>, double>;

// Forward value run, collecting one checkpoint per accepted step.
static void checkpoints(unsigned n_steps, std::vector<cp_type>& cps,
                        std::vector<double>& x_end)
{
  std::vector<double> p(P, P + NP);
  auto sys = make_system<double>(p);

  std::vector<double> x(X0, X0 + NX), xout(NX), xerr(NX);
  cppde::tsit5<double> st;

  cps.assign(n_steps, cp_type());
  for (unsigned s = 0; s < n_steps; ++s) {
    const double t = T0 + static_cast<double>(s) * DT;
    cps[s].capture(st, x, t, DT);
    st.do_step(sys, x, t, xout, DT, xerr);
    st.prepare_dense_output();
    x = xout;
  }
  x_end = x;
}


// ---------------------------------------------------------------------------
//  The same trajectory backwards, written rather than recorded. The step size
//  is frozen here, so this answers the x and theta block and says nothing
//  about t or h.
// ---------------------------------------------------------------------------

static void closed_steps(const std::vector<cp_type>& cps,
                         std::vector<double>& wx, std::vector<double>& wp)
{
  std::vector<double> pv(P, P + NP);
  auto sys = make_system<double>(pv);
  adjoint_terms adj{pv};
  cppde::tsit5<double> st;
  cppde::adjoint::onestep_workspace ws;

  std::vector<double> wphi(NX + NP, 0.0), win(NX, 0.0);
  for (std::size_t s = cps.size(); s-- > 0;) {
    const auto& cp = cps[s];
    cppde::adjoint::apply_onestep_adjoint(
        st, sys, cp.start_state(), cp.t, cp.dt, NX, NX + NP,
        wx.data(), adj, win.data(), wphi.data(), ws);
    wx = win;
  }
  wp.assign(NP, 0.0);
  for (std::size_t j = 0; j < NP; ++j) wp[j] = wphi[NX + j];
}

// ---------------------------------------------------------------------------

static void compare(const char* name, unsigned n_steps, const double* w)
{
  std::vector<double> S, fwd_end;
  forward_steps(n_steps, S, fwd_end);

  std::vector<cp_type> cps;
  std::vector<double>  cp_end;
  checkpoints(n_steps, cps, cp_end);

  std::vector<double> wx(w, w + NX), wp;
  closed_steps(cps, wx, wp);

  std::printf("%-20s", name);
  for (std::size_t i = 0; i < NX; ++i) std::printf(" %.17g", wx[i]);
  std::printf("  |");
  for (std::size_t j = 0; j < NP; ++j) std::printf(" %.17g", wp[j]);
  std::printf("\n");

  // w' S, the contraction the reverse pass computes in one sweep.
  for (unsigned j = 0; j < ND; ++j) {
    double wS = 0.0;
    for (std::size_t i = 0; i < NX; ++i) wS += w[i] * S[i * ND + j];
    const double got = (j < NX) ? wx[j] : wp[j - NX];
    const std::string tag = (j < NX) ? "  dx" + std::to_string(j)
                                     : "  dp" + std::to_string(j - NX);
    close(wS, got, std::string(name) + tag);
  }

  // The step end must survive recomputation under a different scalar type and,
  // beyond one step, without the FSAL carry. A few ulp over a handful of steps:
  // only the double run reaches the BLAS axpy, which fuses where a tape cannot.
  for (std::size_t i = 0; i < NX; ++i) {
    close(fwd_end[i], cp_end[i],
          std::string(name) + " value run x" + std::to_string(i), 1e-15);
  }
}

int main() {
  std::printf("%-20s %-20s %s\n", "case", "w'dx/dx0", "| w'dx/dtheta");

  // Unit seeds are the rows of the Jacobian, one sweep each.
  static const double e0[NX] = {1.0, 0.0, 0.0};
  static const double e1[NX] = {0.0, 1.0, 0.0};
  static const double e2[NX] = {0.0, 0.0, 1.0};
  // What an objective seeds: a reduction over every state at once.
  static const double wm[NX] = {0.6, -1.3, 2.2};

  compare("1 step e0", 1, e0);
  compare("1 step e1", 1, e1);
  compare("1 step e2", 1, e2);
  compare("1 step mixed", 1, wm);

  compare("4 steps mixed", 4, wm);

  compare("2 steps e0", 2, e0);
  compare("2 steps mixed", 2, wm);
  compare("4 steps mixed", 4, wm);

  // A zero seed must produce a zero cotangent.
  static const double z[NX] = {0.0, 0.0, 0.0};
  compare("1 step zero", 1, z);

  std::printf(g_failures == 0 ? "\nOK\n" : "\n%d FAILURES\n", g_failures);
  return g_failures == 0 ? 0 : 1;
}
