// One integration step backwards on tsit5: stage 3a of dev/adjoint-plan.md.
//
// The oracle is the forward mode on the same step. Forward carries the full
// sensitivity S of (x, theta) -> x_out, reverse carries w' S for one w. Both
// differentiate the same discrete step, so they agree to rounding and not to a
// solver tolerance.
//
// Covered: w' S against S' w for the unit seeds and for a mixed one; the
// replayed step end against the forward one, bit for bit; chained steps, where
// the forward run recycles k7 as the next k1 and the replay recomputes it.
//
// Every number is printed at %.17g, so the output is the assertion as well.
//
// Build and run:  dev/cxx/run.sh --reverse-step
//
// Copyright (C) 2026 Simon Beyer

#include <cstdio>
#include <cmath>
#include <string>
#include <utility>
#include <vector>

#include <cppde/cppde.hpp>

using cppde::codual;
using cppde::dual;

static int g_failures = 0;

static void check(bool ok, const std::string& what) {
  if (!ok) { std::printf("FAIL  %s\n", what.c_str()); ++g_failures; }
}

static void close(double a, double b, const std::string& what) {
  const double scale = std::fabs(a) > 1.0 ? std::fabs(a) : 1.0;
  check(std::fabs(a - b) <= 1e-14 * scale,
        what + "  forward " + std::to_string(a) + "  reverse " + std::to_string(b));
}

// ---------------------------------------------------------------------------
//  A small nonlinear model, so no partial is accidentally constant. Parameters
//  are held as the value type, the way the emitted ode_system holds them.
// ---------------------------------------------------------------------------

static constexpr std::size_t NX = 3;
static constexpr std::size_t NP = 4;

template<class V>
struct model {
  std::vector<V> params;

  explicit model(const std::vector<V>& p) : params(p) {}

  void operator()(const std::vector<V>& x, std::vector<V>& dxdt, const V& /*t*/) const {
    const V& p0 = params[0];
    const V& p1 = params[1];
    const V& p2 = params[2];
    const V& p3 = params[3];
    const V  s  = cppde::sqrt(x[1]);
    dxdt[0] = -p0 * x[0] * x[1] + p1 * x[2];
    dxdt[1] =  p0 * x[0] * x[1] - p2 * s;
    dxdt[2] =  p2 * s - p1 * x[2] + p3 * cppde::exp(-x[0]);
  }
};

// The stepper takes a (deriv, jacobian) pair; tsit5 never reads the second.
template<class V>
static std::pair<model<V>, int> make_system(const std::vector<V>& p) {
  return std::make_pair(model<V>(p), 0);
}

static const double X0[NX] = {1.4, 0.9, 0.3};
static const double P [NP] = {0.7, 0.35, 1.1, 0.25};

static const double T0 = 0.2;
static const double DT = 0.13;

// ---------------------------------------------------------------------------
//  Forward reference: n steps under dual, seeded with the identity over
//  (x, theta), so the tangents are the sensitivity matrix.
// ---------------------------------------------------------------------------

static constexpr unsigned ND = NX + NP;
using D = dual<double, ND>;

// S is [NX x ND] row-major, x_end is the step end.
static void forward_steps(unsigned n_steps, std::vector<double>& S,
                          std::vector<double>& x_end)
{
  std::vector<D> p(NP);
  for (std::size_t j = 0; j < NP; ++j) { p[j] = D(P[j]); p[j].diff(NX + j); }
  auto sys = make_system<D>(p);

  std::vector<D> x(NX), xout(NX), xerr(NX);
  for (std::size_t i = 0; i < NX; ++i) { x[i] = D(X0[i]); x[i].diff(i); }

  cppde::tsit5<D> st;
  // The controller does this; driving do_step directly must too, or the FSAL
  // branch copies k7's value into k1 and leaves its tangents stale.
  st.prepare_sensitivities(ND);
  double t = T0;
  for (unsigned s = 0; s < n_steps; ++s) {
    st.do_step(sys, x, t, xout, DT, xerr);
    st.prepare_dense_output();   // acceptance: k7 becomes reusable as the next k1
    x = xout;
    t += DT;
  }

  S.assign(NX * ND, 0.0);
  x_end.assign(NX, 0.0);
  for (std::size_t i = 0; i < NX; ++i) {
    x_end[i] = x[i].x();
    for (unsigned j = 0; j < ND; ++j) S[i * ND + j] = x[i][j];
  }
}

// ---------------------------------------------------------------------------
//  Reverse: replay each step from its checkpoint, newest first, handing the
//  step-start cotangent to the previous step and summing the parameter one.
// ---------------------------------------------------------------------------

using rec_type = cppde::reverse::step_recorder<cppde::tsit5<double>, double>;
using cp_type  = cppde::reverse::step_checkpoint<cppde::tsit5<double>, double>;

// Forward value run, collecting one checkpoint per accepted step.
static void checkpoints(unsigned n_steps, std::vector<cp_type>& cps,
                        std::vector<double>& x_end)
{
  std::vector<double> p(P, P + NP);
  auto sys = make_system<double>(p);

  std::vector<double> x(X0, X0 + NX), xout(NX), xerr(NX);
  cppde::tsit5<double> st;
  double t = T0;

  cps.assign(n_steps, cp_type());
  for (unsigned s = 0; s < n_steps; ++s) {
    cps[s].capture(st, x, t, DT);
    st.do_step(sys, x, t, xout, DT, xerr);
    st.prepare_dense_output();
    x = xout;
    t += DT;
  }
  x_end = x;
}

// wx on entry is the cotangent of the trajectory end; on return that of its
// start. wp accumulates over the steps, as theta is shared by all of them.
static void reverse_steps(const std::vector<cp_type>& cps,
                          std::vector<double>& wx, std::vector<double>& wp,
                          std::vector<double>& replayed_end)
{
  wp.assign(NP, 0.0);
  for (std::size_t s = cps.size(); s-- > 0;) {
    rec_type rec;
    rec.begin();

    std::vector<codual<double>> p(NP);
    for (std::size_t j = 0; j < NP; ++j) p[j] = codual<double>(P[j]);
    rec.independent(p);
    // The functor copies params; that copy carries the slots to read back.
    auto sys = make_system<codual<double>>(p);
    rec.record(sys, cps[s]);

    if (s + 1 == cps.size()) {
      replayed_end.assign(NX, 0.0);
      for (std::size_t i = 0; i < NX; ++i) replayed_end[i] = rec.xout()[i].x();
    }

    rec.seed(wx);
    rec.sweep();
    rec.accumulate(sys.first.params, wp);
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

  std::printf("%-18s", name);
  for (std::size_t i = 0; i < NX; ++i) std::printf(" %.17g", wx[i]);
  std::printf("  |");
  for (std::size_t j = 0; j < NP; ++j) std::printf(" %.17g", wp[j]);
  std::printf("\n");

  // w' S, the contraction the reverse pass computes in one sweep.
  for (unsigned j = 0; j < ND; ++j) {
    double wS = 0.0;
    for (std::size_t i = 0; i < NX; ++i) wS += w[i] * S[i * ND + j];
    const double got = (j < NX) ? wx[j] : wp[j - NX];
    close(wS, got, std::string(name) + (j < NX ? "  dx" : "  dp") +
                   std::to_string(j < NX ? j : j - NX));
  }

  // The step end must survive recomputation under a different scalar type and,
  // beyond one step, without the FSAL carry.
  for (std::size_t i = 0; i < NX; ++i) {
    check(fwd_end[i] == cp_end[i],
          std::string(name) + " value run x" + std::to_string(i));
    check(cp_end[i] == replayed_end[i],
          std::string(name) + " replay x" + std::to_string(i));
  }
}

int main() {
  std::printf("%-18s %-20s %-20s %-20s\n", "case", "w'dx/dx0[0]", "[1]", "[2]");

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

  // Two steps: the second hands its start cotangent to the first.
  compare("2 steps e0", 2, e0);
  compare("2 steps mixed", 2, wm);

  // Four steps, to show the chain does not drift.
  compare("4 steps mixed", 4, wm);

  // A zero seed must produce a zero cotangent.
  static const double z[NX] = {0.0, 0.0, 0.0};
  compare("1 step zero", 1, z);

  std::printf(g_failures == 0 ? "\nOK\n" : "\n%d FAILURES\n", g_failures);
  return g_failures == 0 ? 0 : 1;
}
