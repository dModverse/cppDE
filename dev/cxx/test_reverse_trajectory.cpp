// The trajectory backwards on tsit5: stage 4 of dev/adjoint-plan.md.
//
// The oracle is the forward mode over the same trajectory. Forward carries the
// full sensitivity S of (x0, theta) -> x(t_obs) at every observation time,
// reverse carries the sum over observations of w' S for one w per observation.
// Both differentiate the same discretisation, so they agree to rounding.
//
// Both sides run the same step sequence: the adaptive run in double picks it,
// the checkpoints record it, and the forward reference replays it rather than
// adapting again, which under sensitivities would take different steps. Nothing
// here is compared across two discretisations.
//
// Covered: a seed at one observation and at all of them, on an adaptive run
// through integrate_times_dense, so the checkpoint collector sits in the
// production loop; an observation at the trajectory start, which reaches the
// initial state without passing a step; a zero seed.
//
// Not covered here: the control law, dt_{k+1} = Ctrl(err_k). The trajectory
// recorder reads each step's wt and wdt but does not hand them back, which is
// the frozen path and exactly what the forward sensitivities compute.
//
// Every number is printed at %.17g, so the output is the assertion as well.
//
// Build and run:  dev/cxx/run.sh --reverse-trajectory
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

static void close(double a, double b, const std::string& what, double tol = 1e-13) {
  const double scale = std::fabs(a) > 1.0 ? std::fabs(a) : 1.0;
  check(std::fabs(a - b) <= tol * scale,
        what + "  forward " + std::to_string(a) + "  reverse " + std::to_string(b));
}

// ---------------------------------------------------------------------------
//  A small nonlinear, non-autonomous model, positivity-preserving so the run
//  stays on one smooth branch and the step sequence is set by accuracy alone.
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
    const V  s  = x[1] / (V(1) + x[1]);
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

static const double X0[NX] = {1.4, 0.9, 0.3};
static const double P [NP] = {0.7, 0.35, 1.1, 0.25};

static const double ATOL = 1e-9;
static const double RTOL = 1e-9;
// Deliberately far too large: the controller throws the first attempts away,
// which is the only way the replay of a rejected attempt gets exercised.
static const double HINI = 2.0;

// Set per run by the tolerance sweep at the end.
static double g_atol = ATOL;
static double g_rtol = RTOL;

// Observation times, the first one at the trajectory start.
static const std::vector<double> TIMES = {0.0, 0.37, 1.1, 1.15, 2.0, 3.4, 5.0};

static constexpr unsigned ND = NX + NP;
using D = dual<double, ND>;

using Stepper = cppde::tsit5<double>;
using store_type = cppde::reverse::trajectory_store<Stepper, double>;

// ---------------------------------------------------------------------------
//  The forward value run: the adaptive driver, with the checkpoint collector
//  hooked into it. What comes back is the step sequence and the observed
//  states, and the store holds everything the reverse pass reads.
// ---------------------------------------------------------------------------

static void forward_value(store_type& store, std::vector<double>& x_obs)
{
  std::vector<double> p(P, P + NP);
  auto sys = make_system<double>(p);

  auto controlled = cppde::onestep_controller<Stepper>(g_atol, g_rtol);
  auto dense = cppde::onestep_dense_output<decltype(controlled)>(std::move(controlled));

  std::vector<double> x(X0, X0 + NX);
  std::vector<cppde::detail::FixedEvent<std::vector<double>, double>> fixed;
  std::vector<cppde::detail::RootEvent<std::vector<double>, double>>  root;
  cppde::StepChecker checker(1000000, 1000000);

  x_obs.clear();
  auto obs = [&](const std::vector<double>& xs, const double& t) {
    store.observe(t);
    x_obs.insert(x_obs.end(), xs.begin(), xs.end());
  };
  cppde::reverse::step_collector<decltype(dense), Stepper, double>
      step_obs(store, dense, HINI);

  cppde::integrate_times_dense(dense, sys, x, TIMES.begin(), TIMES.end(), HINI,
                               obs, fixed, root, checker, 1e-8, 1,
                               cppde::detail::no_dt_estimator{}, nullptr,
                               std::ref(step_obs));
}

// ---------------------------------------------------------------------------
//  Forward reference: the recorded step sequence replayed under dual with the
//  identity seed, interpolating at the same observation times. S_obs is
//  [n_obs x NX x ND] row-major.
// ---------------------------------------------------------------------------

static void forward_sens(const store_type& store, std::vector<double>& S_obs,
                         std::vector<double>& x_obs)
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

  S_obs.assign(store.n_obs() * NX * ND, 0.0);
  x_obs.assign(store.n_obs() * NX, 0.0);

  auto record_obs = [&](std::size_t i, const std::vector<D>& xs) {
    for (std::size_t k = 0; k < NX; ++k) {
      x_obs[i * NX + k] = xs[k].x();
      for (unsigned j = 0; j < ND; ++j)
        S_obs[(i * NX + k) * ND + j] = xs[k][j];
    }
  };

  std::size_t next_obs = 0;
  while (next_obs < store.n_obs() && store.obs(next_obs).step == 0) {
    record_obs(next_obs, x);
    ++next_obs;
  }

  std::vector<D> x_interp(NX);
  for (std::size_t k = 0; k < store.n_steps(); ++k) {
    const D t(store.step(k).t), h(store.step(k).dt);
    st.do_step(sys, x, t, xout, h, xerr);
    st.prepare_dense_output();   // acceptance: k7 becomes reusable as the next k1

    while (next_obs < store.n_obs() && store.obs(next_obs).step == k + 1) {
      const D t_end = t + h;
      st.calc_state(D(store.obs(next_obs).t), x_interp, x, t, xout, t_end);
      record_obs(next_obs, x_interp);
      ++next_obs;
    }
    x = xout;
  }
  check(next_obs == store.n_obs(), "every observation reached by the reference");
}

// ---------------------------------------------------------------------------

static std::size_t g_max_nodes = 0;
static double      g_worst = 0.0;   // largest deviation from the dual reference

// One sweep on the recorded grid.
static void sweep_once(const store_type& store, const std::vector<double>& seeds,
                       std::vector<double>& out)
{
  std::vector<double> p(P, P + NP);
  cppde::reverse::trajectory_recorder<Stepper, double> rev;
  rev.sweep(store, p,
            [](const std::vector<codual<double>>& pc) {
              return make_system<codual<double>>(pc);
            },
            seeds);
  if (rev.max_tape_nodes() > g_max_nodes) g_max_nodes = rev.max_tape_nodes();

  out.assign(ND, 0.0);
  for (std::size_t i = 0; i < NX; ++i) out[i] = rev.wx0()[i];
  for (std::size_t j = 0; j < NP; ++j) out[NX + j] = rev.wp()[j];
}

static void compare(const char* name, const store_type& store,
                    const std::vector<double>& S_obs,
                    const std::vector<double>& seeds)
{
  std::vector<double> w;
  sweep_once(store, seeds, w);

  std::printf("%-22s", name);
  for (unsigned j = 0; j < ND; ++j) {
    if (j == NX) std::printf("  |");
    std::printf(" %.17g", w[j]);
  }

  // sum over observations of w' S, the contraction the reverse pass computes.
  double worst = 0.0;
  for (unsigned j = 0; j < ND; ++j) {
    double wS = 0.0;
    for (std::size_t o = 0; o < store.n_obs(); ++o)
      for (std::size_t i = 0; i < NX; ++i)
        wS += seeds[o * NX + i] * S_obs[(o * NX + i) * ND + j];
    const std::string tag = (j < NX) ? "  dx" + std::to_string(j)
                                     : "  dp" + std::to_string(j - NX);
    close(wS, w[j], std::string(name) + tag);

    const double scale = std::fabs(wS) > 1.0 ? std::fabs(wS) : 1.0;
    const double rel = std::fabs(w[j] - wS) / scale;
    if (rel > worst) worst = rel;
  }
  if (worst > g_worst) g_worst = worst;
  std::printf("   dev %.2e\n", worst);
}

int main() {
  store_type store;
  std::vector<double> x_run;
  forward_value(store, x_run);

  std::vector<double> S_obs, x_ref;
  forward_sens(store, S_obs, x_ref);

  std::printf("steps %zu  observations %zu\n", store.n_steps(), store.n_obs());
  // HINI is far too large on purpose, so the controller throws attempts away.
  // Nothing records how many: a rejected attempt is a control decision, and
  // only the accepted grid reaches the store.
  check(store.n_steps() > 4, "the controller took more than a handful of steps");
  check(store.n_obs() == TIMES.size(), "one observation per requested time");
  check(store.obs(0).step == 0, "the first observation is the trajectory start");

  // The replayed sequence must reproduce the driver's own output, or the two
  // sides are not differentiating the same trajectory.
  for (std::size_t k = 0; k < x_run.size() && k < x_ref.size(); ++k)
    close(x_run[k], x_ref[k], "replayed observation " + std::to_string(k), 1e-13);

  const std::size_t n_seed = store.n_obs() * NX;

  // One observation at a time: every column of the trajectory Jacobian.
  for (std::size_t o = 0; o < store.n_obs(); ++o) {
    for (std::size_t i = 0; i < NX; ++i) {
      std::vector<double> seeds(n_seed, 0.0);
      seeds[o * NX + i] = 1.0;
      const std::string nm = "obs" + std::to_string(o) + " e" + std::to_string(i);
      compare(nm.c_str(), store, S_obs, seeds);
    }
  }

  // What an objective seeds: every observation at once, one sweep.
  std::vector<double> all(n_seed);
  for (std::size_t k = 0; k < n_seed; ++k)
    all[k] = 0.3 * static_cast<double>(k % 5) - 0.7;
  compare("all observations", store, S_obs, all);

  std::vector<double> zero(n_seed, 0.0);
  compare("zero", store, S_obs, zero);

  // The point of a checkpoint per step: the tape holds one step, never the
  // trajectory. What grows with the run is the checkpoint store, and for a
  // one-step method that is the state plus two scalars per step.
  const std::size_t node = sizeof(cppde::codual_tape<double>::node);
  std::printf("\ntape %zu nodes at the widest step, %zu bytes; taping the whole "
              "run would be %zu\ncheckpoints %zu bytes over %zu steps\n",
              g_max_nodes, g_max_nodes * node,
              g_max_nodes * store.n_steps() * node,
              store.n_steps() * (NX + 2) * sizeof(double), store.n_steps());

  // The two modes share a grid, so their disagreement is roundoff and not
  // discretisation: it must not grow with the tolerance the way a
  // discretisation term would. Over six decades the step sequence changes
  // completely, which is what makes this more than a repeat of the run above.
  std::printf("\ndeviation from the dual reference at most %.2e at rtol %.0e\n",
              g_worst, RTOL);
  std::printf("\n%-10s %-8s %-12s\n", "rtol", "steps", "deviation");
  double worst_all = 0.0;
  for (const double r : {1e-6, 1e-8, 1e-10, 1e-12}) {
    g_atol = g_rtol = r;
    store_type st2;
    std::vector<double> run2;
    forward_value(st2, run2);

    std::vector<double> S2, x2;
    forward_sens(st2, S2, x2);

    std::vector<double> all2(st2.n_obs() * NX);
    for (std::size_t k = 0; k < all2.size(); ++k)
      all2[k] = 0.3 * static_cast<double>(k % 5) - 0.7;

    std::vector<double> w2;
    sweep_once(st2, all2, w2);

    double worst = 0.0;
    for (unsigned j = 0; j < ND; ++j) {
      double wS = 0.0;
      for (std::size_t o = 0; o < st2.n_obs(); ++o)
        for (std::size_t i = 0; i < NX; ++i)
          wS += all2[o * NX + i] * S2[(o * NX + i) * ND + j];
      const double sc = std::fabs(wS) > 1.0 ? std::fabs(wS) : 1.0;
      const double rel = std::fabs(w2[j] - wS) / sc;
      if (rel > worst) worst = rel;
    }
    std::printf("%-10.0e %-8zu %-12.3e\n", r, st2.n_steps(), worst);
    if (worst > worst_all) worst_all = worst;
  }
  g_atol = ATOL; g_rtol = RTOL;
  check(worst_all < 1e-9,
        "reverse and dual agree on every step sequence, not just one");

  std::printf(g_failures == 0 ? "\nOK\n" : "\n%d FAILURES\n", g_failures);
  return g_failures == 0 ? 0 : 1;
}
