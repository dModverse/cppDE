// The trajectory backwards on all four methods: stage 4 of dev/adjoint-plan.md.
//
// The oracle is the forward mode over the same trajectory, and the sharpness
// comes from both sides running the same step sequence. The adaptive run in
// double picks it and the checkpoints record it, decisions included: the step
// times and sizes, and for a multistep method the order the controller chose and
// the rescale it applied. The reference then replays that sequence under dual
// rather than adapting again, which under sensitivities would take other steps.
//
// What the seed lands on is the trajectory start: the initial state and the
// parameters, for every method. A multistep run builds its Nordsieck history
// out of that state through initialize(), so the reference builds it the same
// way and the reverse sweep collapses the history's cotangent back onto the
// state. That is the map dMod2 asks for, and the map the two sides
// differentiate here.
//
// Covered: bdf, adams, rb4 and tsit5, each on an adaptive run with observations
// interpolated inside the steps, seeded one carry slot at a time and with a
// mixed seed over every observation at once.
//
// Every number is printed at %.17g, so the output is the assertion as well.
//
// Build and run:  dev/cxx/run.sh --reverse-trajectory-methods
//
// Copyright (C) 2026 Simon Beyer

#include <cstdio>
#include <cmath>
#include <functional>
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

static void close(double a, double b, const std::string& what, double tol) {
  const double scale = std::fabs(a) > 1.0 ? std::fabs(a) : 1.0;
  char buf[160];
  std::snprintf(buf, sizeof buf, "  forward %.17g  reverse %.17g  rel %.2e",
                a, b, std::fabs(a - b) / scale);
  check(std::fabs(a - b) <= tol * scale, what + buf);
}

// ---------------------------------------------------------------------------
//  A small nonlinear, non-autonomous model, positivity-preserving so the run
//  stays on one smooth branch, with minus its Jacobian beside it, which is what
//  an iteration matrix is built from.
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
static double P [NP] = {0.9, 1.4, 0.6};

static const double ATOL = 1e-8;
static const double RTOL = 1e-8;
static const double HINI = 0.02;

static const std::vector<double> TIMES = {0.1, 0.37, 0.8, 0.85, 1.4, 2.0};

// Directions: the widest carry any method here can produce, which for Adams is
// its full order range, plus the parameters.
static constexpr unsigned ND_MAX = 13 * NX + NP;
using D = dual<double, ND_MAX>;

// ---------------------------------------------------------------------------
//  The controller and dense-output wrapper each stepper family is driven with.
// ---------------------------------------------------------------------------

template<class S> struct pipeline;

template<cppde::multistep_method M, class V, class J, class R>
struct pipeline<cppde::multistepper<M, V, J, R>> {
  using controller = cppde::multistepper_controller<cppde::multistepper<M, V, J, R>>;
  using dense      = cppde::multistepper_dense_output<controller>;
};
template<class V, class R> struct pipeline<cppde::tsit5<V, R>> {
  using controller = cppde::onestep_controller<cppde::tsit5<V, R>>;
  using dense      = cppde::onestep_dense_output<controller>;
};
template<class V, class R> struct pipeline<cppde::rosenbrock4<V, R>> {
  using controller = cppde::onestep_controller<cppde::rosenbrock4<V, R>>;
  using dense      = cppde::onestep_dense_output<controller>;
};

template<class S> using store_of = cppde::reverse::trajectory_store<S, double>;

template<class S, class = void> struct is_adams : std::false_type {};
template<class S> struct is_adams<S, std::void_t<decltype(S::method)>>
  : std::bool_constant<S::method == cppde::multistep_method::adams> {};

// ---------------------------------------------------------------------------
//  The forward value run: the adaptive driver with the checkpoint collector
//  hooked into it, so the store holds the step sequence and its decisions.
// ---------------------------------------------------------------------------

template<class S>
static void forward_value(store_of<S>& store, std::vector<double>& x_obs)
{
  std::vector<double> p(P, P + NP);
  auto sys = make_system<double>(p);

  typename pipeline<S>::controller ctl(ATOL, RTOL);
  typename pipeline<S>::dense dense(std::move(ctl));

  std::vector<double> x(X0, X0 + NX);
  std::vector<cppde::detail::FixedEvent<std::vector<double>, double>> fixed;
  std::vector<cppde::detail::RootEvent<std::vector<double>, double>>  root;
  cppde::StepChecker checker(1000000, 1000000);

  x_obs.clear();
  auto obs = [&](const std::vector<double>& xs, const double& t) {
    store.observe(t);
    x_obs.insert(x_obs.end(), xs.begin(), xs.end());
  };
  cppde::reverse::step_collector<decltype(dense), S, double> step_obs(store, dense, HINI);

  cppde::integrate_times_dense(dense, sys, x, TIMES.begin(), TIMES.end(), HINI,
                               obs, fixed, root, checker, 1e-8, 1,
                               cppde::detail::no_dt_estimator{}, nullptr,
                               std::ref(step_obs));
}

// ---------------------------------------------------------------------------
//  Forward reference: the recorded sequence replayed under dual, seeded with the
//  identity over the trajectory start's carry and the parameters. Only the
//  decisions come from the checkpoints; the carry itself evolves under dual.
// ---------------------------------------------------------------------------

template<class S>
static void forward_sens(const store_of<S>& store, std::size_t n_carry,
                         unsigned nd, std::vector<double>& S_obs,
                         std::vector<double>& x_ref)
{
  using rstep = typename S::template rebind_value<D>;
  constexpr bool multistep = cppde::reverse::has_step_snapshot<S>::value;

  // Each direction is seeded at the magnitude of what it perturbs. The
  // corrector's error norm takes the maximum over every sensitivity direction,
  // and a unit tangent on a high Nordsieck slot is orders above the slot itself,
  // which stops the iteration early rather than sharpening it. The map is linear
  // in the seed, so S is unscaled again below.
  std::vector<double> scale(nd, 1.0);
  for (std::size_t j = 0; j < NP; ++j)
    scale[n_carry + j] = std::max(std::fabs(P[j]), 1e-2);

  std::vector<D> p(NP);
  for (std::size_t j = 0; j < NP; ++j) {
    p[j] = D(P[j]);
    p[j].diff(static_cast<unsigned>(n_carry + j));
    p[j].d(static_cast<unsigned>(n_carry + j)) = scale[n_carry + j];
  }
  auto sys = make_system<D>(p);

  rstep st;
  if constexpr (multistep) {
    st.set_tolerances(ATOL, RTOL);
    // The corrector's own stopping rule takes the maximum over every
    // sensitivity direction, so under the identity seed it stops early and at
    // high order not at all. Both modes are defined to differentiate the exactly
    // solved equation, so the reference is given the room to solve it.
    st.set_max_corrector_iters(200);
  }
  st.prepare_sensitivities(nd);

  const auto& cp0 = store.step(0);
  std::vector<D> x(NX), xout(NX), xerr(NX), x_interp(NX);

  for (std::size_t i = 0; i < NX; ++i) {
    scale[i] = std::max(std::fabs(cp0.start_state()[i]), 1e-2);
    x[i] = D(cp0.start_state()[i]);
    x[i].diff(static_cast<unsigned>(i));
    x[i].d(static_cast<unsigned>(i)) = scale[i];
  }
  if constexpr (multistep) {
    // The history the run started from, built the way the run built it.
    std::vector<D> f0(NX);
    sys.first(x, f0, D(cp0.t));
    st.initialize(x, D(cp0.t), f0, D(cp0.carry.h));
    st.load_carry(cp0.carry, NX);
  }

  S_obs.assign(store.n_obs() * NX * nd, 0.0);
  x_ref.assign(store.n_obs() * NX, 0.0);
  auto record_obs = [&](std::size_t o, const std::vector<D>& xs) {
    for (std::size_t i = 0; i < NX; ++i) {
      x_ref[o * NX + i] = xs[i].x();
      for (unsigned d = 0; d < nd; ++d)
        S_obs[(o * NX + i) * nd + d] = xs[i][d] / scale[d];
    }
  };

  std::size_t next_obs = 0;
  while (next_obs < store.n_obs() && store.obs(next_obs).step == 0) {
    record_obs(next_obs, x);
    ++next_obs;
  }

  for (std::size_t k = 0; k < store.n_steps(); ++k) {
    const auto& cp = store.step(k);
    // The recorded decisions, not the ones this run would make: load_carry sets
    // the order, the step-size history and the counters the coefficients are
    // built from, and leaves the Nordsieck slots alone, so the duals chain while
    // the bookkeeping follows the run being differentiated.
    if constexpr (multistep) st.load_carry(cp.carry, NX);
    st.do_step(sys, x, D(cp.t), xout, D(cp.dt), xerr);
    if constexpr (multistep) {
      check(st.newton_converged(), "the reference corrector converged");
      // The same tail the reverse replay runs, from the same record: the
      // Nordsieck update, the order the controller chose, and every rescale it
      // applied before the next step, thrown-away attempts included.
      cp.apply_tail(st, sys);
      for (std::size_t i = 0; i < NX; ++i) x[i] = st.zn_mut(0)[i];
    } else {
      st.prepare_dense_output();
    }

    while (next_obs < store.n_obs() && store.obs(next_obs).step == k + 1) {
      const double t_obs = store.obs(next_obs).t;
      if constexpr (multistep) {
        st.eval_dense_into(D(t_obs), x_interp);
      } else {
        const D t0(cp.t), t1(cp.t + cp.dt);
        st.calc_state(D(t_obs), x_interp, x, t0, xout, t1);
      }
      record_obs(next_obs, x_interp);
      ++next_obs;
    }

    if constexpr (!multistep) x = xout;
  }
  check(next_obs == store.n_obs(), "every observation reached by the reference");
}

// Replays the recorded sequence in double, for a finite-difference probe of the
// same prescribed discretisation both modes differentiate.
template<class S>
static void replay_values(const store_of<S>& store, std::vector<double>& x_obs,
                          bool verify_carry = false, const char* name = "")
{
  constexpr bool multistep = cppde::reverse::has_step_snapshot<S>::value;
  std::vector<double> p(P, P + NP);
  auto sys = make_system<double>(p);
  S st;
  if constexpr (multistep) st.set_tolerances(ATOL, RTOL);

  const auto& cp0 = store.step(0);
  std::vector<double> x(NX), xout(NX), xerr(NX), xi(NX);
  for (std::size_t i = 0; i < NX; ++i) x[i] = cp0.start_state()[i];
  if constexpr (multistep) {
    std::vector<double> f0(NX);
    sys.first(x, f0, cp0.t);
    st.initialize(x, cp0.t, f0, cp0.carry.h);
    st.load_carry(cp0.carry, NX);
  }

  x_obs.assign(store.n_obs() * NX, 0.0);
  std::size_t o = 0;
  while (o < store.n_obs() && store.obs(o).step == 0) {
    for (std::size_t i = 0; i < NX; ++i) x_obs[o * NX + i] = x[i];
    ++o;
  }
  for (std::size_t k = 0; k < store.n_steps(); ++k) {
    const auto& cp = store.step(k);
    if constexpr (multistep) st.load_carry(cp.carry, NX);
    st.do_step(sys, x, cp.t, xout, cp.dt, xerr);
    if constexpr (multistep) {
      cp.apply_tail(st, sys);
      for (std::size_t i = 0; i < NX; ++i) x[i] = st.zn_mut(0)[i];
    } else {
      st.prepare_dense_output();
    }
    while (o < store.n_obs() && store.obs(o).step == k + 1) {
      if constexpr (multistep) st.eval_dense_into(store.obs(o).t, xi);
      else st.calc_state(store.obs(o).t, xi, x, cp.t, xout, cp.t + cp.dt);
      for (std::size_t i = 0; i < NX; ++i) x_obs[o * NX + i] = xi[i];
      ++o;
    }
    if constexpr (multistep) {
      // The carry this tail produced has to be the carry the run recorded for
      // the next step, or the checkpoints do not describe the run. It holds
      // across thrown-away attempts too: what they did to the history is part
      // of the tail record.
      if (verify_carry && k + 1 < store.n_steps()) {
        const auto& nxt = store.step(k + 1);
        for (int j = 0; j <= nxt.carry.q; ++j)
          for (std::size_t i = 0; i < NX; ++i)
            close(nxt.zn[static_cast<std::size_t>(j) * NX + i], st.zn(j)[i],
                  std::string(name) + " carry handoff step " + std::to_string(k)
                  + " slot "
                  + std::to_string(j) + "." + std::to_string(i), 1e-9);
      }
    } else {
      x = xout;
    }
  }
}

template<class S>
static void run_method(const char* name, double tol)
{
  store_of<S> store;
  std::vector<double> x_run;
  forward_value<S>(store, x_run);

  // One map for every method: the initial state and the parameters.
  const std::size_t n_carry = NX;
  const unsigned nd = static_cast<unsigned>(n_carry + NP);

  std::vector<double> S_obs, x_ref;
  forward_sens<S>(store, n_carry, nd, S_obs, x_ref);
  for (std::size_t k = 0; k < x_run.size() && k < x_ref.size(); ++k)
    close(x_run[k], x_ref[k], std::string(name) + " observed value " + std::to_string(k), tol);

  std::printf("%-8s steps %3zu  observations %zu  carry %zu\n",
              name, store.n_steps(), store.n_obs(), n_carry);
  if constexpr (cppde::reverse::has_step_snapshot<S>::value) {
    // The carry a step hands on has to be the carry the next one reads, or the
    // cotangents are seeded onto the wrong slots and nothing says so.
    std::size_t bad = 0;
    for (std::size_t k = 0; k + 1 < store.n_steps(); ++k)
      if (store.step(k).q_next != store.step(k + 1).carry.q) ++bad;
    std::printf("  order handoff mismatches: %zu\n", bad);
  }
  check(store.n_steps() > 4, std::string(name) + " took more than a handful of steps");
  check(store.n_obs() == TIMES.size(), std::string(name) + " one observation per time");

  {
    // A third source, independent of both modes: finite differences on the same
    // prescribed sequence. Loose by nature, and only ever a sanity check, since
    // the oracle is the forward mode.
    // Unperturbed first: the carry this replay hands on has to be the carry the
    // run recorded, or the checkpoints do not describe the run.
    std::vector<double> plain;
    replay_values<S>(store, plain, true, name);

    const double h = 1e-6 * std::max(1.0, std::fabs(P[0]));
    std::vector<double> hi, lo;
    const double saved = P[0];
    P[0] = saved + h; replay_values<S>(store, hi);
    P[0] = saved - h; replay_values<S>(store, lo);
    P[0] = saved;
    for (std::size_t k = 0; k < hi.size(); ++k)
      close((hi[k] - lo[k]) / (2 * h), S_obs[k * nd + n_carry],
            std::string(name) + " fd d(obs" + std::to_string(k) + ")/dp0", 1e-5);
  }

  std::vector<double> pv(P, P + NP);
  auto jac_d = jacobian<double>{pv};

  const std::size_t n_seed = store.n_obs() * NX;
  std::size_t max_nodes = 0;
  auto sweep_with = [&](const std::vector<double>& seeds, std::vector<double>& out) {
    cppde::reverse::equation_solver<jacobian<double>, double> solver(jac_d);
    cppde::reverse::trajectory_recorder<S, double> rev;
    rev.sweep(store, pv,
              [](const std::vector<codual<double>>& pc) {
                return make_system<codual<double>>(pc);
              },
              seeds, solver);
    if (rev.max_tape_nodes() > max_nodes) max_nodes = rev.max_tape_nodes();
    // An observation before the first step never passes through one, so the
    // replay has nothing to interpolate there.
    for (std::size_t o = 0; o < store.n_obs(); ++o)
      if (store.obs(o).step > 0)
        for (std::size_t i = 0; i < NX; ++i)
          close(x_run[o * NX + i], rev.replayed_obs()[o * NX + i],
                std::string(name) + " replayed value " + std::to_string(o * NX + i), tol);
    out.assign(nd, 0.0);
    for (std::size_t i = 0; i < NX; ++i) out[i] = rev.wx0()[i];
    check(rev.whistory0().empty(),
          std::string(name) + " the start's history cotangent is collapsed");
    for (std::size_t j = 0; j < NP; ++j) out[n_carry + j] = rev.wp()[j];
  };

  auto compare = [&](const char* what, const std::vector<double>& seeds) {
    std::vector<double> got;
    sweep_with(seeds, got);
    for (unsigned d = 0; d < nd; ++d) {
      double wS = 0.0;
      for (std::size_t o = 0; o < store.n_obs(); ++o)
        for (std::size_t i = 0; i < NX; ++i)
          wS += seeds[o * NX + i] * S_obs[(o * NX + i) * nd + d];
      const std::string tag = (d < n_carry) ? "  dz" + std::to_string(d)
                                            : "  dp" + std::to_string(d - n_carry);
      close(wS, got[d], std::string(name) + " " + what + tag, tol);
    }
  };

  // One observation at a time, so no column of the trajectory Jacobian hides.
  for (std::size_t o = 0; o < store.n_obs(); ++o)
    for (std::size_t i = 0; i < NX; ++i) {
      std::vector<double> e(n_seed, 0.0);
      e[o * NX + i] = 1.0;
      compare(("obs" + std::to_string(o) + " e" + std::to_string(i)).c_str(), e);
    }

  std::vector<double> all(n_seed);
  for (std::size_t k = 0; k < n_seed; ++k)
    all[k] = 0.3 * static_cast<double>(k % 5) - 0.7;
  compare("all", all);

  // The bound the design claims: one step of tape, not one trajectory. Printed
  // per method because what a step costs differs by an order of magnitude
  // between them, rosenbrock4 putting its whole iteration matrix on the tape
  // where an explicit method puts only its stages.
  const std::size_t node = sizeof(cppde::codual_tape<double>::node);
  std::printf("  tape %zu nodes at the widest step, %zu bytes;"
              " the whole run taped would be %zu\n",
              max_nodes, max_nodes * node, max_nodes * store.n_steps() * node);
}

int main() {
  using cppde::multistep_method;
  // The two floors are the oracle's, not the adjoint's. An explicit method's
  // reference reproduces the recorded trajectory exactly, so it stands at 1e-10.
  // An implicit one's re-solves the corrector, and its stopping rule takes the
  // maximum over every sensitivity direction, so it stops elsewhere than the
  // value run did and the trajectory drifts by the residual it stopped on.
  // Where the adjoint itself is checked sharply is the step, in
  // test_reverse_step_multistep.cpp, which runs both methods at 1e-9 over every
  // order, order change, rescale and dense-output seed.
  run_method<cppde::multistepper<multistep_method::bdf, double, cppde::dense_lu_tag>>("bdf", 1e-6);
  run_method<cppde::multistepper<multistep_method::adams, double, cppde::dense_lu_tag>>("adams", 1e-6);
  run_method<cppde::rosenbrock4<double>>("rb4", 1e-10);
  run_method<cppde::tsit5<double>>("tsit5", 1e-10);



  std::printf(g_failures == 0 ? "\nOK\n" : "\n%d FAILURES\n", g_failures);
  return g_failures == 0 ? 0 : 1;
}
