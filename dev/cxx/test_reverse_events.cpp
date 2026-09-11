// The trajectory backwards across events: stage 5 of dev/adjoint-plan.md.
//
// An event is two maps between two steps. The jump, which for a fixed time is
// the reset carried across the discontinuity by the Heun sandwich and for a
// root the same sandwich around a crossing time the implicit function theorem
// supplies; and the restart, which for a multistep method throws the Nordsieck
// history away and builds a new one out of the post-jump state alone. The
// reverse pass replays both. Which events fired, and where the root sat, are
// control decisions and are read back rather than decided again.
//
// The oracle is the forward mode over the same trajectory and the same jumps:
// the saltation functions are templated on the scalar type, so both directions
// run the very same code and cannot drift apart. Both sides also follow the
// recorded step sequence rather than adapting again.
//
// Covered: a fixed-time reset whose value depends on a parameter, a
// root-triggered reset whose value does, on all four methods, so every carry
// shape is exercised; observations before, between and after the jumps,
// including the post-jump observation the engine emits, which is a value and
// not an interpolation; a seed one observation at a time and over all at once.
//
// Every number is printed at %.17g, so the output is the assertion as well.
//
// Build and run:  dev/cxx/run.sh --reverse-events
//
// Copyright (C) 2026 Simon Beyer

#include <cstdio>
#include <cmath>
#include <functional>
#include <string>
#include <utility>
#include <vector>

#include <cppde/cppde.hpp>

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
//  Model: three states, four parameters, non-autonomous so the crossing time
//  carries a derivative through the right-hand side as well as through g.
// ---------------------------------------------------------------------------

static constexpr std::size_t NX = 3;
static constexpr std::size_t NP = 4;

template<class V>
struct model {
  std::vector<V> p;
  void operator()(const std::vector<V>& x, std::vector<V>& d, const V& t) const {
    d[0] = -p[0] * x[0] + p[1] * x[1] * x[2];
    d[1] =  p[0] * x[0] - p[1] * x[1] * x[2] - p[2] * x[1] * x[1];
    d[2] =  p[2] * x[1] * x[1] * cppde::cos(t) - p[3] * x[2];
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
                      J(2,2) =  p[3];
    dfdt[0] = V(0.0);
    dfdt[1] = V(0.0);
    dfdt[2] = -p[2] * x[1] * x[1] * cppde::sin(t);
  }
};

template<class V>
static std::pair<model<V>, jacobian<V>> make_system(const std::vector<V>& p) {
  return std::make_pair(model<V>{p}, jacobian<V>{p});
}

static const double X0[NX] = {1.0, 0.55, 0.15};
static double       P [NP] = {0.9, 1.4, 0.6, 0.35};

static const double ATOL = 1e-9;
static const double RTOL = 1e-9;
static const double HINI = 0.02;

// The fixed reset sits at 1.0, the root condition on the way down before it.
static const double T_FIXED = 1.0;
static const double G_LEVEL = 0.70;

static const std::vector<double> TIMES = {0.0, 0.3, 0.9, 1.0, 1.4, 2.0, 3.0};

static constexpr unsigned ND = NX + NP;
using D = dual<double, ND>;

// ---------------------------------------------------------------------------
//  The events, on whichever scalar type the caller runs.
//
//  Both read a parameter, so both put a cotangent on wp: the fixed reset scales
//  x0 by a function of p3, the root reset adds a multiple of p1 to x2 when x0
//  falls through G_LEVEL. The root carries its partials, without which the
//  engine applies the reset with no transport and the crossing time drops out
//  of the derivative entirely.
// ---------------------------------------------------------------------------

template<class V>
struct event_set {
  std::vector<cppde::detail::FixedEvent<std::vector<V>, V>> fixed;
  std::vector<cppde::detail::RootEvent<std::vector<V>, V>>  root;
};

template<class V>
static event_set<V> make_events(const std::vector<V>& p) {
  using state = std::vector<V>;
  event_set<V> ev;

  cppde::detail::FixedEvent<state, V> f;
  f.time        = V(T_FIXED);
  f.state_index = 0;
  const V a = p[3];
  f.value_func  = [a](const state&, const V&) { return V(0.5) + a; };
  f.method      = cppde::detail::EventMethod::Multiply;
  ev.fixed.push_back(f);

  cppde::detail::RootEvent<state, V> r;
  r.func        = [](const state& x, const V&) { return x[0] - V(G_LEVEL); };
  r.state_index = 2;
  const V b = p[1];
  r.value_func  = [b](const state&, const V&) { return V(0.05) * b; };
  r.method      = cppde::detail::EventMethod::Add;
  r.terminal    = false;
  r.direction   = 0;
  r.dg_dx       = [](const state& x, const V&, state& g) {
    g.assign(x.size(), V(0)); g[0] = V(1);
  };
  r.dg_dt       = [](const state&, const V&) { return V(0); };
  ev.root.push_back(r);
  return ev;
}

// ---------------------------------------------------------------------------
//  Driver plumbing, one shape per stepper family.
// ---------------------------------------------------------------------------

// What the generator emits beside the model for a reverse build.
struct adjoint_terms {
  std::vector<double> p;

  void jac_t_vec(const std::vector<double>& x, const std::vector<double>& lam,
                 const double& t, std::vector<double>& out) const {
    out.assign(NX, 0.0);
    out[0] = (-p[0]) * lam[0] + (p[0]) * lam[1];
    out[1] = (p[1] * x[2]) * lam[0]
           + (-p[1] * x[2] - 2.0 * p[2] * x[1]) * lam[1]
           + (2.0 * p[2] * x[1] * std::cos(t)) * lam[2];
    out[2] = (p[1] * x[1]) * lam[0] + (-p[1] * x[1]) * lam[1]
           + (-p[3]) * lam[2];
  }

  void dfdp_t_vec_axpy(const std::vector<double>& x,
                       const std::vector<double>& lam,
                       const double& t, const double& sc,
                       double* out) const {
    out[NX + 0] += sc * ((-x[0]) * lam[0] + (x[0]) * lam[1]);
    out[NX + 1] += sc * ((x[1] * x[2]) * lam[0] + (-x[1] * x[2]) * lam[1]);
    out[NX + 2] += sc * ((-x[1] * x[1]) * lam[1]
                         + (x[1] * x[1] * std::cos(t)) * lam[2]);
    out[NX + 3] += sc * ((-x[2]) * lam[2]);
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
    out[NX + 3] += sc * ((-v[2]) * lam[2]);
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

// And what it emits for the jump: the derivatives of the event expressions
// themselves. The fixed event multiplies state 0 by 0.5 + p3 at a constant
// time; the root event adds 0.05 * p1 to state 2 where x0 crosses a level.
struct event_adjoint_terms {
  std::vector<double> p;

  void fixed_dh_dx(int, const std::vector<double>&, const double&,
                   std::vector<double>& out) const { out.assign(NX, 0.0); }
  void fixed_dh_dp_axpy(int ev, const std::vector<double>&, const double&,
                        const double& sc, double* out) const {
    if (ev == 0) out[NX + 3] += sc * 1.0;
  }
  void fixed_dtime_dp_axpy(int, const double&, double*) const {}

  void root_dh_dx(int, const std::vector<double>&, const double&,
                  std::vector<double>& out) const { out.assign(NX, 0.0); }
  void root_dh_dp_axpy(int ev, const std::vector<double>&, const double&,
                       const double& sc, double* out) const {
    if (ev == 0) out[NX + 1] += sc * 0.05;
  }
  void root_dg_dp_axpy(int, const std::vector<double>&, const double&,
                       const double&, double*) const {}
};

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

static constexpr std::size_t NPOS = cppde::reverse::event_record<double>::npos;

// ---------------------------------------------------------------------------
//  The forward value run, with both collectors hooked into the dense driver.
// ---------------------------------------------------------------------------

template<class S>
static void forward_value(store_of<S>& store, std::vector<double>& x_obs)
{
  std::vector<double> p(P, P + NP);
  auto sys = make_system<double>(p);
  auto ev  = make_events<double>(p);

  typename pipeline<S>::controller ctl(ATOL, RTOL);
  typename pipeline<S>::dense dense(std::move(ctl));

  std::vector<double> x(X0, X0 + NX);
  cppde::StepChecker checker(1000000, 1000000);

  cppde::reverse::step_collector<decltype(dense), S, double> coll(store, dense, HINI);

  x_obs.clear();
  auto obs = [&](const std::vector<double>& xs, const double& t) {
    store.observe(t);
    x_obs.insert(x_obs.end(), xs.begin(), xs.end());
  };

  cppde::integrate_times_dense(dense, sys, x, TIMES.begin(), TIMES.end(), HINI,
                               obs, ev.fixed, ev.root, checker, 1e-10, 1,
                               cppde::detail::no_dt_estimator{}, nullptr,
                               std::ref(coll), coll.event_observer());
}

// ---------------------------------------------------------------------------
//  Forward reference: the recorded sequence replayed under dual, seeded with
//  the identity over the initial state and the parameters. Every decision comes
//  from the store, every value from the arithmetic.
// ---------------------------------------------------------------------------

template<class S>
static void forward_sens(const store_of<S>& store, std::vector<double>& S_obs,
                         std::vector<double>& x_ref)
{
  using rstep = typename S::template rebind_value<D>;
  constexpr bool multistep = cppde::reverse::has_step_snapshot<S>::value;

  std::vector<double> scale(ND, 1.0);
  std::vector<D> p(NP);
  for (std::size_t j = 0; j < NP; ++j) {
    scale[NX + j] = std::max(std::fabs(P[j]), 1e-2);
    p[j] = D(P[j]);
    p[j].diff(static_cast<unsigned>(NX + j));
    p[j].d(static_cast<unsigned>(NX + j)) = scale[NX + j];
  }
  auto sys = make_system<D>(p);
  auto ev  = make_events<D>(p);

  rstep st;
  if constexpr (multistep) {
    st.set_tolerances(ATOL, RTOL);
    st.set_max_corrector_iters(200);
  }
  st.prepare_sensitivities(ND);

  const auto& cp0 = store.step(0);
  std::vector<D> x(NX), xout(NX), xerr(NX), x_interp(NX), f0(NX);
  for (std::size_t i = 0; i < NX; ++i) {
    scale[i] = std::max(std::fabs(cp0.start_state()[i]), 1e-2);
    x[i] = D(cp0.start_state()[i]);
    x[i].diff(static_cast<unsigned>(i));
    x[i].d(static_cast<unsigned>(i)) = scale[i];
  }

  S_obs.assign(store.n_obs() * NX * ND, 0.0);
  x_ref.assign(store.n_obs() * NX, 0.0);
  auto record_obs = [&](std::size_t o, const std::vector<D>& xs) {
    for (std::size_t i = 0; i < NX; ++i) {
      x_ref[o * NX + i] = xs[i].x();
      for (unsigned d = 0; d < ND; ++d)
        S_obs[(o * NX + i) * ND + d] = xs[i][d] / scale[d];
    }
  };

  // A jump the store recorded, applied on whichever scalar type is running.
  auto apply_jump = [&](std::vector<D>& xa, const std::vector<D>& xb,
                        const cppde::reverse::event_record<double>& e) {
    xa = xb;
    if (e.root) {
      cppde::detail::saltation_root_analytical_batch(xa, xb, D(e.t), sys,
                                                     ev.root, e.triggered);
    } else {
      auto at_surface = [&](std::vector<D>& xs, const D& ts) {
        for (std::size_t i : e.switched)
          cppde::detail::apply_event_action(xs, xs, ts, ev.root[i]);
      };
      cppde::detail::apply_fixed_events_at_time(xa, D(e.t), ev.fixed, sys,
                                                at_surface);
    }
  };

  // The trajectory start, and any reset that sat on it.
  {
    const std::size_t e0 = store.event_before(0);
    if (e0 != NPOS) {
      std::vector<D> xb = x;
      apply_jump(x, xb, store.event(e0));
      for (std::size_t o = 0; o < store.n_obs(); ++o)
        if (store.obs(o).event == e0) record_obs(o, x);
    }
    if constexpr (multistep) {
      sys.first(x, f0, D(cp0.t));
      st.initialize(x, D(cp0.t), f0, D(cp0.carry.h));
      st.load_carry(cp0.carry, NX);
    }
  }

  std::size_t next_obs = 0;
  while (next_obs < store.n_obs() && store.obs(next_obs).step == 0) {
    if (store.obs(next_obs).event == NPOS) record_obs(next_obs, x);
    ++next_obs;
  }

  for (std::size_t k = 0; k < store.n_steps(); ++k) {
    const auto& cp = store.step(k);
    if constexpr (multistep) st.load_carry(cp.carry, NX);
    st.do_step(sys, x, D(cp.t), xout, D(cp.dt), xerr);
    if constexpr (multistep) {
      check(st.newton_converged(), "the reference corrector converged");
      cp.apply_tail(st, sys);
      for (std::size_t i = 0; i < NX; ++i) x[i] = st.zn_mut(0)[i];
    } else {
      st.prepare_dense_output();
    }

    auto interp = [&](double t, std::vector<D>& out) {
      if constexpr (multistep) {
        st.eval_dense_into(D(t), out);
      } else {
        st.calc_state(D(t), out, x, D(cp.t), xout, D(cp.t + cp.dt));
      }
    };

    while (next_obs < store.n_obs() && store.obs(next_obs).step == k + 1) {
      if (store.obs(next_obs).event == NPOS) {
        interp(store.obs(next_obs).t, x_interp);
        record_obs(next_obs, x_interp);
      }
      ++next_obs;
    }

    const std::size_t ei = store.event_before(k + 1);
    if (ei != NPOS) {
      const auto& e = store.event(ei);
      std::vector<D> xb(NX), xa(NX);
      interp(e.t_before, xb);
      apply_jump(xa, xb, e);
      for (std::size_t o = 0; o < store.n_obs(); ++o)
        if (store.obs(o).event == ei) record_obs(o, xa);
      if constexpr (multistep) {
        sys.first(xa, f0, D(e.t));
        st.initialize(xa, D(e.t), f0, D(e.dt_restart));
      } else {
        // The engine restarts the stepper on the far side, which drops the
        // FSAL stage. Reusing it would carry a derivative from before the jump
        // into the step after it.
        st.invalidate_lu();
      }
      x = xa;
    } else if constexpr (!multistep) {
      x = xout;
    }
  }
  check(next_obs == store.n_obs(), "every observation reached by the reference");
}

// ---------------------------------------------------------------------------

template<class S>
static void run_method(const char* name, double tol)
{
  store_of<S> store;
  std::vector<double> x_run;
  forward_value<S>(store, x_run);

  std::vector<double> S_obs, x_ref;
  forward_sens<S>(store, S_obs, x_ref);
  for (std::size_t k = 0; k < x_run.size() && k < x_ref.size(); ++k)
    close(x_run[k], x_ref[k],
          std::string(name) + " observed value " + std::to_string(k), tol);

  std::size_t n_rooted = 0;
  for (std::size_t i = 0; i < store.n_events(); ++i)
    if (store.event(i).root) ++n_rooted;
  std::printf("%-8s steps %3zu  observations %zu  events %zu (%zu rooted)\n",
              name, store.n_steps(), store.n_obs(), store.n_events(), n_rooted);
  check(store.n_events() == 2, std::string(name) + " both jumps recorded");
  check(n_rooted == 1, std::string(name) + " one of them root-triggered");

  std::vector<double> pv(P, P + NP);
  auto jac_d = jacobian<double>{pv};

  const std::size_t n_seed = store.n_obs() * NX;
  auto sweep_with = [&](const std::vector<double>& seeds, std::vector<double>& out) {
    adjoint_terms adj{pv};
    event_adjoint_terms eadj{pv};
    auto sysd = make_system<double>(pv);
    auto evd  = make_events<double>(pv);
    auto jmp  = cppde::adjoint::make_jump_terms(sysd, eadj, evd.fixed, evd.root);
    out.assign(ND, 0.0);
    if constexpr (cppde::reverse::has_step_snapshot<S>::value) {
      cppde::reverse::equation_solver<jacobian<double>, double> solver(jac_d);
      cppde::adjoint::closed_multistep_trajectory<S> tr;
      tr.sweep(store, NX + NP, seeds.data(), adj, solver, jmp);
      for (std::size_t i = 0; i < NX; ++i) out[i] = tr.wx0()[i];
      for (std::size_t j = 0; j < NP; ++j) out[NX + j] = tr.wp()[NX + j];
    } else {
      S st;
      cppde::adjoint::closed_onestep_trajectory<S> tr;
      tr.sweep(store, NX + NP, seeds.data(), sysd, adj, st, jmp);
      for (std::size_t i = 0; i < NX; ++i) out[i] = tr.wx0()[i];
      for (std::size_t j = 0; j < NP; ++j) out[NX + j] = tr.wp()[NX + j];
    }
  };

  auto compare = [&](const std::string& what, const std::vector<double>& seeds) {
    std::vector<double> got;
    sweep_with(seeds, got);
    for (unsigned d = 0; d < ND; ++d) {
      double wS = 0.0;
      for (std::size_t o = 0; o < store.n_obs(); ++o)
        for (std::size_t i = 0; i < NX; ++i)
          wS += seeds[o * NX + i] * S_obs[(o * NX + i) * ND + d];
      const std::string tag = (d < NX) ? "  dx" + std::to_string(d)
                                       : "  dp" + std::to_string(d - NX);
      close(wS, got[d], std::string(name) + " " + what + tag, tol);
    }
  };

  for (std::size_t o = 0; o < store.n_obs(); ++o)
    for (std::size_t i = 0; i < NX; ++i) {
      std::vector<double> e(n_seed, 0.0);
      e[o * NX + i] = 1.0;
      compare("obs" + std::to_string(o) + " e" + std::to_string(i), e);
    }

  std::vector<double> all(n_seed);
  for (std::size_t k = 0; k < n_seed; ++k)
    all[k] = 0.3 * static_cast<double>(k % 5) - 0.7;
  compare("all", all);

}

int main() {
  using cppde::multistep_method;
  // All four, because a jump is the one place where the three carry shapes
  // differ: the multistep methods throw their Nordsieck history away and build
  // a new one, the one-step methods carry the state alone, and rb4 solves its
  // stages against a matrix the restart invalidates.
  run_method<cppde::multistepper<multistep_method::bdf, double, cppde::dense_lu_tag>>("bdf", 1e-6);
  run_method<cppde::multistepper<multistep_method::adams, double, cppde::dense_lu_tag>>("adams", 1e-6);
  run_method<cppde::rosenbrock4<double>>("rb4", 1e-9);
  run_method<cppde::tsit5<double>>("tsit5", 1e-9);

  std::printf(g_failures == 0 ? "\nOK\n" : "\n%d FAILURES\n", g_failures);
  return g_failures == 0 ? 0 : 1;
}
