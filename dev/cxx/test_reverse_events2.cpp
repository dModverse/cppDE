// The trajectory backwards across events, at second order: forward over
// forward against forward over reverse on ONE recorded step sequence.
//
// What this adds over --reverse-events is the order, and over the R tests the
// shared grid: both directions replay the sequence the value run recorded, so a
// residue is algebra and not discretisation. What it adds over
// dev/jump-adjoint-check.cpp is the store. That check holds the jump's
// transpose against saltation_root_analytical_batch, and both sides share the
// convention that the localisation is tangent-blind; it is sharp for the
// transposition and blind to the convention. Here the store supplies the
// pre-jump state, as it does in a real solve.
//
// The root sits on the state whose right-hand side reads the clock, so
// d(g_dot)/dt is not zero. A root on a state with an autonomous right-hand side
// leaves that term at zero and the case goes untested. The fixed event's time
// depends on a parameter and its reset scales that same state, for the same
// reason: the shift then has a tangent, and the right-hand side read at the
// moving event time adds df/dt before and after a reset that changes it.
//
// Both jumps once failed here in the second order alone, the root one by 3.4e+01
// and the fixed one by 1.8e-01, on both steppers and on one step sequence: the
// gaps were algebra, rank one on the event time's parameter, and absent from
// the gradient.
//
// Build and run:  dev/cxx/run.sh --reverse-events2
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
//  Model: three states, four parameters, and a right-hand side that reads the
//  clock in the component the root watches.
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

static const double ATOL = 1e-11;
static const double RTOL = 1e-11;
static const double HINI = 0.02;

static const double T_FIXED = 1.0;
static const double G_LEVEL = 0.30;   // on x2, which reads the clock

static const std::vector<double> TIMES = {0.0, 0.3, 0.9, 1.0, 1.4, 2.0, 3.0};

static constexpr unsigned ND = NX + NP;
using D  = dual<double, ND>;
using D2 = cppde::dual2nd<double, ND>;

// ---------------------------------------------------------------------------
//  Events. The root watches x2 and adds a multiple of p1 to x1.
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

  // The fixed time moves with p3 and sits at T_FIXED in value, so the jump's
  // shift has a tangent and the right-hand side is read at a moving time. The
  // reset scales x2, whose right-hand side reads the clock: the df/dt terms on
  // either side of the reset then differ and do not cancel.
  cppde::detail::FixedEvent<state, V> f;
  f.time        = p[3] * V(T_FIXED / P[3]);
  f.state_index = 2;
  const V a = p[3];
  f.value_func  = [a](const state&, const V&) { return V(0.5) + a; };
  f.method      = cppde::detail::EventMethod::Multiply;
  ev.fixed.push_back(f);

  cppde::detail::RootEvent<state, V> r;
  r.func        = [](const state& x, const V&) { return x[2] - V(G_LEVEL); };
  r.state_index = 1;
  const V b = p[1];
  r.value_func  = [b](const state&, const V&) { return V(0.05) * b; };
  r.method      = cppde::detail::EventMethod::Add;
  r.terminal    = false;
  r.direction   = 0;
  r.dg_dx       = [](const state& x, const V&, state& g) {
    g.assign(x.size(), V(0)); g[2] = V(1);
  };
  r.dg_dt       = [](const state&, const V&) { return V(0); };
  // g_dot = f2, so g_ddot = df2/dt + grad f2 . f, and the first term is what a
  // root on an autonomous component never sees.
  const std::vector<V> pc = p;
  r.g_dot_dot   = [pc](const state& x, const V& t) -> double {
    const double x1 = cppde::ad_traits::scalar_value(x[1]);
    const double x2 = cppde::ad_traits::scalar_value(x[2]);
    const double tv = cppde::ad_traits::scalar_value(t);
    double pv[NP];
    for (std::size_t j = 0; j < NP; ++j)
      pv[j] = cppde::ad_traits::scalar_value(pc[j]);
    const double f1 = pv[0] * cppde::ad_traits::scalar_value(x[0])
                    - pv[1] * x1 * x2 - pv[2] * x1 * x1;
    const double f2 = pv[2] * x1 * x1 * std::cos(tv) - pv[3] * x2;
    return -pv[2] * x1 * x1 * std::sin(tv)
         + 2.0 * pv[2] * x1 * std::cos(tv) * f1
         - pv[3] * f2;
  };
  ev.root.push_back(r);
  return ev;
}

// ---------------------------------------------------------------------------
//  What the generator emits beside the model, on whichever type is running.
// ---------------------------------------------------------------------------

template<class V>
struct adjoint_terms {
  std::vector<V> p;

  void jac_t_vec(const std::vector<V>& x, const std::vector<V>& lam,
                 const V& t, std::vector<V>& out) const {
    out.assign(NX, V(0.0));
    out[0] = (V(0.0) - p[0]) * lam[0] + p[0] * lam[1];
    out[1] = (p[1] * x[2]) * lam[0]
           + (V(0.0) - p[1] * x[2] - V(2.0) * p[2] * x[1]) * lam[1]
           + (V(2.0) * p[2] * x[1] * cppde::cos(t)) * lam[2];
    out[2] = (p[1] * x[1]) * lam[0] + (V(0.0) - p[1] * x[1]) * lam[1]
           + (V(0.0) - p[3]) * lam[2];
  }

  void dfdp_t_vec_axpy(const std::vector<V>& x, const std::vector<V>& lam,
                       const V& t, const V& sc, V* out) const {
    out[NX + 0] += sc * ((V(0.0) - x[0]) * lam[0] + x[0] * lam[1]);
    out[NX + 1] += sc * ((x[1] * x[2]) * lam[0] + (V(0.0) - x[1] * x[2]) * lam[1]);
    out[NX + 2] += sc * ((V(0.0) - x[1] * x[1]) * lam[1]
                         + (x[1] * x[1] * cppde::cos(t)) * lam[2]);
    out[NX + 3] += sc * ((V(0.0) - x[2]) * lam[2]);
  }

  void jvp_x_t_vec(const std::vector<V>& x, const std::vector<V>& v,
                   const std::vector<V>& lam, const V& t,
                   std::vector<V>& out) const {
    out.assign(NX, V(0.0));
    out[1] = (p[1] * v[2]) * lam[0]
           + (V(0.0) - V(2.0) * p[2] * v[1] - p[1] * v[2]) * lam[1]
           + (V(2.0) * p[2] * cppde::cos(t) * v[1]) * lam[2];
    out[2] = (p[1] * v[1]) * lam[0] + (V(0.0) - p[1] * v[1]) * lam[1];
    (void)x;
  }

  void jvp_p_t_vec_axpy(const std::vector<V>& x, const std::vector<V>& v,
                        const std::vector<V>& lam, const V& t,
                        const V& sc, V* out) const {
    const V q = x[2] * v[1] + x[1] * v[2];
    out[NX + 0] += sc * ((V(0.0) - v[0]) * lam[0] + v[0] * lam[1]);
    out[NX + 1] += sc * (q * lam[0] + (V(0.0) - q) * lam[1]);
    out[NX + 2] += sc * ((V(0.0) - V(2.0) * x[1] * v[1]) * lam[1]
                         + (V(2.0) * x[1] * cppde::cos(t) * v[1]) * lam[2]);
    out[NX + 3] += sc * ((V(0.0) - v[2]) * lam[2]);
  }

  void dfdt_x_t_vec(const std::vector<V>& x, const std::vector<V>& lam,
                    const V& t, std::vector<V>& out) const {
    out.assign(NX, V(0.0));
    out[1] = (V(0.0) - V(2.0) * p[2] * x[1] * cppde::sin(t)) * lam[2];
  }

  void dfdt_p_t_vec_axpy(const std::vector<V>& x, const std::vector<V>& lam,
                         const V& t, const V& sc, V* out) const {
    out[NX + 2] += sc * ((V(0.0) - x[1] * x[1] * cppde::sin(t)) * lam[2]);
  }

  V dfdt_dot(const std::vector<V>& x, const std::vector<V>& lam,
             const V& t) const {
    return (V(0.0) - p[2] * x[1] * x[1] * cppde::sin(t)) * lam[2];
  }


};

template<class V>
struct event_adjoint_terms {
  std::vector<V> p;

  void fixed_dh_dx(int, const std::vector<V>&, const V&,
                   std::vector<V>& out) const { out.assign(NX, V(0.0)); }
  void fixed_dh_dp_axpy(int ev, const std::vector<V>&, const V&,
                        const V& sc, V* out) const {
    if (ev == 0) out[NX + 3] += sc;
  }
  void fixed_dtime_dp_axpy(int ev, const V& sc, V* out) const {
    if (ev == 0) out[NX + 3] += sc * V(T_FIXED / P[3]);
  }
  void fixed_dh_dt_axpy(int, const std::vector<V>&, const V&, const V&,
                        V*) const {}

  void root_dh_dx(int, const std::vector<V>&, const V&,
                  std::vector<V>& out) const { out.assign(NX, V(0.0)); }
  void root_dh_dp_axpy(int ev, const std::vector<V>&, const V&,
                       const V& sc, V* out) const {
    if (ev == 0) out[NX + 1] += sc * V(0.05);
  }
  void root_dg_dp_axpy(int, const std::vector<V>&, const V&,
                       const V&, V*) const {}
  void root_dh_dt_axpy(int, const std::vector<V>&, const V&, const V&,
                       V*) const {}

  // grad g_dot for g = x2 - level, so g_dot = f2.
  void root_gdot_dx(int, const std::vector<V>& x, const V& t,
                    std::vector<V>& out) const {
    out.assign(NX, V(0.0));
    out[1] = V(2.0) * p[2] * x[1] * cppde::cos(t);
    out[2] = V(0.0) - p[3];
  }
  void root_gdot_dp_axpy(int, const std::vector<V>& x, const V& t,
                         const V& sc, V* out) const {
    out[NX + 2] += sc * (x[1] * x[1] * cppde::cos(t));
    out[NX + 3] += sc * (V(0.0) - x[2]);
  }
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

static constexpr std::size_t NPOS = cppde::reverse::event_record<D>::npos;

// ---------------------------------------------------------------------------
//  The value run, on the type forward over reverse runs it on: a dual seeded
//  with the identity, so the store's event states carry dx/dtheta.
//
//  One-step methods only. A multistep carry is typed by the scalar it was
//  recorded on, so a dual2nd replay cannot load a dual store's Nordsieck
//  history; the gap this harness is aimed at is the same on all four methods.
// ---------------------------------------------------------------------------

template<class SAD>
static void value_run(cppde::reverse::trajectory_store<SAD, D>& store,
                      std::vector<D>& p, std::vector<double>& scale)
{
  auto sys = make_system<D>(p);
  auto ev  = make_events<D>(p);

  typename pipeline<SAD>::controller ctl(ATOL, RTOL);
  typename pipeline<SAD>::dense dense(std::move(ctl));

  std::vector<D> x(NX);
  for (std::size_t i = 0; i < NX; ++i) {
    scale[i] = std::max(std::fabs(X0[i]), 1e-2);
    x[i] = D(X0[i]);
    x[i].diff(static_cast<unsigned>(i));
    x[i].d(static_cast<unsigned>(i)) = scale[i];
  }

  cppde::StepChecker checker(1000000, 1000000);
  cppde::reverse::step_collector<decltype(dense), SAD, D> coll(store, dense, HINI);
  auto obs = [&](const std::vector<D>&, const D& t) { store.observe(t.x()); };

  cppde::integrate_times_dense(dense, sys, x, TIMES.begin(), TIMES.end(), D(HINI),
                               obs, ev.fixed, ev.root, checker, 1e-12, 1,
                               cppde::detail::no_dt_estimator{}, nullptr,
                               std::ref(coll), coll.event_observer());
}

// ---------------------------------------------------------------------------
//  The oracle: the recorded sequence replayed over a nested dual. Every
//  decision comes from the store, so this is the same discretisation.
// ---------------------------------------------------------------------------

template<class SAD, class S>
static void forward_second(const cppde::reverse::trajectory_store<SAD, D>& store,
                           const std::vector<double>& scale,
                           std::vector<double>& G_obs,
                           std::vector<double>& H_obs)
{
  using S2 = typename S::template rebind_value<D2>;

  std::vector<D2> p(NP);
  for (std::size_t j = 0; j < NP; ++j) {
    p[j].arm_full(ND);
    p[j].scalar() = P[j];
    p[j].d1_at(static_cast<unsigned>(NX + j)) = scale[NX + j];
  }
  auto sys = make_system<D2>(p);
  auto ev  = make_events<D2>(p);

  S2 st;
  st.prepare_sensitivities(ND);

  const auto& cp0 = store.step(0);
  std::vector<D2> x(NX), xout(NX), xerr(NX), x_interp(NX);
  for (std::size_t i = 0; i < NX; ++i) {
    x[i].arm_full(ND);
    x[i].scalar() = cp0.start_state()[i].x();
    x[i].d1_at(static_cast<unsigned>(i)) = scale[i];
  }
  for (auto& v : {&xout, &xerr, &x_interp})
    for (auto& e : *v) e.arm_full(ND);

  G_obs.assign(store.n_obs() * NX * ND, 0.0);
  H_obs.assign(store.n_obs() * NX * ND * ND, 0.0);
  auto record = [&](std::size_t o, const std::vector<D2>& xs) {
    for (std::size_t i = 0; i < NX; ++i)
      for (unsigned a = 0; a < ND; ++a) {
        G_obs[(o * NX + i) * ND + a] = xs[i].d1_at(a) / scale[a];
        for (unsigned b = 0; b < ND; ++b)
          H_obs[((o * NX + i) * ND + a) * ND + b] =
              xs[i].dd_at(a, b) / (scale[a] * scale[b]);
      }
  };

  auto apply_jump = [&](std::vector<D2>& xa, const std::vector<D2>& xb,
                        const cppde::reverse::event_record<D>& e) {
    xa = xb;
    if (e.root) {
      cppde::detail::saltation_root_analytical_batch(xa, xb, D2(e.t), sys,
                                                     ev.root, e.triggered);
    } else {
      auto at_surface = [&](std::vector<D2>& xs, const D2& ts) {
        for (std::size_t i : e.switched)
          cppde::detail::apply_event_action(xs, xs, ts, ev.root[i]);
      };
      cppde::detail::apply_fixed_events_at_time(xa, D2(e.t), ev.fixed, sys,
                                                at_surface);
    }
  };

  {
    const std::size_t e0 = store.event_before(0);
    if (e0 != NPOS) {
      std::vector<D2> xb = x;
      apply_jump(x, xb, store.event(e0));
      for (std::size_t o = 0; o < store.n_obs(); ++o)
        if (store.obs(o).event == e0) record(o, x);
    }
  }

  std::size_t next_obs = 0;
  while (next_obs < store.n_obs() && store.obs(next_obs).step == 0) {
    if (store.obs(next_obs).event == NPOS) record(next_obs, x);
    ++next_obs;
  }

  for (std::size_t k = 0; k < store.n_steps(); ++k) {
    const auto& cp = store.step(k);
    st.do_step(sys, x, D2(cp.t), xout, D2(cp.dt), xerr);
    st.prepare_dense_output();

    auto interp = [&](double t, std::vector<D2>& out) {
      st.calc_state(D2(t), out, x, D2(cp.t), xout, D2(cp.t + cp.dt));
    };

    while (next_obs < store.n_obs() && store.obs(next_obs).step == k + 1) {
      if (store.obs(next_obs).event == NPOS) {
        interp(store.obs(next_obs).t, x_interp);
        record(next_obs, x_interp);
      }
      ++next_obs;
    }

    const std::size_t ei = store.event_before(k + 1);
    if (ei != NPOS) {
      const auto& e = store.event(ei);
      std::vector<D2> xb(NX), xa(NX);
      for (auto& v : {&xb, &xa}) for (auto& q : *v) q.arm_full(ND);
      interp(e.t_before, xb);
      apply_jump(xa, xb, e);
      for (std::size_t o = 0; o < store.n_obs(); ++o)
        if (store.obs(o).event == ei) record(o, xa);
      st.invalidate_lu();
      x = xa;
    } else {
      x = xout;
    }
  }
  check(next_obs == store.n_obs(), "every observation reached by the oracle");
}

// ---------------------------------------------------------------------------

template<class S>
static void run_method(const char* name, double tol)
{
  using SAD = typename S::template rebind_value<D>;

  std::vector<double> scale(ND, 1.0);
  std::vector<D> p(NP);
  for (std::size_t j = 0; j < NP; ++j) {
    scale[NX + j] = std::max(std::fabs(P[j]), 1e-2);
    p[j] = D(P[j]);
    p[j].diff(static_cast<unsigned>(NX + j));
    p[j].d(static_cast<unsigned>(NX + j)) = scale[NX + j];
  }

  cppde::reverse::trajectory_store<SAD, D> store;
  value_run<SAD>(store, p, scale);

  std::size_t n_rooted = 0;
  for (std::size_t i = 0; i < store.n_events(); ++i)
    if (store.event(i).root) ++n_rooted;
  std::printf("%-8s steps %3zu  observations %zu  events %zu (%zu rooted)\n",
              name, store.n_steps(), store.n_obs(), store.n_events(), n_rooted);
  check(store.n_events() >= 2, std::string(name) + " both jumps recorded");
  check(n_rooted >= 1, std::string(name) + " a root among them");

  std::vector<double> G_obs, H_obs;
  forward_second<SAD, S>(store, scale, G_obs, H_obs);

  auto sweep_with = [&](const std::vector<D>& seeds, std::vector<D>& out) {
    adjoint_terms<D> adj{p};
    event_adjoint_terms<D> eadj{p};
    auto sysd = make_system<D>(p);
    auto evd  = make_events<D>(p);
    auto jmp  = cppde::adjoint::make_jump_terms(sysd, eadj, evd.fixed, evd.root);
    out.assign(ND, D(0.0));
    SAD st;
    // The sweep re-runs each step on this stepper, so its buffers need the
    // same tangent slab the forward one got.
    st.prepare_sensitivities(ND);
    cppde::adjoint::closed_onestep_trajectory<SAD> tr;
    tr.sweep(store, NX + NP, seeds.data(), sysd, adj, st, jmp);
    for (std::size_t i = 0; i < NX; ++i) out[i] = tr.wx0()[i] + tr.wp()[i];
    for (std::size_t j = 0; j < NP; ++j) out[NX + j] = tr.wp()[NX + j];
  };

  // One observation and one state at a time, so the contraction is the gradient
  // and Hessian of that row and the comparison is entry by entry.
  const std::size_t n_seed = store.n_obs() * NX;
  double e1 = 0.0, e2 = 0.0, sum2 = 0.0;
  std::string w1, w2;
  for (std::size_t o = 0; o < store.n_obs(); ++o)
    for (std::size_t i = 0; i < NX; ++i) {
      std::vector<D> e(n_seed, D(0.0));
      e[o * NX + i] = D(1.0);
      std::vector<D> got;
      sweep_with(e, got);
      for (unsigned a = 0; a < ND; ++a) {
        const double g1 = G_obs[(o * NX + i) * ND + a];
        const double h1 = got[a].x();
        if (std::fabs(g1 - h1) > e1) {
          e1 = std::fabs(g1 - h1);
          w1 = "obs" + std::to_string(o) + " x" + std::to_string(i) +
               " d" + std::to_string(a);
        }
        for (unsigned b = 0; b < ND; ++b) {
          const double g2 = H_obs[((o * NX + i) * ND + a) * ND + b];
          const double h2 = got[a].d(b) / scale[b];
          sum2 += std::fabs(g2 - h2);
          sum2 += std::fabs(g2 - h2);
          if (std::fabs(g2 - h2) > e2) {
            e2 = std::fabs(g2 - h2);
            w2 = "obs" + std::to_string(o) + " x" + std::to_string(i) +
                 " d" + std::to_string(a) + std::to_string(b);
          }
        }
      }
    }
  std::printf("%-8s first order  max |reverse - forward| = %.3e  at %s\n",
              name, e1, w1.c_str());
  std::printf("%-8s second order max |reverse - forward| = %.3e  at %s\n",
              name, e2, w2.c_str());
  std::printf("%-8s second order sum |reverse - forward| = %.10e\n", name, sum2);
  check(e1 <= tol, std::string(name) + " first order");
  check(e2 <= tol, std::string(name) + " second order");
}

int main() {
  run_method<cppde::tsit5<double>>("tsit5", 1e-6);
  run_method<cppde::rosenbrock4<double>>("rb4", 1e-6);
  std::printf(g_failures == 0 ? "\nOK\n" : "\n%d FAILURES\n", g_failures);
  return g_failures == 0 ? 0 : 1;
}
