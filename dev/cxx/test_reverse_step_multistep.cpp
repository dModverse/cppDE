// One multistep step backwards, BDF and Adams: stage 3b of dev/adjoint-plan.md.
//
// The oracle is the forward mode on the same step. What a multistep method
// carries across a step boundary is the Nordsieck history and not the state
// alone, so the sensitivity is that of the map
//
//   (zn_in, theta) -> zn_out,
//
// forward as the full matrix S and reverse as w' S for one w. Both differentiate
// the same discrete step, so they agree to rounding.
//
// The corrector is where the two differ in construction and must not differ in
// result. Forward iterates and peels the tangents by the implicit function
// theorem at every iteration; reverse puts the solution back and records the
// equation it solves, closing it by a transposed solve against the same matrix.
// Neither differentiates the iterates.
//
// Covered, for both methods: the unit seeds on every Nordsieck slot and a mixed
// one, at order 2 through 4 reached by a real controller warm-up, with the order
// held, raised and lowered, and with a rescale applied.
//
// Every number is printed at %.17g, so the output is the assertion as well.
//
// Build and run:  dev/cxx/run.sh --reverse-step-multistep
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

static void close(double a, double b, const std::string& what, double tol = 1e-9) {
  const double scale = std::fabs(a) > 1.0 ? std::fabs(a) : 1.0;
  char buf[128];
  std::snprintf(buf, sizeof buf, "  forward %.17g  reverse %.17g  rel %.2e",
                a, b, std::fabs(a - b) / scale);
  check(std::fabs(a - b) <= tol * scale, what + buf);
}

// ---------------------------------------------------------------------------
//  A small stiff, nonlinear, non-autonomous model, so the corrector does real
//  work and no partial is accidentally constant.
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

// Minus the Jacobian, which is what the solver's iteration matrix is built from:
// factorize_W only adds the diagonal, W = -J + (1/gamma) I. A generated model
// emits it the same way.
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

// What the code generator emits beside the model for a reverse build: J' lambda
// and (df/dp)' lambda, in plain double. Written out here by hand for the same
// model, which is what makes the closed-form adjoint checkable without a
// generated source.
struct adjoint_terms {
  std::vector<double> p;

  void jac_t_vec(const std::vector<double>& x, const std::vector<double>& lam,
                 const double& t, std::vector<double>& out) const {
    out.assign(NX, 0.0);
    // J(i,j) = d f_i / d x_j, contracted over i.
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
    out[NX + 0] += sc*((-x[0]) * lam[0] + (x[0]) * lam[1]);
    out[NX + 1] += sc*((x[1] * x[2]) * lam[0] + (-x[1] * x[2]) * lam[1]);
    out[NX + 2] += sc*((-x[1] * x[1]) * lam[1]
                + (x[1] * x[1] * std::cos(t)) * lam[2]);
  }
};

template<class V>
static std::pair<model<V>, jacobian<V>> make_system(const std::vector<V>& p) {
  return std::make_pair(model<V>{p}, jacobian<V>{p});
}

static const double X0[NX] = {1.0, 0.35, 0.15};
static const double P [NP] = {0.9, 1.4, 0.6};

static const double T0   = 0.1;
// Loose on purpose. The error norm takes the maximum over every sensitivity
// direction, and a unit seed on a Nordsieck slot is not a physical state, so a
// tight tolerance stalls the reference's corrector rather than sharpening it.
// What the comparison needs is that both sides converge, not that they are
// accurate: they differentiate the same equation either way.
static const double ATOL = 1e-6;
static const double RTOL = 1e-6;

// Directions: every Nordsieck slot a warm-up can produce, plus the parameters.
// Fixed at compile time so the reference runs on a static dual; the slots a case
// does not reach stay zero.
static constexpr unsigned ND_MAX = (12 + 1) * NX + NP;
using D = dual<double, ND_MAX>;

template<cppde::multistep_method M>
using stepper_d = cppde::multistepper<M, double, cppde::dense_lu_tag>;
template<cppde::multistep_method M>
using checkpoint = cppde::reverse::step_checkpoint<stepper_d<M>, double>;

// Iterates the corrector to machine precision, value and tangents together.
// PECE stops at its own dcon, which is only about the solver tolerance, and both
// modes are defined to differentiate the exactly solved equation, so that
// stopping rule would set the comparison's floor rather than rounding. Not for
// the BDF family: there the fixed point is Newton's, and this iteration does not
// converge to it on a stiff step.
template<cppde::multistep_method M, class V, class Sys, class St>
static void polish_corrector(Sys& sys, St& st, double t, std::vector<V>& y)
{
  if constexpr (M == cppde::multistep_method::adams) {
    const std::size_t n = y.size();
    const double gam = static_cast<double>(st.gamma());
    const double h   = static_cast<double>(st.h());
    const V t_new = V(t + h);
    std::vector<V> f(n);
    for (int it = 0; it < 60; ++it) {
      sys.first(y, f, t_new);
      for (std::size_t i = 0; i < n; ++i)
        y[i] = st.zn(0)[i] + gam * f[i] - (gam / h) * st.zn(1)[i];
    }
  } else {
    (void)sys; (void)st; (void)t; (void)y;
  }
}

// ---------------------------------------------------------------------------
//  Warm-up: a real controller run, so the order, the step-size history and the
//  Nordsieck slots are what a solve actually produces rather than set by hand.
// ---------------------------------------------------------------------------

template<cppde::multistep_method M>
static void warm_up(checkpoint<M>& cp, std::vector<double>& x, double& t, double& dt,
                    int n_steps)
{
  std::vector<double> p(P, P + NP);
  auto sys = make_system<double>(p);
  cppde::multistepper_controller<stepper_d<M>> ctl(ATOL, RTOL);

  x.assign(X0, X0 + NX);
  t  = T0;
  dt = 0.01;
  std::vector<double> xout(NX);
  for (int s = 0; s < n_steps; ++s) {
    while (ctl.try_step(sys, x, t, xout, dt) != cppde::success) {}
    x = xout;
  }
  cp.capture(ctl.stepper(), x, t, dt);
}

// ---------------------------------------------------------------------------
//  Forward reference: one step under dual, seeded with the identity over the
//  Nordsieck slots and the parameters, then the same tail the checkpoint's
//  finish() replays. S is [n_carry_out x ND] row-major.
// ---------------------------------------------------------------------------

template<cppde::multistep_method M>
static void forward_step(const checkpoint<M>& cp, int q_next, double eta,
                         unsigned nd, std::vector<double>& S,
                         std::vector<double>& carry_out, double t_interp = 0.0)
{
  const std::size_t n  = cp.n_states;
  const std::size_t nz = static_cast<std::size_t>(cp.carry.q + 1) * n;

  // Each direction is seeded at the magnitude of what it perturbs, not at one.
  // The corrector's error norm takes the maximum over the state and every
  // sensitivity direction, and a unit tangent on a high Nordsieck slot is many
  // orders above the slot itself, which stalls the iteration instead of
  // sharpening it. The map is linear in the seed, so S is unscaled again below.
  std::vector<double> scale(nd, 1.0);
  for (std::size_t k = 0; k < nz; ++k)
    scale[k] = std::max(std::fabs(cp.zn[k]), 1e-2);
  for (std::size_t j = 0; j < NP; ++j)
    scale[nz + j] = std::max(std::fabs(P[j]), 1e-2);

  std::vector<D> p(NP);
  for (std::size_t j = 0; j < NP; ++j) {
    p[j] = D(P[j]);
    p[j].diff(static_cast<unsigned>(nz + j));
    p[j].d(static_cast<unsigned>(nz + j)) = scale[nz + j];
  }
  auto sys = make_system<D>(p);

  cppde::multistepper<M, D, cppde::dense_lu_tag> st;
  st.set_tolerances(ATOL, RTOL);
  st.prepare_sensitivities(nd);
  st.load_carry(cp.carry, n);
  for (int j = 0; j <= cp.carry.q; ++j) {
    auto& slot = st.zn_mut(j);
    for (std::size_t i = 0; i < n; ++i) {
      const std::size_t k = static_cast<std::size_t>(j) * n + i;
      slot[i] = D(cp.zn[k]);
      slot[i].diff(static_cast<unsigned>(k));
      slot[i].d(static_cast<unsigned>(k)) = scale[k];
    }
  }

  std::vector<D> x(n), xout(n), xerr(n);
  for (std::size_t i = 0; i < n; ++i) x[i] = st.zn_mut(0)[i];

  st.do_step(sys, x, cp.t, xout, cp.dt, xerr);
  check(st.newton_converged(), "the forward reference corrector converged");
  polish_corrector<M>(sys, st, cp.t, xout);
  st.replay_outputs(xout, xout, xerr);
  st.complete_step();
  st.set_tn_current(cp.t + cp.dt);
  st.prepare_dense_output();
  if (q_next > cp.carry.q) st.save_acor_to_zn_qmax();
  if (q_next != cp.carry.q) st.set_order_for_next_step(q_next);
  if (std::abs(eta - 1.0) > 1e-14) st.rescale(eta);

  if (t_interp > 0.0) {
    // The dense output of this step, which is what an observation inside it
    // reaches. Seeded like a carry slot so the same comparison covers it.
    std::vector<D> xi(n);
    st.eval_dense_into(D(t_interp), xi);
    S.assign(n * nd, 0.0);
    carry_out.assign(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
      carry_out[i] = xi[i].x();
      for (unsigned d = 0; d < nd; ++d) S[i * nd + d] = xi[i][d] / scale[d];
    }
    return;
  }

  const int q_out = st.current_order();
  const std::size_t n_out = static_cast<std::size_t>(q_out + 1) * n;
  S.assign(n_out * nd, 0.0);
  carry_out.assign(n_out, 0.0);
  for (int j = 0; j <= q_out; ++j) {
    const auto& slot = st.zn(j);
    for (std::size_t i = 0; i < n; ++i) {
      const std::size_t r = static_cast<std::size_t>(j) * n + i;
      carry_out[r] = slot[i].x();
      for (unsigned d = 0; d < nd; ++d) S[r * nd + d] = slot[i][d] / scale[d];
    }
  }
}

// ---------------------------------------------------------------------------
//  Reverse: the same step replayed under codual, seeded on the carry out.
// ---------------------------------------------------------------------------

template<cppde::multistep_method M>
static void reverse_step(const checkpoint<M>& cp, const std::vector<double>& w,
                         std::vector<double>& wz, std::vector<double>& wp,
                         std::vector<double>& carry_out, double t_interp = 0.0)
{
  using C = codual<double>;
  const std::size_t n = cp.n_states;

  cppde::reverse::step_recorder<stepper_d<M>, double> rec;
  rec.begin();

  std::vector<C> p(NP);
  for (std::size_t j = 0; j < NP; ++j) p[j] = C(P[j]);
  rec.independent(p);
  auto sys = make_system<C>(p);

  rec.load(cp, cp.dt);

  // The equation's matrix at the solution, in plain doubles: the forward run's
  // own factorisation belongs to its iteration and is stale by design. Scaled by
  // gamma, since factorize_W builds W = (1/gamma) I - J and the residual's
  // derivative in y is gamma * W.
  std::vector<double> pv(P, P + NP);
  auto jac_d = jacobian<double>{pv};
  cppde::reverse::equation_solver<jacobian<double>, double> solver(jac_d);

  rec.attempt_implicit(sys, rec.dt_in(), cp.y, solver);
  // The equation is res = (y - zn0) + rl1*zn1 - gamma*f, so its derivative in y
  // is gamma * W, which is what the scale argument carries.
  const double gamma = rec.implicit_gamma();
  solver.prepare(cp.y, rec.implicit_t_new(), 1.0 / gamma, gamma);

  if (t_interp > 0.0) {
    std::vector<C> xi;
    rec.interpolate(t_interp, xi);
    carry_out.assign(xi.size(), 0.0);
    for (std::size_t i = 0; i < xi.size(); ++i) {
      carry_out[i] = xi[i].x();
      xi[i].seed(i < w.size() ? w[i] : 0.0);
    }
  } else {
    const auto& co = rec.carry_out();
    carry_out.assign(co.size(), 0.0);
    for (std::size_t i = 0; i < co.size(); ++i) {
      carry_out[i] = co[i].x();
      co[i].seed(i < w.size() ? w[i] : 0.0);
    }
  }

  rec.sweep();

  wz.assign(static_cast<std::size_t>(cp.carry.q + 1) * n, 0.0);
  for (std::size_t i = 0; i < n; ++i) wz[i] = rec.wx()[i];
  for (std::size_t i = 0; i < rec.whistory().size(); ++i)
    wz[n + i] = rec.whistory()[i];

  wp.assign(NP, 0.0);
  rec.accumulate(p, wp);
}

// ---------------------------------------------------------------------------
//  The same step adjoint, written rather than recorded.
//
//  The operators come out of the stepper itself, so this is not a second
//  statement of what a step does; only the transposes and the implicit
//  function theorem are stated here.
// ---------------------------------------------------------------------------

template<cppde::multistep_method M>
static void closed_step(const checkpoint<M>& cp, const std::vector<double>& w,
                        std::vector<double>& wz, std::vector<double>& wp)
{
  const std::size_t n = cp.n_states;

  // The probe never evaluates the right-hand side: the tail only reads it for
  // an order-one restart, which a step-level checkpoint does not record.
  struct null_rhs {
    void operator()(const std::vector<double>&, std::vector<double>& d,
                    const double&) const { d.assign(d.size(), 0.0); }
  };
  struct null_sys { null_rhs first; } nsys;

  // Held across calls, as a sweep holds them: constructing a probe stepper
  // allocates every slot it could need and would otherwise dominate.
  static cppde::adjoint::multistep_probe<stepper_d<M>> probe;
  static cppde::adjoint::multistep_operators<stepper_d<M>> ops;
  static cppde::adjoint::multistep_workspace ws;
  // The tail's own key: what the recorded ops will do to the history.
  const double tail_key = cp.q_next + 1e3 * cp.eta + 1e6 * cp.ops.size();
  probe.build(cp.carry, cp.dt, tail_key,
              [&](stepper_d<M>& pr) { cp.apply_tail(pr, nsys); }, ops);

  std::vector<double> pv(P, P + NP);
  auto jac_d = jacobian<double>{pv};
  cppde::reverse::equation_solver<jacobian<double>, double> solver(jac_d);
  const double t_new = cp.t + ops.h;
  solver.prepare(cp.y, t_new, 1.0 / ops.gamma, ops.gamma);

  adjoint_terms adj{pv};
  std::vector<double> w_out(static_cast<std::size_t>(ops.q_out + 1) * n, 0.0);
  for (std::size_t i = 0; i < w_out.size() && i < w.size(); ++i) w_out[i] = w[i];

  wz.assign(static_cast<std::size_t>(cp.carry.q + 1) * n, 0.0);
  std::vector<double> wphi(NX + NP, 0.0);
  cppde::adjoint::apply_multistep_adjoint<stepper_d<M>>(
      ops, n, NX + NP, cp.y, t_new, w_out.data(), solver, adj,
      wz.data(), wphi.data(), ws);

  wp.assign(NP, 0.0);
  for (std::size_t j = 0; j < NP; ++j) wp[j] = wphi[NX + j];
}

// ---------------------------------------------------------------------------

template<cppde::multistep_method M>
static void compare(const char* name, const checkpoint<M>& cp_in, int q_next,
                    double eta, const std::vector<double>& w,
                    double t_interp = 0.0)
{
  checkpoint<M> cp = cp_in;
  cp.q_next  = q_next;
  cp.eta     = eta;

  const std::size_t n  = cp.n_states;
  const std::size_t nz = static_cast<std::size_t>(cp.carry.q + 1) * n;
  const unsigned    nd = static_cast<unsigned>(nz + NP);

  std::vector<double> S, fwd_carry;
  forward_step<M>(cp, q_next, eta, nd, S, fwd_carry, t_interp);

  std::vector<double> wz, wp, rev_carry;
  reverse_step<M>(cp, w, wz, wp, rev_carry, t_interp);

  // The written adjoint, against the same reference. Not for an observation
  // inside the step: the dense-output adjoint is its own piece and comes with
  // the trajectory, not with the step.
  std::vector<double> cz, cp_par;
  const bool closed = (t_interp <= 0.0);
  if (closed) closed_step<M>(cp, w, cz, cp_par);

  check(fwd_carry.size() == rev_carry.size(),
        std::string(name) + " same carry width");
  for (std::size_t i = 0; i < fwd_carry.size() && i < rev_carry.size(); ++i)
    close(fwd_carry[i], rev_carry[i], std::string(name) + " carry " + std::to_string(i));

  std::printf("%-30s q %d -> %d  eta %.3g  |", name, cp.carry.q, q_next, eta);
  for (std::size_t i = 0; i < wz.size(); ++i) std::printf(" %.17g", wz[i]);
  std::printf("  |");
  for (std::size_t j = 0; j < NP; ++j) std::printf(" %.17g", wp[j]);
  std::printf("\n");

  for (unsigned d = 0; d < nd; ++d) {
    double wS = 0.0;
    for (std::size_t r = 0; r < fwd_carry.size(); ++r)
      wS += (r < w.size() ? w[r] : 0.0) * S[r * nd + d];
    const double got = (d < nz) ? wz[d] : wp[d - nz];
    const std::string tag = (d < nz) ? "  dz" + std::to_string(d)
                                     : "  dp" + std::to_string(d - nz);
    close(wS, got, std::string(name) + tag);
    if (closed) {
      const double gotc = (d < nz) ? cz[d] : cp_par[d - nz];
      close(wS, gotc, std::string(name) + " closed" + tag);
    }
  }
}

template<cppde::multistep_method M>
static void run_method(const char* method)
{
  for (const int warm : {3, 8, 20, 45, 80}) {
    checkpoint<M> cp;
    std::vector<double> x;
    double t = 0.0, dt = 0.0;
    warm_up<M>(cp, x, t, dt, warm);

    const std::size_t n  = cp.n_states;
    const std::size_t nz = static_cast<std::size_t>(cp.carry.q + 1) * n;

    // The state the corrector converged to at this step, from the forward run.
    {
      std::vector<double> p(P, P + NP);
      auto sys = make_system<double>(p);
      stepper_d<M> st;
      st.set_tolerances(ATOL, RTOL);
      st.load_carry(cp.carry, n);
      for (int j = 0; j <= cp.carry.q; ++j)
        for (std::size_t i = 0; i < n; ++i)
          st.zn_mut(j)[i] = cp.zn[static_cast<std::size_t>(j) * n + i];
      std::vector<double> xin(cp.zn.begin(), cp.zn.begin() + n), xo(n), xe(n);
      st.do_step(sys, xin, cp.t, xo, cp.dt, xe);
      check(st.newton_converged(),
            std::string(method) + " the captured step converged");
      polish_corrector<M>(sys, st, cp.t, xo);
      cp.y.assign(xo.begin(), xo.end());
    }

    std::vector<double> wmix(static_cast<std::size_t>(cp.carry.q + 2) * n);
    for (std::size_t k = 0; k < wmix.size(); ++k)
      wmix[k] = 0.4 * static_cast<double>(k % 3) - 0.5;

    const std::string tag = std::string(method) + " warm" + std::to_string(warm);
    check(nz >= 2 * n, tag + " reached order two or higher");

    // Every Nordsieck slot on its own, so no column of the carry Jacobian hides.
    for (std::size_t r = 0; r < nz && r < 4 * n; ++r) {
      std::vector<double> e(nz + n, 0.0);
      e[r] = 1.0;
      compare<M>((tag + " e" + std::to_string(r)).c_str(), cp, cp.carry.q, 1.0, e);
    }

    // Seeded through the dense output rather than on the carry, which is how an
    // observation inside a step reaches the step.
    compare<M>((tag + " interp mid").c_str(), cp, cp.carry.q, 1.0, wmix,
               cp.t + 0.4 * cp.dt);
    compare<M>((tag + " interp end").c_str(), cp, cp.carry.q, 1.0, wmix,
               cp.t + cp.dt);

    compare<M>((tag + " mixed").c_str(),   cp, cp.carry.q,     1.0,  wmix);
    compare<M>((tag + " rescale").c_str(), cp, cp.carry.q,     0.83, wmix);
    compare<M>((tag + " order up").c_str(), cp, cp.carry.q + 1, 1.0, wmix);
    if (cp.carry.q > 1)
      compare<M>((tag + " order down").c_str(), cp, cp.carry.q - 1, 1.0, wmix);
  }
}

int main() {
  std::printf("%-30s %-22s %s\n", "case", "order and rescale",
              "| w' dzn_out/dzn_in | w' dzn_out/dtheta");

  run_method<cppde::multistep_method::bdf>("bdf");
  run_method<cppde::multistep_method::adams>("adams");

  std::printf(g_failures == 0 ? "\nOK\n" : "\n%d FAILURES\n", g_failures);
  return g_failures == 0 ? 0 : 1;
}
