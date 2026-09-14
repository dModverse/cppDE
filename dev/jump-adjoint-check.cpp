// A root jump, and nothing else. The forward sandwich over a nested dual is
// the oracle; the backward sweep over a flat one is what is tested. No solver,
// no store, no event engine, so a disagreement is the jump's transpose and
// cannot be anything else.
#include <cppde/cppde.hpp>
#include <cppde/cppde_adjoint_step.hpp>

#include <cstdio>
#include <cmath>
#include <type_traits>
#include <vector>

static constexpr unsigned ND = 5;   // x0, x1, p0, p1, p2
static constexpr std::size_t NS = 2, NP = 3;

// Values. x1 sits on the surface, which is where the engine localises it.
static const double X0 = 1.20, X1 = 0.55, P0 = 0.70, P1 = 0.35, P2 = 0.40;
static const double GC = 0.55;             // root: g = x1 - GC
static const double TD = 0.30;             // how hard f1 reads the clock
static const double TE = 1.10;             // event time
static const double W[NS] = {0.83, -0.42}; // the cotangent the sweep starts on
// It depends on theta too, the way one handed down from a later step does.
static const double DW[NS][ND] = {{ 0.11, -0.23,  0.31,  0.07, -0.19},
                                  {-0.05,  0.17, -0.09,  0.26,  0.13}};

// f0 = -p0 x0 + p1 x1
// f1 =  p0 x0 - p1 x1 - p2 x1^2 - TD t x1
template<class T>
struct Sys {
  std::vector<T> p;
  template<class St, class Tm>
  void first(const St& x, St& f, const Tm& t) const {
    f[0] = -p[0] * x[0] + p[1] * x[1];
    f[1] = p[0] * x[0] - p[1] * x[1] - p[2] * x[1] * x[1] - T(TD) * T(t) * x[1];
  }
};

struct no_scope { };
template<bool On>
using maybe_scope = std::conditional_t<On, cppde::dual_arena::scope, no_scope>;

// Two shapes of the same contraction. Generated = false arms its output and
// opens no scope; Generated = true is what the codegen emits, an unarmed
// assign followed by the model's own arena scope. The buffer binding is the
// one thing a hand-written oracle gets wrong by being careful.
template<class T, bool Generated>
struct Adj {
  const Sys<T>* sys;
  mutable int calls = 0;

  void jac_t_vec(const std::vector<T>& x, const std::vector<T>& lam, double t,
                 std::vector<T>& out) const {
    const auto& p = sys->p;
    auto& A = cppde::dual_arena::arena();
    const std::size_t top_in = A.slabs_[A.active_].top;
    out.assign(NS, T(0.0));
    if constexpr (!Generated)
      for (auto& o : out) cppde::ad_traits::arm_tangents(o);
    maybe_scope<Generated> _rhs_arena_scope;
    out[0] = -p[0] * lam[0] + p[0] * lam[1];
    out[1] = p[1] * lam[0]
           + (T(0.0) - p[1] - T(2.0) * p[2] * x[1] - T(TD) * T(t)) * lam[1];
    if constexpr (Generated) {
      // A tangent below the top this call entered at is storage the caller has
      // since handed to someone else.
      const char* base = (const char*) A.slabs_[A.active_].data;
      std::printf("  jac_t_vec #%d  top at entry %zu  out0 @%td  out1 @%td\n",
                  ++calls, top_in,
                  (const char*) &out[0].d(0) - base,
                  (const char*) &out[1].d(0) - base);
    }
  }
  template<class S>
  void dfdp_t_vec_axpy(const std::vector<T>& x, const std::vector<T>& lam,
                       double, const S& sc, T* w) const {
    maybe_scope<Generated> _rhs_arena_scope;
    w[0] += sc * (T(0.0) - x[0] * lam[0] + x[0] * lam[1]);
    w[1] += sc * (x[1] * lam[0] - x[1] * lam[1]);
    w[2] += sc * (T(0.0) - x[1] * x[1] * lam[1]);
  }
  // f1 reads the clock through -TD t x1, so df/dt is not zero.
  T dfdt_dot(const std::vector<T>& x, const std::vector<T>& lam,
             const T&) const {
    return (T(0.0) - T(TD)) * x[1] * lam[1];
  }
};

// h = p2, so dh/dx is zero and dh/dp is the unit on p2.
template<class T>
struct EvAdj {
  const Sys<T>* sys;
  void root_dh_dx(int, const std::vector<T>&, const T&, std::vector<T>& o) const {
    for (std::size_t i = 0; i < NS; ++i) o[i] = T(0.0);
  }
  void root_dh_dp_axpy(int, const std::vector<T>&, const T&, const T& sc,
                       T* w) const { w[2] += sc; }
  void root_dg_dp_axpy(int, const std::vector<T>&, double, const T&,
                       T*) const {}
  // h = p2 reads no clock.
  void root_dh_dt_axpy(int, const std::vector<T>&, const T&, const T&,
                       T*) const {}
  // g = x1 - GC, so g_dot = f1 and grad g_dot is row 1 of J on the state and
  // of df/dp on the parameters. Written out, not assembled, as the model does.
  void root_gdot_dx(int, const std::vector<T>& x, double t,
                    std::vector<T>& o) const {
    const auto& p = sys->p;
    o.assign(NS, T(0.0));
    o[0] = p[0];
    o[1] = T(0.0) - p[1] - T(2.0) * p[2] * x[1] - T(TD) * T(t);
  }
  void root_gdot_dp_axpy(int, const std::vector<T>& x, double, const T& sc,
                         T* w) const {
    w[0] += sc * x[0];
    w[1] += sc * (T(0.0) - x[1]);
    w[2] += sc * (T(0.0) - x[1] * x[1]);
  }
};

template<class T>
std::vector<cppde::detail::RootEvent<std::vector<T>, T> > make_events() {
  cppde::detail::RootEvent<std::vector<T>, T> e;
  e.func = [](const std::vector<T>& x, const T&) -> T {
    T r; r = x[1] - T(GC); return r;
  };
  e.state_index = 0;
  e.value_func = [](const std::vector<T>&, const T&) -> T { return T(P2); };
  e.method = cppde::detail::EventMethod::Add;
  e.terminal = false;
  e.dg_dx = [](const std::vector<T>&, const T&, std::vector<T>& g) {
    g[0] = T(0.0); g[1] = T(1.0);
  };
  e.dg_dt = [](const std::vector<T>&, const T&) -> T { return T(0.0); };
  // G_tt along the flow. Present, so root_ift_curvature takes the analytical
  // branch and allocates nothing, which is what a model without forcings does.
  // On the finite-difference branch the temporaries push the arena top past
  // the window this check is about.
  e.g_dot_dot = [](const std::vector<T>& x, const T& t) -> double {
    const double x0 = cppde::ad_traits::scalar_value(x[0]);
    const double x1 = cppde::ad_traits::scalar_value(x[1]);
    const double tv = cppde::ad_traits::scalar_value(t);
    const double f0 = -P0 * x0 + P1 * x1;
    const double f1 = P0 * x0 - P1 * x1 - P2 * x1 * x1 - TD * tv * x1;
    // g_ddot = df1/dt + grad f1 . f
    return -TD * x1 + P0 * f0 + (-P1 - 2.0 * P2 * x1 - TD * tv) * f1;
  };
  return {e};
}

static const char* const NAMES[ND] = {"x0", "x1", "p0", "p1", "p2"};

static void report(const char* what,
                   const double* grad_ff, const double (*hess_ff)[ND],
                   const double* grad_rv, const double (*hess_rv)[ND]) {
  std::printf("\n=== %s ===\n", what);
  std::printf("gradient        forward        reverse           diff\n");
  double gmax = 0.0;
  for (unsigned a = 0; a < ND; ++a) {
    std::printf("  %-3s %14.8f %14.8f %14.2e\n", NAMES[a], grad_ff[a],
                grad_rv[a], grad_rv[a] - grad_ff[a]);
    gmax = std::fmax(gmax, std::fabs(grad_rv[a] - grad_ff[a]));
  }
  std::printf("\nHessian, reverse minus forward\n      ");
  for (unsigned b = 0; b < ND; ++b) std::printf("%12s", NAMES[b]);
  std::printf("\n");
  double hmax = 0.0, href = 0.0;
  for (unsigned a = 0; a < ND; ++a) {
    std::printf("  %-3s ", NAMES[a]);
    for (unsigned b = 0; b < ND; ++b) {
      const double d = hess_rv[a][b] - hess_ff[a][b];
      std::printf("%12.6f", d);
      hmax = std::fmax(hmax, std::fabs(d));
      href = std::fmax(href, std::fabs(hess_ff[a][b]));
    }
    std::printf("\n");
  }
  std::printf("\ngrad %.3e   hess %.3e (relative %.3e)\n",
              gmax, hmax, hmax / (href > 0 ? href : 1.0));
}

// The sandwich forward over a flat dual, then its transpose. Generated picks
// which buffer shape the contractions present to the sweep.
template<bool Generated>
static void run_reverse(const std::vector<cppde::detail::TriggeredEvent>& trig,
                        double* grad_rv, double (*hess_rv)[ND]) {
  using T1 = cppde::dual<double, 0>;
  Sys<T1> sys;
  sys.p.resize(NP);
  const double pv[NP] = {P0, P1, P2};
  std::vector<T1> xb(NS);
  const double xv[NS] = {X0, X1};
  for (std::size_t i = 0; i < NS; ++i) {
    xb[i].diff((unsigned)i, ND); xb[i].x() = xv[i];
  }
  for (std::size_t j = 0; j < NP; ++j) {
    sys.p[j].diff((unsigned)(NS + j), ND); sys.p[j].x() = pv[j];
  }
  auto ev = make_events<T1>();
  T1 amt = sys.p[2];
  ev[0].value_func = [amt](const std::vector<T1>&, const T1&) -> T1 { return amt; };

  std::vector<T1> xa(NS);
  for (std::size_t i = 0; i < NS; ++i) { xa[i].arm(); xa[i] = xb[i]; }
  T1 te; te.arm(); te.x() = TE;
  cppde::detail::saltation_root_analytical_batch(xa, xb, te, sys, ev, trig);

  cppde::dual_arena::width_scope _w(ND);
  Adj<T1, Generated> adj{&sys};
  EvAdj<T1> eadj{&sys};
  cppde::adjoint::jump_workspace<T1> ws;
  std::vector<T1> wout(NS), win(NS), wp(NP);
  cppde::adjoint::zero_armed(wout, NS);
  cppde::adjoint::zero_armed(win, NS);
  cppde::adjoint::zero_armed(wp, NP);
  for (std::size_t i = 0; i < NS; ++i) {
    wout[i].diff(0u, ND); wout[i].x() = W[i];
    for (unsigned b = 0; b < ND; ++b) wout[i].d(b) = DW[i][b];
  }

  cppde::adjoint::apply_root_jump_adjoint(
      xb, xa, TE, ev, trig, sys, eadj, adj, NS,
      wout.data(), win.data(), wp.data(), ws);

  for (std::size_t i = 0; i < NS; ++i) {
    grad_rv[i] = win[i].x();
    for (unsigned b = 0; b < ND; ++b) hess_rv[i][b] = win[i].d(b);
  }
  for (std::size_t j = 0; j < NP; ++j) {
    grad_rv[NS + j] = wp[j].x();
    for (unsigned b = 0; b < ND; ++b) hess_rv[NS + j][b] = wp[j].d(b);
  }
}

int main() {
  const std::vector<cppde::detail::TriggeredEvent> trig{{0, -1.0, 1.0}};

  // ---- oracle: the sandwich over a nested dual -----------------------------
  using T2 = cppde::dual2nd<double, 0>;
  double grad_ff[ND] = {0}, hess_ff[ND][ND] = {{0}};
  {
    Sys<T2> sys;
    sys.p.resize(NP);
    const double pv[NP] = {P0, P1, P2};
    std::vector<T2> xb(NS);
    const double xv[NS] = {X0, X1};
    for (std::size_t i = 0; i < NS; ++i) {
      xb[i].arm_full(ND); xb[i].scalar() = xv[i]; xb[i].d1_at((unsigned)i) = 1.0;
    }
    for (std::size_t j = 0; j < NP; ++j) {
      sys.p[j].arm_full(ND); sys.p[j].scalar() = pv[j];
      sys.p[j].d1_at((unsigned)(NS + j)) = 1.0;
    }
    // The event's height is the parameter p2, seeded like the others.
    auto ev = make_events<T2>();
    T2 amt = sys.p[2];
    ev[0].value_func = [amt](const std::vector<T2>&, const T2&) -> T2 { return amt; };

    std::vector<T2> xo(NS);
    for (std::size_t i = 0; i < NS; ++i) { xo[i].arm_full(ND); xo[i] = xb[i]; }
    T2 te; te.arm_full(ND); te.scalar() = TE;
    cppde::detail::saltation_root_analytical_batch(xo, xb, te, sys, ev, trig);

    // What the sweep returns is g_a = sum_i w_i dx_i/dtheta_a with w an input,
    // so its tangent in direction b carries dw/dtheta_b against dx/dtheta_a and
    // w against the second derivative. a is the route, b the direction.
    for (unsigned a = 0; a < ND; ++a) {
      for (std::size_t i = 0; i < NS; ++i) grad_ff[a] += W[i] * xo[i].d1_at(a);
      for (unsigned b = 0; b < ND; ++b)
        for (std::size_t i = 0; i < NS; ++i)
          hess_ff[a][b] += DW[i][b] * xo[i].d1_at(a) + W[i] * xo[i].dd_at(a, b);
    }
  }

  // ---- test: the same sandwich forward, then its transpose -----------------
  double grad_rv[ND] = {0}, hess_rv[ND][ND] = {{0}};
  double grad_gn[ND] = {0}, hess_gn[ND][ND] = {{0}};
  run_reverse<false>(trig, grad_rv, hess_rv);
  std::printf("arena, under the buffer shape the codegen emits\n");
  run_reverse<true>(trig, grad_gn, hess_gn);
  std::printf("\n");

  report("output buffers armed by the caller", grad_ff, hess_ff,
         grad_rv, hess_rv);
  report("output buffers bound the way the codegen binds them", grad_ff,
         hess_ff, grad_gn, hess_gn);

  std::printf("\nHessian, forward\n      ");
  for (unsigned b = 0; b < ND; ++b) std::printf("%12s", NAMES[b]);
  std::printf("\n");
  for (unsigned a = 0; a < ND; ++a) {
    std::printf("  %-3s ", NAMES[a]);
    for (unsigned b = 0; b < ND; ++b) std::printf("%12.6f", hess_ff[a][b]);
    std::printf("\n");
  }
  return 0;
}
