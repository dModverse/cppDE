// A root jump, and nothing else. The forward sandwich over a nested dual is
// the oracle; the backward sweep over a flat one is what is tested. No solver,
// no store, no event engine, so a disagreement is the jump's transpose and
// cannot be anything else.
#include <cppde/cppde.hpp>
#include <cppde/cppde_adjoint_step.hpp>

#include <cstdio>
#include <cmath>
#include <vector>

static constexpr unsigned ND = 5;   // x0, x1, p0, p1, p2
static constexpr std::size_t NS = 2, NP = 3;

// Values. x1 sits on the surface, which is where the engine localises it.
static const double X0 = 1.20, X1 = 0.55, P0 = 0.70, P1 = 0.35, P2 = 0.40;
static const double GC = 0.55;             // root: g = x1 - GC
static const double TE = 1.10;             // event time
static const double W[NS] = {0.83, -0.42}; // the cotangent the sweep starts on
// It depends on theta too, the way one handed down from a later step does.
static const double DW[NS][ND] = {{ 0.11, -0.23,  0.31,  0.07, -0.19},
                                  {-0.05,  0.17, -0.09,  0.26,  0.13}};

// f0 = -p0 x0 + p1 x1
// f1 =  p0 x0 - p1 x1 - p2 x1^2
template<class T>
struct Sys {
  std::vector<T> p;
  template<class St, class Tm>
  void first(const St& x, St& f, const Tm&) const {
    f[0] = -p[0] * x[0] + p[1] * x[1];
    f[1] = p[0] * x[0] - p[1] * x[1] - p[2] * x[1] * x[1];
  }
};

template<class T>
struct Adj {
  const Sys<T>* sys;
  void jac_t_vec(const std::vector<T>& x, const std::vector<T>& lam, double,
                 std::vector<T>& out) const {
    const auto& p = sys->p;
    out.assign(NS, T(0.0));
    for (auto& o : out) cppde::ad_traits::arm_tangents(o);
    out[0] = -p[0] * lam[0] + p[0] * lam[1];
    out[1] = p[1] * lam[0] + (T(0.0) - p[1] - T(2.0) * p[2] * x[1]) * lam[1];
  }
  template<class S>
  void dfdp_t_vec_axpy(const std::vector<T>& x, const std::vector<T>& lam,
                       double, const S& sc, T* w) const {
    w[0] += sc * (T(0.0) - x[0] * lam[0] + x[0] * lam[1]);
    w[1] += sc * (x[1] * lam[0] - x[1] * lam[1]);
    w[2] += sc * (T(0.0) - x[1] * x[1] * lam[1]);
  }
};

// h = p2, so dh/dx is zero and dh/dp is the unit on p2.
template<class T>
struct EvAdj {
  void root_dh_dx(int, const std::vector<T>&, double, std::vector<T>& o) const {
    for (std::size_t i = 0; i < NS; ++i) o[i] = T(0.0);
  }
  void root_dh_dp_axpy(int, const std::vector<T>&, double, const T& sc,
                       T* w) const { w[2] += sc; }
  void root_dg_dp_axpy(int, const std::vector<T>&, double, const T&,
                       T*) const {}
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
  return {e};
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
  using T1 = cppde::dual<double, 0>;
  double grad_rv[ND] = {0}, hess_rv[ND][ND] = {{0}};
  {
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
    Adj<T1> adj{&sys};
    EvAdj<T1> eadj;
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

  const char* nm[ND] = {"x0", "x1", "p0", "p1", "p2"};
  std::printf("gradient        forward        reverse           diff\n");
  double gmax = 0.0;
  for (unsigned a = 0; a < ND; ++a) {
    std::printf("  %-3s %14.8f %14.8f %14.2e\n", nm[a], grad_ff[a], grad_rv[a],
                grad_rv[a] - grad_ff[a]);
    gmax = std::fmax(gmax, std::fabs(grad_rv[a] - grad_ff[a]));
  }
  std::printf("\nHessian, reverse minus forward\n      ");
  for (unsigned b = 0; b < ND; ++b) std::printf("%12s", nm[b]);
  std::printf("\n");
  double hmax = 0.0, href = 0.0;
  for (unsigned a = 0; a < ND; ++a) {
    std::printf("  %-3s ", nm[a]);
    for (unsigned b = 0; b < ND; ++b) {
      const double d = hess_rv[a][b] - hess_ff[a][b];
      std::printf("%12.6f", d);
      hmax = std::fmax(hmax, std::fabs(d));
      href = std::fmax(href, std::fabs(hess_ff[a][b]));
    }
    std::printf("\n");
  }
  std::printf("\nHessian, forward\n      ");
  for (unsigned b = 0; b < ND; ++b) std::printf("%12s", nm[b]);
  std::printf("\n");
  for (unsigned a = 0; a < ND; ++a) {
    std::printf("  %-3s ", nm[a]);
    for (unsigned b = 0; b < ND; ++b) std::printf("%12.6f", hess_ff[a][b]);
    std::printf("\n");
  }
  std::printf("\ngrad %.3e   hess %.3e (relative %.3e)\n",
              gmax, hmax, hmax / (href > 0 ? href : 1.0));
  return 0;
}
