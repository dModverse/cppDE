/**
 * @file cppde_adjoint_step.hpp
 * @brief The adjoint of one accepted step, written out per stepper family.
 *
 * The accepted grid is frozen, so a step is a fixed map and its adjoint is
 * stated. What a multistep step does to the Nordsieck history splits in three:
 *
 *   zn_pred = A zn_in                      rescale and Pascal shift
 *   res(y, zn_pred, theta) = 0             the corrector, closed by the IFT
 *   zn_out  = B zn_pred + c acor           the tail, acor = y - zn_pred[0]
 *
 * A, B and c are read off the stepper by running its own routines on a unit
 * slot. The model supplies `jac_t_vec` (J' lambda) and `dfdp_t_vec_axpy`
 * ((df/dp)' lambda, scaled and added). See vignette("Methods"), section
 * "The step map and its adjoint".
 *
 * Copyright (C) 2026 Simon Beyer
 */

#ifndef CPPDE_ADJOINT_STEP_HPP
#define CPPDE_ADJOINT_STEP_HPP

#include <cmath>
#include <cstddef>
#include <limits>
#include <type_traits>
#include <vector>

#include <cppde/cppde_dual_math.hpp>
#include <cppde/cppde_events.hpp>
#include <cppde/cppde_saltation.hpp>
#include <cppde/cppde_profiler.hpp>
#include <cppde/cppde_reverse_step.hpp>

// The hot loops here run over buffers that are always distinct. Saying so is
// what lets them vectorise.
#if defined(__GNUC__) || defined(__clang__)
#  define CPPDE_RESTRICT __restrict__
#elif defined(_MSC_VER)
#  define CPPDE_RESTRICT __restrict
#else
#  define CPPDE_RESTRICT
#endif

namespace cppde {
namespace adjoint {

// ---------------------------------------------------------------------------
//  The forward loop observes a time before the step bracket at the bracket
//  start instead, which happens after an event restart. Clamping reproduces
//  that branch; inside the bracket, which is every other case, it does nothing.
// ---------------------------------------------------------------------------
template<class Checkpoint>
inline double clamp_to_step(const Checkpoint& cp, double t) {
  const double a = cp.t, b = cp.t + cp.dt;
  const double lo = a < b ? a : b;
  const double hi = a < b ? b : a;
  return t < lo ? lo : (t > hi ? hi : t);
}

// ---------------------------------------------------------------------------
//  n zeros with their tangent storage bound, so a callee that opens its own
//  dual_arena::scope writes into them in place. Sized to exactly what the
//  callee writes: a later growth copy-constructs the new elements unarmed.
// ---------------------------------------------------------------------------
template<class T>
inline void zero_armed(std::vector<T>& v, std::size_t n) {
  v.assign(n, T(0.0));
  for (std::size_t i = 0; i < n; ++i) cppde::ad_traits::arm_tangents(v[i]);
}

// ---------------------------------------------------------------------------
//  A dense operator over Nordsieck slots, row-major [n_out x n_in].
//
//  Small by construction: the order never exceeds the stepper's own maximum,
//  so this is at most 14 by 14 doubles.
// ---------------------------------------------------------------------------
struct slot_operator {
  int rows = 0, cols = 0;
  std::vector<double> a;

  void resize(int r, int c) { rows = r; cols = c; a.assign(static_cast<std::size_t>(r) * c, 0.0); }
  double&       operator()(int i, int j)       { return a[static_cast<std::size_t>(i) * cols + j]; }
  double        operator()(int i, int j) const { return a[static_cast<std::size_t>(i) * cols + j]; }

  /// out[j*n + i] += sum_k this(k, j) * w[k*n + i], the transposed apply.
  /// The two arrays are always distinct buffers, which the compiler cannot see
  /// and which decides whether the inner loop vectorises.
  template<class T>
  void apply_transposed(const T* CPPDE_RESTRICT w, std::size_t n,
                        T* CPPDE_RESTRICT out) const {
    for (int j = 0; j < cols; ++j)
      for (int k = 0; k < rows; ++k) {
        const double m = (*this)(k, j);
        if (m == 0.0) continue;
        const T* CPPDE_RESTRICT wk = w + static_cast<std::size_t>(k) * n;
        T* CPPDE_RESTRICT oj = out + static_cast<std::size_t>(j) * n;
        for (std::size_t i = 0; i < n; ++i) oj[i] += m * wk[i];
      }
  }
};

// ---------------------------------------------------------------------------
//  The three operators of one multistep step, read off the stepper.
//
//  `probe` is a one-state copy of the stepper the run used. It is driven
//  through the very same routines with a unit slot in place of the state, so
//  every column is what the forward step would have done to that slot.
// ---------------------------------------------------------------------------
template<class Stepper>
struct multistep_operators {
  slot_operator A;      ///< zn_in -> zn_pred
  slot_operator B;      ///< zn_pred -> zn_out
  std::vector<double> c;///< the acor column of the tail
  double rl1 = 0.0;     ///< the residual's coefficient on zn_pred[1]
  double gamma = 0.0;   ///< and on f
  double h = 0.0;
  double err_scale = 1.0;  ///< what acor has to be multiplied by to be the error
  int q_in = 1, q_out = 1;
};

// ---------------------------------------------------------------------------
//  Building the operators. They act on the slot index identically for every
//  state component, so a probe whose states are the slots gives every column
//  in one pass: put the identity in and read the result back as a matrix.
// ---------------------------------------------------------------------------

/// A: what the rescale and the Pascal shift do, plus the step's coefficients.
template<class Stepper, class Carry>
void probe_pre(Stepper& probe, const Carry& carry, double dt, slot_operator& A,
               double& rl1, double& gamma, double& h, double* err_scale = nullptr)
{
  const int q = carry.q;
  const std::size_t m = static_cast<std::size_t>(q + 1);

  probe.load_carry(carry, m);
  for (int j = 0; j <= Stepper::max_order; ++j) {
    auto& slot = probe.zn_mut(j);
    for (std::size_t i = 0; i < m; ++i)
      slot[i] = (j == static_cast<int>(i)) ? 1.0 : 0.0;
  }
  rl1   = static_cast<double>(probe.replay_predict(m, dt));
  gamma = static_cast<double>(probe.gamma());
  h     = static_cast<double>(probe.h());
  // Set by the same call that set the coefficients, and only meaningful after.
  if (err_scale) *err_scale = static_cast<double>(probe.error_constant());

  A.resize(q + 1, q + 1);
  for (int j = 0; j <= q; ++j)
    for (int k = 0; k <= q; ++k) A(j, k) = static_cast<double>(probe.zn(j)[k]);
}

/// B and the acor column c: what the tail does. Component q+1 holds acor.
template<class Stepper, class Carry, class TailFn>
void probe_tail(Stepper& probe, const Carry& carry, double dt, TailFn&& tail,
                slot_operator& B, std::vector<double>& c, int& q_out)
{
  const int q = carry.q;
  const std::size_t m = static_cast<std::size_t>(q + 2);

  probe.load_carry(carry, m);
  for (int j = 0; j <= Stepper::max_order; ++j) probe.zn_mut(j).assign(m, 0.0);
  probe.replay_predict(m, dt);
  // The tail reads the predicted array, so the identity goes in after the
  // shift, not before it.
  for (int j = 0; j <= Stepper::max_order; ++j) {
    auto& slot = probe.zn_mut(j);
    slot.assign(m, 0.0);
    if (j <= q) slot[static_cast<std::size_t>(j)] = 1.0;
  }

  // acor = y - zn_pred[0]. Zero for the columns of B, one for the last
  // component, which is therefore the acor column.
  std::vector<typename Stepper::value_type> y(m, 0.0), out(m), err(m);
  y[0] = 1.0;              // component 0 holds slot 0's unit
  y[m - 1] = 1.0;
  probe.replay_outputs(y, out, err);
  tail(probe);

  q_out = probe.current_order();
  B.resize(q_out + 1, q + 1);
  c.assign(static_cast<std::size_t>(q_out + 1), 0.0);
  for (int j = 0; j <= q_out; ++j) {
    const auto& slot = probe.zn(j);
    for (int k = 0; k <= q; ++k) B(j, k) = static_cast<double>(slot[k]);
    c[static_cast<std::size_t>(j)] = static_cast<double>(slot[m - 1]);
  }
}

/// The two probe steppers, kept across a sweep. Constructing one allocates
/// every Nordsieck slot it could ever need, which on a step budget is the
/// whole cost of reading the operators off.
template<class Stepper>
struct multistep_probe {
  // The operators are polynomials in the step-size history, which is a double
  // whatever the run's scalar type is. So the probe is a double stepper even
  // under a sweep that propagates tangents, and the matrices come out double.
  using probe_stepper = typename Stepper::template rebind_value<double>;
  probe_stepper pre, tail_probe;

  // What the operators depend on; a step whose key matches reuses them. qwait
  // is left out: it reaches only the controller's order-change constants,
  // never l, gamma or the Nordsieck shift.
  struct key {
    int q = -1, L = 0;
    bool started = false;
    double dt = 0.0, h = 0.0, hscale = 0.0, eta = 0.0;
    double tail = 0.0;
    std::vector<double> tau;
    // The recorded tail itself: two tails of the same length can rescale by
    // different factors, and the operators differ with them.
    std::vector<double> ops;
    bool operator==(const key& o) const {
      return q == o.q && L == o.L && started == o.started &&
             dt == o.dt && h == o.h && hscale == o.hscale && eta == o.eta &&
             tail == o.tail && tau == o.tau && ops == o.ops;
    }
  };
  key last, pending;
  bool valid = false;

  /// Whether the operators already in hand are this step's. The key it built
  /// to answer stays, so the rebuild that may follow does not build it twice.
  template<class Carry>
  bool matches(const Carry& carry, double dt, double tail_key,
               const history_log* ops = nullptr) {
    pending.ops.clear();
    if (ops)
      for (const history_entry& e : *ops) {
        pending.ops.push_back(static_cast<double>(e.op));
        pending.ops.push_back(e.value);
      }
    pending.q = carry.q; pending.L = carry.L;
    pending.started = carry.nst > 0;
    pending.dt = dt; pending.h = carry.h;
    pending.hscale = carry.hscale; pending.eta = carry.eta;
    pending.tail = tail_key;
    pending.tau.assign(carry.tau.begin(), carry.tau.begin() + (carry.q + 1));
    return valid && pending == last;
  }

  /// Read A, B and c off the stepper. `tail` applies whatever the run recorded
  /// between this step's acceptance and the next one's. Call after a matches()
  /// that said no, whose key this commits.
  template<class Carry, class TailFn>
  void rebuild(const Carry& carry, double dt, TailFn&& tail,
               multistep_operators<Stepper>& ops)
  {
    ops.q_in = carry.q;
    probe_pre(pre, carry, dt, ops.A, ops.rl1, ops.gamma, ops.h, &ops.err_scale);
    probe_tail(tail_probe, carry, dt, tail, ops.B, ops.c, ops.q_out);
    last = pending;
    valid = true;
  }

  /// Both at once, for a caller with nothing to time between them.
  template<class Carry, class TailFn>
  void build(const Carry& carry, double dt, double tail_key, TailFn&& tail,
             multistep_operators<Stepper>& ops)
  {
    if (!matches(carry, dt, tail_key)) rebuild(carry, dt, tail, ops);
  }
};

/// One-shot form, for a caller that has no sweep to hang a probe on.
template<class Stepper, class Carry, class TailFn>
multistep_operators<Stepper> build_multistep_operators(
    const Carry& carry, double dt, TailFn&& tail)
{
  multistep_operators<Stepper> ops;
  multistep_probe<Stepper> probe;
  probe.build(carry, dt, 0.0, tail, ops);
  return ops;
}

// ---------------------------------------------------------------------------
//  The adjoint of one step: w_pred = B' w_out, mu = (I - gamma J)^-T w_y for
//  the corrector, w_in = A' w_pred. See vignette("Methods"), "The step map and
//  its adjoint". `solver`'s transposed apply includes the gamma scale.
// ---------------------------------------------------------------------------
/// What the step above hands down, passed back through the tail onto the
/// predicted history and onto acor. Separate from the step adjoint itself
/// because a trajectory adds its observations to the same two.
template<class Stepper, class T>
void carry_into_pred(const multistep_operators<Stepper>& ops, std::size_t n,
                     const T* CPPDE_RESTRICT w_out,
                     T* CPPDE_RESTRICT w_pred,
                     T* CPPDE_RESTRICT w_acor)
{
  ops.B.apply_transposed(w_out, n, w_pred);
  for (int k = 0; k <= ops.q_out; ++k) {
    const double m = ops.c[static_cast<std::size_t>(k)];
    if (m == 0.0) continue;
    const T* CPPDE_RESTRICT wk = w_out + static_cast<std::size_t>(k) * n;
    for (std::size_t i = 0; i < n; ++i) w_acor[i] += m * wk[i];
  }
}

/// The buffers one step adjoint needs, kept across a sweep so a step allocates
/// nothing.
template<class T = double>
struct multistep_workspace {
  std::vector<T> w_pred, w_acor, mu, x, q;
};

/// The step adjoint's first half: the right-hand side of its transposed solve,
/// left in ws.mu. Split from the second so a caller can time the solve apart.
template<class Stepper, class T>
void multistep_adjoint_rhs(const multistep_operators<Stepper>& ops,
                           std::size_t n, multistep_workspace<T>& ws,
                           T* w_out_state = nullptr)
{
  T* CPPDE_RESTRICT w_pred = ws.w_pred.data();
  const T* CPPDE_RESTRICT w_acor = ws.w_acor.data();

  ws.mu = ws.w_acor;
  for (std::size_t i = 0; i < n; ++i) w_pred[i] -= w_acor[i];
  if (w_out_state) for (std::size_t i = 0; i < n; ++i) w_out_state[i] = w_acor[i];
  (void)ops;
}

/// The second half, on a ws.mu the caller has already solved with.
template<class Stepper, class AdjTerms, class T>
void multistep_adjoint_finish(const multistep_operators<Stepper>& ops,
                              std::size_t n, std::size_t n_phi,
                              const std::vector<T>& y, double t_new,
                              const AdjTerms& adj,
                              T* w_in, T* w_theta,
                              multistep_workspace<T>& ws)
{
  T* CPPDE_RESTRICT w_pred = ws.w_pred.data();
  const T* CPPDE_RESTRICT mu = ws.mu.data();

  for (std::size_t i = 0; i < n; ++i) w_pred[i] += mu[i];
  if (ops.q_in >= 1) {
    const double rl1 = ops.rl1;
    for (std::size_t i = 0; i < n; ++i) w_pred[n + i] -= rl1 * mu[i];
  }

  if (w_theta) adj.dfdp_t_vec_axpy(y, ws.mu, t_new, ops.gamma, w_theta);
  ops.A.apply_transposed(ws.w_pred.data(), n, w_in);
  (void)n_phi;
}

/// The step adjoint proper, on a predicted-history cotangent that the caller
/// has already assembled.
template<class Stepper, class Solver, class AdjTerms, class T>
void apply_multistep_adjoint_pre(const multistep_operators<Stepper>& ops,
                                 std::size_t n, std::size_t n_phi,
                                 const std::vector<T>& y, double t_new,
                                 Solver& solver, const AdjTerms& adj,
                                 T* w_in, T* w_theta,
                                 multistep_workspace<T>& ws,
                                 T* w_out_state = nullptr)
{
  multistep_adjoint_rhs(ops, n, ws, w_out_state);
  solver.transposed(ws.mu);
  multistep_adjoint_finish(ops, n, n_phi, y, t_new, adj, w_in, w_theta, ws);
}

template<class Stepper, class Solver, class AdjTerms, class T>
void apply_multistep_adjoint(const multistep_operators<Stepper>& ops,
                             std::size_t n, std::size_t n_phi,
                             const std::vector<T>& y, double t_new,
                             const T* w_out,
                             Solver& solver, const AdjTerms& adj,
                             T* w_in, T* w_theta,
                             multistep_workspace<T>& ws,
                             T* w_out_state = nullptr)
{
  const std::size_t nz_in = static_cast<std::size_t>(ops.q_in + 1) * n;
  std::vector<T>& w_pred = ws.w_pred;
  std::vector<T>& w_acor = ws.w_acor;
  w_pred.assign(nz_in, T(0.0));
  w_acor.assign(n, T(0.0));
  carry_into_pred(ops, n, w_out, w_pred.data(), w_acor.data());
  apply_multistep_adjoint_pre(ops, n, n_phi, y, t_new, solver, adj,
                              w_in, w_theta, ws, w_out_state);
}

/// What a trajectory needs to take a jump apart: the right-hand side, the
/// model's derivatives of the event expressions, and the events themselves.
/// A model without any passes `no_jumps`, and the boundary compiles out.
struct no_jumps { static constexpr bool active = false; };

template<class System, class EvAdj, class FixedEvents, class RootEvents>
struct jump_terms {
  static constexpr bool active = true;
  System& sys;
  const EvAdj& eadj;
  const FixedEvents& fixed;
  const RootEvents& root;
};

template<class System, class EvAdj, class FixedEvents, class RootEvents>
jump_terms<System, EvAdj, FixedEvents, RootEvents>
make_jump_terms(System& s, const EvAdj& e, const FixedEvents& f,
                const RootEvents& r) { return {s, e, f, r}; }

/// The restart, transposed. initialize() builds the whole Nordsieck history out
/// of one state, zn[0] = x and zn[1] = h f(x, t) with the rest zero, so the
/// history's cotangent collapses onto that state and nothing is left above it.
/// The trajectory start and every event boundary apply the same map.
template<class AdjTerms, class T>
void collapse_restart(const std::vector<T>& w_carry, std::size_t n,
                      const std::vector<T>& x0, double t0, double h0,
                      const AdjTerms& adj, T* w_state, T* w_theta,
                      std::vector<T>& mu, std::vector<T>& jv)
{
  for (std::size_t i = 0; i < n; ++i)
    w_state[i] = (i < w_carry.size()) ? w_carry[i] : T(0.0);
  if (w_carry.size() < 2 * n) return;
  mu.assign(w_carry.begin() + static_cast<std::ptrdiff_t>(n),
            w_carry.begin() + static_cast<std::ptrdiff_t>(2 * n));
  adj.jac_t_vec(x0, mu, t0, jv);
  for (std::size_t i = 0; i < n; ++i) w_state[i] += h0 * jv[i];
  adj.dfdp_t_vec_axpy(x0, mu, t0, h0, w_theta);
}

// ---------------------------------------------------------------------------
//  The adjoint of one jump. In value the Heun sandwich collapses to the
//  saltation dx = R'(dx_before + f_before s) - f_after s, with s the event
//  time's differential. See vignette("Methods"), "The adjoint of an event".
// ---------------------------------------------------------------------------
template<class T = double>
struct jump_workspace {
  std::vector<T> fb, fa, g, wy;
  std::vector<std::vector<T> > path;
  // The reset chain on one sandwich's surface, rebuilt per event.
  std::vector<std::vector<T> > micro;
};

/// One reset, transposed. `w_z` is the cotangent on the state it wrote, `w_y`
/// the one on the state it read, which starts as a copy of `w_z`.
///     Replace   z[k] = h(y, t)
///     Add       z[k] = y[k] + h(y, t)
///     Multiply  z[k] = y[k] * h(y, t)
/// `w_time` takes what the reset contributes to the event time's cotangent.
template<class GradX, class GradP, class GradT, class T>
void reset_transpose(int k, cppde::detail::EventMethod method, const T& h,
                     int idx,
                     const std::vector<T>& y, const T& t, std::size_t n,
                     const T* w_z, T* w_y, T* w_theta,
                     GradX&& dh_dx, GradP&& dh_dp_axpy, std::vector<T>& g,
                     GradT&& dh_dt_axpy, T& w_time)
{
  if (k < 0) return;
  const T wk = w_z[k];
  T c = T(1.0);
  using cppde::detail::EventMethod;
  switch (method) {
    case EventMethod::Replace:  w_y[k] -= wk; break;
    case EventMethod::Add:      break;
    case EventMethod::Multiply: w_y[k] += wk * (h - T(1.0)); c = y[k]; break;
  }
  // A zero cotangent contributes nothing, and asking is not free of the
  // scalar type: under tangents the value alone does not say so.
  zero_armed(g, n);
  dh_dx(idx, y, t, g);
  for (std::size_t i = 0; i < n; ++i) w_y[i] += wk * c * g[i];
  dh_dp_axpy(idx, y, t, wk * c, w_theta);
  // A reset that reads the clock moves with the event time, so dh/dt belongs
  // on the scalar the shift is scaled by.
  dh_dt_axpy(idx, y, t, wk * c, &w_time);
}

/// A batch of root events, which share one surface and one dt*. f reads the
/// modes `mb` on the near side of the surface and `ma` on the far side.
template<class System, class RootEvents, class EvAdj, class AdjTerms, class T>
void apply_root_jump_adjoint(const std::vector<T>& x_before,
                             const std::vector<T>& x_after,
                             double t, const RootEvents& root_events,
                             const std::vector<cppde::detail::TriggeredEvent>& triggered,
                             System& sys, const EvAdj& eadj,
                             const AdjTerms& adj, std::size_t n,
                             const T* w_out, T* w_in, T* w_theta,
                             jump_workspace<T>& ws,
                             const std::vector<signed char>& mb = {},
                             const std::vector<signed char>& ma = {})
{
  using cppde::detail::use_switch_modes;
  use_switch_modes(mb);
  zero_armed(ws.fb, n);
  sys.first(x_before, ws.fb, t);

  // Which event dt* came from, picked the way the forward run picks it: the
  // first triggered, non-terminal one whose gradients are there.
  std::size_t src = triggered.size();
  T g_dot = T(0.0);
  for (std::size_t j = 0; j < triggered.size(); ++j) {
    const auto& evt = root_events[triggered[j].index];
    if (evt.terminal) continue;
    if (evt.dg_dx && evt.dg_dt) {
      zero_armed(ws.g, n);
      evt.dg_dx(x_before, t, ws.g);
      T gd = evt.dg_dt(x_before, t);
      for (std::size_t i = 0; i < n; ++i) gd += ws.g[i] * ws.fb[i];
      if (std::abs(ad_traits::scalar_value(gd)) >= 1e-15) { src = j; g_dot = gd; }
      break;
    }
  }
  const bool shifted = (src < triggered.size());

  const T t_ad = T(t);

  // The shift, s = -(g - value(g)) / g_dot: zero in value, dt*/dv in the
  // tangent. It reads off the stored state because the localisation is
  // tangent-blind, so that state holds dx/dv at a fixed time and g contracts
  // it to the event time's. Zero without gradients, which collapses everything
  // below to the bare reset transpose.
  T s = T(0.0);
  if (shifted) {
    const T gv = root_events[triggered[src].index].func(x_before, t_ad);
    s = (T(ad_traits::scalar_value(gv)) - gv) / g_dot;
  }
  const T t_s = t_ad + s;

  // Forward: x_e = x_b + f_b s, x_* = x_b + (f_b + f_e) s / 2, x_a = R(x_*),
  // x_k = x_a - f_a s = the stored state, x_out = x_a - (f_a + f_k) s / 2.
  // Transposed below; see vignette("Methods"), "Transposing the root sandwich".
  use_switch_modes(ma);
  zero_armed(ws.fa, n);
  sys.first(x_after, ws.fa, t);
  use_switch_modes(mb);
  std::vector<T> xe, xst, xa, xk, fe, fa, fk;
  std::vector<T> wfa, wfk, wfb, wfe, wxb, wv, wa;
  for (auto* v : {&xe, &xst, &xa, &xk, &fe, &fa, &fk,
                  &wfa, &wfk, &wfb, &wfe, &wxb, &wv, &wa})
    zero_armed(*v, n);

  // The forward sandwich, rebuilt as saltation_root_analytical_batch builds it.
  // The transpose below reverses that straight-line program assignment by
  // assignment: pruning on the value of s would drop an s^2 factor, which
  // leaves the gradient and stays in the Hessian.
  for (std::size_t i = 0; i < n; ++i) xe[i] = x_before[i] + ws.fb[i] * s;
  sys.first(xe, fe, t_s);
  for (std::size_t i = 0; i < n; ++i)
    xst[i] = x_before[i] + T(0.5) * (ws.fb[i] + fe[i]) * s;
  for (std::size_t i = 0; i < n; ++i) xa[i] = xst[i];
  for (std::size_t j = 0; j < triggered.size(); ++j) {
    const auto& evt = root_events[triggered[j].index];
    if (!evt.terminal)
      cppde::detail::apply_event_action(xa, xa, t_s, evt);
  }
  use_switch_modes(ma);
  sys.first(xa, fa, t_s);
  for (std::size_t i = 0; i < n; ++i) xk[i] = xa[i] - fa[i] * s;
  sys.first(xk, fk, t_ad);

  //  x_out = x_a - (f_a + f_k) s / 2
  T w_s = T(0.0);
  for (std::size_t i = 0; i < n; ++i) {
    wa[i]  = w_out[i];
    wfa[i] = (T(0.0) - T(0.5)) * s * w_out[i];
    wfk[i] = wfa[i];
    w_s   -= T(0.5) * (fa[i] + fk[i]) * w_out[i];
  }

  //  f_k = f(x_k, t)
  adj.jac_t_vec(xk, wfk, t, wv);
  adj.dfdp_t_vec_axpy(xk, wfk, t, T(1.0), w_theta);

  //  x_k = x_a - f_a s
  for (std::size_t i = 0; i < n; ++i) {
    wa[i]  += wv[i];
    wfa[i] -= s * wv[i];
    w_s    -= fa[i] * wv[i];
  }

  //  f_a = f(x_a, t_s), so t_s takes df/dt
  adj.jac_t_vec(xa, wfa, t, wv);
  for (std::size_t i = 0; i < n; ++i) wa[i] += wv[i];
  adj.dfdp_t_vec_axpy(xa, wfa, t, T(1.0), w_theta);
  w_s += adj.dfdt_dot(xa, wfa, t_s);

  //  x_a = R(x_*), every reset reading the surface state the forward reset
  for (std::size_t i = 0; i < n; ++i) w_in[i] = wa[i];
  for (std::size_t j = triggered.size(); j-- > 0;) {
    const std::size_t idx = triggered[j].index;
    const auto& evt = root_events[idx];
    if (evt.terminal) continue;
    const T h = (evt.state_index >= 0 && evt.value_func)
                   ? evt.value_func(xst, t_s) : T(0.0);
    for (std::size_t i = 0; i < n; ++i) wa[i] = w_in[i];
    reset_transpose(
        evt.state_index, evt.method, h, static_cast<int>(idx), xst, t_s, n,
        wa.data(), w_in, w_theta,
        [&](int e, const std::vector<T>& y, const T& tt, std::vector<T>& o)
          { eadj.root_dh_dx(e, y, tt, o); },
        [&](int e, const std::vector<T>& y, const T& tt, const T& sc, T* o)
          { eadj.root_dh_dp_axpy(e, y, tt, sc, o); },
        ws.g,
        [&](int e, const std::vector<T>& y, const T& tt, const T& sc, T* o)
          { eadj.root_dh_dt_axpy(e, y, tt, sc, o); },
        w_s);
  }

  //  x_* = x_b + (f_b + f_e) s / 2
  use_switch_modes(mb);
  for (std::size_t i = 0; i < n; ++i) {
    wxb[i] = w_in[i];
    wfb[i] = T(0.5) * s * w_in[i];
    wfe[i] = wfb[i];
    w_s   += T(0.5) * (ws.fb[i] + fe[i]) * w_in[i];
  }

  //  f_e = f(x_e, t_s), the same pair the other way
  adj.jac_t_vec(xe, wfe, t, wv);
  adj.dfdp_t_vec_axpy(xe, wfe, t, T(1.0), w_theta);
  w_s += adj.dfdt_dot(xe, wfe, t_s);

  //  x_e = x_b + f_b s
  for (std::size_t i = 0; i < n; ++i) {
    wxb[i] += wv[i];
    wfb[i] += s * wv[i];
    w_s    += ws.fb[i] * wv[i];
  }
  for (std::size_t i = 0; i < n; ++i) w_in[i] = wxb[i];

  if (!shifted) return;

  // Where the shift comes from: s = -(grad g . dx + dg/dp . dp) / g_dot, and
  // the forward run corrects it, so ds/ds_lin is 1 + 2 c2 s, one in value and
  // a tangent in the second order.
  const std::size_t idx = triggered[src].index;
  const double c2 = cppde::detail::root_ift_curvature(
      x_before, t_ad, sys, ws.fb, root_events[idx],
      ad_traits::scalar_value(g_dot));
  const T c = (T(0.0) - w_s / g_dot) * (T(1.0) + T(2.0 * c2) * s);
  zero_armed(ws.g, n);
  root_events[idx].dg_dx(x_before, t, ws.g);
  for (std::size_t i = 0; i < n; ++i) w_in[i] += c * ws.g[i];
  eadj.root_dg_dp_axpy(static_cast<int>(idx), x_before, t, c, w_theta);

  // g_dot reads the state too, so the quotient has a second half:
  //
  //   ds/dx = -grad g / g_dot - (s / g_dot) grad g_dot,
  //
  // which contains g, zero in value: absent from the gradient and full
  // strength in the Hessian. The model differentiates g_dot whole, because
  // J' grad g alone holds only for a g linear in x and blind to the clock.
  const T cs = c * s;
  zero_armed(wv, n);
  eadj.root_gdot_dx(static_cast<int>(idx), x_before, t, wv);
  for (std::size_t i = 0; i < n; ++i) w_in[i] += cs * wv[i];
  eadj.root_gdot_dp_axpy(static_cast<int>(idx), x_before, t, cs, w_theta);

  //  f_b = f(x_b, t). grad g_dot already includes the route through f_b that
  //  the shift takes, so only the two Heun legs feed this one.
  adj.jac_t_vec(x_before, wfb, t, wv);
  for (std::size_t i = 0; i < n; ++i) w_in[i] += wv[i];
  adj.dfdp_t_vec_axpy(x_before, wfb, t, T(1.0), w_theta);
}

/// A restart after a stall, transposed. No event fired; the state followed
/// f(t_before, x_before) from t_before to t, past a jump of f in t.
template<class AdjTerms, class T>
void apply_crossing_adjoint(const std::vector<T>& x_before, double t_before,
                            double t, const AdjTerms& adj, std::size_t n,
                            const T* w_after, T* w_before, T* w_theta,
                            jump_workspace<T>& ws)
{
  const double d = t - t_before;
  zero_armed(ws.wy, n);
  for (std::size_t i = 0; i < n; ++i) ws.wy[i] = w_after[i];
  zero_armed(ws.g, n);
  adj.jac_t_vec(x_before, ws.wy, t_before, ws.g);
  for (std::size_t i = 0; i < n; ++i) w_before[i] = w_after[i] + d * ws.g[i];
  adj.dfdp_t_vec_axpy(x_before, ws.wy, t_before, d, w_theta);
}

/// The fixed events at one time, each its own sandwich, applied in order. The
/// last of them applies the root resets a jump switched on.
template<class System, class FixedEvents, class RootEvents, class EvAdj,
         class AdjTerms, class T>
void apply_fixed_jump_adjoint(const std::vector<T>& x_before,
                              double t, const FixedEvents& fixed_events,
                              const RootEvents& root_events,
                              const std::vector<std::size_t>& switched,
                              System& sys, const EvAdj& eadj,
                              const AdjTerms& adj, std::size_t n,
                              const T* w_out, T* w_in, T* w_theta,
                              jump_workspace<T>& ws,
                              const std::vector<signed char>& mb = {},
                              const std::vector<signed char>& ma = {})
{
  using cppde::detail::use_switch_modes;
  // Which of them fire here, and which is last: the same test the engine makes.
  std::vector<int> fired;
  for (std::size_t j = 0; j < fixed_events.size(); ++j)
    if (std::abs(ad_traits::scalar_value(fixed_events[j].time) - t) < 1e-14)
      fired.push_back(static_cast<int>(j));

  // The path the forward took, which the store does not keep: one whole
  // sandwich per event, the switched resets on the surface of the last. A bare
  // reset agrees in value and misses dx_after/dt* = f_before - f_after.
  const std::size_t nf = fired.size();
  std::vector<T> tmp;
  zero_armed(tmp, n);
  // The modes change on the surface of the last of them, as in the forward run.
  auto at_surface = [&](std::vector<T>& xs, const T& te) {
    for (std::size_t q = 0; q < switched.size(); ++q) {
      tmp = xs;
      cppde::detail::apply_event_action(xs, tmp, te, root_events[switched[q]]);
    }
    use_switch_modes(ma);
  };
  ws.path.assign(nf + 1, std::vector<T>());
  ws.path[0] = x_before;
  for (std::size_t j = 0; j < nf; ++j) {
    use_switch_modes(mb);
    zero_armed(ws.path[j + 1], n);
    ws.path[j + 1] = ws.path[j];
    const auto& evt = fixed_events[fired[j]];
    if constexpr (std::is_arithmetic<T>::value) {
      cppde::detail::apply_event_action_fixed(ws.path[j + 1], ws.path[j], evt);
      if (j + 1 == nf) at_surface(ws.path[j + 1], evt.time);
    } else if (j + 1 == nf) {
      cppde::detail::saltation_fixed_analytical(ws.path[j + 1], ws.path[j], sys,
                                                evt, at_surface);
    } else {
      cppde::detail::saltation_fixed_analytical(ws.path[j + 1], ws.path[j], sys,
                                                evt);
    }
  }

  // Backwards through the same chain: the root path's sandwich with the event
  // time's residual tau in place of the shift. Products of two factors that
  // vanish in value are dropped.
  std::vector<T> w(w_out, w_out + n);
  std::vector<T> f1, f2, g1, g2, xe, xk, jv;
  zero_armed(f1, n); zero_armed(f2, n);
  zero_armed(g1, n); zero_armed(g2, n);
  zero_armed(xe, n); zero_armed(xk, n); zero_armed(jv, n);
  ws.wy.assign(n, T(0.0));
  zero_armed(ws.wy, n);

  for (std::size_t j = nf; j-- > 0;) {
    const auto& evt = fixed_events[fired[j]];
    const std::vector<T>& y = ws.path[j];
    const T tau = evt.time - T(ad_traits::scalar_value(evt.time));

    // The half that reaches the surface, replayed so the resets are
    // transposed where the forward applied them.
    use_switch_modes(mb);
    sys.first(y, f1, t);
    for (std::size_t i = 0; i < n; ++i) xe[i] = y[i] + f1[i] * tau;
    sys.first(xe, f2, evt.time);
    const std::size_t n_res = 1 + ((j + 1 == nf) ? switched.size() : 0);
    ws.micro.assign(n_res + 1, std::vector<T>());
    zero_armed(ws.micro[0], n);
    for (std::size_t i = 0; i < n; ++i)
      ws.micro[0][i] = y[i] + T(0.5) * (f1[i] + f2[i]) * tau;
    zero_armed(ws.micro[1], n);
    ws.micro[1] = ws.micro[0];
    cppde::detail::apply_event_action_fixed(ws.micro[1], ws.micro[0], evt);
    for (std::size_t q = 1; q < n_res; ++q) {
      zero_armed(ws.micro[q + 1], n);
      ws.micro[q + 1] = ws.micro[q];
      cppde::detail::apply_event_action(ws.micro[q + 1], ws.micro[q], evt.time,
                                        root_events[switched[q - 1]]);
    }
    const std::vector<T>& xa = ws.micro[n_res];

    // Leaving the surface: x_out = x_a - (g1 + g2) tau / 2, g1 at the event
    // time and g2 at the grid time.
    if (j + 1 == nf) use_switch_modes(ma);
    sys.first(xa, g1, evt.time);
    for (std::size_t i = 0; i < n; ++i) xk[i] = xa[i] - g1[i] * tau;
    sys.first(xk, g2, t);
    for (std::size_t i = 0; i < n; ++i) ws.wy[i] = w[i];
    adj.jac_t_vec(xa, ws.wy, t, jv);
    adj.dfdp_t_vec_axpy(xa, ws.wy, t, T(0.0) - tau, w_theta);
    T w_s = T(0.0);
    for (std::size_t i = 0; i < n; ++i)
      w_s += T(0.5) * (tau * (g1[i] * jv[i]) - (g1[i] + g2[i]) * w[i]);
    // df/dt of g1 at the event time, cotangent -tau w / 2.
    for (std::size_t i = 0; i < n; ++i) ws.wy[i] = (T(0.0) - T(0.5) * tau) * w[i];
    w_s += adj.dfdt_dot(xa, ws.wy, evt.time);
    for (std::size_t i = 0; i < n; ++i) w[i] -= tau * jv[i];

    // The resets on that surface, switched ones first because they came last.
    for (std::size_t q = n_res; q-- > 1;) {
      const auto& re = root_events[switched[q - 1]];
      const std::vector<T>& yy = ws.micro[q];
      const T h = (re.state_index >= 0 && re.value_func)
                     ? re.value_func(yy, evt.time) : T(0.0);
      ws.wy = w;
      reset_transpose(
          re.state_index, re.method, h, static_cast<int>(switched[q - 1]), yy,
          evt.time, n, w.data(), ws.wy.data(), w_theta,
          [&](int e, const std::vector<T>& yv, const T& tt, std::vector<T>& o)
            { eadj.root_dh_dx(e, yv, tt, o); },
          [&](int e, const std::vector<T>& yv, const T& tt, const T& sc, T* o)
            { eadj.root_dh_dp_axpy(e, yv, tt, sc, o); },
          ws.g,
          [&](int e, const std::vector<T>& yv, const T& tt, const T& sc, T* o)
            { eadj.root_dh_dt_axpy(e, yv, tt, sc, o); },
          w_s);
      w.swap(ws.wy);
    }
    {
      const std::vector<T>& yy = ws.micro[0];
      const T h = (evt.state_index >= 0 && evt.value_func)
                     ? evt.value_func(yy, evt.time) : T(0.0);
      ws.wy = w;
      reset_transpose(
          evt.state_index, evt.method, h, fired[j], yy, evt.time, n,
          w.data(), ws.wy.data(), w_theta,
          [&](int e, const std::vector<T>& yv, const T& tt, std::vector<T>& o)
            { eadj.fixed_dh_dx(e, yv, tt, o); },
          [&](int e, const std::vector<T>& yv, const T& tt, const T& sc, T* o)
            { eadj.fixed_dh_dp_axpy(e, yv, tt, sc, o); },
          ws.g,
          [&](int e, const std::vector<T>& yv, const T& tt, const T& sc, T* o)
            { eadj.fixed_dh_dt_axpy(e, yv, tt, sc, o); },
          w_s);
      w.swap(ws.wy);
    }

    // Entering it: x_* = x_b + (f1 + f2) tau / 2, the same pair the other way.
    use_switch_modes(mb);
    for (std::size_t i = 0; i < n; ++i) ws.wy[i] = w[i];
    adj.jac_t_vec(y, ws.wy, t, jv);
    adj.dfdp_t_vec_axpy(xe, ws.wy, t, tau, w_theta);
    for (std::size_t i = 0; i < n; ++i)
      w_s += T(0.5) * ((f1[i] + f2[i]) * w[i] + tau * (f1[i] * jv[i]));
    // df/dt of f2 at the event time, cotangent tau w / 2.
    for (std::size_t i = 0; i < n; ++i) ws.wy[i] = T(0.5) * tau * w[i];
    w_s += adj.dfdt_dot(xe, ws.wy, evt.time);
    for (std::size_t i = 0; i < n; ++i) w[i] += tau * jv[i];

    // kappa scales the event time's derivative, which reads the parameters.
    eadj.fixed_dtime_dp_axpy(fired[j], w_s, w_theta);
  }

  for (std::size_t i = 0; i < n; ++i) w_in[i] = w[i];
}

/// One seed inside a step: the time and the row of seeds.
template<class T>
struct obs_ref { double t; const T* w; };

/// The most substeps a checked step is swept in, the most pieces for the
/// extrapolated flow.
constexpr int default_max_sub = 1024;
/// From this many substeps on, a halving that does not shrink the difference
/// ends the refinement.
constexpr int stall_substeps = 32;

/// The cotangent on the state entering a jump that was read off the dense
/// output of the step below it, waiting for that step's sweep.
template<class T>
struct jump_handover {
  std::vector<T> w;
  double t = 0.0;
  bool armed = false;

  /// Hands `w_before` on to the step below, or onto the initial state `w_x`
  /// when the jump is at the run's own start.
  template<class Event>
  void take(const Event& e, const std::vector<T>& w_before, std::vector<T>& w_x) {
    if (e.after_step > 0) {
      w = w_before; t = e.t_before; armed = true;
    } else {
      for (std::size_t i = 0; i < w_x.size(); ++i) w_x[i] += w_before[i];
    }
  }
};

/// The seeds step `k` holds into `out`: the handover a jump left, then the
/// observations inside it. One a jump produced is a value and is seeded at the
/// boundary instead. Returns whether anything observes inside the step.
template<class Store, class T>
bool step_seeds(const Store& store, std::size_t k, const T* seeds, std::size_t n,
                std::size_t& next_obs, jump_handover<T>& h, std::vector<obs_ref<T>>& out)
{
  std::size_t obs_lo = next_obs;
  while (obs_lo > 0 && store.obs(obs_lo - 1).step == k + 1) --obs_lo;
  const bool observed = obs_lo < next_obs || h.armed;
  out.clear();
  if (h.armed) { out.push_back({h.t, h.w.data()}); h.armed = false; }
  for (std::size_t o = obs_lo; o < next_obs; ++o) {
    if (store.obs(o).event < store.n_events()) continue;
    out.push_back({store.obs(o).t, seeds + o * n});
  }
  next_obs = obs_lo;
  return observed;
}

/// The observations before the first step, which are the initial state itself,
/// except the value a jump at the start produced, which that jump seeded.
template<class Store, class T>
void seed_initial_state(const Store& store, std::size_t next_obs, const T* seeds,
                        std::size_t n, std::vector<T>& w_x)
{
  while (next_obs > 0) {
    --next_obs;
    if (store.obs(next_obs).event < store.n_events()) continue;
    const T* w = seeds + next_obs * n;
    for (std::size_t i = 0; i < n; ++i) w_x[i] += w[i];
  }
}

/// The adjoint of the intervention `ei`: `w_after`, the cotangent on the state
/// it produced, gets the observations of that value and goes through the jump
/// into `w_before` on the state entering it, the share of theta into `wp`.
template<class Store, class AdjTerms, class Jumps, class T>
void event_adjoint(const Store& store, std::size_t ei, const T* seeds, std::size_t n,
                   const AdjTerms& adj, const Jumps& jumps, std::vector<T>& w_after,
                   std::vector<T>& w_before, T* wp, jump_workspace<T>& ws)
{
  const auto& e = store.event(ei);
  for (std::size_t o = 0; o < store.n_obs(); ++o)
    if (store.obs(o).event == ei)
      for (std::size_t i = 0; i < n; ++i) w_after[i] += seeds[o * n + i];
  w_before.assign(n, T(0.0));
  if (e.crossing)
    apply_crossing_adjoint(e.x_before, e.t_before, e.t, adj, n,
                           w_after.data(), w_before.data(), wp, ws);
  else if (e.root)
    apply_root_jump_adjoint(e.x_before, e.x_after, e.t, jumps.root,
                            e.triggered, jumps.sys, jumps.eadj, adj, n,
                            w_after.data(), w_before.data(), wp, ws,
                            e.modes_before, e.modes_after);
  else
    apply_fixed_jump_adjoint(e.x_before, e.t, jumps.fixed, jumps.root,
                             e.switched, jumps.sys, jumps.eadj, adj, n,
                             w_after.data(), w_before.data(), wp, ws,
                             e.modes_before, e.modes_after);
}

// ---------------------------------------------------------------------------
//  A run that took no step, which is a run over a single time. Every
//  observation there is the initial state itself, or the value a jump at that
//  time produced, which reaches the initial state and the parameters through
//  that jump's adjoint. The flow in between has length zero and is the
//  identity.
// ---------------------------------------------------------------------------
template<class Store, class AdjTerms, class Jumps, class T>
void sweep_without_steps(const Store& store, std::size_t n, const T* seeds,
                         const AdjTerms& adj, const Jumps& jumps,
                         std::vector<T>& w_x, std::vector<T>& w_theta,
                         jump_workspace<T>& ws)
{
  seed_initial_state(store, store.n_obs(), seeds, n, w_x);
  if constexpr (Jumps::active) {
    std::vector<T> w_after, w_before;
    for (std::size_t ei = 0; ei < store.n_events(); ++ei) {
      w_after.assign(n, T(0.0));
      event_adjoint(store, ei, seeds, n, adj, jumps, w_after, w_before,
                    w_theta.data(), ws);
      for (std::size_t i = 0; i < n; ++i) w_x[i] += w_before[i];
    }
  } else {
    (void)adj; (void)jumps; (void)w_theta; (void)ws;
  }
}

// ---------------------------------------------------------------------------
//  A whole multistep trajectory backwards, step adjoints and jump adjoints in
//  reverse order. An observation inside a step takes its row
//  d x_interp / d (zn_pred, acor) from the probe's interpolant, so an observed
//  step rebuilds the operators rather than taking the cached ones.
// ---------------------------------------------------------------------------
template<class Stepper>
class closed_multistep_trajectory {
public:
  // The sweep runs in the stepper's own scalar type. In plain double that is
  // the gradient; over a dual it is the gradient and its directional
  // derivatives, which is forward over reverse.
  using scalar_type = typename Stepper::value_type;

  /// Whether the sweep also keeps lambda and the refinement indicator per step.
  void trace_lambda(bool on) { m_trace = on; }

  /// One seed column. `seeds` is [n_obs x n_states] row-major.
  template<class Store, class AdjTerms, class Solver, class Jumps = no_jumps>
  void sweep(const Store& store, std::size_t n_phi, const scalar_type* seeds,
             const AdjTerms& adj, Solver& solver, const Jumps& jumps = Jumps())
  {
    using T = scalar_type;
    const std::size_t n = store.n_states();
    const std::size_t n_steps = store.n_steps();

    zero_armed(m_wp, n_phi);
    m_wx.assign(n, T(0.0));
    m_whist.clear();
    m_lam.assign(m_trace ? n_steps * n : 0u, T(0.0));
    m_eta.assign(m_trace ? n_steps : 0u, T(0.0));
    m_ref_m.assign(m_trace ? n_steps : 0u, 1);
    if (n_steps == 0) {
      sweep_without_steps(store, n, seeds, adj, jumps, m_wx, m_wp, m_jws);
      return;
    }

    // The cotangent the step above hands down, on its own carry. Across an
    // intervention it is instead a cotangent on a point inside the step below,
    // because the restart threw the carry away.
    std::vector<T> w_carry, w_after, w_before;
    jump_handover<T> handover;
    std::size_t next_obs = store.n_obs();

    std::vector<T> w_in, jtv;
    zero_armed(jtv, n);

    for (std::size_t k = n_steps; k-- > 0;) {
      const auto& cp = store.step(k);
      store.use_step_modes(k);

      // Does anything observe inside this step? Then the probe has to hold
      // this step's own interpolant, not the one it kept from another.
      const bool observed = step_seeds(store, k, seeds, n, next_obs, handover, m_obs);
      if (observed) m_probe.valid = false;

      T eta_k = T(0.0);
      w_in.clear();
      ms_step_adjoint(cp, observed, m_obs.data(), m_obs.size(), w_carry, n, n_phi,
                      adj, solver, w_in, m_wp.data(), m_trace ? &eta_k : nullptr);

      // lambda is what this step hands the one below it, on the state slot.
      if (m_trace) {
        for (std::size_t i = 0; i < n; ++i) m_lam[k * n + i] = w_in[i];
        m_eta[k] = eta_k;
      }
      w_carry.swap(w_in);

      // The intervention this step was entered through, if any. Two maps in the
      // order the forward run applied them, so swept the other way round: the
      // restart, which collapses the carry onto the state the jump ended on,
      // and the jump itself.
      if constexpr (Jumps::active) {
        const std::size_t ei = store.event_before(k);
        if (ei >= store.n_events()) continue;
        const auto& e = store.event(ei);
        w_after.assign(n, T(0.0));
        m_ws.x.assign(e.x_after.begin(), e.x_after.end());
        collapse_restart(w_carry, n, m_ws.x, e.t,
                         e.restart ? e.dt_restart
                                   : static_cast<double>(cp.carry.h),
                         adj, w_after.data(), m_wp.data(), m_ws.mu, jtv);
        event_adjoint(store, ei, seeds, n, adj, jumps, w_after, w_before,
                      m_wp.data(), m_jws);
        w_carry.clear();
        handover.take(e, w_before, m_wx);
      }
    }

    // The trajectory start, which is the same restart with no jump under it.
    if (!w_carry.empty()) {
      const auto& cp0 = store.step(0);
      store.use_step_modes(0);
      w_after.assign(n, T(0.0));
      m_ws.x.assign(cp0.start_state(), cp0.start_state() + n);
      collapse_restart(w_carry, n, m_ws.x, cp0.t,
                       static_cast<double>(cp0.carry.h), adj, w_after.data(),
                       m_wp.data(), m_ws.mu, jtv);
      for (std::size_t i = 0; i < n; ++i) m_wx[i] += w_after[i];
    }
    seed_initial_state(store, next_obs, seeds, n, m_wx);
  }

  /// The sweep that takes the flow between the forward run's grid points
  /// instead of the multistep scheme: each step's interval from the state the
  /// step started at, swept by `flow`. Only the state's cotangent passes from
  /// one step to the next, as the exact flow depends on nothing else.
  template<class Store, class System, class AdjTerms, class Flow, class Jumps = no_jumps>
  void sweep_flow(const Store& store, std::size_t n_phi, const scalar_type* seeds,
                  System& sys, const AdjTerms& adj, Flow& flow,
                  const Jumps& jumps = Jumps())
  {
    using T = scalar_type;
    const std::size_t n = store.n_states();
    const std::size_t n_steps = store.n_steps();
    zero_armed(m_wp, n_phi);
    m_wx.assign(n, T(0.0));
    m_whist.clear();
    m_lam.assign(m_trace ? n_steps * n : 0u, T(0.0));
    // The flow has no error estimate of the multistep scheme to weigh.
    m_eta.assign(m_trace ? n_steps : 0u, T(std::numeric_limits<double>::quiet_NaN()));
    m_ref_m.assign(m_trace ? n_steps : 0u, 1);
    flow.begin_flow(n_phi);
    if (n_steps == 0) {
      sweep_without_steps(store, n, seeds, adj, jumps, m_wx, m_wp, m_jws);
      return;
    }
    std::size_t next_obs = store.n_obs();
    std::vector<T> w_out(n, T(0.0)), w_in(n, T(0.0)), x0, w_after, w_before;
    jump_handover<T> handover;
    for (std::size_t k = n_steps; k-- > 0;) {
      const auto& cp = store.step(k);
      store.use_step_modes(k);
      step_seeds(store, k, seeds, n, next_obs, handover, m_obs);
      x0.assign(cp.start_state(), cp.start_state() + n);
      const int mk = flow.flow_interval(sys, adj, n, n_phi, x0, cp.t, cp.dt,
                                        m_obs, w_out, w_in);
      if (m_trace) {
        for (std::size_t i = 0; i < n; ++i) m_lam[k * n + i] = w_in[i];
        m_ref_m[k] = mk;
      }
      w_out.swap(w_in);

      // The intervention this step was entered through, on the state alone.
      if constexpr (Jumps::active) {
        const std::size_t ei = store.event_before(k);
        if (ei < store.n_events()) {
          w_after = w_out;
          event_adjoint(store, ei, seeds, n, adj, jumps, w_after, w_before,
                        m_wp.data(), m_jws);
          w_out.assign(n, T(0.0));
          handover.take(store.event(ei), w_before, m_wx);
        }
      }
    }
    for (std::size_t i = 0; i < n; ++i) m_wx[i] += w_out[i];
    seed_initial_state(store, next_obs, seeds, n, m_wx);
    const auto& fp = flow.wp();
    for (std::size_t i = 0; i < n_phi; ++i) m_wp[i] += fp[i];
  }

  /// Under trace_lambda: the substeps each step's interval was swept in.
  const std::vector<int>& refine_m() const { return m_ref_m; }

private:
  using T_ = scalar_type;

  /// The adjoint of one step from its checkpoint: `w_carry` the cotangent on
  /// the carry the step hands on, empty where an intervention handed a point
  /// inside the step instead; `w_in` receives the one on the carry it was
  /// entered with, `wp` the step's share of the gradient, `eta` the cotangent
  /// on acor against the step's error estimate.
  template<class CP, class AdjTerms, class Solver>
  void ms_step_adjoint(const CP& cp, bool observed, const obs_ref<T_>* ob, std::size_t nob,
                       std::vector<T_>& w_carry, std::size_t n, std::size_t n_phi,
                       const AdjTerms& adj, Solver& solver, std::vector<T_>& w_in,
                       T_* wp, T_* eta)
  {
    using T = T_;
    multistep_operators<Stepper>& ops = m_ops;
    // An observed step needs its own interpolant, not one kept from another.
    if (observed || nob > 0) m_probe.valid = false;
    const double tail_key =
        cp.q_next + 1e3 * cp.eta + 1e6 * static_cast<double>(cp.ops.size());
    // Timed on the rebuild alone, so the call count is the miss count.
    if (!m_probe.matches(cp.carry, cp.dt, tail_key,
                         cp.ops_recorded ? &cp.ops : nullptr)) {
      auto _tp = m_prof.timer(cppde::prof_cat::rev_operators);
      m_probe.rebuild(cp.carry, cp.dt,
                      [&](typename multistep_probe<Stepper>::probe_stepper& pr)
                        { cp.apply_tail(pr, m_null); }, ops);
    }

    const std::size_t nz_in = static_cast<std::size_t>(ops.q_in + 1) * n;
    const std::size_t nz_out = static_cast<std::size_t>(ops.q_out + 1) * n;
    m_ws.w_pred.assign(nz_in, T(0.0));
    m_ws.w_acor.assign(n, T(0.0));

    // A cotangent at a time inside the step, through its own interpolant. The
    // probe's states are the slots plus one for acor, so the row it writes is
    // d x_interp / d (zn_pred[0..q], acor).
    for (std::size_t o = 0; o < nob; ++o) {
      auto _tp = m_prof.timer(cppde::prof_cat::rev_interp);
      const T* w = ob[o].w;
      m_dense_row.assign(static_cast<std::size_t>(ops.q_in + 2), 0.0);
      m_probe.tail_probe.eval_dense_into(clamp_to_step(cp, ob[o].t), m_dense_row);
      for (int j = 0; j <= ops.q_in; ++j)
        for (std::size_t i = 0; i < n; ++i)
          m_ws.w_pred[static_cast<std::size_t>(j) * n + i] +=
              m_dense_row[static_cast<std::size_t>(j)] * w[i];
      const std::size_t acor_slot = m_dense_row.size() - 1;
      for (std::size_t i = 0; i < n; ++i)
        m_ws.w_acor[i] += m_dense_row[acor_slot] * w[i];
    }

    // What the step above handed back, on its carry.
    if (!w_carry.empty()) {
      auto _tp = m_prof.timer(cppde::prof_cat::rev_adjoint);
      w_carry.resize(nz_out, T(0.0));
      carry_into_pred(ops, n, w_carry.data(), m_ws.w_pred.data(),
                      m_ws.w_acor.data());
    }

    const double t_new = cp.t + ops.h;
    solver.prepare(cp.y, t_new, 1.0 / ops.gamma, ops.gamma);

    w_in.assign(nz_in, T(0.0));
    if (eta) m_wout.assign(n, T(0.0));
    // Timed either side of the solve, which reports itself.
    { auto _tp = m_prof.timer(cppde::prof_cat::rev_adjoint);
      multistep_adjoint_rhs(ops, n, m_ws, eta ? m_wout.data() : nullptr); }
    solver.transposed(m_ws.mu);
    { auto _tp = m_prof.timer(cppde::prof_cat::rev_adjoint);
      multistep_adjoint_finish(ops, n, n_phi, cp.y, t_new, adj, w_in.data(),
                               wp, m_ws); }

    // eta is the cotangent on acor against the step's own error estimate,
    // which for a corrector method is acor scaled by the order's constant.
    if (eta) {
      T e = T(0.0);
      for (std::size_t i = 0; i < n; ++i) {
        T pred0 = T(0.0);
        for (int j = 0; j <= ops.q_in; ++j)
          pred0 += ops.A(0, j) * cp.zn[static_cast<std::size_t>(j) * n + i];
        e += m_wout[i] * (cp.y[i] - pred0);
      }
      *eta = e * ops.err_scale;
    }
  }

public:
  /// Per-category timings of the sweep, to stderr. Compiled away without
  /// CPPDE_PROFILE. The transposed algebra reports itself, from the solver.
  void report_profile() const { m_prof.report("cppDE written adjoint"); }

  const std::vector<scalar_type>& wx0() const { return m_wx; }
  const std::vector<scalar_type>& whistory0() const { return m_whist; }
  const std::vector<scalar_type>& wp() const { return m_wp; }

  /// Under trace_lambda: [n_steps, n_states] step-major, and one per step.
  const std::vector<scalar_type>& lambda() const { return m_lam; }
  const std::vector<scalar_type>& eta() const { return m_eta; }

private:
  // The tail reads the right-hand side only for an order-one restart, which
  // belongs to a boundary and not to a step.
  struct null_rhs {
    void operator()(const std::vector<double>&, std::vector<double>& d,
                    const double&) const { d.assign(d.size(), 0.0); }
  };
  struct null_sys { null_rhs first; };

  null_sys m_null;
  cppde::profiler m_prof;
  jump_workspace<scalar_type> m_jws;
  multistep_probe<Stepper> m_probe;
  multistep_workspace<scalar_type> m_ws;
  // The step adjoint's buffers.
  multistep_operators<Stepper> m_ops;
  std::vector<double> m_dense_row;
  std::vector<obs_ref<scalar_type>> m_obs;
  std::vector<int> m_ref_m;
  bool m_trace = false;
  std::vector<scalar_type> m_wx, m_whist, m_wp, m_lam, m_eta, m_wout;
};

// ---------------------------------------------------------------------------
//  The adjoint of one explicit Runge-Kutta step, stages recomputed by one
//  forward run: m_i = h b_i w + h sum_{j>i} a_ji u_j with u_j = J(X_j)' m_j,
//  w_x = w + sum_i u_i. See vignette("Methods"), "The step map and its adjoint".
// ---------------------------------------------------------------------------
template<class T = double>
struct onestep_workspace {
  std::vector<T> xout, xerr, m, u, x_stage, q;
  std::vector<std::vector<T> > mm;
};

/// The backward recursion alone, on a stepper that has just run the step. A
/// trajectory runs the step itself, because its interpolant reads the stages
/// before the recursion consumes them.
template<class Stepper, class AdjTerms, class T>
void onestep_recurse(Stepper& st,
                     const T* x0, double t, double dt,
                     std::size_t n, std::size_t n_phi,
                     const T* w_out,
                     const AdjTerms& adj,
                     T* w_in, T* w_theta,
                     onestep_workspace<T>& ws,
                     T* w_out_state = nullptr)
{
  constexpr int S = Stepper::n_stages_used;

  if (w_out_state) for (std::size_t i = 0; i < n; ++i) w_out_state[i] = w_out[i];

  if (ws.mm.size() < static_cast<std::size_t>(S) + 1) ws.mm.resize(S + 1);
  std::vector<std::vector<T> >& U = ws.mm;
  for (int i = 1; i <= S; ++i) U[i].assign(n, T(0.0));
  zero_armed(ws.u, n);

  for (std::size_t i = 0; i < n; ++i) w_in[i] = w_out[i];

  // Newest stage first: u_j is needed by every earlier stage and by nothing
  // later, so one pass suffices.
  for (int i = S; i >= 1; --i) {
    ws.m.assign(n, T(0.0));
    const double bi = Stepper::tableau_b(i);
    for (std::size_t k = 0; k < n; ++k) ws.m[k] = dt * bi * w_out[k];
    for (int j = i + 1; j <= S; ++j) {
      const double aji = Stepper::tableau_a(j, i);
      if (aji == 0.0) continue;
      for (std::size_t k = 0; k < n; ++k) ws.m[k] += dt * aji * U[j][k];
    }

    // The stage state, rebuilt from the step start and the stage derivatives.
    ws.x_stage.assign(x0, x0 + n);
    for (int j = 1; j < i; ++j) {
      const double aij = Stepper::tableau_a(i, j);
      if (aij == 0.0) continue;
      const auto& kj = st.stage_k(j);
      for (std::size_t k = 0; k < n; ++k) ws.x_stage[k] += dt * aij * kj[k];
    }
    const double ti = t + Stepper::tableau_c(i) * dt;

    adj.jac_t_vec(ws.x_stage, ws.m, ti, ws.u);
    U[i] = ws.u;
    for (std::size_t k = 0; k < n; ++k) w_in[k] += ws.u[k];

    adj.dfdp_t_vec_axpy(ws.x_stage, ws.m, ti, 1.0, w_theta);
  }
}

/// Step and recursion in one, for a caller with one step to adjoint.
template<class Stepper, class System, class AdjTerms, class T>
void apply_onestep_adjoint(Stepper& st, System& sys,
                           const T* x0, double t, double dt,
                           std::size_t n, std::size_t n_phi,
                           const T* w_out,
                           const AdjTerms& adj,
                           T* w_in, T* w_theta,
                           onestep_workspace<T>& ws,
                           T* w_out_state = nullptr)
{
  zero_armed(ws.xout, n);
  zero_armed(ws.xerr, n);
  std::vector<T> x(x0, x0 + n);
  st.do_step(sys, x, t, ws.xout, dt, ws.xerr);
  onestep_recurse(st, x0, t, dt, n, n_phi, w_out, adj, w_in, w_theta, ws,
                  w_out_state);
}

// ---------------------------------------------------------------------------
//  The adjoint of one Rosenbrock step. dg = W^-1 (dr + (dJ) g), so besides the
//  transposed solves against the step's own W it needs lambda contracted with
//  d(J g). See vignette("Methods"), section "The step map and its adjoint".
// ---------------------------------------------------------------------------
template<class T = double>
struct rosenbrock_workspace {
  std::vector<T> xout, xerr, x, lam, wx, wD, jv, xs, M;
  std::vector<std::vector<T> > mu, wX;
};

// Whether the generated terms take sum_s lam_s' J v_s over the Jacobian
// pattern in one pullback, rather than one contraction per stage.
template<class A, class = void> struct has_jvp_bundle : std::false_type {};
template<class A>
struct has_jvp_bundle<A, std::void_t<decltype(A::jvp_bundle_rows())>>
: std::integral_constant<bool, A::has_jvp_bundle> {};

// Whether the generated terms take the right-hand side's pullback onto the
// state and theta in one pass.
template<class A, class = void> struct has_vjp_fused : std::false_type {};
template<class A>
struct has_vjp_fused<A, std::void_t<decltype(&A::vjp_t_vec)>> : std::true_type {};

// out_x = J(x)' lam, out_p[n + k] += sc (f_p(x)' lam)_k.
template<class AdjTerms, class V, class Tm, class T>
inline void rhs_pullback(const AdjTerms& adj, const V& x, const V& lam, const Tm& t,
                         double sc, std::vector<T>& out_x, T* out_p)
{
  if constexpr (has_vjp_fused<AdjTerms>::value) {
    adj.vjp_t_vec(x, lam, t, sc, out_x, out_p);
  } else {
    adj.jac_t_vec(x, lam, t, out_x);
    adj.dfdp_t_vec_axpy(x, lam, t, sc, out_p);
  }
}

// The bilinear form M over the Jacobian's pattern pulled back onto the state,
// out_x, and added onto theta, out_p, in one pass where the terms fuse them.
template<class A, class = void> struct has_jvp_mat_fused : std::false_type {};
template<class A>
struct has_jvp_mat_fused<A, std::void_t<decltype(&A::jvp_t_mat)>> : std::true_type {};

template<class AdjTerms, class V, class Tm, class T>
inline void bundle_pullback(const AdjTerms& adj, const V& x, const V& M, const Tm& t,
                            std::vector<T>& out_x, T* out_p)
{
  if constexpr (has_jvp_mat_fused<AdjTerms>::value) {
    adj.jvp_t_mat(x, M, t, 1.0, out_x, out_p);
  } else {
    adj.jvp_x_t_mat(x, M, t, out_x);
    adj.jvp_p_t_mat_axpy(x, M, t, 1.0, out_p);
  }
}

template<class Stepper, class System, class AdjTerms, class T>
void apply_rosenbrock_adjoint(Stepper& st, System& sys,
                              const T* x0, double t, double dt,
                              std::size_t n, std::size_t n_phi,
                              const T* w_out,
                              const AdjTerms& adj,
                              T* w_in, T* w_theta,
                              rosenbrock_workspace<T>& ws,
                              const T* mu_seed = nullptr,
                              T* w_out_state = nullptr,
                              T* w_in_emb = nullptr, T* w_theta_emb = nullptr)
{
  constexpr int S = Stepper::n_stages_used;   // five stages, then the error

  // The step once forward, to recover the stages and the factorisation. An
  // accepted step rebuilds the Jacobian at its own start, so this is the same
  // matrix and the same stages the run produced.
  zero_armed(ws.xout, n);
  zero_armed(ws.xerr, n);
  zero_armed(ws.jv, n);
  ws.x.assign(x0, x0 + n);
  st.do_step(sys, ws.x, t, ws.xout, dt, ws.xerr);

  if (w_out_state) for (std::size_t i = 0; i < n; ++i) w_out_state[i] = w_out[i];

  if (ws.mu.size() < static_cast<std::size_t>(S) + 2) ws.mu.resize(S + 2);
  if (ws.wX.size() < static_cast<std::size_t>(S) + 2) ws.wX.resize(S + 2);
  for (int i = 1; i <= S + 1; ++i) {
    ws.mu[i].assign(n, T(0.0));
    ws.wX[i].assign(n, T(0.0));
  }
  ws.wx.assign(n, T(0.0));
  ws.wD.assign(n, T(0.0));
  bool wD_touched = false;

  // An observation inside the step reaches the stages directly: the continuous
  // extension is linear in them.
  if (mu_seed)
    for (int i = 1; i <= S; ++i)
      for (std::size_t k = 0; k < n; ++k)
        ws.mu[i][k] += mu_seed[static_cast<std::size_t>(i - 1) * n + k];

  // The stage state X_i, rebuilt from the step start and the stage vectors.
  // The error stage reads the last stage's state plus that stage's own vector,
  // so it takes the last stage's combination and not a row of its own.
  auto stage_state = [&](int i) {
    const int row = (i == S + 1) ? S : i;
    ws.xs.assign(x0, x0 + n);
    for (int j = 1; j < row; ++j) {
      const double a = Stepper::stage_a(row, j);
      if (a == 0.0) continue;
      const auto& gj = st.stage_g(j);
      for (std::size_t k = 0; k < n; ++k) ws.xs[k] += a * gj[k];
    }
    if (i == S + 1) {
      const auto& gs = st.stage_g(S);
      for (std::size_t k = 0; k < n; ++k) ws.xs[k] += gs[k];
    }
  };

  // One solved stage: the matrix term, the right-hand side's couplings, and
  // the stage's own evaluation of f.
  auto sweep_stage_into = [&](int i, const std::vector<T>& gi, T* wXi, T* w_theta) {
    ws.lam = ws.mu[i];
    st.stage_solve_transposed(ws.lam);

    if constexpr (has_jvp_bundle<AdjTerms>::value) {
      const int* r = AdjTerms::jvp_bundle_rows();
      const int* c = AdjTerms::jvp_bundle_cols();
      for (int k = 0; k < AdjTerms::jvp_bundle_size; ++k) ws.M[k] += ws.lam[r[k]] * gi[c[k]];
    } else {
      adj.jvp_x_t_vec(ws.x, gi, ws.lam, t, ws.jv);
      for (std::size_t k = 0; k < n; ++k) ws.wx[k] += ws.jv[k];
      adj.jvp_p_t_vec_axpy(ws.x, gi, ws.lam, t, 1.0, w_theta);
    }

    for (int j = 1; j < i; ++j) {
      const double c = Stepper::stage_c(i, j) / dt;
      if (c == 0.0) continue;
      for (std::size_t k = 0; k < n; ++k) ws.mu[j][k] += c * ws.lam[k];
    }
    const double d = Stepper::stage_d(i);
    if (d != 0.0) {
      for (std::size_t k = 0; k < n; ++k) ws.wD[k] += dt * d * ws.lam[k];
      wD_touched = true;
    }

    const double ti = t + Stepper::stage_node(i) * dt;
    if (i == 1) {
      rhs_pullback(adj, ws.x, ws.lam, ti, 1.0, ws.jv, w_theta);
      for (std::size_t k = 0; k < n; ++k) ws.wx[k] += ws.jv[k];
      return;
    }
    stage_state(i);
    rhs_pullback(adj, ws.xs, ws.lam, ti, 1.0, ws.jv, w_theta);
    for (std::size_t k = 0; k < n; ++k) wXi[k] += ws.jv[k];
  };

  // The cotangent on X_i, spread onto the step start and the stages it reads.
  auto spread_stage_state = [&](int i, const T* wXi) {
    for (std::size_t k = 0; k < n; ++k) ws.wx[k] += wXi[k];
    for (int j = 1; j < i; ++j) {
      const double a = Stepper::stage_a(i, j);
      if (a == 0.0) continue;
      for (std::size_t k = 0; k < n; ++k) ws.mu[j][k] += a * wXi[k];
    }
  };

  // The stages backwards from the seeds on X_6 and on e, onto the step start
  // and theta.
  auto backward = [&](bool on_state, T* wth) {
    if constexpr (has_jvp_bundle<AdjTerms>::value)
      ws.M.assign(AdjTerms::jvp_bundle_size, T(0.0));
    for (std::size_t k = 0; k < n; ++k) {
      if (on_state) ws.wX[S + 1][k] += w_out[k];
      ws.mu[S + 1][k] += w_out[k];
    }
    // The error solve, whose stage vector is the error estimate itself.
    sweep_stage_into(S + 1, ws.xerr, ws.wX[S + 1].data(), wth);
    // X_6 = X_5 + g_5
    for (std::size_t k = 0; k < n; ++k) {
      ws.wX[S][k] += ws.wX[S + 1][k];
      ws.mu[S][k] += ws.wX[S + 1][k];
    }
    for (int i = S; i >= 1; --i) {
      sweep_stage_into(i, st.stage_g(i), ws.wX[i].data(), wth);
      if (i >= 2) spread_stage_state(i, ws.wX[i].data());
    }
    if constexpr (has_jvp_bundle<AdjTerms>::value) {
      bundle_pullback(adj, ws.x, ws.M, t, ws.jv, wth);
      for (std::size_t k = 0; k < n; ++k) ws.wx[k] += ws.jv[k];
    }
    // D = df/dt(x, t), which the Jacobian evaluation filled and every early
    // stage reads.
    if (wD_touched) {
      adj.dfdt_x_t_vec(ws.x, ws.wD, t, ws.jv);
      for (std::size_t k = 0; k < n; ++k) ws.wx[k] += ws.jv[k];
      adj.dfdt_p_t_vec_axpy(ws.x, ws.wD, t, 1.0, wth);
    }
  };
  // x_out = X_6 + e
  backward(true, w_theta);
  for (std::size_t k = 0; k < n; ++k) w_in[k] = ws.wx[k];

  // The embedded solution is X_6 alone, so the two adjoints differ by the
  // response to the seed on e: the local error of this step's adjoint, at
  // the embedded order, on the same stages and factorisation.
  if (w_in_emb) {
    for (int i = 1; i <= S + 1; ++i) {
      ws.mu[i].assign(n, T(0.0));
      ws.wX[i].assign(n, T(0.0));
    }
    ws.wx.assign(n, T(0.0));
    ws.wD.assign(n, T(0.0));
    wD_touched = false;
    backward(false, w_theta_emb);
    for (std::size_t k = 0; k < n; ++k) w_in_emb[k] = ws.wx[k];
  }
  (void)n_phi;
}

// Whether a one-step method's stages are linear solves it can hand out, rather
// than explicit combinations its tableau already describes. That is what decides
// how its adjoint moves through a step.
template<class S, class = void> struct has_stage_vectors : std::false_type {};
template<class S>
struct has_stage_vectors<S, std::void_t<decltype(std::declval<const S&>().stage_g(1))>>
: std::true_type {};

// ---------------------------------------------------------------------------
//  A whole trajectory backwards on a one-step method, which passes only the
//  state across a step boundary. An observation inside a step goes through the
//  continuous extension, whose weights the stepper hands out.
// ---------------------------------------------------------------------------
template<class Stepper>
class closed_onestep_trajectory {
public:
  // As in the multistep form: the sweep runs in the stepper's scalar type.
  using scalar_type = typename Stepper::value_type;

  /// Whether the sweep also keeps lambda and the refinement indicator per step.
  void trace_lambda(bool on) { m_trace = on; }

  /// Hold each step's adjoint to the local error test of CVODES' backward
  /// problem: lambda to rtol, atol, the share of the gradient to rtol, gradtol
  /// when gradtol > 0; rtol 0 takes the steps as they are.
  void set_refine(double rtol, double atol, double gradtol = 0.0,
                  int max_sub = default_max_sub) {
    m_refine_rtol = rtol; m_refine_atol = atol; m_refine_gradtol = gradtol;
    m_refine_max = max_sub;
  }
  /// The steps of the last sweep that met the test at no number of substeps.
  std::size_t refine_failures() const { return m_ref_fail; }
  /// Under trace_lambda: the substeps each step was swept in.
  const std::vector<int>& refine_m() const { return m_ref_m; }

  /// One seed column. `seeds` is [n_obs x n_states] row-major.
  template<class Store, class System, class AdjTerms, class Jumps = no_jumps>
  void sweep(const Store& store, std::size_t n_phi, const scalar_type* seeds,
             System& sys, const AdjTerms& adj, Stepper& st,
             const Jumps& jumps = Jumps())
  {
    using T = scalar_type;
    const std::size_t n = store.n_states();
    const std::size_t n_steps = store.n_steps();

    zero_armed(m_wp, n_phi);
    m_wx.assign(n, T(0.0));
    m_whist.clear();
    m_lam.assign(m_trace ? n_steps * n : 0u, T(0.0));
    m_eta.assign(m_trace ? n_steps : 0u, T(0.0));
    m_ref_m.assign(m_trace ? n_steps : 0u, 1);
    m_ref_fail = 0;
    if (n_steps == 0) {
      sweep_without_steps(store, n, seeds, adj, jumps, m_wx, m_wp, m_jws);
      return;
    }

    std::size_t next_obs = store.n_obs();
    std::vector<T> w_out(n, T(0.0)), w_in(n, T(0.0)), x0, w_after, w_before;
    jump_handover<T> handover;

    for (std::size_t k = n_steps; k-- > 0;) {
      const auto& cp = store.step(k);
      store.use_step_modes(k);
      x0.assign(cp.start_state(), cp.start_state() + n);
      step_seeds(store, k, seeds, n, next_obs, handover, m_obs);

      T eta_k = T(0.0);
      if (m_refine_rtol > 0.0) {
        refined_step(sys, adj, st, n, n_phi, x0, cp.t, cp.dt, w_out, w_in, eta_k, k);
      } else {
        w_in.assign(n, T(0.0));
        step_adjoint(sys, adj, st, n, n_phi, x0.data(), cp.t, cp.dt,
                     m_obs.data(), m_obs.size(), w_out.data(), w_in.data(),
                     m_wp.data(), &eta_k);
      }

      // lambda is what this step hands the one below it. eta is the cotangent
      // at the step end against the step's own embedded error, which a
      // one-step method reports already scaled.
      if (m_trace) {
        for (std::size_t i = 0; i < n; ++i) m_lam[k * n + i] = w_in[i];
        m_eta[k] = eta_k;
      }
      w_out.swap(w_in);

      // The intervention this step was entered through. A one-step method
      // passes on only the state, so there is no restart to collapse: what the
      // step hands down is already the cotangent on the state the jump left.
      if constexpr (Jumps::active) {
        const std::size_t ei = store.event_before(k);
        if (ei < store.n_events()) {
          w_after = w_out;
          event_adjoint(store, ei, seeds, n, adj, jumps, w_after, w_before,
                        m_wp.data(), m_jws);
          w_out.assign(n, T(0.0));
          handover.take(store.event(ei), w_before, m_wx);
        }
      }
    }

    for (std::size_t i = 0; i < n; ++i) m_wx[i] += w_out[i];
    seed_initial_state(store, next_obs, seeds, n, m_wx);
  }

  void report_profile() const { m_prof.report("cppDE written adjoint"); }

  const std::vector<scalar_type>& wx0() const { return m_wx; }
  const std::vector<scalar_type>& whistory0() const { return m_whist; }
  const std::vector<scalar_type>& wp() const { return m_wp; }

  /// Under trace_lambda: [n_steps, n_states] step-major, and one per step.
  const std::vector<scalar_type>& lambda() const { return m_lam; }
  const std::vector<scalar_type>& eta() const { return m_eta; }

private:
  using obs_ref = adjoint::obs_ref<scalar_type>;

  /// The adjoint of the step from `xs` at `t0` over `dt`, its seeds included:
  /// `w_end` is the cotangent at the step's end, `w_in` receives the one at
  /// its start, `wp` the step's share of the gradient. `eta`, if given, gets
  /// the cotangent at the end against the step's embedded error.
  template<class System, class AdjTerms>
  void step_adjoint(System& sys, const AdjTerms& adj, Stepper& st, std::size_t n,
                    std::size_t n_phi, const scalar_type* xs, double t0, double dt,
                    const obs_ref* ob, std::size_t nob, const scalar_type* w_end,
                    scalar_type* w_in, scalar_type* wp, scalar_type* eta,
                    scalar_type* w_in_emb = nullptr, scalar_type* wp_emb = nullptr)
  {
    using T = scalar_type;
    constexpr bool rosen = has_stage_vectors<Stepper>::value;
    m_x0.assign(xs, xs + n);
    // An explicit interpolant reads f at both ends of the step, so its stages
    // must exist before an observation is seeded; a Rosenbrock interpolant is
    // linear in its stage vectors and needs only the weights.
    if constexpr (!rosen) {
      zero_armed(m_ws.xout, n);
      zero_armed(m_ws.xerr, n);
      st.do_step(sys, m_x0, t0, m_ws.xout, dt, m_ws.xerr);
    }
    m_wstart.assign(n, T(0.0));
    m_wend.assign(w_end, w_end + n);
    if constexpr (rosen) m_mu_seed.assign(stage_slots(n), T(0.0));
    zero_armed(m_jv, n);
    const double lo = dt > 0 ? t0 : t0 + dt, hi = dt > 0 ? t0 + dt : t0;
    for (std::size_t o = 0; o < nob; ++o) {
      auto _tp = m_prof.timer(cppde::prof_cat::rev_interp);
      const double t_obs = ob[o].t < lo ? lo : (ob[o].t > hi ? hi : ob[o].t);
      const T* w = ob[o].w;
      const double s = (t_obs - t0) / dt;
      double cw[4];
      Stepper::dense_weights(s, cw);
      for (std::size_t i = 0; i < n; ++i) {
        m_wstart[i] += cw[0] * w[i];
        m_wend[i]   += cw[1] * w[i];
      }
      if constexpr (rosen) {
        // The other two weights sit on the continuous extension's own two
        // vectors, each a fixed combination of the stages.
        for (int j = 1; j <= Stepper::n_stages_used; ++j) {
          const double c = Stepper::dense_stage_weight(3, j) * cw[2]
                         + Stepper::dense_stage_weight(4, j) * cw[3];
          if (c == 0.0) continue;
          T* mj = m_mu_seed.data() + static_cast<std::size_t>(j - 1) * n;
          for (std::size_t i = 0; i < n; ++i) mj[i] += c * w[i];
        }
      } else {
        // k1 = f(x_old, t) and k7 = f(x_new, t + dt), so the two derivative
        // weights land on the states at the two ends and on theta.
        seed_through_rhs(adj, m_x0, t0, w, dt * cw[2], n, n_phi,
                         m_wstart.data(), m_jv, wp);
        seed_through_rhs(adj, m_ws.xout, t0 + dt, w, dt * cw[3], n,
                         n_phi, m_wend.data(), m_jv, wp);
      }
    }
    { auto _tp = m_prof.timer(cppde::prof_cat::rev_adjoint);
      if constexpr (rosen)
        apply_rosenbrock_adjoint(st, sys, m_x0.data(), t0, dt, n, n_phi,
                                 m_wend.data(), adj, w_in, wp, m_rws,
                                 m_mu_seed.data(), static_cast<scalar_type*>(nullptr),
                                 w_in_emb, wp_emb);
      else
        onestep_recurse(st, m_x0.data(), t0, dt, n, n_phi, m_wend.data(),
                        adj, w_in, wp, m_ws); }
    for (std::size_t i = 0; i < n; ++i) w_in[i] += m_wstart[i];
    if (eta) {
      const std::vector<T>& xe = rosen ? m_rws.xerr : m_ws.xerr;
      T e = T(0.0);
      for (std::size_t i = 0; i < n; ++i) e += m_wend[i] * xe[i];
      *eta = e;
    }
  }

  /// The step's adjoint, checked per component of its share of the gradient:
  /// by the embedded estimate where the method has one, then by halving, the
  /// step taken as m substeps of the same method from its own start, m
  /// doubled until the share agrees with that of m / 2 after Richardson
  /// extrapolation at the method's order. The finest result is kept.
  template<class System, class AdjTerms>
  void refined_step(System& sys, const AdjTerms& adj, Stepper& st, std::size_t n,
                    std::size_t n_phi, const std::vector<scalar_type>& x0, double t0,
                    double dt, const std::vector<scalar_type>& w_end,
                    std::vector<scalar_type>& w_in, scalar_type& eta, std::size_t k)
  {
    using T = scalar_type;
    using ad_traits::scalar_value;
    constexpr bool rosen = has_stage_vectors<Stepper>::value;
    zero_armed(m_pc, n_phi);
    w_in.assign(n, T(0.0));
    const double rt = m_refine_rtol;
    // The local error test of CVODES' backward problem: the WRMS norm of the
    // cotangent handed down against rtol |lambda| + atol, and under a gradient
    // tolerance that of the step's share against rtol |running sum| + gradtol.
    auto measure = [&](const std::vector<T>& pa, const std::vector<T>& pb,
                       const std::vector<T>& share, const std::vector<T>& la,
                       const std::vector<T>& lb, const std::vector<T>& lam, bool diff) {
      auto dif = [&](const std::vector<T>& u, const std::vector<T>& v, std::size_t i) {
        return diff ? std::abs(static_cast<double>(scalar_value(u[i]) - scalar_value(v[i])))
                    : std::abs(static_cast<double>(scalar_value(u[i])));
      };
      double sl = 0.0, sq = 0.0;
      for (std::size_t i = 0; i < n; ++i) {
        const double w = rt * std::abs(static_cast<double>(scalar_value(lam[i]))) + m_refine_atol;
        if (w > 0.0) { const double r = dif(la, lb, i) / w; sl += r * r; }
      }
      double e = n > 0 ? std::sqrt(sl / static_cast<double>(n)) : 0.0;
      if (m_refine_gradtol > 0.0 && n_phi > n) {
        for (std::size_t i = n; i < n_phi; ++i) {
          const double q = static_cast<double>(scalar_value(m_wp[i]) + scalar_value(share[i]));
          const double r = dif(pa, pb, i) / (rt * std::abs(q) + m_refine_gradtol);
          sq += r * r;
        }
        e = std::max(e, std::sqrt(sq / static_cast<double>(n_phi - n)));
      }
      return e;
    };
    double est = 0.0;
    if constexpr (rosen) {
      // The embedded estimate first, at no new step: a step that meets the
      // tolerance at the lower order is taken as it is.
      zero_armed(m_pe, n_phi);
      m_we.assign(n, T(0.0));
      step_adjoint(sys, adj, st, n, n_phi, x0.data(), t0, dt, m_obs.data(),
                   m_obs.size(), w_end.data(), w_in.data(), m_pc.data(), &eta,
                   m_we.data(), m_pe.data());
      est = measure(m_pe, m_pe, m_pc, m_we, m_we, w_in, false);
      if (est <= 1.0) {
        for (std::size_t i = 0; i < n_phi; ++i) m_wp[i] += m_pc[i];
        if (m_trace && k < m_ref_m.size()) m_ref_m[k] = 1;
        return;
      }
    } else {
      step_adjoint(sys, adj, st, n, n_phi, x0.data(), t0, dt, m_obs.data(),
                   m_obs.size(), w_end.data(), w_in.data(), m_pc.data(), &eta);
    }
    // Richardson's factor 2^q - 1 only as far as two successive differences
    // show the order: the first halving is taken at face value.
    const double R = std::pow(2.0, static_cast<double>(Stepper::stepper_order)) - 1.0;
    double d_prev = 0.0;
    int m = 1;
    while (m < m_refine_max) {
      m *= 2;
      const double h = dt / m;
      // The substeps forward, from the step's own start.
      m_xs.resize(static_cast<std::size_t>(m + 1));
      m_xs[0] = x0;
      for (int j = 0; j < m; ++j) {
        zero_armed(m_xs[static_cast<std::size_t>(j + 1)], n);
        zero_armed(m_xe, n);
        st.do_step(sys, m_xs[static_cast<std::size_t>(j)], t0 + j * h,
                   m_xs[static_cast<std::size_t>(j + 1)], h, m_xe);
      }
      // And back, each with the seeds that fall into it: those before the
      // step's start into the first substep, those past its end into the last,
      // in the direction of integration.
      zero_armed(m_pf, n_phi);
      m_wa = w_end;
      const double dir = dt > 0 ? 1.0 : -1.0;
      for (int j = m - 1; j >= 0; --j) {
        const double a = t0 + j * h, b = t0 + (j + 1) * h;
        m_sub.clear();
        for (const obs_ref& o : m_obs)
          if ((j == 0 || dir * (o.t - a) > 0) && (j == m - 1 || dir * (o.t - b) <= 0))
            m_sub.push_back(o);
        m_wb.assign(n, T(0.0));
        step_adjoint(sys, adj, st, n, n_phi, m_xs[static_cast<std::size_t>(j)].data(),
                     a, h, m_sub.data(), m_sub.size(), m_wa.data(), m_wb.data(),
                     m_pf.data(), nullptr);
        m_wa.swap(m_wb);
      }
      // The finer result against the coarser, per component.
      const double d = measure(m_pf, m_pc, m_pf, m_wa, w_in, m_wa, true);
      const double r = (m > 2 && d > 0.0) ? std::min(R, std::max(1.0, d_prev / d - 1.0)) : 1.0;
      est = d / r;
      // A halving that no longer shrinks the difference at all meets an error
      // the substeps do not reduce, and more of them buy nothing.
      const bool stalled = m >= stall_substeps && d >= d_prev;
      d_prev = d;
      m_pc.swap(m_pf);
      w_in.swap(m_wa);
      if (est <= 1.0 || stalled) break;
    }
    for (std::size_t i = 0; i < n_phi; ++i) m_wp[i] += m_pc[i];
    if (est > 1.0) ++m_ref_fail;
    if (m_trace && k < m_ref_m.size()) m_ref_m[k] = m;
  }

  /// A cotangent `w` scaled by `c` on f(x, t), put onto x and theta.
  template<class AdjTerms>
  void seed_through_rhs(const AdjTerms& adj, const std::vector<scalar_type>& x,
                        double t, const scalar_type* w, double c, std::size_t n,
                        std::size_t n_phi, scalar_type* w_x,
                        std::vector<scalar_type>& jv, scalar_type* wp)
  {
    if (c == 0.0) return;
    m_ws.m.assign(n, scalar_type(0.0));
    for (std::size_t i = 0; i < n; ++i) m_ws.m[i] = c * w[i];
    adj.jac_t_vec(x, m_ws.m, t, jv);
    for (std::size_t i = 0; i < n; ++i) w_x[i] += jv[i];
    adj.dfdp_t_vec_axpy(x, m_ws.m, t, 1.0, wp);
    (void)n_phi;
  }

  /// How many stage cotangents a step has, when it has any.
  static std::size_t stage_slots(std::size_t n) {
    if constexpr (has_stage_vectors<Stepper>::value)
      return static_cast<std::size_t>(Stepper::n_stages_used) * n;
    else
      return 0u;
  }

  cppde::profiler m_prof;
  jump_workspace<scalar_type> m_jws;
  onestep_workspace<scalar_type> m_ws;
  rosenbrock_workspace<scalar_type> m_rws;
  bool m_trace = false;
  std::vector<scalar_type> m_wx, m_whist, m_wp, m_lam, m_eta;
  // The step adjoint's buffers, and the halving check's.
  std::vector<obs_ref> m_obs, m_sub;
  std::vector<scalar_type> m_x0, m_wstart, m_wend, m_mu_seed, m_jv, m_xe,
                           m_pc, m_pf, m_wa, m_wb, m_pe, m_we;
  std::vector<std::vector<scalar_type>> m_xs;
  double m_refine_rtol = 0.0;
  double m_refine_atol = 0.0, m_refine_gradtol = 0.0;
  int m_refine_max = default_max_sub;
  std::size_t m_ref_fail = 0;
  std::vector<int> m_ref_m;
};

}  // namespace adjoint
}  // namespace cppde

#endif  // CPPDE_ADJOINT_STEP_HPP
