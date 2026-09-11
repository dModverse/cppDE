/**
 * @file cppde_adjoint_step.hpp
 * @brief The adjoint of one accepted step, written rather than recorded.
 *
 * The reverse mode used to obtain a step's adjoint by re-running the step on a
 * tape type and sweeping the tape. On a stiff model that cost about three
 * hundred right-hand-side evaluations per step, for a step whose forward
 * variant costs twenty-seven. It bought generality that is not needed: the
 * accepted grid is frozen, so a step is a fixed map and its adjoint can be
 * stated.
 *
 * What a step does to the Nordsieck history splits in three:
 *
 *   zn_pred = A zn_in                      rescale and Pascal shift
 *   res(y, zn_pred, theta) = 0             the corrector, closed by the IFT
 *   zn_out  = B zn_pred + c acor           the tail, acor = y - zn_pred[0]
 *
 * A, B and c act on the slot index alone, with plain-double coefficients that
 * follow from the carry. They are therefore the same small matrices for every
 * state component, and they are obtained from the stepper itself: the same
 * routines the forward step runs, applied to a unit slot on a one-state copy.
 * Nothing here restates what the stepper does, so nothing here can drift from
 * it.
 *
 * The model supplies two contractions, which the code generator emits:
 * `jac_t_vec` for J' lambda, which sizes its own output, and
 * `dfdp_t_vec_axpy` for (df/dp)' lambda, scaled and added into what the caller
 * already holds. So nothing here has to know the model's dimensions.
 *
 * Copyright (C) 2026 Simon Beyer
 */

#ifndef CPPDE_ADJOINT_STEP_HPP
#define CPPDE_ADJOINT_STEP_HPP

#include <cmath>
#include <cstddef>
#include <type_traits>
#include <vector>

#include <cppde/cppde_dual_math.hpp>
#include <cppde/cppde_events.hpp>
#include <cppde/cppde_profiler.hpp>

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
//  n zeros with their tangent storage bound.
//
//  Everything the sweep hands to the model or to the stepper goes through
//  this. Those open a dual_arena::scope of their own, and a tangent allocated
//  inside one is reclaimed when it pops; a buffer that already has one bound
//  is written in place instead. In plain double it is an assign and nothing
//  more. Sizing to exactly what the callee writes matters: a later growth
//  copy-constructs the new elements unarmed.
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
//  Building the operators.
//
//  Every one of them acts on the slot index with coefficients that follow from
//  the carry, identically for each state component. So a probe whose states
//  are the slots themselves gives every column in one pass: put the identity
//  into the [slot x component] array and read the result back as a matrix.
//  One pass instead of one per column, which is what makes reading the
//  operators off the stepper affordable.
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

/// B and the acor column c: what the tail does. Component q+1 carries acor.
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
  // under a sweep that carries tangents, and the matrices come out double.
  using probe_stepper = typename Stepper::template rebind_value<double>;
  probe_stepper pre, tail_probe;

  // What the operators depend on. Long stretches of a run hold the order and
  // the step size, and then every step asks for the same matrices; comparing
  // is a handful of doubles against reading them off again.
  //
  // The carry's qwait is deliberately not here. It counts down to the next
  // order decision and so moves almost every step, but it reaches only the
  // error constants for going up or down an order, never l, gamma or the
  // Nordsieck shift. Those constants belong to the controller, which is not
  // differentiated. Keeping it in would cost every hit in the cache.
  struct key {
    int q = -1, L = 0;
    bool started = false;
    double dt = 0.0, h = 0.0, hscale = 0.0, eta = 0.0;
    double tail = 0.0;
    std::vector<double> tau;
    bool operator==(const key& o) const {
      return q == o.q && L == o.L && started == o.started &&
             dt == o.dt && h == o.h && hscale == o.hscale && eta == o.eta &&
             tail == o.tail && tau == o.tau;
    }
  };
  key last, pending;
  bool valid = false;

  /// Whether the operators already in hand are this step's. The key it built
  /// to answer stays, so the rebuild that may follow does not build it twice.
  template<class Carry>
  bool matches(const Carry& carry, double dt, double tail_key) {
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
//  The adjoint of one step.
//
//    w_pred = B' w_out,  w_acor = c' w_out
//    acor = y - zn_pred[0]        ->  w_y = w_acor,  w_pred[0] -= w_acor
//    res  = (y - zn_pred0) + rl1 zn_pred1 - gamma f(y, t_new) = 0
//         ->  mu = (I - gamma J)^-T w_y
//             w_pred[0] += mu,  w_pred[1] -= rl1 mu
//             w_theta   += gamma (df/dp)' mu
//    w_in = A' w_pred
//
//  `solver` is the same equation_solver the tape path uses: prepared on the
//  step's own Jacobian, its transposed apply already carries the gamma scale.
// ---------------------------------------------------------------------------
/// What the step above hands down, carried back through the tail onto the
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

  adj.dfdp_t_vec_axpy(y, ws.mu, t_new, ops.gamma, w_theta);
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
//  The adjoint of one jump.
//
//  The engine carries a discontinuity across on a Heun sandwich: a forward
//  shift to the event surface, the resets, a backward shift to the grid time.
//  The shift is by the event time's own residual, whose value is subtracted off
//  before it is used, so it is numerically zero and only its derivative
//  survives. In value the jump is therefore a jump, and the sandwich collapses
//  to the classical saltation:
//
//    dx = R'(dx_before + f_before s) - f_after s
//
//  with R the resets and s the residual's differential. For a root event
//  s = -(grad g . dx + dg/dp . dp) / g_dot, for a fixed one it is the
//  differential of the event's time. The second-order correction of the root's
//  dt* multiplies dt* itself and drops with it.
//
//  So a jump's adjoint needs the right-hand side at the two ends and the model's
//  own derivatives of the event expressions. It needs no Jacobian: every term
//  the shifts would have contributed carries a factor that is zero.
// ---------------------------------------------------------------------------
template<class T = double>
struct jump_workspace {
  std::vector<T> fb, fa, g, wy;
  std::vector<std::vector<T> > path;
};

/// One reset, transposed. `w_z` is the cotangent on the state it wrote, `w_y`
/// the one on the state it read, which starts as a copy of `w_z`.
///     Replace   z[k] = h(y, t)
///     Add       z[k] = y[k] + h(y, t)
///     Multiply  z[k] = y[k] * h(y, t)
template<class GradX, class GradP, class T>
void reset_transpose(int k, cppde::detail::EventMethod method, const T& h,
                     int idx,
                     const std::vector<T>& y, double t, std::size_t n,
                     const T* w_z, T* w_y, T* w_theta,
                     GradX&& dh_dx, GradP&& dh_dp_axpy, std::vector<T>& g)
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
}

/// A batch of root events, which share one surface and one dt*.
template<class System, class RootEvents, class EvAdj, class T>
void apply_root_jump_adjoint(const std::vector<T>& x_before,
                             const std::vector<T>& x_after,
                             double t, const RootEvents& root_events,
                             const std::vector<cppde::detail::TriggeredEvent>& triggered,
                             System& sys, const EvAdj& eadj, std::size_t n,
                             const T* w_out, T* w_in, T* w_theta,
                             jump_workspace<T>& ws)
{
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

  T w_s = T(0.0);
  if (shifted) {
    zero_armed(ws.fa, n);
    sys.first(x_after, ws.fa, t);
    for (std::size_t i = 0; i < n; ++i) w_s -= ws.fa[i] * w_out[i];
  }

  // The resets, every one reading the same pre-jump state.
  for (std::size_t i = 0; i < n; ++i) w_in[i] = w_out[i];
  for (std::size_t j = triggered.size(); j-- > 0;) {
    const std::size_t idx = triggered[j].index;
    const auto& evt = root_events[idx];
    if (evt.terminal) continue;
    const T h = (evt.state_index >= 0 && evt.value_func)
                   ? evt.value_func(x_before, t) : T(0.0);
    reset_transpose(
        evt.state_index, evt.method, h, static_cast<int>(idx), x_before, t, n,
        w_out, w_in, w_theta,
        [&](int e, const std::vector<T>& y, double tt, std::vector<T>& o)
          { eadj.root_dh_dx(e, y, tt, o); },
        [&](int e, const std::vector<T>& y, double tt, const T& sc, T* o)
          { eadj.root_dh_dp_axpy(e, y, tt, sc, o); },
        ws.g);
  }

  if (!shifted) return;
  for (std::size_t i = 0; i < n; ++i) w_s += ws.fb[i] * w_in[i];

  // s = -(grad g . dx + dg/dp . dp) / g_dot
  const std::size_t idx = triggered[src].index;
  const T c = -w_s / g_dot;
  zero_armed(ws.g, n);
  root_events[idx].dg_dx(x_before, t, ws.g);
  for (std::size_t i = 0; i < n; ++i) w_in[i] += c * ws.g[i];
  eadj.root_dg_dp_axpy(static_cast<int>(idx), x_before, t, c, w_theta);
}

/// The fixed events at one time, each its own sandwich, applied in order. The
/// last of them carries the root resets a jump switched on.
template<class System, class FixedEvents, class RootEvents, class EvAdj, class T>
void apply_fixed_jump_adjoint(const std::vector<T>& x_before,
                              double t, const FixedEvents& fixed_events,
                              const RootEvents& root_events,
                              const std::vector<std::size_t>& switched,
                              System& sys, const EvAdj& eadj, std::size_t n,
                              const T* w_out, T* w_in, T* w_theta,
                              jump_workspace<T>& ws)
{
  // Which of them fire here, and which is last: the same test the engine makes.
  std::vector<int> fired;
  for (std::size_t j = 0; j < fixed_events.size(); ++j)
    if (std::abs(ad_traits::scalar_value(fixed_events[j].time) - t) < 1e-14)
      fired.push_back(static_cast<int>(j));

  // The value path through the resets, which the store does not keep: a jump of
  // several events runs one sandwich each, and every one reads what the last
  // left. In value a sandwich is its reset, so this is the reset chain.
  const std::size_t n_steps = fired.size() + switched.size();
  ws.path.assign(n_steps + 1, std::vector<T>());
  ws.path[0] = x_before;
  std::size_t s = 0;
  for (std::size_t j = 0; j < fired.size(); ++j, ++s) {
    zero_armed(ws.path[s + 1], n);
    ws.path[s + 1] = ws.path[s];
    cppde::detail::apply_event_action_fixed(ws.path[s + 1], ws.path[s],
                                            fixed_events[fired[j]]);
  }
  for (std::size_t j = 0; j < switched.size(); ++j, ++s) {
    zero_armed(ws.path[s + 1], n);
    ws.path[s + 1] = ws.path[s];
    cppde::detail::apply_event_action(ws.path[s + 1], ws.path[s], T(t),
                                      root_events[switched[j]]);
  }

  // Backwards through the same chain. The switched resets ride on the last
  // sandwich's surface, so they sit inside its shift rather than beside it.
  std::vector<T> w(w_out, w_out + n);
  ws.wy.assign(n, T(0.0));
  for (std::size_t j = switched.size(); j-- > 0;) {
    const std::size_t p = fired.size() + j;
    const auto& evt = root_events[switched[j]];
    const T h = (evt.state_index >= 0 && evt.value_func)
                   ? evt.value_func(ws.path[p], t) : T(0.0);
    ws.wy = w;
    reset_transpose(
        evt.state_index, evt.method, h, static_cast<int>(switched[j]),
        ws.path[p], t, n, w.data(), ws.wy.data(), w_theta,
        [&](int e, const std::vector<T>& y, double tt, std::vector<T>& o)
          { eadj.root_dh_dx(e, y, tt, o); },
        [&](int e, const std::vector<T>& y, double tt, const T& sc, T* o)
          { eadj.root_dh_dp_axpy(e, y, tt, sc, o); },
        ws.g);
    w.swap(ws.wy);
  }

  for (std::size_t j = fired.size(); j-- > 0;) {
    const auto& evt = fixed_events[fired[j]];
    const std::vector<T>& y = ws.path[j];
    const std::vector<T>& z = ws.path[j + 1];

    zero_armed(ws.fb, n); zero_armed(ws.fa, n);
    sys.first(y, ws.fb, t);
    // The last sandwich ends on the surface the switched resets left.
    const std::vector<T>& zz =
        (j + 1 == fired.size()) ? ws.path[n_steps] : z;
    sys.first(zz, ws.fa, t);

    T w_s = T(0.0);
    for (std::size_t i = 0; i < n; ++i) w_s -= ws.fa[i] * w[i];

    const T h = (evt.state_index >= 0 && evt.value_func)
                   ? evt.value_func(y, evt.time) : T(0.0);
    ws.wy = w;
    reset_transpose(
        evt.state_index, evt.method, h, fired[j], y, t, n,
        w.data(), ws.wy.data(), w_theta,
        [&](int e, const std::vector<T>& yy, double tt, std::vector<T>& o)
          { eadj.fixed_dh_dx(e, yy, tt, o); },
        [&](int e, const std::vector<T>& yy, double tt, const T& sc, T* o)
          { eadj.fixed_dh_dp_axpy(e, yy, tt, sc, o); },
        ws.g);
    w.swap(ws.wy);

    for (std::size_t i = 0; i < n; ++i) w_s += ws.fb[i] * w[i];
    eadj.fixed_dtime_dp_axpy(fired[j], w_s, w_theta);
  }

  for (std::size_t i = 0; i < n; ++i) w_in[i] = w[i];
}

// ---------------------------------------------------------------------------
//  A whole trajectory backwards, without a tape.
//
//  The same store the recorder walks, the same order, the same outputs. What
//  changes is what happens inside a step: the operators are read off the
//  stepper and the adjoint is applied, rather than the step being re-run on a
//  tape type and the tape swept.
//
//  Observations inside a step reach it through the dense output, and their row
//  is read off the same probe: after the tail the probe carries a valid
//  interpolant over B, so evaluating it at the observation time gives
//  d x_interp / d (zn_pred, acor) directly. The probe's interpolant belongs to
//  the step whose tail last ran on it, so a step that is observed rebuilds
//  rather than taking the cached operators.
//
//  Events are not here yet: a jump is its own map between two steps and comes
//  with the saltation adjoint.
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
    if (n_steps == 0) return;

    // The cotangent the step above hands down, on its own carry. Across an
    // intervention it is instead a cotangent on a point inside the step below,
    // because the restart threw the carry away.
    std::vector<T> w_carry, pending, w_after, w_before;
    bool pending_interp = false;
    double pending_t = 0.0;
    std::size_t next_obs = store.n_obs();

    multistep_operators<Stepper> ops;
    // A dense row is one entry per Nordsieck slot, a contraction one per
    // state. Two lengths, two buffers, and the row is the interpolant's own
    // coefficients, which are double whatever the sweep carries.
    std::vector<double> dense_row;
    std::vector<T> w_in, jtv;
    zero_armed(jtv, n);

    for (std::size_t k = n_steps; k-- > 0;) {
      const auto& cp = store.step(k);

      // Does anything observe inside this step? Then the probe has to carry
      // this step's own interpolant, not the one it kept from another.
      std::size_t obs_lo = next_obs;
      while (obs_lo > 0 && store.obs(obs_lo - 1).step == k + 1) --obs_lo;
      const bool observed = (obs_lo < next_obs) || pending_interp;
      if (observed) m_probe.valid = false;

      const double tail_key =
          cp.q_next + 1e3 * cp.eta + 1e6 * static_cast<double>(cp.ops.size());
      // Timed on the rebuild alone, so the call count is the miss count.
      if (!m_probe.matches(cp.carry, cp.dt, tail_key)) {
        auto _tp = m_prof.timer(cppde::prof_cat::rev_operators);
        m_probe.rebuild(cp.carry, cp.dt,
                        [&](typename multistep_probe<Stepper>::probe_stepper& pr)
                          { cp.apply_tail(pr, m_null); }, ops);
      }

      const std::size_t nz_in = static_cast<std::size_t>(ops.q_in + 1) * n;
      const std::size_t nz_out = static_cast<std::size_t>(ops.q_out + 1) * n;
      m_ws.w_pred.assign(nz_in, T(0.0));
      m_ws.w_acor.assign(n, T(0.0));

      // A cotangent at a time inside this step, through its own interpolant.
      // The probe's states are the slots plus one for acor, so the row it
      // writes is d x_interp / d (zn_pred[0..q], acor).
      auto seed_dense = [&](double t_obs, const T* w) {
        auto _tp = m_prof.timer(cppde::prof_cat::rev_interp);
        dense_row.assign(static_cast<std::size_t>(ops.q_in + 2), 0.0);
        m_probe.tail_probe.eval_dense_into(clamp_to_step(cp, t_obs), dense_row);
        for (int j = 0; j <= ops.q_in; ++j)
          for (std::size_t i = 0; i < n; ++i)
            m_ws.w_pred[static_cast<std::size_t>(j) * n + i] +=
                dense_row[static_cast<std::size_t>(j)] * w[i];
        const std::size_t acor_slot = dense_row.size() - 1;
        for (std::size_t i = 0; i < n; ++i)
          m_ws.w_acor[i] += dense_row[acor_slot] * w[i];
      };

      // What the step above handed back. Ordinarily its carry, read off this
      // step's end; across an intervention a cotangent on a point inside this
      // step, and then the carry is not read at all.
      if (pending_interp) {
        seed_dense(pending_t, pending.data());
        pending_interp = false;
      } else if (!w_carry.empty()) {
        auto _tp = m_prof.timer(cppde::prof_cat::rev_adjoint);
        w_carry.resize(nz_out, T(0.0));
        carry_into_pred(ops, n, w_carry.data(), m_ws.w_pred.data(),
                        m_ws.w_acor.data());
      }

      // The observations this step carries. One a jump produced is a value and
      // is seeded at the boundary instead.
      for (std::size_t o = obs_lo; o < next_obs; ++o) {
        if (store.obs(o).event < store.n_events()) continue;
        seed_dense(store.obs(o).t, seeds + o * n);
      }
      next_obs = obs_lo;

      const double t_new = cp.t + ops.h;
      solver.prepare(cp.y, t_new, 1.0 / ops.gamma, ops.gamma);

      w_in.assign(nz_in, T(0.0));
      if (m_trace) m_wout.assign(n, T(0.0));
      // Timed either side of the solve, which reports itself.
      { auto _tp = m_prof.timer(cppde::prof_cat::rev_adjoint);
        multistep_adjoint_rhs(ops, n, m_ws,
                              m_trace ? m_wout.data() : nullptr); }
      solver.transposed(m_ws.mu);
      { auto _tp = m_prof.timer(cppde::prof_cat::rev_adjoint);
        multistep_adjoint_finish(ops, n, n_phi, cp.y, t_new, adj, w_in.data(),
                                 m_wp.data(), m_ws); }

      // lambda is what this step hands the one below it, on the state slot.
      // eta is that cotangent against the step's own error estimate, which for
      // a corrector method is acor scaled by the order's error constant.
      if (m_trace) {
        for (std::size_t i = 0; i < n; ++i) m_lam[k * n + i] = w_in[i];
        T e = T(0.0);
        for (std::size_t i = 0; i < n; ++i) {
          T pred0 = T(0.0);
          for (int j = 0; j <= ops.q_in; ++j)
            pred0 += ops.A(0, j) * cp.zn[static_cast<std::size_t>(j) * n + i];
          e += m_wout[i] * (cp.y[i] - pred0);
        }
        m_eta[k] = e * ops.err_scale;
      }
      w_carry.swap(w_in);

      // The intervention this step was entered through, if any. Two maps in the
      // order the forward run applied them, so swept the other way round: the
      // restart, which collapses the carry onto the state the jump ended on,
      // and the jump itself.
      if constexpr (Jumps::active) {
        const std::size_t ei = store.event_before(k);
        if (ei < store.n_events()) {
          const auto& e = store.event(ei);
          w_after.assign(n, T(0.0));
          m_ws.x.assign(e.x_after.begin(), e.x_after.end());
          collapse_restart(w_carry, n, m_ws.x, e.t,
                           e.restart ? e.dt_restart
                                     : static_cast<double>(cp.carry.h),
                           adj, w_after.data(), m_wp.data(), m_ws.mu, jtv);

          // Observations the jump produced are values, not interpolations.
          for (std::size_t o = 0; o < store.n_obs(); ++o)
            if (store.obs(o).event == ei)
              for (std::size_t i = 0; i < n; ++i)
                w_after[i] += seeds[o * n + i];

          w_before.assign(n, T(0.0));
          if (e.root)
            apply_root_jump_adjoint(e.x_before, e.x_after, e.t, jumps.root,
                                    e.triggered, jumps.sys, jumps.eadj, n,
                                    w_after.data(), w_before.data(),
                                    m_wp.data(), m_jws);
          else
            apply_fixed_jump_adjoint(e.x_before, e.t, jumps.fixed, jumps.root,
                                     e.switched, jumps.sys, jumps.eadj, n,
                                     w_after.data(), w_before.data(),
                                     m_wp.data(), m_jws);

          // The state entering the jump was read off the dense output of the
          // step below it. At the run's own start there is no such step and it
          // is the initial state's.
          w_carry.clear();
          if (e.after_step > 0) {
            pending = w_before;
            pending_t = e.t_before;
            pending_interp = true;
          } else {
            for (std::size_t i = 0; i < n; ++i) m_wx[i] += w_before[i];
          }
        }
      }
    }

    // The trajectory start, which is the same restart with no jump under it.
    if (!w_carry.empty()) {
      const auto& cp0 = store.step(0);
      w_after.assign(n, T(0.0));
      m_ws.x.assign(cp0.start_state(), cp0.start_state() + n);
      collapse_restart(w_carry, n, m_ws.x, cp0.t,
                       static_cast<double>(cp0.carry.h), adj, w_after.data(),
                       m_wp.data(), m_ws.mu, jtv);
      for (std::size_t i = 0; i < n; ++i) m_wx[i] += w_after[i];
    }

    // Anything observed before the first step is the initial state itself.
    while (next_obs > 0) {
      --next_obs;
      const T* w = seeds + next_obs * n;
      for (std::size_t i = 0; i < n; ++i) m_wx[i] += w[i];
    }
  }

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
  bool m_trace = false;
  std::vector<scalar_type> m_wx, m_whist, m_wp, m_lam, m_eta, m_wout;
};

// ---------------------------------------------------------------------------
//  The adjoint of one explicit Runge-Kutta step.
//
//    X_i    = x + h sum_{j<i} a_ij k_j,   k_i = f(X_i, t + c_i h)
//    x_out  = x + h sum_i b_i k_i
//
//  and backwards, with u_j = J(X_j, t_j)' m_j,
//
//    m_i    = h b_i w + h sum_{j>i} a_ji u_j
//    w_x    = w + sum_i u_i
//    w_th  += sum_i (df/dp)(X_i, t_i)' m_i
//
//  The stage states are not checkpointed, so the step is run forward once in
//  plain double to recover them. That is the same work the forward step did,
//  and it is why an explicit method's adjoint costs about twice its step
//  rather than the three hundred right-hand sides a tape costs.
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
//  The adjoint of one Rosenbrock step.
//
//  Six solves against one matrix W = I/(gamma h) - J(x, t):
//
//    X_i   = x + sum_{j<i} a_ij g_j,     F_i = f(X_i, t + node_i h)
//    g_i   = W^-1 ( F_i + h d_i D + sum_{j<i} (c_ij/h) g_j )
//    X_6   = X_5 + g_5,   e = W^-1 ( F_6 + sum_j (c_6j/h) g_j )
//    x_out = X_6 + e
//
//  with F_1 = f(x, t) and D = df/dt(x, t).
//
//  W is built from the Jacobian, so a stage depends on x and theta through the
//  matrix as well as through its right-hand side: g = W^-1 r gives
//  dg = W^-1 (dr + (dJ) g), whose adjoint is the contraction of lambda with the
//  derivative of J g. That is a second derivative, which a tape supplied
//  silently and a written adjoint has to ask the model for.
//
//  The transposed solves go through the stepper's own factorisation, which is
//  the one the stages solved against. A corrector method may not do that, its
//  iteration matrix being stale by design; a Rosenbrock stage is a direct solve
//  and the matrix it used is the matrix its derivative needs.
// ---------------------------------------------------------------------------
template<class T = double>
struct rosenbrock_workspace {
  std::vector<T> xout, xerr, x, lam, wx, wD, jv, xs;
  std::vector<std::vector<T> > mu, wX;
};

template<class Stepper, class System, class AdjTerms, class T>
void apply_rosenbrock_adjoint(Stepper& st, System& sys,
                              const T* x0, double t, double dt,
                              std::size_t n, std::size_t n_phi,
                              const T* w_out,
                              const AdjTerms& adj,
                              T* w_in, T* w_theta,
                              rosenbrock_workspace<T>& ws,
                              const T* mu_seed = nullptr,
                              T* w_out_state = nullptr)
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
  auto sweep_stage = [&](int i, const std::vector<T>& gi, T* wXi) {
    ws.lam = ws.mu[i];
    st.stage_solve_transposed(ws.lam);

    adj.jvp_x_t_vec(ws.x, gi, ws.lam, t, ws.jv);
    for (std::size_t k = 0; k < n; ++k) ws.wx[k] += ws.jv[k];
    adj.jvp_p_t_vec_axpy(ws.x, gi, ws.lam, t, 1.0, w_theta);

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
      adj.jac_t_vec(ws.x, ws.lam, ti, ws.jv);
      for (std::size_t k = 0; k < n; ++k) ws.wx[k] += ws.jv[k];
      adj.dfdp_t_vec_axpy(ws.x, ws.lam, ti, 1.0, w_theta);
      return;
    }
    stage_state(i);
    adj.jac_t_vec(ws.xs, ws.lam, ti, ws.jv);
    for (std::size_t k = 0; k < n; ++k) wXi[k] += ws.jv[k];
    adj.dfdp_t_vec_axpy(ws.xs, ws.lam, ti, 1.0, w_theta);
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

  // x_out = X_6 + e
  for (std::size_t k = 0; k < n; ++k) {
    ws.wX[S + 1][k] += w_out[k];
    ws.mu[S + 1][k] += w_out[k];
  }
  // The error solve, whose stage vector is the error estimate itself.
  sweep_stage(S + 1, ws.xerr, ws.wX[S + 1].data());
  // X_6 = X_5 + g_5
  for (std::size_t k = 0; k < n; ++k) {
    ws.wX[S][k] += ws.wX[S + 1][k];
    ws.mu[S][k] += ws.wX[S + 1][k];
  }

  for (int i = S; i >= 1; --i) {
    sweep_stage(i, st.stage_g(i), ws.wX[i].data());
    if (i >= 2) spread_stage_state(i, ws.wX[i].data());
  }

  // D = df/dt(x, t), which the Jacobian evaluation filled and every early
  // stage reads.
  if (wD_touched) {
    adj.dfdt_x_t_vec(ws.x, ws.wD, t, ws.jv);
    for (std::size_t k = 0; k < n; ++k) ws.wx[k] += ws.jv[k];
    adj.dfdt_p_t_vec_axpy(ws.x, ws.wD, t, 1.0, w_theta);
  }

  for (std::size_t k = 0; k < n; ++k) w_in[k] = ws.wx[k];
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
//  A whole trajectory backwards on a one-step method, without a tape.
//
//  The same store, the same order, the same outputs as the multistep form, and
//  a simpler shape: a one-step method carries only the state across a step
//  boundary, so there is no history to collapse and no start boundary.
//
//  An observation inside a step reaches it through the continuous extension,
//  which for tsit5 is a Hermite cubic over (x_old, x_new, h k1, h k7) and whose
//  weights the stepper hands out. Two of those four are right-hand sides, so an
//  observation puts a cotangent on f at both ends of the step, and that is one
//  J' and one (df/dp)' contraction apiece.
// ---------------------------------------------------------------------------
template<class Stepper>
class closed_onestep_trajectory {
public:
  // As in the multistep form: the sweep runs in the stepper's scalar type.
  using scalar_type = typename Stepper::value_type;

  /// Whether the sweep also keeps lambda and the refinement indicator per step.
  void trace_lambda(bool on) { m_trace = on; }

  /// One seed column. `seeds` is [n_obs x n_states] row-major.
  template<class Store, class System, class AdjTerms, class Jumps = no_jumps>
  void sweep(const Store& store, std::size_t n_phi, const scalar_type* seeds,
             System& sys, const AdjTerms& adj, Stepper& st,
             const Jumps& jumps = Jumps())
  {
    using T = scalar_type;
    constexpr bool rosen = has_stage_vectors<Stepper>::value;
    const std::size_t n = store.n_states();
    const std::size_t n_steps = store.n_steps();

    zero_armed(m_wp, n_phi);
    m_wx.assign(n, T(0.0));
    m_whist.clear();
    m_lam.assign(m_trace ? n_steps * n : 0u, T(0.0));
    m_eta.assign(m_trace ? n_steps : 0u, T(0.0));
    if (n_steps == 0) return;

    std::size_t next_obs = store.n_obs();
    std::vector<T> w_out(n, T(0.0)), w_in(n, T(0.0)), w_start(n, T(0.0)), jv,
                   x0, mu_seed, pending, w_after, w_before;
    zero_armed(jv, n);
    bool pending_interp = false;
    double pending_t = 0.0;

    for (std::size_t k = n_steps; k-- > 0;) {
      const auto& cp = store.step(k);

      std::size_t obs_lo = next_obs;
      while (obs_lo > 0 && store.obs(obs_lo - 1).step == k + 1) --obs_lo;

      x0.assign(cp.start_state(), cp.start_state() + n);

      // An explicit method's interpolant reads the right-hand side at both ends
      // of the step, so its stages have to be there before an observation can
      // be seeded. A Rosenbrock interpolant is linear in its own stage vectors
      // with constant weights, and needs nothing but the weights.
      if constexpr (!rosen) {
        zero_armed(m_ws.xout, n);
        zero_armed(m_ws.xerr, n);
        st.do_step(sys, x0, cp.t, m_ws.xout, cp.dt, m_ws.xerr);
      }

      w_start.assign(n, T(0.0));
      if constexpr (rosen) mu_seed.assign(stage_slots(n), T(0.0));

      // A cotangent at a time inside this step, through its interpolant.
      auto seed_at = [&](double t_raw, const T* w) {
        auto _tp = m_prof.timer(cppde::prof_cat::rev_interp);
        const double t_obs = clamp_to_step(cp, t_raw);
        const double s = (t_obs - cp.t) / cp.dt;
        double cw[4];
        Stepper::dense_weights(s, cw);

        for (std::size_t i = 0; i < n; ++i) {
          w_start[i] += cw[0] * w[i];
          w_out[i]   += cw[1] * w[i];
        }
        if constexpr (rosen) {
          // The other two weights sit on the continuous extension's own two
          // vectors, each a fixed combination of the stages.
          for (int j = 1; j <= Stepper::n_stages_used; ++j) {
            const double c = Stepper::dense_stage_weight(3, j) * cw[2]
                           + Stepper::dense_stage_weight(4, j) * cw[3];
            if (c == 0.0) continue;
            T* mj = mu_seed.data() + static_cast<std::size_t>(j - 1) * n;
            for (std::size_t i = 0; i < n; ++i) mj[i] += c * w[i];
          }
        } else {
          // k1 = f(x_old, t) and k7 = f(x_new, t + dt), so the two derivative
          // weights land on the states at the two ends and on theta.
          seed_through_rhs(adj, x0, cp.t, w, cp.dt * cw[2], n, n_phi,
                           w_start.data(), jv);
          seed_through_rhs(adj, m_ws.xout, cp.t + cp.dt, w, cp.dt * cw[3], n,
                           n_phi, w_out.data(), jv);
        }
      };

      // What the step above handed back, where a jump put it inside this step.
      if (pending_interp) { seed_at(pending_t, pending.data()); pending_interp = false; }
      for (std::size_t o = obs_lo; o < next_obs; ++o) {
        if (store.obs(o).event < store.n_events()) continue;
        seed_at(store.obs(o).t, seeds + o * n);
      }
      next_obs = obs_lo;

      w_in.assign(n, T(0.0));
      { auto _tp = m_prof.timer(cppde::prof_cat::rev_adjoint);
        if constexpr (rosen)
          apply_rosenbrock_adjoint(st, sys, x0.data(), cp.t, cp.dt, n, n_phi,
                                   w_out.data(), adj, w_in.data(), m_wp.data(),
                                   m_rws, mu_seed.data());
        else
          onestep_recurse(st, x0.data(), cp.t, cp.dt, n, n_phi, w_out.data(),
                          adj, w_in.data(), m_wp.data(), m_ws); }
      for (std::size_t i = 0; i < n; ++i) w_in[i] += w_start[i];

      // lambda is what this step hands the one below it. eta is the cotangent
      // at the step end against the step's own embedded error, which a
      // one-step method reports already scaled.
      if (m_trace) {
        for (std::size_t i = 0; i < n; ++i) m_lam[k * n + i] = w_in[i];
        const std::vector<T>& xe = rosen ? m_rws.xerr : m_ws.xerr;
        T e = T(0.0);
        for (std::size_t i = 0; i < n; ++i) e += w_out[i] * xe[i];
        m_eta[k] = e;
      }
      w_out.swap(w_in);

      // The intervention this step was entered through. A one-step method
      // carries only the state, so there is no restart to collapse: what the
      // step hands down is already the cotangent on the state the jump left.
      if constexpr (Jumps::active) {
        const std::size_t ei = store.event_before(k);
        if (ei < store.n_events()) {
          const auto& e = store.event(ei);
          w_after = w_out;
          for (std::size_t o = 0; o < store.n_obs(); ++o)
            if (store.obs(o).event == ei)
              for (std::size_t i = 0; i < n; ++i)
                w_after[i] += seeds[o * n + i];

          w_before.assign(n, T(0.0));
          if (e.root)
            apply_root_jump_adjoint(e.x_before, e.x_after, e.t, jumps.root,
                                    e.triggered, jumps.sys, jumps.eadj, n,
                                    w_after.data(), w_before.data(),
                                    m_wp.data(), m_jws);
          else
            apply_fixed_jump_adjoint(e.x_before, e.t, jumps.fixed, jumps.root,
                                     e.switched, jumps.sys, jumps.eadj, n,
                                     w_after.data(), w_before.data(),
                                     m_wp.data(), m_jws);

          w_out.assign(n, T(0.0));
          if (e.after_step > 0) {
            pending = w_before;
            pending_t = e.t_before;
            pending_interp = true;
          } else {
            for (std::size_t i = 0; i < n; ++i) m_wx[i] += w_before[i];
          }
        }
      }
    }

    for (std::size_t i = 0; i < n; ++i) m_wx[i] += w_out[i];

    // Anything observed before the first step is the initial state itself.
    while (next_obs > 0) {
      --next_obs;
      const T* w = seeds + next_obs * n;
      for (std::size_t i = 0; i < n; ++i) m_wx[i] += w[i];
    }
  }

  void report_profile() const { m_prof.report("cppDE written adjoint"); }

  const std::vector<scalar_type>& wx0() const { return m_wx; }
  const std::vector<scalar_type>& whistory0() const { return m_whist; }
  const std::vector<scalar_type>& wp() const { return m_wp; }

  /// Under trace_lambda: [n_steps, n_states] step-major, and one per step.
  const std::vector<scalar_type>& lambda() const { return m_lam; }
  const std::vector<scalar_type>& eta() const { return m_eta; }

private:
  /// A cotangent `w` scaled by `c` on f(x, t), put onto x and theta.
  template<class AdjTerms>
  void seed_through_rhs(const AdjTerms& adj, const std::vector<scalar_type>& x,
                        double t, const scalar_type* w, double c, std::size_t n,
                        std::size_t n_phi, scalar_type* w_x,
                        std::vector<scalar_type>& jv)
  {
    if (c == 0.0) return;
    m_ws.m.assign(n, scalar_type(0.0));
    for (std::size_t i = 0; i < n; ++i) m_ws.m[i] = c * w[i];
    adj.jac_t_vec(x, m_ws.m, t, jv);
    for (std::size_t i = 0; i < n; ++i) w_x[i] += jv[i];
    adj.dfdp_t_vec_axpy(x, m_ws.m, t, 1.0, m_wp.data());
    (void)n_phi;
  }

  /// How many stage cotangents a step carries, when it carries any.
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
};

}  // namespace adjoint
}  // namespace cppde

#endif  // CPPDE_ADJOINT_STEP_HPP
