/**
 * @file cppde_adjoint_step.hpp
 * @brief The adjoint of one accepted step, written out per stepper family.
 *
 * The accepted grid is frozen, so a step is a fixed map with an adjoint that is
 * stated. What a multistep step does to the Nordsieck history splits in three:
 *
 *   zn_pred = A zn_in                      rescale and Pascal shift
 *   acor = N(zn_pred, theta)               the k chord-Newton iterations taken
 *   zn_out  = B zn_pred + c acor           the tail
 *
 * A, B and c are read off the stepper by running the stepper routines on a unit
 * slot. The model supplies `jac_t_vec` (J' lambda) and `dfdp_t_vec_axpy`
 * ((df/dp)' lambda, scaled and added). See "Adjoint of the step map" in
 * https://dmodverse.github.io/cppDE/articles/methods-reverse.html.
 *
 * Copyright (C) 2026 Simon Beyer
 */

#ifndef CPPDE_ADJOINT_STEP_HPP
#define CPPDE_ADJOINT_STEP_HPP

#include <algorithm>
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
#include <cppde/cppde_jump_adjoint.hpp>

// The hot loops here run over buffers that are always distinct. Saying so is
// what lets them vectorize.
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
  /// and which decides whether the inner loop vectorizes.
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

// ---------------------------------------------------------------------------
//  The adjoint of one step: w_pred = B' w_out, the corrector iterations
//  transposed in reverse order against the frozen iteration matrix, then
//  w_in = A' w_pred. See "Adjoint of the step map" in
//  https://dmodverse.github.io/cppDE/articles/methods-reverse.html.
// ---------------------------------------------------------------------------
/// What the step above hands down, passed back through the tail onto the
/// predicted history and onto acor. Separate from the step adjoint itself
/// because a trajectory adds the observations to the same two.
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
  std::vector<T> w_pred, w_acor, x, rho, abar, g, mu, rho_sum, g_sum;
  std::vector<double> wt, x_val;
};

/// The stopping rule of the converged corrector adjoint below: rtol 0 leaves
/// the iterations the forward run took.
struct corrector_control {
  double rtol = 0.0, atol = 0.0;
  int max_iter = 10;   ///< iterations per matrix before the next one is tried
  /// The share of the backward local error test the corrector may take, as
  /// CVODES lets the Newton iteration take a tenth of the forward test.
  double coef = 0.1;
};

/// Steps of the converged corrector adjoint, iterations, steps moved to the
/// matrix at the accepted state, and steps that missed the test there.
struct corrector_stats {
  std::size_t steps = 0, iterations = 0, fresh = 0, unconverged = 0;
};

/// The corrector's adjoint. The forward iteration from a_0 = 0 is
///   a_{m+1} = a_m + s M^-1 [f(zp + a_m) - (rl1 zp1 + a_m) / gamma],
/// acor = a_k, M frozen. With abar_k = w_acor, for m = k-1 .. 0:
///   rho = s M^-T abar_{m+1},  abar_m = abar_{m+1} + J(y_m)' rho - rho / gamma,
/// and J(y_m)' rho, -(rl1/gamma) rho and (df/dtheta)(y_m)' rho go to zp, zp1
/// and theta. `iter` holds y_0 .. y_{k-1}, slot-major.
template<class Stepper, class Solver, class AdjTerms, class T>
void multistep_corrector_adjoint(const multistep_operators<Stepper>& ops,
                                 std::size_t n, const T* iter, int k_iter,
                                 double s_iter, double t_new,
                                 Solver& solver, const AdjTerms& adj,
                                 T* w_theta, multistep_workspace<T>& ws)
{
  T* CPPDE_RESTRICT w_pred = ws.w_pred.data();
  ws.abar = ws.w_acor;
  const double inv_gamma = 1.0 / ops.gamma;
  const double c_zn1 = ops.rl1 * inv_gamma;
  for (int m = k_iter - 1; m >= 0; --m) {
    ws.rho = ws.abar;
    solver.transposed(ws.rho);
    if (s_iter != 1.0) for (T& v : ws.rho) v *= s_iter;
    ws.x.assign(iter + static_cast<std::size_t>(m) * n,
                iter + static_cast<std::size_t>(m + 1) * n);
    zero_armed(ws.g, n);
    adj.jac_t_vec(ws.x, ws.rho, T(t_new), ws.g);
    for (std::size_t i = 0; i < n; ++i) {
      ws.abar[i] += ws.g[i] - inv_gamma * ws.rho[i];
      w_pred[i] += ws.g[i];
    }
    if (ops.q_in >= 1)
      for (std::size_t i = 0; i < n; ++i) w_pred[n + i] -= c_zn1 * ws.rho[i];
    if (w_theta) adj.dfdp_t_vec_axpy(ws.x, ws.rho, T(t_new), T(1.0), w_theta);
  }
}

/// The corrector's adjoint with the corrector solved at the accepted state y:
/// rho solves (I/gamma - J(y))' rho = w_acor, and J(y)' rho, -(rl1/gamma) rho
/// and (df/dtheta)(y)' rho go to zp, zp1 and theta. rho sums the iteration above
/// with J at y, against the frozen matrix and, where abar misses `ctl` within
/// max_iter, against the matrix at y and gamma prepared under `key`.
template<class Stepper, class Solver, class AdjTerms, class T>
void multistep_corrector_adjoint_converged(const multistep_operators<Stepper>& ops,
                                           std::size_t n, const T* y, double s_iter,
                                           double t_new, Solver& solver,
                                           const AdjTerms& adj, T* w_theta,
                                           multistep_workspace<T>& ws,
                                           const corrector_control& ctl,
                                           std::size_t key, corrector_stats& stats)
{
  using ad_traits::scalar_value;
  T* CPPDE_RESTRICT w_pred = ws.w_pred.data();
  const double inv_gamma = 1.0 / ops.gamma;
  const double c_zn1 = ops.rl1 * inv_gamma;
  ws.abar = ws.w_acor;
  ws.x.assign(y, y + n);
  ws.wt.resize(n);
  for (std::size_t i = 0; i < n; ++i) {
    const double w = ctl.rtol * std::abs(static_cast<double>(scalar_value(w_pred[i]))) + ctl.atol;
    ws.wt[i] = w > 0.0 ? 1.0 / w : 0.0;
  }
  zero_armed(ws.rho_sum, n);
  zero_armed(ws.g_sum, n);
  auto norm = [&]() {
    double sq = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
      const double r = static_cast<double>(scalar_value(ws.abar[i])) * ws.wt[i];
      sq += r * r;
    }
    return n > 0 ? std::sqrt(sq / static_cast<double>(n)) : 0.0;
  };
  ++stats.steps;
  bool fresh = false;
  double s = s_iter, a = norm(), a_prev = 0.0;
  int it = 0;
  while (a > ctl.coef) {
    // A rate that cannot reach the test within the budget, or no budget left,
    // moves to the matrix at y; there the budget ends the iteration.
    const bool slow = it >= 2 && a * std::pow(std::min(a / a_prev, 1.0),
                                              ctl.max_iter - it) > ctl.coef;
    if (it >= ctl.max_iter || slow) {
      if (fresh) { ++stats.unconverged; break; }
      ws.x_val.resize(n);
      for (std::size_t i = 0; i < n; ++i)
        ws.x_val[i] = static_cast<double>(scalar_value(y[i]));
      solver.prepare_fresh(ws.x_val, t_new, inv_gamma, key);
      fresh = true; s = 1.0; it = 0;
      ++stats.fresh;
    }
    ws.rho = ws.abar;
    if (fresh) solver.fresh_transposed(ws.rho); else solver.transposed(ws.rho);
    if (s != 1.0) for (T& v : ws.rho) v *= s;
    zero_armed(ws.g, n);
    adj.jac_t_vec(ws.x, ws.rho, T(t_new), ws.g);
    for (std::size_t i = 0; i < n; ++i) {
      ws.abar[i] += ws.g[i] - inv_gamma * ws.rho[i];
      ws.rho_sum[i] += ws.rho[i];
      ws.g_sum[i] += ws.g[i];
    }
    a_prev = a;
    a = norm();
    ++it;
    ++stats.iterations;
  }
  for (std::size_t i = 0; i < n; ++i) w_pred[i] += ws.g_sum[i];
  if (ops.q_in >= 1)
    for (std::size_t i = 0; i < n; ++i) w_pred[n + i] -= c_zn1 * ws.rho_sum[i];
  if (w_theta) adj.dfdp_t_vec_axpy(ws.x, ws.rho_sum, T(t_new), T(1.0), w_theta);
}

/// One step adjoint from the cotangent on the Nordsieck carry the step hands on, for a
/// caller holding the operators, the checkpoint and the prepared matrix.
template<class Stepper, class CP, class Solver, class AdjTerms, class T>
void apply_multistep_adjoint(const multistep_operators<Stepper>& ops,
                             const CP& cp, std::size_t n, double t_new,
                             const T* w_out, Solver& solver, const AdjTerms& adj,
                             T* w_in, T* w_theta, multistep_workspace<T>& ws)
{
  const std::size_t nz_in = static_cast<std::size_t>(ops.q_in + 1) * n;
  ws.w_pred.assign(nz_in, T(0.0));
  ws.w_acor.assign(n, T(0.0));
  carry_into_pred(ops, n, w_out, ws.w_pred.data(), ws.w_acor.data());
  multistep_corrector_adjoint(ops, n, cp.iter.data(), cp.k_iter, cp.s_iter,
                              t_new, solver, adj, w_theta, ws);
  ops.A.apply_transposed(ws.w_pred.data(), n, w_in);
}

// ---------------------------------------------------------------------------
//  A whole multistep trajectory backwards, step adjoints and jump adjoints in
//  reverse order. An observation inside a step takes the row
//  d x_interp / d (zn_pred, acor) from the probe's interpolant, so an observed
//  step rebuilds the operators rather than taking the cached ones.
// ---------------------------------------------------------------------------
template<class Stepper>
class closed_multistep_trajectory {
public:
  // The sweep runs in the stepper's own scalar type. In plain double that is
  // the gradient; over a dual it is the gradient and the directional
  // derivatives, which is forward over reverse.
  using scalar_type = typename Stepper::value_type;

  /// Whether the sweep also keeps lambda and the refinement indicator per step.
  void trace_lambda(bool on) { m_trace = on; }

  /// Solve the corrector of every step to convergence, to a tenth of the
  /// backward local error test against rtol |lambda| + atol; rtol 0, the
  /// default, transposes the iterations the forward run took.
  void set_corrector(double rtol, double atol) {
    m_corr.rtol = rtol; m_corr.atol = atol;
  }
  /// What the converged corrector adjoint did since the trajectory was made.
  const corrector_stats& corrector() const { return m_corr_stats; }

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
    solver.forget();

    // The cotangent the step above hands down, on the Nordsieck carry of that step.
    // Across an intervention it is instead a cotangent on a point inside the step
    // below, because the restart threw the Nordsieck carry away.
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
      ms_step_adjoint(store, k, observed, m_obs.data(), m_obs.size(), w_carry, n, n_phi,
                      adj, solver, w_in, m_wp.data(), m_trace ? &eta_k : nullptr);

      // lambda is what this step hands the one below it, on the state slot.
      if (m_trace) {
        for (std::size_t i = 0; i < n; ++i) m_lam[k * n + i] = w_in[i];
        m_eta[k] = eta_k;
      }
      w_carry.swap(w_in);

      // The intervention this step was entered through, if any. Two maps in the
      // order the forward run applied them, so swept the other way round: the
      // restart, which collapses the Nordsieck carry onto the state the jump ended on,
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

  /// Every seed column at once, step by step from the last: each step's
  /// operators and factorization are built once and serve all columns, which
  /// keep their own Nordsieck cotangent, handover, observations and cotangents. `seeds[c]` is
  /// column c, [n_obs x n_states] row-major; afterwards wx0(c) and wp(c) hold
  /// the cotangents of that column. Not with trace_lambda, which takes sweep() per column.
  template<class Store, class AdjTerms, class Solver, class Jumps = no_jumps>
  void sweep_many(const Store& store, std::size_t n_phi, const std::vector<const scalar_type*>& seeds,
                  const AdjTerms& adj, Solver& solver, const Jumps& jumps = Jumps())
  {
    using T = scalar_type;
    const std::size_t n = store.n_states();
    const std::size_t n_steps = store.n_steps();
    const std::size_t S = seeds.size();
    m_cols.resize(S);
    for (column& c : m_cols) {
      zero_armed(c.wp, n_phi);
      c.wx.assign(n, T(0.0));
      c.w_carry.clear();
      c.handover = jump_handover<T>();
      c.next_obs = store.n_obs();
    }
    if (n_steps == 0) {
      for (std::size_t s = 0; s < S; ++s)
        sweep_without_steps(store, n, seeds[s], adj, jumps, m_cols[s].wx, m_cols[s].wp, m_jws);
      return;
    }
    std::vector<T> w_in, w_after, w_before, jtv;
    zero_armed(jtv, n);
    solver.forget();
    for (std::size_t k = n_steps; k-- > 0;) {
      const auto& cp = store.step(k);
      store.use_step_modes(k);
      // the observations inside the step are the same for every column
      bool observed = false;
      for (std::size_t s = 0; s < S; ++s)
        observed = step_seeds(store, k, seeds[s], n, m_cols[s].next_obs, m_cols[s].handover, m_cols[s].obs) || observed;
      ms_step_prepare(store, k, observed, solver);
      for (std::size_t s = 0; s < S; ++s) {
        column& c = m_cols[s];
        w_in.clear();
        ms_step_column(cp, c.obs.data(), c.obs.size(), c.w_carry, n, n_phi, adj, solver, w_in,
                       c.wp.data(), nullptr);
        c.w_carry.swap(w_in);
      }
      if constexpr (Jumps::active) {
        const std::size_t ei = store.event_before(k);
        if (ei >= store.n_events()) continue;
        const auto& e = store.event(ei);
        for (std::size_t s = 0; s < S; ++s) {
          column& c = m_cols[s];
          w_after.assign(n, T(0.0));
          m_ws.x.assign(e.x_after.begin(), e.x_after.end());
          collapse_restart(c.w_carry, n, m_ws.x, e.t,
                           e.restart ? e.dt_restart
                                     : static_cast<double>(cp.carry.h),
                           adj, w_after.data(), c.wp.data(), m_ws.mu, jtv);
          event_adjoint(store, ei, seeds[s], n, adj, jumps, w_after, w_before,
                        c.wp.data(), m_jws);
          c.w_carry.clear();
          c.handover.take(e, w_before, c.wx);
        }
      }
    }
    // The trajectory start, which is the same restart with no jump under it.
    const auto& cp0 = store.step(0);
    store.use_step_modes(0);
    for (std::size_t s = 0; s < S; ++s) {
      column& c = m_cols[s];
      if (!c.w_carry.empty()) {
        w_after.assign(n, T(0.0));
        m_ws.x.assign(cp0.start_state(), cp0.start_state() + n);
        collapse_restart(c.w_carry, n, m_ws.x, cp0.t,
                         static_cast<double>(cp0.carry.h), adj, w_after.data(),
                         c.wp.data(), m_ws.mu, jtv);
        for (std::size_t i = 0; i < n; ++i) c.wx[i] += w_after[i];
      }
      seed_initial_state(store, c.next_obs, seeds[s], n, c.wx);
    }
  }

  /// The cotangents of column s after sweep_many().
  const std::vector<scalar_type>& wx0(std::size_t s) const { return m_cols[s].wx; }
  const std::vector<scalar_type>& wp(std::size_t s) const { return m_cols[s].wp; }

  /// The sweep that takes the flow between the forward run's grid points
  /// instead of the multistep scheme: each step's interval from the state the
  /// step started at, swept by `flow`. Only the state's cotangent passes from
  /// one step to the next, as the exact flow depends on nothing else.
  template<class Store, class System, class AdjTerms, class Flow, class Jumps = no_jumps>
  void sweep_flow(const Store& store, std::size_t n_phi, const scalar_type* seeds,
                  System& sys, const AdjTerms& adj, Flow& flow,
                  const Jumps& jumps = Jumps())
  {
    sweep_flow_many(store, n_phi, {seeds}, sys, adj, flow, jumps);
    m_wx = m_cols[0].wx;
    m_wp = m_cols[0].wp;
  }

  /// sweep_flow() for every seed column at once, the chains and factorizations
  /// of each interval shared; afterwards wx0(c) and wp(c) hold the cotangents
  /// of column c. Under trace_lambda lambda() is that of column 0.
  template<class Store, class System, class AdjTerms, class Flow, class Jumps = no_jumps>
  void sweep_flow_many(const Store& store, std::size_t n_phi,
                       const std::vector<const scalar_type*>& seeds,
                       System& sys, const AdjTerms& adj, Flow& flow,
                       const Jumps& jumps = Jumps())
  {
    using T = scalar_type;
    const std::size_t n = store.n_states();
    const std::size_t n_steps = store.n_steps();
    const std::size_t S = seeds.size();
    m_cols.resize(S);
    for (column& c : m_cols) {
      zero_armed(c.wp, n_phi);
      c.wx.assign(n, T(0.0));
      c.w_state.assign(n, T(0.0));
      c.handover = jump_handover<T>();
      c.next_obs = store.n_obs();
    }
    m_whist.clear();
    m_lam.assign(m_trace ? n_steps * n : 0u, T(0.0));
    // The flow has no error estimate of the multistep scheme to weigh.
    m_eta.assign(m_trace ? n_steps : 0u, T(std::numeric_limits<double>::quiet_NaN()));
    m_ref_m.assign(m_trace ? n_steps : 0u, 1);
    flow.begin_flow(n_phi, S);
    if (n_steps == 0) {
      for (std::size_t s = 0; s < S; ++s)
        sweep_without_steps(store, n, seeds[s], adj, jumps, m_cols[s].wx, m_cols[s].wp, m_jws);
      return;
    }
    // Per column, the observations in the step and the cotangent on the state at the grid point reached.
    std::vector<const std::vector<obs_ref<T>>*> obs(S);
    std::vector<std::vector<T>*> w(S);
    for (std::size_t s = 0; s < S; ++s) { obs[s] = &m_cols[s].obs; w[s] = &m_cols[s].w_state; }
    std::vector<T> x0, w_after, w_before;
    for (std::size_t k = n_steps; k-- > 0;) {
      const auto& cp = store.step(k);
      store.use_step_modes(k);
      for (std::size_t s = 0; s < S; ++s)
        step_seeds(store, k, seeds[s], n, m_cols[s].next_obs, m_cols[s].handover, m_cols[s].obs);
      x0.assign(cp.start_state(), cp.start_state() + n);
      // The end state the run computed, where a step follows without an event.
      const double* xh_end =
          (k + 1 < n_steps && store.event_before(k + 1) >= store.n_events())
          ? store.step(k + 1).start_state() : nullptr;
      const int mk = flow.flow_interval(sys, adj, n, n_phi, x0, cp.t, cp.dt, obs, w, xh_end);
      if (m_trace) {
        for (std::size_t i = 0; i < n; ++i) m_lam[k * n + i] = m_cols[0].w_state[i];
        m_ref_m[k] = mk;
      }

      // The intervention this step was entered through, on the state alone.
      if constexpr (Jumps::active) {
        const std::size_t ei = store.event_before(k);
        if (ei >= store.n_events()) continue;
        for (std::size_t s = 0; s < S; ++s) {
          column& c = m_cols[s];
          w_after = c.w_state;
          event_adjoint(store, ei, seeds[s], n, adj, jumps, w_after, w_before,
                        c.wp.data(), m_jws);
          c.w_state.assign(n, T(0.0));
          c.handover.take(store.event(ei), w_before, c.wx);
        }
      }
    }
    for (std::size_t s = 0; s < S; ++s) {
      column& c = m_cols[s];
      for (std::size_t i = 0; i < n; ++i) c.wx[i] += c.w_state[i];
      seed_initial_state(store, c.next_obs, seeds[s], n, c.wx);
      const auto& fp = flow.wp(s);
      for (std::size_t i = 0; i < n_phi; ++i) c.wp[i] += fp[i];
    }
  }

  /// Under trace_lambda: the substeps each step's interval was swept in.
  const std::vector<int>& refine_m() const { return m_ref_m; }

private:
  using T_ = scalar_type;

  /// The adjoint of one step from the step checkpoint: `w_carry` the cotangent on
  /// the Nordsieck carry the step hands on, empty where an intervention handed a
  /// point inside the step instead; `w_in` receives the one on the Nordsieck
  /// carry the step was entered with, `wp` the step's share of the gradient, `eta` the cotangent
  /// on acor against the step's error estimate.
  template<class Store, class AdjTerms, class Solver>
  void ms_step_adjoint(const Store& store, std::size_t k, bool observed,
                       const obs_ref<T_>* ob, std::size_t nob,
                       std::vector<T_>& w_carry, std::size_t n, std::size_t n_phi,
                       const AdjTerms& adj, Solver& solver, std::vector<T_>& w_in,
                       T_* wp, T_* eta)
  {
    const auto& cp = store.step(k);
    ms_step_prepare(store, k, observed || nob > 0, solver);
    ms_step_column(cp, ob, nob, w_carry, n, n_phi, adj, solver, w_in, wp, eta);
  }

  /// What a step's adjoint shares between seed columns: the step operators, and
  /// the factorized iteration matrix the corrector ran against, built under the
  /// switch modes it was evaluated under.
  template<class Store, class Solver>
  void ms_step_prepare(const Store& store, std::size_t k, bool observed, Solver& solver)
  {
    const auto& cp = store.step(k);
    multistep_operators<Stepper>& ops = m_ops;
    m_step = k;
    // An observed step needs a fresh interpolant, not one kept from another.
    if (observed) m_probe.valid = false;
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
    if (cp.k_iter == 0) return;
    const auto& M = store.matrix(cp.mat);
    store.use_step_modes(M.modes_step);
    solver.prepare_frozen(M.x, M.t, M.inv_gamma_p, cp.mat);
    store.use_step_modes(k);
  }

  /// The step adjoint of one seed column on the operators and factorization
  /// ms_step_prepare() built for the step.
  template<class CP, class AdjTerms, class Solver>
  void ms_step_column(const CP& cp, const obs_ref<T_>* ob, std::size_t nob,
                      std::vector<T_>& w_carry, std::size_t n, std::size_t n_phi,
                      const AdjTerms& adj, Solver& solver, std::vector<T_>& w_in,
                      T_* wp, T_* eta)
  {
    using T = T_;
    multistep_operators<Stepper>& ops = m_ops;
    const std::size_t nz_in = static_cast<std::size_t>(ops.q_in + 1) * n;
    const std::size_t nz_out = static_cast<std::size_t>(ops.q_out + 1) * n;
    m_ws.w_pred.assign(nz_in, T(0.0));
    m_ws.w_acor.assign(n, T(0.0));

    // A cotangent at a time inside the step, through the interpolant of the step. The
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

    // What the step above handed back, on the Nordsieck carry.
    if (!w_carry.empty()) {
      auto _tp = m_prof.timer(cppde::prof_cat::rev_adjoint);
      w_carry.resize(nz_out, T(0.0));
      carry_into_pred(ops, n, w_carry.data(), m_ws.w_pred.data(),
                      m_ws.w_acor.data());
    }

    const double t_new = cp.t + ops.h;

    // eta is the cotangent on acor against the step's own error estimate,
    // which for a corrector method is acor scaled by the order's constant;
    // the first iterate is the prediction.
    if (eta) {
      T e = T(0.0);
      for (std::size_t i = 0; i < n; ++i)
        e += m_ws.w_acor[i] * (cp.y[i] - cp.iter[i]);
      *eta = e * ops.err_scale;
    }

    w_in.assign(nz_in, T(0.0));
    { auto _tp = m_prof.timer(cppde::prof_cat::rev_adjoint);
      if (m_corr.rtol > 0.0 && cp.k_iter > 0)
        multistep_corrector_adjoint_converged(ops, n, cp.y.data(), cp.s_iter, t_new,
                                              solver, adj, wp, m_ws, m_corr, m_step,
                                              m_corr_stats);
      else
        multistep_corrector_adjoint(ops, n, cp.iter.data(), cp.k_iter, cp.s_iter,
                                    t_new, solver, adj, wp, m_ws);
      ops.A.apply_transposed(m_ws.w_pred.data(), n, w_in.data()); }
    (void)n_phi;
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
  corrector_control m_corr;
  corrector_stats m_corr_stats;
  std::size_t m_step = 0;
  std::vector<scalar_type> m_wx, m_whist, m_wp, m_lam, m_eta;
  // per seed column of sweep_many()
  struct column {
    // w_state: the cotangent on the state at the grid point, in sweep_flow_many()
    std::vector<scalar_type> wx, wp, w_carry, w_state;
    jump_handover<scalar_type> handover;
    std::size_t next_obs = 0;
    std::vector<obs_ref<scalar_type>> obs;
  };
  std::vector<column> m_cols;
};

// ---------------------------------------------------------------------------
//  The adjoint of one explicit Runge-Kutta step, stages recomputed by one
//  forward run: m_i = h b_i w + h sum_{j>i} a_ji u_j with u_j = J(X_j)' m_j,
//  w_x = w + sum_i u_i. See "Adjoint of the step map" in
//  https://dmodverse.github.io/cppDE/articles/methods-reverse.html.
// ---------------------------------------------------------------------------
template<class T = double>
struct onestep_workspace {
  std::vector<T> xout, xerr, m, u, x_stage, q;
  std::vector<std::vector<T> > mm;
};

/// The backward recursion alone, on a stepper that has just run the step. A
/// trajectory runs the step itself, because the trajectory interpolant reads the stages
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
//  d(J g). See "Adjoint of the step map" in
//  https://dmodverse.github.io/cppDE/articles/methods-reverse.html.
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

  // The step once forward, to recover the stages and the factorization. An
  // accepted step rebuilds the Jacobian at the step start, so this is the same
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
  // so it takes the last stage's combination and not a separate row.
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
  // the embedded order, on the same stages and factorization.
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
// than explicit combinations the tableau already describes. That is what decides
// how the adjoint moves through a step.
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
  /// Under refine: the first-order change of the seeded functional from the
  /// computed trajectory to the one the last sweep differentiated, the sum of
  /// the seeds and cotangents times the defects of the substepped maps.
  double correction() const { return m_corr; }

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
    m_corr = 0.0;
    if (n_steps == 0) {
      sweep_without_steps(store, n, seeds, adj, jumps, m_wx, m_wp, m_jws);
      return;
    }

    std::size_t next_obs = store.n_obs();
    std::vector<T> w_out(n, T(0.0)), w_in(n, T(0.0)), x0, w_after, w_before;
    jump_handover<T> handover;
    std::vector<double> xh_end;

    for (std::size_t k = n_steps; k-- > 0;) {
      const auto& cp = store.step(k);
      store.use_step_modes(k);
      x0.assign(cp.start_state(), cp.start_state() + n);
      step_seeds(store, k, seeds, n, next_obs, handover, m_obs);

      T eta_k = T(0.0);
      if (m_refine_rtol > 0.0) {
        // The end state the run computed; after an event the cotangent at the
        // step end is zero and the state the jump started from is not this one.
        const bool has_end = k + 1 < n_steps && store.event_before(k + 1) >= store.n_events();
        xh_end.clear();
        if (has_end) {
          const T* xe = store.step(k + 1).start_state();
          for (std::size_t i = 0; i < n; ++i)
            xh_end.push_back(static_cast<double>(ad_traits::scalar_value(xe[i])));
        }
        refined_step(sys, adj, st, n, n_phi, x0, cp.t, cp.dt, w_out, w_in, eta_k, k,
                     has_end ? xh_end.data() : nullptr);
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

  /// The adjoint of the step from `xs` at `t0` over `dt`, seeds included:
  /// `w_end` is the cotangent at the step's end, `w_in` receives the one at
  /// the step start, `wp` the step's share of the gradient. `eta`, if given, gets
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
    // An explicit interpolant reads f at both ends of the step, so the stages
    // must exist before an observation is seeded; a Rosenbrock interpolant is
    // linear in the stage vectors and needs only the weights.
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

  /// The step's adjoint, checked per component of the step's share of the gradient:
  /// by the embedded estimate where the method has one, then by halving, the
  /// step taken as m substeps of the same method from the step start, m
  /// doubled until the share agrees with that of m / 2 after Richardson
  /// extrapolation at the method's order. The finest result is kept.
  template<class System, class AdjTerms>
  void refined_step(System& sys, const AdjTerms& adj, Stepper& st, std::size_t n,
                    std::size_t n_phi, const std::vector<scalar_type>& x0, double t0,
                    double dt, const std::vector<scalar_type>& w_end,
                    std::vector<scalar_type>& w_in, scalar_type& eta, std::size_t k,
                    const double* xh_end)
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
    double corr = 0.0;
    int m = 1;
    const double dir = dt > 0 ? 1.0 : -1.0;
    // The seeds of substep j of m: those before the step's start into the
    // first, those past the step end into the last, in the direction of integration.
    auto in_substep = [&](const obs_ref& o, int j, int m_, double a, double b) {
      return (j == 0 || dir * (o.t - a) > 0) && (j == m_ - 1 || dir * (o.t - b) <= 0);
    };
    while (m < m_refine_max) {
      m *= 2;
      const double h = dt / m;
      // The substeps forward, from the step's own start, and the value change
      // they make: the cotangents times the defects at the step end and at the
      // observations, read off each substep's dense output.
      double c = 0.0;
      m_xs.resize(static_cast<std::size_t>(m + 1));
      m_xs[0] = x0;
      for (int j = 0; j < m; ++j) {
        zero_armed(m_xs[static_cast<std::size_t>(j + 1)], n);
        zero_armed(m_xe, n);
        st.do_step(sys, m_xs[static_cast<std::size_t>(j)], t0 + j * h,
                   m_xs[static_cast<std::size_t>(j + 1)], h, m_xe);
        const double a = t0 + j * h, b = t0 + (j + 1) * h;
        for (const obs_ref& o : m_obs) {
          if (!o.xh || !in_substep(o, j, m, a, b)) continue;
          dense_value(sys, st, m_xs[static_cast<std::size_t>(j)],
                      m_xs[static_cast<std::size_t>(j + 1)], a, h, o.t, n);
          for (std::size_t i = 0; i < n; ++i)
            c += static_cast<double>(scalar_value(o.w[i]))
               * (static_cast<double>(scalar_value(m_xd[i])) - o.xh[i]);
        }
      }
      if (xh_end)
        for (std::size_t i = 0; i < n; ++i)
          c += static_cast<double>(scalar_value(w_end[i]))
             * (static_cast<double>(scalar_value(m_xs[static_cast<std::size_t>(m)][i])) - xh_end[i]);
      // And back, each with the seeds that fall into it.
      zero_armed(m_pf, n_phi);
      m_wa = w_end;
      for (int j = m - 1; j >= 0; --j) {
        const double a = t0 + j * h, b = t0 + (j + 1) * h;
        m_sub.clear();
        for (const obs_ref& o : m_obs)
          if (in_substep(o, j, m, a, b)) m_sub.push_back(o);
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
      corr = c;
      if (est <= 1.0 || stalled) break;
    }
    m_corr += corr;
    for (std::size_t i = 0; i < n_phi; ++i) m_wp[i] += m_pc[i];
    if (est > 1.0) ++m_ref_fail;
    if (m_trace && k < m_ref_m.size()) m_ref_m[k] = m;
  }

  /// The dense output of the substep from `xa` at `a` over `h` to `xb`, at t
  /// clamped to the substep, into m_xd: the interpolant the adjoint seeds
  /// through. Rosenbrock reads the stages of the substep just taken; an
  /// explicit method uses f at both ends, so no stage of `st` is touched.
  template<class System>
  void dense_value(System& sys, Stepper& st, const std::vector<scalar_type>& xa,
                   const std::vector<scalar_type>& xb, double a, double h, double t,
                   std::size_t n)
  {
    const double lo = h > 0 ? a : a + h, hi = h > 0 ? a + h : a;
    const double tc = t < lo ? lo : (t > hi ? hi : t);
    zero_armed(m_xd, n);
    if constexpr (has_stage_vectors<Stepper>::value) {
      st.prepare_dense_output();
      st.calc_state(tc, m_xd, xa, a, xb, a + h);
    } else {
      double cw[4];
      Stepper::dense_weights((tc - a) / h, cw);
      zero_armed(m_fa, n);
      zero_armed(m_fb, n);
      sys.first(xa, m_fa, a);
      sys.first(xb, m_fb, a + h);
      for (std::size_t i = 0; i < n; ++i)
        m_xd[i] = cw[0] * xa[i] + cw[1] * xb[i] + (h * cw[2]) * m_fa[i]
                + (h * cw[3]) * m_fb[i];
    }
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
                           m_pc, m_pf, m_wa, m_wb, m_pe, m_we, m_xd, m_fa, m_fb;
  double m_corr = 0.0;
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
