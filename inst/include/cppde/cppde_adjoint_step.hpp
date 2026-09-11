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
 * `jac_t_vec` for J' lambda and `dfdp_t_vec` for (df/dp)' lambda. Both size
 * their own output, so nothing here has to know the model's dimensions to hand
 * one a buffer.
 *
 * Copyright (C) 2026 Simon Beyer
 */

#ifndef CPPDE_ADJOINT_STEP_HPP
#define CPPDE_ADJOINT_STEP_HPP

#include <cstddef>
#include <vector>

#include <cppde/cppde_profiler.hpp>

namespace cppde {
namespace adjoint {

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
  void apply_transposed(const double* w, std::size_t n, double* out) const {
    for (int j = 0; j < cols; ++j)
      for (int k = 0; k < rows; ++k) {
        const double m = (*this)(k, j);
        if (m == 0.0) continue;
        const double* wk = w + static_cast<std::size_t>(k) * n;
        double* oj = out + static_cast<std::size_t>(j) * n;
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
  Stepper pre, tail_probe;

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
template<class Stepper>
void carry_into_pred(const multistep_operators<Stepper>& ops, std::size_t n,
                     const double* w_out, double* w_pred, double* w_acor)
{
  ops.B.apply_transposed(w_out, n, w_pred);
  for (int k = 0; k <= ops.q_out; ++k) {
    const double m = ops.c[static_cast<std::size_t>(k)];
    if (m == 0.0) continue;
    const double* wk = w_out + static_cast<std::size_t>(k) * n;
    for (std::size_t i = 0; i < n; ++i) w_acor[i] += m * wk[i];
  }
}

/// The buffers one step adjoint needs, kept across a sweep so a step allocates
/// nothing.
struct multistep_workspace {
  std::vector<double> w_pred, w_acor, mu, x, q;
};

/// The step adjoint proper, on a predicted-history cotangent that the caller
/// has already assembled.
template<class Stepper, class Solver, class AdjTerms>
void apply_multistep_adjoint_pre(const multistep_operators<Stepper>& ops,
                                 std::size_t n, std::size_t n_phi,
                                 const double* y, double t_new,
                                 Solver& solver, const AdjTerms& adj,
                                 double* w_in, double* w_theta,
                                 multistep_workspace& ws,
                                 double* w_out_state = nullptr)
{
  std::vector<double>& w_pred = ws.w_pred;
  std::vector<double>& w_acor = ws.w_acor;

  std::vector<double>& mu = ws.mu;
  mu = w_acor;
  for (std::size_t i = 0; i < n; ++i) w_pred[i] -= w_acor[i];

  if (w_out_state) for (std::size_t i = 0; i < n; ++i) w_out_state[i] = mu[i];

  solver.transposed(mu);

  for (std::size_t i = 0; i < n; ++i) w_pred[i] += mu[i];
  if (ops.q_in >= 1)
    for (std::size_t i = 0; i < n; ++i) w_pred[n + i] -= ops.rl1 * mu[i];

  {
    ws.x.assign(y, y + n);
    ws.q.assign(n_phi, 0.0);
    adj.dfdp_t_vec(ws.x, mu, t_new, ws.q);
    for (std::size_t k = 0; k < ws.q.size(); ++k) w_theta[k] += ops.gamma * ws.q[k];
  }

  ops.A.apply_transposed(w_pred.data(), n, w_in);
}

template<class Stepper, class Solver, class AdjTerms>
void apply_multistep_adjoint(const multistep_operators<Stepper>& ops,
                             std::size_t n, std::size_t n_phi,
                             const double* y, double t_new,
                             const double* w_out,
                             Solver& solver, const AdjTerms& adj,
                             double* w_in, double* w_theta,
                             multistep_workspace& ws,
                             double* w_out_state = nullptr)
{
  const std::size_t nz_in = static_cast<std::size_t>(ops.q_in + 1) * n;
  std::vector<double>& w_pred = ws.w_pred;
  std::vector<double>& w_acor = ws.w_acor;
  w_pred.assign(nz_in, 0.0);
  w_acor.assign(n, 0.0);
  carry_into_pred(ops, n, w_out, w_pred.data(), w_acor.data());
  apply_multistep_adjoint_pre(ops, n, n_phi, y, t_new, solver, adj,
                              w_in, w_theta, ws, w_out_state);
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
  /// Whether the sweep also keeps lambda and the refinement indicator per step.
  void trace_lambda(bool on) { m_trace = on; }

  /// One seed column. `seeds` is [n_obs x n_states] row-major.
  template<class Store, class AdjTerms, class Solver>
  void sweep(const Store& store, std::size_t n_phi, const double* seeds,
             const AdjTerms& adj, Solver& solver)
  {
    const std::size_t n = store.n_states();
    const std::size_t n_steps = store.n_steps();

    m_wp.assign(n_phi, 0.0);
    m_wx.assign(n, 0.0);
    m_whist.clear();
    m_lam.assign(m_trace ? n_steps * n : 0u, 0.0);
    m_eta.assign(m_trace ? n_steps : 0u, 0.0);
    if (n_steps == 0) return;

    // The cotangent the step above hands down, on its own carry.
    std::vector<double> w_carry;
    std::size_t next_obs = store.n_obs();

    multistep_operators<Stepper> ops;
    // A dense row is one entry per Nordsieck slot, a contraction one per
    // state. Two lengths, two buffers.
    std::vector<double> dense_row, w_in, jtv;

    for (std::size_t k = n_steps; k-- > 0;) {
      const auto& cp = store.step(k);

      // Does anything observe inside this step? Then the probe has to carry
      // this step's own interpolant, not the one it kept from another.
      std::size_t obs_lo = next_obs;
      while (obs_lo > 0 && store.obs(obs_lo - 1).step == k + 1) --obs_lo;
      const bool observed = (obs_lo < next_obs);
      if (observed) m_probe.valid = false;

      const double tail_key =
          cp.q_next + 1e3 * cp.eta + 1e6 * static_cast<double>(cp.ops.size());
      // Timed on the rebuild alone, so the call count is the miss count.
      if (!m_probe.matches(cp.carry, cp.dt, tail_key)) {
        auto _tp = m_prof.timer(cppde::prof_cat::rev_operators);
        m_probe.rebuild(cp.carry, cp.dt,
                        [&](Stepper& pr) { cp.apply_tail(pr, m_null); }, ops);
      }

      const std::size_t nz_in = static_cast<std::size_t>(ops.q_in + 1) * n;
      const std::size_t nz_out = static_cast<std::size_t>(ops.q_out + 1) * n;
      m_ws.w_pred.assign(nz_in, 0.0);
      m_ws.w_acor.assign(n, 0.0);

      if (!w_carry.empty()) {
        auto _tp = m_prof.timer(cppde::prof_cat::rev_adjoint);
        w_carry.resize(nz_out, 0.0);
        carry_into_pred(ops, n, w_carry.data(), m_ws.w_pred.data(),
                        m_ws.w_acor.data());
      }

      // The observations this step carries, through its own interpolant.
      for (std::size_t o = obs_lo; o < next_obs; ++o) {
        auto _tp = m_prof.timer(cppde::prof_cat::rev_interp);
        // The probe's states are the slots plus one for acor, so the row it
        // writes is d x_interp / d (zn_pred[0..q], acor).
        dense_row.assign(static_cast<std::size_t>(ops.q_in + 2), 0.0);
        m_probe.tail_probe.eval_dense_into(store.obs(o).t, dense_row);
        const double* w = seeds + o * n;
        for (int j = 0; j <= ops.q_in; ++j)
          for (std::size_t i = 0; i < n; ++i)
            m_ws.w_pred[static_cast<std::size_t>(j) * n + i] +=
                dense_row[static_cast<std::size_t>(j)] * w[i];
        const std::size_t acor_slot = dense_row.size() - 1;
        for (std::size_t i = 0; i < n; ++i)
          m_ws.w_acor[i] += dense_row[acor_slot] * w[i];
      }
      next_obs = obs_lo;

      const double t_new = cp.t + ops.h;
      solver.prepare(cp.y, t_new, 1.0 / ops.gamma, ops.gamma);

      w_in.assign(nz_in, 0.0);
      if (m_trace) m_wout.assign(n, 0.0);
      { auto _tp = m_prof.timer(cppde::prof_cat::rev_adjoint);
        apply_multistep_adjoint_pre(ops, n, n_phi, cp.y.data(), t_new, solver,
                                    adj, w_in.data(), m_wp.data(), m_ws,
                                    m_trace ? m_wout.data() : nullptr); }

      // lambda is what this step hands the one below it, on the state slot.
      // eta is that cotangent against the step's own error estimate, which for
      // a corrector method is acor scaled by the order's error constant.
      if (m_trace) {
        for (std::size_t i = 0; i < n; ++i) m_lam[k * n + i] = w_in[i];
        double e = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
          double pred0 = 0.0;
          for (int j = 0; j <= ops.q_in; ++j)
            pred0 += ops.A(0, j) * cp.zn[static_cast<std::size_t>(j) * n + i];
          e += m_wout[i] * (cp.y[i] - pred0);
        }
        m_eta[k] = e * ops.err_scale;
      }
      w_carry.swap(w_in);
    }

    // The trajectory start. initialize() builds the whole history out of one
    // state, zn[0] = x0 and zn[1] = h f(x0, t0) with the rest zero, so the
    // history's cotangent collapses onto that state and nothing is left above
    // it. The same map an event restart applies, which is why the boundary is
    // not a special case of the step but its own.
    const auto& cp0 = store.step(0);
    for (std::size_t i = 0; i < n && i < w_carry.size(); ++i) m_wx[i] = w_carry[i];
    if (w_carry.size() >= 2 * n) {
      const double h0 = static_cast<double>(cp0.carry.h);
      m_ws.x.assign(cp0.start_state(), cp0.start_state() + n);
      m_ws.mu.assign(w_carry.begin() + static_cast<std::ptrdiff_t>(n),
                     w_carry.begin() + static_cast<std::ptrdiff_t>(2 * n));
      jtv.assign(n, 0.0);
      adj.jac_t_vec(m_ws.x, m_ws.mu, cp0.t, jtv);
      for (std::size_t i = 0; i < n; ++i) m_wx[i] += h0 * jtv[i];
      m_ws.q.assign(n_phi, 0.0);
      adj.dfdp_t_vec(m_ws.x, m_ws.mu, cp0.t, m_ws.q);
      for (std::size_t k = 0; k < n_phi; ++k) m_wp[k] += h0 * m_ws.q[k];
    }

    // Anything observed before the first step is the initial state itself.
    while (next_obs > 0) {
      --next_obs;
      const double* w = seeds + next_obs * n;
      for (std::size_t i = 0; i < n; ++i) m_wx[i] += w[i];
    }
  }

  /// Per-category timings of the sweep, to stderr. Compiled away without
  /// CPPDE_PROFILE. The transposed algebra reports itself, from the solver.
  void report_profile() const { m_prof.report("cppDE written adjoint"); }

  const std::vector<double>& wx0() const { return m_wx; }
  const std::vector<double>& whistory0() const { return m_whist; }
  const std::vector<double>& wp() const { return m_wp; }

  /// Under trace_lambda: [n_steps, n_states] step-major, and one per step.
  const std::vector<double>& lambda() const { return m_lam; }
  const std::vector<double>& eta() const { return m_eta; }

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
  multistep_probe<Stepper> m_probe;
  multistep_workspace m_ws;
  bool m_trace = false;
  std::vector<double> m_wx, m_whist, m_wp, m_lam, m_eta, m_wout;
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
struct onestep_workspace {
  std::vector<double> xout, xerr, m, u, x_stage, q;
  std::vector<std::vector<double> > mm;
};

template<class Stepper, class System, class AdjTerms>
void apply_onestep_adjoint(Stepper& st, System& sys,
                           const double* x0, double t, double dt,
                           std::size_t n, std::size_t n_phi,
                           const double* w_out,
                           const AdjTerms& adj,
                           double* w_in, double* w_theta,
                           onestep_workspace& ws,
                           double* w_out_state = nullptr)
{
  constexpr int S = Stepper::n_stages_used;

  ws.xout.assign(n, 0.0);
  ws.xerr.assign(n, 0.0);
  std::vector<double> x(x0, x0 + n);
  st.do_step(sys, x, t, ws.xout, dt, ws.xerr);

  if (w_out_state) for (std::size_t i = 0; i < n; ++i) w_out_state[i] = w_out[i];

  if (ws.mm.size() < static_cast<std::size_t>(S) + 1) ws.mm.resize(S + 1);
  std::vector<std::vector<double> >& U = ws.mm;
  for (int i = 1; i <= S; ++i) U[i].assign(n, 0.0);

  for (std::size_t i = 0; i < n; ++i) w_in[i] = w_out[i];

  // Newest stage first: u_j is needed by every earlier stage and by nothing
  // later, so one pass suffices.
  for (int i = S; i >= 1; --i) {
    ws.m.assign(n, 0.0);
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

    ws.q.assign(n_phi, 0.0);
    adj.dfdp_t_vec(ws.x_stage, ws.m, ti, ws.q);
    for (std::size_t k = 0; k < ws.q.size(); ++k) w_theta[k] += ws.q[k];
  }
}

}  // namespace adjoint
}  // namespace cppde

#endif  // CPPDE_ADJOINT_STEP_HPP
