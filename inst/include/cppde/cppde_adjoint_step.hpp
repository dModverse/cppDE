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
 * `jac_t_vec` for J' lambda and `dfdp_t_vec` for (df/dp)' lambda.
 *
 * Copyright (C) 2026 Simon Beyer
 */

#ifndef CPPDE_ADJOINT_STEP_HPP
#define CPPDE_ADJOINT_STEP_HPP

#include <cstddef>
#include <vector>

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
               double& rl1, double& gamma, double& h)
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
  struct key {
    int q = -1, L = 0, qwait = 0;
    bool started = false;
    double dt = 0.0, h = 0.0, hscale = 0.0, eta = 0.0;
    double tail = 0.0;
    std::vector<double> tau;
    bool operator==(const key& o) const {
      return q == o.q && L == o.L && qwait == o.qwait && started == o.started &&
             dt == o.dt && h == o.h && hscale == o.hscale && eta == o.eta &&
             tail == o.tail && tau == o.tau;
    }
  };
  key last;
  bool valid = false;

  /// Build A, B and c for one step. `tail` applies whatever the run recorded
  /// between this step's acceptance and the next one's; `tail_key` stands for
  /// what it will do, so a repeat is recognised without running it.
  template<class Carry, class TailFn>
  void build(const Carry& carry, double dt, double tail_key, TailFn&& tail,
             multistep_operators<Stepper>& ops)
  {
    key k;
    k.q = carry.q; k.L = carry.L; k.qwait = carry.qwait;
    k.started = carry.nst > 0;
    k.dt = dt; k.h = carry.h; k.hscale = carry.hscale; k.eta = carry.eta;
    k.tail = tail_key;
    k.tau.assign(carry.tau.begin(), carry.tau.begin() + (carry.q + 1));
    if (valid && k == last) return;

    ops.q_in = carry.q;
    probe_pre(pre, carry, dt, ops.A, ops.rl1, ops.gamma, ops.h);
    probe_tail(tail_probe, carry, dt, tail, ops.B, ops.c, ops.q_out);
    last = k;
    valid = true;
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
/// The buffers one step adjoint needs, kept across a sweep so a step allocates
/// nothing.
struct multistep_workspace {
  std::vector<double> w_pred, w_acor, mu, x, q;
};

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

  ops.B.apply_transposed(w_out, n, w_pred.data());
  for (int k = 0; k <= ops.q_out; ++k) {
    const double m = ops.c[static_cast<std::size_t>(k)];
    if (m == 0.0) continue;
    const double* wk = w_out + static_cast<std::size_t>(k) * n;
    for (std::size_t i = 0; i < n; ++i) w_acor[i] += m * wk[i];
  }

  // acor = y - zn_pred[0]
  std::vector<double>& mu = ws.mu;
  mu = w_acor;
  for (std::size_t i = 0; i < n; ++i) w_pred[i] -= w_acor[i];

  // The step end's own cotangent, before the equation drives y's to zero. The
  // refinement indicator reads it, as wout() does on the tape path.
  if (w_out_state) for (std::size_t i = 0; i < n; ++i) w_out_state[i] = mu[i];

  solver.transposed(mu);

  for (std::size_t i = 0; i < n; ++i) w_pred[i] += mu[i];
  if (ops.q_in >= 1)
    for (std::size_t i = 0; i < n; ++i) w_pred[n + i] -= ops.rl1 * mu[i];

  // (df/dp)' mu, scaled by gamma, onto the flat parameter cotangent.
  {
    ws.x.assign(y, y + n);
    ws.q.assign(n_phi, 0.0);
    adj.dfdp_t_vec(ws.x, mu, t_new, ws.q);
    for (std::size_t k = 0; k < ws.q.size(); ++k) w_theta[k] += ops.gamma * ws.q[k];
  }

  ops.A.apply_transposed(w_pred.data(), n, w_in);
}

}  // namespace adjoint
}  // namespace cppde

#endif  // CPPDE_ADJOINT_STEP_HPP
