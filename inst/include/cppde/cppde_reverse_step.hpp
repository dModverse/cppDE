/*
 One integration step backwards: checkpoint, replay under codual, sweep.

 A step is the map (x, t, h, theta) -> (x_out, xerr), where x, t and h are what
 the step reads that an earlier step produced. The forward run stores them per
 accepted step; the reverse run replays the step under cppde::codual, and one
 sweep turns a cotangent on x_out into cotangents on all four. Tape and replay
 live one step at a time, so the memory bound is a step, not a trajectory.

 h is a tape independent, not a constant. The controller derives it from the
 previous step's error estimate, so the adjoint runs through the step-size and
 order control rather than around it; xerr is recorded for the same reason and
 is what the trajectory sweep seeds the control law through.

 Replay rather than hand-written adjoint equations: the sweep must differentiate
 the stepper's own arithmetic, which no codegen emits.

 The tape follows the branch the forward run took. Acceptance, order, iteration
 and rebuild counts are piecewise constant and stay control decisions; what is
 smooth within one control path is differentiated.

 Copyright (C) 2026 Simon Beyer
 */

#ifndef CPPDE_REVERSE_STEP_HPP
#define CPPDE_REVERSE_STEP_HPP

#include <cstddef>
#include <type_traits>
#include <vector>

#include <cppde/cppde_ad_traits.hpp>
#include <cppde/cppde_codual.hpp>
#include <cppde/cppde_codual_math.hpp>
#include <cppde/cppde_lu.hpp>
#include <cppde/cppde_multistepper.hpp>
#include <cppde/cppde_rosenbrock4.hpp>
#include <cppde/cppde_tsit5.hpp>

namespace cppde {
namespace reverse {

// ============================================================================
//  step_checkpoint<Stepper, T>
//
//  Declared, never defined: every stepper family supplies its own. The contract,
//  in the order step_recorder calls it:
//
//    std::size_t n() const     states in the step-start state
//    double t, dt              the two the controller chose; the recorder makes
//                              them independents, they are not stepper-specific
//
//    capture(const Stepper&, const std::vector<Value>& x, double t, double dt)
//        Reads the live forward stepper. Values only, whatever it integrated in.
//
//    load(RStepper&, std::vector<codual<T>>& x,
//                    std::vector<codual<T>>& history) const
//        Writes the checkpoint into a replay stepper of the reverse value type
//        and registers everything crossing the step boundary as a tape
//        independent: the step-start state in x, stepper-internal carry in
//        history. Their adjoints after the sweep are what the previous step
//        consumes.
//
//    finish(RStepper&, const std::vector<codual<T>>& xout,
//                       std::vector<codual<T>>& carry_out)
//        The step's tail: whatever the controller does to the carry after
//        accepting, replayed from the decisions the forward run took, and the
//        carry the next step reads written into carry_out. Optional; without it
//        the carry out is the step end, which is a one-step method's whole
//        carry.
//
//  history is what carries the multistepper, whose carry is the Nordsieck array
//  and not the state alone. Empty for a one-step method.
// ============================================================================

template<class Stepper, class T> struct step_checkpoint;

// ----------------------------------------------------------------------------
//  tsit5: an explicit one-step method carries nothing across a step boundary.
//
//  FSAL is an optimisation, not a dependence. The recycled k1 is f(x, t) at the
//  checkpointed x, so the replay recomputes it bit for bit. Storing k7 would save
//  one right-hand-side call and cost the tape the dependence of k1 on x.
// ----------------------------------------------------------------------------

template<class Stepper, class T>
struct onestep_checkpoint {
  using scalar_type  = T;
  using stepper_type = Stepper;

  std::vector<T> x;         // step-start state
  double         t  = 0.0;
  double         dt = 0.0;

  std::size_t n() const { return x.size(); }

  template<class Value>
  void capture(const stepper_type& /*st*/, const std::vector<Value>& x_in,
               double t_in, double dt_in)
  {
    x.resize(x_in.size());
    for (std::size_t i = 0; i < x_in.size(); ++i)
      x[i] = static_cast<T>(ad_traits::scalar_value(x_in[i]));
    t  = t_in;
    dt = dt_in;
  }

  template<class RStepper>
  void load(RStepper& /*rst*/,
            std::vector<codual<T>>& x_rev,
            std::vector<codual<T>>& history) const
  {
    x_rev.assign(x.size(), codual<T>());
    for (std::size_t i = 0; i < x.size(); ++i) {
      x_rev[i] = codual<T>(x[i]);
      x_rev[i].independent();
    }
    history.clear();
  }
};

template<class Value, class Resizer, class T>
struct step_checkpoint<cppde::tsit5<Value, Resizer>, T>
  : onestep_checkpoint<cppde::tsit5<Value, Resizer>, T> {};

// ----------------------------------------------------------------------------
//  rosenbrock4: a one-step method too, so the same carry. Its stages are linear
//  solves rather than explicit combinations, and the replay recovers each stage
//  value from the same factorisation it sweeps with, so none of them is stored.
// ----------------------------------------------------------------------------

template<class Value, class Resizer, class T>
struct step_checkpoint<cppde::rosenbrock4<Value, Resizer>, T>
  : onestep_checkpoint<cppde::rosenbrock4<Value, Resizer>, T>
{
  // Its interpolant is built from the stages rather than read off them, so the
  // snapshot has to be taken before anything is interpolated. tsit5 needs no
  // such call, and must not get one: it would arm the FSAL recycle, and the next
  // replayed step starts from its own checkpoint rather than from this one.
  template<class RStepper>
  void finish(RStepper& rst, const std::vector<codual<T>>& xout,
              std::vector<codual<T>>& carry_out) const
  {
    rst.prepare_dense_output();
    carry_out = xout;
  }
};

// ----------------------------------------------------------------------------
//  multistepper: the carry is the Nordsieck history, not the state alone.
//
//  What the forward run has to leave behind is that history, the scalars that
//  pin its meaning, and the state the corrector converged to, which the replay
//  puts back rather than iterating for. The tail is the controller's: the rank-1
//  Nordsieck update, the order it chose and the rescale it applied, all replayed
//  from the decisions rather than recomputed, since those are control decisions.
// ----------------------------------------------------------------------------

template<multistep_method Method, class Value, class JacobianPattern,
         class Resizer, class T>
struct step_checkpoint<cppde::multistepper<Method, Value, JacobianPattern, Resizer>, T> {
  using scalar_type  = T;
  using stepper_type = cppde::multistepper<Method, Value, JacobianPattern, Resizer>;
  using carry_type   = typename stepper_type::carry;

  static constexpr int max_order = stepper_type::max_order;

  carry_type     carry;
  std::vector<T> zn;       // (carry.q + 1) slots of n_states, slot-major
  std::vector<T> y;        // the state the corrector converged to
  std::size_t    n_states = 0;
  double         t  = 0.0;
  double         dt = 0.0;
  int            q_next = 1;    // order the controller picked for the next step
  double         eta    = 1.0;  // and the rescale it applied

  std::size_t n() const { return n_states; }

  // The slots beyond the state, whose cotangents the previous step receives on
  // its own carry out.
  std::size_t n_history() const {
    return static_cast<std::size_t>(carry.q) * n_states;
  }

  void capture(const stepper_type& st, const std::vector<Value>& x_in,
               double t_in, double dt_in)
  {
    st.save_carry(carry);
    n_states = x_in.size();
    zn.assign(static_cast<std::size_t>(carry.q + 1) * n_states, T());
    for (int j = 0; j <= carry.q; ++j) {
      const auto& slot = st.zn(j);
      for (std::size_t i = 0; i < n_states; ++i)
        zn[static_cast<std::size_t>(j) * n_states + i] =
            static_cast<T>(ad_traits::scalar_value(slot[i]));
    }
    t  = t_in;
    dt = dt_in;
  }

  template<class RStepper>
  void load(RStepper& rst,
            std::vector<codual<T>>& x_rev,
            std::vector<codual<T>>& history) const
  {
    rst.load_carry(carry, n_states);
    // Everything the previous replayed step left behind names nodes of a tape
    // that has been rewound, so a value read but not written this step would
    // carry a dependence on an unrelated node: a wrong derivative with a right
    // value, which nothing downstream notices.
    x_rev.assign(n_states, codual<T>());
    history.assign(n_history(), codual<T>());
    for (int j = 0; j <= carry.q; ++j) {
      auto& slot = rst.zn_mut(j);
      for (std::size_t i = 0; i < n_states; ++i) {
        const std::size_t k = static_cast<std::size_t>(j) * n_states + i;
        slot[i] = codual<T>(zn[k]);
        slot[i].independent();
        if (j == 0) x_rev[i] = slot[i];
        else        history[k - n_states] = slot[i];
      }
    }
  }

  template<class RStepper>
  void finish(RStepper& rst, const std::vector<codual<T>>& /*xout*/,
              std::vector<codual<T>>& carry_out) const
  {
    rst.complete_step();
    // The Nordsieck interpolant is anchored at tn_current, which the controller
    // sets here and the replay has to as well.
    rst.set_tn_current(t + dt);
    rst.prepare_dense_output();
    // An order increase reads the top slot, which the controller fills with the
    // accumulated correction first. Without it the new slot is whatever the
    // previous replayed step left there.
    if (q_next > carry.q) rst.save_acor_to_zn_qmax();
    if (q_next != carry.q) rst.set_order_for_next_step(q_next);
    if (std::abs(eta - 1.0) > 1e-14) rst.rescale(static_cast<T>(eta));

    const int q_out = rst.current_order();
    carry_out.assign(static_cast<std::size_t>(q_out + 1) * n_states, codual<T>());
    for (int j = 0; j <= q_out; ++j) {
      const auto& slot = rst.zn(j);
      for (std::size_t i = 0; i < n_states; ++i)
        carry_out[static_cast<std::size_t>(j) * n_states + i] = slot[i];
    }
  }
};

// ============================================================================
//  equation_solver
//
//  The matrix an implicit method's equations are linearised against, factorised
//  once per reverse step and used in both directions: forward to recover a value
//  a stage solved for, transposed to carry a cotangent back through it.
//
//  Fresh rather than the forward run's own, which belongs to its iteration and is
//  stale by design; MSBP and MSBJ are exactly that staleness. res_scale is what
//  the residual's derivative in its solution is as a multiple of W: the step size
//  for a corrector written in Nordsieck form, one for a Rosenbrock stage.
// ============================================================================

template<class JacFunc, class T = double, bool Sparse = false>
class equation_solver {
public:
  explicit equation_solver(JacFunc& jac) : m_jac(&jac) {}

  void prepare(const std::vector<T>& x, T t, T inv_gamma_dt, T res_scale = T(1)) {
    m_lu.resize(x);
    m_lu.call_jacobian(*m_jac, const_cast<std::vector<T>&>(x), t);
    m_lu.factorize_W(x.size(), inv_gamma_dt);
    m_scale = res_scale;
  }

  void forward(std::vector<T>& b) { m_lu.solve(b); }

  void transposed(std::vector<T>& b) {
    m_lu.solve_transposed(b);
    if (m_scale != T(1)) for (T& v : b) v /= m_scale;
  }

private:
  JacFunc*             m_jac;
  cppde::lu_W<T, Sparse> m_lu;
  T                    m_scale = T(1);
};

// ============================================================================
//  replay_interpolate
//
//  Dense output on the replayed step, so an observation between step ends can
//  be seeded. It records onto the same tape as the step itself.
//
//  The one-step families take the endpoints and their times as arguments; the
//  multistepper evaluates its own Nordsieck snapshot and takes only the time.
//  The branch is on which of the two a stepper offers.
// ============================================================================

// Which of the three replay shapes a stepper offers. An explicit one-step method
// has neither: its step is do_step and nothing about it is implicit.
template<class S, class = void> struct has_corrector_replay : std::false_type {};
template<class S>
struct has_corrector_replay<S, std::void_t<decltype(std::declval<S&>().replay_outputs(
    std::declval<const typename S::state_type&>(),
    std::declval<typename S::state_type&>(),
    std::declval<typename S::state_type&>()))>
> : std::true_type {};

template<class S, class = void> struct has_stage_solves : std::false_type {};
template<class S>
struct has_stage_solves<S, std::void_t<decltype(std::declval<const S&>().replay_inv_gamma_dt(
    std::declval<typename S::time_type>()))>
> : std::true_type {};

// Whether the carry has to be read before the step rather than after it, which
// is the multistepper, whose history the step mutates in place.
template<class S, class = void> struct has_step_snapshot : std::false_type {};
template<class S>
struct has_step_snapshot<S, std::void_t<decltype(std::declval<S&>().set_step_snapshot(
    std::declval<typename S::step_snapshot>()))>
> : std::true_type {};

template<class Checkpoint, class RStepper, class State, class = void>
struct has_finish : std::false_type {};

template<class Checkpoint, class RStepper, class State>
struct has_finish<Checkpoint, RStepper, State,
    std::void_t<decltype(std::declval<const Checkpoint&>().finish(
        std::declval<RStepper&>(), std::declval<const State&>(),
        std::declval<State&>()))>
> : std::true_type {};

template<class RStepper, class Time, class State, class = void>
struct has_endpoint_calc_state : std::false_type {};

template<class RStepper, class Time, class State>
struct has_endpoint_calc_state<RStepper, Time, State,
    std::void_t<decltype(std::declval<RStepper&>().calc_state(
        std::declval<Time>(), std::declval<State&>(),
        std::declval<const State&>(), std::declval<Time>(),
        std::declval<const State&>(), std::declval<Time>()))>
> : std::true_type {};

template<class RStepper, class Time, class State>
inline void replay_interpolate(RStepper& st, Time t, State& x,
                               const State& x_old, Time t_old,
                               const State& x_new, Time t_new)
{
  if constexpr (has_endpoint_calc_state<RStepper, Time, State>::value) {
    st.calc_state(t, x, x_old, t_old, x_new, t_new);
  } else {
    (void)x_old; (void)t_old; (void)x_new; (void)t_new;
    st.eval_dense_into(t, x);
  }
}

// ============================================================================
//  step_recorder<Stepper, T>
//
//  Drives one reverse step. Stepper is the forward type; the replay runs on the
//  same method rebound to codual<T>. The checkpoint holds all stepper specifics.
//
//    begin()                drops the previous step's tape and adjoints
//    independent(v)         registers parameters, and whatever else the system
//                           functor reads, as tape inputs
//    record(sys, cp)        replays the step, leaving the stepper live
//    load(cp, dt_in)        the first half of record(): checkpoint in, tape
//                           inputs registered, nothing stepped
//    attempt(sys, dt)       one step at a size the caller derived. The last one
//                           wins. Two halves rather than one call because the
//                           control law turns the size the controller offered
//                           into the size the step took, and the attempts it
//                           threw away on the way are part of that dependence.
//    seed(w) / seed(i, w)   cotangent of the step end; repeated seeds add
//    seed_err(w)            cotangent of the embedded error estimate, the way
//                           the control law reaches back into the step
//    interpolate(t, x)      dense output inside the step, onto the same tape
//    sweep()                one backwards pass
//    wx(), wt(), wdt()      cotangents of the step start, its time and its size
//    accumulate(v, out)     adds the cotangents of registered inputs onto out
//
//  Seeding sits after record() so a caller may record more on top of the step;
//  stage 4 seeds through calc_state, whose interpolation is on the same tape.
//
//  tape_stepsize(false) drops t and h from the independents, which is the
//  frozen path: it computes exactly what the forward sensitivities do, and is
//  the verification switch, not a shipped mode.
// ============================================================================

template<class Stepper, class T = double>
class step_recorder {
public:
  using scalar_type     = T;
  using rev_type        = codual<T>;
  using rev_stepper     = typename Stepper::template rebind_value<rev_type>;
  using checkpoint_type = step_checkpoint<Stepper, T>;
  using tape_type       = codual_tape<T>;

  // One implicit equation on the tape: the solution it determines, the residual
  // that stands in for solving it, and the mark just above the residual's nodes.
  struct solve_point {
    std::vector<rev_type> y, res;
    std::size_t           mark = 0;
  };

  tape_type& tape() const { return codual_tape_for<T>(); }

  // Drops the previous step's nodes and adjoints. Parameter cotangents are read
  // out per step, so nothing has to survive this.
  void begin() {
    tape().rewind();
    m_n_implicit = 0;
    m_solve_fn = nullptr;
  }

  // Registers each element as a tape input. Use the system functor's own
  // parameter copy: that is what the right-hand side reads.
  void independent(std::vector<rev_type>& v) {
    for (std::size_t i = 0; i < v.size(); ++i) v[i].independent();
  }

  // Whether the step's time and size go on the tape. On by default: the
  // adjoint is meant to run through the step-size control.
  void tape_stepsize(bool on) { m_tape_stepsize = on; }

  // Replays the step. The stepper is left holding the recorded stages so
  // calc_state can be seeded on top.
  template<class RSys>
  void record(RSys& sys, const checkpoint_type& cp) {
    load(cp, static_cast<T>(cp.dt));
    attempt(sys, m_dt);
  }

  // Checkpoint in, tape inputs registered, nothing stepped. dt_in is the size
  // the step was entered with, which is the checkpoint's own size unless the
  // controller threw attempts away before it.
  void load(const checkpoint_type& cp, const T& dt_in) {
    m_cp_finish = &cp;
    cp.load(m_stepper, m_x, m_history);
    m_t  = rev_type(static_cast<T>(cp.t));
    m_dt = rev_type(dt_in);
    if (m_tape_stepsize) { m_t.independent(); m_dt.independent(); }
    const std::size_t n = m_x.size();
    m_xout.assign(n, rev_type());
    m_xerr.assign(n, rev_type());
  }

  // One attempt from the loaded step start, onto the tape. The last one wins:
  // its end state and error estimate are the step's.
  template<class RSys>
  void attempt(RSys& sys, const rev_type& dt) {
    m_stepper.do_step(sys, m_x, m_t, m_xout, dt, m_xerr);
    m_dt_used = dt;
    m_tend    = m_t + dt;
    close_step(m_cp_finish);
  }

  // Registers the solution of an implicit equation as tape inputs and returns the
  // point to record its residual into. A stepper that would have solved calls
  // this instead, writes the residual, and calls close_implicit().
  solve_point& open_implicit(const std::vector<T>& y_star) {
    // Grown once and then reused: a Rosenbrock step opens six of these, and a
    // fresh pair of vectors per equation per step would be the reverse pass's
    // largest source of allocation.
    if (m_n_implicit == m_implicit.size()) m_implicit.emplace_back();
    solve_point& sp = m_implicit[m_n_implicit++];
    sp.y.assign(y_star.size(), rev_type());
    for (std::size_t i = 0; i < y_star.size(); ++i) {
      sp.y[i] = rev_type(y_star[i]);
      sp.y[i].independent();
    }
    sp.res.assign(y_star.size(), rev_type());
    return sp;
  }

  // Marks the tape above the equation just recorded. Everything that reads its
  // solution afterwards lands above the mark, which is what lets the sweep
  // reach the equation with the solution's cotangent already complete.
  void close_implicit() { m_implicit[m_n_implicit - 1].mark = tape().size(); }

  // What a stepper's replay calls in place of each linear solve. On entry the
  // vector holds the right-hand side; the value the stage solved for comes back
  // from the same factorisation the sweep uses transposed, so nothing about the
  // stages has to be checkpointed.
  template<class Solver>
  struct staged_sink {
    step_recorder* rec;
    Solver*        solver;

    void begin(std::vector<rev_type>& v) {
      std::vector<T>& b = rec->solve_scratch();
      b.assign(v.size(), T());
      for (std::size_t i = 0; i < v.size(); ++i) b[i] = v[i].x();
      solver->forward(b);
      v = rec->open_implicit(b).y;
    }
    std::vector<rev_type>& residual() { return rec->live_implicit().res; }
    void done() { rec->close_implicit(); }
  };

  // A step whose stages are linear solves against one shared matrix. solver
  // carries that matrix factorised, forward() and transposed() against it.
  template<class RSys, class Solver>
  void attempt_staged(RSys& sys, const rev_type& dt, Solver& solver)
  {
    const std::size_t n = m_x.size();
    m_xout.assign(n, rev_type());
    m_xerr.assign(n, rev_type());
    staged_sink<Solver> sink{this, &solver};
    m_stepper.replay_step(sys, m_x, m_t, dt, m_xout, m_xerr, sink);

    m_solve_ctx = &solver;
    m_solve_fn  = [](void* ctx, std::vector<T>& b) {
      static_cast<Solver*>(ctx)->transposed(b);
    };

    m_dt_used = dt;
    m_tend    = m_t + dt;
    close_step(m_cp_finish);
  }

  // One attempt of an implicit method. The corrector is not iterated: y is the
  // value the forward run converged to and the equation it solves goes on the
  // tape instead, which the sweep closes by a transposed solve. solve_t is that
  // solve, taking a solution's cotangent to the seed for its equation.
  template<class RSys, class SolveT>
  void attempt_implicit(RSys& sys, const rev_type& dt,
                        const std::vector<T>& y_star, SolveT& solve_t)
  {
    const std::size_t n = m_x.size();
    solve_point& sp = open_implicit(y_star);
    m_stepper.replay_residual(sys, m_x, m_t, dt, sp.y, sp.res);
    close_implicit();
    // Read before the tail runs: the controller's rescale moves the stepper's
    // step size, and the equation's matrix belongs to the step that was taken.
    m_impl_gamma = static_cast<T>(m_stepper.gamma());
    m_impl_t_new = m_t.x() + static_cast<T>(m_stepper.h());
    m_xout.assign(n, rev_type());
    m_xerr.assign(n, rev_type());
    m_stepper.replay_outputs(m_implicit.back().y, m_xout, m_xerr);

    set_solver(solve_t);

    m_dt_used = dt;
    m_tend    = m_t + dt;
    close_step(m_cp_finish);
  }

  solve_point&    live_implicit()  { return m_implicit[m_n_implicit - 1]; }
  std::vector<T>& solve_scratch() { return m_solve_scratch; }

  // The transposed solve the sweep runs at every equation.
  template<class SolveT>
  void set_solver(SolveT& solve_t) {
    m_solve_ctx = &solve_t;
    m_solve_fn  = [](void* ctx, std::vector<T>& b) {
      static_cast<SolveT*>(ctx)->transposed(b);
    };
  }

  // The carry the next step reads. Equal to the step end for a one-step method,
  // and the whole Nordsieck history for the multistepper.
  const std::vector<rev_type>& carry_out() const { return m_carry_out; }

  // Cotangent of the carry: the state part on the step end, the rest on what
  // the checkpoint's finish() wrote beyond it.
  void seed_carry(const std::vector<T>& wx, const std::vector<T>& whist) {
    for (std::size_t i = 0; i < m_carry_out.size(); ++i) {
      const T w = (i < wx.size()) ? wx[i]
                : (i - wx.size() < whist.size()) ? whist[i - wx.size()] : T();
      m_carry_out[i].seed(w);
    }
  }

  // Dense output at a time inside the step, recorded on top of it. The step end
  // is t + h as an expression, so an interpolated seed reaches the step size the
  // same way the step end does.
  void interpolate(const T& t, std::vector<rev_type>& x) {
    x.assign(m_x.size(), rev_type());
    replay_interpolate(m_stepper, rev_type(t), x, m_x, m_t, m_xout, m_tend);
  }


  // Cotangent of the step end. Repeated seeds accumulate, so several outputs
  // reduce onto one sweep.
  void seed(const std::vector<T>& w) {
    for (std::size_t i = 0; i < m_xout.size() && i < w.size(); ++i)
      m_xout[i].seed(w[i]);
  }
  void seed(std::size_t i, const T& w) { m_xout[i].seed(w); }

  // Cotangent of the embedded error estimate. The controller reads xerr to pick
  // the next step size, so this is where that dependence re-enters the step.
  void seed_err(const std::vector<T>& w) {
    for (std::size_t i = 0; i < m_xerr.size() && i < w.size(); ++i)
      m_xerr[i].seed(w[i]);
  }

  void sweep() {
    // Newest equation first. Everything above one is swept before its solve
    // runs, so the solution's cotangent is complete; the solve's result seeds
    // the equation, and the next segment carries it into the one below.
    std::size_t hi = tape().size();
    for (std::size_t k = m_n_implicit; k-- > 0;) {
      solve_point& sp = m_implicit[k];
      tape().reverse(hi, sp.mark);
      m_solve_buf.assign(sp.y.size(), T());
      for (std::size_t i = 0; i < sp.y.size(); ++i) m_solve_buf[i] = sp.y[i].adjoint();
      m_solve_fn(m_solve_ctx, m_solve_buf);
      for (std::size_t i = 0; i < sp.res.size(); ++i) sp.res[i].seed(-m_solve_buf[i]);
      hi = sp.mark;
    }
    tape().reverse(hi, 0);
    m_wx.assign(m_x.size(), T());
    for (std::size_t i = 0; i < m_x.size(); ++i) m_wx[i] = m_x[i].adjoint();
    m_whistory.assign(m_history.size(), T());
    for (std::size_t i = 0; i < m_history.size(); ++i)
      m_whistory[i] = m_history[i].adjoint();
    m_wt  = m_t.adjoint();
    m_wdt = m_dt.adjoint();
  }

  // Valid after sweep(): the cotangents the previous step receives. wdt is what
  // the control law consumes, wt what the step before it adds to its own.
  const std::vector<T>& wx()        const { return m_wx; }
  const std::vector<T>& whistory()  const { return m_whistory; }
  const T&              wt()        const { return m_wt; }
  const T&              wdt()       const { return m_wdt; }

  // out += cotangents of v. Parameters are shared across steps, so their
  // cotangent is a sum over the trajectory.
  void accumulate(const std::vector<rev_type>& v, std::vector<T>& out) const {
    if (out.size() < v.size()) out.resize(v.size(), T());
    for (std::size_t i = 0; i < v.size(); ++i) out[i] += v[i].adjoint();
  }

  // The replayed step end, for comparison against the forward run's.
  const std::vector<rev_type>& xout() const { return m_xout; }
  const std::vector<rev_type>& xerr() const { return m_xerr; }
  const std::vector<rev_type>& xin()  const { return m_x; }
  rev_stepper&                 stepper()    { return m_stepper; }

  // The step's tape inputs and the size it actually took, for a caller that
  // builds the control law on top.
  const rev_type& t()       const { return m_t; }
  const rev_type& dt_in()   const { return m_dt; }
  const rev_type& dt_used() const { return m_dt_used; }
  const rev_type& t_end()   const { return m_tend; }

  // Where and with what scaling the corrector's equation has to be linearised,
  // valid after attempt_implicit and unaffected by the step's tail.
  const T& implicit_gamma() const { return m_impl_gamma; }
  const T& implicit_t_new() const { return m_impl_t_new; }

private:
  // The step's tail, where the stepper has one. Without it the carry out is the
  // step end itself, which a one-step method's next step reads unchanged.
  void close_step(const checkpoint_type* cp) {
    using state = std::vector<rev_type>;
    if constexpr (has_finish<checkpoint_type, rev_stepper, state>::value) {
      cp->finish(m_stepper, m_xout, m_carry_out);
    } else {
      (void)cp;
      m_carry_out = m_xout;
    }
  }

  rev_stepper            m_stepper;
  const checkpoint_type* m_cp_finish = nullptr;
  std::vector<rev_type>  m_x, m_history, m_xout, m_xerr, m_carry_out;
  std::vector<solve_point> m_implicit;      // one per equation, in tape order
  std::size_t              m_n_implicit = 0; // how many of them this step opened
  std::vector<T>           m_solve_buf, m_solve_scratch;
  T                      m_impl_gamma{}, m_impl_t_new{};
  void*                  m_solve_ctx = nullptr;
  void                 (*m_solve_fn)(void*, std::vector<T>&) = nullptr;
  rev_type               m_t, m_dt, m_dt_used, m_tend;
  std::vector<T>         m_wx, m_whistory;
  T                      m_wt{}, m_wdt{};
  bool                   m_tape_stepsize = true;
};

}  // namespace reverse
}  // namespace cppde

#endif  // CPPDE_REVERSE_STEP_HPP
