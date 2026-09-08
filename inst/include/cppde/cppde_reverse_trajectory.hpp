/*
 The trajectory backwards: checkpoint store, reverse loop, seeds at the
 observation times.

 The forward run integrates in double and drops one checkpoint per accepted
 step, plus a note of which observation fell where. The reverse run walks the
 steps newest first, replays each one under cppde::codual, seeds it with the
 cotangent the later step handed back and with whatever observation lands
 inside it, and sweeps. Memory is one step's tape plus the checkpoints.

 What comes out is the cotangent of the trajectory start and, summed over every
 step, the cotangent of the parameters. The parameter accumulator has no state
 dimension: it is the quadrature the continuous adjoint writes as an integral.

 Observations sit at interpolated times, not at step ends, so the dense output
 is on the tape with the step that carries it. An observation before the first
 step reaches the initial state directly.

 The step-size chain, dt_{k+1} = Ctrl(err_k), is on the tape as well. A step
 reads four things from the one before it, not one: the state, the time, the
 size the controller offered, and the PI memory err_old. The attempts the
 controller threw away are replayed too, because they are what turned the
 offered size into the size the step took.

 control_chain(false) drops the last three and keeps only the state. That is the
 frozen path: it computes exactly what the forward sensitivities compute, which
 makes it the oracle, and the difference between the two settings is the
 step-size term itself.

 Copyright (C) 2026 Simon Beyer
 */

#ifndef CPPDE_REVERSE_TRAJECTORY_HPP
#define CPPDE_REVERSE_TRAJECTORY_HPP

#include <cstddef>
#include <vector>

#include <cppde/cppde_onestep_controller.hpp>
#include <cppde/cppde_reverse_step.hpp>

namespace cppde {
namespace reverse {

// ============================================================================
//  The controller state a step is entered with, and the constants of its law.
//
//  Everything the PI controller reads that an earlier step wrote. dt_in is not
//  the checkpoint's size: attempts thrown away shrink it on the way, and how
//  many is `rejected`. Their sizes are not stored because the replay derives
//  them, which is the same reason the checkpoint holds no stage values.
//
//  This is the one-step law. The multistepper controls order as well and gets
//  its own pair with stage 3b.
// ============================================================================

struct control_state {
  double   dt_in      = 0.0;
  double   err_old    = 1.0;
  bool     first_step = false;   // no usable err_old
  unsigned rejected   = 0;       // attempts before the accepted one
};

struct control_params {
  double atol = 1e-6, rtol = 1e-6;
  double alpha = 0.0, beta = 0.0, safety = 0.9;
  double min_factor = 0.2, max_factor = 5.0;
  double order = 4.0;            // error order of the embedded pair

  // Reads them off a live one-step controller, so the replay cannot drift from
  // the settings the forward run used.
  template<class Controller>
  static control_params of(const Controller& c) {
    control_params p;
    p.atol = c.atol();   p.rtol = c.rtol();
    p.alpha = c.alpha(); p.beta = c.beta(); p.safety = c.safety();
    p.min_factor = c.min_factor(); p.max_factor = c.max_factor();
    p.order = static_cast<double>(Controller::order);
    return p;
  }
};

// ============================================================================
//  trajectory_store<Stepper, T>
//
//  What the forward run leaves behind. Filled through two calls, both cheap
//  enough to sit in the integration loop:
//
//    capture(stepper, x, t, dt)   after an accepted step, with the step start
//    observe(t)                   at every observer call
//
//  observe() records the number of steps taken before it, which is what ties an
//  observation to the step whose bracket contains it. Zero means the initial
//  state, before any step ran.
// ============================================================================

template<class Stepper, class T = double>
class trajectory_store {
public:
  using checkpoint_type = step_checkpoint<Stepper, T>;

  struct observation {
    double      t;
    std::size_t step;   // accepted steps before it; 0 is the trajectory start
  };

  void clear() { m_steps.clear(); m_controls.clear(); m_obs.clear(); }

  void reserve(std::size_t n_steps, std::size_t n_obs) {
    m_steps.reserve(n_steps);
    m_controls.reserve(n_steps);
    m_obs.reserve(n_obs);
  }

  // No default control state: an unset dt_in would replay the step at size zero
  // once the chain is on, and silently. step_collector is the caller.
  template<class Value>
  void capture(const Stepper& st, const std::vector<Value>& x,
               double t, double dt, const control_state& cs)
  {
    m_steps.emplace_back();
    m_steps.back().capture(st, x, t, dt);
    m_controls.push_back(cs);
  }

  // A checkpoint the collector filled itself, for a family whose carry cannot be
  // read off the stepper after the step.
  void push(const checkpoint_type& cp, const control_state& cs) {
    m_steps.push_back(cp);
    m_controls.push_back(cs);
  }

  void observe(double t) { m_obs.push_back(observation{t, m_steps.size()}); }

  std::size_t n_steps() const { return m_steps.size(); }
  std::size_t n_obs()   const { return m_obs.size(); }
  std::size_t n_states() const {
    return m_steps.empty() ? 0u : m_steps.front().n();
  }

  const checkpoint_type& step(std::size_t k)    const { return m_steps[k]; }
  const control_state&   control(std::size_t k) const { return m_controls[k]; }
  const observation&     obs(std::size_t i)     const { return m_obs[i]; }

  control_params&       params()       { return m_params; }
  const control_params& params() const { return m_params; }

private:
  std::vector<checkpoint_type> m_steps;
  std::vector<control_state>   m_controls;
  std::vector<observation>     m_obs;
  control_params               m_params;
};

// ============================================================================
//  step_collector
//
//  The step observer the dense driver calls, and the only place that knows how
//  to read a live stepper. It carries the controller state forward by one step:
//  what the controller holds after an accepted step is what the next one enters
//  with, and the first step's is handed in at construction.
//
//  An event restart resets the controller between two firings, so the carried
//  state would be wrong across one. That is stage 5, where the restart cuts the
//  chain anyway.
// ============================================================================

template<class DenseStepper, class Stepper, class T = double>
class step_collector {
public:
  step_collector(trajectory_store<Stepper, T>& store, DenseStepper& st, double dt0)
    : m_store(store), m_st(st)
  {
    if constexpr (!snapshot_family) {
      m_next.dt_in      = dt0;
      m_next.err_old    = 1.0;
      m_next.first_step = true;
      m_store.params() = control_params::of(st.controlled_stepper());
    } else {
      (void)dt0;
      st.controlled_stepper().stepper().set_step_snapshot(
          [this](double t, double h) { this->snapshot(t, h); });
    }
  }

  void operator()() {
    auto& ctl = m_st.controlled_stepper();

    if constexpr (snapshot_family) {
      // The carry came from the last snapshot before this acceptance; what the
      // step ended on and what the controller then chose come from here.
      m_pending.y.assign(m_st.current_state().begin(), m_st.current_state().end());
      m_pending.q_next = ctl.stepper().current_order();
      m_pending.eta    = static_cast<double>(ctl.stepper().hscale()) / m_pending.dt;
      control_state cs;
      cs.rejected = static_cast<unsigned>(ctl.n_rejected() - m_rejected_seen);
      m_rejected_seen = ctl.n_rejected();
      m_store.push(m_pending, cs);
      return;
    } else {
      const double t0 = m_st.previous_time();

      control_state cs = m_next;
      cs.rejected = static_cast<unsigned>(ctl.n_rejected() - m_rejected_seen);
      m_rejected_seen = ctl.n_rejected();

      // dt_old() and not current_time() - previous_time(): the controller
      // advances t by dt, and fl(t + dt) - t is not dt. The replay has to step
      // the size the forward run stepped, not a rounded version of it.
      m_store.capture(ctl.stepper(), m_st.previous_state(), t0, ctl.dt_old(), cs);

      m_next.dt_in      = m_st.current_time_step();
      m_next.err_old    = ctl.last_error();
      m_next.first_step = ctl.first_step();
    }
  }

private:
  static constexpr bool snapshot_family = has_step_snapshot<Stepper>::value;

  void snapshot(double t, double h) {
    if constexpr (snapshot_family) {
      auto& st = m_st.controlled_stepper().stepper();
      m_pending.capture(st, st.zn(0), t, h);
    } else {
      (void)t; (void)h;
    }
  }

  trajectory_store<Stepper, T>& m_store;
  DenseStepper&                 m_st;
  typename trajectory_store<Stepper, T>::checkpoint_type m_pending;
  control_state                 m_next;
  int                           m_rejected_seen = 0;
};

// ============================================================================
//  trajectory_recorder<Stepper, T>
//
//  One backwards pass over a filled store.
//
//    sweep(store, params, make_sys, seeds)
//
//  make_sys is called once per step with the step's own codual parameters,
//  already registered as tape inputs, and returns the system the replay
//  integrates. It has to be the same system the forward run used.
//
//  seeds is [n_obs, n_states] row-major, the cotangent of each observed state.
//  For a least-squares objective it is the residual weighted by the error
//  model, which is where the chain above the solver hands over.
//
//  After the sweep, wx0() is the cotangent of the trajectory start and wp()
//  that of the parameters.
// ============================================================================

template<class Stepper, class T = double>
class trajectory_recorder {
public:
  using store_type    = trajectory_store<Stepper, T>;
  using recorder_type = step_recorder<Stepper, T>;
  using rev_type      = codual<T>;

  // An explicit method asks for no equation solver, so the four-argument form
  // hands in one that is never called.
  struct no_solver {
    void forward(std::vector<T>&)    {}
    void transposed(std::vector<T>&) {}
  };

  template<class MakeSys>
  void sweep(const store_type& store, const std::vector<T>& params,
             MakeSys make_sys, const std::vector<T>& seeds)
  {
    no_solver none;
    sweep(store, params, make_sys, seeds, none);
  }

  template<class MakeSys, class Solver>
  void sweep(const store_type& store, const std::vector<T>& params,
             MakeSys make_sys, const std::vector<T>& seeds, Solver& solver)
  {
    const std::size_t n_x = store.n_states();
    const std::size_t n_p = params.size();

    m_wx.assign(n_x, T());
    m_whist.clear();
    m_wp.assign(n_p, T());
    m_wt.assign(store.n_steps(), T());
    m_wdt.assign(store.n_steps(), T());
    m_max_nodes = 0;
    m_x_obs.assign(store.n_obs() * n_x, T());

    std::size_t next_obs = store.n_obs();
    std::vector<rev_type> p, x_interp;

    // What the step after this one hands back through the control law.
    T w_dt_in = T(), w_t = T(), w_err_old = T();

    for (std::size_t k = store.n_steps(); k-- > 0;) {
      m_rec.begin();

      p.assign(n_p, rev_type());
      for (std::size_t j = 0; j < n_p; ++j) p[j] = rev_type(params[j]);
      m_rec.independent(p);

      auto sys = make_sys(p);
      const control_state& cs = store.control(k);

      rev_type err_old, dt_next, t_next, err_old_next;
      const bool chained = replay_one(sys, store, k, cs, solver, err_old,
                                      dt_next, t_next, err_old_next);

      // Observations inside this step, seeded through the dense output. They
      // are recorded on the step's own tape, so one sweep carries both.
      while (next_obs > 0 && store.obs(next_obs - 1).step == k + 1) {
        --next_obs;
        m_rec.interpolate(clamp_to_step(store.step(k), store.obs(next_obs).t),
                          x_interp);
        const T* w = seeds.data() + next_obs * n_x;
        for (std::size_t i = 0; i < n_x && i < x_interp.size(); ++i) {
          m_x_obs[next_obs * n_x + i] = x_interp[i].x();
          x_interp[i].seed(w[i]);
        }
      }

      // What the later step handed back: its carry, which is the step end for a
      // one-step method and the whole Nordsieck history for the multistepper,
      // and the three the control law carries on top of it.
      m_rec.seed_carry(m_wx, m_whist);
      if (chained) {
        dt_next.seed(w_dt_in);
        t_next.seed(w_t);
        err_old_next.seed(w_err_old);
      }
      m_rec.sweep();
      if (chained) {
        w_dt_in   = m_rec.wdt();
        w_t       = m_rec.wt();
        w_err_old = err_old.adjoint();
      }

      // Parameters are shared by every step, so their cotangent is a sum. The
      // copies inside the system name the same tape slots as p.
      m_rec.accumulate(p, m_wp);
      m_wt[k]  = m_rec.wt();
      m_wdt[k] = m_rec.wdt();
      m_wx     = m_rec.wx();
      m_whist  = m_rec.whistory();

      const std::size_t used = m_rec.tape().size();
      if (used > m_max_nodes) m_max_nodes = used;
    }

    // Whatever was observed before the first step is the initial state itself.
    while (next_obs > 0) {
      --next_obs;
      const T* w = seeds.data() + next_obs * n_x;
      for (std::size_t i = 0; i < n_x; ++i) m_wx[i] += w[i];
    }
  }

  // Whether the step-size chain across step boundaries goes on the tape. On is
  // the shipped state; off drops it and computes what the forward sensitivities
  // compute, which is the oracle.
  void control_chain(bool on) { m_control_chain = on; }

  const std::vector<T>& wx0() const { return m_wx; }
  const std::vector<T>& wp()  const { return m_wp; }

  // The rest of the trajectory start's carry, which for a multistep method is
  // the Nordsieck slots above the state and for a one-step method is empty.
  const std::vector<T>& whistory0() const { return m_whist; }

  // What the replay itself observed, [n_obs, n_states] row-major. It has to be
  // the forward run's own output: a replay that lands elsewhere differentiates
  // elsewhere, and nothing else in the sweep would say so.
  const std::vector<T>& replayed_obs() const { return m_x_obs; }

  // Per step, the cotangents the control law carries: the step's own time and
  // the size it was entered with.
  const std::vector<T>& wt()  const { return m_wt; }
  const std::vector<T>& wdt() const { return m_wdt; }

  // The largest tape any single step needed. The bound is a step, not a
  // trajectory: begin() rewinds before each replay and keeps the capacity, so
  // the tape allocates once and is reused for every step after the first.
  std::size_t max_tape_nodes() const { return m_max_nodes; }

  recorder_type& step_rec() { return m_rec; }

private:
  // Dispatches the step onto the shape its method has: an implicit corrector,
  // a set of linear stage solves, or neither. Returns whether the control law
  // went on the tape with it, which only the explicit one-step path does today.
  template<class RSys, class Solver>
  bool replay_one(RSys& sys, const store_type& store, std::size_t k,
                  const control_state& cs, Solver& solver, rev_type& err_old,
                  rev_type& dt_next, rev_type& t_next, rev_type& err_old_next)
  {
    using rstep = typename recorder_type::rev_stepper;
    const auto& cp = store.step(k);

    if constexpr (has_stage_solves<rstep>::value) {
      m_rec.load(cp, static_cast<T>(cp.dt));
      solver.prepare(cp.x, static_cast<T>(cp.t),
                     m_rec.stepper().replay_inv_gamma_dt(cp.dt));
      m_rec.attempt_staged(sys, m_rec.dt_in(), solver);
      return false;
    } else if constexpr (has_corrector_replay<rstep>::value) {
      m_rec.load(cp, static_cast<T>(cp.dt));
      m_rec.attempt_implicit(sys, m_rec.dt_in(), cp.y, solver);
      // The corrector is written in Nordsieck form, so its derivative in the
      // solution is gamma times the matrix the solver factorised.
      const T gamma = m_rec.implicit_gamma();
      solver.prepare(cp.y, m_rec.implicit_t_new(), T(1) / gamma, gamma);
      return false;
    } else if (m_control_chain) {
      replay_controlled(sys, store, k, cs, err_old, dt_next, t_next,
                        err_old_next);
      return true;
    } else {
      m_rec.record(sys, cp);
      return false;
    }
  }

  // One step with its control law on the tape: the attempts the controller threw
  // away, the accepted one, and the three values the next step reads. The law is
  // onestep_controller's, called through the same functions it calls itself.
  template<class RSys>
  void replay_controlled(RSys& sys, const store_type& store, std::size_t k,
                         const control_state& cs, rev_type& err_old,
                         rev_type& dt_next, rev_type& t_next,
                         rev_type& err_old_next)
  {
    const control_params& par = store.params();

    m_rec.load(store.step(k), static_cast<T>(cs.dt_in));
    err_old = rev_type(static_cast<T>(cs.err_old));
    err_old.independent();

    rev_type dt = m_rec.dt_in();
    for (unsigned j = 0; j < cs.rejected; ++j) {
      m_rec.attempt(sys, dt);
      dt = dt * onestep_detail::reject_factor(step_error(par), par.order,
                                              par.safety, par.min_factor);
    }
    m_rec.attempt(sys, dt);

    const rev_type err = step_error(par);
    // pi_form and cap_at_one are what update_stepsize reads off m_first_step and
    // m_last_rejected: a step with a thrown-away attempt behind it has the flag.
    const rev_type factor = onestep_detail::accept_factor(
        err, err_old, !(cs.first_step || cs.rejected > 0), cs.rejected > 0,
        par.order, par.alpha, par.beta, par.safety, par.min_factor,
        par.max_factor);

    t_next       = m_rec.t_end();
    dt_next      = m_rec.dt_used() * factor;
    err_old_next = cppde::max(rev_type(T(0.01)), err);
  }

  // The controller's error norm on the replayed attempt, floored the way
  // try_step floors it before it reaches the control law.
  rev_type step_error(const control_params& par) const {
    const rev_type e = onestep_detail::wrms_state(
        m_rec.xout(), m_rec.xin(), m_rec.xerr(), par.atol, par.rtol,
        [](const rev_type& v) { return v; });
    return cppde::max(e, T(1e-15));
  }

  // The forward loop observes a time before the step bracket at the bracket
  // start instead, which happens after an event restart. Clamping reproduces
  // that branch; inside the bracket, which is every other case, it does nothing.
  template<class Checkpoint>
  static T clamp_to_step(const Checkpoint& cp, double t) {
    const double a = cp.t, b = cp.t + cp.dt;
    const double lo = a < b ? a : b;
    const double hi = a < b ? b : a;
    return static_cast<T>(t < lo ? lo : (t > hi ? hi : t));
  }

  recorder_type  m_rec;
  std::vector<T> m_wx, m_whist, m_wp, m_wt, m_wdt, m_x_obs;
  std::size_t    m_max_nodes = 0;
  bool           m_control_chain = true;
};

}  // namespace reverse
}  // namespace cppde

#endif  // CPPDE_REVERSE_TRAJECTORY_HPP
