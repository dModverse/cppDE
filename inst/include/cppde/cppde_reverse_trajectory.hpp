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

 The step grid is the forward run's and is not differentiated. A step reads one
 thing from the one before it, the state; the size it took is a constant read
 off its checkpoint. This is Bock's internal numerical differentiation: the
 nominal run adapts freely, and the derivative is taken of the scheme that run
 actually applied.

 Differentiating the controller instead puts spurious derivatives of the time
 steps on the tape, which make the discrete adjoint inconsistent with the
 adjoint ODE. cppDE did that behind a switch until 2026-09-09; the term was
 measured at O(tol) and the switch is gone. See dev/adjoint-plan.md.

 Copyright (C) 2026 Simon Beyer
 */

#ifndef CPPDE_REVERSE_TRAJECTORY_HPP
#define CPPDE_REVERSE_TRAJECTORY_HPP

#include <cstddef>
#include <vector>

#include <cppde/cppde_event_engine.hpp>
#include <cppde/cppde_onestep_controller.hpp>
#include <cppde/cppde_profiler.hpp>
#include <cppde/cppde_adjoint_step.hpp>
#include <cppde/cppde_reverse_step.hpp>
#include <cppde/cppde_saltation.hpp>

namespace cppde {
namespace reverse {

// ============================================================================
//  event_record<T>
//
//  One intervention between two steps, as the reverse mode needs it: which
//  events fired, the state on both sides of the jump, and the size the stepper
//  was restarted with.
//
//  Which events fired and when the root sat is a control decision and is read
//  back, not decided again. What is replayed is the arithmetic: the saltation
//  sandwich for a root event, the reset for a fixed one, and the restart, which
//  for a multistep method rebuilds the whole Nordsieck history out of the
//  post-jump state and so is the only thing the carry then depends on.
// ============================================================================

template<class T>
struct event_record {
  static constexpr std::size_t npos = static_cast<std::size_t>(-1);

  std::size_t after_step = 0;    // accepted steps before it
  double      t          = 0.0;  // where the jump sits
  double      t_before   = 0.0;  // where the state entering it was read
  bool        root       = false;
  bool        restart    = false;
  double      dt_restart = 0.0;

  std::vector<T> x_before, x_after;
  std::vector<cppde::detail::TriggeredEvent> triggered;
  // Root conditions a fixed jump switched on, in the order the engine applied
  // them on the event surface. Empty unless a model has both kinds.
  std::vector<std::size_t> switched;
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
    // The event whose post-jump state this is, rather than an interpolation
    // inside the step above. npos for every ordinary observation.
    std::size_t event = event_record<T>::npos;
  };

  void clear() {
    m_steps.clear(); m_obs.clear(); m_events.clear();
  }

  void reserve(std::size_t n_steps, std::size_t n_obs) {
    m_steps.reserve(n_steps);
    m_obs.reserve(n_obs);
  }

  template<class Value>
  void capture(const Stepper& st, const std::vector<Value>& x,
               double t, double dt)
  {
    auto _tp = m_prof.timer(cppde::prof_cat::rev_checkpoint);
    m_steps.emplace_back();
    m_steps.back().capture(st, x, t, dt);
  }

  // Where the store's own time went. Empty and free unless CPPDE_PROFILE.
  const cppde::profiler& prof() const { return m_prof; }

  // A checkpoint the collector filled itself, for a family whose carry cannot be
  // read off the stepper after the step.
  void push(const checkpoint_type& cp) { m_steps.push_back(cp); }

  // The post-jump state is observed at the time the jump sits at, before or
  // after the event note depending on which site fired it, so the marking runs
  // in both directions: push_event scans backwards over what already stands at
  // that time, observe() forwards while the window is open. A root event also
  // observes just before the jump, at t minus a whisker, which is why the
  // backward scan stops on the first time that differs.
  void observe(double t) {
    std::size_t ev = event_record<T>::npos;
    if (m_open_event != event_record<T>::npos) {
      if (t == m_open_t) ev = m_open_event;
      else               m_open_event = event_record<T>::npos;
    }
    m_obs.push_back(observation{t, m_steps.size(), ev});
  }

  void push_event(const event_record<T>& e) {
    m_events.push_back(e);
    const std::size_t idx = m_events.size() - 1;
    for (std::size_t i = m_obs.size(); i-- > 0;) {
      if (m_obs[i].t != e.t) break;
      m_obs[i].event = idx;
    }
    m_open_event = idx;
    m_open_t     = e.t;
  }

  std::size_t n_events() const { return m_events.size(); }
  const event_record<T>& event(std::size_t i) const { return m_events[i]; }

  // The intervention a step is entered through, npos where the step reads the
  // one before it directly. Linear over a list that is empty on most models.
  std::size_t event_before(std::size_t step) const {
    for (std::size_t i = 0; i < m_events.size(); ++i)
      if (m_events[i].after_step == step) return i;
    return event_record<T>::npos;
  }

  std::size_t n_steps() const { return m_steps.size(); }
  std::size_t n_obs()   const { return m_obs.size(); }
  std::size_t n_states() const {
    return m_steps.empty() ? 0u : m_steps.front().n();
  }

  const checkpoint_type& step(std::size_t k)    const { return m_steps[k]; }
  // The collector fills a checkpoint's tail record only once the step after it
  // has been accepted, which is when the operations in between are complete.
  checkpoint_type&       step_mut(std::size_t k)       { return m_steps[k]; }
  const observation&     obs(std::size_t i)     const { return m_obs[i]; }

private:
  std::vector<checkpoint_type> m_steps;
  std::vector<observation>     m_obs;
  std::vector<event_record<T>> m_events;
  std::size_t                  m_open_event = event_record<T>::npos;
  double                       m_open_t     = 0.0;
  cppde::profiler              m_prof;
};

// ============================================================================
//  step_collector
//
//  The step observer the dense driver calls, and the only place that knows how
//  to read a live stepper. The controller's own state is not recorded: the grid
//  it produced is a constant to the reverse pass, so only the step and what it
//  ended on matter.
// ============================================================================

template<class DenseStepper, class Stepper, class T = double>
class step_collector {
public:
  step_collector(trajectory_store<Stepper, T>& store, DenseStepper& st, double dt0)
    : m_store(store), m_st(st)
  {
    (void)dt0;   // kept for call-site compatibility; the grid needs no seed
    if constexpr (snapshot_family) {
      st.controlled_stepper().stepper().set_step_snapshot(
          [this](double t, double h) { this->snapshot(t, h); });
      st.controlled_stepper().stepper().set_history_log(&m_hlog);
    }
  }

  // The engine's event observer. Handed to integrate_times_dense beside the
  // step one; both write into the same store.
  std::function<void(const cppde::detail::event_note<
                       typename DenseStepper::state_type>&)>
  event_observer() {
    return [this](const auto& e) { this->on_event(e); };
  }

  void operator()() {
    auto& ctl = m_st.controlled_stepper();

    if constexpr (snapshot_family) {
      // The carry came from the last snapshot before this acceptance; what the
      // step ended on and what the controller then chose come from here.
      m_pending.y.assign(m_st.current_state().begin(), m_st.current_state().end());
      m_pending.q_next = ctl.stepper().current_order();
      m_pending.eta    = static_cast<double>(ctl.stepper().hscale()) / m_pending.dt;
      hand_over_history();
      m_store.push(m_pending);
      return;
    } else {
      // dt_old() and not current_time() - previous_time(): the controller
      // advances t by dt, and fl(t + dt) - t is not dt. The replay has to step
      // the size the forward run stepped, not a rounded version of it.
      m_store.capture(ctl.stepper(), m_st.previous_state(),
                      m_st.previous_time(), ctl.dt_old());
    }
  }

private:
  static constexpr bool snapshot_family = has_step_snapshot<Stepper>::value;

  // The history operations logged since the previous acceptance are the
  // previous step's tail, then whatever attempts this step threw away, then
  // this step's own completion and tail. complete_step writes the cut: what
  // stands before it is what the previous checkpoint has to replay, what stands
  // after it opens the next window.
  void hand_over_history() {
    if constexpr (snapshot_family) {
      std::size_t cut = m_hlog.size();
      for (std::size_t i = m_hlog.size(); i-- > 0;)
        if (m_hlog[i].op == history_op::complete) { cut = i; break; }
      if (m_store.n_steps() > 0) {
        auto& prev = m_store.step_mut(m_store.n_steps() - 1);
        prev.ops.assign(m_hlog.begin(), m_hlog.begin() + cut);
        prev.ops_recorded = true;
      }
      m_hlog.erase(m_hlog.begin(),
                   m_hlog.begin() + (cut < m_hlog.size() ? cut + 1 : cut));
    }
  }

public:
  // One intervention, turned into a record. The store does the marking of the
  // observation the jump produced, so an observer only ever needs the store.
  template<class Note>
  void on_event(const Note& e) {
    event_record<T> r;
    r.after_step = m_store.n_steps();
    r.t          = e.t;
    r.t_before   = e.t_before;
    r.root       = e.root;
    r.restart    = e.restart;
    r.dt_restart = e.dt_restart;
    if (e.x_before) copy_state(*e.x_before, r.x_before);
    if (e.x_after)  copy_state(*e.x_after,  r.x_after);
    if (e.triggered) r.triggered = *e.triggered;
    if (e.switched)  r.switched  = *e.switched;
    m_store.push_event(r);
  }

  template<class State>
  static void copy_state(const State& in, std::vector<T>& out) {
    out.resize(in.size());
    for (std::size_t i = 0; i < in.size(); ++i)
      out[i] = static_cast<T>(ad_traits::scalar_value(in[i]));
  }

private:
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
  history_log                   m_hlog;
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
  using rev_state     = std::vector<rev_type>;

  // The events a model carries, rebuilt on the reverse scalar type. A jump's
  // value and, in dMod2, its time depend on the parameters, so they are built
  // from the step's own codual parameter copies and their cotangents land on wp
  // like the right-hand side's.
  struct event_set {
    std::vector<cppde::detail::FixedEvent<rev_state, rev_type>> fixed;
    std::vector<cppde::detail::RootEvent<rev_state, rev_type>>  root;
  };

  // A model without events hands in this, and no event replay ever runs.
  struct no_events {
    event_set operator()(const rev_state&) const { return event_set{}; }
  };

  // An explicit method asks for no equation solver, so the shorter forms hand
  // in one that is never called.
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
    no_events none;
    sweep(store, params, make_sys, none, seeds, solver);
  }

  // An explicit method with events: no equation to solve, jumps to replay.
  template<class MakeSys, class MakeEvents>
  void sweep(const store_type& store, const std::vector<T>& params,
             MakeSys make_sys, MakeEvents make_events, const std::vector<T>& seeds)
  {
    no_solver none;
    sweep(store, params, make_sys, make_events, seeds, none);
  }

  template<class MakeSys, class MakeEvents, class Solver>
  void sweep(const store_type& store, const std::vector<T>& params,
             MakeSys make_sys, MakeEvents make_events,
             const std::vector<T>& seeds, Solver& solver)
  {
    const std::size_t n_x = store.n_states();
    const std::size_t n_p = params.size();

    m_wx.assign(n_x, T());
    m_whist.clear();
    m_wp.assign(n_p, T());
    m_wt.assign(store.n_steps(), T());
    m_wdt.assign(store.n_steps(), T());
    m_lam.assign(m_trace_lambda ? store.n_steps() * n_x : 0u, T());
    m_eta.assign(m_trace_lambda ? store.n_steps() : 0u, T());
    m_max_nodes = 0;
    m_pending.clear();
    m_pending_interp = false;
    m_x_obs.assign(store.n_obs() * n_x, T());

    std::size_t next_obs = store.n_obs();
    std::vector<rev_type> p, x_interp;

    for (std::size_t k = store.n_steps(); k-- > 0;) {
      m_rec.begin();

      p.assign(n_p, rev_type());
      for (std::size_t j = 0; j < n_p; ++j) p[j] = rev_type(params[j]);
      m_rec.independent(p);

      auto sys = make_sys(p);
      { auto _tp = m_prof.timer(cppde::prof_cat::rev_replay);
        replay_one(sys, store, k, solver); }

      // Observations inside this step, seeded through the dense output. They
      // are recorded on the step's own tape, so one sweep carries both. One
      // that an event produced is not an interpolation and is seeded there.
      while (next_obs > 0 && store.obs(next_obs - 1).step == k + 1) {
        --next_obs;
        if (store.obs(next_obs).event != event_record<T>::npos) continue;
        auto _tp = m_prof.timer(cppde::prof_cat::rev_interp);
        m_rec.interpolate(static_cast<T>(cppde::adjoint::clamp_to_step(
                              store.step(k), store.obs(next_obs).t)),
                          x_interp);
        const T* w = seeds.data() + next_obs * n_x;
        for (std::size_t i = 0; i < n_x && i < x_interp.size(); ++i) {
          m_x_obs[next_obs * n_x + i] = x_interp[i].x();
          x_interp[i].seed(w[i]);
        }
      }

      // What the later step handed back. Ordinarily that is its carry, read off
      // this step's end; across an intervention it is a cotangent the event
      // replay put on a point inside this step, and the carry is not read at
      // all, because the restart threw it away.
      if (m_pending_interp) {
        m_rec.interpolate(static_cast<T>(cppde::adjoint::clamp_to_step(
                              store.step(k), m_pending_t)), x_interp);
        for (std::size_t i = 0; i < x_interp.size() && i < m_pending.size(); ++i)
          x_interp[i].seed(m_pending[i]);
        m_pending_interp = false;
      } else {
        m_rec.seed_carry(m_wx, m_whist);
      }
      { auto _tp = m_prof.timer(cppde::prof_cat::rev_sweep);
        m_rec.sweep(); }

      // Parameters are shared by every step, so their cotangent is a sum. The
      // copies inside the system name the same tape slots as p.
      m_rec.accumulate(p, m_wp);
      m_wt[k]  = m_rec.wt();
      m_wdt[k] = m_rec.wdt();
      m_wx     = m_rec.wx();
      m_whist  = m_rec.whistory();
      if (m_trace_lambda) {
        for (std::size_t i = 0; i < n_x && i < m_wx.size(); ++i)
          m_lam[k * n_x + i] = m_wx[i];
        // The refinement indicator: lambda at the step end against the step's
        // own error estimate. wout() and not xout().adjoint(), whose adjoint the
        // sweep zeroes for a corrector method.
        T e = T();
        const std::vector<T>& lo = m_rec.wout();
        const std::vector<rev_type>& xe = m_rec.xerr();
        for (std::size_t i = 0; i < xe.size() && i < lo.size(); ++i)
          e += lo[i] * xe[i].x();
        m_eta[k] = e * m_rec.error_scale();
      }

      const std::size_t used = m_rec.tape().size();
      if (used > m_max_nodes) m_max_nodes = used;

      // The boundary this step was entered through. For a multistep method the
      // carry there was built by initialize() out of one state, at the run's
      // start and again after every event, so the cotangent of the whole
      // Nordsieck array reduces to a cotangent of that state. An intervention
      // then carries it back across the jump.
      const std::size_t ei = store.event_before(k);
      if (ei != event_record<T>::npos || k == 0) {
        const event_record<T>* e =
            (ei != event_record<T>::npos) ? &store.event(ei) : nullptr;
        replay_boundary(store, k, ei, e, params, make_sys, make_events, seeds);
      }
    }

    // Whatever was observed before the first step is the initial state itself.
    while (next_obs > 0) {
      --next_obs;
      if (store.obs(next_obs).event != event_record<T>::npos) continue;
      const T* w = seeds.data() + next_obs * n_x;
      for (std::size_t i = 0; i < n_x; ++i) m_wx[i] += w[i];
    }
  }

  const std::vector<T>& wx0() const { return m_wx; }
  const std::vector<T>& wp()  const { return m_wp; }

  // The rest of the trajectory start's carry, which for a multistep method is
  // the Nordsieck slots above the state and for a one-step method is empty.
  const std::vector<T>& whistory0() const { return m_whist; }

  // What the replay itself observed, [n_obs, n_states] row-major. It has to be
  // the forward run's own output: a replay that lands elsewhere differentiates
  // elsewhere, and nothing else in the sweep would say so.
  const std::vector<T>& replayed_obs() const { return m_x_obs; }

  // Per step, the cotangents of its own time and of the size it was entered
  // with. Diagnostic: not chained through the control law, see replay_one.
  const std::vector<T>& wt()  const { return m_wt; }
  const std::vector<T>& wdt() const { return m_wdt; }

  // lambda on the forward run's grid, [n_steps, n_states] row-major, row k the
  // state cotangent at step k's start. Off by default: it is another store's
  // worth of doubles and only the diagnostics read it.
  void trace_lambda(bool on) { m_trace_lambda = on; }
  const std::vector<T>& lambda() const { return m_lam; }

  // Per step, lambda^T e_k. This and not |dJ/dh_k| h_k is the refinement
  // indicator: wdt() is the transport derivative, lambda^T f, and is the size of
  // J rather than the size of the error.
  const std::vector<T>& eta() const { return m_eta; }

  // The largest tape any single step needed. The bound is a step, not a
  // trajectory: begin() rewinds before each replay and keeps the capacity, so
  // the tape allocates once and is reused for every step after the first.
  std::size_t max_tape_nodes() const { return m_max_nodes; }

  recorder_type& step_rec() { return m_rec; }

private:
  // Dispatches the step onto the shape its method has: an implicit corrector,
  // a set of linear stage solves, or neither.
  //
  // The controller stays off the tape in every shape. h_k is held at what the
  // forward run stepped and seeded nowhere, which is Bock's IND: differentiate
  // the scheme the adaptive decisions produced, not the decisions.
  template<class RSys, class Solver>
  void replay_one(RSys& sys, const store_type& store, std::size_t k,
                  Solver& solver)
  {
    using rstep = typename recorder_type::rev_stepper;
    const auto& cp = store.step(k);

    if constexpr (has_stage_solves<rstep>::value) {
      m_rec.load(cp, static_cast<T>(cp.dt));
      solver.prepare(cp.x, static_cast<T>(cp.t),
                     m_rec.stepper().replay_inv_gamma_dt(cp.dt));
      m_rec.attempt_staged(sys, m_rec.dt_in(), solver);
    } else if constexpr (has_corrector_replay<rstep>::value) {
      m_rec.load(cp, static_cast<T>(cp.dt));
      m_rec.attempt_implicit(sys, m_rec.dt_in(), cp.y, solver);
      // The corrector is written in Nordsieck form, so its derivative in the
      // solution is gamma times the matrix the solver factorised.
      const T gamma = m_rec.implicit_gamma();
      solver.prepare(cp.y, m_rec.implicit_t_new(), T(1) / gamma, gamma);
    } else {
      m_rec.record(sys, cp);
    }
  }

  // ------------------------------------------------------------------------
  //  The boundary a step was entered through, backwards.
  //
  //  Two maps, in the order the forward run applied them and so swept in
  //  reverse: the restart, which for a multistep method builds the whole
  //  Nordsieck history out of one state and therefore collapses the carry's
  //  cotangent onto that state, and the jump, which carries it back across the
  //  discontinuity. Both go on a tape of their own, rewound after; together
  //  they are smaller than one step.
  //
  //  Where there is no event and the step is not the first, nothing runs: the
  //  carry chains straight into the step before.
  // ------------------------------------------------------------------------
  template<class MakeSys, class MakeEvents>
  void replay_boundary(const store_type& store, std::size_t k, std::size_t ei,
                       const event_record<T>* e, const std::vector<T>& params,
                       MakeSys make_sys, MakeEvents make_events,
                       const std::vector<T>& seeds)
  {
    constexpr bool multistep = has_step_snapshot<Stepper>::value;
    const std::size_t n_x = store.n_states();
    const std::size_t n_p = params.size();
    const bool has_jump   = (e != nullptr);
    if constexpr (!multistep) { if (!has_jump) return; }

    codual_tape<T>& tp = codual_tape_for<T>();
    tp.rewind();

    std::vector<rev_type> p(n_p);
    for (std::size_t j = 0; j < n_p; ++j) {
      p[j] = rev_type(params[j]);
      p[j].independent();
    }
    auto sys = make_sys(p);

    // The state the restart reads, which is the state the jump ends on where
    // there is one and the step's own start where there is not.
    std::vector<rev_type> xb, xa(n_x);
    if (has_jump) {
      auto ev = make_events(p);
      xb.assign(n_x, rev_type());
      for (std::size_t i = 0; i < n_x; ++i) {
        xb[i] = rev_type(e->x_before[i]);
        xb[i].independent();
      }
      xa = xb;
      apply_jump(xa, xb, *e, sys, ev);
    } else {
      const T* x0 = store.step(k).start_state();
      for (std::size_t i = 0; i < n_x; ++i) {
        xa[i] = rev_type(x0[i]);
        xa[i].independent();
      }
    }

    // The restart. Its time and size are what the forward run used: the first
    // checkpoint's own at the trajectory start, the size the engine
    // re-estimated after an event, both control decisions.
    if constexpr (multistep) {
      const double t_r  = has_jump ? e->t : store.step(k).t;
      const double dt_r = (has_jump && e->restart) ? e->dt_restart
                                                   : store.step(k).carry.h;
      auto& rst = m_rec.stepper();
      std::vector<rev_type> f0(n_x);
      sys.first(xa, f0, rev_type(t_r));
      rst.initialize(xa, rev_type(t_r), f0, rev_type(dt_r));
      // Slot by slot, the state part off wx and the rest off whistory, which is
      // how the checkpoint above registered them.
      const std::size_t n_slot = 1 + (n_x ? m_whist.size() / n_x : 0);
      for (std::size_t j = 0; j < n_slot; ++j) {
        const auto& slot = rst.zn(static_cast<int>(j));
        for (std::size_t i = 0; i < n_x; ++i) {
          const std::size_t idx = j * n_x + i;
          slot[i].seed(idx < n_x ? m_wx[i] : m_whist[idx - n_x]);
        }
      }
    } else {
      for (std::size_t i = 0; i < n_x && i < m_wx.size(); ++i) xa[i].seed(m_wx[i]);
    }

    // Observations the jump produced. They are values, not interpolations, so
    // nothing above has seeded them.
    if (has_jump) {
      for (std::size_t o = 0; o < store.n_obs(); ++o) {
        if (store.obs(o).event != ei) continue;
        for (std::size_t i = 0; i < n_x; ++i) {
          m_x_obs[o * n_x + i] = xa[i].x();
          xa[i].seed(seeds[o * n_x + i]);
        }
      }
    }

    tp.reverse();

    m_wx.assign(n_x, T());
    if (has_jump) {
      for (std::size_t i = 0; i < n_x; ++i) m_wx[i] = xb[i].adjoint();
      // The state entering the jump was read off the dense output of the step
      // below it, so that is where its cotangent goes. At the run's own start
      // there is no such step and it is the initial state's.
      m_pending_interp = (e->after_step > 0);
      m_pending_t      = e->t_before;
      m_pending        = m_wx;
    } else {
      for (std::size_t i = 0; i < n_x; ++i) m_wx[i] = xa[i].adjoint();
    }
    m_whist.clear();
    m_rec.accumulate(p, m_wp);

    const std::size_t used = tp.size();
    if (used > m_max_nodes) m_max_nodes = used;
  }

  // The reset itself, on the reverse scalar type. Which events fired is read
  // back; what they did is replayed through the very functions the forward run
  // called, so the two cannot drift apart.
  template<class RSys, class Events>
  static void apply_jump(std::vector<rev_type>& x,
                         const std::vector<rev_type>& x_before,
                         const event_record<T>& e, RSys& sys, Events& ev)
  {
    if (e.root) {
      cppde::detail::saltation_root_analytical_batch(
          x, x_before, rev_type(e.t), sys, ev.root, e.triggered);
    } else {
      // Root conditions the jump switched on belong on the same surface, which
      // is what the engine's at_surface hook does forwards.
      auto at_surface = [&](std::vector<rev_type>& xs, const rev_type& ts) {
        for (std::size_t i : e.switched)
          cppde::detail::apply_event_action(xs, xs, ts, ev.root[i]);
      };
      cppde::detail::apply_fixed_events_at_time(x, rev_type(e.t), ev.fixed,
                                                sys, at_surface);
    }
  }

  recorder_type  m_rec;
  std::vector<T> m_wx, m_whist, m_wp, m_wt, m_wdt, m_x_obs, m_lam, m_eta;
  // What an intervention handed back: a cotangent on a point inside the step
  // below it rather than on that step's carry.
  std::vector<T> m_pending;
  double         m_pending_t = 0.0;
  bool           m_pending_interp = false;
  std::size_t    m_max_nodes = 0;
  bool           m_trace_lambda = false;
  cppde::profiler m_prof;

public:
  // Per-category timings of the sweep and of the store, to stderr. Compiled
  // away without CPPDE_PROFILE.
  void report_profile(const store_type& store) const {
    m_prof.report("cppDE reverse sweep");
    store.prof().report("cppDE reverse checkpoints");
  }
};

}  // namespace reverse
}  // namespace cppde

#endif  // CPPDE_REVERSE_TRAJECTORY_HPP
