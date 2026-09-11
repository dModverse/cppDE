/*
 The trajectory backwards: checkpoint store, reverse loop, seeds at the
 observation times.

 The forward run integrates in double and drops one checkpoint per accepted
 step, plus a note of which observation fell where and of every intervention.
 That is what this header holds. The backward walk over it is written rather
 than recorded and lives in cppde_adjoint_step.hpp.

 What comes out there is the cotangent of the trajectory start and, summed over
 every step, the cotangent of the parameters. The parameter accumulator has no
 state dimension: it is the quadrature the continuous adjoint writes as an
 integral.

 Observations sit at interpolated times, not at step ends, so a step's own
 continuous extension carries them. An observation before the first step reaches
 the initial state directly.

 The step grid is the forward run's and is not differentiated. A step reads one
 thing from the one before it, the state; the size it took is a constant read
 off its checkpoint. This is Bock's internal numerical differentiation: the
 nominal run adapts freely, and the derivative is taken of the scheme that run
 actually applied.

 Differentiating the controller instead puts spurious derivatives of the time
 steps into the chain rule, which make the discrete adjoint inconsistent with
 the adjoint ODE. cppDE did that behind a switch until 2026-09-09; the term was
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
      // The grid is not differentiated, so the times a checkpoint keeps are
      // values whatever the run integrates in.
      m_store.capture(ctl.stepper(), m_st.previous_state(),
                      ad_traits::scalar_value(m_st.previous_time()),
                      ad_traits::scalar_value(ctl.dt_old()));
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
      out[i] = ad_traits::store_as<T>(in[i]);
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

}  // namespace reverse
}  // namespace cppde

#endif  // CPPDE_REVERSE_TRAJECTORY_HPP
