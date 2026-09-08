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
#include <vector>

#include <cppde/cppde_ad_traits.hpp>
#include <cppde/cppde_codual.hpp>
#include <cppde/cppde_codual_math.hpp>
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

template<class Value, class Resizer, class T>
struct step_checkpoint<cppde::tsit5<Value, Resizer>, T> {
  using scalar_type  = T;
  using stepper_type = cppde::tsit5<Value, Resizer>;

  std::vector<T> x;         // step-start state
  double         t  = 0.0;
  double         dt = 0.0;

  std::size_t n() const { return x.size(); }

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
//    seed(w) / seed(i, w)   cotangent of the step end; repeated seeds add
//    seed_err(w)            cotangent of the embedded error estimate, the way
//                           the control law reaches back into the step
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

  tape_type& tape() const { return codual_tape_for<T>(); }

  // Drops the previous step's nodes and adjoints. Parameter cotangents are read
  // out per step, so nothing has to survive this.
  void begin() { tape().rewind(); }

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
    cp.load(m_stepper, m_x, m_history);
    m_t  = rev_type(static_cast<T>(cp.t));
    m_dt = rev_type(static_cast<T>(cp.dt));
    if (m_tape_stepsize) { m_t.independent(); m_dt.independent(); }
    const std::size_t n = m_x.size();
    m_xout.assign(n, rev_type());
    m_xerr.assign(n, rev_type());
    m_stepper.do_step(sys, m_x, m_t, m_xout, m_dt, m_xerr);
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
    tape().reverse();
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
  rev_stepper&                 stepper()    { return m_stepper; }

private:
  rev_stepper            m_stepper;
  std::vector<rev_type>  m_x, m_history, m_xout, m_xerr;
  rev_type               m_t, m_dt;
  std::vector<T>         m_wx, m_whistory;
  T                      m_wt{}, m_wdt{};
  bool                   m_tape_stepsize = true;
};

}  // namespace reverse
}  // namespace cppde

#endif  // CPPDE_REVERSE_STEP_HPP
