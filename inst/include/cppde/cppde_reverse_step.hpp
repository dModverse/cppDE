/*
 What one integration step leaves behind, and the matrix the step adjoint solves
 against.

 A step is the map (x, t, h, theta) -> (x_out, xerr), where x, t and h are what
 the step reads that an earlier step produced. The forward run stores them per
 accepted step; the adjoint of that map is written in cppde_adjoint_step.hpp.

 The controller is not differentiated: the size is read off the checkpoint and
 the grid is a constant (Bock's internal numerical differentiation, see
 "Internal numerical differentiation" in
 https://dmodverse.github.io/cppDE/articles/methods-reverse-checks.html).

 Acceptance, order, iteration and rebuild counts are piecewise constant and stay
 control decisions; what is smooth within one control path is differentiated.

 Copyright (C) 2026 Simon Beyer
 */

#ifndef CPPDE_REVERSE_STEP_HPP
#define CPPDE_REVERSE_STEP_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <limits>
#include <type_traits>
#include <vector>

#include <cppde/cppde_ad_traits.hpp>
#include <cppde/cppde_lu.hpp>
#include <cppde/cppde_multistepper.hpp>
#include <cppde/cppde_profiler.hpp>
#include <cppde/cppde_rosenbrock4.hpp>
#include <cppde/cppde_tsit5.hpp>

namespace cppde {
namespace reverse {

// ============================================================================
//  step_checkpoint<Stepper, T>
//
//  Declared, never defined: every stepper family supplies a specialization. The contract:
//
//    std::size_t n() const     states in the step-start state
//    double t, dt              the two the controller chose
//    const T* start_state()    the state the step begins from, flat
//
//    capture(const Stepper&, const std::vector<Value>& x, double t, double dt)
//        Reads the live forward stepper. Values only, whatever it integrated in.
//
//    apply_tail(St&, RSys&) const
//        Whatever the controller does to the state passed between this step's
//        acceptance and the next one's, replayed from the record the forward run
//        left rather than decided again. The system comes with it because one of
//        those operations, the order-1 restart, reseeds the history from the
//        right-hand side. A one-step method has no tail: it passes on the state.
//
//  The multistepper passes on the Nordsieck carry and not the state alone,
//  which is what makes the two families differ here at all.
// ============================================================================

template<class Stepper, class T> struct step_checkpoint;

// Whether the state passed on has to be read before the step rather than after it, which
// is the multistepper, whose history the step mutates in place.
template<class S, class = void> struct has_step_snapshot : std::false_type {};
template<class S>
struct has_step_snapshot<S, std::void_t<decltype(std::declval<S&>().set_step_snapshot(
    std::declval<typename S::step_snapshot>()))>
> : std::true_type {};

// ----------------------------------------------------------------------------
//  tsit5: an explicit one-step method keeps nothing across a step boundary.
//
//  FSAL is an optimization, not a dependence: the recycled k1 is f(x, t) at the
//  checkpointed x, so the replay recomputes it bit for bit and k7 is not stored.
// ----------------------------------------------------------------------------

template<class Stepper, class T>
struct onestep_checkpoint {
  using scalar_type  = T;
  using stepper_type = Stepper;

  std::vector<T> x;         // step-start state
  double         t  = 0.0;
  double         dt = 0.0;

  std::size_t n() const { return x.size(); }

  // The state the step starts from, which for a one-step method is all a step
  // passes on. Flat so both families answer the same question the same way.
  const T* start_state() const { return x.data(); }

  template<class Value>
  void capture(const stepper_type& /*st*/, const std::vector<Value>& x_in,
               double t_in, double dt_in, int /*slots*/ = -1)
  {
    x.resize(x_in.size());
    for (std::size_t i = 0; i < x_in.size(); ++i)
      x[i] = ad_traits::store_as<T>(x_in[i]);
    t  = t_in;
    dt = dt_in;
  }

};

template<class Value, class T>
struct step_checkpoint<cppde::tsit5<Value>, T>
  : onestep_checkpoint<cppde::tsit5<Value>, T> {};

// ----------------------------------------------------------------------------
//  rosenbrock4: a one-step method too, so the same checkpoint. The stages are linear
//  solves rather than explicit combinations, and the replay recovers each stage
//  value from the same factorization it sweeps with, so none of them is stored.
// ----------------------------------------------------------------------------

template<class Value, class JacobianPattern, class T>
struct step_checkpoint<cppde::rosenbrock4<Value, JacobianPattern>, T>
  : onestep_checkpoint<cppde::rosenbrock4<Value, JacobianPattern>, T>
{
  // The rb4 interpolant is built from the stages rather than read off them, so the
  // snapshot has to be taken before anything is interpolated. tsit5 needs no
  // such call, and must not get one: it would arm the FSAL recycle, and the next
  // replayed step starts from the checkpoint of that step rather than from this one.
};

// ----------------------------------------------------------------------------
//  multistepper: the Nordsieck carry passes from step to step, not the state alone.
//
//  What the forward run has to leave behind is that history, the scalars that
//  pin the meaning of the history, the state the corrector ended on and the
//  iterations it took to get there. The tail is the controller's: the rank-1
//  Nordsieck update, the order it chose and the rescale it applied, all replayed
//  from the decisions rather than recomputed, since those are control decisions.
// ----------------------------------------------------------------------------

template<multistep_method Method, class Value, class JacobianPattern, class T>
struct step_checkpoint<cppde::multistepper<Method, Value, JacobianPattern>, T> {
  using scalar_type  = T;
  using stepper_type = cppde::multistepper<Method, Value, JacobianPattern>;
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

  // What the controller did to the history between this step's acceptance and
  // the next one's: the tail of this step, and every attempt the next step threw away.
  // A thrown-away attempt is not free, it rescales the history the accepted one
  // is then entered with, so without this record the reverse chain hands the
  // next step a Nordsieck carry at the wrong scale: right value, wrong gradient.
  // Not recorded means no run filled it and the two fields above describe the
  // tail on their own, which is how the step-level tests build a checkpoint.
  history_log ops;
  bool        ops_recorded = false;

  // The corrector iterations the accepted attempt took: the k_iter iterates f
  // was evaluated at (slot-major, the prediction first), the scale of each
  // solve, and the index of the iteration matrix in the store.
  std::vector<T> iter;
  int            k_iter = 0;
  double         s_iter = 1.0;
  std::size_t    mat    = 0;

  std::size_t n() const { return n_states; }

  // Nordsieck slot 0, which is the state the step starts from; empty where the
  // store kept no history for this step.
  const T* start_state() const { return zn.data(); }

  // `slots` of the history kept, all of them where negative.
  void capture(const stepper_type& st, const std::vector<Value>& x_in,
               double t_in, double dt_in, int slots = -1)
  {
    st.save_carry(carry);
    n_states = x_in.size();
    const int keep = slots < 0 ? carry.q + 1 : std::min(slots, carry.q + 1);
    zn.assign(static_cast<std::size_t>(keep) * n_states, T());
    for (int j = 0; j < keep; ++j) {
      const auto& slot = st.zn(j);
      for (std::size_t i = 0; i < n_states; ++i)
        zn[static_cast<std::size_t>(j) * n_states + i] =
            ad_traits::store_as<T>(slot[i]);
    }
    t  = t_in;
    dt = dt_in;
  }


  // The tail on a stepper of any value type: the reverse replay and a forward
  // reference over the recorded sequence both need it. Everything in it is a
  // decision the forward run took, replayed rather than taken again.
  template<class St, class RSys>
  void apply_tail(St& st, RSys& sys) const
  {
    st.complete_step();
    // The Nordsieck interpolant is anchored at tn_current, which the controller
    // sets here and the replay has to as well.
    st.set_tn_current(t + dt);
    st.prepare_dense_output();
    if (ops_recorded) {
      for (const history_entry& e : ops) {
        switch (e.op) {
          case history_op::rescale:
            st.rescale(static_cast<typename St::time_type>(e.value)); break;
          case history_op::order:
            st.set_order_for_next_step(static_cast<int>(e.value)); break;
          // An order increase reads the top slot, which the controller fills
          // with the accumulated correction first.
          case history_op::save_acor:  st.save_acor_to_zn_qmax(); break;
          case history_op::hscale:     st.set_hscale(e.value); break;
          case history_op::reload_zn1: st.reload_zn1_from_f(sys.first); break;
          case history_op::complete:   break;   // the cut, applied above
        }
      }
    } else {
      if (q_next > carry.q) st.save_acor_to_zn_qmax();
      if (q_next != carry.q) st.set_order_for_next_step(q_next);
      if (std::abs(eta - 1.0) > 1e-14)
        st.rescale(static_cast<typename St::time_type>(eta));
    }
  }

};

// ============================================================================
//  equation_solver
//
//  The iteration matrix of a multistep corrector, (1/gamma_p) I - J(x, t), as
//  the forward run factorized it: rebuilt from the point, time and gamma_p the
//  run recorded, on values alone, and factorized once for every step that
//  iterated against it. The sweep solves with the transpose.
// ============================================================================

template<class JacFunc, class T = double, bool Sparse = false>
class equation_solver {
public:
  explicit equation_solver(JacFunc& jac) : m_jac(&jac) { m_lu.set_frozen(true); }

  // The matrix under `key`; a no-op when that key is the one held.
  void prepare_frozen(const std::vector<double>& x, double t, double inv_gamma_p,
                      std::size_t key) {
    if (m_held && key == m_key) return;
    const std::size_t n = x.size();
    m_x.assign(x.begin(), x.end());
    m_lu.resize(m_x);
    { auto _tp = m_prof.timer(cppde::prof_cat::jac_eval);
      m_lu.call_jacobian(*m_jac, m_x, T(t)); }
    { auto _tp = m_prof.timer(cppde::prof_cat::w_build);
      m_lu.build_W(n, T(inv_gamma_p)); }
    { auto _tp = m_prof.timer(cppde::prof_cat::lu_factor);
      m_lu.factorize_built_W(); }
    m_held = true;
    m_key = key;
  }

  // Drops the held matrix, so the next prepare_frozen() builds the matrix anew.
  void forget() { m_held = false; }

  void transposed(std::vector<T>& b) {
    auto _tp = m_prof.timer(cppde::prof_cat::rev_solve);
    m_lu.solve_transposed(b);
  }

  // Per-category timings to stderr. Compiled away without CPPDE_PROFILE.
  void report_profile() const { m_prof.report("cppDE reverse linear algebra"); }

private:
  JacFunc*               m_jac;
  cppde::lu_W<T, Sparse> m_lu;
  std::vector<T>         m_x;
  bool                   m_held = false;
  std::size_t            m_key = 0;
  cppde::profiler        m_prof;
};

}  // namespace reverse
}  // namespace cppde

#endif  // CPPDE_REVERSE_STEP_HPP
