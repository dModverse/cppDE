/*
 What one integration step leaves behind, and the matrix its adjoint solves
 against.

 A step is the map (x, t, h, theta) -> (x_out, xerr), where x, t and h are what
 the step reads that an earlier step produced. The forward run stores them per
 accepted step; the adjoint of that map is written in cppde_adjoint_step.hpp.

 The controller is not differentiated: the size is read off the checkpoint and
 the grid is a constant (Bock's internal numerical differentiation, see
 vignette("Methods"), "Internal numerical differentiation").

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
//  Declared, never defined: every stepper family supplies its own. The contract:
//
//    std::size_t n() const     states in the step-start state
//    double t, dt              the two the controller chose
//    const T* start_state()    the state the step begins from, flat
//
//    capture(const Stepper&, const std::vector<Value>& x, double t, double dt)
//        Reads the live forward stepper. Values only, whatever it integrated in.
//
//    apply_tail(St&, RSys&) const
//        Whatever the controller does to the carry between this step's
//        acceptance and the next one's, replayed from the record the forward run
//        left rather than decided again. The system comes with it because one of
//        those operations, the order-1 restart, reseeds the history from the
//        right-hand side. A one-step method has no tail: its carry is the state.
//
//  The multistepper's carry is the Nordsieck array and not the state alone,
//  which is what makes the two families differ here at all.
// ============================================================================

template<class Stepper, class T> struct step_checkpoint;

// Whether the carry has to be read before the step rather than after it, which
// is the multistepper, whose history the step mutates in place.
template<class S, class = void> struct has_step_snapshot : std::false_type {};
template<class S>
struct has_step_snapshot<S, std::void_t<decltype(std::declval<S&>().set_step_snapshot(
    std::declval<typename S::step_snapshot>()))>
> : std::true_type {};

// Whether a stepper scales its embedded estimate into a local error. The
// multistepper's xerr is the raw correction and the error is that times tq[2],
// so steps of different order are otherwise not comparable.
template<class S, class = void>
struct has_error_constant : std::false_type {};

template<class S>
struct has_error_constant<S, std::void_t<decltype(std::declval<const S&>().error_constant())>>
: std::true_type {};

// ----------------------------------------------------------------------------
//  tsit5: an explicit one-step method keeps nothing across a step boundary.
//
//  FSAL is an optimisation, not a dependence: the recycled k1 is f(x, t) at the
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

  // The state the step starts from, which for a one-step method is the whole
  // carry. Flat so both families answer the same question the same way.
  const T* start_state() const { return x.data(); }

  template<class Value>
  void capture(const stepper_type& /*st*/, const std::vector<Value>& x_in,
               double t_in, double dt_in)
  {
    x.resize(x_in.size());
    for (std::size_t i = 0; i < x_in.size(); ++i)
      x[i] = ad_traits::store_as<T>(x_in[i]);
    t  = t_in;
    dt = dt_in;
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

  // What the controller did to the history between this step's acceptance and
  // the next one's: its own tail, and every attempt the next step threw away.
  // A thrown-away attempt is not free, it rescales the history the accepted one
  // is then entered with, so without this record the reverse chain hands the
  // next step a carry at the wrong scale: right value, wrong gradient.
  // Not recorded means no run filled it and the two fields above describe the
  // tail on their own, which is how the step-level tests build a checkpoint.
  history_log ops;
  bool        ops_recorded = false;

  std::size_t n() const { return n_states; }

  // Nordsieck slot 0, which is the state the step starts from.
  const T* start_state() const { return zn.data(); }

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
//  The matrix an implicit method's equations are linearised against, factorised
//  once per reverse step and used in both directions: forward to recover a value
//  a stage solved for, transposed to pass a cotangent back through it.
//
//  Fresh rather than the forward run's own, which belongs to its iteration and is
//  stale by design; MSBP and MSBJ are exactly that staleness. res_scale is what
//  the residual's derivative in its solution is as a multiple of W: the step size
//  for a corrector written in Nordsieck form, one for a Rosenbrock stage.
// ============================================================================

template<class JacFunc, class T = double, bool Sparse = false>
class equation_solver {
  // A transposed solve in double may go through the factorisation of a step
  // the sweep has already passed and be refined against this step's matrix;
  // any other use factorises this step's matrix.
  static constexpr bool refinable = std::is_same<T, double>::value;

public:
  explicit equation_solver(JacFunc& jac) : m_jac(&jac) {}

  void prepare(const std::vector<T>& x, T t, T inv_gamma_dt, T res_scale = T(1)) {
    m_lu.resize(x);
    m_n = x.size();
    { auto _tp = m_prof.timer(cppde::prof_cat::jac_eval);
      m_lu.call_jacobian(*m_jac, const_cast<std::vector<T>&>(x), t); }
    m_c = inv_gamma_dt;
    m_scale = res_scale;
    m_pending = true;
    if constexpr (!refinable) factorise();
  }

  void forward(std::vector<T>& b) {
    if (m_pending) factorise();
    auto _tp = m_prof.timer(cppde::prof_cat::rev_solve);
    m_lu.solve(b);
  }

  void transposed(std::vector<T>& b) {
    bool done = false;
    if constexpr (refinable) {
      if (m_pending && m_held) done = refine(b);
    }
    if (!done) {
      if (m_pending) factorise();
      auto _tp = m_prof.timer(cppde::prof_cat::rev_solve);
      m_lu.solve_transposed(b);
    }
    if (m_scale != T(1)) for (T& v : b) v /= m_scale;
  }

  // Per-category timings to stderr: the per-step Jacobian and factorisation
  // against the transposed solves, which is what a step's linear algebra is.
  // Compiled away without CPPDE_PROFILE.
  void report_profile() const {
    m_prof.report("cppDE reverse linear algebra");
#ifdef CPPDE_PROFILE
    std::fprintf(stderr, "  refinement: ratio %.2f, budget %d, %d tried, %d accepted, %d passes\n",
                 m_ratio_seen, m_budget, m_n_try, m_n_ok, m_n_pass);
#endif
  }

private:
  void factorise() {
    { auto _tp = m_prof.timer(cppde::prof_cat::w_build);
      m_lu.build_W(m_n, m_c); }
    { auto _tp = m_prof.timer(cppde::prof_cat::lu_factor);
      m_lu.factorize_built_W(); }
    m_pending = false;
    m_held = true;
  }

  // W_k^T x = b by refinement on the held factorisation of a swept step's W_h,
  // x <- x + W_h^{-T} (b - W_k^T x), until the correction is at rounding level.
  // False, b untouched, when it contracts too slowly to beat a factorisation.
  bool refine(std::vector<T>& b) {
    const int budget = iteration_budget();
    m_budget = budget;
    if (budget < min_budget) return false;
    // After a refinement that did not pay, the next few steps factorise
    // directly, twice as many after each further miss.
    if (m_skip > 0) { --m_skip; return false; }
    ++m_n_try;
    auto _tp = m_prof.timer(cppde::prof_cat::rev_solve);
    const std::size_t n = b.size();
    m_b = b;
    m_x = b;
    m_lu.solve_transposed(m_x);
    double prev = 0.0;
    for (int it = 0; it < budget; ++it) {
      m_lu.negJ_transposed_apply(m_x, m_r);
      double xmax = 0.0;
      for (std::size_t i = 0; i < n; ++i) {
        m_r[i] = m_b[i] - (m_c * m_x[i] + m_r[i]);
        xmax = std::max(xmax, std::abs(m_x[i]));
      }
      m_lu.solve_transposed(m_r);
      ++m_n_pass;
      double dmax = 0.0;
      for (std::size_t i = 0; i < n; ++i) {
        m_x[i] += m_r[i];
        dmax = std::max(dmax, std::abs(m_r[i]));
      }
      if (dmax <= 4.0 * eps * xmax ||
          (it > 0 && dmax <= stall_tol * xmax && dmax > stall_rate * prev)) {
        b.swap(m_x);
        ++m_n_ok;
        m_backoff = 0;
        return true;
      }
      if (it > 0 && dmax > rate_max * prev) break;
      prev = dmax;
    }
    m_backoff = std::min(2 * m_backoff + 1, max_backoff);
    m_skip = m_backoff;
    return false;
  }

  // As many passes as one factorisation costs in flops, a pass being a pair of
  // triangular solves and a product with J. Fixed by the factors, so a sweep
  // is reproducible.
  int iteration_budget() {
    m_ratio_seen = m_lu.refactor_solve_ratio();
    const double k = m_ratio_seen / pass_cost - 1.0;
    return static_cast<int>(std::min(static_cast<double>(max_iter), std::floor(k)));
  }

  static constexpr double eps = std::numeric_limits<double>::epsilon();
  static constexpr double rate_max = 0.5;
  // A correction this small that no longer halves has reached its floor.
  static constexpr double stall_tol = 1e-12;
  static constexpr double stall_rate = 0.5;
  // One pass, the solve plus the product with J, in units of one solve.
  static constexpr double pass_cost = 1.5;
  static constexpr int max_iter = 12;
  static constexpr int max_backoff = 31;
  static constexpr int min_budget = 3;

  JacFunc*             m_jac;
  cppde::lu_W<T, Sparse> m_lu;
  T                    m_scale = T(1);
  T                    m_c = T(0);
  bool                 m_pending = false;
  bool                 m_held = false;
  std::size_t          m_n = 0;
  std::vector<double>  m_b, m_x, m_r;
  int                  m_budget = 0, m_n_try = 0, m_n_ok = 0, m_n_pass = 0;
  int                  m_skip = 0, m_backoff = 0;
  double               m_ratio_seen = 0.0;
  cppde::profiler      m_prof;
};

}  // namespace reverse
}  // namespace cppde

#endif  // CPPDE_REVERSE_STEP_HPP
