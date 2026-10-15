/**
 * @file cppde_sliding.hpp
 * @brief State switches that slide along their surfaces (Filippov).
 *
 * Where the branches of the switches in a set S push the solution onto their
 * surfaces g_k = 0, the solution stays on them and f is the mix
 * sum_c w_c(a) f_c over the combinations c of their modes, w_c the product of
 * a_k or 1 - a_k, with the weights a chosen so that r'_k = dg_k/dt +
 * grad g_k . f is zero for every k in S. For one switch
 * a = r'_off / (r'_off - r'_on). A slide ends where its weight reaches 0 or 1,
 * a root of its own, and the solution leaves to that side; f is continuous
 * there. a is formed in the model's scalar type, so dual sensitivities take
 * its dependence on x and the parameters.
 *
 * Copyright (C) 2026 Simon Beyer
 */

#ifndef CPPDE_SLIDING_HPP
#define CPPDE_SLIDING_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include <cppde/cppde_events.hpp>
#include <cppde/cppde_saltation.hpp>
#include <cppde/cppde_switch_modes.hpp>
#include <cppde/cppde_switch_surface.hpp>
#include <cppde/cppde_types.hpp>

namespace cppde {
namespace detail {

// The model's right-hand side with the switch roots it slides along.
template<class System, class State, class Time>
struct sliding_core {
  using V = typename State::value_type;
  using Roots = std::vector<RootEvent<State, Time>>;

  explicit sliding_core(const typename System::first_type& f) : rhs(f) {}

  typename System::first_type rhs;
  Roots roots;                 // the model's roots, then two exits per switch
  std::vector<int> root_of;    // the root of each switch
  std::vector<State> fm;       // f of each combination of the slides
  State f_tmp;

  // Whether the last Jacobian wrote entries outside the model's pattern, and
  // whether a slide met a Jacobian it could not write exactly.
  bool wrote_full = false;
  bool inexact = false;

  void eval(const State& x, State& f, const Time& t) {
    f.assign(x.size(), V(0.0));
    for (auto& v : f) ad_traits::arm_tangents(v);
    rhs(x, f, t);
  }

  // f of every combination of the slides S, into fm; the modes are restored.
  void combos(const std::vector<int>& S, const State& x, const Time& t) {
    signed char* md = switch_mode_sink().m;
    const std::size_t M = std::size_t(1) << S.size();
    fm.resize(M);
    for (std::size_t c = 0; c < M; ++c) {
      for (std::size_t i = 0; i < S.size(); ++i) md[S[i]] = ((c >> i) & 1u) ? 1 : 0;
      eval(x, fm[c], t);
    }
    for (int k : S) md[k] = sliding_mode;
  }

  // The weights a of the slides S at (x, t): the mix of the combinations in fm
  // under which r' of every root in S is zero. Newton in double, then in the
  // model's scalar so a takes its derivatives; false where it does not
  // converge. A weight may leave [0, 1] by a step before its exit is located.
  // `rr` receives r' of root k under combination c.
  bool solve(const std::vector<int>& S, const State& x, const Time& t, std::vector<V>& a,
             std::vector<std::vector<V>>* rr = nullptr) {
    const std::size_t n = S.size(), M = std::size_t(1) << n;
    combos(S, x, t);
    std::vector<std::vector<V>> r(n, std::vector<V>(M));
    double terms = 0.0;
    for (std::size_t k = 0; k < n; ++k)
      for (std::size_t c = 0; c < M; ++c)
        r[k][c] = root_rate(roots[root_of[S[k]]], x, t, fm[c], terms);
    auto newton = [&](auto& w) {
      using W = std::decay_t<decltype(w[0])>;
      std::vector<W> R(n, W(0.0)), A(n * n, W(0.0));
      for (std::size_t k = 0; k < n; ++k)
        for (std::size_t c = 0; c < M; ++c) {
          const W rk = W(r[k][c]);
          R[k] = R[k] + rk * combo_weight(w, c);
          for (std::size_t j = 0; j < n; ++j)
            A[k * n + j] = A[k * n + j] + rk * combo_weight(w, c, static_cast<int>(j));
        }
      small_solve(A, R, [](const auto& v) { return cppde::detail::scalar_value(v); });
      double step = 0.0;
      for (std::size_t j = 0; j < n; ++j) {
        w[j] = w[j] - R[j];
        step = std::max(step, std::abs(scalar_value(R[j])));
      }
      return step;
    };
    std::vector<double> ad(n, 0.5);
    double step = 1.0;
    for (int it = 0; it < 60 && step > 1e-15; ++it) step = newton_double(r, ad);
    for (double v : ad)
      if (!std::isfinite(v)) return false;
    if (step > 1e-10) return false;
    a.assign(ad.begin(), ad.end());
    for (int it = 0; it < 3; ++it) newton(a);
    if (rr) *rr = r;
    return true;
  }

  // One Newton step of the weights in double.
  static double newton_double(const std::vector<std::vector<V>>& r, std::vector<double>& w) {
    const std::size_t n = w.size(), M = std::size_t(1) << n;
    std::vector<double> R(n, 0.0), A(n * n, 0.0);
    for (std::size_t k = 0; k < n; ++k)
      for (std::size_t c = 0; c < M; ++c) {
        const double rk = scalar_value(r[k][c]);
        R[k] += rk * combo_weight(w, c);
        for (std::size_t j = 0; j < n; ++j)
          A[k * n + j] += rk * combo_weight(w, c, static_cast<int>(j));
      }
    small_solve(A, R, [](const auto& v) { return cppde::detail::scalar_value(v); });
    double step = 0.0;
    for (std::size_t j = 0; j < n; ++j) {
      w[j] -= R[j];
      step = std::max(step, std::abs(R[j]));
    }
    return step;
  }

  // f = sum_c w_c(a) fm[c].
  void mix(const std::vector<V>& a, State& f) const {
    for (std::size_t i = 0; i < f.size(); ++i) {
      V s(0.0);
      for (std::size_t c = 0; c < fm.size(); ++c) s = s + combo_weight(a, c) * fm[c][i];
      f[i] = s;
    }
  }

  static void no_mix() {
    throw std::runtime_error(
        "the state switches that slide at once have no Filippov mix of their "
        "branches that keeps the solution on all their surfaces");
  }

  // f where the switches in S slide.
  void slide(const std::vector<int>& S, const State& x, State& dxdt, const Time& t) {
    std::vector<V> a;
    if (!solve(S, x, t, a)) no_mix();
    mix(a, dxdt);
  }

  // r' of switch k with its mode held at `on`, the other slides mixed, and the
  // size of its terms.
  V rate(int k, bool on, const State& x, const Time& t, double& terms) {
    signed char& m = switch_mode_sink().m[k];
    const signed char keep = m;
    m = on ? 1 : 0;
    const std::vector<int> S = sliding_switches();
    if (S.empty()) {
      eval(x, f_tmp, t);
    } else {
      std::vector<V> a;
      if (!solve(S, x, t, a)) {
        m = keep;
        no_mix();
      }
      f_tmp.assign(x.size(), V(0.0));
      mix(a, f_tmp);
    }
    m = keep;
    return root_rate(roots[root_of[k]], x, t, f_tmp, terms);
  }

  // The weight of slide k at (x, t); k slides.
  V weight(int k, const State& x, const Time& t) {
    const std::vector<int> S = sliding_switches();
    std::vector<V> a;
    if (!solve(S, x, t, a)) no_mix();
    for (std::size_t i = 0; i < S.size(); ++i)
      if (S[i] == k) return a[i];
    return V(0.0);
  }
};

// The pair (right-hand side, Jacobian) the steppers take, sliding where
// switches slide. A dense Jacobian is that of the slide,
// sum_c w_c J_c + sum_k df/da_k grad a_k', grad a by central differences,
// since the reverse mode reads it as exact; a sparse one keeps the combination
// with the largest weight, enough for the Newton iteration.
template<class System, class State, class Time>
struct sliding_system {
  using core_type = sliding_core<System, State, Time>;
  using V = typename State::value_type;

  struct first_type {
    std::shared_ptr<core_type> c;
    void operator()(const State& x, State& dxdt, const Time& t) const {
      const std::vector<int> S = sliding_switches();
      if (S.empty()) {
        c->rhs(x, dxdt, t);
        return;
      }
      c->slide(S, x, dxdt, t);
    }
  };

  struct second_type {
    typename System::second_type jac;
    std::shared_ptr<core_type> c;

    template<class Matrix>
    void operator()(const State& x, Matrix& J, const Time& t, State& dfdt) {
      const std::vector<int> S = sliding_switches();
      if (S.empty()) {
        if (c->wrote_full) {
          clear(J);
          c->wrote_full = false;
        }
        jac(x, J, t, dfdt);
        return;
      }
      std::vector<V> a;
      if (!c->solve(S, x, t, a)) core_type::no_mix();
      const std::size_t M = c->fm.size(), n = x.size(), ns = S.size();
      signed char* md = switch_mode_sink().m;
      if constexpr (std::is_same_v<Matrix, dense_matrix<V>>) {
        // df/da_k from the combinations at x, before the differences overwrite them.
        std::vector<State> dfda(ns, State(n));
        for (std::size_t k = 0; k < ns; ++k)
          for (std::size_t i = 0; i < n; ++i) {
            V s(0.0);
            for (std::size_t cc = 0; cc < M; ++cc)
              s = s + combo_weight(a, cc, static_cast<int>(k)) * c->fm[cc][i];
            dfda[k][i] = s;
          }
        for (int i = 0; i < static_cast<int>(n); ++i)
          for (int j = 0; j < static_cast<int>(n); ++j) J(i, j) = V(0.0);
        for (std::size_t i = 0; i < n; ++i) dfdt[i] = V(0.0);
        Matrix Jc(static_cast<int>(n), static_cast<int>(n));
        State dtc(n);
        for (std::size_t cc = 0; cc < M; ++cc) {
          for (std::size_t q = 0; q < ns; ++q) md[S[q]] = ((cc >> q) & 1u) ? 1 : 0;
          jac(x, Jc, t, dtc);
          const V w = combo_weight(a, cc);
          for (int i = 0; i < static_cast<int>(n); ++i) {
            for (int j = 0; j < static_cast<int>(n); ++j) J(i, j) = J(i, j) + w * Jc(i, j);
            dfdt[i] = dfdt[i] + w * dtc[i];
          }
        }
        for (int k : S) md[k] = sliding_mode;
        // grad a and da/dt by central differences of the weights.
        State xh = x;
        std::vector<V> up, dn;
        for (std::size_t j = 0; j < n; ++j) {
          const double h = 1e-6 * std::max(1.0, std::abs(scalar_value(x[j])));
          xh[j] = x[j] + V(h);
          if (!c->solve(S, xh, t, up)) core_type::no_mix();
          xh[j] = x[j] - V(h);
          if (!c->solve(S, xh, t, dn)) core_type::no_mix();
          xh[j] = x[j];
          // J holds -J: subtract df/da_k grad a_k'.
          for (std::size_t k = 0; k < ns; ++k) {
            const V g = (up[k] - dn[k]) / V(2.0 * h);
            for (std::size_t i = 0; i < n; ++i)
              J(static_cast<int>(i), static_cast<int>(j)) =
                  J(static_cast<int>(i), static_cast<int>(j)) - dfda[k][i] * g;
          }
        }
        const double ht = 1e-6 * std::max(1.0, std::abs(scalar_value(t)));
        if (!c->solve(S, x, t + Time(ht), up) || !c->solve(S, x, t - Time(ht), dn))
          core_type::no_mix();
        for (std::size_t k = 0; k < ns; ++k) {
          const V da = (up[k] - dn[k]) / V(2.0 * ht);
          for (std::size_t i = 0; i < n; ++i) dfdt[i] = dfdt[i] + dfda[k][i] * da;
        }
        c->wrote_full = true;
      } else {
        std::size_t best = 0;
        double wbest = -1.0;
        for (std::size_t cc = 0; cc < M; ++cc) {
          const double w = scalar_value(combo_weight(a, cc));
          if (w > wbest) { wbest = w; best = cc; }
        }
        for (std::size_t q = 0; q < ns; ++q) md[S[q]] = ((best >> q) & 1u) ? 1 : 0;
        jac(x, J, t, dfdt);
        for (int k : S) md[k] = sliding_mode;
        c->inexact = true;
      }
    }

    template<class Matrix>
    static void clear(Matrix& J) {
      if constexpr (std::is_same_v<Matrix, dense_matrix<V>>)
        for (int i = 0; i < J.rows(); ++i)
          for (int j = 0; j < J.cols(); ++j) J(i, j) = V(0.0);
    }
  };

  first_type first;
  second_type second;

  const typename core_type::Roots& roots() const { return first.c->roots; }

  // Signs of r' of switch k with its mode off and on, the other slides mixed.
  void rate_signs(int k, const State& x, const Time& t, int& off, int& on) const {
    double terms = 0.0;
    off = rate_sign(scalar_value(first.c->rate(k, false, x, t, terms)), terms);
    on = rate_sign(scalar_value(first.c->rate(k, true, x, t, terms)), terms);
  }

  // Whether the slides as they stand have a mix at (x, t), weights in [0, 1].
  bool mixable(const State& x, const Time& t) const {
    const std::vector<int> S = sliding_switches();
    std::vector<V> a;
    if (S.empty()) return true;
    if (!first.c->solve(S, x, t, a)) return false;
    for (const V& v : a) {
      const double w = scalar_value(v);
      if (!(w >= -1e-9 && w <= 1.0 + 1e-9)) return false;
    }
    return true;
  }
};

template<class S> struct is_sliding_system : std::false_type {};
template<class System, class State, class Time>
struct is_sliding_system<sliding_system<System, State, Time>> : std::true_type {};

// `system` with the switch roots of `root`, and two exit roots per switch
// appended to them, negative while it slides: a_k - 1 and -a_k.
template<class State, class Time, class System>
sliding_system<System, State, Time>
make_sliding_system(const System& system, const std::vector<RootEvent<State, Time>>& root) {
  using S = sliding_system<System, State, Time>;
  using V = typename State::value_type;
  auto c = std::make_shared<typename S::core_type>(system.first);
  c->roots = root;
  for (std::size_t i = 0; i < root.size(); ++i) {
    const int k = root[i].mode;
    if (k < 0) continue;
    if (static_cast<int>(c->root_of.size()) <= k) c->root_of.resize(k + 1, -1);
    c->root_of[k] = static_cast<int>(i);
  }
  typename S::core_type* cp = c.get();
  for (int k = 0; k < static_cast<int>(c->root_of.size()); ++k) {
    for (int on = 1; on >= 0; --on) {
      RootEvent<State, Time> ex;
      auto g = [cp, k, on](const State& x, const Time& t) -> V {
        if (!is_sliding(k)) return V(1.0);
        const V a = cp->weight(k, x, t);
        return on ? a - V(1.0) : V(0.0) - a;
      };
      ex.func = g;
      // The exit time moves with the parameters, which the second order takes
      // through a saltation: central differences of g in x and t, in the
      // model's scalar type so they have their own sensitivities.
      ex.dg_dx = [g](const State& x, const Time& t, State& out) {
        State xh = x;
        for (std::size_t i = 0; i < x.size(); ++i) {
          const double h = 1e-6 * std::max(1.0, std::abs(scalar_value(x[i])));
          xh[i] = x[i] + V(h);
          const V up = g(xh, t);
          xh[i] = x[i] - V(h);
          const V dn = g(xh, t);
          xh[i] = x[i];
          out[i] = (up - dn) / V(2.0 * h);
        }
      };
      ex.dg_dt = [g](const State& x, const Time& t) -> V {
        const double h = 1e-6 * std::max(1.0, std::abs(scalar_value(t)));
        return (g(x, t + Time(h)) - g(x, t - Time(h))) / V(2.0 * h);
      };
      ex.state_index = -1;
      ex.value_func = nullptr;
      ex.method = EventMethod::Replace;
      ex.slide_exit = k;
      ex.exit_on = on == 1;
      c->roots.push_back(std::move(ex));
    }
  }
  return S{typename S::first_type{c}, typename S::second_type{system.second, c}};
}

}  // namespace detail
}  // namespace cppde

#endif  // CPPDE_SLIDING_HPP
