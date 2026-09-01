/*
 Arithmetic operators, math functions, and comparisons for cppde::dual<T, N>,
 plus the scalar (arithmetic-type) overloads of the same math functions.

 All defined in namespace cppde so generated codegen output `cppde::exp(x)`
 etc resolves via ADL or qualified call. Math functions use the `using std::fn`
 idiom to dispatch to std (T=double) or to user-namespace overloads (T=mpfr).

 Convention:
 - Result tangent for unary y = f(x):    y.tan[i] = f'(x.val) * x.tan[i]
 - Result tangent for binary y = f(a,b): y.tan[i] = f_a * a.tan[i] + f_b * b.tan[i]
 - Comparisons fall back to .x().

 Copyright (C) 2026 Simon Beyer
 */

#ifndef CPPDE_DUAL_MATH_HPP
#define CPPDE_DUAL_MATH_HPP

#include <cppde/cppde_dual.hpp>
#include <cppde/cppde_ad_traits.hpp>

#include <cmath>
#include <cstdlib>
#include <type_traits>

namespace cppde {

// =============================================================================
// Scalar (arithmetic-type) overloads.
//
// The codegen writes every math call as cppde::<fn>, so one emitted expression
// compiles for double, dual and dual2nd. The AD overloads are constrained to AD
// operands, which leaves an all-numeric call, such as a folded derivative
// constant, without a candidate. These forward to <cmath>.
// =============================================================================
namespace detail {
template<class... Ts>
inline constexpr bool all_arithmetic_v = (std::is_arithmetic_v<Ts> && ...);
}  // namespace detail

// Return types follow <cmath>: floating-point arguments keep their own type,
// integral ones promote to double.
#define CPPDE_DEFINE_SCALAR_UNARY(NAME)                                       \
  template<class T, std::enable_if_t<std::is_arithmetic_v<T>, int> = 0>       \
  inline auto NAME(T x) { return std::NAME(x); }

CPPDE_DEFINE_SCALAR_UNARY(exp)
CPPDE_DEFINE_SCALAR_UNARY(log)
CPPDE_DEFINE_SCALAR_UNARY(sqrt)
CPPDE_DEFINE_SCALAR_UNARY(sin)
CPPDE_DEFINE_SCALAR_UNARY(cos)
CPPDE_DEFINE_SCALAR_UNARY(tan)
CPPDE_DEFINE_SCALAR_UNARY(asin)
CPPDE_DEFINE_SCALAR_UNARY(acos)
CPPDE_DEFINE_SCALAR_UNARY(atan)
CPPDE_DEFINE_SCALAR_UNARY(sinh)
CPPDE_DEFINE_SCALAR_UNARY(cosh)
CPPDE_DEFINE_SCALAR_UNARY(tanh)
CPPDE_DEFINE_SCALAR_UNARY(asinh)
CPPDE_DEFINE_SCALAR_UNARY(acosh)
CPPDE_DEFINE_SCALAR_UNARY(atanh)

#undef CPPDE_DEFINE_SCALAR_UNARY

// abs: std::abs has no unsigned overload (and its integral one lives in
// <cstdlib>), so route unsigned types around it instead of letting the call
// go ambiguous.
template<class T, std::enable_if_t<std::is_arithmetic_v<T>, int> = 0>
inline auto abs(T x) {
  if constexpr (std::is_unsigned_v<T>) return x;
  else return std::abs(x);
}

template<class A, class B,
         std::enable_if_t<detail::all_arithmetic_v<A, B>, int> = 0>
inline auto pow(A a, B b) { return std::pow(a, b); }

// min / max return by value on the common type: std::min / std::max return a
// reference to a parameter, which would dangle for a mixed-type call that has
// to materialise a converted temporary.
template<class A, class B,
         std::enable_if_t<detail::all_arithmetic_v<A, B>, int> = 0>
inline std::common_type_t<A, B> min(A a, B b) {
  using C = std::common_type_t<A, B>;
  return (static_cast<C>(a) < static_cast<C>(b)) ? static_cast<C>(a)
                                                 : static_cast<C>(b);
}

template<class A, class B,
         std::enable_if_t<detail::all_arithmetic_v<A, B>, int> = 0>
inline std::common_type_t<A, B> max(A a, B b) {
  using C = std::common_type_t<A, B>;
  return (static_cast<C>(a) < static_cast<C>(b)) ? static_cast<C>(b)
                                                 : static_cast<C>(a);
}

// =============================================================================
// Eager against ET routing. The expression overlay covers every dual over a
// non-AD scalar, heap and static-N alike. Eager stays only where T is itself an
// AD type, that is the outer layer of a nested dual2nd.
// =============================================================================
namespace detail {
template<class T, unsigned /*N*/>
struct eager_dual_active
  : std::bool_constant<ad_traits::is_ad<T>::value> {};
} // namespace detail

#define CPPDE_EAGER_GATE(T, N) \
  std::enable_if_t<::cppde::detail::eager_dual_active<T, N>::value, int> = 0

// =============================================================================
// Arithmetic operators (dual op dual)
// =============================================================================

template<class T, unsigned N, CPPDE_EAGER_GATE(T, N)>
inline dual<T, N> operator+(const dual<T, N>& a, const dual<T, N>& b) {
  dual<T, N> r(a.x() + b.x());
  if (!a.depend() && !b.depend()) return r;
  if (a.depend() && b.depend()) {
    r.set_depend_from(a, b);
    for (unsigned i = 0; i < r.loop_size(); ++i) r[i] = a[i] + b[i];
  } else if (a.depend()) {
    r.set_depend_from(a);
    for (unsigned i = 0; i < r.loop_size(); ++i) r[i] = a[i];
  } else {
    r.set_depend_from(b);
    for (unsigned i = 0; i < r.loop_size(); ++i) r[i] = b[i];
  }
  return r;
}

template<class T, unsigned N, CPPDE_EAGER_GATE(T, N)>
inline dual<T, N> operator-(const dual<T, N>& a, const dual<T, N>& b) {
  dual<T, N> r(a.x() - b.x());
  if (!a.depend() && !b.depend()) return r;
  if (a.depend() && b.depend()) {
    r.set_depend_from(a, b);
    for (unsigned i = 0; i < r.loop_size(); ++i) r[i] = a[i] - b[i];
  } else if (a.depend()) {
    r.set_depend_from(a);
    for (unsigned i = 0; i < r.loop_size(); ++i) r[i] = a[i];
  } else {
    r.set_depend_from(b);
    for (unsigned i = 0; i < r.loop_size(); ++i) r[i] = -b[i];
  }
  return r;
}

template<class T, unsigned N, CPPDE_EAGER_GATE(T, N)>
inline dual<T, N> operator*(const dual<T, N>& a, const dual<T, N>& b) {
  dual<T, N> r(a.x() * b.x());
  if (!a.depend() && !b.depend()) return r;
  if (a.depend() && b.depend()) {
    r.set_depend_from(a, b);
    for (unsigned i = 0; i < r.loop_size(); ++i) r[i] = a[i] * b.x() + a.x() * b[i];
  } else if (a.depend()) {
    r.set_depend_from(a);
    for (unsigned i = 0; i < r.loop_size(); ++i) r[i] = a[i] * b.x();
  } else {
    r.set_depend_from(b);
    for (unsigned i = 0; i < r.loop_size(); ++i) r[i] = a.x() * b[i];
  }
  return r;
}

template<class T, unsigned N, CPPDE_EAGER_GATE(T, N)>
inline dual<T, N> operator/(const dual<T, N>& a, const dual<T, N>& b) {
  T inv = T(1) / b.x();
  dual<T, N> r(a.x() * inv);
  if (!a.depend() && !b.depend()) return r;
  if (a.depend() && b.depend()) {
    r.set_depend_from(a, b);
    for (unsigned i = 0; i < r.loop_size(); ++i)
      r[i] = (a[i] - r.x() * b[i]) * inv;
  } else if (a.depend()) {
    r.set_depend_from(a);
    for (unsigned i = 0; i < r.loop_size(); ++i) r[i] = a[i] * inv;
  } else {
    r.set_depend_from(b);
    for (unsigned i = 0; i < r.loop_size(); ++i) r[i] = -r.x() * b[i] * inv;
  }
  return r;
}

// Unary
template<class T, unsigned N, CPPDE_EAGER_GATE(T, N)>
inline dual<T, N> operator+(const dual<T, N>& a) { return a; }

template<class T, unsigned N, CPPDE_EAGER_GATE(T, N)>
inline dual<T, N> operator-(const dual<T, N>& a) {
  dual<T, N> r;
  r.x() = -a.x();
  if (a.depend()) {
    r.set_depend_from(a);
    for (unsigned i = 0; i < r.loop_size(); ++i) r[i] = -a[i];
  }
  return r;
}

// =============================================================================
// Arithmetic operators (dual op scalar U / U op dual)
// =============================================================================

template<class T, unsigned N, class U,
         std::enable_if_t<std::is_arithmetic_v<U>
                          && detail::eager_dual_active<T, N>::value, int> = 0>
inline dual<T, N> operator+(const dual<T, N>& a, const U& b) {
  dual<T, N> r;
  r.x() = a.x() + static_cast<T>(b);
  if (a.depend()) {
    r.set_depend_from(a);
    for (unsigned i = 0; i < r.loop_size(); ++i) r[i] = a[i];
  }
  return r;
}
template<class T, unsigned N, class U,
         std::enable_if_t<std::is_arithmetic_v<U>
                          && detail::eager_dual_active<T, N>::value, int> = 0>
inline dual<T, N> operator+(const U& a, const dual<T, N>& b) { return b + a; }

template<class T, unsigned N, class U,
         std::enable_if_t<std::is_arithmetic_v<U>
                          && detail::eager_dual_active<T, N>::value, int> = 0>
inline dual<T, N> operator-(const dual<T, N>& a, const U& b) {
  dual<T, N> r;
  r.x() = a.x() - static_cast<T>(b);
  if (a.depend()) {
    r.set_depend_from(a);
    for (unsigned i = 0; i < r.loop_size(); ++i) r[i] = a[i];
  }
  return r;
}
template<class T, unsigned N, class U,
         std::enable_if_t<std::is_arithmetic_v<U>
                          && detail::eager_dual_active<T, N>::value, int> = 0>
inline dual<T, N> operator-(const U& a, const dual<T, N>& b) {
  dual<T, N> r;
  r.x() = static_cast<T>(a) - b.x();
  if (b.depend()) {
    r.set_depend_from(b);
    for (unsigned i = 0; i < r.loop_size(); ++i) r[i] = -b[i];
  }
  return r;
}

template<class T, unsigned N, class U,
         std::enable_if_t<std::is_arithmetic_v<U>
                          && detail::eager_dual_active<T, N>::value, int> = 0>
inline dual<T, N> operator*(const dual<T, N>& a, const U& b) {
  dual<T, N> r;
  T tb = static_cast<T>(b);
  r.x() = a.x() * tb;
  if (a.depend()) {
    r.set_depend_from(a);
    for (unsigned i = 0; i < r.loop_size(); ++i) r[i] = a[i] * tb;
  }
  return r;
}
template<class T, unsigned N, class U,
         std::enable_if_t<std::is_arithmetic_v<U>
                          && detail::eager_dual_active<T, N>::value, int> = 0>
inline dual<T, N> operator*(const U& a, const dual<T, N>& b) { return b * a; }

template<class T, unsigned N, class U,
         std::enable_if_t<std::is_arithmetic_v<U>
                          && detail::eager_dual_active<T, N>::value, int> = 0>
inline dual<T, N> operator/(const dual<T, N>& a, const U& b) {
  dual<T, N> r;
  T inv = T(1) / static_cast<T>(b);
  r.x() = a.x() * inv;
  if (a.depend()) {
    r.set_depend_from(a);
    for (unsigned i = 0; i < r.loop_size(); ++i) r[i] = a[i] * inv;
  }
  return r;
}
template<class T, unsigned N, class U,
         std::enable_if_t<std::is_arithmetic_v<U>
                          && detail::eager_dual_active<T, N>::value, int> = 0>
inline dual<T, N> operator/(const U& a, const dual<T, N>& b) {
  dual<T, N> r;
  T inv = T(1) / b.x();
  r.x() = static_cast<T>(a) * inv;
  if (b.depend()) {
    r.set_depend_from(b);
    for (unsigned i = 0; i < r.loop_size(); ++i) r[i] = -r.x() * b[i] * inv;
  }
  return r;
}

// =============================================================================
// Compound assignment (member declarations resolve here)
// =============================================================================

template<class T, unsigned N>
inline dual<T, N>& dual<T, N>::operator+=(const dual<T, N>& o) {
  *this = *this + o; return *this;
}
template<class T, unsigned N>
inline dual<T, N>& dual<T, N>::operator-=(const dual<T, N>& o) {
  *this = *this - o; return *this;
}
template<class T, unsigned N>
inline dual<T, N>& dual<T, N>::operator*=(const dual<T, N>& o) {
  *this = *this * o; return *this;
}
template<class T, unsigned N>
inline dual<T, N>& dual<T, N>::operator/=(const dual<T, N>& o) {
  *this = *this / o; return *this;
}

template<class T>
inline dual<T, 0>& dual<T, 0>::operator+=(const dual<T, 0>& o) {
  *this = *this + o; return *this;
}
template<class T>
inline dual<T, 0>& dual<T, 0>::operator-=(const dual<T, 0>& o) {
  *this = *this - o; return *this;
}
template<class T>
inline dual<T, 0>& dual<T, 0>::operator*=(const dual<T, 0>& o) {
  *this = *this * o; return *this;
}
template<class T>
inline dual<T, 0>& dual<T, 0>::operator/=(const dual<T, 0>& o) {
  *this = *this / o; return *this;
}

// =============================================================================
// Math: unary functions  y = f(x), y.tan[i] = f'(x.val) * x.tan[i]
// =============================================================================

// VAL_EXPR and DERIV_EXPR may invoke any std math function (e.g. asin's
// derivative needs sqrt). Bringing the whole std-math suite into scope keeps
// ADL working for both T=double (resolves to std::) and T=user-namespace
// (resolves to that namespace via ADL).
#define CPPDE_DEFINE_UNARY(NAME, VAL_EXPR, DERIV_EXPR)                       \
  template<class T, unsigned N, CPPDE_EAGER_GATE(T, N)>                      \
  inline dual<T, N> NAME(const dual<T, N>& a) {                               \
    using std::exp;   using std::log;   using std::sqrt;                      \
    using std::sin;   using std::cos;   using std::tan;                       \
    using std::asin;  using std::acos;  using std::atan;                      \
    using std::sinh;  using std::cosh;  using std::tanh;                      \
    using std::asinh; using std::acosh; using std::atanh;                     \
    dual<T, N> r;                                                             \
    const T xv = a.x();                                                       \
    r.x() = (VAL_EXPR);                                                       \
    if (a.depend()) {                                                         \
      const T fp = (DERIV_EXPR);                                              \
      r.set_depend_from(a);                                                   \
      for (unsigned i = 0; i < r.loop_size(); ++i) r[i] = fp * a[i];               \
    }                                                                         \
    return r;                                                                 \
  }

CPPDE_DEFINE_UNARY(exp,  exp(xv),    r.x())
CPPDE_DEFINE_UNARY(log,  log(xv),    T(1) / xv)
CPPDE_DEFINE_UNARY(sqrt, sqrt(xv),   T(1) / (T(2) * r.x()))

// trig
template<class T, unsigned N, CPPDE_EAGER_GATE(T, N)>
inline dual<T, N> sin(const dual<T, N>& a) {
  using std::sin; using std::cos;
  dual<T, N> r;
  const T xv = a.x();
  r.x() = sin(xv);
  if (a.depend()) {
    const T fp = cos(xv);
    r.set_depend_from(a);
    for (unsigned i = 0; i < r.loop_size(); ++i) r[i] = fp * a[i];
  }
  return r;
}
template<class T, unsigned N, CPPDE_EAGER_GATE(T, N)>
inline dual<T, N> cos(const dual<T, N>& a) {
  using std::sin; using std::cos;
  dual<T, N> r;
  const T xv = a.x();
  r.x() = cos(xv);
  if (a.depend()) {
    const T fp = -sin(xv);
    r.set_depend_from(a);
    for (unsigned i = 0; i < r.loop_size(); ++i) r[i] = fp * a[i];
  }
  return r;
}
template<class T, unsigned N, CPPDE_EAGER_GATE(T, N)>
inline dual<T, N> tan(const dual<T, N>& a) {
  using std::tan;
  dual<T, N> r;
  const T xv = a.x();
  const T tv = tan(xv);
  r.x() = tv;
  if (a.depend()) {
    const T fp = T(1) + tv * tv;  // sec^2 = 1 + tan^2
    r.set_depend_from(a);
    for (unsigned i = 0; i < r.loop_size(); ++i) r[i] = fp * a[i];
  }
  return r;
}

CPPDE_DEFINE_UNARY(asin, asin(xv),  T(1) / sqrt(T(1) - xv * xv))
CPPDE_DEFINE_UNARY(acos, acos(xv), -T(1) / sqrt(T(1) - xv * xv))
CPPDE_DEFINE_UNARY(atan, atan(xv),  T(1) / (T(1) + xv * xv))

// hyperbolic
template<class T, unsigned N, CPPDE_EAGER_GATE(T, N)>
inline dual<T, N> sinh(const dual<T, N>& a) {
  using std::sinh; using std::cosh;
  dual<T, N> r;
  const T xv = a.x();
  r.x() = sinh(xv);
  if (a.depend()) {
    const T fp = cosh(xv);
    r.set_depend_from(a);
    for (unsigned i = 0; i < r.loop_size(); ++i) r[i] = fp * a[i];
  }
  return r;
}
template<class T, unsigned N, CPPDE_EAGER_GATE(T, N)>
inline dual<T, N> cosh(const dual<T, N>& a) {
  using std::sinh; using std::cosh;
  dual<T, N> r;
  const T xv = a.x();
  r.x() = cosh(xv);
  if (a.depend()) {
    const T fp = sinh(xv);
    r.set_depend_from(a);
    for (unsigned i = 0; i < r.loop_size(); ++i) r[i] = fp * a[i];
  }
  return r;
}
template<class T, unsigned N, CPPDE_EAGER_GATE(T, N)>
inline dual<T, N> tanh(const dual<T, N>& a) {
  using std::tanh;
  dual<T, N> r;
  const T xv = a.x();
  const T tv = tanh(xv);
  r.x() = tv;
  if (a.depend()) {
    const T fp = T(1) - tv * tv;
    r.set_depend_from(a);
    for (unsigned i = 0; i < r.loop_size(); ++i) r[i] = fp * a[i];
  }
  return r;
}

CPPDE_DEFINE_UNARY(asinh, asinh(xv), T(1) / sqrt(xv * xv + T(1)))
CPPDE_DEFINE_UNARY(acosh, acosh(xv), T(1) / sqrt(xv * xv - T(1)))
CPPDE_DEFINE_UNARY(atanh, atanh(xv), T(1) / (T(1) - xv * xv))

#undef CPPDE_DEFINE_UNARY

// =============================================================================
// abs: piecewise linear, derivative sign(x); at x=0 we return 0
// =============================================================================
template<class T, unsigned N, CPPDE_EAGER_GATE(T, N)>
inline dual<T, N> abs(const dual<T, N>& a) {
  using std::abs;
  dual<T, N> r;
  const T xv = a.x();
  r.x() = abs(xv);
  if (a.depend()) {
    const T fp = (xv > T(0)) ? T(1) : ((xv < T(0)) ? T(-1) : T(0));
    r.set_depend_from(a);
    for (unsigned i = 0; i < r.loop_size(); ++i) r[i] = fp * a[i];
  }
  return r;
}

// =============================================================================
// pow: dual^dual, dual^scalar, scalar^dual
// y = a^b ⇒  dy = b * a^(b-1) * da + log(a) * a^b * db
// =============================================================================
template<class T, unsigned N, CPPDE_EAGER_GATE(T, N)>
inline dual<T, N> pow(const dual<T, N>& a, const dual<T, N>& b) {
  using std::pow; using std::log;
  dual<T, N> r;
  const T av = a.x();
  const T bv = b.x();
  const T y  = pow(av, bv);
  r.x() = y;
  if (!a.depend() && !b.depend()) return r;
  const T fa = (a.depend()) ? bv * pow(av, bv - T(1)) : T(0);
  const T fb = (b.depend()) ? log(av) * y              : T(0);
  if (a.depend() && b.depend()) {
    r.set_depend_from(a, b);
    for (unsigned i = 0; i < r.loop_size(); ++i) r[i] = fa * a[i] + fb * b[i];
  } else if (a.depend()) {
    r.set_depend_from(a);
    for (unsigned i = 0; i < r.loop_size(); ++i) r[i] = fa * a[i];
  } else {
    r.set_depend_from(b);
    for (unsigned i = 0; i < r.loop_size(); ++i) r[i] = fb * b[i];
  }
  return r;
}
template<class T, unsigned N, class U,
         std::enable_if_t<std::is_arithmetic_v<U>
                          && detail::eager_dual_active<T, N>::value, int> = 0>
inline dual<T, N> pow(const dual<T, N>& a, const U& b) {
  using std::pow;
  dual<T, N> r;
  const T av = a.x();
  const T bv = static_cast<T>(b);
  r.x() = pow(av, bv);
  if (a.depend()) {
    const T fa = bv * pow(av, bv - T(1));
    r.set_depend_from(a);
    for (unsigned i = 0; i < r.loop_size(); ++i) r[i] = fa * a[i];
  }
  return r;
}
template<class T, unsigned N, class U,
         std::enable_if_t<std::is_arithmetic_v<U>
                          && detail::eager_dual_active<T, N>::value, int> = 0>
inline dual<T, N> pow(const U& a, const dual<T, N>& b) {
  using std::pow; using std::log;
  dual<T, N> r;
  const T av = static_cast<T>(a);
  const T bv = b.x();
  const T y  = pow(av, bv);
  r.x() = y;
  if (b.depend()) {
    const T fb = log(av) * y;
    r.set_depend_from(b);
    for (unsigned i = 0; i < r.loop_size(); ++i) r[i] = fb * b[i];
  }
  return r;
}

// =============================================================================
// min / max / clamp: select-by-value, derivatives = winner's derivatives
// =============================================================================
template<class T, unsigned N>
inline dual<T, N> min(const dual<T, N>& a, const dual<T, N>& b) {
  return (a.x() < b.x()) ? a : b;
}
template<class T, unsigned N>
inline dual<T, N> max(const dual<T, N>& a, const dual<T, N>& b) {
  return (a.x() < b.x()) ? b : a;
}
template<class T, unsigned N, class U,
         std::enable_if_t<std::is_arithmetic_v<U>, int> = 0>
inline dual<T, N> min(const dual<T, N>& a, const U& b) {
  return (a.x() < static_cast<T>(b)) ? a : dual<T, N>(b);
}
template<class T, unsigned N, class U,
         std::enable_if_t<std::is_arithmetic_v<U>, int> = 0>
inline dual<T, N> min(const U& a, const dual<T, N>& b) { return min(b, a); }
template<class T, unsigned N, class U,
         std::enable_if_t<std::is_arithmetic_v<U>, int> = 0>
inline dual<T, N> max(const dual<T, N>& a, const U& b) {
  return (a.x() < static_cast<T>(b)) ? dual<T, N>(b) : a;
}
template<class T, unsigned N, class U,
         std::enable_if_t<std::is_arithmetic_v<U>, int> = 0>
inline dual<T, N> max(const U& a, const dual<T, N>& b) { return max(b, a); }

template<class T, unsigned N, class L, class H>
inline dual<T, N> clamp(const dual<T, N>& a, const L& lo, const H& hi) {
  return min(max(a, lo), hi);
}

// =============================================================================
// Comparisons (always on .x())
// =============================================================================
template<class T, unsigned N>
inline bool operator==(const dual<T, N>& a, const dual<T, N>& b) { return a.x() == b.x(); }
template<class T, unsigned N>
inline bool operator!=(const dual<T, N>& a, const dual<T, N>& b) { return a.x() != b.x(); }
template<class T, unsigned N>
inline bool operator< (const dual<T, N>& a, const dual<T, N>& b) { return a.x() <  b.x(); }
template<class T, unsigned N>
inline bool operator<=(const dual<T, N>& a, const dual<T, N>& b) { return a.x() <= b.x(); }
template<class T, unsigned N>
inline bool operator> (const dual<T, N>& a, const dual<T, N>& b) { return a.x() >  b.x(); }
template<class T, unsigned N>
inline bool operator>=(const dual<T, N>& a, const dual<T, N>& b) { return a.x() >= b.x(); }

template<class T, unsigned N, class U,
         std::enable_if_t<std::is_arithmetic_v<U>, int> = 0>
inline bool operator==(const dual<T, N>& a, const U& b) { return a.x() == static_cast<T>(b); }
template<class T, unsigned N, class U,
         std::enable_if_t<std::is_arithmetic_v<U>, int> = 0>
inline bool operator==(const U& a, const dual<T, N>& b) { return static_cast<T>(a) == b.x(); }
template<class T, unsigned N, class U,
         std::enable_if_t<std::is_arithmetic_v<U>, int> = 0>
inline bool operator!=(const dual<T, N>& a, const U& b) { return a.x() != static_cast<T>(b); }
template<class T, unsigned N, class U,
         std::enable_if_t<std::is_arithmetic_v<U>, int> = 0>
inline bool operator!=(const U& a, const dual<T, N>& b) { return static_cast<T>(a) != b.x(); }
template<class T, unsigned N, class U,
         std::enable_if_t<std::is_arithmetic_v<U>, int> = 0>
inline bool operator< (const dual<T, N>& a, const U& b) { return a.x() <  static_cast<T>(b); }
template<class T, unsigned N, class U,
         std::enable_if_t<std::is_arithmetic_v<U>, int> = 0>
inline bool operator< (const U& a, const dual<T, N>& b) { return static_cast<T>(a) <  b.x(); }
template<class T, unsigned N, class U,
         std::enable_if_t<std::is_arithmetic_v<U>, int> = 0>
inline bool operator<=(const dual<T, N>& a, const U& b) { return a.x() <= static_cast<T>(b); }
template<class T, unsigned N, class U,
         std::enable_if_t<std::is_arithmetic_v<U>, int> = 0>
inline bool operator<=(const U& a, const dual<T, N>& b) { return static_cast<T>(a) <= b.x(); }
template<class T, unsigned N, class U,
         std::enable_if_t<std::is_arithmetic_v<U>, int> = 0>
inline bool operator> (const dual<T, N>& a, const U& b) { return a.x() >  static_cast<T>(b); }
template<class T, unsigned N, class U,
         std::enable_if_t<std::is_arithmetic_v<U>, int> = 0>
inline bool operator> (const U& a, const dual<T, N>& b) { return static_cast<T>(a) >  b.x(); }
template<class T, unsigned N, class U,
         std::enable_if_t<std::is_arithmetic_v<U>, int> = 0>
inline bool operator>=(const dual<T, N>& a, const U& b) { return a.x() >= static_cast<T>(b); }
template<class T, unsigned N, class U,
         std::enable_if_t<std::is_arithmetic_v<U>, int> = 0>
inline bool operator>=(const U& a, const dual<T, N>& b) { return static_cast<T>(a) >= b.x(); }

// =============================================================================
// max_abs_all_levels(v): recursive absolute-value sweep over a value and all of
// its tangent slots, so a nested dual is covered down to the second-order
// tangents.
// =============================================================================

template<class T>
inline std::enable_if_t<std::is_arithmetic_v<T>, double>
max_abs_all_levels(const T& v) { return std::abs(static_cast<double>(v)); }

template<class T, unsigned N>
inline double max_abs_all_levels(const dual<T, N>& v) {
  auto& vm = const_cast<dual<T, N>&>(v);
  double m = max_abs_all_levels(vm.x());
  unsigned nd = vm.size();
  for (unsigned i = 0; i < nd; ++i)
    m = std::max(m, max_abs_all_levels(vm.d(i)));
  return m;
}

// max_abs_all_levels_vec(v): vector wrapper. Used by integrate_times to
// evaluate the equilibrate stop condition across both states and all AD
// levels of dxdt.
template<class State>
inline double max_abs_all_levels_vec(const State& v) {
  double m = 0.0;
  for (std::size_t i = 0; i < v.size(); ++i)
    m = std::max(m, max_abs_all_levels(v[i]));
  return m;
}

} // namespace cppde

#endif // CPPDE_DUAL_MATH_HPP
