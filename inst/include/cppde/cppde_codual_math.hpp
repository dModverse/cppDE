/*
 Arithmetic operators, math functions and comparisons for cppde::codual<T>.

 Mirrors the surface of cppde_dual_math.hpp, so a model body compiles unchanged
 under either type. Where dual multiplies the local partial into every tangent,
 codual records it and the sweep does the multiplication.

 Convention:
 - Unary  y = f(x):    record(x.slot, f'(x.val))
 - Binary y = f(a, b): record(a.slot, f_a, b.slot, f_b)
 - An operand without dependence contributes no operand slot, which keeps
   constants off the tape rather than recording a zero partial.
 - A single operand whose partial is one needs no node either: the result names
   the operand's own slot, and the sweep is spared an edge it would multiply by
   one. That covers every shift by a constant, and a stepper does many.
 - Comparisons fall back to .x(), as they do for dual.

 No eager gate: these are the only definitions of these operators for codual.

 Copyright (C) 2026 Simon Beyer
 */

#ifndef CPPDE_CODUAL_MATH_HPP
#define CPPDE_CODUAL_MATH_HPP

#include <cppde/cppde_codual.hpp>

#include <cmath>
#include <type_traits>

namespace cppde {

namespace detail {

// The two recording shapes every operator below reduces to, with the depend()
// case analysis in one place.
template<class T>
inline codual<T> codual_unary(const T& value, const codual<T>& a, const T& fa) {
  codual_tape<T>& tp = codual_tape_for<T>();
  if (!tp.live(a.slot())) return codual<T>(value, codual<T>::none);
  return codual<T>(value, tp.record(a.slot(), fa));
}

// y = x + c: the derivative is the identity, so the result is the same tape
// node under another value. A node with one operand and partial one is exactly
// that identity, so recording it would add work and no information.
template<class T>
inline codual<T> codual_shift(const T& value, const codual<T>& a) {
  return codual<T>(value, a.slot());
}

template<class T>
inline codual<T> codual_binary(const T& value,
                               const codual<T>& a, const T& fa,
                               const codual<T>& b, const T& fb) {
  codual_tape<T>& tp = codual_tape_for<T>();
  // Resolved once and handed to the tape, rather than tested here and resolved
  // again inside record(). One node covers all three live combinations: a dead
  // operand is `nolocal`, and the sweep skips it.
  const unsigned la = tp.local(a.slot());
  const unsigned lb = tp.local(b.slot());
  if (la == codual_tape<T>::nolocal && lb == codual_tape<T>::nolocal)
    return codual<T>(value, codual<T>::none);
  return codual<T>(value, tp.record_local(la, fa, lb, fb));
}

}  // namespace detail

// =============================================================================
// Arithmetic operators (codual op codual)
// =============================================================================

template<class T>
inline codual<T> operator+(const codual<T>& a, const codual<T>& b) {
  return detail::codual_binary(a.x() + b.x(), a, T(1), b, T(1));
}

template<class T>
inline codual<T> operator-(const codual<T>& a, const codual<T>& b) {
  return detail::codual_binary(a.x() - b.x(), a, T(1), b, T(-1));
}

template<class T>
inline codual<T> operator*(const codual<T>& a, const codual<T>& b) {
  return detail::codual_binary(a.x() * b.x(), a, b.x(), b, a.x());
}

template<class T>
inline codual<T> operator/(const codual<T>& a, const codual<T>& b) {
  const T inv = T(1) / b.x();
  const T y   = a.x() * inv;
  return detail::codual_binary(y, a, inv, b, -y * inv);
}

// Unary
template<class T>
inline codual<T> operator+(const codual<T>& a) { return a; }

template<class T>
inline codual<T> operator-(const codual<T>& a) {
  return detail::codual_unary(-a.x(), a, T(-1));
}

// =============================================================================
// Mixed codual/scalar operators
// =============================================================================

#define CPPDE_CODUAL_MIXED(OP, VAL_AS, FA_AS, VAL_SA, FA_SA)                  \
  template<class T, class U,                                                  \
           std::enable_if_t<std::is_arithmetic_v<U>, int> = 0>                \
  inline codual<T> operator OP (const codual<T>& a, const U& s) {             \
    const T sv = static_cast<T>(s); (void) sv;                                \
    return detail::codual_unary((VAL_AS), a, (FA_AS));                        \
  }                                                                           \
  template<class T, class U,                                                  \
           std::enable_if_t<std::is_arithmetic_v<U>, int> = 0>                \
  inline codual<T> operator OP (const U& s, const codual<T>& a) {             \
    const T sv = static_cast<T>(s); (void) sv;                                \
    return detail::codual_unary((VAL_SA), a, (FA_SA));                        \
  }

CPPDE_CODUAL_MIXED(*, a.x() * sv, sv,    sv * a.x(), sv)
CPPDE_CODUAL_MIXED(/, a.x() / sv, T(1) / sv,
                      sv / a.x(), -sv / (a.x() * a.x()))

#undef CPPDE_CODUAL_MIXED

// + and - against a scalar are shifts and carry the operand's slot instead of a
// node. Only s - a turns the derivative over, and that one records.
template<class T, class U, std::enable_if_t<std::is_arithmetic_v<U>, int> = 0>
inline codual<T> operator+(const codual<T>& a, const U& s) {
  return detail::codual_shift(a.x() + static_cast<T>(s), a);
}
template<class T, class U, std::enable_if_t<std::is_arithmetic_v<U>, int> = 0>
inline codual<T> operator+(const U& s, const codual<T>& a) {
  return detail::codual_shift(static_cast<T>(s) + a.x(), a);
}
template<class T, class U, std::enable_if_t<std::is_arithmetic_v<U>, int> = 0>
inline codual<T> operator-(const codual<T>& a, const U& s) {
  return detail::codual_shift(a.x() - static_cast<T>(s), a);
}
template<class T, class U, std::enable_if_t<std::is_arithmetic_v<U>, int> = 0>
inline codual<T> operator-(const U& s, const codual<T>& a) {
  return detail::codual_unary(static_cast<T>(s) - a.x(), a, T(-1));
}

// =============================================================================
// Compound assignment
//
// Routed through the free operators, so the recording rules are stated once.
// The scalar += and -= in the class body keep the slot: d(x + c)/dx is one.
// =============================================================================

template<class T>
inline codual<T>& codual<T>::operator+=(const codual<T>& o) { return *this = *this + o; }
template<class T>
inline codual<T>& codual<T>::operator-=(const codual<T>& o) { return *this = *this - o; }
template<class T>
inline codual<T>& codual<T>::operator*=(const codual<T>& o) { return *this = *this * o; }
template<class T>
inline codual<T>& codual<T>::operator/=(const codual<T>& o) { return *this = *this / o; }

template<class T>
template<class U, std::enable_if_t<std::is_arithmetic_v<U>, int>>
inline codual<T>& codual<T>::operator*=(const U& v) { return *this = *this * v; }
template<class T>
template<class U, std::enable_if_t<std::is_arithmetic_v<U>, int>>
inline codual<T>& codual<T>::operator/=(const U& v) { return *this = *this / v; }

// =============================================================================
// Math: unary functions
//
// VAL_EXPR and DERIV_EXPR may invoke any std math function, so the whole suite
// is in scope, keeping ADL working for T = double and for a user namespace.
// =============================================================================

#define CPPDE_CODUAL_UNARY(NAME, VAL_EXPR, DERIV_EXPR)                        \
  template<class T>                                                           \
  inline codual<T> NAME(const codual<T>& a) {                                 \
    using std::exp;   using std::log;   using std::sqrt;                      \
    using std::sin;   using std::cos;   using std::tan;                       \
    using std::asin;  using std::acos;  using std::atan;                      \
    using std::sinh;  using std::cosh;  using std::tanh;                      \
    using std::asinh; using std::acosh; using std::atanh;                     \
    const T xv = a.x();                                                       \
    const T yv = (VAL_EXPR);                                                  \
    if (!a.depend()) return codual<T>(yv, codual<T>::none);                   \
    return detail::codual_unary(yv, a, (DERIV_EXPR));                         \
  }

CPPDE_CODUAL_UNARY(exp,   exp(xv),   yv)
CPPDE_CODUAL_UNARY(log,   log(xv),   T(1) / xv)
CPPDE_CODUAL_UNARY(sqrt,  sqrt(xv),  T(1) / (T(2) * yv))
CPPDE_CODUAL_UNARY(sin,   sin(xv),   cos(xv))
CPPDE_CODUAL_UNARY(cos,   cos(xv),  -sin(xv))
CPPDE_CODUAL_UNARY(tan,   tan(xv),   T(1) + yv * yv)
CPPDE_CODUAL_UNARY(asin,  asin(xv),  T(1) / sqrt(T(1) - xv * xv))
CPPDE_CODUAL_UNARY(acos,  acos(xv), -T(1) / sqrt(T(1) - xv * xv))
CPPDE_CODUAL_UNARY(atan,  atan(xv),  T(1) / (T(1) + xv * xv))
CPPDE_CODUAL_UNARY(sinh,  sinh(xv),  cosh(xv))
CPPDE_CODUAL_UNARY(cosh,  cosh(xv),  sinh(xv))
CPPDE_CODUAL_UNARY(tanh,  tanh(xv),  T(1) - yv * yv)
CPPDE_CODUAL_UNARY(asinh, asinh(xv), T(1) / sqrt(xv * xv + T(1)))
CPPDE_CODUAL_UNARY(acosh, acosh(xv), T(1) / sqrt(xv * xv - T(1)))
CPPDE_CODUAL_UNARY(atanh, atanh(xv), T(1) / (T(1) - xv * xv))

#undef CPPDE_CODUAL_UNARY

// Non-differentiable at 0, where the subgradient 0 is taken, as dual does.
template<class T>
inline codual<T> abs(const codual<T>& a) {
  using std::abs;
  const T xv = a.x();
  const T fp = (xv > T(0)) ? T(1) : ((xv < T(0)) ? T(-1) : T(0));
  return detail::codual_unary(abs(xv), a, fp);
}

// =============================================================================
// Math: pow
//
// y = a^b  =>  dy = b * a^(b-1) * da + log(a) * a^b * db
// The b partial is formed only when b depends on something: log(a) is undefined
// for a <= 0.
// =============================================================================

template<class T>
inline codual<T> pow(const codual<T>& a, const codual<T>& b) {
  using std::pow; using std::log;
  const T av = a.x(), bv = b.x();
  const T y  = pow(av, bv);
  const T fa = bv * pow(av, bv - T(1));
  const T fb = b.depend() ? log(av) * y : T(0);
  return detail::codual_binary(y, a, fa, b, fb);
}

template<class T, class U,
         std::enable_if_t<std::is_arithmetic_v<U>, int> = 0>
inline codual<T> pow(const codual<T>& a, const U& b) {
  using std::pow;
  const T bv = static_cast<T>(b);
  return detail::codual_unary(pow(a.x(), bv), a, bv * pow(a.x(), bv - T(1)));
}

template<class T, class U,
         std::enable_if_t<std::is_arithmetic_v<U>, int> = 0>
inline codual<T> pow(const U& a, const codual<T>& b) {
  using std::pow; using std::log;
  const T av = static_cast<T>(a);
  const T y  = pow(av, b.x());
  return detail::codual_unary(y, b, log(av) * y);
}

// =============================================================================
// Math: selection
//
// A branch selects one operand whole, so the result already names that
// operand's tape slot and nothing is recorded.
// =============================================================================

template<class T>
inline codual<T> min(const codual<T>& a, const codual<T>& b) {
  return (a.x() < b.x()) ? a : b;
}
template<class T>
inline codual<T> max(const codual<T>& a, const codual<T>& b) {
  return (a.x() < b.x()) ? b : a;
}
template<class T, class U, std::enable_if_t<std::is_arithmetic_v<U>, int> = 0>
inline codual<T> min(const codual<T>& a, const U& b) {
  return (a.x() < static_cast<T>(b)) ? a : codual<T>(b);
}
template<class T, class U, std::enable_if_t<std::is_arithmetic_v<U>, int> = 0>
inline codual<T> min(const U& a, const codual<T>& b) { return min(b, a); }
template<class T, class U, std::enable_if_t<std::is_arithmetic_v<U>, int> = 0>
inline codual<T> max(const codual<T>& a, const U& b) {
  return (a.x() < static_cast<T>(b)) ? codual<T>(b) : a;
}
template<class T, class U, std::enable_if_t<std::is_arithmetic_v<U>, int> = 0>
inline codual<T> max(const U& a, const codual<T>& b) { return max(b, a); }

template<class T, class L, class H>
inline codual<T> clamp(const codual<T>& a, const L& lo, const H& hi) {
  return min(max(a, lo), hi);
}

// =============================================================================
// select
//
// The conditional operator as a call, the spelling the piecewise codegen emits.
// The branch is a control decision, so the chosen operand is returned whole and
// nothing is recorded. Both branches are evaluated, as cppde_scalar_ops.hpp
// requires.
// =============================================================================

template<class T>
inline codual<T> select(bool cond, const codual<T>& a, const codual<T>& b) {
  return cond ? a : b;
}
template<class T, class U, std::enable_if_t<std::is_arithmetic_v<U>, int> = 0>
inline codual<T> select(bool cond, const codual<T>& a, const U& b) {
  return cond ? a : codual<T>(b);
}
template<class T, class U, std::enable_if_t<std::is_arithmetic_v<U>, int> = 0>
inline codual<T> select(bool cond, const U& a, const codual<T>& b) {
  return cond ? codual<T>(a) : b;
}

// =============================================================================
// Comparisons
// =============================================================================

#define CPPDE_CODUAL_CMP(OP)                                                  \
  template<class T>                                                           \
  inline bool operator OP (const codual<T>& a, const codual<T>& b) {          \
    return a.x() OP b.x();                                                    \
  }                                                                           \
  template<class T, class U,                                                  \
           std::enable_if_t<std::is_arithmetic_v<U>, int> = 0>                \
  inline bool operator OP (const codual<T>& a, const U& b) {                  \
    return a.x() OP static_cast<T>(b);                                        \
  }                                                                           \
  template<class T, class U,                                                  \
           std::enable_if_t<std::is_arithmetic_v<U>, int> = 0>                \
  inline bool operator OP (const U& a, const codual<T>& b) {                  \
    return static_cast<T>(a) OP b.x();                                        \
  }

CPPDE_CODUAL_CMP(==)
CPPDE_CODUAL_CMP(!=)
CPPDE_CODUAL_CMP(<)
CPPDE_CODUAL_CMP(<=)
CPPDE_CODUAL_CMP(>)
CPPDE_CODUAL_CMP(>=)

#undef CPPDE_CODUAL_CMP

}  // namespace cppde

#endif  // CPPDE_CODUAL_MATH_HPP
