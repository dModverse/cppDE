/*
 Scalar building blocks that every generated model needs, with or without AD.

 value_of(x) is the uniform value accessor. The AD types and the expression
 templates add their own overloads in their own headers, so a single spelling
 reaches the innermost scalar of anything a code generator can emit.

 select(cond, a, b) replaces the conditional operator in generated code. A
 ternary cannot type check when one branch is an expression-template node and
 the other a literal, because neither converts to the other; a function call
 converts both operands to a common node instead. Both arguments are evaluated,
 so a guarded branch has to be safe to evaluate.

 Copyright (C) 2026 Simon Beyer
 */

#ifndef CPPDE_SCALAR_OPS_HPP
#define CPPDE_SCALAR_OPS_HPP

#include <type_traits>

namespace cppde {

template<class S, std::enable_if_t<std::is_arithmetic<S>::value, int> = 0>
constexpr S value_of(const S& s) { return s; }

template<class A, class B,
         std::enable_if_t<std::is_arithmetic<A>::value
                          && std::is_arithmetic<B>::value, int> = 0>
constexpr std::common_type_t<A, B> select(bool cond, A a, B b) {
  using R = std::common_type_t<A, B>;
  return cond ? static_cast<R>(a) : static_cast<R>(b);
}

} // namespace cppde

#endif // CPPDE_SCALAR_OPS_HPP
