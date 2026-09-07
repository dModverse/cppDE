/*
 cppde::codual<T>: first-order reverse-mode AD type.

 A codual is a value plus a slot on the thread-local cppde::codual_tape<T>. It
 carries no tangent storage, so its width is independent of the number of
 parameters; the derivatives come out of one backwards sweep over the tape.

 Slot `codual_tape<T>::none` marks a value that depends on no input. Operations
 on such values record nothing, which is what `dual::depend()` does for the
 tangent loop.

 No expression-template overlay: cppde_dual_expr.hpp collapses the temporaries
 that a tape needs.

 This header defines the data class only. Arithmetic, math functions and
 comparisons live in cppde_codual_math.hpp.

 Copyright (C) 2026 Simon Beyer
 */

#ifndef CPPDE_CODUAL_HPP
#define CPPDE_CODUAL_HPP

#include <cppde/cppde_codual_tape.hpp>
#include <cppde/cppde_scalar_ops.hpp>

#include <type_traits>

namespace cppde {

template<class T = double>
class codual {
public:
  using value_type = T;
  using tape_type  = codual_tape<T>;

  static constexpr unsigned none = tape_type::none;

  // -- constructors -----------------------------------------------------------
  codual() : val_(), slot_(none) {}

  template<class U,
           std::enable_if_t<std::is_convertible_v<U, T>, int> = 0>
  codual(const U& v) : val_(static_cast<T>(v)), slot_(none) {}

  codual(const T& v, unsigned slot) : val_(v), slot_(slot) {}

  // A copy names the same tape node: the value was recorded once and copying it
  // creates no new dependence.
  codual(const codual&)            = default;
  codual& operator=(const codual&) = default;

  // Assigning a plain number drops the dependence, as dual::operator=(const U&)
  // clears its tangents.
  template<class U,
           std::enable_if_t<std::is_convertible_v<U, T>, int> = 0>
  codual& operator=(const U& v) {
    val_  = static_cast<T>(v);
    slot_ = none;
    return *this;
  }

  // Compound assignment. Defined out of line in cppde_codual_math.hpp via the
  // free operators, so every recording rule lives in one place.
  codual& operator+=(const codual& o);
  codual& operator-=(const codual& o);
  codual& operator*=(const codual& o);
  codual& operator/=(const codual& o);
  template<class U, std::enable_if_t<std::is_arithmetic_v<U>, int> = 0>
  codual& operator+=(const U& v) { val_ += static_cast<T>(v); return *this; }
  template<class U, std::enable_if_t<std::is_arithmetic_v<U>, int> = 0>
  codual& operator-=(const U& v) { val_ -= static_cast<T>(v); return *this; }
  template<class U, std::enable_if_t<std::is_arithmetic_v<U>, int> = 0>
  codual& operator*=(const U& v);
  template<class U, std::enable_if_t<std::is_arithmetic_v<U>, int> = 0>
  codual& operator/=(const U& v);

  // -- accessors --------------------------------------------------------------
  const T& x()   const { return val_; }
  T&       x()         { return val_; }
  const T& val() const { return val_; }   // .x() alias
  T&       val()       { return val_; }

  unsigned slot()   const { return slot_; }
  void     set_slot(unsigned s) { slot_ = s; }
  bool     depend() const { return slot_ != none; }

  // -- tape interaction -------------------------------------------------------

  // Registers this value as an input. Its adjoint after the sweep is the
  // derivative of the seeded output with respect to it.
  codual& independent() {
    slot_ = codual_tape_for<T>().independent();
    return *this;
  }

  // Seeds this value's adjoint before the sweep. Repeated seeds accumulate.
  void seed(const T& w) const { codual_tape_for<T>().seed(slot_, w); }

  // Valid only after codual_tape<T>::reverse(). A value with no dependence has
  // a zero derivative by construction.
  T adjoint() const {
    if (slot_ == none) return T();
    return codual_tape_for<T>().adjoint(slot_);
  }

private:
  T        val_;
  unsigned slot_;
};

// value_of for the reverse type: peel to the innermost scalar, one spelling
// across the AD types.
template<class T>
inline auto value_of(const codual<T>& d) { return value_of(d.x()); }

}  // namespace cppde

#endif  // CPPDE_CODUAL_HPP
