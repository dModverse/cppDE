/*
 cppde::dual2nd<T>: second-order forward-mode AD type.

 dual2nd publicly inherits from dual<dual<T>>, so the storage layout and
 accessor surface (.x(), .d(j), .size(), .depend(), .diff(), operator[]) serve
 the generic AD machinery and the LU's IFT recursion.

 The distinct type is for dispatch. The dual2nd operators (cppde_dual2nd_expr.hpp,
 cppde_dual2nd_math.hpp) match dual2nd exactly and fill only the lower Hessian
 triangle (j <= i), mirrored to the upper; the eager nested-dual operators in
 cppde_dual_math.hpp do not match dual2nd.

   d1_at(i)    -> outer.tan_[i].x()                      gradient
   dd_at(i, j) -> outer.tan_[max(i,j)].tan_[min(i,j)]    canonical lower triangle

 The Hessian is stored as the full m x m inner block, m the tangent width. The LU does not read
 outer.val_.tan_.

 Copyright (C) 2026 Simon Beyer
 */

#ifndef CPPDE_DUAL2ND_HPP
#define CPPDE_DUAL2ND_HPP

#include <cppde/cppde_dual.hpp>

#include <cassert>
#include <cstddef>
#include <type_traits>

namespace cppde {

// Forward-declare the ET CRTP base so dual2nd can declare templated
// assignment / construction from any Expr2nd<D> tree. Definitions live in
// cppde_dual2nd_expr.hpp.
namespace dual2nd_expr {
  template<class Derived> struct Expr2nd;
}

// =============================================================================
// dual2nd<T>: distinct type, public-inheritance refinement of nested dual.
// All storage and standard accessors are inherited; the new type adds
// triangular convenience accessors used by math primitives.
// =============================================================================
template<class T = double>
class dual2nd : public dual<dual<T>> {
public:
  using base         = dual<dual<T>>;
  using inner_dual_t = dual<T>;
  using value_type   = inner_dual_t;

    // Inherit the base ctors and assignments: default, converting from an
    // arithmetic type, copy and move, and the compound operators for plain dual
    // and arithmetic operands.
  using base::base;
  using base::operator=;
  using base::operator+=;
  using base::operator-=;
  using base::operator*=;
  using base::operator/=;

  // Default ctor must be redeclared explicitly because using-base-ctors
  // does not import the implicit default ctor in all compiler versions.
  dual2nd() : base() {}

    // Copies and moves never bind the tangents of the value layer, which repeat
    // the gradient of the outer layer: an entry bound before an arena scope and
    // assigned inside it stays valid after the scope.
  dual2nd(const dual2nd& o) : base() { copy_from(o); }
  dual2nd(dual2nd&& o) noexcept : base(std::move(static_cast<base&>(o))) {}
  dual2nd& operator=(const dual2nd& o) {
    if (this != &o) copy_from(o);
    return *this;
  }
  dual2nd& operator=(dual2nd&& o) noexcept {
    if (this != &o) move_from(o);
    return *this;
  }

    // The scalar of v into the value layer, and the tangents of v where that layer
    // is bound; binds nothing.
  void set_value_layer(const inner_dual_t& v) {
    inner_dual_t& own = base::x();
    if (own.depend()) own = v;
    else own.x() = v.x();
  }

    // ---------------------------------------------------------------------------
    // Triangular accessors. The gradient lives in outer.tan_[i].x(), which is
    // what d1_at(i) reads. dd_at(i, j) is the canonical lower-triangle Hessian
    // slot and canonicalizes an unordered index pair.
    // ---------------------------------------------------------------------------

    // First-order tangent df/dtheta_i, the tan_-side copy. The mutable overload
    // writes into outer.tan_[i] and needs armed storage; the const one goes
    // through base::d(i), which is bounds-safe on an unarmed operand.
  T& d1_at(unsigned i) {
    return base::operator[](i).x();
  }
  const T& d1_at(unsigned i) const {
    return static_cast<const base*>(this)->d(i).x();
  }

    // Second-order tangent d^2f/dtheta_i dtheta_j, stored at
    // outer.tan_[max(i,j)].tan_[min(i,j)]. Both overloads canonicalize the pair,
    // the mutable one needs armed storage.
  T& dd_at(unsigned i, unsigned j) {
    const unsigned r = (j <= i) ? i : j;
    const unsigned c = (j <= i) ? j : i;
    return base::operator[](r)[c];
  }
  const T& dd_at(unsigned i, unsigned j) const {
    const unsigned r = (j <= i) ? i : j;
    const unsigned c = (j <= i) ? j : i;
    return static_cast<const base*>(this)->d(r).d(c);
  }

    // Raw (i, j) inner-tangent slot, without canonicalization. mirror_upper_dd
    // copies the lower triangle into the upper one with it, so the LU's IFT
    // extraction finds a populated [j][i] cell.
  T& dd_raw(unsigned i, unsigned j) {
    return base::operator[](i)[j];
  }
  const T& dd_raw(unsigned i, unsigned j) const {
    return base::operator[](i)[j];
  }

    // Arm the outer tangent slots so d1_at and dd_at write into allocated
    // storage: afterwards base::depend() and every base::operator[](i).depend()
    // hold. The outer.val_ layer stays unarmed, the LU does not read it.
  void arm_full(unsigned m) {
    if (!this->depend()) base::set_depend_size(m);
    for (unsigned i = 0; i < m; ++i) {
      auto& ti = base::operator[](i);
      if (!ti.depend()) ti.set_depend_size(m);
    }
  }

  // The tangent width, the loop bound of the math primitives.
  unsigned width() const { return this->size(); }

  // Innermost scalar f(theta). Equivalent to base::x().x(), peeling both AD
  // layers. Math primitives read/write this directly to avoid going through
  // the nested-dual eager operators (which would re-evaluate the value layer
  // as a dual<T> expression instead of a single T).
  T& scalar()             { return base::x().x(); }
  const T& scalar() const { return base::x().x(); }

    // -- Expression-template assignment / construction -------------------------
    // Materializes an Expr2nd<D> tree into this dual2nd's storage, defined in
    // cppde_dual2nd_expr.hpp. Deduction matches CRTP derivations only, so the
    // dual2nd-to-dual2nd and scalar paths fall back to the base operators.
  template<class D>
  inline dual2nd& operator=(const dual2nd_expr::Expr2nd<D>& e);

  template<class D>
  inline dual2nd(const dual2nd_expr::Expr2nd<D>& e);

  // Compound assignment from Expr2nd: route through (*this op other) which
  // builds a BinExpr2 and materializes into *this. The eager nested base
  // operator+=(const dual<dual<T>>&) is preserved for plain dual2nd
  // operands via inheritance.
  template<class D>
  inline dual2nd& operator+=(const dual2nd_expr::Expr2nd<D>& e);
  template<class D>
  inline dual2nd& operator-=(const dual2nd_expr::Expr2nd<D>& e);
  template<class D>
  inline dual2nd& operator*=(const dual2nd_expr::Expr2nd<D>& e);
  template<class D>
  inline dual2nd& operator/=(const dual2nd_expr::Expr2nd<D>& e);

    // The same with a dual2nd operand, with the Hessian cross terms the inherited
    // nested-dual operators omit. Templates, so that an operand that only
    // converts to dual2nd keeps the inherited overload.
  template<class U, std::enable_if_t<std::is_same_v<U, dual2nd>, int> = 0>
  inline dual2nd& operator+=(const U& o);
  template<class U, std::enable_if_t<std::is_same_v<U, dual2nd>, int> = 0>
  inline dual2nd& operator-=(const U& o);
  template<class U, std::enable_if_t<std::is_same_v<U, dual2nd>, int> = 0>
  inline dual2nd& operator*=(const U& o);
  template<class U, std::enable_if_t<std::is_same_v<U, dual2nd>, int> = 0>
  inline dual2nd& operator/=(const U& o);

private:
    // The base copy with the value layer through set_value_layer.
  void copy_from(const dual2nd& o) {
    set_value_layer(static_cast<const base&>(o).x());
    base& b = *this;
    const base& ob = o;
    const unsigned m = ob.size();
    if (m == 0) {
      for (unsigned i = 0; i < b.size(); ++i) b[i] = inner_dual_t();
      return;
    }
    if (b.size() == 0) b.set_depend_size(m);
    assert(b.size() == m && "dual2nd: tangent-vector size mismatch");
    for (unsigned i = 0; i < m; ++i) b[i] = ob[i];
  }

    // As copy_from, but an unbound outer layer takes the buffer of o and a
    // bound one takes the rows of o.
  void move_from(dual2nd& o) noexcept {
    set_value_layer(static_cast<const base&>(o).x());
    base& b = *this;
    base& ob = o;
    const unsigned m = ob.size();
    if (m == 0) {
      for (unsigned i = 0; i < b.size(); ++i) b[i] = inner_dual_t();
      return;
    }
    if (b.size() > 0) {
      for (unsigned i = 0; i < m; ++i) b[i] = std::move(ob[i]);
      return;
    }
    b.rebind_storage(&ob[0], m);
    ob.rebind_storage(nullptr, 0);
  }
};

// The first-order dual of v, the scalar and inline gradient outer.tan_[k].x(),
// written into the tangent storage dst already has, so that a buffer filled on
// every solve allocates on the first fill only.
template<class S>
inline void assign_first_order(cppde::dual<S>& dst, const cppde::dual2nd<S>& v) {
  const unsigned m = v.size();
  if (!v.depend() || m == 0) { dst = v.scalar(); return; }
  dst.set_depend_size(m);
  dst.x() = v.scalar();
  for (unsigned k = 0; k < m; ++k) dst[k] = v.d1_at(k);
}

// value_of on a dual2nd: the innermost scalar, matching the first-order
// overload so an emitter has one spelling for both orders.
template<class T>
inline auto value_of(const dual2nd<T>& d) { return value_of(d.scalar()); }

} // namespace cppde

#endif // CPPDE_DUAL2ND_HPP
