/*
 cppde::dual2nd<T, N>: second-order forward-mode AD type.

 dual2nd publicly inherits from dual<dual<T, N>, N>. Inheritance preserves the
 storage layout and the full accessor surface (.x(), .d(j), .size(), .depend(),
 .diff(), operator[]) so that the generic AD machinery (cppde_ad_lu.hpp,
 cppde_newton.hpp, cppde_integrate_times.hpp, the LU IFT recursion in
 cppde_ad_traits.hpp) continues to work unchanged: a dual2nd is-a nested
 dual.

 The reason for a distinct type is dispatch. Math primitives in
 cppde_dual2nd_math.hpp are templated on dual2nd<T, N> specifically and
 exploit Hessian symmetry by computing only the lower triangle (j <= i) and
 mirroring to the upper. Function-template argument deduction is strict: the
 dual2nd-specific operators match dual2nd exactly and do NOT match the base
 dual<dual<T, N>, N> via subclass slicing. Conversely, the eager nested-dual
 operators in cppde_dual_math.hpp deduce on dual<dual<T, N>, N> and do NOT
 match dual2nd. The two operator sets are unambiguous.

 Convenience accessors d1_at(i) / dd_at(i, j) translate to the underlying
 nested-dual storage:
   - d1_at(i) -> outer.tan_[i].x()  (gradient slot, also redundantly stored
                                     in outer.val_.tan_[i] for LU-IFT
                                     correctness; primitives mirror writes)
   - dd_at(i, j) -> outer.tan_[max(i,j)].tan_[min(i,j)]  (canonical lower
                                     triangle, fed by symmetric computation)

 Storage layout for now matches dual<dual<T, N>, N> (full N x N inner-tangent
 block, gradient redundantly stored across both outer-val and outer-tan).
 The compute saving comes from math primitives only filling the lower
 triangle and mirroring at the end of each operation. Storage compaction to
 a packed N(N+1)/2 Hessian is a follow-up that requires reworking the LU
 IFT extraction pipeline (the recursive layer-by-layer model assumes .x()
 peels exactly one AD layer; a packed dual2nd peels two).

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

// ----------------------------------------------------------------------------
// Triangular indexing helpers. tri_idx canonicalises (i, j) to the lower
// triangle and returns the packed offset; tri_size(N) is the total entry
// count in the lower triangle. Used by math primitives that iterate the
// symmetric Hessian in canonical order (j <= i).
// ----------------------------------------------------------------------------
namespace detail {
constexpr inline unsigned tri_idx(unsigned i, unsigned j) noexcept {
  return (j <= i) ? (i * (i + 1u) / 2u + j) : (j * (j + 1u) / 2u + i);
}
constexpr inline unsigned tri_size(unsigned n) noexcept {
  return n * (n + 1u) / 2u;
}
} // namespace detail

// =============================================================================
// dual2nd<T, N>: distinct type, public-inheritance refinement of nested dual.
// All storage and standard accessors are inherited; the new type adds
// triangular convenience accessors used by math primitives.
// =============================================================================
template<class T = double, unsigned N = 0>
class dual2nd : public dual<dual<T, N>, N> {
public:
  using base         = dual<dual<T, N>, N>;
  using inner_dual_t = dual<T, N>;
  using value_type   = inner_dual_t;
  static constexpr unsigned static_size = N;

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

  // Implicit copy / move are generated and copy/move the base subobject.
  // No new data members, so default copy/move are correct.

    // ---------------------------------------------------------------------------
    // Triangular accessors. The gradient lives in outer.tan_[i].x(), which is
    // what d1_at(i) reads. dd_at(i, j) is the canonical lower-triangle Hessian
    // slot and canonicalises an unordered index pair.
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
    // outer.tan_[max(i,j)].tan_[min(i,j)]. Both overloads canonicalise the pair,
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

    // Raw (i, j) inner-tangent slot, without canonicalisation. mirror_upper_dd
    // copies the lower triangle into the upper one with it, so the LU's IFT
    // extraction finds a populated [j][i] cell.
  T& dd_raw(unsigned i, unsigned j) {
    return base::operator[](i)[j];
  }
  const T& dd_raw(unsigned i, unsigned j) const {
    return base::operator[](i)[j];
  }

    // No-op. The LU reads the gradient from the inline outer.tan_[k].x() slot
    // through first_order_view, so there is nothing to mirror. Kept for the call
    // sites that still invoke it.
  void sync_d1_redundant() noexcept {}

    // Arm the outer tangent slots so d1_at and dd_at write into allocated
    // storage: afterwards base::depend() and every base::operator[](i).depend()
    // hold. The outer.val_ layer stays unarmed, the LU does not read it.
  void arm_full(unsigned m) {
    if constexpr (N > 0) {
      (void)m;
      if (!this->depend()) base::set_depend();
      for (unsigned i = 0; i < N; ++i) {
        auto& ti = base::operator[](i);
        if (!ti.depend()) ti.set_depend();
      }
    } else {
      if (!this->depend()) base::set_depend_size(m);
      for (unsigned i = 0; i < m; ++i) {
        auto& ti = base::operator[](i);
        if (!ti.depend()) ti.set_depend_size(m);
      }
    }
  }

  // Width helper. For static-N this is the compile-time N; for dynamic-N
  // it is the runtime size. Math primitives need this as the loop bound.
  unsigned width() const {
    if constexpr (N > 0) return N;
    else                  return this->size();
  }

  // Innermost scalar f(theta). Equivalent to base::x().x(), peeling both AD
  // layers. Math primitives read/write this directly to avoid going through
  // the nested-dual eager operators (which would re-evaluate the value layer
  // as a dual<T,N> expression instead of a single T).
  T& scalar()             { return base::x().x(); }
  const T& scalar() const { return base::x().x(); }

    // -- Expression-template assignment / construction -------------------------
    // Materialises an Expr2nd<D> tree into this dual2nd's storage, defined in
    // cppde_dual2nd_expr.hpp. Deduction matches CRTP derivations only, so the
    // dual2nd-to-dual2nd and scalar paths fall back to the base operators.
  template<class D>
  inline dual2nd& operator=(const dual2nd_expr::Expr2nd<D>& e);

  template<class D>
  inline dual2nd(const dual2nd_expr::Expr2nd<D>& e);

  // Compound assignment from Expr2nd: route through (*this op other) which
  // builds a BinExpr2 and materialises into *this. The eager nested base
  // operator+=(const dual<dual<T,N>,N>&) is preserved for plain dual2nd
  // operands via inheritance.
  template<class D>
  inline dual2nd& operator+=(const dual2nd_expr::Expr2nd<D>& e);
  template<class D>
  inline dual2nd& operator-=(const dual2nd_expr::Expr2nd<D>& e);
  template<class D>
  inline dual2nd& operator*=(const dual2nd_expr::Expr2nd<D>& e);
  template<class D>
  inline dual2nd& operator/=(const dual2nd_expr::Expr2nd<D>& e);
};

// ----------------------------------------------------------------------------
// Synthesise a first-order dual<S, N> from a dual2nd's scalar + inline gradient
// (outer.tan_[k].x()), bypassing the redundant val_tan_block. Used by the LU
// dual2nd dispatch to extract the value layer without requiring val_tan to be
// kept in sync via sync_d1_redundant.
// ----------------------------------------------------------------------------
template<class S, unsigned N>
inline cppde::dual<S, N> first_order_view(const cppde::dual2nd<S, N>& v) {
  cppde::dual<S, N> r;
  r.x() = v.scalar();
  if (v.depend()) {
    if constexpr (N > 0) {
      r.set_depend_size();
      for (unsigned k = 0; k < N; ++k) r[k] = v.d1_at(k);
    } else {
      const unsigned m = v.size();
      if (m > 0) {
        r.set_depend_size(m);
        for (unsigned k = 0; k < m; ++k) r[k] = v.d1_at(k);
      }
    }
  }
  return r;
}

// value_of on a dual2nd: the innermost scalar, matching the first-order
// overload so an emitter has one spelling for both orders.
template<class T, unsigned N>
inline auto value_of(const dual2nd<T, N>& d) { return value_of(d.scalar()); }

} // namespace cppde

#endif // CPPDE_DUAL2ND_HPP
