/*
 Centralized AD type traits for cppDE.

 is_ad / is_dual2nd / inner_type / scalar_type / scalar_value plus the bulk
 helpers max_deriv_size and bulk_inject_results. Specialized for
 cppde::dual<T> and cppde::dual2nd<T>.

 The bulk helpers are written generically over any AD type that exposes the
 accessor surface: `.x()`, `.d(j)`, `.size()`, `.depend()`,
 `.diff(idx[, n])`, `operator[](j)`.

 Copyright (C) 2026 Simon Beyer
 */

#ifndef CPPDE_AD_TRAITS_HPP
#define CPPDE_AD_TRAITS_HPP

#include <type_traits>
#include <vector>

#include <cppde/cppde_types.hpp>

// Forward declaration for the dual backend so traits can specialize without
// pulling in the full headers from every call site.
namespace cppde {
  template<class T> class dual;
  template<class T> class dual2nd;
}

namespace cppde {
namespace ad_traits {

// ============================================================================
//  is_ad<T>
// ============================================================================

template<class T> struct is_ad : std::false_type {};
template<class T> struct is_ad<cppde::dual<T>>    : std::true_type {};
template<class T> struct is_ad<cppde::dual2nd<T>> : std::true_type {};

// is_dual2nd<T>: matches only cppde::dual2nd<S>, not the base class dual<dual<S>>. The
// LU, slab and multistepper paths dispatch on it to extraction routines that
// read the gradient from outer.tan_[k].x().
template<class T> struct is_dual2nd : std::false_type {};
template<class T> struct is_dual2nd<cppde::dual2nd<T>> : std::true_type {};

// ============================================================================
//  inner_type<T>: strip one layer of AD wrapping.
//
//  dual2nd peels to dual<T>, which holds value and gradient; the LU IFT
//  recursion runs dual2nd<T> -> dual<T> -> T.
// ============================================================================

template<class T> struct inner_type                    { using type = T; };
template<class T> struct inner_type<cppde::dual<T>>    { using type = T; };
template<class T> struct inner_type<cppde::dual2nd<T>> { using type = cppde::dual<T>; };
template<class T> using inner_type_t = typename inner_type<T>::type;

// ============================================================================
//  scalar_type<T>: recursively unwrap to innermost non-AD scalar
// ============================================================================

template<class T> struct scalar_type                    { using type = T; };
template<class T> struct scalar_type<cppde::dual<T>>     : scalar_type<T> {};
template<class T> struct scalar_type<cppde::dual2nd<T>>  : scalar_type<T> {};
template<class T> using scalar_type_t = typename scalar_type<T>::type;

// ============================================================================
//  scalar_value(v): extract innermost double from any (nested) AD type
// ============================================================================

template<class T>
inline typename std::enable_if<std::is_arithmetic<T>::value, double>::type
scalar_value(const T& v) { return static_cast<double>(v); }

template<class T>
inline double scalar_value(const cppde::dual<T>& v) {
  return scalar_value(v.x());
}

template<class T>
inline double scalar_value(const cppde::dual2nd<T>& v) {
  return scalar_value(v.x());
}

// ============================================================================
//  store_as<T>(v): the value in the type a store keeps it in
//
//  A reverse trajectory checkpoints states in the scalar type of the run. In plain
//  double that is the value alone; over the same type the run integrates in it
//  is the whole number, tangents included, which is what forward over reverse
//  needs.
// ============================================================================

// ============================================================================
//  arm_tangents(v): bind v's tangent storage before a callee's arena scope
//
//  A generated model body opens a private dual_arena::scope, so a tangent
//  it allocates for one of the caller's numbers dies when that scope pops. The
//  stepper avoids this by holding the stepper buffers in a tangent slab; a caller that
//  hands the model a caller-owned buffer arms it instead. Nothing to do for a
//  plain scalar.
// ============================================================================

template<class T>
inline void arm_tangents(T&) {}

template<class S>
inline void arm_tangents(cppde::dual<S>& v) { v.arm(); }

// ============================================================================
//  arm_outputs(a, b, x, p): binds the tangent storage of a and b, as wide as x
//  and p are. An entry that already has storage is left unchanged.
// ============================================================================

template<class T>
inline unsigned tangent_count(const T&) { return 0; }

template<class S>
inline unsigned tangent_count(const cppde::dual<S>& v) { return v.size(); }

template<class T>
inline void arm_width(T&, unsigned) {}

template<class T> struct is_heap_dual : std::false_type {};
template<class S> struct is_heap_dual<cppde::dual<S>> : std::true_type {};

// A second-order dual also binds the rows of the second derivatives, as wide as
// the first.
template<class S>
inline void arm_width(cppde::dual<S>& v, unsigned w) {
  if (w == 0) return;
  if (v.size() == 0) {
    v.set_depend_size(w);
    for (unsigned i = 0; i < w; ++i) v[i] = S();
  }
  if constexpr (is_heap_dual<S>::value)
    for (unsigned i = 0; i < v.size(); ++i) arm_width(v[i], w);
}

// dual2nd derives from the nested dual, so the templates for that base lose to
// the generic ones by their conversion; these take it exactly and bind both
// layers zeroed, so that an entry the body does not write reads as zero.
template<class S>
inline unsigned tangent_count(const cppde::dual2nd<S>& v) { return v.width(); }

template<class S>
inline void arm_width(cppde::dual2nd<S>& v, unsigned w) {
  auto& b = static_cast<typename cppde::dual2nd<S>::base&>(v);
  arm_width(b, w);
}

template<class T>
inline void arm_outputs(std::vector<T>& a, std::vector<T>& b,
                        const std::vector<T>& x, const std::vector<T>& p) {
  unsigned w = 0;
  for (const T& v : x) { const unsigned c = tangent_count(v); if (c > w) w = c; }
  for (const T& v : p) { const unsigned c = tangent_count(v); if (c > w) w = c; }
  for (T& v : a) arm_width(v, w);
  for (T& v : b) arm_width(v, w);
}

template<class T, class V>
inline T store_as(const V& v) {
  if constexpr (std::is_same<T, V>::value) return v;
  else return static_cast<T>(scalar_value(v));
}

// ============================================================================
//  step_coef<TimeArg>: the type the steppers combine their stages in.
//
//  The step size is a constant of the map being differentiated, the grid being
//  frozen, so it has no tangent and the stage weights are plain doubles
//  whatever the state is integrated in.
// ============================================================================

template<class TimeArg> using step_coef_t = double;

template<class TimeArg>
inline double step_coef_of(const TimeArg& dt) { return scalar_value(dt); }

// ============================================================================
//  Bulk extraction / injection helpers (generic over any AD type with the
//  standard accessor surface).
// ============================================================================

// Max number of active derivative directions in a vector
template<class AD,
         std::enable_if_t<is_ad<AD>::value, int> = 0>
inline unsigned max_deriv_size(const std::vector<AD>& v)
{
  unsigned mx = 0;
  for (std::size_t i = 0; i < v.size(); ++i) {
    unsigned sz = const_cast<AD&>(v[i]).size();
    if (sz > mx) mx = sz;
  }
  return mx;
}

// Max number of active derivative directions in a dense_matrix
template<class AD,
         std::enable_if_t<is_ad<AD>::value, int> = 0>
inline unsigned max_deriv_size(const dense_matrix<AD>& M)
{
  unsigned mx = 0;
  for (std::size_t k = 0; k < M.data.size(); ++k) {
    unsigned sz = const_cast<AD&>(M.data[k]).size();
    if (sz > mx) mx = sz;
  }
  return mx;
}

// Max number of active derivative directions in a csc_matrix. The sparse
// solver's dual specialization reads the entries' derivatives contiguously and
// never asks; the nested-dual one keeps the matrix and does.
template<class AD,
         std::enable_if_t<is_ad<AD>::value, int> = 0>
inline unsigned max_deriv_size(const csc_matrix<AD>& M)
{
  unsigned mx = 0;
  for (std::size_t k = 0; k < M.Ax.size(); ++k) {
    unsigned sz = const_cast<AD&>(M.Ax[k]).size();
    if (sz > mx) mx = sz;
  }
  return mx;
}

// The value layer of b from v; for a dual2nd through set_value_layer, which
// binds no tangents there.
template<class AD, class V>
inline void set_value(AD& b, const V& v) { b.x() = v; }
template<class S, class V>
inline void set_value(cppde::dual2nd<S>& b, const V& v) { b.set_value_layer(v); }

// One pass from a column-major n by n_derivs scratch into b, activating
// dependence where it is not yet set and overwriting every tangent slot.
template<class AD,
         std::enable_if_t<is_ad<AD>::value, int> = 0>
inline void bulk_inject_results(
    std::vector<AD>& b,
    const std::vector<inner_type_t<AD>>& b_val,
    const std::vector<inner_type_t<AD>>& dx_all,
    int n, unsigned n_derivs)
{
  for (int i = 0; i < n; ++i) {
    set_value(b[i], b_val[i]);
    if (n_derivs > 0) {
      // diff(0, n_derivs) sizes the tangent buffer; every slot is overwritten below.
      if (!b[i].depend()) b[i].diff(0, n_derivs);
      for (unsigned j = 0; j < n_derivs; ++j)
        b[i][j] = dx_all[j * n + i];
    }
  }
}

} // namespace ad_traits

namespace detail {

// The names the rest of the library reaches for.
using ad_traits::scalar_value;

} // namespace detail
} // namespace cppde

#endif // CPPDE_AD_TRAITS_HPP
