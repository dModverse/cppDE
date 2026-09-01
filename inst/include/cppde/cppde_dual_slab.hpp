/*
 Contiguous tangent storage for std::vector<dual<T, 0>>.

 Replaces the per-element arena allocations of dual<T, 0>::tan_ with a single
 [n_rows × n_cols] block owned by the multistepper / controller. Each
 dual<T,0>::tan_ then points into one row of the block (set via
 rebind_storage). Subsequent dual = expr materialisations hit the in-place
 reuse branch (size_ matches), so the hot path makes zero arena allocations.

 The slab is sized once per solve via prepare_sensitivities(n_sens) and never
 grown afterward, which keeps the embedded tan_ pointers stable.

 For non-dynamic-dual T (double, static-N dual<T,N!=0>, …) tangent_slab is
 specialised as an empty stub so multistepper<double, …> instances pay no
 size or codegen cost for the slab machinery.

 Copyright (C) 2026 Simon Beyer
 */

#ifndef CPPDE_DUAL_SLAB_HPP
#define CPPDE_DUAL_SLAB_HPP

#include <cassert>
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <type_traits>
#include <vector>

#include <cppde/cppde_dual.hpp>
#include <cppde/cppde_types.hpp>

namespace cppde {
// Forward declaration: dual2nd<S, N> is defined in cppde_dual2nd.hpp;
// the slab specialisation below references it for nested-AD slabbing.
template<class S, unsigned N> class dual2nd;

namespace detail {

// =============================================================================
//  is_dynamic_dual<T>
//
//  True for cppde::dual<S, N> with a non-AD scalar S, that is every dual the
//  slab can back. A nested dual is excluded: its inner tangents would need one
//  buffer per outer tangent, which this layout does not model.
// =============================================================================

template<class T>
struct is_dynamic_dual : std::false_type {};

template<class S, unsigned N>
struct is_dynamic_dual<cppde::dual<S, N>>
  : std::bool_constant<!std::is_class_v<S> || std::is_arithmetic_v<S>>
{};

// dual2nd<S, N> is slab-backed through the two-level specialisation below. In
// external mode nordsieck_block owns the storage and hands each slot its
// pointers; in own-storage mode the slab allocates its blocks itself.
template<class S, unsigned N>
struct is_dynamic_dual<cppde::dual2nd<S, N>>
  : std::bool_constant<!std::is_class_v<S> || std::is_arithmetic_v<S>>
{};

// Trait: is T a dual2nd? Used by tangent_slab to dispatch to the nested
// specialisation.
template<class T>
struct is_dual2nd : std::false_type {};
template<class S, unsigned N>
struct is_dual2nd<cppde::dual2nd<S, N>> : std::true_type {};

// =============================================================================
//  tangent_slab<T>: contiguous [n_rows × n_cols] storage for dual<S,0>.
//
//  Primary template handles dual<S,0>; non-dynamic-dual T uses the empty
//  stub specialisation below.
// =============================================================================

template<class T, bool = is_dynamic_dual<T>::value>
class tangent_slab; // primary template (active path)

template<class T>
class tangent_slab<T, true> {
public:
  using inner_type = typename T::value_type;  // dual<S,0>::value_type == S

  tangent_slab() = default;
  tangent_slab(const tangent_slab&)            = delete;
  tangent_slab& operator=(const tangent_slab&) = delete;
  tangent_slab(tangent_slab&&) noexcept        = default;
  tangent_slab& operator=(tangent_slab&&) noexcept = default;

  unsigned n_rows() const noexcept { return n_rows_; }
  unsigned n_cols() const noexcept { return n_cols_; }
  bool     primed() const noexcept { return n_cols_ != 0; }

  // Raw access to the contiguous tangent block. Layout is row-major:
  // tangents of v[i] live at base + i * n_cols + [0, n_cols). Free helpers
  // (vec_axpy_with_slab / vec_scale_with_slab) call BLAS over this whole
  // block; per-element pointers are still reachable via v[i].
  inner_type*       tangent_data()       noexcept { return base_(); }
  const inner_type* tangent_data() const noexcept { return base_(); }
  std::size_t       tangent_size() const noexcept {
    return static_cast<std::size_t>(n_rows_) * n_cols_;
  }

  // Size the slab and rebind every v[i].tan_ into row i, idempotent when the
  // shape already matches. A static-N dual pins n_cols to N even for a smaller
  // active width: the loops run to N and a shorter row would read the next one.
  void prime(std::vector<T>& v, unsigned n_rows, unsigned n_cols) {
    assert(static_cast<unsigned>(v.size()) == n_rows
           && "tangent_slab::prime: vector size mismatch");
    if constexpr (T::static_size > 0) {
      n_cols = T::static_size;
    }
    if (external_ == nullptr && n_rows == n_rows_ && n_cols == n_cols_) {
      rebind_only(v);
      return;
    }
    n_rows_ = n_rows;
    n_cols_ = n_cols;
    external_ = nullptr;
    storage_.assign(static_cast<std::size_t>(n_rows) * n_cols, inner_type{});
    rebind_only(v);
  }

  // Bind v[i].tan_ into a caller-owned buffer instead of allocating; the owner
  // keeps it alive. nordsieck_block uses this to give all K slots slices of one
  // contiguous block, so BLAS-3 can work across the whole of it.
  void prime_external(std::vector<T>& v, inner_type* base,
                      unsigned n_rows, unsigned n_cols) {
    assert(static_cast<unsigned>(v.size()) == n_rows
           && "tangent_slab::prime_external: vector size mismatch");
    if constexpr (T::static_size > 0) {
      n_cols = T::static_size;
    }
    storage_.clear();
    storage_.shrink_to_fit();
    external_ = base;
    n_rows_   = n_rows;
    n_cols_   = n_cols;
    rebind_only(v);
  }

  // Re-bind v[i].tan_ pointers without resizing storage. Used after a
  // std::vector::resize that may have moved the dual elements (their tan_
  // pointers remain valid as long as storage_ hasn't moved).
  void rebind_only(std::vector<T>& v) {
    if (n_cols_ == 0) return;  // not yet primed; nothing to bind
    assert(static_cast<unsigned>(v.size()) == n_rows_
           && "tangent_slab::rebind_only: vector size changed");
    inner_type* base = base_();
    for (unsigned i = 0; i < n_rows_; ++i) {
      v[i].rebind_storage(base + static_cast<std::size_t>(i) * n_cols_, n_cols_);
    }
  }

private:
  inner_type* base_() noexcept {
    return external_ ? external_ : storage_.data();
  }
  const inner_type* base_() const noexcept {
    return external_ ? external_ : storage_.data();
  }

  std::vector<inner_type> storage_;
  inner_type*             external_ = nullptr;  // non-null = external-storage mode
  unsigned                n_rows_ = 0;
  unsigned                n_cols_ = 0;
};

// Empty stub for non-dynamic-dual T: zero size, all methods no-op.
template<class T>
class tangent_slab<T, false> {
public:
  unsigned n_rows() const noexcept { return 0; }
  unsigned n_cols() const noexcept { return 0; }
  bool     primed() const noexcept { return false; }

  void prime(std::vector<T>&, unsigned, unsigned) noexcept {}
  void rebind_only(std::vector<T>&) noexcept {}
};

// =============================================================================
//  tangent_slab specialisation for cppde::dual2nd<S, N>: two blocks
//
//     outer_storage_:  n_rows * N        dual<S, N>   (outer.tan_)
//     hess_storage_:   n_rows * N * N    S            (Hessian rows)
//
//  Per state i, v[i].base.tan_ points at outer_storage_ + i*N, and
//  outer_storage_[i*N + k].tan_ at hess_storage_ + i*N*N + k*N.
template<class S, unsigned N>
class tangent_slab<cppde::dual2nd<S, N>, true> {
public:
  using value_type   = cppde::dual2nd<S, N>;
  using outer_inner_t = cppde::dual<S, N>;
  using inner_type   = S;

  tangent_slab() = default;
  tangent_slab(const tangent_slab&)            = delete;
  tangent_slab& operator=(const tangent_slab&) = delete;
  tangent_slab(tangent_slab&&) noexcept        = default;
  tangent_slab& operator=(tangent_slab&&) noexcept = default;

  unsigned n_rows() const noexcept { return n_rows_; }
  unsigned n_cols() const noexcept { return n_cols_; }
  bool     primed() const noexcept { return n_cols_ != 0; }

  // hess_data: n_rows * N * N S, the second-order Hessian rows. Flat double
  // array suitable for BLAS daxpy / dscal.
  S*       hess_data()       noexcept { return hess_base_(); }
  const S* hess_data() const noexcept { return hess_base_(); }
  std::size_t hess_size() const noexcept {
    return static_cast<std::size_t>(n_rows_) * n_cols_ * n_cols_;
  }

  // tangent_data: returns the OUTER block (dual<S,N>*), matching the
  // signature T::value_type* expected by the multistepper's nordsieck-
  // evaluation generic loop.
  outer_inner_t*       tangent_data()       noexcept { return outer_base_(); }
  const outer_inner_t* tangent_data() const noexcept { return outer_base_(); }
  std::size_t tangent_size() const noexcept {
    return static_cast<std::size_t>(n_rows_) * n_cols_;
  }

  // 4-arg prime_external (compatibility): owns the inner hess block.
  void prime_external(std::vector<value_type>& v, outer_inner_t* outer_base,
                      unsigned n_rows, unsigned n_cols) {
    if constexpr (N > 0) n_cols = N;
    n_rows_ = n_rows;
    n_cols_ = n_cols;
    external_outer_ = outer_base;
    external_hess_  = nullptr;
    hess_storage_.assign(static_cast<std::size_t>(n_rows) * n_cols * n_cols, S(0));
    outer_storage_.clear();
    outer_storage_.shrink_to_fit();
    rebind_only(v);
  }

  // 6-arg prime_external: nordsieck_block<dual2nd> binds a slot's slab onto
  // two slices (outer, hess) of the unified K-slot blocks. (val_tan_base
  // is accepted for ABI compatibility but ignored.)
  void prime_external(std::vector<value_type>& v,
                      outer_inner_t* outer_base,
                      S* /*val_tan_base*/,
                      S* hess_base,
                      unsigned n_rows, unsigned n_cols) {
    if constexpr (N > 0) n_cols = N;
    n_rows_ = n_rows;
    n_cols_ = n_cols;
    external_outer_ = outer_base;
    external_hess_  = hess_base;
    outer_storage_.clear();
    outer_storage_.shrink_to_fit();
    hess_storage_.clear();
    hess_storage_.shrink_to_fit();
    rebind_only(v);
  }

  void prime(std::vector<value_type>& v, unsigned n_rows, unsigned n_cols) {
    assert(static_cast<unsigned>(v.size()) == n_rows
           && "tangent_slab<dual2nd>::prime: vector size mismatch");
    if constexpr (N > 0) {
      n_cols = N;
    }
    if (external_outer_ == nullptr && external_hess_ == nullptr
        && n_rows == n_rows_ && n_cols == n_cols_) {
      rebind_only(v);
      return;
    }
    n_rows_ = n_rows;
    n_cols_ = n_cols;
    external_outer_ = nullptr;
    external_hess_  = nullptr;
    const std::size_t outer_total = static_cast<std::size_t>(n_rows) * n_cols;
    const std::size_t hess_total  = outer_total * n_cols;
    outer_storage_.clear();
    outer_storage_.resize(outer_total);
    hess_storage_.assign(hess_total, S(0));
    rebind_only(v);
  }

  // Re-bind without resizing.
  void rebind_only(std::vector<value_type>& v) {
    if (n_cols_ == 0) return;
    assert(static_cast<unsigned>(v.size()) == n_rows_
           && "tangent_slab<dual2nd>::rebind_only: vector size changed");
    outer_inner_t* outer_base = outer_base_();
    S*             hess_base  = hess_base_();
    for (unsigned i = 0; i < n_rows_; ++i) {
      auto& di = v[i];
      di.rebind_storage(outer_base + static_cast<std::size_t>(i) * n_cols_,
                        n_cols_);
      for (unsigned k = 0; k < n_cols_; ++k) {
        outer_inner_t& ok = outer_base[static_cast<std::size_t>(i) * n_cols_ + k];
        ok.rebind_storage(hess_base
                            + static_cast<std::size_t>(i) * n_cols_ * n_cols_
                            + static_cast<std::size_t>(k) * n_cols_,
                          n_cols_);
      }
    }
  }

private:
  outer_inner_t* outer_base_() noexcept {
    return external_outer_ ? external_outer_ : outer_storage_.data();
  }
  const outer_inner_t* outer_base_() const noexcept {
    return external_outer_ ? external_outer_ : outer_storage_.data();
  }
  S* hess_base_() noexcept {
    return external_hess_ ? external_hess_ : hess_storage_.data();
  }
  const S* hess_base_() const noexcept {
    return external_hess_ ? external_hess_ : hess_storage_.data();
  }

  std::vector<outer_inner_t> outer_storage_;
  std::vector<S>             hess_storage_;
  outer_inner_t*             external_outer_ = nullptr;
  S*                         external_hess_  = nullptr;
  unsigned n_rows_ = 0;
  unsigned n_cols_ = 0;
};

} // namespace detail

// =============================================================================
//  Slab-aware AXPY / SCAL on std::vector<dual<S, 0>>.
//
//  With the slab primed the tangent buffers are slices of one contiguous block,
//  so the tangent half collapses into a single BLAS call instead of n
//  expression evaluations. Without a slab these forward to vec_axpy / vec_scale.
// =============================================================================

// `with_hess = false` leaves the dual2nd Hessian layer untouched, for a
// caller that combines many slots into one BLAS-3 call.  No effect for
// any other T.
template<class T>
inline void vec_axpy_with_slab(
    std::vector<T>& y, detail::tangent_slab<T>& y_slab,
    double alpha,
    const std::vector<T>& x, const detail::tangent_slab<T>& x_slab,
    bool with_hess = true)
{
  if constexpr (detail::is_dual2nd<T>::value) {
    // BLAS-3 hybrid: scalar and inline d1 per element, one daxpy over the
    // Hessian block. arm_full restores the outer depend_ flags, which a
    // base::operator=(scalar) may have reset.
    using S_inner = typename detail::tangent_slab<T>::inner_type;
    const std::size_t n = y.size();
    if (y_slab.primed() && x_slab.primed()) {
      const S_inner a_s = static_cast<S_inner>(alpha);
      const unsigned m = y_slab.n_cols();
      for (std::size_t i = 0; i < n; ++i) {
        y[i].arm_full(m);  // ensure depend_ is true throughout
        y[i].scalar() += a_s * x[i].scalar();
        for (unsigned k = 0; k < m; ++k)
          y[i].d1_at(k) += a_s * x[i].d1_at(k);
      }
      const std::size_t h_total = with_hess ? y_slab.hess_size() : 0;
      if (h_total > 0) {
        if constexpr (std::is_same_v<S_inner, double>) {
          int len = static_cast<int>(h_total);
          int inc = 1;
          double a = alpha;
          F77_CALL(daxpy)(&len, &a,
                          const_cast<double*>(x_slab.hess_data()), &inc,
                          y_slab.hess_data(), &inc);
        } else {
          S_inner* yh = y_slab.hess_data();
          const S_inner* xh = x_slab.hess_data();
          for (std::size_t k = 0; k < h_total; ++k) yh[k] += a_s * xh[k];
        }
      }
    } else {
      const T a_t = T(static_cast<S_inner>(alpha));
      for (std::size_t i = 0; i < n; ++i) y[i] += a_t * x[i];
    }
    return;
  }
  if constexpr (detail::is_dynamic_dual<T>::value) {
    using S = typename T::value_type;
    const std::size_t n = y.size();
    if (y_slab.primed() && x_slab.primed()) {
      // BLAS hot path: per-element val + one BLAS pass over tangent block.
      // The two passes are disjoint (val ≠ tangents), no double-counting.
      for (std::size_t i = 0; i < n; ++i)
        y[i].x() += static_cast<S>(alpha) * x[i].x();
      assert(y_slab.tangent_size() == x_slab.tangent_size()
             && "vec_axpy_with_slab: slab size mismatch");
      const std::size_t total = y_slab.tangent_size();
      if (total > 0) {
        if constexpr (std::is_same_v<S, double>) {
          int len = static_cast<int>(total);
          int inc = 1;
          double a = alpha;
          F77_CALL(daxpy)(&len, &a,
                          const_cast<double*>(x_slab.tangent_data()), &inc,
                          y_slab.tangent_data(), &inc);
        } else {
          S* yp = y_slab.tangent_data();
          const S* xp = x_slab.tangent_data();
          const S a_s = static_cast<S>(alpha);
          for (std::size_t k = 0; k < total; ++k)
            yp[k] += a_s * xp[k];
        }
      }
    } else {
      // Slab unprimed (sensitivities-off solve, e.g. M=0 reparam path):
      // per-element ET handles val AND tangents in one shot. Wrap alpha
      // in a dual so dual::operator+=(Expr<>) fires.
      const T a_t = T(static_cast<S>(alpha));
      for (std::size_t i = 0; i < n; ++i) y[i] += a_t * x[i];
    }
  } else {
    // Non-dual scalar (double): ordinary vec_axpy.
    vec_axpy(y, alpha, x);
  }
}

// Slab-aware vector zero. It must not assign T(0) to the elements: that runs
// dual<T,0>::operator=(const U&), which drops tan_ and the slab binding with
// it. Without a slab it multiplies by zero, keeping any allocated buffer.
template<class T>
inline void vec_zero_with_slab(
    std::vector<T>& y, detail::tangent_slab<T>& y_slab)
{
  if constexpr (detail::is_dual2nd<T>::value) {
    // BLAS-3 hybrid: per-element scalar + inline d1 zero + memset on hess
    // slab block. (val_tan_block dropped; LU reads inline_d1 via
    // first_order_view.)
    using S_inner = typename detail::tangent_slab<T>::inner_type;
    const std::size_t n = y.size();
    if (y_slab.primed()) {
      const unsigned m = y_slab.n_cols();
      for (std::size_t i = 0; i < n; ++i) {
        y[i].arm_full(m);
        y[i].scalar() = S_inner(0);
        for (unsigned k = 0; k < m; ++k) y[i].d1_at(k) = S_inner(0);
      }
      const std::size_t h_total = y_slab.hess_size();
      if (h_total > 0)
        std::memset(y_slab.hess_data(), 0, h_total * sizeof(S_inner));
    } else {
      const S_inner zero = S_inner(0);
      for (std::size_t i = 0; i < n; ++i) y[i] *= zero;
    }
    return;
  }
  if constexpr (detail::is_dynamic_dual<T>::value) {
    using S = typename T::value_type;
    const std::size_t n = y.size();
    for (std::size_t i = 0; i < n; ++i) y[i].x() = S(0);
    // memset only over trivially copyable tangents: a dual tangent element
    // carries a tan_ pointer that must survive.
    if constexpr (std::is_trivially_copyable_v<S>) {
      if (y_slab.primed()) {
        const std::size_t total = y_slab.tangent_size();
        if (total > 0)
          std::memset(y_slab.tangent_data(), 0, total * sizeof(S));
        return;
      }
    }
    // Slab unprimed: zero values via *=0 to preserve any arena tan_
    // buffers (so they stay live for the next vec_axpy materialisation).
    const S zero = S(0);
    for (std::size_t i = 0; i < n; ++i) y[i] *= zero;
  } else {
    vec_zero(y);
  }
}

// Slab-aware vector copy: replaces the per-element dual::operator=
// loop (which iterates val_ + N tangents per dual) with a flat
// std::memcpy over the contiguous slab block. Falls back to plain
// per-element assignment for non-dynamic-dual T or unprimed slabs.
template<class T>
inline void vec_copy_with_slab(
    std::vector<T>& y, detail::tangent_slab<T>& y_slab,
    const std::vector<T>& x, const detail::tangent_slab<T>& x_slab)
{
  if constexpr (detail::is_dual2nd<T>::value) {
    // BLAS-3 hybrid: per-element scalar + inline d1 copy + memcpy on hess.
    // (sync_d1_redundant dropped: LU reads inline_d1 via first_order_view.)
    using S_inner = typename detail::tangent_slab<T>::inner_type;
    const std::size_t n = y.size();
    if (y_slab.primed() && x_slab.primed()) {
      const unsigned m = y_slab.n_cols();
      for (std::size_t i = 0; i < n; ++i) {
        y[i].arm_full(m);
        y[i].scalar() = x[i].scalar();
        for (unsigned k = 0; k < m; ++k) y[i].d1_at(k) = x[i].d1_at(k);
      }
      const std::size_t h_total = y_slab.hess_size();
      if (h_total > 0)
        std::memcpy(y_slab.hess_data(), x_slab.hess_data(),
                    h_total * sizeof(S_inner));
    } else {
      for (std::size_t i = 0; i < n; ++i) y[i] = x[i];
    }
    return;
  }
  if constexpr (detail::is_dynamic_dual<T>::value) {
    using S = typename T::value_type;
    const std::size_t n = y.size();
    for (std::size_t i = 0; i < n; ++i)
      y[i].x() = x[i].x();
    // memcpy only over trivially copyable tangents: a dual tangent element
    // carries a tan_ pointer that must not be aliased into y.
    if constexpr (std::is_trivially_copyable_v<S>) {
      if (y_slab.primed() && x_slab.primed()) {
        assert(y_slab.tangent_size() == x_slab.tangent_size()
               && "vec_copy_with_slab: slab size mismatch");
        const std::size_t total = y_slab.tangent_size();
        if (total > 0)
          std::memcpy(y_slab.tangent_data(), x_slab.tangent_data(),
                      total * sizeof(S));
        return;
      }
    }
    // Slab unprimed: per-element copy through dual::operator=.
    for (std::size_t i = 0; i < n; ++i) y[i] = x[i];
  } else {
    // Non-dual / static-N: std::vector copy assignment is fine.
    y = x;
  }
}

template<class T>
inline void vec_scale_with_slab(
    std::vector<T>& y, detail::tangent_slab<T>& y_slab, double alpha)
{
  if constexpr (detail::is_dual2nd<T>::value) {
    // BLAS-3 hybrid: per-element scalar + inline d1 scale + dscal on hess
    // block. (sync_d1_redundant dropped: LU reads inline_d1 via
    // first_order_view.)
    using S_inner = typename detail::tangent_slab<T>::inner_type;
    const std::size_t n = y.size();
    if (y_slab.primed()) {
      const S_inner a_s = static_cast<S_inner>(alpha);
      const unsigned m = y_slab.n_cols();
      for (std::size_t i = 0; i < n; ++i) {
        y[i].arm_full(m);
        y[i].scalar() *= a_s;
        for (unsigned k = 0; k < m; ++k) y[i].d1_at(k) *= a_s;
      }
      const std::size_t h_total = y_slab.hess_size();
      if (h_total > 0) {
        if constexpr (std::is_same_v<S_inner, double>) {
          int len = static_cast<int>(h_total);
          int inc = 1;
          double a = alpha;
          F77_CALL(dscal)(&len, &a, y_slab.hess_data(), &inc);
        } else {
          S_inner* yh = y_slab.hess_data();
          for (std::size_t k = 0; k < h_total; ++k) yh[k] *= a_s;
        }
      }
    } else {
      const S_inner a_s = static_cast<S_inner>(alpha);
      for (std::size_t i = 0; i < n; ++i) y[i] *= a_s;
    }
    return;
  }
  if constexpr (detail::is_dynamic_dual<T>::value) {
    using S = typename T::value_type;
    const std::size_t n = y.size();
    if (y_slab.primed()) {
      // BLAS hot path: per-element val * scalar + BLAS dscal on tangent block.
      for (std::size_t i = 0; i < n; ++i)
        y[i].x() *= static_cast<S>(alpha);
      const std::size_t total = y_slab.tangent_size();
      if (total > 0) {
        if constexpr (std::is_same_v<S, double>) {
          int len = static_cast<int>(total);
          int inc = 1;
          double a = alpha;
          F77_CALL(dscal)(&len, &a, y_slab.tangent_data(), &inc);
        } else {
          S* yp = y_slab.tangent_data();
          const S a_s = static_cast<S>(alpha);
          for (std::size_t k = 0; k < total; ++k)
            yp[k] *= a_s;
        }
      }
    } else {
      // Sensitivities off: per-element scalar *= goes through the inline
      // dual<T,0>::operator*=(U) path (no allocations).
      const S a_s = static_cast<S>(alpha);
      for (std::size_t i = 0; i < n; ++i) y[i] *= a_s;
    }
  } else {
    vec_scale(y, alpha);
  }
}

} // namespace cppde

#endif // CPPDE_DUAL_SLAB_HPP
