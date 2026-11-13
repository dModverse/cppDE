/*
 cppde::dual<T>: first-order forward-mode AD type with a runtime tangent width.

 Storage is `T* tan_;` pointing either at a row of an externally-owned
 cppde::detail::tangent_slab block (set via rebind_storage), or at a buffer
 obtained from the thread-local cppde::dual_arena. Per-RHS-eval allocation
 cost for temporaries is a bump-pointer increment; rolled back LIFO via a
 cppde::dual_arena::scope guard at the RHS body level. Storage layout in
 std::vector<dual<T>> is therefore SoA at the tangent level (rows of one
 contiguous slab), enabling BLAS daxpy/dscal across the tangent block.

 Provides the standard accessor surface used by codegen output:
 `.x()`, `.diff(idx, n)`, `.d(j)`, `.size()`, `.depend()`, `operator[]`.

 This header defines the data class only. Arithmetic operators and math
 functions (eager path, gated to nested-AD T only) live in
 cppde_dual_math.hpp; the expression-template overlay (active for non-AD T)
 lives in cppde_dual_expr.hpp.

 Copyright (C) 2026 Simon Beyer
 */

#ifndef CPPDE_DUAL_HPP
#define CPPDE_DUAL_HPP

#include <cppde/cppde_dual_arena.hpp>
#include <cppde/cppde_scalar_ops.hpp>

#include <cassert>
#include <cstddef>
#include <type_traits>

namespace cppde {

// Forward declaration of the ET CRTP base (cppde_dual_expr.hpp). A templated
// assignment from any ET node lets `dxdt[i] = a*b + c*d` materialize the whole
// tree into dxdt[i] without intermediate arena allocations.
namespace dual_expr {
  template<class Derived> struct Expr;
}

// -----------------------------------------------------------------------------
// Internal helper: allocate T[n] from the TLS arena.
//
// A trivially destructible and trivially default-constructible type gets a bare
// bump. A non-trivial default ctor still has to run, or a nested dual would find
// a stale tan_ and skip the allocation. A real dtor adds dtor tracking.
namespace detail {

template<class T>
inline T* arena_alloc_t(std::size_t n) {
  if constexpr (std::is_trivially_destructible_v<T>
                && std::is_trivially_default_constructible_v<T>) {
    return dual_arena::arena().alloc_trivial<T>(n);
  } else if constexpr (std::is_trivially_destructible_v<T>) {
    T* arr = dual_arena::arena().alloc_trivial<T>(n);
    for (std::size_t i = 0; i < n; ++i) ::new (static_cast<void*>(arr + i)) T();
    return arr;
  } else {
    return dual_arena::arena().template alloc<T>(n);
  }
}

} // namespace detail

// =============================================================================
// dual<T>: the tangent width is a runtime size_, set by the first diff(),
// set_depend_size() or rebind_storage().
// =============================================================================
template<class T = double>
class dual {
  T*       tan_;
  unsigned size_;
  T        val_;

public:
  using value_type   = T;

  // -- constructors -----------------------------------------------------------
  dual() : tan_(nullptr), size_(0), val_() {}

  template<class U,
           std::enable_if_t<std::is_convertible_v<U, T>, int> = 0>
  dual(const U& v) : tan_(nullptr), size_(0), val_(static_cast<T>(v)) {}

  dual(const dual& o) : tan_(nullptr), size_(0), val_(o.val_) {
    if (o.size_ > 0) {
      size_ = o.size_;
      tan_  = detail::arena_alloc_t<T>(size_);
      for (unsigned i = 0; i < size_; ++i) tan_[i] = o.tan_[i];
    }
  }

  dual& operator=(const dual& o) {
    if (this == &o) return *this;
    val_ = o.val_;
    if (o.size_ > 0) {
      if (size_ == 0) {
        size_ = o.size_;
        tan_  = detail::arena_alloc_t<T>(size_);
      }
      assert(size_ == o.size_ && "dual: tangent-vector size mismatch");
      for (unsigned i = 0; i < size_; ++i) tan_[i] = o.tan_[i];
    } else if (size_ > 0) {
      for (unsigned i = 0; i < size_; ++i) tan_[i] = T();
    }
    return *this;
  }

  template<class U,
           std::enable_if_t<std::is_convertible_v<U, T>, int> = 0>
  dual& operator=(const U& v) {
    val_ = static_cast<T>(v);
    // Keep the buffer and zero the tangent values. Dropping to size_ = 0 would
    // be equivalent downstream, but a slab-bound dual would lose the binding on
    // every `dual = scalar`, and the bound paths take the BLAS route anyway.
    if (size_ > 0) {
      for (unsigned i = 0; i < size_; ++i) tan_[i] = T();
    }
    return *this;
  }

  // Move steals the arena pointer instead of allocating and copying, one
  // allocation per top-level expression. Both moves are declared explicitly
  // because the copy assignment operator suppresses the implicit ones.
  dual(dual&& o) noexcept : tan_(o.tan_), size_(o.size_), val_(std::move(o.val_)) {
    o.tan_ = nullptr;
    o.size_ = 0;
  }
  dual& operator=(dual&& o) noexcept {
    if (this == &o) return *this;
    val_ = std::move(o.val_);
    if (o.size_ > 0) {
      if (size_ == 0) {
        // Steal: cheap path, no allocation, no copy.
        tan_   = o.tan_;
        size_  = o.size_;
        o.tan_  = nullptr;
        o.size_ = 0;
      } else {
        // Pre-existing tangent buffer: copy values into it (preserves
        // address stability for callers that hold references).
        assert(size_ == o.size_ && "dual: tangent-vector size mismatch in move-assign");
        for (unsigned i = 0; i < size_; ++i) tan_[i] = std::move(o.tan_[i]);
      }
    } else if (size_ > 0) {
      // A non-depend right-hand side zeroes the tangents of an existing buffer,
      // as copy assignment does, so no stale derivative survives.
      for (unsigned i = 0; i < size_; ++i) tan_[i] = T();
    }
    return *this;
  }

  // -- expression-template assignment / construction --------------------------
  // Materializes an Expr<D> tree straight into this dual: one scalar pass for
  // the value, one arena allocation or in-place reuse, one fused tangent loop.
  // Definitions live in cppde_dual_expr.hpp, once Expr<D> is complete.
  template<class D>
  dual& operator=(const dual_expr::Expr<D>& e);

  template<class D>
  dual(const dual_expr::Expr<D>& e);

  // In-place compound assignment from an Expr<D>. Without these overloads the
  // compiler routes through operator+=(const dual&) and synthesizes a temporary
  // that allocates a tangent buffer per element, the dominant cost in axpy.
  template<class D>
  dual& operator+=(const dual_expr::Expr<D>& e);
  template<class D>
  dual& operator-=(const dual_expr::Expr<D>& e);
  template<class D>
  dual& operator*=(const dual_expr::Expr<D>& e);
  template<class D>
  dual& operator/=(const dual_expr::Expr<D>& e);

  // -- accessors --------------------------------------------------------------
  const T& x()   const { return val_; }
  T&       x()         { return val_; }
  const T& val() const { return val_; }
  T&       val()       { return val_; }

  unsigned size()   const { return size_; }
  bool     depend() const { return size_ != 0; }
  // Loop bound inside the arithmetic operators, after set_depend_from*().
  unsigned loop_size() const { return size_; }

  // Mutable + const overloads, both returning a reference.
  // Out-of-bounds returns ref to a thread-local zero.
  T& d(unsigned j) {
    if (j >= size_) {
      thread_local T zero{};
      zero = T();
      return zero;
    }
    return tan_[j];
  }
  const T& d(unsigned j) const {
    if (j >= size_) {
      thread_local T zero{};
      zero = T();
      return zero;
    }
    return tan_[j];
  }

  const T& operator[](unsigned j) const { return tan_[j]; }
  T&       operator[](unsigned j)       { return tan_[j]; }

  // -- seeding ----------------------------------------------------------------
  T& diff(unsigned idx, unsigned n) {
    assert(idx < n && "diff(idx, n): idx out of bounds");
    if (size_ == 0) {
      size_ = n;
      tan_  = detail::arena_alloc_t<T>(n);
    } else {
      assert(size_ == n && "diff(idx, n): n must match existing size");
    }
    for (unsigned i = 0; i < size_; ++i) tan_[i] = T();
    tan_[idx] = T(1);
    return tan_[idx];
  }

  // Allocate-without-init helper used by operators (fills tangents themselves).
  void set_depend_size(unsigned n) {
    if (size_ == 0) {
      size_ = n;
      tan_  = detail::arena_alloc_t<T>(n);
    } else {
      assert(size_ == n && "dual: tangent size mismatch");
    }
  }

  // Bind a zeroed tangent buffer now, without claiming a dependence, so a write
  // inside a callee's own dual_arena::scope reuses it instead of allocating there.
  // A dual with no width set takes the caller's declared one, if any.
  void arm() {
    if (size_ == 0) {
      const unsigned w = dual_arena::default_tangent_width();
      if (w == 0) return;
      size_ = w;
      tan_  = detail::arena_alloc_t<T>(size_);
    }
    if (tan_ == nullptr) tan_ = detail::arena_alloc_t<T>(size_);
    for (unsigned i = 0; i < size_; ++i) tan_[i] = T();
  }

  // Non-allocating bind onto an externally owned buffer, typically a row of
  // tangent_slab. The dual never frees tan_, so the owner has to keep it alive.
  void rebind_storage(T* p, unsigned n) noexcept {
    tan_  = p;
    size_ = n;
  }
  void set_depend_from(const dual& a) {
    set_depend_size(a.size_);
  }
  void set_depend_from(const dual& a, const dual& b) {
    assert(a.size_ == b.size_ && "operands have different tangent sizes");
    set_depend_size(a.size_);
  }

  // Compound assignment operators (free + / - / * / / defined in math header).
  dual& operator+=(const dual& o);
  dual& operator-=(const dual& o);
  dual& operator*=(const dual& o);
  dual& operator/=(const dual& o);
  template<class U,
           std::enable_if_t<std::is_arithmetic_v<U>, int> = 0>
  dual& operator+=(const U& v) { val_ += static_cast<T>(v); return *this; }
  template<class U,
           std::enable_if_t<std::is_arithmetic_v<U>, int> = 0>
  dual& operator-=(const U& v) { val_ -= static_cast<T>(v); return *this; }
  template<class U,
           std::enable_if_t<std::is_arithmetic_v<U>, int> = 0>
  dual& operator*=(const U& v) {
    val_ *= static_cast<T>(v);
    if (size_) for (unsigned i = 0; i < size_; ++i) tan_[i] *= static_cast<T>(v);
    return *this;
  }
  template<class U,
           std::enable_if_t<std::is_arithmetic_v<U>, int> = 0>
  dual& operator/=(const U& v) {
    val_ /= static_cast<T>(v);
    if (size_) for (unsigned i = 0; i < size_; ++i) tan_[i] /= static_cast<T>(v);
    return *this;
  }
};

// =============================================================================
// value_of for the first-order dual: peel to the innermost scalar. Nested T
// recurses through the same overload set.
// =============================================================================
template<class T>
inline auto value_of(const dual<T>& d) { return value_of(d.x()); }

// dual2nd<T>, a refinement of dual<dual<T>>, lives in cppde_dual2nd.hpp.

} // namespace cppde

#endif // CPPDE_DUAL_HPP
