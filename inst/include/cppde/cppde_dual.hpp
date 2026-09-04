/*
 cppde::dual<T, N>: first-order forward-mode AD type.

 Both static-N (N > 0) and dynamic-N (N == 0) specs share a single storage
 strategy: `T* tan_;` pointing either at a row of an externally-owned
 cppde::detail::tangent_slab block (set via rebind_storage), or at a buffer
 obtained from the thread-local cppde::dual_arena. Per-RHS-eval allocation
 cost for temporaries is a bump-pointer increment; rolled back LIFO via a
 cppde::dual_arena::scope guard at the RHS body level. Storage layout in
 std::vector<dual<T, N>> is therefore SoA at the tangent level (rows of one
 contiguous slab), enabling BLAS daxpy/dscal across the tangent block.

 The only structural difference between the two specs is that the static-N
 primary template has no runtime size_ field: N is a compile-time constant
 used as the loop bound (constexpr-foldable for unrolling / SIMD).

 Provides the standard accessor surface used by codegen output:
 `.x()`, `.diff(idx)`, `.diff(idx, n)`, `.d(j)`, `.size()`, `.depend()`,
 `operator[]`.

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
// assignment from any ET node lets `dxdt[i] = a*b + c*d` materialise the whole
// tree into dxdt[i] without intermediate arena allocations.
namespace dual_expr {
  template<class Derived> struct Expr;
}

// -----------------------------------------------------------------------------
// Internal helper: allocate T[n] from the TLS arena.
//
// A trivially destructible and trivially default-constructible type gets a bare
// bump. A non-trivial default ctor still has to run, or a nested dual would find
// a stale tan_ and skip its own allocation. A real dtor adds dtor tracking.
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
// Primary template: static N, the compile-time tangent width (N == 0 below).
//
// tan_ points either into a tangent_slab row after rebind_storage, or into an
// arena buffer for temporaries. N is the loop bound, not inline storage, which
// is what lets this specialisation share the dynamic one's SoA layout.
// =============================================================================
template<class T = double, unsigned N = 0>
class dual {
  static_assert(N > 0, "primary template requires N > 0; N == 0 is specialized");

  T    val_;
  T*   tan_;
  bool depend_;

public:
  using value_type   = T;
  static constexpr unsigned static_size = N;

  // -- constructors -----------------------------------------------------------
  dual() : val_(), tan_(nullptr), depend_(false) {}

  template<class U,
           std::enable_if_t<std::is_convertible_v<U, T>, int> = 0>
  dual(const U& v) : val_(static_cast<T>(v)), tan_(nullptr), depend_(false) {}

  dual(const dual& o) : val_(o.val_), tan_(nullptr), depend_(false) {
    if (o.depend_) {
      depend_ = true;
      tan_    = detail::arena_alloc_t<T>(N);
      for (unsigned i = 0; i < N; ++i) tan_[i] = o.tan_[i];
    }
  }

  dual& operator=(const dual& o) {
    if (this == &o) return *this;
    val_ = o.val_;
    if (o.depend_) {
      if (tan_ == nullptr) tan_ = detail::arena_alloc_t<T>(N);
      depend_ = true;
      for (unsigned i = 0; i < N; ++i) tan_[i] = o.tan_[i];
    } else {
      // Non-depend source: zero the tangents but keep the binding, so a later
      // .diff() or ET assignment reuses the buffer. depend_ goes to false, as in
      // operator=(const U&).
      if (tan_ != nullptr) {
        for (unsigned i = 0; i < N; ++i) tan_[i] = T();
      }
      depend_ = false;
    }
    return *this;
  }

  template<class U,
           std::enable_if_t<std::is_convertible_v<U, T>, int> = 0>
  dual& operator=(const U& v) {
    val_ = static_cast<T>(v);
    // Keep the tan_ buffer and zero its values, as dual<T,0> does. Without it a
    // slab-bound dual would lose its binding on every `dual = scalar`, which is
    // what the codegen does when it re-seeds the state.
    if (tan_ != nullptr) {
      for (unsigned i = 0; i < N; ++i) tan_[i] = T();
    }
    depend_ = false;
    return *this;
  }

  // Move steals tan_ instead of allocating and copying, which is what makes
  // `dxdt[i] = a + b + c` cheap. Both moves are declared explicitly because the
  // copy assignment operator suppresses the implicit ones.
  dual(dual&& o) noexcept
    : val_(std::move(o.val_)), tan_(o.tan_), depend_(o.depend_)
  {
    o.tan_    = nullptr;
    o.depend_ = false;
  }
  dual& operator=(dual&& o) noexcept {
    if (this == &o) return *this;
    val_ = std::move(o.val_);
    if (o.depend_) {
      if (tan_ == nullptr) {
        // Steal: cheap path, no allocation, no copy.
        tan_      = o.tan_;
        depend_   = true;
        o.tan_    = nullptr;
        o.depend_ = false;
      } else {
        // Pre-existing tangent buffer (likely slab-bound): copy values into
        // it. This preserves address stability for callers that hold
        // references to our tangent storage (slab is the canonical case).
        for (unsigned i = 0; i < N; ++i) tan_[i] = std::move(o.tan_[i]);
        depend_ = true;
      }
    } else if (tan_ != nullptr) {
      // RHS is non-depend, we have a buffer: zero the tangents (mirror
      // copy-assignment behaviour to avoid stale derivative values).
      for (unsigned i = 0; i < N; ++i) tan_[i] = T();
      depend_ = false;
    } else {
      depend_ = false;
    }
    return *this;
  }

  // -- accessors --------------------------------------------------------------
  const T& x()   const { return val_; }
  T&       x()         { return val_; }
  const T& val() const { return val_; }   // .x() alias
  T&       val()       { return val_; }

  unsigned size()   const { return depend_ ? N : 0u; }
  bool     depend() const { return depend_; }
  // Loop bound for use INSIDE arithmetic operators after set_depend_from*().
  // Returns the compile-time constant N (constexpr-foldable for unrolling /
  // SIMD vectorisation). The runtime-checked size() is for callers that
  // need to honour the depend_ state.
  static constexpr unsigned loop_size() { return N; }

  // Returns mutable reference, even when out-of-bounds or !depend
  // (returns ref to a thread-local zero). Codegen sometimes does
  // `dual.d(j) = value` after .diff(idx) seeding.
  T& d(unsigned j) {
    if (!depend_ || j >= N) {
      thread_local T zero{};
      zero = T();
      return zero;
    }
    return tan_[j];
  }
  const T& d(unsigned j) const {
    if (!depend_ || j >= N) {
      thread_local T zero{};
      zero = T();
      return zero;
    }
    return tan_[j];
  }

  // Direct tangent access (caller must have called .diff() / set_depend_size()
  // first to ensure tan_ is allocated; otherwise this dereferences nullptr).
  const T& operator[](unsigned j) const { return tan_[j]; }
  T&       operator[](unsigned j)       { return tan_[j]; }

  // -- seeding ---------------------------------------------------------------
  // Activate this dual as the idx-th independent variable: tan_[idx]=1, rest=0.
  T& diff(unsigned idx) {
    assert(idx < N && "diff(idx): idx out of compile-time bound N");
    if (tan_ == nullptr) tan_ = detail::arena_alloc_t<T>(N);
    depend_ = true;
    for (unsigned i = 0; i < N; ++i) tan_[i] = T();
    tan_[idx] = T(1);
    return tan_[idx];
  }
  // 2-arg form (kept for API uniformity with the dynamic spec): n must
  // equal N at compile time.
  T& diff(unsigned idx, unsigned n) {
    assert(n == N && "diff(idx, n): n must equal compile-time N");
    (void)n;
    return diff(idx);
  }

  // Internal: mark this dual as having an active tangent vector without
  // initializing values (operators / ET fill the tangents themselves). If
  // tan_ has not been bound yet, allocate a fresh buffer from the arena.
  void set_depend() {
    if (tan_ == nullptr) tan_ = detail::arena_alloc_t<T>(N);
    depend_ = true;
  }

  // Allocate-without-init helper used by the ET path (mirrors the dual<T,0>
  // member of the same name). The 0-arg form is the natural one for static-N
  // (size is the template parameter); the 1-arg form keeps API parity with
  // the dynamic spec for shared ET callers: n must equal N.
  void set_depend_size() { set_depend(); }
  void set_depend_size(unsigned n) {
    assert(n == N && "dual<T,N>::set_depend_size(n): n must equal compile-time N");
    (void)n;
    set_depend();
  }

  // Non-allocating bind onto an externally owned buffer, typically a row of
  // tangent_slab. The dual never frees tan_, so the owner has to outlive it.
  // depend_ goes to true: a slab-bound dual is always active.
  void rebind_storage(T* p, unsigned n) noexcept {
    assert(n == N && "dual<T,N>::rebind_storage: n must equal compile-time N");
    (void)n;
    tan_    = p;
    depend_ = true;
  }

  void set_depend_from(const dual&)              { set_depend(); }
  void set_depend_from(const dual&, const dual&) { set_depend(); }

  // Compound assignment operators (defined inline via free + / - / * / /).
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
    if (depend_) for (unsigned i = 0; i < N; ++i) tan_[i] *= static_cast<T>(v);
    return *this;
  }
  template<class U,
           std::enable_if_t<std::is_arithmetic_v<U>, int> = 0>
  dual& operator/=(const U& v) {
    val_ /= static_cast<T>(v);
    if (depend_) for (unsigned i = 0; i < N; ++i) tan_[i] /= static_cast<T>(v);
    return *this;
  }

  // -- expression-template assignment / construction --------------------------
  // Materialises an Expr<D> tree straight into tan_: value first, then one
  // allocation if not yet bound, then the tangent loop over the constant N.
  // Definitions live in cppde_dual_expr.hpp, once Expr<D> is complete.
  template<class D>
  dual& operator=(const dual_expr::Expr<D>& e);

  template<class D>
  dual(const dual_expr::Expr<D>& e);

  template<class D>
  dual& operator+=(const dual_expr::Expr<D>& e);
  template<class D>
  dual& operator-=(const dual_expr::Expr<D>& e);
  template<class D>
  dual& operator*=(const dual_expr::Expr<D>& e);
  template<class D>
  dual& operator/=(const dual_expr::Expr<D>& e);
};

// =============================================================================
// Partial specialization: dynamic N (arena-backed storage).
// =============================================================================
template<class T>
class dual<T, 0> {
  T*       tan_;
  unsigned size_;
  T        val_;

public:
  using value_type   = T;
  static constexpr unsigned static_size = 0;

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
      assert(size_ == o.size_ && "dual<T,0>: tangent-vector size mismatch");
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
    // be equivalent downstream, but a slab-bound dual would lose its binding on
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
        assert(size_ == o.size_ && "dual<T,0>: tangent-vector size mismatch in move-assign");
        for (unsigned i = 0; i < size_; ++i) tan_[i] = std::move(o.tan_[i]);
      }
    } else if (size_ > 0) {
      // RHS is non-depend, we have a buffer: zero the tangents (mirror
      // the copy-assignment behaviour to avoid stale derivative values).
      for (unsigned i = 0; i < size_; ++i) tan_[i] = T();
    }
    return *this;
  }

  // -- expression-template assignment / construction --------------------------
  // Materialises an Expr<D> tree straight into this dual: one scalar pass for
  // the value, one arena allocation or in-place reuse, one fused tangent loop.
  // Definitions live in cppde_dual_expr.hpp, once Expr<D> is complete.
  template<class D>
  dual& operator=(const dual_expr::Expr<D>& e);

  template<class D>
  dual(const dual_expr::Expr<D>& e);

  // In-place compound assignment from an Expr<D>. Without these overloads the
  // compiler routes through operator+=(const dual&) and synthesises a temporary
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
  // Mirrors dual<T,N>::loop_size() for the dynamic spec: see comment there.
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
      assert(size_ == n && "dual<T,0>: tangent size mismatch");
    }
  }

  // Non-allocating bind onto an externally owned buffer, typically a row of
  // tangent_slab. As in the arena-backed case the dual never frees tan_, so the
  // owner has to keep it alive.
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
template<class T, unsigned N>
inline auto value_of(const dual<T, N>& d) { return value_of(d.x()); }

// dual2nd<T, N> is now a distinct class (cppde_dual2nd.hpp): a public-
// inheritance refinement of dual<dual<T, N>, N> with hand-derived symmetric
// math primitives. The previous typedef alias has been removed; downstream
// code should #include <cppde/cppde_dual2nd.hpp> when needed.

} // namespace cppde

#endif // CPPDE_DUAL_HPP
