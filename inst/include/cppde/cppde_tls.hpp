/*
 Thread-local scratch storage, kept in the dual arena of the thread: a solve's
 outermost arena scope frees it, and no thread_local with a destructor pins the
 shared object against dyn.unload().

 Copyright (C) 2026 Simon Beyer
 */

#ifndef CPPDE_TLS_HPP
#define CPPDE_TLS_HPP

#include <vector>

#include "cppde_dual_arena.hpp"

namespace cppde {
namespace detail {

// Per-thread scratch. Slot must be distinct per call site; callers assign()
// or resize() before use, inside a dual_arena::scope.
template <int Slot>
inline std::vector<double>& tls_scratch_f64() {
  return dual_arena::arena().template scratch<Slot>();
}

}  // namespace detail
}  // namespace cppde

#endif  // CPPDE_TLS_HPP
