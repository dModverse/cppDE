/*
 Thread-local scratch storage.

 A thread_local with a non-trivial destructor registers __cxa_thread_atexit,
 which makes glibc refuse to unmap the shared object -- dyn.unload() becomes a
 no-op and every model rebuild leaks its .so.  A thread_local pointer registers
 nothing; the buffer behind it is leaked on purpose (one per thread per slot).

 Copyright (C) 2026 Simon Beyer
 */

#ifndef CPPDE_TLS_HPP
#define CPPDE_TLS_HPP

#include <vector>

namespace cppde {
namespace detail {

// Per-thread scratch. Slot must be distinct per call site; callers assign()
// or resize() before use.
template <int Slot>
inline std::vector<double>& tls_scratch_f64() {
  thread_local std::vector<double>* p = nullptr;
  if (p == nullptr) p = new std::vector<double>();  // leaked, see above
  return *p;
}

}  // namespace detail
}  // namespace cppde

#endif  // CPPDE_TLS_HPP
