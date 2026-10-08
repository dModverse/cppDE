/*
 Modes of the state switches of the right-hand side.

 A comparison of the model that reads a state is held as a mode, one boolean
 per switch, which the generated right-hand side, its Jacobian and every
 derivative code read instead of the comparison. The integration loops locate
 the root of the comparison's argument and set the mode there, see
 cppde_event_engine.hpp; between two roots f is smooth in x.

 The generated solve hands its modes over per thread, since a batch runs
 several solves at once. Plain types and free of R.

 Copyright (C) 2026 Simon Beyer
 */

#ifndef CPPDE_SWITCH_MODES_HPP
#define CPPDE_SWITCH_MODES_HPP

#include <cstddef>
#include <vector>

namespace cppde {

struct switch_modes {
  signed char* m = nullptr;
  std::size_t  n = 0;
};

// Where the generated code reads the modes of the running solve. Trivially
// destructible, see cppde_tls.hpp.
inline switch_modes& switch_mode_sink() {
  thread_local switch_modes s;
  return s;
}

struct switch_mode_scope {
  switch_modes prev;
  switch_mode_scope(signed char* m, std::size_t n) : prev(switch_mode_sink()) {
    switch_mode_sink() = switch_modes{m, n};
  }
  ~switch_mode_scope() { switch_mode_sink() = prev; }
  switch_mode_scope(const switch_mode_scope&) = delete;
  switch_mode_scope& operator=(const switch_mode_scope&) = delete;
};

// Mode k of the running solve.
inline bool switch_mode(int k) { return switch_mode_sink().m[k] != 0; }

namespace detail {

// A mode holds where its argument is positive, or zero for a non-strict
// comparison.
inline bool mode_holds(double g, bool closed) {
  return g > 0.0 || (closed && g == 0.0);
}

inline void set_switch_mode(int k, bool on) {
  switch_mode_sink().m[k] = on ? 1 : 0;
}

// The modes as they stand, empty for a model without any.
inline std::vector<signed char> switch_mode_snapshot() {
  const switch_modes& s = switch_mode_sink();
  return std::vector<signed char>(s.m, s.m + s.n);
}

// Puts back a snapshot; an empty one leaves the modes alone.
inline void use_switch_modes(const std::vector<signed char>& v) {
  const switch_modes& s = switch_mode_sink();
  for (std::size_t i = 0; i < v.size() && i < s.n; ++i) s.m[i] = v[i];
}

}  // namespace detail
}  // namespace cppde

#endif  // CPPDE_SWITCH_MODES_HPP
