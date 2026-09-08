/*
 Operation tape for cppde::codual<T>, the reverse-mode AD type.

 One node per recorded operation, holding the local partial derivatives and the
 slots of its operands. Every expression reduces to unary and binary nodes, so
 the node is fixed width and the tape is a flat array.

 Three invariants carry the design. An operand slot is always smaller than the
 slot it feeds, so the reverse sweep is one descending loop and never a graph
 traversal. The tape is thread-local and rewound per step, so the memory bound
 is one step rather than one trajectory. And slots are monotone across rewinds:
 the base advances instead of the indices being reused, so a value left over in
 a reused buffer names a slot the tape no longer owns and reads as a constant.
 Without that, a buffer the replay does not write would carry a dependence on an
 unrelated node, which is a wrong derivative with a right value.

 Storage follows cppde_tls.hpp: a thread_local pointer, never a thread_local
 object, whose destructor would register __cxa_thread_atexit and make
 dyn.unload() leak every rebuilt model .so.

 Copyright (C) 2026 Simon Beyer
 */

#ifndef CPPDE_CODUAL_TAPE_HPP
#define CPPDE_CODUAL_TAPE_HPP

#include <cstddef>
#include <vector>

namespace cppde {

// =============================================================================
//  codual_tape<T>
//
//  Slots index both the node array and the adjoint array; independent variables
//  get a node with no operands so the two stay in step.
// =============================================================================
template<class T = double>
class codual_tape {
public:
  using value_type = T;

  // Slot of a value that carries no dependence. Operations against it skip the
  // corresponding accumulation instead of adding a zero.
  static constexpr std::size_t none = static_cast<std::size_t>(-1);

  struct node {
    unsigned a, b;   // operand positions in this tape, `nolocal` when absent
    T        pa, pb; // local partials with respect to those operands
  };

  static constexpr unsigned nolocal = static_cast<unsigned>(-1);

  // Whether a slot belongs to the tape as it stands. Everything older is a
  // constant, which is what makes a reused buffer safe without clearing it.
  bool live(std::size_t slot) const {
    return slot != none && slot >= base_ && slot - base_ < nodes_.size();
  }

  // -- recording --------------------------------------------------------------

  std::size_t independent() {
    nodes_.push_back(node{nolocal, nolocal, T(), T()});
    return base_ + nodes_.size() - 1u;
  }

  std::size_t record(std::size_t a, const T& pa) {
    nodes_.push_back(node{local(a), nolocal, pa, T()});
    return base_ + nodes_.size() - 1u;
  }

  std::size_t record(std::size_t a, const T& pa, std::size_t b, const T& pb) {
    nodes_.push_back(node{local(a), local(b), pa, pb});
    return base_ + nodes_.size() - 1u;
  }

  std::size_t size() const { return nodes_.size(); }

  // -- sweeping ---------------------------------------------------------------

  // Clears the adjoints and sizes them to the tape. Call once before seeding.
  void prepare() { adj_.assign(nodes_.size(), T()); }

  // Adds w onto the adjoint of one slot. Repeated seeds accumulate.
  void seed(std::size_t slot, const T& w) {
    if (!live(slot)) return;
    if (adj_.size() < nodes_.size()) adj_.resize(nodes_.size(), T());
    adj_[slot - base_] = adj_[slot - base_] + w;
  }

  // Single backwards pass. Operand slots are strictly smaller than the node
  // they feed, so one descending loop suffices.
  void reverse() { reverse(nodes_.size(), 0); }

  // The nodes in [lo, hi), newest first. An implicit equation interrupts the
  // sweep at its own nodes: what its inputs receive is not a chain rule but a
  // transposed solve, which the caller does between the two calls. Everything
  // that reads the solution has to be recorded above hi, or its share of the
  // solution's cotangent is not there yet when the solve runs.
  void reverse(std::size_t hi, std::size_t lo) {
    if (adj_.size() < nodes_.size()) adj_.resize(nodes_.size(), T());
    if (hi > nodes_.size()) hi = nodes_.size();
    for (std::size_t i = hi; i-- > lo;) {
      const T& w = adj_[i];
      if (w == T()) continue;
      const node& n = nodes_[i];
      if (n.a != nolocal) adj_[n.a] = adj_[n.a] + n.pa * w;
      if (n.b != nolocal) adj_[n.b] = adj_[n.b] + n.pb * w;
    }
  }

  T adjoint(std::size_t slot) const {
    return live(slot) ? adj_[slot - base_] : T();
  }

  // -- lifetime ---------------------------------------------------------------

  // Drops the nodes and adjoints, keeps both capacities, and moves the slot base
  // past everything just dropped so no old slot can name a new node.
  void rewind() {
    base_ += nodes_.size();
    nodes_.clear();
    adj_.clear();
  }

  // RAII rewind to the state at construction, LIFO like dual_arena::scope.
  class scope {
    codual_tape& t_;
    std::size_t  mark_;
  public:
    explicit scope(codual_tape& t) : t_(t), mark_(t.nodes_.size()) {}
    scope(const scope&)            = delete;
    scope& operator=(const scope&) = delete;
    ~scope() {
      t_.nodes_.resize(mark_);
      if (t_.adj_.size() > mark_) t_.adj_.resize(mark_);
    }
  };

private:
  unsigned local(std::size_t slot) const {
    return live(slot) ? static_cast<unsigned>(slot - base_) : nolocal;
  }

  std::vector<node> nodes_;
  std::vector<T>    adj_;
  std::size_t       base_ = 0;
};

// The tape every codual<T> in this thread records onto.
template<class T>
inline codual_tape<T>& codual_tape_for() {
  thread_local codual_tape<T>* p = nullptr;
  if (p == nullptr) p = new codual_tape<T>();  // leaked on purpose, see header
  return *p;
}

}  // namespace cppde

#endif  // CPPDE_CODUAL_TAPE_HPP
