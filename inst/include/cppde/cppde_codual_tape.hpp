/*
 Operation tape for cppde::codual<T>, the reverse-mode AD type.

 One node per recorded operation, holding the local partial derivatives and the
 slots of its operands. Every expression reduces to unary and binary nodes, so
 the node is fixed width and the tape is a flat array.

 Two invariants carry the design. An operand slot is always smaller than the
 slot it feeds, so the reverse sweep is one descending loop and never a graph
 traversal. The tape is thread-local and rewound per step, so the memory bound
 is one step rather than one trajectory.

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
  static constexpr unsigned none = static_cast<unsigned>(-1);

  struct node {
    unsigned a, b;   // operand slots, `none` when absent
    T        pa, pb; // local partials with respect to those operands
  };

  // -- recording --------------------------------------------------------------

  unsigned independent() {
    nodes_.push_back(node{none, none, T(), T()});
    return static_cast<unsigned>(nodes_.size() - 1u);
  }

  unsigned record(unsigned a, const T& pa) {
    nodes_.push_back(node{a, none, pa, T()});
    return static_cast<unsigned>(nodes_.size() - 1u);
  }

  unsigned record(unsigned a, const T& pa, unsigned b, const T& pb) {
    nodes_.push_back(node{a, b, pa, pb});
    return static_cast<unsigned>(nodes_.size() - 1u);
  }

  std::size_t size() const { return nodes_.size(); }

  // -- sweeping ---------------------------------------------------------------

  // Clears the adjoints and sizes them to the tape. Call once before seeding.
  void prepare() { adj_.assign(nodes_.size(), T()); }

  // Adds w onto the adjoint of one slot. Repeated seeds accumulate.
  void seed(unsigned slot, const T& w) {
    if (slot == none) return;
    if (adj_.size() < nodes_.size()) adj_.resize(nodes_.size(), T());
    adj_[slot] = adj_[slot] + w;
  }

  // Single backwards pass. Operand slots are strictly smaller than the node
  // they feed, so one descending loop suffices.
  void reverse() {
    if (adj_.size() < nodes_.size()) adj_.resize(nodes_.size(), T());
    for (std::size_t i = nodes_.size(); i-- > 0;) {
      const T& w = adj_[i];
      if (w == T()) continue;
      const node& n = nodes_[i];
      if (n.a != none) adj_[n.a] = adj_[n.a] + n.pa * w;
      if (n.b != none) adj_[n.b] = adj_[n.b] + n.pb * w;
    }
  }

  const T& adjoint(unsigned slot) const { return adj_[slot]; }

  // -- lifetime ---------------------------------------------------------------

  // Drops the nodes and adjoints, keeps both capacities.
  void rewind() {
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
  std::vector<node> nodes_;
  std::vector<T>    adj_;
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
