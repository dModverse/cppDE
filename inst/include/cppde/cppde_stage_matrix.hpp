/*
 stage_matrix<T, S>: contiguous tangent storage for S RK stages.

 Owns one [(S) × n_rows × n_cols] inner_type buffer; provides S
 tangent_slab views (each one [n_rows × n_cols] block) bound via
 prime_external() onto the corresponding contiguous slice.

 Storage layout matches nordsieck_block: column-major over (stage j,
 n_rows*n_cols flat), so the unified buffer is effectively a
 (n_rows*n_cols) × S matrix with leading dimension n_rows*n_cols and
 column j = stage j. Per-slot views (slab(j)) keep the row-major
 (state, deriv) layout of a tangent_slab.

 Empty stub for non-dynamic-dual T so non-AD steppers pay no size or
 codegen cost.

 Copyright (C) 2026 Simon Beyer
 */

#ifndef CPPDE_STAGE_MATRIX_HPP
#define CPPDE_STAGE_MATRIX_HPP

#include <array>
#include <cassert>
#include <cstddef>
#include <vector>

#include <cppde/cppde_dual_slab.hpp>

namespace cppde {
namespace detail {

template<class T, unsigned S, bool = is_dynamic_dual<T>::value>
class stage_matrix;

template<class T, unsigned S>
class stage_matrix<T, S, true> {
public:
  using inner_type = typename T::value_type;

  stage_matrix() = default;
  stage_matrix(const stage_matrix&)            = delete;
  stage_matrix& operator=(const stage_matrix&) = delete;
  stage_matrix(stage_matrix&&) noexcept        = default;
  stage_matrix& operator=(stage_matrix&&) noexcept = default;

  tangent_slab<T>&       slab(unsigned j)       noexcept { assert(j < S); return slabs_[j]; }
  const tangent_slab<T>& slab(unsigned j) const noexcept { assert(j < S); return slabs_[j]; }

  void prime(std::array<std::vector<T>*, S> facades,
             unsigned n_rows, unsigned n_cols) {
    const bool reshape = (n_rows != n_rows_ || n_cols != n_cols_);
    if (reshape) {
      n_rows_ = n_rows;
      n_cols_ = n_cols;
      storage_.assign(slot_stride() * S, inner_type{});
    }
    for (unsigned j = 0; j < S; ++j) {
      if (facades[j] == nullptr) continue;
      slabs_[j].prime_external(*facades[j],
                               storage_.data() + j * slot_stride(),
                               n_rows_, n_cols_);
    }
  }

private:
  std::size_t slot_stride() const noexcept {
    return static_cast<std::size_t>(n_rows_) * n_cols_;
  }

  std::vector<inner_type>          storage_;
  std::array<tangent_slab<T>, S>   slabs_;
  unsigned                         n_rows_ = 0;
  unsigned                         n_cols_ = 0;
};

template<class T, unsigned S>
class stage_matrix<T, S, false> {
public:
  void prime(std::array<std::vector<T>*, S>, unsigned, unsigned) noexcept {}

  tangent_slab<T>&       slab(unsigned)       noexcept { return stub_; }
  const tangent_slab<T>& slab(unsigned) const noexcept { return stub_; }

private:
  mutable tangent_slab<T> stub_;
};

} // namespace detail
} // namespace cppde

#endif // CPPDE_STAGE_MATRIX_HPP
