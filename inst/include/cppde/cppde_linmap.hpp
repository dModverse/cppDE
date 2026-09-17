/*
 Constant linear maps of generated models.

 A map C holds the numeric coefficients of long linear sums over the state
 vector as static tables. Generated code evaluates C x once per call and
 assembles Jacobian rows as J += G C. The element type of the vectors is a
 template parameter, so the same code serves double, dual and dual2nd.

 Copyright (C) 2026 Simon Beyer
 */

#ifndef CPPDE_LINMAP_HPP
#define CPPDE_LINMAP_HPP

#include <type_traits>
#include <utility>
#include <vector>

namespace cppde {

// Compressed rows over static arrays.
struct linmap_csr {
  int rows;
  int cols;
  const int* ptr;     // rows + 1
  const int* idx;     // ptr[rows]
  const double* val;  // ptr[rows]
};

// Row-major dense map over a static array.
struct linmap_dense {
  int rows;
  int cols;
  const double* val;  // rows * cols
};

namespace linmap_detail {

template<class V>
using elem_t = std::decay_t<decltype(std::declval<V&>()[0])>;

template<class M>
using mat_elem_t = std::decay_t<decltype(std::declval<M&>()(0, 0))>;

} // namespace linmap_detail

// sum_j C[r, j] * x[j]
template<class X>
linmap_detail::elem_t<const X> row_dot(const linmap_csr& m, int r, const X& x) {
  linmap_detail::elem_t<const X> s(0.0);
  for (int k = m.ptr[r]; k < m.ptr[r + 1]; ++k) s += m.val[k] * x[m.idx[k]];
  return s;
}

template<class X>
linmap_detail::elem_t<const X> row_dot(const linmap_dense& m, int r, const X& x) {
  linmap_detail::elem_t<const X> s(0.0);
  const double* v = m.val + static_cast<long>(r) * m.cols;
  for (int j = 0; j < m.cols; ++j) s += v[j] * x[j];
  return s;
}

// y = C x
template<class Map, class X, class Y>
void apply(const Map& m, const X& x, Y& y) {
  for (int r = 0; r < m.rows; ++r) y[r] = row_dot(m, r, x);
}

// y[j] += a * C[r, j]
template<class A, class Y>
void axpy_row_dense(const linmap_csr& m, int r, const A& a, Y& y) {
  const linmap_detail::elem_t<Y> s(a);
  for (int k = m.ptr[r]; k < m.ptr[r + 1]; ++k) y[m.idx[k]] += s * m.val[k];
}

template<class A, class Y>
void axpy_row_dense(const linmap_dense& m, int r, const A& a, Y& y) {
  const linmap_detail::elem_t<Y> s(a);
  const double* v = m.val + static_cast<long>(r) * m.cols;
  for (int j = 0; j < m.cols; ++j) y[j] += s * v[j];
}

// J(i, j) += a * C[r, j]
template<class A, class Mat>
void axpy_row_dense(const linmap_csr& m, int r, const A& a, Mat& J, int i) {
  const linmap_detail::mat_elem_t<Mat> s(a);
  for (int k = m.ptr[r]; k < m.ptr[r + 1]; ++k) J(i, m.idx[k]) += s * m.val[k];
}

template<class A, class Mat>
void axpy_row_dense(const linmap_dense& m, int r, const A& a, Mat& J, int i) {
  const linmap_detail::mat_elem_t<Mat> s(a);
  const double* v = m.val + static_cast<long>(r) * m.cols;
  for (int j = 0; j < m.cols; ++j) J(i, j) += s * v[j];
}

// ax[dst[k]] += a * (k-th stored coefficient of row r)
template<class A, class V>
void axpy_row_idx(const linmap_csr& m, int r, const A& a, V& ax, const int* dst) {
  const linmap_detail::elem_t<V> s(a);
  const int k0 = m.ptr[r];
  for (int k = k0; k < m.ptr[r + 1]; ++k) ax[dst[k - k0]] += s * m.val[k];
}

template<class A, class V>
void axpy_row_idx(const linmap_dense& m, int r, const A& a, V& ax, const int* dst) {
  const linmap_detail::elem_t<V> s(a);
  const double* v = m.val + static_cast<long>(r) * m.cols;
  for (int j = 0; j < m.cols; ++j) ax[dst[j]] += s * v[j];
}

// val[k] = pool[id[k]]
template<class Id>
std::vector<double> linmap_values(const double* pool, const Id* id, int nnz) {
  std::vector<double> val(nnz);
  for (int k = 0; k < nnz; ++k) val[k] = pool[id[k]];
  return val;
}

// y += C^T w
template<class Map, class W, class Y>
void apply_t_add(const Map& m, const W& w, Y& y) {
  for (int r = 0; r < m.rows; ++r) axpy_row_dense(m, r, w[r], y);
}

} // namespace cppde

#endif // CPPDE_LINMAP_HPP
