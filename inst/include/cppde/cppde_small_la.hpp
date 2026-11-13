/*
 Dense kernels of the step path in plain loops: the LU factorization with partial
 pivoting and the solves with it, and the products of the corrector and Nordsieck updates.

 OpenBLAS takes the work buffers of the level-2 and level-3 routines from an
 allocator behind one process-wide lock. Called once per step from the threads of
 a batch, that lock serializes them; these loops take none. All matrices are
 column-major.

 Copyright (C) 2026 Simon Beyer
 */

#ifndef CPPDE_SMALL_LA_HPP
#define CPPDE_SMALL_LA_HPP

#include <cmath>
#include <cstddef>
#include <utility>

namespace cppde {
namespace small_la {

// Systems up to this many unknowns are factorized and solved here, larger ones by
// LAPACK, whose blocked kernels outweigh the lock there.
constexpr int max_n = 64;

// LU factorization P A = L U of the n x n matrix `a` in place, L unit lower, as
// dgetf2 does; `ipiv[k]` is the row swapped with row k, 0-based. Returns 0, or
// k + 1 for the first zero pivot k, whose column is then left as it is.
inline int getrf(int n, double* a, int* ipiv)
{
  int info = 0;
  for (int k = 0; k < n; ++k) {
    double* ck = a + static_cast<std::size_t>(k) * n;
    int p = k;
    double amax = std::fabs(ck[k]);
    for (int i = k + 1; i < n; ++i) {
      const double v = std::fabs(ck[i]);
      if (v > amax) { amax = v; p = i; }
    }
    ipiv[k] = p;
    if (ck[p] == 0.0) {
      if (info == 0) info = k + 1;
      continue;
    }
    if (p != k)
      for (int j = 0; j < n; ++j)
        std::swap(a[k + static_cast<std::size_t>(j) * n], a[p + static_cast<std::size_t>(j) * n]);
    const double inv = 1.0 / ck[k];
    for (int i = k + 1; i < n; ++i) ck[i] *= inv;
    for (int j = k + 1; j < n; ++j) {
      double* cj = a + static_cast<std::size_t>(j) * n;
      const double f = cj[k];
      if (f == 0.0) continue;
      for (int i = k + 1; i < n; ++i) cj[i] -= ck[i] * f;
    }
  }
  return info;
}

// Solves A X = B (`trans` 'N') or A' X = B ('T') with the factors of getrf(); B is
// n x nrhs with leading dimension ldb and is overwritten by X.
inline void getrs(char trans, int n, int nrhs, const double* a, const int* ipiv,
                  double* b, int ldb)
{
  for (int c = 0; c < nrhs; ++c) {
    double* x = b + static_cast<std::size_t>(c) * ldb;
    if (trans == 'N') {
      for (int k = 0; k < n; ++k)
        if (ipiv[k] != k) std::swap(x[k], x[ipiv[k]]);
      for (int j = 0; j < n; ++j) {
        const double xj = x[j];
        if (xj == 0.0) continue;
        const double* cj = a + static_cast<std::size_t>(j) * n;
        for (int i = j + 1; i < n; ++i) x[i] -= cj[i] * xj;
      }
      for (int j = n - 1; j >= 0; --j) {
        const double* cj = a + static_cast<std::size_t>(j) * n;
        x[j] /= cj[j];
        const double xj = x[j];
        if (xj == 0.0) continue;
        for (int i = 0; i < j; ++i) x[i] -= cj[i] * xj;
      }
    } else {
      for (int j = 0; j < n; ++j) {
        const double* cj = a + static_cast<std::size_t>(j) * n;
        double s = x[j];
        for (int i = 0; i < j; ++i) s -= cj[i] * x[i];
        x[j] = s / cj[j];
      }
      for (int j = n - 1; j >= 0; --j) {
        const double* cj = a + static_cast<std::size_t>(j) * n;
        double s = x[j];
        for (int i = j + 1; i < n; ++i) s -= cj[i] * x[i];
        x[j] = s;
      }
      for (int k = n - 1; k >= 0; --k)
        if (ipiv[k] != k) std::swap(x[k], x[ipiv[k]]);
    }
  }
}

// y -= A x for the m x n matrix A with leading dimension lda.
inline void gemv_sub(int m, int n, const double* A, int lda, const double* x, double* y)
{
  for (int j = 0; j < n; ++j) {
    const double xj = x[j];
    if (xj == 0.0) continue;
    const double* cj = A + static_cast<std::size_t>(j) * lda;
    for (int i = 0; i < m; ++i) y[i] -= cj[i] * xj;
  }
}

// y -= A' x for the m x n matrix A with leading dimension lda; y has n entries.
inline void gemv_t_sub(int m, int n, const double* A, int lda, const double* x, double* y)
{
  for (int j = 0; j < n; ++j) {
    const double* cj = A + static_cast<std::size_t>(j) * lda;
    double s = 0.0;
    for (int i = 0; i < m; ++i) s += cj[i] * x[i];
    y[j] -= s;
  }
}

// C -= A B for A m x k (lda), B k x nrhs (ldb), C m x nrhs (ldc).
inline void gemm_sub(int m, int nrhs, int k, const double* A, int lda,
                     const double* B, int ldb, double* C, int ldc)
{
  for (int c = 0; c < nrhs; ++c)
    gemv_sub(m, k, A, lda, B + static_cast<std::size_t>(c) * ldb,
             C + static_cast<std::size_t>(c) * ldc);
}

// B := B P for B m x n (ldb) and P n x n unit lower triangular (ldp). Column j of
// the product reads the columns j and above, so ascending j works in place.
inline void trmm_right_lower_unit(int m, int n, const double* P, int ldp, double* B, int ldb)
{
  for (int j = 0; j < n; ++j) {
    double* bj = B + static_cast<std::size_t>(j) * ldb;
    for (int k = j + 1; k < n; ++k) {
      const double pkj = P[k + static_cast<std::size_t>(j) * ldp];
      if (pkj == 0.0) continue;
      const double* bk = B + static_cast<std::size_t>(k) * ldb;
      for (int i = 0; i < m; ++i) bj[i] += bk[i] * pkj;
    }
  }
}

// A += x y' for A m x n (lda).
inline void ger(int m, int n, const double* x, const double* y, double* A, int lda)
{
  for (int j = 0; j < n; ++j) {
    const double yj = y[j];
    if (yj == 0.0) continue;
    double* cj = A + static_cast<std::size_t>(j) * lda;
    for (int i = 0; i < m; ++i) cj[i] += x[i] * yj;
  }
}

}  // namespace small_la
}  // namespace cppde

#endif
