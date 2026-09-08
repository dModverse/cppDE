// The sparse solve and its transpose against a dense reference.
//
// The reverse mode uses one factorisation in both directions: forward to
// recover a value a stage solved for, transposed to carry a cotangent back
// through it. The forward half has been exercised by every stiff solve the
// package has ever run; the transposed half had not been, and the first stiff
// model to lean on it disagreed with the dense path by four orders of
// magnitude while the dense path was exact.
//
// The matrix here is deliberately asymmetric, so W and W' are different
// questions and an implementation that answers one for the other is caught.
//
// Build and run:  dev/cxx/run.sh --sparse-transpose
//
// Copyright (C) 2026 Simon Beyer

#include <cstdio>
#include <cmath>
#include <string>
#include <vector>

#include <cppde/cppde.hpp>

static int g_failures = 0;

static void close(double a, double b, const std::string& what, double tol) {
  const double scale = std::fabs(a) > 1.0 ? std::fabs(a) : 1.0;
  if (std::fabs(a - b) > tol * scale) {
    std::printf("FAIL  %s  dense %.17g  sparse %.17g  rel %.2e\n",
                what.c_str(), a, b, std::fabs(a - b) / scale);
    ++g_failures;
  }
}

int main() {
#if !defined(KLU)
  std::printf("skipped: built without KLU\n");
  return 0;
#else
  // A 5x5 asymmetric matrix with a stiff diagonal, the shape an iteration
  // matrix I/gamma - J has on a stiff model.
  constexpr int n = 5;
  double dense[n * n] = {0};   // column-major, as LAPACK wants
  auto at = [&](int i, int j) -> double& { return dense[i + n * j]; };
  at(0,0) =  1.0e4; at(0,1) = -2.0;   at(0,3) =  0.5;
  at(1,0) = -3.0;   at(1,1) =  7.0;   at(1,2) =  1.5;
  at(2,1) =  0.25;  at(2,2) =  1.0e3; at(2,4) = -4.0;
  at(3,0) =  2.0;   at(3,3) =  9.0;   at(3,4) =  0.75;
  at(4,2) = -1.25;  at(4,4) =  2.0e2;

  // The same matrix in compressed sparse column.
  cppde::csc_matrix<double> W;
  W.n = n;
  W.Ap.assign(n + 1, 0);
  for (int j = 0; j < n; ++j) {
    W.Ap[j] = static_cast<int>(W.Ai.size());
    for (int i = 0; i < n; ++i)
      if (at(i, j) != 0.0) { W.Ai.push_back(i); W.Ax.push_back(at(i, j)); }
  }
  W.Ap[n] = static_cast<int>(W.Ai.size());
  W.nnz = static_cast<int>(W.Ai.size());

  const double rhs[n] = {1.0, -2.0, 3.0, 0.5, -1.5};

  // Dense reference, both directions, via the same LAPACK path the dense
  // stepper uses.
  std::vector<double> ref_fwd(rhs, rhs + n), ref_tr(rhs, rhs + n);
  {
    std::vector<double> A(dense, dense + n * n);
    std::vector<int> ipiv(n);
    int N = n, one = 1, info = 0;
    F77_CALL(dgetrf)(&N, &N, A.data(), &N, ipiv.data(), &info);
    if (info != 0) { std::printf("FAIL  dense factorisation\n"); return 1; }
    char nt = 'N', tr = 'T';
    F77_CALL(dgetrs)(&nt, &N, &one, A.data(), &N, ipiv.data(),
                     ref_fwd.data(), &N, &info FCONE);
    F77_CALL(dgetrs)(&tr, &N, &one, A.data(), &N, ipiv.data(),
                     ref_tr.data(), &N, &info FCONE);
  }

  cppde::ad_lu::sparse_lu_solver<double> lu;
  lu.analyze_pattern(W);
  lu.factorize(W);

  std::vector<double> got_fwd(rhs, rhs + n), got_tr(rhs, rhs + n);
  lu.solve(got_fwd);
  lu.solve_transposed(got_tr);

  for (int i = 0; i < n; ++i) {
    close(ref_fwd[i], got_fwd[i], "solve " + std::to_string(i), 1e-12);
    close(ref_tr[i],  got_tr[i],  "solve_transposed " + std::to_string(i), 1e-12);
  }

  std::printf("W  x = b   ");
  for (int i = 0; i < n; ++i) std::printf(" %.17g", got_fwd[i]);
  std::printf("\nW' x = b   ");
  for (int i = 0; i < n; ++i) std::printf(" %.17g", got_tr[i]);
  std::printf("\n");

  std::printf(g_failures == 0 ? "\nOK\n" : "\n%d FAILURES\n", g_failures);
  return g_failures == 0 ? 0 : 1;
#endif
}
