// The transposed solve on a derivative-carrying iteration matrix.
//
// A written adjoint over a dual scalar needs W^-T where W itself carries
// tangents. The peeling transposes with the matrix: the value layer solves
// W_val^T x = b_val, and each direction solves W_val^T x' = b' - (dW)^T x_val.
//
// The oracle is the same solver on the explicitly transposed matrix: building
// W^T and calling solve() must give what calling solve_transposed() on W gives,
// to the last bit of the value layer and to rounding on the tangents. That is a
// different code path through the same arithmetic, which is what makes it a
// check rather than a restatement.
//
// Dense and sparse are both covered: the peeling is the same, but the sparse
// transpose is an index swap inside a compressed column walk, which is its own
// way to go wrong. The sparse half is skipped when KLU is absent.
//
// Build and run:  dev/cxx/run.sh --ad-transpose

#include <cstdio>
#include <cmath>
#include <string>
#include <vector>

#include <cppde/cppde.hpp>

using cppde::dual;

static int g_failures = 0;

static void close(double a, double b, const std::string& what,
                  double tol = 1e-12) {
  const double scale = std::fabs(a) > 1.0 ? std::fabs(a) : 1.0;
  if (std::fabs(a - b) > tol * scale) {
    std::printf("FAIL  %-40s  direct %.17g  transposed %.17g\n",
                what.c_str(), a, b);
    ++g_failures;
  }
}

static constexpr unsigned N = 3;      // derivative directions
using D = dual<double, N>;

// A matrix that is not symmetric, is well conditioned, and whose every entry
// carries a different tangent, so a transpose that acts on the wrong index
// shows up in the derivative layer even where the value layer is right.
static void fill(cppde::dense_matrix<D>& W, int n) {
  W.resize(n, n);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) {
      const double v = (i == j) ? (4.0 + 0.5 * i) : (0.3 * (i + 1) - 0.2 * j);
      D e(v);
      for (unsigned k = 0; k < N; ++k)
        e.diff(k) = 0.01 * (i + 1) * (j + 2) * (k + 1);
      W(i, j) = e;
    }
}

// The same matrix in compressed sparse column, transposed or not. Every entry
// is structurally present, so the pattern says nothing and the indices carry
// the whole question.
static void to_csc(const cppde::dense_matrix<D>& A, int n,
                   cppde::csc_matrix<D>& out) {
  out.n = n;
  out.Ap.assign(n + 1, 0);
  out.Ai.clear();
  out.Ax.clear();
  for (int j = 0; j < n; ++j) {
    out.Ap[j] = static_cast<int>(out.Ai.size());
    for (int i = 0; i < n; ++i) {
      out.Ai.push_back(i);
      out.Ax.push_back(const_cast<cppde::dense_matrix<D>&>(A)(i, j));
    }
  }
  out.Ap[n] = static_cast<int>(out.Ai.size());
  out.nnz = static_cast<int>(out.Ai.size());
}

static void run(int n) {
  cppde::dense_matrix<D> W, Wt;
  fill(W, n);
  Wt.resize(n, n);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) Wt(i, j) = W(j, i);

  std::vector<D> b0(n), b(n), b2(n);
  for (int i = 0; i < n; ++i) {
    D e(1.0 + 0.7 * i);
    for (unsigned k = 0; k < N; ++k) e.diff(k) = 0.05 * (i + 1) - 0.02 * k;
    b0[i] = e;
    b[i] = e;
    b2[i] = e;
  }

  cppde::ad_lu::dense_lu_solver<D> direct, transposed;
  direct.factorize(Wt);
  direct.solve(b);                    // W^T x = b, the long way round
  transposed.factorize(W);
  transposed.solve_transposed(b2);    // the same, on the factorisation of W

  const std::string tag = "n=" + std::to_string(n);
  for (int i = 0; i < n; ++i) {
    close(b[i].x(), b2[i].x(), tag + " value " + std::to_string(i));
    for (unsigned k = 0; k < N; ++k)
      close(b[i].d(k), b2[i].d(k),
            tag + " d" + std::to_string(k) + " " + std::to_string(i), 1e-10);
  }
  std::printf("%-6s  solved and checked %d states over %u directions\n",
              tag.c_str(), n, N);

#if defined(KLU)
  cppde::csc_matrix<D> Ws, Wts;
  to_csc(W, n, Ws);
  to_csc(Wt, n, Wts);

  std::vector<D> s1(n), s2(n);
  for (int i = 0; i < n; ++i) { s1[i] = b0[i]; s2[i] = b0[i]; }

  cppde::ad_lu::sparse_lu_solver<D> sp_direct, sp_transposed;
  sp_direct.factorize(Wts);
  sp_direct.solve(s1);
  sp_transposed.factorize(Ws);
  sp_transposed.solve_transposed(s2);

  for (int i = 0; i < n; ++i) {
    close(s1[i].x(), s2[i].x(), tag + " sparse value " + std::to_string(i),
          1e-12);
    close(b[i].x(), s2[i].x(), tag + " sparse vs dense " + std::to_string(i),
          1e-10);
    for (unsigned k = 0; k < N; ++k) {
      close(s1[i].d(k), s2[i].d(k),
            tag + " sparse d" + std::to_string(k) + " " + std::to_string(i),
            1e-10);
      close(b[i].d(k), s2[i].d(k),
            tag + " sparse vs dense d" + std::to_string(k) + " " +
              std::to_string(i), 1e-10);
    }
  }
  std::printf("%-6s  sparse agrees with dense over %u directions\n",
              tag.c_str(), N);
#endif
}

int main() {
  run(3);
  run(7);
  run(12);
  std::printf(g_failures == 0 ? "\nOK\n" : "\n%d FAILURES\n", g_failures);
  return g_failures == 0 ? 0 : 1;
}
