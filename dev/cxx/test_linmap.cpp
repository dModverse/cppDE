/*
 Constant linear maps (cppde_linmap.hpp) against written-out sums.

 Every helper runs over double, dual<double, 3> and dual2nd<double, 3>, for a
 compressed and a dense map with the same coefficients. The reference is the
 plain loop over the full coefficient matrix in the same element type, so the
 results must agree exactly. Values, first and some second derivatives are
 printed at %.17g.

 Build and run:  dev/cxx/run.sh --linmap
 */

#include <cstdio>
#include <string>
#include <vector>

#include <cppde/cppde.hpp>
#include <cppde/cppde_linmap.hpp>

using cppde::dual;
using cppde::dual2nd;

static int g_failures = 0;

static void check(bool ok, const std::string& what) {
  if (!ok) { std::printf("FAIL  %s\n", what.c_str()); ++g_failures; }
}

// 3 x 4 with a zero row and zeros inside rows.
static const int R = 3, C = 4;
static const double FULL[R * C] = {
  0.5, 0.0, -1.25, 2.0,
  0.0, 0.0, 0.0, 0.0,
  1e-3, 3.0, 0.0, -7.5,
};
static const int PTR[] = {0, 3, 3, 6};
static const int IDX[] = {0, 2, 3, 0, 1, 3};
static const double VAL[] = {0.5, -1.25, 2.0, 1e-3, 3.0, -7.5};
static const cppde::linmap_csr CSR{R, C, PTR, IDX, VAL};
static const cppde::linmap_dense DENSE{R, C, FULL};

static double val(double x) { return x; }
static double tan1(double, int) { return 0.0; }
static double val(const dual<double, 3>& x) { return x.x(); }
static double tan1(const dual<double, 3>& x, int k) { return x.d(k); }
static double val(const dual2nd<double, 3>& x) { return x.x().x(); }
static double tan1(const dual2nd<double, 3>& x, int k) { return x.d1_at(k); }
static double tan2(double, int) { return 0.0; }
static double tan2(const dual<double, 3>&, int) { return 0.0; }
static double tan2(const dual2nd<double, 3>& x, int k) { return x.dd_at(k, (k + 1) % 3); }

template<class T>
static bool same(const T& a, const T& b) {
  if (val(a) != val(b)) return false;
  for (int k = 0; k < 3; ++k)
    if (tan1(a, k) != tan1(b, k) || tan2(a, k) != tan2(b, k)) return false;
  return true;
}

template<class T>
static void emit(const char* tag, const std::vector<T>& y) {
  std::printf("%-24s", tag);
  for (const T& v : y) {
    std::printf(" %.17g", val(v));
    for (int k = 0; k < 3; ++k) std::printf(" %.17g", tan1(v, k));
    for (int k = 0; k < 3; ++k) std::printf(" %.17g", tan2(v, k));
  }
  std::printf("\n");
}

template<class T> static void seed(T&, int) {}
static void seed(dual<double, 3>& x, int k) { x.diff(k % 3); }
static void seed(dual2nd<double, 3>& x, int k) { x.x().diff(k % 3); x.diff(k % 3); }

template<class T>
static std::vector<T> inputs(int n, double base) {
  std::vector<T> x(n);
  for (int j = 0; j < n; ++j) {
    x[j] = T(base + 0.37 * j - 0.11 * j * j);
    seed(x[j], j);
  }
  return x;
}

template<class T, class Map>
static void run(const char* type, const char* kind, const Map& m) {
  const std::string pre = std::string(type) + "/" + kind + " ";
  std::vector<T> x = inputs<T>(C, 0.8);
  std::vector<T> w = inputs<T>(R, -0.3);

  // y = C x
  std::vector<T> y(R), yr(R);
  cppde::apply(m, x, y);
  for (int r = 0; r < R; ++r) {
    T s(0.0);
    for (int j = 0; j < C; ++j) s += FULL[r * C + j] * x[j];
    yr[r] = s;
    check(same(y[r], yr[r]), pre + "apply row " + std::to_string(r));
    check(same(cppde::row_dot(m, r, x), yr[r]), pre + "row_dot " + std::to_string(r));
  }
  emit((pre + "apply").c_str(), y);

  // z += C^T w
  std::vector<T> z = inputs<T>(C, 2.0), zr = z;
  cppde::apply_t_add(m, w, z);
  for (int j = 0; j < C; ++j) {
    for (int r = 0; r < R; ++r) zr[j] += w[r] * FULL[r * C + j];
    check(same(z[j], zr[j]), pre + "apply_t_add " + std::to_string(j));
  }
  emit((pre + "apply_t_add").c_str(), z);

  // J(1, .) += a * C[2, .], a an expression
  cppde::dense_matrix<T> J(R, C);
  cppde::dense_matrix<T> Jr(R, C);
  cppde::axpy_row_dense(m, 2, -(w[0] * w[1]), J, 1);
  T a = -(w[0] * w[1]);
  for (int j = 0; j < C; ++j) {
    Jr(1, j) += a * FULL[2 * C + j];
    for (int i = 0; i < R; ++i)
      check(same(J(i, j), Jr(i, j)), pre + "axpy_row_dense J " + std::to_string(i));
  }
  std::vector<T> jrow(C);
  for (int j = 0; j < C; ++j) jrow[j] = J(1, j);
  emit((pre + "axpy_row_dense J").c_str(), jrow);

  // ax[dst] += a * row 0, dst over the stored coefficients of the row
  const bool dense = kind[0] == 'd';
  const int stored = dense ? C : PTR[1] - PTR[0];
  std::vector<int> dst(stored);
  for (int k = 0; k < stored; ++k) dst[k] = 2 * k + 1;
  std::vector<T> ax(2 * C + 1, T(0.0)), axr = ax;
  cppde::axpy_row_idx(m, 0, a, ax, dst.data());
  for (int k = 0; k < stored; ++k) {
    const int j = dense ? k : IDX[k];
    axr[dst[k]] += a * FULL[j];
  }
  for (size_t k = 0; k < ax.size(); ++k)
    check(same(ax[k], axr[k]), pre + "axpy_row_idx " + std::to_string(k));
  emit((pre + "axpy_row_idx").c_str(), ax);

  // a zero row leaves y alone
  std::vector<T> y1 = x;
  cppde::axpy_row_dense(m, 1, a, y1);
  for (int j = 0; j < C; ++j) check(same(y1[j], x[j]), pre + "zero row " + std::to_string(j));
}

int main() {
  run<double>("double", "csr", CSR);
  run<double>("double", "dense", DENSE);
  run<dual<double, 3>>("dual", "csr", CSR);
  run<dual<double, 3>>("dual", "dense", DENSE);
  run<dual2nd<double, 3>>("dual2nd", "csr", CSR);
  run<dual2nd<double, 3>>("dual2nd", "dense", DENSE);

  // a double map applied to a pointer target of an AD type
  std::vector<dual<double, 3>> w = inputs<dual<double, 3>>(R, 1.5);
  std::vector<dual<double, 3>> out(C + 2, dual<double, 3>(0.0));
  cppde::apply_t_add(CSR, w, out);
  dual<double, 3>* p = out.data();
  cppde::axpy_row_dense(CSR, 0, w[2], p);
  emit("dual pointer target", out);

  // coefficients through a pool
  static const double POOL[] = {-7.5, 0.5, 3.0, -1.25, 1e-3, 2.0};
  static const unsigned short VID[] = {1, 3, 5, 4, 2, 0};
  std::vector<double> pooled = cppde::linmap_values(POOL, VID, 6);
  for (int k = 0; k < 6; ++k) check(pooled[k] == VAL[k], "linmap_values " + std::to_string(k));

  if (g_failures) {
    std::printf("%d failure(s)\n", g_failures);
    return 1;
  }
  std::printf("all linmap checks passed\n");
  return 0;
}
