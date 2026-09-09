/*
 Where a reverse step spends its time.

 A seeded solve does three things a value solve does not: it records the
 arithmetic on a tape, it walks that tape backwards, and it pays for whatever
 the recording type costs over a double. This separates the three so an
 optimisation can be aimed rather than guessed at.

 Reported as multiples of the same right-hand side evaluated in plain double,
 which is the only unit that survives a change of machine.

   build: dev/cxx/run.sh --bench-codual
 */

#include <cppde/cppde.hpp>
#include <cppde/cppde_codual.hpp>
#include <cppde/cppde_codual_math.hpp>

#include <chrono>
#include <cstdio>
#include <vector>

using cppde::codual;
using cppde::dual;

// A right-hand side with the shape a real model has: a few products, a
// division, and a transcendental, over a handful of states.
template<class T>
static void rhs(const T* x, const T* p, T* out) {
  const T s = x[0] + x[1] + x[2];
  out[0] = -p[0] * x[0] + p[1] * x[1] * x[2] / (T(1) + s);
  out[1] =  p[0] * x[0] - p[1] * x[1] * x[2] - p[2] * x[1] * x[1];
  out[2] =  p[2] * x[1] * x[1] * cppde::exp(-p[3] * x[2]);
}

static const double X0[3] = {1.0, 0.4, 0.1};
static const double P0[4] = {0.7, 0.35, 1.1, 0.25};

template<class F>
static double timeit(F&& f, int reps = 5) {
  double best = 1e300;
  for (int r = 0; r < reps; ++r) {
    const auto t0 = std::chrono::steady_clock::now();
    f();
    const auto t1 = std::chrono::steady_clock::now();
    const double dt = std::chrono::duration<double>(t1 - t0).count();
    if (dt < best) best = dt;
  }
  return best;
}

int main() {
  const int N = 200000;   // right-hand side evaluations per measurement
  double sink = 0.0;

  // --- plain double, the unit ------------------------------------------------
  const double t_dbl = timeit([&] {
    double x[3] = {X0[0], X0[1], X0[2]}, o[3];
    for (int i = 0; i < N; ++i) {
      x[0] = X0[0] + 1e-12 * i;
      rhs<double>(x, P0, o);
      sink += o[0] + o[1] + o[2];
    }
  });

  // --- forward dual over 4 parameters, for scale -----------------------------
  using D4 = dual<double, 4>;
  const double t_dual = timeit([&] {
    D4 x[3], p[4], o[3];
    for (unsigned j = 0; j < 4; ++j) { p[j] = D4(P0[j]); p[j].diff(j); }
    for (int i = 0; i < N; ++i) {
      for (int k = 0; k < 3; ++k) x[k] = D4(X0[k]);
      x[0] = D4(X0[0] + 1e-12 * i);
      rhs<D4>(x, p, o);
      sink += o[0].x() + o[1].x() + o[2].x();
    }
  });

  // --- codual: recording only ------------------------------------------------
  using C = codual<double>;
  cppde::codual_tape<double>& tp = cppde::codual_tape_for<double>();
  std::size_t nodes_per_eval = 0;

  const double t_rec = timeit([&] {
    C x[3], p[4], o[3];
    for (int i = 0; i < N; ++i) {
      tp.rewind();
      for (int k = 0; k < 3; ++k) { x[k] = C(X0[k]); x[k].independent(); }
      for (int j = 0; j < 4; ++j) { p[j] = C(P0[j]); p[j].independent(); }
      x[0] = C(X0[0] + 1e-12 * i); x[0].independent();
      rhs<C>(x, p, o);
      nodes_per_eval = tp.size();
      sink += o[0].x() + o[1].x() + o[2].x();
    }
  });

  // --- codual: recording plus the backward sweep -----------------------------
  const double t_sweep = timeit([&] {
    C x[3], p[4], o[3];
    for (int i = 0; i < N; ++i) {
      tp.rewind();
      for (int k = 0; k < 3; ++k) { x[k] = C(X0[k]); x[k].independent(); }
      for (int j = 0; j < 4; ++j) { p[j] = C(P0[j]); p[j].independent(); }
      x[0] = C(X0[0] + 1e-12 * i); x[0].independent();
      rhs<C>(x, p, o);
      tp.prepare();
      for (int k = 0; k < 3; ++k) o[k].seed(1.0);
      tp.reverse();
      sink += x[0].adjoint() + p[0].adjoint();
    }
  });

  // --- the tape alone, no codual layer and no thread-local lookup ----------
  // Same node count, written through a reference the compiler keeps in a
  // register. What it leaves out is exactly what the operator layer and the
  // per-operation lookup cost.
  const std::size_t NPE = nodes_per_eval;
  const double t_raw = timeit([&] {
    for (int i = 0; i < N; ++i) {
      tp.rewind();
      const std::size_t s0 = tp.independent();
      std::size_t s1 = tp.independent();
      for (std::size_t k = 2; k < NPE; ++k)
        s1 = tp.record(s0, 1.0, s1, 1.0);
      sink += (double)s1;
    }
  });

  // --- one thread-local lookup per node, nothing else ----------------------
  const double t_tls = timeit([&] {
    for (int i = 0; i < N; ++i)
      for (std::size_t k = 0; k < NPE; ++k)
        sink += (double)cppde::codual_tape_for<double>().size();
  });

  std::printf("right-hand side evaluations: %d, tape nodes each: %zu\n\n",
              N, nodes_per_eval);
  std::printf("%-28s %-12s %-8s\n", "", "seconds", "x double");
  std::printf("%-28s %-12.5f %-8.2f\n", "plain double",       t_dbl,   1.0);
  std::printf("%-28s %-12.5f %-8.2f\n", "forward dual, 4 par", t_dual,  t_dual  / t_dbl);
  std::printf("%-28s %-12.5f %-8.2f\n", "codual, record only", t_rec,   t_rec   / t_dbl);
  std::printf("%-28s %-12.5f %-8.2f\n", "codual, record+sweep", t_sweep, t_sweep / t_dbl);
  std::printf("%-28s %-12.5f %-8.2f\n", "tape writes alone",   t_raw, t_raw / t_dbl);
  std::printf("%-28s %-12.5f %-8.2f\n", "tls lookups alone",   t_tls, t_tls / t_dbl);

  const double per = 1e9 / (double)N / (double)nodes_per_eval;
  std::printf("\nper tape node, nanoseconds\n");
  std::printf("  recording, all in      %6.2f\n", (t_rec - t_dbl) * per);
  std::printf("  tape write alone       %6.2f\n", t_raw * per);
  std::printf("  thread-local lookup    %6.2f\n", t_tls * per);
  std::printf("  the backward sweep     %6.2f\n", (t_sweep - t_rec) * per);

  if (sink == 12345.6789) std::printf(" ");   // keep the work
  return 0;
}
