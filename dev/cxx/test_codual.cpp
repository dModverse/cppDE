// Direct C++ coverage for the reverse-AD scalar type cppde::codual.
//
// The oracle is cppde::dual: both compute the same gradient of the same
// expression, one forward and one backward, so they must agree to rounding and
// not to a solver tolerance. Every number is printed at %.17g, so the OUTPUT is
// the assertion as well: build against two revisions and diff.
//
// Build and run:  dev/cxx/run.sh --codual
//
// Copyright (C) 2026 Simon Beyer

#include <cstdio>
#include <cmath>
#include <string>
#include <vector>

#include <cppde/cppde.hpp>

using cppde::dual;
using cppde::codual;

static int g_failures = 0;

static void check(bool ok, const std::string& what) {
  if (!ok) { std::printf("FAIL  %s\n", what.c_str()); ++g_failures; }
}

// Three independents throughout, so one driver covers every arity; a function
// that needs fewer simply ignores the rest and must report a zero for them.
static constexpr unsigned NV = 3;

// Forward gradient: seed the identity, read the tangents off the result.
template<class F>
static std::vector<double> grad_forward(F f, const double* x) {
  using D = dual<double, NV>;
  D v[NV];
  for (unsigned i = 0; i < NV; ++i) { v[i] = D(x[i]); v[i].diff(i); }
  D y = f(v[0], v[1], v[2]);
  std::vector<double> g(NV, 0.0);
  for (unsigned i = 0; i < NV; ++i) g[i] = y.depend() ? y[i] : 0.0;
  return g;
}

// Reverse gradient: record, seed the output, sweep, read the inputs.
template<class F>
static std::vector<double> grad_reverse(F f, const double* x, std::size_t* nodes) {
  using C = codual<double>;
  cppde::codual_tape<double>& tp = cppde::codual_tape_for<double>();
  tp.rewind();

  C v[NV];
  for (unsigned i = 0; i < NV; ++i) { v[i] = C(x[i]); v[i].independent(); }
  C y = f(v[0], v[1], v[2]);

  tp.prepare();
  y.seed(1.0);
  tp.reverse();

  if (nodes != nullptr) *nodes = tp.size();
  std::vector<double> g(NV, 0.0);
  for (unsigned i = 0; i < NV; ++i) g[i] = v[i].adjoint();
  return g;
}

template<class F>
static void compare(const char* name, F f, double x0, double x1, double x2) {
  const double x[NV] = {x0, x1, x2};
  std::size_t nodes = 0;
  const std::vector<double> gf = grad_forward(f, x);
  const std::vector<double> gr = grad_reverse(f, x, &nodes);

  std::printf("%-14s", name);
  for (unsigned i = 0; i < NV; ++i) std::printf(" %.17g", gr[i]);
  std::printf("   nodes %zu\n", nodes);

  for (unsigned i = 0; i < NV; ++i) {
    const double scale = std::fabs(gf[i]) > 1.0 ? std::fabs(gf[i]) : 1.0;
    check(std::fabs(gf[i] - gr[i]) <= 1e-14 * scale,
          std::string(name) + " d" + std::to_string(i) +
          "  forward " + std::to_string(gf[i]) +
          "  reverse " + std::to_string(gr[i]));
  }
}

#define CMP(NAME, EXPR, X0, X1, X2) \
  compare(NAME, [](auto a, auto b, auto c) { (void)a; (void)b; (void)c; return EXPR; }, X0, X1, X2)

int main() {
  std::printf("%-14s %-20s %-20s %-20s\n", "expr", "d/da", "d/db", "d/dc");

  // Arithmetic, including the mixed codual/scalar overloads.
  CMP("a+b",        a + b,                         1.3, -0.7,  2.1);
  CMP("a-b",        a - b,                         1.3, -0.7,  2.1);
  CMP("a*b",        a * b,                         1.3, -0.7,  2.1);
  CMP("a/b",        a / b,                         1.3, -0.7,  2.1);
  CMP("-a",         -a,                            1.3, -0.7,  2.1);
  CMP("a*b*c",      a * b * c,                     1.3, -0.7,  2.1);
  CMP("2.5*a",      2.5 * a,                       1.3, -0.7,  2.1);
  CMP("a/2.5",      a / 2.5,                       1.3, -0.7,  2.1);
  CMP("2.5/a",      2.5 / a,                       1.3, -0.7,  2.1);
  CMP("3.0-a",      3.0 - a,                       1.3, -0.7,  2.1);

  // Unary math. Arguments are chosen inside each function's domain.
  CMP("exp",        cppde::exp(a),                 0.4,  0.3,  0.6);
  CMP("log",        cppde::log(a),                 0.4,  0.3,  0.6);
  CMP("sqrt",       cppde::sqrt(a),                0.4,  0.3,  0.6);
  CMP("sin",        cppde::sin(a),                 0.4,  0.3,  0.6);
  CMP("cos",        cppde::cos(a),                 0.4,  0.3,  0.6);
  CMP("tan",        cppde::tan(a),                 0.4,  0.3,  0.6);
  CMP("asin",       cppde::asin(a),                0.4,  0.3,  0.6);
  CMP("acos",       cppde::acos(a),                0.4,  0.3,  0.6);
  CMP("atan",       cppde::atan(a),                0.4,  0.3,  0.6);
  CMP("sinh",       cppde::sinh(a),                0.4,  0.3,  0.6);
  CMP("cosh",       cppde::cosh(a),                0.4,  0.3,  0.6);
  CMP("tanh",       cppde::tanh(a),                0.4,  0.3,  0.6);
  CMP("asinh",      cppde::asinh(a),               0.4,  0.3,  0.6);
  CMP("acosh",      cppde::acosh(a),               1.7,  0.3,  0.6);
  CMP("atanh",      cppde::atanh(a),               0.4,  0.3,  0.6);
  CMP("abs+",       cppde::abs(a),                 0.4,  0.3,  0.6);
  CMP("abs-",       cppde::abs(a),                -0.4,  0.3,  0.6);

  // pow in all three shapes.
  CMP("pow(a,b)",   cppde::pow(a, b),              1.4,  2.3,  0.6);
  CMP("pow(a,3.0)", cppde::pow(a, 3.0),            1.4,  2.3,  0.6);
  CMP("pow(2.0,a)", cppde::pow(2.0, a),            1.4,  2.3,  0.6);

  // Selection returns one operand whole, so the other gets a zero.
  CMP("min",        cppde::min(a, b),              1.4,  2.3,  0.6);
  CMP("max",        cppde::max(a, b),              1.4,  2.3,  0.6);
  CMP("clamp",      cppde::clamp(a, 0.0, 1.0),     1.4,  2.3,  0.6);

  // Compositions, where a shared subexpression is recorded once and its
  // adjoint accumulates from both consumers.
  CMP("shared",     cppde::exp(a * b) * (a * b),   0.7,  1.1,  0.6);
  CMP("chain",      cppde::log(cppde::sqrt(a * a + b * b) + c), 0.7, 1.1, 0.6);
  CMP("mm-kinetics", (a * b) / (c + b),            2.0,  0.5,  1.5);

  // A value built only from constants records nothing and reports no
  // derivative, which is what keeps forcings and fixed parameters off the tape.
  {
    cppde::codual_tape<double>& tp = cppde::codual_tape_for<double>();
    tp.rewind();
    codual<double> k(3.0);
    codual<double> y = cppde::exp(k * 2.0) + 1.0;
    check(tp.size() == 0, "constant expression records no node");
    check(!y.depend(), "constant expression carries no dependence");
    check(y.adjoint() == 0.0, "constant expression has a zero adjoint");
    std::printf("%-14s %.17g   nodes %zu\n", "const", y.x(), tp.size());
  }

  // Repeated seeds accumulate, which is how observation times feed one sweep.
  {
    cppde::codual_tape<double>& tp = cppde::codual_tape_for<double>();
    tp.rewind();
    codual<double> a(2.0); a.independent();
    codual<double> y = a * a;          // dy/da = 2a = 4
    tp.prepare();
    y.seed(1.0);
    y.seed(2.0);                       // total weight 3
    tp.reverse();
    check(std::fabs(a.adjoint() - 12.0) <= 1e-14, "repeated seeds accumulate");
    std::printf("%-14s %.17g\n", "seed x3", a.adjoint());
  }

  // A shift by a constant is the identity in the derivative, so it records
  // nothing and the result names its operand's node. What matters is that the
  // gradient is unchanged; the node count is the saving.
  {
    cppde::codual_tape<double>& tp = cppde::codual_tape_for<double>();
    tp.rewind();
    codual<double> a(2.0); a.independent();
    const std::size_t before = tp.size();
    codual<double> y = ((a + 3.0) - 1.0) + 0.5;    // value 4.5, derivative 1
    check(tp.size() == before, "a shift by a constant records no node");
    check(std::fabs(y.x() - 4.5) <= 1e-15, "the shift still moves the value");
    tp.prepare();
    y.seed(1.0);
    tp.reverse();
    check(std::fabs(a.adjoint() - 1.0) <= 1e-15, "the shift leaves the derivative at one");

    // The turned shift is not one: s - a has derivative -1 and records.
    tp.rewind();
    codual<double> b(2.0); b.independent();
    const std::size_t one_input = tp.size();
    codual<double> z = 5.0 - b;
    check(tp.size() == one_input + 1, "the turned shift records one node");
    tp.prepare();
    z.seed(1.0);
    tp.reverse();
    check(std::fabs(b.adjoint() + 1.0) <= 1e-15, "and turns the derivative");
    std::printf("%-14s %.17g\n", "shift", y.x());
  }

  // scope restores the node count, so a nested recording leaves no residue.
  {
    cppde::codual_tape<double>& tp = cppde::codual_tape_for<double>();
    tp.rewind();
    codual<double> a(2.0); a.independent();
    const std::size_t before = tp.size();
    {
      cppde::codual_tape<double>::scope guard(tp);
      codual<double> t = cppde::exp(a) * a;
      (void) t;
      check(tp.size() > before, "nested recording grows the tape");
    }
    check(tp.size() == before, "scope restores the node count");
    std::printf("%-14s %zu\n", "scope", tp.size());
  }

  std::printf(g_failures == 0 ? "\nOK\n" : "\n%d FAILURES\n", g_failures);
  return g_failures == 0 ? 0 : 1;
}
