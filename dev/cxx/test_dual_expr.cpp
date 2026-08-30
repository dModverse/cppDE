// Direct C++ coverage for the forward-AD expression templates.
//
// Every computed number is printed at %.17g, so the OUTPUT is the assertion:
// build against two revisions and diff. That is a stronger oracle than the R
// suite's 1e-8/1e-10 tolerances, which would swallow a small perturbation.
//
// Build and run:  dev/cxx/run.sh
//
// Copyright (C) 2026 Simon Beyer

#include <cstdio>
#include <cmath>
#include <string>
#include <type_traits>

#include <cppde/cppde.hpp>

using cppde::dual;
namespace ex = cppde::expr;

static int g_failures = 0;

static void check(bool ok, const std::string& what) {
  if (!ok) { std::printf("FAIL  %s\n", what.c_str()); ++g_failures; }
}

static void close_to(double got, double want, const std::string& what) {
  const double scale = std::fabs(want) > 1.0 ? std::fabs(want) : 1.0;
  check(std::fabs(got - want) <= 1e-12 * scale,
        what + "  got " + std::to_string(got) + " want " + std::to_string(want));
}

// Print value and every tangent, so the diff covers all of them.
template <class D>
static void emit(const char* tag, const D& d, unsigned n) {
  std::printf("%-28s %.17g", tag, d.x());
  for (unsigned j = 0; j < n; ++j) std::printf(" %.17g", d[j]);
  std::printf("\n");
}

// =============================================================================
// Structural assertions: node identity and size.
// =============================================================================

using D3 = dual<double, 3>;
using LD = ex::DualLeaf<double, 3>;
using LS = ex::ScalarLeaf<double>;
using Mul = ex::BinExpr<LD, LD, ex::MulOp>;

static_assert(std::is_same<decltype(std::declval<D3>() * std::declval<D3>()),
                           Mul>::value,
              "dual * dual must build BinExpr<DualLeaf, DualLeaf, MulOp>");
static_assert(std::is_same<decltype(std::declval<D3>() * 2.0),
                           ex::BinExpr<LD, LS, ex::MulOp>>::value,
              "dual * scalar must align the scalar leaf to the dual's T");
static_assert(std::is_same<decltype(2.0 * std::declval<D3>()),
                           ex::BinExpr<LS, LD, ex::MulOp>>::value,
              "scalar * dual");
static_assert(std::is_same<decltype(std::declval<D3>() * std::declval<D3>()
                                    + std::declval<D3>()),
                           ex::BinExpr<Mul, LD, ex::AddOp>>::value,
              "composite + dual nests the composite directly");

// Children are stored by value (see cppde_dual_expr.hpp).
static_assert(!std::is_reference<decltype(Mul::l_)>::value, "leaf by value");
static_assert(!std::is_reference<
                  decltype(ex::BinExpr<Mul, Mul, ex::AddOp>::l_)>::value,
              "composite by value");

// =============================================================================
// Numerical checks, run for both materialiser families:
// dual<T,0> (runtime slab) and dual<T,N> (compile-time slab) are separate code
// paths in cppde_dual_expr.hpp and must both be covered.
// =============================================================================

template <class D>
static void seed(D& d, unsigned idx, unsigned n) {
  if constexpr (std::is_same<D, dual<double, 0>>::value) d.diff(idx, n);
  else                                                   d.diff(idx, n);
}

template <class D>
static void arithmetic_suite(const char* label, unsigned N) {
  D a(2.5), b(1.5), c(3.25);
  seed(a, 0, N); seed(b, 1 % N, N); seed(c, 2 % N, N);

  // The eight operand combinations, each op.
  D r;
  r = a * b;        emit((std::string(label) + " a*b").c_str(), r, N);
  close_to(r.x(), 3.75, "a*b value");
  r = a * 2.0;      emit((std::string(label) + " a*s").c_str(), r, N);
  r = 2.0 * a;      emit((std::string(label) + " s*a").c_str(), r, N);
  r = a * b + c;    emit((std::string(label) + " a*b+c").c_str(), r, N);
  r = c + a * b;    emit((std::string(label) + " c+a*b").c_str(), r, N);
  r = a * b - 2.0;  emit((std::string(label) + " a*b-s").c_str(), r, N);
  r = 2.0 - a * b;  emit((std::string(label) + " s-a*b").c_str(), r, N);
  r = (a * b) / (b * c);
  emit((std::string(label) + " (a*b)/(b*c)").c_str(), r, N);

  r = -a;           emit((std::string(label) + " -a").c_str(), r, N);
  r = +a;           emit((std::string(label) + " +a").c_str(), r, N);

  // Every math function, at arguments inside its domain (UBSan would flag
  // otherwise: log/sqrt need x>0, asin/acos/atanh need |x|<1, acosh needs x>1).
  D u(0.5), v(2.0);
  seed(u, 0, N); seed(v, 1 % N, N);
  r = exp(u);   emit((std::string(label) + " exp(u)").c_str(), r, N);
  r = log(v);   emit((std::string(label) + " log(v)").c_str(), r, N);
  r = sqrt(v);  emit((std::string(label) + " sqrt(v)").c_str(), r, N);
  r = sin(u);   emit((std::string(label) + " sin(u)").c_str(), r, N);
  r = cos(u);   emit((std::string(label) + " cos(u)").c_str(), r, N);
  r = tan(u);   emit((std::string(label) + " tan(u)").c_str(), r, N);
  r = asin(u);  emit((std::string(label) + " asin(u)").c_str(), r, N);
  r = acos(u);  emit((std::string(label) + " acos(u)").c_str(), r, N);
  r = atan(u);  emit((std::string(label) + " atan(u)").c_str(), r, N);
  r = sinh(u);  emit((std::string(label) + " sinh(u)").c_str(), r, N);
  r = cosh(u);  emit((std::string(label) + " cosh(u)").c_str(), r, N);
  r = tanh(u);  emit((std::string(label) + " tanh(u)").c_str(), r, N);
  r = asinh(u); emit((std::string(label) + " asinh(u)").c_str(), r, N);
  r = acosh(v); emit((std::string(label) + " acosh(v)").c_str(), r, N);
  r = atanh(u); emit((std::string(label) + " atanh(u)").c_str(), r, N);
  r = abs(-u);  emit((std::string(label) + " abs(-u)").c_str(), r, N);

  r = pow(v, u);   emit((std::string(label) + " pow(v,u)").c_str(), r, N);
  r = pow(v, 3.0); emit((std::string(label) + " pow(v,3)").c_str(), r, N);
  r = pow(3.0, u); emit((std::string(label) + " pow(3,u)").c_str(), r, N);

  // A deep composite mixing every node kind.
  r = ((a * b + c / v) - exp(u)) * sqrt(v) + pow(a, 2.0) - (-c);
  emit((std::string(label) + " deep").c_str(), r, N);
}

// =============================================================================
// Aliasing. Every pattern below is lifted from production code and runs for
// several iterations, because the failure mode this guards against is
// "first call correct, later calls drift".
// =============================================================================

template <class D>
static void aliasing_suite(const char* label, unsigned N) {
  // cppde_dual_math.hpp: *this = *this + o  (and the other compound forms)
  {
    D a(1.25), o(0.5);
    seed(a, 0, N); seed(o, 1 % N, N);
    for (int k = 0; k < 6; ++k) {
      a = a + o;
      emit((std::string(label) + " alias a=a+o").c_str(), a, N);
    }
    close_to(a.x(), 1.25 + 6 * 0.5, "alias a=a+o value");
  }
  { D a(4.0), o(0.5); seed(a, 0, N); seed(o, 1 % N, N);
    for (int k = 0; k < 6; ++k) { a = a - o;
      emit((std::string(label) + " alias a=a-o").c_str(), a, N); } }
  { D a(1.5), o(1.1); seed(a, 0, N); seed(o, 1 % N, N);
    for (int k = 0; k < 6; ++k) { a = a * o;
      emit((std::string(label) + " alias a=a*o").c_str(), a, N); } }
  { D a(1.0); seed(a, 0, N);
    for (int k = 0; k < 6; ++k) { a += D(0.25);
      emit((std::string(label) + " alias a+=c").c_str(), a, N); } }
  { D a(8.0); seed(a, 0, N);
    for (int k = 0; k < 6; ++k) { a *= D(0.5);
      emit((std::string(label) + " alias a*=c").c_str(), a, N); } }

  // cppde_multistepper.hpp:1947   x[i] = x[i] * s + zj[i];
  { D x(2.0), zj(0.75); seed(x, 0, N); seed(zj, 1 % N, N);
    const double s = 0.9;
    for (int k = 0; k < 6; ++k) { x = x * s + zj;
      emit((std::string(label) + " alias x=x*s+zj").c_str(), x, N); } }

  // cppde_newton.hpp:229
  { D ftemp(1.0), zn1(0.5), acor(0.25), tempv(0.0);
    seed(ftemp, 0, N); seed(zn1, 1 % N, N); seed(acor, 2 % N, N);
    const double rl1 = 0.8, inv_gamma = 1.3;
    for (int k = 0; k < 6; ++k) {
      tempv = ftemp - (rl1 * zn1 + acor) * inv_gamma;
      emit((std::string(label) + " newton tempv").c_str(), tempv, N); } }

  // cppde_rosenbrock4.hpp:501 -- the deepest nest in the codebase
  { D xo(1.0), xn(0.5), c3(0.25), c4(0.125), x(0.0);
    seed(xo, 0, N); seed(xn, 1 % N, N); seed(c3, 2 % N, N);
    const double s = 0.4, s1 = 0.6;
    for (int k = 0; k < 6; ++k) {
      x = xo * s1 + s * (xn + s1 * (c3 + s * c4));
      emit((std::string(label) + " rb4 dense-out").c_str(), x, N); } }
}

// =============================================================================
// Lifetime stress: materialise a deep tree, then clobber the stack and re-read
// a value captured beforehand. Under -fsanitize=address with
// detect_stack_use_after_scope=1 a node that outlived its operands becomes a
// hard failure instead of a plausible-looking wrong number.
// =============================================================================

template <class D>
static void lifetime_stress(const char* label, unsigned N) {
  D a(1.5), b(2.5), c(0.75), r;
  seed(a, 0, N); seed(b, 1 % N, N); seed(c, 2 % N, N);
  r = ((a * b + c) * (a - c) + exp(c) / (b + 1.0)) - pow(a, 2.0);
  const double kept = r.x();
  { volatile double scratch[512]; for (int i = 0; i < 512; ++i) scratch[i] = i * 3.5;
    (void)scratch[0]; }
  close_to(r.x(), kept, "value survived a stack clobber");
  emit((std::string(label) + " lifetime").c_str(), r, N);
}

int main() {
  std::printf("# cppde expression-template harness\n");

  arithmetic_suite<dual<double, 0>>("dyn", 3);
  arithmetic_suite<dual<double, 3>>("st3", 3);
  arithmetic_suite<dual<double, 1>>("st1", 1);

  aliasing_suite<dual<double, 0>>("dyn", 3);
  aliasing_suite<dual<double, 3>>("st3", 3);

  lifetime_stress<dual<double, 0>>("dyn", 3);
  lifetime_stress<dual<double, 3>>("st3", 3);

  std::printf("# failures: %d\n", g_failures);
  return g_failures == 0 ? 0 : 1;
}
