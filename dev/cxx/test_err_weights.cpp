/*
 The lambda weight carrier, stage 9 of dev/adjoint-plan.md.

 Checked here rather than through R because the interesting cases are the ones
 an R-level test cannot reach cheaply: what a query does inside a broken
 interval, what a size mismatch does, and that the sink is null outside a scope
 so an unweighted run pays nothing and cannot inherit another solve's weights.
 */

#include <cppde/cppde_err_weights.hpp>
#include <cstdio>
#include <cmath>
int fails = 0;
static void chk(bool ok, const char* m) { if (!ok) { std::printf("FAIL %s\n", m); ++fails; } }
int main() {
  cppde::err_weights w;
  chk(w.empty(), "empty by default");

  // three samples, two states
  w.set({0.0, 1.0, 2.0}, {1,10, 2,20, 4,40}, 2);
  chk(!w.empty(), "not empty after set");
  std::vector<double> o;
  w.at(0.5, o);  chk(std::fabs(o[0]-1.5)<1e-12 && std::fabs(o[1]-15)<1e-12, "linear midpoint");
  w.at(-1.0, o); chk(o[0]==1 && o[1]==10, "clamped left");
  w.at(9.0, o);  chk(o[0]==4 && o[1]==40, "clamped right");
  w.at(2.0, o);  chk(o[0]==4 && o[1]==40, "right endpoint");

  // a break at sample 1 must stop the interpolant spanning [0,1]
  cppde::err_weights b;
  b.set({0.0, 1.0, 2.0}, {1,10, 2,20, 4,40}, 2, {1});
  b.at(0.2, o); chk(o[0]==1, "break: nearer left sample");
  b.at(0.9, o); chk(o[0]==2, "break: nearer right sample");
  b.at(1.5, o); chk(std::fabs(o[0]-3.0)<1e-12, "unbroken interval still blends");

  // mismatched sizes are dropped whole, not half-read
  cppde::err_weights bad;
  bad.set({0.0, 1.0}, {1,2,3}, 2);
  chk(bad.empty(), "size mismatch drops");

  // the floor lifts small components toward the largest
  cppde::err_weights f;
  f.set({0.0}, {1.0, 0.001}, 2);
  f.floor(0.1);
  f.at(0.0, o); chk(std::fabs(o[1]-0.1)<1e-12, "floor lifts");

  // the sink is null until a scope is open, and restores after
  chk(cppde::err_weight_sink() == nullptr, "sink starts null");
  {
    cppde::err_weight_scope s(w);
    chk(cppde::err_weight_sink() == &w, "scope sets sink");
    std::vector<double> e{1.0, 0.0};
    const double q = cppde::detail::weighted_error(e, 0.0, [](double v){ return v; });
    // lambda(0) = (1, 10), e = (1, 0) -> |1| / gradtol(1e-6)
    chk(std::fabs(q - 1e6) < 1e-3, "weighted_error contracts");
  }
  chk(cppde::err_weight_sink() == nullptr, "scope restores sink");

  std::vector<double> e{1.0, 0.0};
  chk(cppde::detail::weighted_error(e, 0.0, [](double v){ return v; }) == 0.0,
      "no sink means no term");

  std::printf(fails == 0 ? "OK\n" : "%d FAILURES\n", fails);
  return fails == 0 ? 0 : 1;
}
