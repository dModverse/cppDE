/*
 What a reverse solve actually costs in memory.

 The trajectory harness reports the checkpoint store as
 n_steps * (n_x + 2) * sizeof(double), which is the payload and not the
 footprint: every checkpoint holds its own std::vector, so each step is a
 separate allocation with its own header and its own vector object. On a small
 model that overhead is larger than the data.

 This counts what is really asked of the allocator, by replacing the global
 operator new. Reported against the payload, which is the number a design
 discussion would use.

   build: dev/cxx/run.sh --bench-revmem
 */

#include <cppde/cppde.hpp>
#include <cppde/cppde_codual.hpp>
#include <cppde/cppde_codual_math.hpp>
#include <cppde/cppde_reverse_trajectory.hpp>

#include <cstdio>
#include <cstdlib>
#include <new>
#include <vector>

namespace {
std::size_t g_bytes = 0, g_calls = 0, g_live = 0, g_peak = 0;
bool g_on = false;
}

void* operator new(std::size_t n) {
  void* p = std::malloc(n ? n : 1);
  if (!p) throw std::bad_alloc();
  if (g_on) {
    g_bytes += n; ++g_calls; g_live += n;
    if (g_live > g_peak) g_peak = g_live;
  }
  return p;
}
void operator delete(void* p) noexcept { std::free(p); }
void operator delete(void* p, std::size_t n) noexcept {
  if (g_on && g_live >= n) g_live -= n;
  std::free(p);
}

using cppde::codual;

template<class T>
static std::pair<cppde::vector_t<T>, int> make_system_dummy() { return {}; }

int main() {
  using Stepper = cppde::tsit5<double>;
  using Store   = cppde::reverse::trajectory_store<Stepper, double>;

  const int NX_LIST[] = {2, 5, 20, 100};
  const std::size_t N_STEPS = 2000;

  std::printf("A store of %zu accepted steps, states across the columns.\n\n",
              N_STEPS);
  std::printf("%-8s %-14s %-14s %-10s %-12s\n",
              "n_x", "payload KiB", "asked KiB", "overhead", "allocations");

  for (int nx : NX_LIST) {
    Stepper st;
    std::vector<double> x((std::size_t)nx, 1.0);

    g_bytes = g_calls = g_live = g_peak = 0;
    g_on = true;
    {
      Store store;
      for (std::size_t k = 0; k < N_STEPS; ++k)
        store.capture(st, x, 0.01 * (double)k, 0.01);
      // Keep it alive to the end of the scope, which is what a real solve does.
      if (store.n_steps() != N_STEPS) std::printf("  (unexpected step count)\n");
    }
    g_on = false;

    const double payload = (double)N_STEPS * ((double)nx + 2.0) * sizeof(double);
    std::printf("%-8d %-14.1f %-14.1f %-10.2f %-12zu\n",
                nx, payload / 1024.0, (double)g_peak / 1024.0,
                (double)g_peak / payload, g_calls);
  }

  std::printf("\n`asked` is the peak the allocator was holding for the store.\n"
              "One allocation per step means the per-allocation header and the\n"
              "vector object itself are paid n_steps times over.\n");
  return 0;
}
