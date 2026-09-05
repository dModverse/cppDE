/*
 Keep BLAS single-threaded where a threaded BLAS would break the process.

 Two hazards share one lever.

 A generated model is compiled with -fopenmp and runs on libgomp, while the
 default threading layer of Intel MKL's single dynamic library is libiomp5.  Two
 live OpenMP runtimes in one process is undefined behaviour.  It shows up as
 silently wrong sensitivity blocks, after which the step controller sees a
 diverging error estimate and gives up, rather than as a solve that returns a
 wrong answer.  (Confirmed by MKL_THREADING_LAYER=GNU and =SEQUENTIAL both
 making it go away, while MKL_NUM_THREADS=2 does not.)

 A threaded BLAS also does not survive fork().  Only the forking thread reaches
 the child, and the first BLAS call there large enough to thread spins forever
 on a team whose workers stayed behind.  It needs a warm pool in the parent, so
 it misses small tests and hits real fits.

 Neither is reachable through the environment: a threading runtime reads
 OMP_NUM_THREADS when its library is loaded, before any R code runs, and
 setenv() afterwards is neither effective nor thread-safe.  The runtime entry
 points probed below are the only lever left.

 Pinning to one thread costs cppDE nothing: parallelism lives one level up,
 where solveODEBatch() runs whole conditions concurrently, and the per-solve
 BLAS calls are far too small to thread usefully anyway.

 Every pin is scoped.  The thread count is process-global state that belongs to
 the caller, so a solve restores what was set before it and the fork guard
 restores it in the parent as soon as fork() has returned.

 Copyright (C) 2026 Simon Beyer
 */

#ifndef CPPDE_BLAS_THREADS_HPP
#define CPPDE_BLAS_THREADS_HPP

#ifdef _WIN32
#include <windows.h>
// windows.h defines TRUE and FALSE as int macros, which shadow R's Rboolean
// constants for every translation unit that includes this header.
#undef TRUE
#undef FALSE
#else
#include <dlfcn.h>
#include <pthread.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

namespace cppde {
namespace detail {

// The most vendors that can be resolved at once.  Under a dispatcher every
// entry is an alias onto one counter, which costs a few redundant writes.
#define CPPDE_BLAS_MAX_VENDORS 4

// One BLAS implementation's runtime thread-count entry points.  `saved` is
// scratch for the fork guard, which has no stack frame to keep it on.
struct blas_vendor {
  const char* name = nullptr;
  void (*set)(int) = nullptr;
  int  (*get)()    = nullptr;
  int  saved       = -1;
};

inline void* find_symbol(const char* symbol) {
#ifdef _WIN32
  // The BLAS lives in a DLL, so the process image's own export table never
  // carries these.  Ask the modules that can, then the image as a fallback.
  static const char* modules[] = {"mkl_rt.dll", "mkl_rt.2.dll", "flexiblas.dll",
                                  "libopenblas.dll", "Rblas.dll", "libblas.dll",
                                  nullptr};
  for (int i = 0; modules[i] != nullptr; ++i) {
    HMODULE h = GetModuleHandleA(modules[i]);
    if (h == nullptr) continue;
    FARPROC p = GetProcAddress(h, symbol);
    if (p != nullptr) return reinterpret_cast<void*>(p);
  }
  return reinterpret_cast<void*>(GetProcAddress(GetModuleHandleA(nullptr), symbol));
#else
  return dlsym(RTLD_DEFAULT, symbol);
#endif
}

// Resolved once per process.  Runtime APIs only, see the file header on why the
// environment variables are not an option.
struct blas_thread_api {
  blas_vendor vendor[CPPDE_BLAS_MAX_VENDORS];
  int n = 0;

  void add(const char* nm, const char* set_sym, const char* get_sym) {
    if (n >= CPPDE_BLAS_MAX_VENDORS) return;
    void* s = find_symbol(set_sym);
    if (s == nullptr) return;
    vendor[n].name = nm;
    vendor[n].set  = reinterpret_cast<void (*)(int)>(s);
    vendor[n].get  = reinterpret_cast<int (*)()>(find_symbol(get_sym));
    ++n;
  }

  blas_thread_api() {
    // FlexiBLAS first.  It dispatches to whichever backend is loaded and
    // exports the other vendors' names as aliases onto its own counter, so
    // wherever it is in play its own name is the honest one to report.
    add("FlexiBLAS", "flexiblas_set_num_threads",  "flexiblas_get_num_threads");
    add("MKL",       "MKL_Set_Num_Threads",        "MKL_Get_Max_Threads");
    add("OpenBLAS",  "openblas_set_num_threads",   "openblas_get_num_threads");
    add("BLIS",      "bli_thread_set_num_threads", "bli_thread_get_num_threads");
  }
};

inline blas_thread_api& blas_api() {
  static blas_thread_api api;
  return api;
}

// OpenMP's own thread count.  A BLAS that threads through OpenMP takes its team
// width from here and not from the vendor counter, so pinning the vendor alone
// leaves a forked child spinning.  Fork guard only: a solve must not touch this
// lever, cppDE's batch parallelism is that same runtime.
struct omp_thread_api {
  void (*set)(int) = nullptr;
  int  (*get)()    = nullptr;
  int  saved       = -1;

  omp_thread_api() {
    set = reinterpret_cast<void (*)(int)>(find_symbol("omp_set_num_threads"));
    get = reinterpret_cast<int (*)()>(find_symbol("omp_get_max_threads"));
  }
};

inline omp_thread_api& omp_api() {
  static omp_thread_api api;
  return api;
}

// The BLAS whose thread count cppDE steers, or nullptr if nothing resolved.
inline const char* blas_name() {
  const blas_thread_api& api = blas_api();
  return (api.n > 0) ? api.vendor[0].name : nullptr;
}

// Current thread count, or -1 when no getter resolved.
inline int blas_threads() {
  const blas_thread_api& api = blas_api();
  for (int i = 0; i < api.n; ++i)
    if (api.vendor[i].get != nullptr) return api.vendor[i].get();
  return -1;
}

// Set BLAS to one thread without recording the previous value.  Used where a
// restore would be unsafe: the thread count is process-global, so workers
// inside a parallel region must not race each other on putting it back.
// Writing the same value from every worker is harmless.
inline void ensure_single_thread_blas() {
  blas_thread_api& api = blas_api();
  for (int i = 0; i < api.n; ++i) api.vendor[i].set(1);
}

  // RAII: pin BLAS to one thread and restore the caller's setting on exit.
  //
  // A no-op inside an OpenMP region: the batch entry installs one scope around
  // the whole region instead. In a region that is not ours there is no serial
  // phase to hook, so the setting is applied without being restored.
class single_thread_blas_scope {
public:
  single_thread_blas_scope() {
#ifdef _OPENMP
    if (omp_in_parallel()) { ensure_single_thread_blas(); return; }
#endif
    blas_thread_api& api = blas_api();
    m_n = api.n;
    for (int i = 0; i < m_n; ++i) {
      if (api.vendor[i].get != nullptr) m_saved[i] = api.vendor[i].get();
      api.vendor[i].set(1);
    }
  }

  ~single_thread_blas_scope() {
    blas_thread_api& api = blas_api();
    for (int i = 0; i < m_n; ++i)
      if (m_saved[i] > 1) api.vendor[i].set(m_saved[i]);
  }

  single_thread_blas_scope(const single_thread_blas_scope&)            = delete;
  single_thread_blas_scope& operator=(const single_thread_blas_scope&) = delete;

private:
  // -1 means "nothing to restore": no runtime API, no getter, or we were
  // inside a parallel region.
  int m_saved[CPPDE_BLAS_MAX_VENDORS] = {-1, -1, -1, -1};
  int m_n = 0;
};

  // Pin BLAS to one thread for the width of a fork() and restore the caller's
  // count in the parent, so a child never meets a pool whose workers did not
  // come along.  Installed once per process, from R_init_cppDE().
inline bool install_blas_fork_guard() {
#ifndef _WIN32
  blas_api();  // resolve here: dlsym() has no business running inside a handler
  omp_api();
  ::pthread_atfork(
      +[]() {
        blas_thread_api& api = blas_api();
        for (int i = 0; i < api.n; ++i) {
          api.vendor[i].saved = (api.vendor[i].get != nullptr) ? api.vendor[i].get() : -1;
          api.vendor[i].set(1);
        }
        omp_thread_api& omp = omp_api();
        if (omp.set == nullptr) return;
        omp.saved = (omp.get != nullptr) ? omp.get() : -1;
        omp.set(1);
      },
      +[]() {
        blas_thread_api& api = blas_api();
        for (int i = 0; i < api.n; ++i)
          if (api.vendor[i].saved > 1) api.vendor[i].set(api.vendor[i].saved);
        omp_thread_api& omp = omp_api();
        if (omp.set != nullptr && omp.saved > 1) omp.set(omp.saved);
      },
      // No child handler: only async-signal-safe calls are allowed there and a
      // vendor setter is not one, OpenBLAS's joins the threads that are gone.
      // The child inherits the pin from the parent's memory image.
      nullptr);
  return true;
#else
  return false;  // no fork() on Windows
#endif
}

}  // namespace detail
}  // namespace cppde

#endif  // CPPDE_BLAS_THREADS_HPP
