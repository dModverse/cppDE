/*
 Keep BLAS single-threaded for the duration of a solve.

 A generated model is compiled with -fopenmp and therefore runs on libgomp.
 A threaded BLAS brings its own OpenMP runtime along: the default threading
 layer of Intel MKL's single dynamic library is libiomp5, so a solve inside an
 R session linked against MKL ends up with two OpenMP runtimes live in one
 process, which is undefined behaviour.  It shows up as silently wrong
 sensitivity blocks, the step controller then sees a diverging error
 estimate and gives up, rather than the solve returning a wrong answer.
 (Confirmed by MKL_THREADING_LAYER=GNU and =SEQUENTIAL both making it go away,
 while MKL_NUM_THREADS=2 does not.)

 Pinning BLAS to one thread avoids the second runtime entirely.  cppDE gives up
 nothing by it: parallelism lives one level up, where solveODEBatch() runs whole
 conditions concurrently, and the per-solve BLAS calls are far too small to
 thread usefully anyway.

 The pin is scoped, not permanent, the thread count is process-global state
 that belongs to the caller, so a solve restores whatever was set before it.

 Copyright (C) 2026 Simon Beyer
 */

#ifndef CPPDE_BLAS_THREADS_HPP
#define CPPDE_BLAS_THREADS_HPP

// --- Runtime detection of the BLAS thread-control entry points ---
#ifdef _WIN32
#include <windows.h>
#else
#include <dlfcn.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

namespace cppde {
namespace detail {

// Resolved once per process.  Runtime APIs only: the environment variables are
// read by the BLAS at its own init, far too late to help here, and writing them
// with setenv() is neither effective nor thread-safe.
struct blas_thread_api {
  void (*mkl_set)(int)   = nullptr;
  int  (*mkl_get)()      = nullptr;
  void (*oblas_set)(int) = nullptr;
  int  (*oblas_get)()    = nullptr;

  blas_thread_api() {
#ifdef _WIN32
    HMODULE h = GetModuleHandle(NULL);
    mkl_set = (void(*)(int))GetProcAddress(h, "MKL_Set_Num_Threads");
    mkl_get = (int(*)())   GetProcAddress(h, "MKL_Get_Max_Threads");
#else
    mkl_set   = (void(*)(int))dlsym(RTLD_DEFAULT, "MKL_Set_Num_Threads");
    mkl_get   = (int(*)())    dlsym(RTLD_DEFAULT, "MKL_Get_Max_Threads");
    oblas_set = (void(*)(int))dlsym(RTLD_DEFAULT, "openblas_set_num_threads");
    oblas_get = (int(*)())    dlsym(RTLD_DEFAULT, "openblas_get_num_threads");
#endif
  }
};

inline const blas_thread_api& blas_api() {
  static const blas_thread_api api;
  return api;
}

// Set BLAS to one thread without recording the previous value.  Used where a
// restore would be unsafe: the thread count is process-global, so workers
// inside a parallel region must not race each other on putting it back.
// Writing the same value from every worker is harmless.
inline void ensure_single_thread_blas() {
  const blas_thread_api& api = blas_api();
  if (api.mkl_set)   api.mkl_set(1);
  if (api.oblas_set) api.oblas_set(1);
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
    const blas_thread_api& api = blas_api();
    if (api.mkl_set)   { if (api.mkl_get)   m_mkl   = api.mkl_get();   api.mkl_set(1); }
    if (api.oblas_set) { if (api.oblas_get) m_oblas = api.oblas_get(); api.oblas_set(1); }
  }

  ~single_thread_blas_scope() {
    const blas_thread_api& api = blas_api();
    if (m_mkl   > 1 && api.mkl_set)   api.mkl_set(m_mkl);
    if (m_oblas > 1 && api.oblas_set) api.oblas_set(m_oblas);
  }

  single_thread_blas_scope(const single_thread_blas_scope&)            = delete;
  single_thread_blas_scope& operator=(const single_thread_blas_scope&) = delete;

private:
  // -1 means "nothing to restore": no runtime API, no getter, or we were
  // inside a parallel region.
  int m_mkl   = -1;
  int m_oblas = -1;
};

}  // namespace detail
}  // namespace cppde

#endif  // CPPDE_BLAS_THREADS_HPP
