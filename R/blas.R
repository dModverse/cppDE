#' Report the BLAS Fork Guard
#'
#' @description
#' A threaded BLAS deadlocks in a forked child: its worker threads do not come
#' along, and the first call large enough to thread waits forever on a lock they
#' held. cppDE pins BLAS to one thread for the width of every `fork()` and
#' restores the previous count in the parent, so [parallel::mclapply()] and
#' anything built on it are safe without the caller wrapping them. This reports
#' what the guard found.
#'
#' @details
#' Nothing is pinned at load time: outside a `fork()` the thread count is the
#' caller's own. A child keeps the single thread it inherited, so BLAS inside a
#' forked worker is serial.
#'
#' `api` is `NA` on a BLAS with no runtime thread-control entry point, and the
#' deadlock is then still reachable; the fallback is `OMP_NUM_THREADS=1` in the
#' environment that starts R, because setting it from R is too late. `guard` is
#' `FALSE` on Windows, which has no `fork()`.
#'
#' Attaching the package prints this in one line. The guard itself is installed
#' when the package DLL loads, so it is in place either way. Set
#' `options(cppDE.quiet = TRUE)` beforehand to suppress the line, or attach with
#' [suppressPackageStartupMessages()].
#'
#' @return List with `api` (the BLAS whose thread count cppDE steers, `NA` when
#'   no entry point resolved), `threads` (its current thread count, `NA` without
#'   a getter) and `guard` (whether the `fork()` handler is installed).
#' @seealso [batchAvailable()]
#' @export
forkGuard <- function() .Call(C_cppde_blas_info)
