#' Report the BLAS Fork Guard
#'
#' @description
#' Reports the guard cppDE installs so that a threaded BLAS does not deadlock in
#' a forked child, as under [parallel::mclapply()]. The guard sets BLAS to one
#' thread for the width of every `fork()` and restores the previous count in
#' the parent.
#'
#' @details
#' BLAS also runs on one thread during each solve. Outside a solve and a
#' `fork()` the thread count is the caller's own; a forked child keeps the
#' single thread it inherited.
#'
#' `api` is `NA` on a BLAS with no runtime thread-control entry point, and the
#' deadlock is then still reachable; the fallback is `OMP_NUM_THREADS=1` in the
#' environment that starts R. `guard` is `FALSE` on Windows, which has no
#' `fork()`.
#'
#' Attaching the package prints this in one line. Set
#' `options(cppDE.quiet = TRUE)` beforehand to suppress the line, or attach with
#' [suppressPackageStartupMessages()].
#'
#' @return List with `api` (the BLAS whose thread count cppDE steers, `NA` when
#'   no entry point resolved), `threads` (its current thread count, `NA` without
#'   a getter) and `guard` (whether the `fork()` handler is installed).
#' @seealso [batchAvailable()]
#' @example inst/examples/forkGuard.R
#' @export
forkGuard <- function() .Call(C_cppde_blas_info)
