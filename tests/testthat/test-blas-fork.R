# The BLAS fork guard: BLAS pinned to one thread for the width of every fork()
# and restored in the parent. A regression here deadlocks a forked worker rather
# than returning a wrong answer, so the fork case runs in its own process under
# a timeout and asserts on the exit status.

skip_on_cran()

reported_threads <- function() forkGuard()$threads

# A separate R process with this session's library path, combined output back.
run_r <- function(code, timeout = 180) {
  script <- tempfile(fileext = ".R")
  writeLines(c(paste0(".libPaths(", paste(deparse(.libPaths()), collapse = ""), ")"),
               code), script)
  on.exit(unlink(script), add = TRUE)
  suppressWarnings(
    system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(script)),
            stdout = TRUE, stderr = TRUE, timeout = timeout))
}

test_that("forkGuard reports the three fields", {
  g <- forkGuard()
  expect_named(g, c("api", "threads", "guard"))
  expect_type(g$api, "character")
  expect_type(g$threads, "integer")
  expect_type(g$guard, "logical")
  expect_length(g$guard, 1L)
})

test_that("the guard is installed where there is a fork to guard", {
  skip_on_os("windows")
  expect_true(forkGuard()$guard)
})

test_that("Windows reports no guard, having no fork", {
  skip_if(.Platform$OS.type != "windows", "not Windows")
  expect_false(forkGuard()$guard)
})

test_that("a fork leaves the BLAS thread count where it was", {
  skip_on_os("windows")
  before <- reported_threads()
  skip_if(is.na(before), "no BLAS thread-count getter resolved")

  invisible(parallel::mclapply(1:2, function(i) i, mc.cores = 2))

  expect_identical(reported_threads(), before)
})

test_that("a forked worker survives a BLAS call the parent warmed", {
  skip_on_os("windows")
  n <- reported_threads()
  skip_if(is.na(n) || n < 2L, "BLAS reports a single thread here")

  out <- run_r(c(
    'library(cppDE)',
    'x <- matrix(rnorm(600 * 600), 600); invisible(x %*% x)',
    'r <- parallel::mclapply(1:2, function(i) {',
    '  m <- matrix(rnorm(200 * 200), 200); sum(m %*% m)',
    '}, mc.cores = 2)',
    'if (length(unlist(r)) != 2L) stop("worker returned nothing")'))

  expect_identical(attr(out, "status"), NULL)
})

test_that("the startup line comes on attach, not on a namespace load", {
  skip_on_os("windows")
  n <- reported_threads()
  skip_if(is.na(n) || n < 2L, "nothing would be reported")

  expect_true(any(grepl("^cppDE:", run_r('library(cppDE)'))))
  expect_false(any(grepl("^cppDE:", run_r('loadNamespace("cppDE")'))))
  expect_false(any(grepl("^cppDE:",
                         run_r('options(cppDE.quiet = TRUE); library(cppDE)'))))
})
