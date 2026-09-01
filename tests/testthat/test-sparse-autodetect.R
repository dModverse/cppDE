# Pins the dense/sparse auto-detection thresholds (>= 8 states, Jacobian
# density <= 0.4), which live in `decide_sparse()` in
# inst/python/codegen_cppODE.py and are shared by both backends.
#
# The models below have an exactly known structural density, so a threshold
# change shows up as a failure rather than as a slow benchmark. `compile =
# FALSE` is enough: the decision is made in codegen and recorded on the
# returned handle as attr(m, "sparse").

skip_on_cran()

# Without KLU the R layer downgrades an auto decision to dense before codegen
# ever sees it (R/cppODE.R, R/cvode.R), so the sparse expectations below would
# fail for a reason that has nothing to do with the thresholds.
skip_if_not(isTRUE(cppDE:::cvodeConfig$klu_available),
            "KLU not available; auto-detection is forced to dense")

lu_of <- function(mod) if (isTRUE(attr(mod, "sparse"))) "sparse" else "dense"

# Both backends must reach the same decision -- they share decide_sparse(), and
# a divergence is exactly the drift the shared helper exists to prevent.
expect_lu <- function(rhs, expected, label, ...) {
  a <- cppODE(rhs, deriv = FALSE, compile = FALSE, verbose = FALSE,
              outdir = tempdir(), modelname = paste0("autodet_cpp_", label), ...)
  b <- cvode(rhs, deriv = FALSE, compile = FALSE, verbose = FALSE,
             outdir = tempdir(), modelname = paste0("autodet_cv_", label), ...)
  expect_identical(lu_of(a), expected, info = paste("cppDE:", label))
  expect_identical(lu_of(b), expected, info = paste("CVODE:", label))
}

# -- model builders with an exactly known Jacobian density ---------------------

# Bidiagonal chain: dx_i = k*x_{i-1} - k*x_i.  nnz = 2n-1, density (2n-1)/n^2.
chain_rhs <- function(n) {
  eqs <- vapply(seq_len(n), function(i)
    if (i == 1L) "-k*x1" else sprintf("k*x%d - k*x%d", i - 1L, i), "")
  stats::setNames(eqs, paste0("x", seq_len(n)))
}

# Each row depends on exactly `w` states (cyclically), so nnz = n*w and the
# density is exactly w/n -- which is what pins the 0.4 cutoff.
band_rhs <- function(n, w) {
  eqs <- vapply(seq_len(n), function(i) {
    idx <- ((i - 1L + seq_len(w) - 1L) %% n) + 1L
    paste0("-k*(", paste0("x", idx, collapse = " + "), ")")
  }, "")
  stats::setNames(eqs, paste0("x", seq_len(n)))
}

# -- the state-count threshold ------------------------------------------------

test_that("auto-detection stays dense below 8 states", {
  # n = 7, density 13/49 = 0.265: sparse enough, but too few states.
  expect_lu(chain_rhs(7), "dense", "chain7")
})

test_that("auto-detection switches to sparse from 8 states", {
  # n = 8, density 15/64 = 0.234.
  expect_lu(chain_rhs(8), "sparse", "chain8")
})

# -- the density threshold ----------------------------------------------------

test_that("auto-detection stays dense for a dense Jacobian", {
  # n = 12, every row depends on every state: density 1.0.
  expect_lu(band_rhs(12, 12), "dense", "full12")
})

test_that("auto-detection accepts density exactly at the 0.4 cutoff", {
  # n = 10, 4 states per row: density exactly 0.4, and the rule is `<= 0.4`.
  expect_lu(band_rhs(10, 4), "sparse", "band10w4")
})

test_that("auto-detection rejects density just above the cutoff", {
  # n = 10, 5 states per row: density 0.5.
  expect_lu(band_rhs(10, 5), "dense", "band10w5")
})

# -- explicit pinning still wins ----------------------------------------------

test_that("an explicit sparse/dense argument overrides the thresholds", {
  # Small and dense -- auto would say dense.
  expect_lu(chain_rhs(3), "sparse", "pin_sparse", sparse = TRUE)
  # Large and sparse -- auto would say sparse.
  expect_lu(chain_rhs(40), "dense", "pin_dense", sparse = FALSE)
})

# -- explicit methods have no Jacobian ----------------------------------------

test_that("an explicit method stays dense however sparse the system looks", {
  # tsit5 skips Jacobian generation, so the nnz count is 0. Read naively that
  # is a perfectly sparse matrix and, above 8 states, would select KLU for a
  # system that has no iteration matrix at all.
  mod <- cppODE(chain_rhs(40), deriv = FALSE, method = "tsit5", compile = FALSE,
                verbose = FALSE, outdir = tempdir(), modelname = "autodet_tsit5")
  expect_identical(lu_of(mod), "dense")
})


# --------------------------------------------------------------------------
# Sparse + first-order sensitivities
#
# A solve keeps BLAS single-threaded throughout, because a threaded MKL brings
# libiomp5 alongside libgomp and corrupts the AD tangent blocks. The pin has to
# be independent of the dense LU, which a sparse model never reaches.

test_that("sparse Jacobian and dense Jacobian agree on first-order sensitivities", {
  # Chain of 8 states: Jacobian is bidiagonal, so auto-detection picks sparse.
  n   <- 8L
  nms <- paste0("x", seq_len(n))
  rhs <- setNames(c("-k1*x1",
                    paste0("k", seq_len(n - 1L), "*x", seq_len(n - 1L), " - k",
                           seq(2L, n), "*x", seq(2L, n))[seq_len(n - 1L)]),
                  nms)
  rhs[[n]] <- paste0("k", n - 1L, "*x", n - 1L)

  parms <- c(setNames(c(10, rep(0, n - 1L)), nms),
             setNames(seq(0.7, by = 0.3, length.out = n - 1L),
                      paste0("k", seq_len(n - 1L))))
  tt <- c(0, 0.5, 1, 2, 4, 8)

  sp <- cppODE(rhs, modelname = "sens_sparse", deriv = TRUE, sparse = TRUE,
               outdir = tempdir(), verbose = FALSE)
  # Sparse first, before anything has run a dense LU in this session.
  out_sp <- solveODE(sp, tt, parms)

  dn <- cppODE(rhs, modelname = "sens_dense", deriv = TRUE, sparse = FALSE,
               outdir = tempdir(), verbose = FALSE)
  out_dn <- solveODE(dn, tt, parms)

  expect_identical(out_sp$time, tt)
  expect_identical(out_sp$diagnostics$return_code, 0L)
  expect_equal(out_sp$variable, out_dn$variable, tolerance = 1e-6)
  expect_equal(out_sp$sens1, out_dn$sens1, tolerance = 1e-6)
})

test_that("an incomplete integration is an error, not partial results", {
  m <- cppODE(c(A = "-k*A"), modelname = "sparse_onfailure", deriv = FALSE,
              outdir = tempdir(), verbose = FALSE)
  tt <- c(0, 1, 2)
  p  <- c(A = 1, k = 1)

  expect_error(solveODE(m, tt, p, maxsteps = 2L), "did not complete")
  partial <- expect_warning(solveODE(m, tt, p, maxsteps = 2L, onFailure = "warn"),
                            "did not complete")
  expect_lt(length(partial$time), length(tt))
  expect_silent(solveODE(m, tt, p, maxsteps = 2L, onFailure = "silent"))
})
