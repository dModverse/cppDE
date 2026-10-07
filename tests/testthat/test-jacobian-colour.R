# The coloured Jacobian (CPPDE_JAC = "colour") against the entry-wise one, for
# both backends, forward sensitivities and the reverse sweep.

skip_on_cran()

n_col <- 12
eqns_col <- setNames(vapply(seq_len(n_col), function(i) if (i == 1) "-k1*x1 + k2*x2" else
  sprintf("k1*x%d*x%d/(1 + x%d) - k2*x%d", i - 1, i - 1, i - 1, i), ""), paste0("x", seq_len(n_col)))
pars_col <- c(setNames(c(1, rep(0, n_col - 1)), names(eqns_col)), k1 = 0.7, k2 = 0.3)
times_col <- seq(0, 10, 1)

build_col <- function(strategy, backend, name, ...) {
  old <- Sys.getenv("CPPDE_JAC", NA)
  on.exit(if (is.na(old)) Sys.unsetenv("CPPDE_JAC") else Sys.setenv(CPPDE_JAC = old))
  Sys.setenv(CPPDE_JAC = strategy)
  backend(eqns_col, modelname = paste0("jac_", name, "_", strategy), compile = FALSE, ...)
}

mods_col <- list(
  cpp_c = build_col("colour", cppODE, "cpp"), cpp_e = build_col("entries", cppODE, "cpp"),
  rev_c = build_col("colour", cppODE, "rev", derivMode = "reverse"),
  rev_e = build_col("entries", cppODE, "rev", derivMode = "reverse"))
if (isTRUE(cvodeConfig$available))
  mods_col <- c(mods_col, list(cv_c = build_col("colour", cvode, "cv"), cv_e = build_col("entries", cvode, "cv")))
do.call(compile, c(unname(mods_col), list(output = "test_jacobian_colour", cores = test_cores())))

test_that("coloured and entry-wise Jacobians give the same forward sensitivities", {
  a <- solveODE(mods_col$cpp_c, times_col, pars_col, abstol = 1e-10, reltol = 1e-10)
  b <- solveODE(mods_col$cpp_e, times_col, pars_col, abstol = 1e-10, reltol = 1e-10)
  expect_equal(a$variable, b$variable, tolerance = 1e-8)
  expect_equal(a$tangent, b$tangent, tolerance = 1e-8)
})

test_that("coloured and entry-wise Jacobians give the same reverse gradient", {
  W <- matrix(1, length(times_col), n_col)
  a <- solveODE(mods_col$rev_c, times_col, pars_col, abstol = 1e-10, reltol = 1e-10, cotangent = W)
  b <- solveODE(mods_col$rev_e, times_col, pars_col, abstol = 1e-10, reltol = 1e-10, cotangent = W)
  expect_length(a$cotangent, n_col + 2)
  expect_equal(a$cotangent, b$cotangent, tolerance = 1e-8)
})

test_that("the coloured CVODE Jacobian matches the entry-wise one", {
  skip_if_not(isTRUE(cvodeConfig$available), "CVODE not available")
  a <- solveODE(mods_col$cv_c, times_col, pars_col, abstol = 1e-10, reltol = 1e-10)
  b <- solveODE(mods_col$cv_e, times_col, pars_col, abstol = 1e-10, reltol = 1e-10)
  expect_equal(a$variable, b$variable, tolerance = 1e-8)
  expect_equal(a$tangent, b$tangent, tolerance = 1e-8)
})
