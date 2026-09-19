# Test that all three backends accept rhs/eqns with zero parameters.
# Pure-state systems and literal-only expressions must not break codegen,
# AD slab sizing, or sensitivity output shape.

skip_on_cran()

# The native and cppFUN models, compiled into one shared object.
mod_nd   <- cppODE(c(x = "-x"), modelname = "noparm_cpp_nd", deriv = FALSE,
                   compile = FALSE)
mod_d    <- cppODE(c(x = "-x"), modelname = "noparm_cpp_d", deriv = TRUE,
                   compile = FALSE)
# convenient = FALSE so we can pass an explicit (n_obs, 0) matrix; the
# convenient wrapper has no way to express n_obs when there are no vars.
fun_lit  <- cppFUN(c(y = "5"), compile = FALSE, modelname = "noparm_fun_lit",
                   convenient = FALSE)
fun_dual <- cppFUN(c(y = "2*x + 3"), compile = FALSE,
                   modelname = "noparm_fun_dual", derivMode = "forward")
fun_d2   <- cppFUN(c(y = "x^2 + 3*x"), compile = FALSE,
                   modelname = "noparm_fun_dual_d2",
                   derivMode = "forward", deriv2 = TRUE)
compile(mod_nd, mod_d, fun_lit, fun_dual, fun_d2,
        output = "test_no_parameters", cores = 1)

if (isTRUE(cvodeConfig$available)) {
  cv_nd <- cvode(c(x = "-x"), modelname = "noparm_cv_nd", compile = FALSE)
  cv_d  <- cvode(c(x = "-x"), modelname = "noparm_cv_d", deriv = TRUE,
                 compile = FALSE)
  compile(cv_nd, cv_d, output = "test_no_parameters_cvode", cores = 1)
}

# -- cppDE: pure decay, no parameters ----------------------------------------

test_that("cppDE compiles and solves with zero parameters (deriv = FALSE)", {
  mod <- mod_nd

  expect_equal(attr(mod, "parameters"), character(0))

  tvec <- seq(0, 2, by = 0.5)
  res  <- solveODE(mod, times = tvec, parms = c(x = 1),
                   abstol = 1e-10, reltol = 1e-10)

  expect_equal(as.numeric(res$variable[, "x"]),
               exp(-tvec), tolerance = 1e-8)
})

test_that("cppDE deriv = TRUE with zero parameters seeds initial-state sens", {
  mod <- mod_d

  expect_equal(attr(mod, "parameters"), character(0))
  expect_equal(attr(mod, "dimNames")$sens, "x")

  tvec <- seq(0, 2, by = 0.5)
  res  <- solveODE(mod, times = tvec, parms = c(x = 1),
                   abstol = 1e-10, reltol = 1e-10)

  expect_equal(dim(res$sens1), c(length(tvec), 1L, 1L))
  # dx(t)/dx0 = exp(-t) for dx/dt = -x
  expect_equal(as.numeric(res$sens1[, 1, 1]),
               exp(-tvec), tolerance = 1e-8)
})

# -- CVODE: same model, same expectations -------------------------------------

test_that("CVODE compiles and solves with zero parameters (deriv = FALSE)", {
  skip_if_not(isTRUE(cvodeConfig$available), "CVODE backend not available")

  mod <- cv_nd

  expect_equal(attr(mod, "parameters"), character(0))

  tvec <- seq(0, 2, by = 0.5)
  res  <- solveODE(mod, times = tvec, parms = c(x = 1),
                   abstol = 1e-10, reltol = 1e-10)

  expect_equal(as.numeric(res$variable[, "x"]),
               exp(-tvec), tolerance = 1e-8)
})

test_that("CVODE deriv = TRUE with zero parameters seeds initial-state sens", {
  skip_if_not(isTRUE(cvodeConfig$available), "CVODE backend not available")

  mod <- cv_d

  expect_equal(attr(mod, "parameters"), character(0))
  expect_equal(attr(mod, "dimNames")$sens, "x")

  tvec <- seq(0, 2, by = 0.5)
  res  <- solveODE(mod, times = tvec, parms = c(x = 1),
                   abstol = 1e-10, reltol = 1e-10)

  expect_equal(dim(res$sens1), c(length(tvec), 1L, 1L))
  expect_equal(as.numeric(res$sens1[, 1, 1]),
               exp(-tvec), tolerance = 1e-8)
})

# -- cppFUN: literal-only equations (no variables, no parameters) -------------

test_that("cppFUN accepts literal-only equations", {
  obj <- fun_lit

  expect_equal(attr(obj, "variables"), character(0))
  expect_null(attr(obj, "parameters"))

  res <- obj$func(matrix(numeric(0), nrow = 3L, ncol = 0L))
  expect_true(is.matrix(res))
  expect_equal(dim(res), c(3L, 1L))
  expect_equal(as.numeric(res[, "y"]), rep(5, 3))
})

# -- cppFUN: state variables only, no parameters -------------------------------

test_that("cppFUN dual mode evaluates with zero parameters", {
  obj <- fun_dual

  expect_equal(attr(obj, "variables"), "x")
  expect_null(attr(obj, "parameters"))

  res <- obj$func(x = c(1, 2, 3))
  expect_equal(as.numeric(res[, "y"]), c(5, 7, 9))

  expect_true(!is.null(obj$jac))
  n_obs <- 3L
  dX <- array(0, c(n_obs, 1L, 1L), list(NULL, "x", "x"))
  dX[, "x", "x"] <- 1
  dP <- matrix(numeric(0), nrow = 0L, ncol = 1L,
               dimnames = list(NULL, "x"))

  jac <- obj$jac(x = c(1, 2, 3), dX = dX, dP = dP)
  # dy/dx = 2 along the lone theta = "x"
  expect_equal(as.numeric(jac[, "y", "x"]), rep(2, n_obs))

  # Combined evaluate() returns y and dy in one nested-dual pass.
  ev <- obj$evaluate(x = c(1, 2, 3), dX = dX, dP = dP)
  expect_equal(as.numeric(ev$y[, "y"]),         c(5, 7, 9))
  expect_equal(as.numeric(ev$dy[, "y", "x"]),   rep(2, n_obs))
})

test_that("cppFUN dual mode supports deriv2 with zero parameters", {
  obj <- fun_d2

  # raw hess (identity seed) at x = c(1, 2)
  hess <- obj$hess(x = c(1, 2))
  expect_equal(dim(hess), c(2L, 1L, 1L, 1L))
  # d2y/dx2 = 2
  expect_equal(as.numeric(hess[, "y", "x", "x"]), c(2, 2))

  # explicit identity seed via dX should give the same.
  n_obs <- 2L
  dX <- array(0, c(n_obs, 1L, 1L), list(NULL, "x", "x"))
  dX[, "x", "x"] <- 1
  hess2 <- obj$hess(x = c(1, 2), dX = dX)
  expect_equal(hess, hess2)
})

test_that("cppFUN forward mode produces correct jac/hess with zero parameters", {
  obj <- fun_d2

  expect_equal(attr(obj, "variables"), "x")
  expect_null(attr(obj, "parameters"))

  jac  <- obj$jac(x = c(1, 2))
  hess <- obj$hess(x = c(1, 2))

  expect_equal(dim(jac),  c(2L, 1L, 1L))
  expect_equal(dim(hess), c(2L, 1L, 1L, 1L))
  # dy/dx = 2x + 3
  expect_equal(as.numeric(jac[, "y", "x"]), c(5, 7))
  # d2y/dx2 = 2
  expect_equal(as.numeric(hess[, "y", "x", "x"]), c(2, 2))
})
