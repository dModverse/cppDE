# The saltation of a switch: taken however slowly its root crosses, and once
# per switch where two cross at the same time but move apart with the
# parameters.

skip_on_cran()

library(cppDE)

sa_rhs <- c(x = "c", y = "piecewise(1, x > 1, 0)", z = "piecewise(1, time > ts, 0)",
            u = "e", v = "piecewise(1, u > e, 0)")
sa <- function(nm, ...) cppODE(sa_rhs, modelname = nm, compile = FALSE, ...)
sa_mods <- list(d2 = sa("sa_d2", deriv2 = TRUE), fr = sa("sa_fr", derivMode = "forward-reverse"))
has_cvode <- isTRUE(cvodeConfig$available)
if (has_cvode)
  sa_mods$cv <- cvode(sa_rhs, modelname = "sa_cv", deriv = TRUE, compile = FALSE)
do.call(compile, c(unname(sa_mods), list(output = "test_switch_saltation",
                                         cores = test_cores())))

sa_times <- c(0, 0.5, 1.5, 2)
sa_p <- c(x = 0.5, y = 0, z = 0, u = 0, v = 0, c = 0.5, ts = 1, e = 1e-17)
sa_on <- as.numeric(sa_times > 1)
sa_solve <- function(m, ...)
  solveODE(m, sa_times, sa_p, abstol = 1e-26, reltol = 1e-10, roottol = 1e-12, ...)

sa_pars <- c("x", "c", "ts")
sa_dy <- cbind(x = 2 * sa_on, c = 2 * sa_on, ts = 0)
sa_dz <- cbind(x = 0, c = 0, ts = -sa_on)
sa_d2y <- array(0, c(4, 3, 3), list(NULL, sa_pars, sa_pars))
sa_d2y[, "x", "c"] <- sa_d2y[, "c", "x"] <- -4 * sa_on
sa_d2y[, "c", "c"] <- -8 * sa_on

test_that("two switches at one time take a saltation each", {
  out <- sa_solve(sa_mods$d2)
  expect_equal(unname(out$variable[, "y"]), pmax(sa_times - 1, 0), tolerance = 1e-9)
  expect_equal(unname(out$variable[, "z"]), pmax(sa_times - 1, 0), tolerance = 1e-9)
  expect_equal(unname(out$tangent[, "y", sa_pars]), unname(sa_dy), tolerance = 1e-8)
  expect_equal(unname(out$tangent[, "z", sa_pars]), unname(sa_dz), tolerance = 1e-8)
  expect_equal(unname(out$hessian[, "y", sa_pars, sa_pars]), unname(sa_d2y),
               tolerance = 1e-7)
  expect_equal(max(abs(out$hessian[, "z", sa_pars, sa_pars])), 0, tolerance = 1e-8)

  W <- array(0, c(4, 5, 1), list(NULL, names(sa_rhs), NULL))
  W[, "y", 1] <- c(0.3, -1.2, 0.8, 1.1)
  W[, "z", 1] <- c(0.5, 0.7, -0.9, 0.4)
  fr <- sa_solve(sa_mods$fr, cotangent = W)
  expect_equal(unname(fr$cotangent[sa_pars, 1]),
               unname(colSums(sa_dy * W[, "y", 1] + sa_dz * W[, "z", 1])), tolerance = 1e-8)
  expect_equal(unname(fr$curvature[sa_pars, sa_pars, 1]),
               unname(apply(sa_d2y * W[, "y", 1], c(2, 3), sum)), tolerance = 1e-7)
})

test_that("a slow crossing takes its saltation", {
  out <- sa_solve(sa_mods$d2)
  expect_equal(unname(out$variable[, "v"]), pmax(sa_times - 1, 0), tolerance = 1e-9)
  expect_equal(unname(out$tangent[, "v", "u"]) * sa_p[["e"]], sa_on, tolerance = 1e-8)
  skip_if_not(has_cvode, "CVODE not available")
  cv <- sa_solve(sa_mods$cv)
  expect_equal(unname(cv$tangent[, "y", sa_pars]), unname(sa_dy), tolerance = 1e-7)
  expect_equal(unname(cv$tangent[, "z", sa_pars]), unname(sa_dz), tolerance = 1e-7)
  expect_equal(unname(cv$tangent[, "v", "u"]) * sa_p[["e"]], sa_on, tolerance = 1e-7)
})
