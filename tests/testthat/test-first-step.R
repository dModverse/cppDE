# The first step of a solve that starts late under a tight tolerance. A
# sensitivity starts at zero with a rate far above atol, which pushes the
# bounds of the multistep methods' estimate across each other; the step they
# took then fell below a tick of t at t0 = 1e4, t + h rounded back to t, and
# the stepper could not move. The estimate lives in
# inst/include/cppde/cppde_utils.hpp.

skip_on_cran()

methods <- c("bdf", "adams", "rb4", "tsit5")
models <- lapply(setNames(nm = methods), function(m)
  cppODE(c(x = "-k*x"), method = m, modelname = paste0("fstep_", m),
         deriv = TRUE, includeTimeZero = FALSE, compile = FALSE))
do.call(compile, c(unname(models), list(output = "test_first_step", cores = 1)))

parms <- c(x = 1e3, k = 0.1)
exact <- function(t) 1e3 * exp(-0.1 * t)

test_that("a solve starting late under a tight tolerance takes its first step", {
  for (m in methods) for (t0 in c(0, 280, 1e4)) {
    out <- solveODE(models[[m]], t0 + c(0, 0.5, 70), parms,
                    abstol = 1e-12, reltol = 1e-12, onFailure = "stop")
    expect_equal(unname(out$variable[, "x"]), exact(c(0, 0.5, 70)),
                 tolerance = 1e-8, info = paste(m, t0))
    # the sensitivity of x in k, which starts at zero and carries the rate
    expect_equal(unname(out$tangent[3, "x", "k"]), -70 * exact(70),
                 tolerance = 1e-6, info = paste(m, t0))
  }
})

test_that("the first step never falls below four ticks of the time variable", {
  for (m in c("bdf", "adams")) for (t0 in c(1e4, 1e6)) {
    out <- solveODE(models[[m]], t0 + c(0, 1), parms,
                    abstol = 1e-14, reltol = 1e-14, onFailure = "stop")
    expect_equal(unname(out$variable[2, "x"]), exact(1), tolerance = 1e-8,
                 info = paste(m, t0))
  }
})
