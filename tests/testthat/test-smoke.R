## Runs wherever code can be generated, CRAN included: one model through the
## whole pipeline, against the closed form.
test_that("a decay model and an observation function compile and solve", {
  skip_if_not(codegenAvailable())

  m <- cppODE(c(A = "-k*A"), modelname = "smoke_ode", compile = FALSE)
  g <- cppFUN(c(y = "s*A"), parameters = "s", modelname = "smoke_fun",
              compile = FALSE)
  compile(m, g, output = "smoke_all")

  times <- c(0, 1, 2)
  res <- solveODE(m, times, c(A = 2, k = 0.5))
  expect_equal(res$variable[, "A"], 2 * exp(-0.5 * times), tolerance = 1e-5)
  expect_equal(unname(res$tangent[, "A", "k"]), -2 * times * exp(-0.5 * times),
               tolerance = 1e-5)

  ev <- g$evaluate(A = 3, s = 2)
  expect_equal(unname(ev$y[1, "y"]), 6)
  expect_equal(unname(ev$tangent[1, "y", ]), c(2, 3))
})
