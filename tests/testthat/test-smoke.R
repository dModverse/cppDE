## Runs wherever code can be generated, CRAN included: one model through the
## whole pipeline, against the closed form.
test_that("a decay model and an observation function compile and solve", {
  skip_if_not(codegenAvailable())

  m <- cppODE(c(A = "-k*A"), modelname = "smoke_ode", compile = FALSE)
  g <- cppFUN(c(y = "s*A"), parameters = "s", modelname = "smoke_fun",
              compile = FALSE)
  # compile() hands the build its own flags and restores whatever was set before.
  old <- Sys.getenv("PKG_CFLAGS", unset = NA)
  Sys.setenv(PKG_CFLAGS = "-DCPPDE_SMOKE")
  on.exit(if (is.na(old)) Sys.unsetenv("PKG_CFLAGS") else Sys.setenv(PKG_CFLAGS = old))
  compile(m, g, output = "smoke_all")
  expect_identical(Sys.getenv("PKG_CFLAGS"), "-DCPPDE_SMOKE")

  times <- c(0, 1, 2)
  res <- solveODE(m, times, c(A = 2, k = 0.5))
  expect_equal(res$variable[, "A"], 2 * exp(-0.5 * times), tolerance = 1e-5)
  expect_equal(unname(res$tangent[, "A", "k"]), -2 * times * exp(-0.5 * times),
               tolerance = 1e-5)

  ev <- g$evaluate(A = 3, s = 2)
  expect_equal(unname(ev$y[1, "y"]), 6)
  expect_equal(unname(ev$tangent[1, "y", ]), c(2, 3))
})
