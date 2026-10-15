# Expressions read as R reads them, checked against R evaluating the same
# strings, or rejected with a message that names the construct.

skip_on_cran()

library(cppDE)

er_eqns <- c(and_or  = "(x > a & x < b) + 2*(x < a | x > b)",
             ifelse  = "ifelse(x > a, a*x, b)",
             trunc   = "trunc(x)",
             round   = "round(x)",
             rnd_par = "round(a*x)")
er_mods <- list(
  fun = cppFUN(er_eqns, variables = "x", parameters = c("a", "b"),
               modelname = "er_fun", convenient = TRUE),
  kinks = cppFUN(c(sg = "sign(a)", hv = "Heaviside(a)", ab = "abs(a)*a",
                   mx = "max(a, b)", mn = "min(a, b)", mx0 = "max(0, a)",
                   mn0 = "min(0, a)"),
                 variables = NULL, parameters = c("a", "b"), deriv2 = TRUE,
                 derivMode = c("forward", "forward-reverse"), modelname = "er_kinks",
                 convenient = FALSE),
  # A pulse written with the single `&` of R.
  pulse = cppODE(c(x = "(time > t1 & time < t2)"), modelname = "er_pulse",
                 deriv = TRUE, compile = FALSE)
)
do.call(compile, c(unname(er_mods), list(output = "test_expressions_r",
                                         cores = test_cores())))

test_that("&, |, ifelse(), trunc() and round() give what R gives", {
  xs <- c(-2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5, -1.7, 1.7, 0.2, 2.2,
          0.49999999999999994)
  a <- 0.5
  b <- 2
  out <- er_mods$fun$func(x = xs, a = a, b = b)
  for (nm in names(er_eqns)) {
    ref <- eval(parse(text = er_eqns[[nm]]), list(x = xs, a = a, b = b))
    expect_equal(unname(out[, nm]), as.numeric(ref), info = nm)
  }
})

test_that("a pulse written with & ends", {
  times <- c(0, 0.5, 1.5, 2.5, 3)
  p <- c(x = 0, t1 = 1, t2 = 2)
  out <- solveODE(er_mods$pulse, times, p, abstol = 1e-10, reltol = 1e-10)
  inside <- times > p[["t1"]] & times < p[["t2"]]
  after <- times > p[["t2"]]
  expect_equal(unname(out$variable[, "x"]), pmin(pmax(times - 1, 0), 1),
               tolerance = 1e-9)
  expect_equal(unname(out$tangent[, "x", "t1"]), -as.numeric(inside | after),
               tolerance = 1e-8)
  expect_equal(unname(out$tangent[, "x", "t2"]), as.numeric(after),
               tolerance = 1e-8)
})

test_that("syntax R reads otherwise, or not at all, is rejected by name", {
  bad <- c("1 < x && x < 2 < a" = "chained comparison",
           "a < x < b" = "chained comparison",
           "x %% a" = "'%%' is not supported; write a - floor",
           "x %/% a" = "'%/%' is not supported; write floor",
           "x % a" = "'%' is not supported",
           "x // a" = "'//' is not supported",
           "x<-a" = "'<-' is an assignment in R",
           "~x" = "'~' is not supported",
           "a if x > 1 else b" = "Python keyword used as a symbol name",
           "foo(x)" = "unknown function 'foo'")
  for (e in names(bad))
    expect_error(cppFUN(c(y = e), variables = "x", parameters = c("a", "b"),
                        modelname = "er_bad"),
                 bad[[e]], fixed = TRUE, info = e)
})

test_that("cppFUN() differentiates kinks alike forward and reverse", {
  f <- er_mods$kinks
  no <- 7L
  W <- array(0, c(1L, no, no))
  for (o in seq_len(no)) W[1L, o, o] <- 1
  for (P in list(c(a = 0, b = 0), c(a = 0.3, b = -0.2), c(a = -0.3, b = 0.3))) {
    J <- f$jac(NULL, P)[1, , ]
    H <- f$hess(NULL, P)[1, , , ]
    r <- f$vjp(NULL, P, W, tangentP = diag(2))
    expect_equal(unname(t(r$cotangentP)), unname(J), info = P)
    for (o in seq_len(no))
      expect_equal(unname(r$curvatureP[, o, ]), unname(H[o, , ]), info = c(P, o))
  }
  # At the switch point the derivative is that of the branch the value takes:
  # a tie goes to the first argument of max() and the second of min(), with a
  # constant argument last.
  J0 <- unname(f$jac(NULL, c(a = 0, b = 0))[1, , ])
  expect_equal(J0, rbind(c(0, 0), c(0, 0), c(0, 0), c(1, 0), c(0, 1), c(1, 0),
                         c(0, 0)))
})
