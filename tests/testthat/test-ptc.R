# Steady states by pseudo-transient continuation (ptc), against closed forms.

skip_on_cran()

# -- Shared setup --------------------------------------------------------------

mods <- list(
  pd = cppFUN(c(A = "k_in - k_out * A"), variables = "A",
              parameters = c("k_in", "k_out"), modelname = "ptc_pd", compile = FALSE),
  ab = cppFUN(c(A = "-k * A + km * B", B = "k * A - km * B"), variables = c("A", "B"),
              parameters = c("k", "km"), modelname = "ptc_ab", compile = FALSE),
  bi = cppFUN(c(X = "k0 + k1 * X^2 / (K^2 + X^2) - d * X"), variables = "X",
              parameters = c("k0", "k1", "K", "d"), modelname = "ptc_bi", compile = FALSE),
  sq = cppFUN(c(x = "x^2 - a"), variables = "x", parameters = "a",
              modelname = "ptc_sq", compile = FALSE),
  ti = cppFUN(c(F = "kf * (1 - F)", S = "ks * F - kd * S"), variables = c("F", "S"),
              parameters = c("kf", "ks", "kd"), modelname = "ptc_ti", compile = FALSE),
  sc = cppFUN(c(small = "ks - d * small", big = "kb - d * big", fix = "0"),
              variables = c("small", "big", "fix"), parameters = c("ks", "kb", "d"),
              modelname = "ptc_sc", compile = FALSE),
  rt = cppFUN(c(slow = "ks - kslow * slow", fast = "kf - kfast * fast"),
              variables = c("slow", "fast"), parameters = c("ks", "kslow", "kf", "kfast"),
              modelname = "ptc_rt", compile = FALSE),
  fo = cppFUN(c(x = "a * (x - 1) - (y - 1)", y = "(x - 1) + a * (y - 1)"), variables = c("x", "y"),
              parameters = "a", modelname = "ptc_fo", compile = FALSE)
)
do.call(compile, c(unname(mods), list(output = "test_ptc", cores = test_cores())))

# -- Flow ---------------------------------------------------------------------------

test_that("ptc reaches production-degradation steady state", {
  r <- ptc(mods$pd, x = c(A = 1), parms = c(k_in = 2, k_out = 0.5))
  expect_true(r$converged)
  expect_equal(r$x[["A"]], 4, tolerance = 1e-10)
})

test_that("ptc keeps a conserved total", {
  r <- ptc(mods$ab, x = c(A = 1, B = 1), parms = c(k = 2, km = 0.5),
           C = rbind(tot = c(A = 1, B = 1)), total = 5)
  expect_true(r$converged)
  expect_equal(unname(r$x[c("A", "B")]), c(1, 4), tolerance = 1e-10)
})

test_that("ptc follows the flow into the basin it starts in", {
  p  <- c(k0 = 0.02, k1 = 1, K = 1, d = 0.4)
  f  <- function(x) p[["k0"]] + p[["k1"]] * x^2 / (p[["K"]]^2 + x^2) - p[["d"]] * x
  lo <- uniroot(f, c(0.01, 0.1), tol = 1e-14)$root
  mi <- uniroot(f, c(0.1, 1), tol = 1e-14)$root
  hi <- uniroot(f, c(1, 3), tol = 1e-14)$root
  expect_equal(ptc(mods$bi, c(X = 0.8 * mi), p)$x[["X"]], lo, tolerance = 1e-8)
  expect_equal(ptc(mods$bi, c(X = 1.2 * mi), p)$x[["X"]], hi, tolerance = 1e-8)
})

test_that("ptc leaves an unstable root to the side it came from", {
  p  <- c(k0 = 0.02, k1 = 1, K = 1, d = 0.4)
  f  <- function(x) p[["k0"]] + p[["k1"]] * x^2 / (p[["K"]]^2 + x^2) - p[["d"]] * x
  lo <- uniroot(f, c(0.01, 0.1), tol = 1e-14)$root
  mi <- uniroot(f, c(0.1, 1), tol = 1e-14)$root
  hi <- uniroot(f, c(1, 3), tol = 1e-14)$root
  for (s in c(0.99, 0.999)) {
    r <- ptc(mods$bi, c(X = s * mi), p)
    expect_true(r$converged)
    expect_equal(r$x[["X"]], lo, tolerance = 1e-8)
  }
  for (s in c(1.001, 1.01, 1.1)) expect_equal(ptc(mods$bi, c(X = s * mi), p)$x[["X"]], hi, tolerance = 1e-8)
  r <- ptc(mods$bi, c(X = 1.01 * mi), p, controls = list(stability = FALSE))
  expect_equal(r$x[["X"]], mi, tolerance = 1e-8)
})

test_that("ptc accepts an unstable root whose growth is below stabilityTol", {
  p <- c(a = 1e-3)
  expect_false(ptc(mods$fo, c(x = 1.01, y = 0.99), p)$converged)
  r <- ptc(mods$fo, c(x = 1.01, y = 0.99), p, controls = list(stabilityTol = 2e-3))
  expect_true(r$converged)
  expect_equal(unname(r$x[c("x", "y")]), c(1, 1), tolerance = 1e-8)
})

test_that("ptc brings a tiny state to the relative tolerance from any start", {
  # a tiny state meets the residual test long before it is accurate
  p <- c(kf = 4, ks = 3e-14, kd = 0.03)
  for (x0 in c(1, 1e-3)) {
    r <- ptc(mods$ti, x = c(F = x0, S = x0), parms = p)
    expect_true(r$converged)
    expect_lt(abs(r$x[["S"]] / 1e-12 - 1), 1e-9)
  }
})

test_that("ptc resolves states fifteen orders apart and holds unsolved ones", {
  r <- ptc(mods$sc, x = c(small = 1, big = 1, fix = 3), parms = c(ks = 1e-10, kb = 1e5, d = 0.1),
           solve = c("small", "big"))
  expect_true(r$converged)
  expect_equal(unname(r$x[c("small", "big")]), c(1e-9, 1e6), tolerance = 1e-10)
  expect_identical(r$x[["fix"]], 3)
})

test_that("ptc resolves rates sixteen orders apart", {
  r <- ptc(mods$rt, x = c(slow = 3, fast = 3),
           parms = c(ks = 1e-16, kslow = 1e-16, kf = 1, kfast = 1))
  expect_true(r$converged)
  expect_equal(unname(r$x), c(1, 1), tolerance = 1e-10)
})

# -- Plain equations ---------------------------------------------------------------

test_that("ptc solves plain equations with flow = FALSE", {
  # The root is unstable read as a flow and found read as an equation.
  r <- ptc(mods$sq, x = c(x = 5), parms = c(a = 4), flow = FALSE)
  expect_true(r$converged)
  expect_equal(r$x[["x"]], 2, tolerance = 1e-10)
  r <- ptc(mods$sq, x = c(x = 1), parms = c(a = 4), flow = FALSE, positive = FALSE)
  expect_equal(r$x[["x"]], 2, tolerance = 1e-10)
})

test_that("ptc reports what it could not do", {
  r <- ptc(mods$sq, x = c(x = 1), parms = c(a = -1), flow = FALSE, positive = FALSE,
           controls = list(maxit = 20))
  expect_false(r$converged)
  expect_match(r$message, "no convergence in 20 iterations")
  expect_error(ptc(mods$pd, x = c(A = 1), parms = c(k_in = 1)), "misses k_out")
  expect_error(ptc(mods$pd, x = c(A = 1), parms = c(k_in = 1, k_out = 1),
                   controls = list(tol = 1)), "unknown controls: tol")
})

test_that("ptc takes reltol and abstol only", {
  p <- c(k_in = 2, k_out = 0.5)
  r <- ptc(mods$pd, x = c(A = 1), parms = p, controls = list(reltol = 1e-8, abstol = 1e-12))
  expect_equal(r$x[["A"]], 4, tolerance = 1e-8)
  expect_error(ptc(mods$pd, x = c(A = 1), parms = p, controls = list(rtol = 1e-8)),
               "unknown controls: rtol")
})
