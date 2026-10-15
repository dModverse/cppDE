# Kinks and steps of a state: abs(), min() and max() switch as held
# comparisons, and floor(), ceiling(), round() and trunc() of a state or of
# the time are rejected, since no root locates their jumps.

skip_on_cran()

library(cppDE)

kk <- function(rhs, nm, ...) cppODE(rhs, modelname = nm, compile = FALSE, ...)
kk_rhs <- c(x = "c", y = "piecewise(1, (x > 1) > 0.5, 0)",
            u1 = "a + b*max(u1 - 1, 0)",
            u2 = "a + b*(u2 - 1 - min(u2 - 1, 0))",
            u3 = "a + b*(u3 - 1 + abs(u3 - 1))/2")
kk_mods <- list(
  steps = kk(c(w = "round(k) + ceiling(k) - trunc(-k)"), "kk_steps", deriv2 = TRUE),
  d2    = kk(kk_rhs, "kk_d2", deriv2 = TRUE),
  fr    = kk(kk_rhs, "kk_fr", derivMode = "forward-reverse")
)
do.call(compile, c(unname(kk_mods), list(output = "test_switch_kinks",
                                         cores = test_cores())))
has_cvode <- isTRUE(cvodeConfig$available)
if (has_cvode) {
  kk_cv <- list(
    fwd = cvode(kk_rhs, modelname = "kk_cv", deriv = TRUE, compile = FALSE),
    rev = cvode(kk_rhs, modelname = "kk_cv_rev", derivMode = "reverse", compile = FALSE)
  )
  do.call(compile, c(unname(kk_cv), list(output = "test_switch_kinks_cv",
                                         cores = test_cores())))
}

tight <- function(m, times, p, ...)
  solveODE(m, times, p, abstol = 1e-10, reltol = 1e-10, roottol = 1e-10, ...)

test_that("steps of a parameter are constant, with derivative 0", {
  times <- c(0, 1, 2)
  out <- tight(kk_mods$steps, times, c(w = 0, k = 2.5))
  expect_equal(unname(out$variable[, "w"]), (2 + 3 + 2) * times, tolerance = 1e-12)
  expect_equal(max(abs(out$tangent[, "w", "k"])), 0)
  expect_equal(max(abs(out$hessian)), 0)
})

test_that("steps of a state or of the time are rejected", {
  for (rhs in c("1 + floor(x)", "ceiling(2*x)", "round(x)", "trunc(x - k)",
                "x - floor(x/k)*k", "floor(time/k)"))
    expect_error(kk(c(x = rhs), "kk_bad"), "of a state or of the time are not supported",
                 fixed = TRUE, info = rhs)
})

# -- A comparison of a comparison ---------------------------------------------

kk_times <- c(0, 0.5, 1.5, 2.5, 3)
kk_p <- c(x = 0.2, y = 0, c = 0.5, u1 = 0.4, u2 = 0.3, u3 = 0.2, a = 0.6, b = 0.8)
kk_pars <- c("x", "c")
kk_closed <- function() {
  x0 <- kk_p[["x"]]
  c <- kk_p[["c"]]
  tau <- (1 - x0) / c
  on <- kk_times > tau
  list(y = pmax(kk_times - tau, 0),
       dy = cbind(x = on / c, c = on * (1 - x0) / c^2),
       d2y = array(c(0 * on, -on / c^2, -on / c^2, -2 * on * (1 - x0) / c^3),
                   c(length(kk_times), 2, 2)))
}
kk_W <- array(c(0.3, -1.2, 0.8, 1.1, -0.4), c(length(kk_times), 1, 1))
kk_seed <- function(s) {
  W <- array(0, c(length(kk_times), length(kk_rhs), 1),
             list(NULL, names(kk_rhs), NULL))
  W[, s, 1] <- kk_W[, 1, 1]
  W
}

test_that("a condition on a held comparison switches with its mode", {
  ref <- kk_closed()
  out <- tight(kk_mods$d2, kk_times, kk_p)
  expect_equal(unname(out$variable[, "y"]), ref$y, tolerance = 1e-9)
  expect_equal(unname(out$tangent[, "y", kk_pars]), unname(ref$dy), tolerance = 1e-8)
  expect_equal(unname(out$hessian[, "y", kk_pars, kk_pars]), ref$d2y, tolerance = 1e-7)

  fr <- tight(kk_mods$fr, kk_times, kk_p, cotangent = kk_seed("y"))
  expect_equal(unname(fr$cotangent[kk_pars, 1]), unname(colSums(ref$dy * kk_W[, 1, 1])),
               tolerance = 1e-8)
  expect_equal(unname(fr$curvature[kk_pars, kk_pars, 1]),
               apply(ref$d2y * kk_W[, 1, 1], c(2, 3), sum), tolerance = 1e-7)
})

test_that("cvode() switches on a comparison of a comparison", {
  skip_if_not(has_cvode, "CVODE not available")
  ref <- kk_closed()
  out <- tight(kk_cv$fwd, kk_times, kk_p)
  expect_equal(unname(out$variable[, "y"]), ref$y, tolerance = 1e-9)
  expect_equal(unname(out$tangent[, "y", kk_pars]), unname(ref$dy), tolerance = 1e-7)
  rev <- tight(kk_cv$rev, kk_times, kk_p, cotangent = kk_seed("y"))
  expect_equal(unname(rev$cotangent[kk_pars, 1]), unname(colSums(ref$dy * kk_W[, 1, 1])),
               tolerance = 1e-7)
})

# -- abs(), min() and max() of a state ----------------------------------------

kk_u <- function(q, t = kk_times) {
  tau <- (1 - q[["u"]]) / q[["a"]]
  ifelse(t < tau, q[["u"]] + q[["a"]] * t,
         1 + q[["a"]] / q[["b"]] * (exp(q[["b"]] * (t - tau)) - 1))
}
kk_q <- function(s) c(u = kk_p[[s]], a = 0.6, b = 0.8)
kk_d1 <- function(q, n, h = 1e-6) {
  up <- dn <- q
  up[n] <- up[n] + h
  dn[n] <- dn[n] - h
  (kk_u(up) - kk_u(dn)) / (2 * h)
}
kk_d2 <- function(q, n, m, h = 1e-4) {
  up <- dn <- q
  up[m] <- up[m] + h
  dn[m] <- dn[m] - h
  (kk_d1(up, n) - kk_d1(dn, n)) / (2 * h)
}

test_that("abs(), min() and max() of a state have the closed-form derivatives", {
  out <- tight(kk_mods$d2, kk_times, kk_p)
  for (s in c("u1", "u2", "u3")) {
    pars <- c(u = s, a = "a", b = "b")
    q <- kk_q(s)
    expect_equal(unname(out$variable[, s]), kk_u(q), tolerance = 1e-8, info = s)
    for (n in names(pars)) {
      expect_equal(unname(out$tangent[, s, pars[[n]]]), kk_d1(q, n),
                   tolerance = 1e-7, info = paste(s, n))
      for (m in names(pars))
        expect_equal(unname(out$hessian[, s, pars[[n]], pars[[m]]]), kk_d2(q, n, m),
                     tolerance = 1e-5, info = paste(s, n, m))
    }
  }
})

test_that("their reverse gradient and curvature match the closed form", {
  W <- kk_seed("u3")
  fr <- tight(kk_mods$fr, kk_times, kk_p, cotangent = W)
  pars <- c(u = "u3", a = "a", b = "b")
  q <- kk_q("u3")
  g <- sapply(names(pars), function(n) sum(kk_d1(q, n) * kk_W[, 1, 1]))
  H <- outer(names(pars), names(pars),
             Vectorize(function(n, m) sum(kk_d2(q, n, m) * kk_W[, 1, 1])))
  expect_equal(unname(fr$cotangent[pars, 1]), unname(g), tolerance = 1e-7)
  expect_equal(unname(fr$curvature[pars, pars, 1]), H, tolerance = 1e-5)
  skip_if_not(has_cvode, "CVODE not available")
  out <- tight(kk_cv$fwd, kk_times, kk_p)
  for (s in c("u1", "u2", "u3"))
    expect_equal(unname(out$tangent[, s, c(s, "a", "b")]),
                 unname(sapply(names(pars), function(n) kk_d1(kk_q(s), n))),
                 tolerance = 1e-6,
                 info = s)
  rev <- tight(kk_cv$rev, kk_times, kk_p, cotangent = W)
  expect_equal(unname(rev$cotangent[pars, 1]), unname(g), tolerance = 1e-6)
})
