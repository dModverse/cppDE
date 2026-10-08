# Switches on a state: a comparison, Heaviside() or sign() whose argument
# reads a state is held as a mode and switched at the located root of its
# argument, so the sensitivities take the jump of f there.

skip_on_cran()

library(cppDE)

# f jumps from k to 0 where x reaches 1.
sw_rhs <- c(x = "piecewise(k, x < 1, 0)")
# A relay: the force flips sign whenever x crosses zero, in both directions.
relay <- c(x = "y", y = "piecewise(-k, x > 0, k) - c*y")
# f is continuous across x = 1 and its derivative in x is not.
kink <- c(x = "piecewise(-a*x, x > 1, -a - b*(x - 1))")
# Heaviside() and sign() of a state, both switching where x passes 1.
hs <- c(x = "c", y = "k*sign(x - 1)", z = "Heaviside(1 - x)")
# Jumps that move the state across the switch and back below it.
sw_events <- data.frame(var = "x", time = c(1, 3), value = c(1, 0.5),
                        method = c("add", "replace"))

sw <- function(rhs, nm, ...) cppODE(rhs, modelname = nm, compile = FALSE, ...)
mods <- list(
  ss_d2        = sw(sw_rhs, "ss_d2", deriv2 = TRUE),
  ss_adams     = sw(sw_rhs, "ss_adams", method = "adams"),
  ss_rb4       = sw(sw_rhs, "ss_rb4", method = "rb4"),
  ss_tsit5     = sw(sw_rhs, "ss_tsit5", method = "tsit5"),
  ss_rb4_grid  = sw(sw_rhs, "ss_rb4_grid", method = "rb4", useDenseOutput = FALSE),
  ss_rev_bdf   = sw(sw_rhs, "ss_rev_bdf", derivMode = "reverse"),
  ss_rev_adams = sw(sw_rhs, "ss_rev_adams", method = "adams", derivMode = "reverse"),
  ss_rev_rb4   = sw(sw_rhs, "ss_rev_rb4", method = "rb4", derivMode = "reverse"),
  ss_rev_tsit5 = sw(sw_rhs, "ss_rev_tsit5", method = "tsit5", derivMode = "reverse"),
  ss_fr        = sw(sw_rhs, "ss_fr", derivMode = "forward-reverse"),
  rl_plain     = sw(relay, "rl_plain", deriv = FALSE),
  rl_d2        = sw(relay, "rl_d2", deriv2 = TRUE),
  rl_rb4       = sw(relay, "rl_rb4", method = "rb4"),
  rl_rev       = sw(relay, "rl_rev", derivMode = "reverse"),
  rl_rev_rb4   = sw(relay, "rl_rev_rb4", method = "rb4", derivMode = "reverse"),
  rl_fr        = sw(relay, "rl_fr", derivMode = "forward-reverse"),
  kink_d2      = sw(kink, "kink_d2", deriv2 = TRUE),
  hs           = sw(hs, "hs_fwd"),
  ev           = sw(sw_rhs, "ss_ev", events = sw_events),
  slide        = sw(c(x = "piecewise(1, x < 1, -1)"), "ss_slide", deriv = FALSE)
)
has_cvode <- isTRUE(cvodeConfig$available)
if (has_cvode) {
  mods$ss_cv      <- cvode(sw_rhs, modelname = "ss_cv", deriv = TRUE, compile = FALSE)
  mods$ss_cv_rev  <- cvode(sw_rhs, modelname = "ss_cv_rev", derivMode = "reverse",
                           compile = FALSE)
  mods$rl_cv      <- cvode(relay, modelname = "rl_cv", deriv = TRUE, compile = FALSE)
  mods$rl_cv_rev  <- cvode(relay, modelname = "rl_cv_rev", derivMode = "reverse",
                           compile = FALSE)
}
do.call(compile, c(unname(mods), list(output = "test_state_switch",
                                      cores = test_cores())))

tight <- function(m, times, p, ...)
  solveODE(m, times, p, abstol = 1e-10, reltol = 1e-10, roottol = 1e-10, ...)

# -- The switch at x = 1 ------------------------------------------------------

sw_times <- seq(0, 4, by = 0.5)
sw_p <- c(x = 0.2, k = 0.5)
sw_tau <- (1 - sw_p[["x"]]) / sw_p[["k"]]
sw_before <- sw_times < sw_tau
sw_x <- pmin(sw_p[["x"]] + sw_p[["k"]] * sw_times, 1)
sw_dx <- cbind(x = as.numeric(sw_before), k = sw_times * sw_before)

test_that("the switch is located and the tangents take its jump", {
  for (nm in c("ss_d2", "ss_adams", "ss_rb4", "ss_tsit5", "ss_rb4_grid")) {
    out <- tight(mods[[nm]], sw_times, sw_p)
    expect_equal(nrow(out$variable), length(sw_times), info = nm)
    expect_equal(unname(out$variable[, "x"]), sw_x, tolerance = 1e-9, info = nm)
    expect_equal(unname(out$tangent[, "x", c("x", "k")]), unname(sw_dx),
                 tolerance = 1e-8, info = nm)
  }
})

test_that("the second-order tangents are those of the switched solution", {
  out <- tight(mods$ss_d2, sw_times, sw_p)
  expect_equal(max(abs(out$hessian)), 0, tolerance = 1e-9)
})

test_that("the reverse gradient takes the jump on every method", {
  set.seed(3)
  W <- array(rnorm(length(sw_times)), c(length(sw_times), 1, 1))
  ref <- colSums(sw_dx * W[, 1, 1])
  for (nm in c("ss_rev_bdf", "ss_rev_adams", "ss_rev_rb4", "ss_rev_tsit5", "ss_fr")) {
    out <- tight(mods[[nm]], sw_times, sw_p, cotangent = W)
    expect_equal(unname(out$variable[, "x"]), sw_x, tolerance = 1e-9, info = nm)
    expect_equal(unname(out$cotangent[c("x", "k"), 1]), unname(ref),
                 tolerance = 1e-8, info = nm)
  }
  out <- tight(mods$ss_fr, sw_times, sw_p, cotangent = W)
  expect_equal(max(abs(out$curvature)), 0, tolerance = 1e-9)
})

test_that("cvode() locates the switch, forward and reverse", {
  skip_if_not(has_cvode, "CVODE not available")
  out <- tight(mods$ss_cv, sw_times, sw_p)
  expect_equal(nrow(out$variable), length(sw_times))
  expect_equal(unname(out$variable[, "x"]), sw_x, tolerance = 1e-9)
  expect_equal(unname(out$tangent[, "x", c("x", "k")]), unname(sw_dx),
               tolerance = 1e-8)

  set.seed(3)
  W <- array(rnorm(length(sw_times)), c(length(sw_times), 1, 1))
  rev <- tight(mods$ss_cv_rev, sw_times, sw_p, cotangent = W)
  expect_equal(unname(rev$cotangent[c("x", "k"), 1]),
               unname(colSums(sw_dx * W[, 1, 1])), tolerance = 1e-7)
})

# -- Re-entry -----------------------------------------------------------------

rl_times <- seq(0, 10, by = 0.25)
rl_p <- c(x = 1, y = 0, k = 1, c = 0.1)

# Central differences of tight solves without derivatives.
rl_fd <- function(h = 1e-6) {
  run <- function(p) solveODE(mods$rl_plain, rl_times, p, abstol = 1e-12,
                              reltol = 1e-12, roottol = 1e-12)$variable
  out <- lapply(names(rl_p), function(n) {
    up <- dn <- rl_p
    up[n] <- up[n] + h
    dn[n] <- dn[n] - h
    (run(up) - run(dn)) / (2 * h)
  })
  array(unlist(out), c(length(rl_times), 2, length(rl_p)))
}

test_that("a relay that switches back and forth matches finite differences", {
  fd <- rl_fd()
  d2 <- tight(mods$rl_d2, rl_times, rl_p)
  # Several crossings in both directions.
  expect_gt(sum(diff(sign(d2$variable[, "x"])) != 0), 2)
  expect_equal(unname(d2$tangent[, , names(rl_p)]), fd, tolerance = 1e-6)
  rb4 <- tight(mods$rl_rb4, rl_times, rl_p)
  expect_equal(unname(rb4$tangent[, , names(rl_p)]), fd, tolerance = 1e-6)
})

test_that("the relay's second-order tangents match differences of the first", {
  h <- 1e-5
  tan <- function(n, s) {
    p <- rl_p
    p[n] <- p[n] + s * h
    tight(mods$rl_d2, rl_times, p)$tangent[, , names(rl_p)]
  }
  H <- tight(mods$rl_d2, rl_times, rl_p)$hessian[, , names(rl_p), names(rl_p)]
  for (n in names(rl_p))
    expect_equal(unname(H[, , , n]), unname((tan(n, 1) - tan(n, -1)) / (2 * h)),
                 tolerance = 1e-5, info = n)
})

test_that("the relay's reverse gradient and curvature match the forward mode", {
  d2 <- tight(mods$rl_d2, rl_times, rl_p)
  set.seed(5)
  W <- array(rnorm(length(d2$variable)), c(dim(d2$variable), 1))
  grad <- apply(d2$tangent * as.vector(W[, , 1]), 3, sum)[names(rl_p)]
  curv <- apply(d2$hessian * as.vector(W[, , 1]), c(3, 4), sum)
  for (nm in c("rl_rev", "rl_rev_rb4", "rl_fr")) {
    out <- tight(mods[[nm]], rl_times, rl_p, cotangent = W)
    expect_equal(unname(out$cotangent[names(rl_p), 1]), unname(grad),
                 tolerance = 1e-7, info = nm)
  }
  fr <- tight(mods$rl_fr, rl_times, rl_p, cotangent = W)
  expect_equal(unname(fr$curvature[names(rl_p), names(rl_p), 1]),
               unname(curv[names(rl_p), names(rl_p)]), tolerance = 1e-6)
})

test_that("cvode() follows the relay forward and backward", {
  skip_if_not(has_cvode, "CVODE not available")
  fd <- rl_fd()
  fwd <- tight(mods$rl_cv, rl_times, rl_p)
  expect_equal(nrow(fwd$variable), length(rl_times))
  expect_equal(unname(fwd$tangent[, , names(rl_p)]), fd, tolerance = 1e-6)

  set.seed(5)
  W <- array(rnorm(length(fwd$variable)), c(dim(fwd$variable), 1))
  rev <- tight(mods$rl_cv_rev, rl_times, rl_p, cotangent = W)
  grad <- apply(fwd$tangent * as.vector(W[, , 1]), 3, sum)[names(rl_p)]
  expect_equal(unname(rev$cotangent[names(rl_p), 1]), unname(grad),
               tolerance = 1e-6)
})

# -- Continuous and discontinuous f -------------------------------------------

test_that("a kink, continuous in f, has the closed-form first and second order", {
  times <- seq(0, 4, by = 0.5)
  p <- c(x = 2, a = 0.5, b = 2)
  closed <- function(q) {
    tau <- log(q[["x"]]) / q[["a"]]
    r <- q[["a"]] / q[["b"]]
    ifelse(times < tau, q[["x"]] * exp(-q[["a"]] * times),
           1 - r + r * exp(-q[["b"]] * (times - tau)))
  }
  d1 <- function(q, n, h = 1e-6) {
    up <- dn <- q
    up[n] <- up[n] + h
    dn[n] <- dn[n] - h
    (closed(up) - closed(dn)) / (2 * h)
  }
  out <- tight(mods$kink_d2, times, p)
  expect_equal(unname(out$variable[, "x"]), closed(p), tolerance = 1e-9)
  for (n in names(p)) {
    expect_equal(unname(out$tangent[, "x", n]), d1(p, n), tolerance = 1e-7, info = n)
    for (m in names(p)) {
      h <- 1e-4
      up <- dn <- p
      up[m] <- up[m] + h
      dn[m] <- dn[m] - h
      ref <- (d1(up, n) - d1(dn, n)) / (2 * h)
      expect_equal(unname(out$hessian[, "x", n, m]), ref, tolerance = 1e-5,
                   info = paste(n, m))
    }
  }
})

test_that("Heaviside() and sign() of a state switch where their argument does", {
  times <- seq(0, 3, by = 0.5)
  p <- c(x = 0.4, y = 0, z = 0, c = 0.5, k = 2)
  tau <- (1 - p[["x"]]) / p[["c"]]
  before <- times < tau
  out <- tight(mods$hs, times, p)
  expect_equal(unname(out$variable[, "y"]),
               ifelse(before, -p[["k"]] * times, p[["k"]] * (times - 2 * tau)),
               tolerance = 1e-9)
  expect_equal(unname(out$variable[, "z"]), pmin(times, tau), tolerance = 1e-9)
  # tau moves with x0 and c; y and z follow it past the switch.
  dtau <- c(x = -1 / p[["c"]], c = -(1 - p[["x"]]) / p[["c"]]^2)
  for (n in names(dtau)) {
    expect_equal(unname(out$tangent[, "y", n]),
                 ifelse(before, 0, -2 * p[["k"]] * dtau[[n]]), tolerance = 1e-8, info = n)
    expect_equal(unname(out$tangent[, "z", n]),
                 ifelse(before, 0, dtau[[n]]), tolerance = 1e-8, info = n)
  }
})

test_that("a jump across the switch sets the mode on its far side", {
  times <- seq(0, 5, by = 0.5)
  p <- c(x = 0.2, k = 0.4)
  out <- tight(mods$ev, times, p)
  # x rises to 0.6 at t = 1, jumps to 1.6 and holds; at t = 3 it drops to 0.5,
  # rises again and stops at 1 from t = 3 + 0.5 / k on.
  t2 <- 3 + 0.5 / p[["k"]]
  x <- ifelse(times < 1, 0.2 + 0.4 * times,
              ifelse(times < 3, 1.6, pmin(0.5 + 0.4 * (times - 3), 1)))
  dk <- ifelse(times < 1, times, ifelse(times < 3, 1,
                                        ifelse(times < t2, times - 3, 0)))
  expect_equal(unname(out$variable[, "x"]), x, tolerance = 1e-9)
  expect_equal(unname(out$tangent[, "x", "k"]), dk, tolerance = 1e-8)
  expect_equal(unname(out$tangent[, "x", "x"]), as.numeric(times < 3), tolerance = 1e-8)
})

# -- Batches ------------------------------------------------------------------

test_that("a batch switches every condition where its own root lies", {
  conds <- list(early = list(parms = c(x = 0.2, k = 0.5)),
                late = list(parms = c(x = 0.5, k = 0.2)),
                above = list(parms = c(x = 1.2, k = 0.5)))
  bat <- solveODEBatch(mods$ss_d2, conds, times = sw_times, abstol = 1e-10,
                       reltol = 1e-10, roottol = 1e-10, cores = 2)
  for (cn in names(conds)) {
    ser <- tight(mods$ss_d2, sw_times, conds[[cn]]$parms)
    expect_identical(bat[[cn]]$variable, ser$variable, info = cn)
    expect_identical(bat[[cn]]$tangent, ser$tangent, info = cn)
    p <- conds[[cn]]$parms
    tau <- (1 - p[["x"]]) / p[["k"]]
    x <- if (p[["x"]] >= 1) rep(p[["x"]], length(sw_times))
         else pmin(p[["x"]] + p[["k"]] * sw_times, 1)
    expect_equal(unname(bat[[cn]]$variable[, "x"]), x, tolerance = 1e-9, info = cn)
    expect_equal(unname(bat[[cn]]$tangent[, "x", "k"]),
                 sw_times * (sw_times < tau), tolerance = 1e-8, info = cn)
  }
  skip_if_not(has_cvode, "CVODE not available")
  bat <- solveODEBatch(mods$ss_cv, conds, times = sw_times, abstol = 1e-10,
                       reltol = 1e-10, roottol = 1e-10, cores = 2)
  for (cn in names(conds)) {
    ser <- tight(mods$ss_cv, sw_times, conds[[cn]]$parms)
    expect_identical(bat[[cn]]$variable, ser$variable, info = cn)
    expect_identical(bat[[cn]]$tangent, ser$tangent, info = cn)
  }
})

test_that("a solution that slides along the surface stops with a message", {
  # Both branches point at x = 1, so the mode would have to flip without end.
  expect_error(tight(mods$slide, sw_times, c(x = 0.5)), "slides along")
})
