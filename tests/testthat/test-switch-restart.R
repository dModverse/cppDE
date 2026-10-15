# A switch whose root sits at zero where the integration starts or restarts:
# at t0, at the time of a fixed event, and after its own crossing. The mode
# follows the side the solution leaves to, the sensitivities take the jump of f
# once, and a solution that both branches push onto the surface slides.

skip_on_cran()

library(cppDE)

rs_rhs <- c(x = "piecewise(1, time > ts, 0)")
rs_events <- data.frame(var = "x", time = 1, value = 0.5, method = "add")
rs_slide <- c(x = "piecewise(1, x < 1, -1)")

rs <- function(nm, ...)
  cppODE(rs_rhs, modelname = nm, events = rs_events, compile = FALSE, ...)
rs_mods <- list(
  d2       = rs("rs_d2", deriv2 = TRUE),
  rb4_grid = rs("rs_rb4_grid", method = "rb4", useDenseOutput = FALSE),
  rev      = rs("rs_rev", derivMode = "reverse"),
  fr       = rs("rs_fr", derivMode = "forward-reverse")
)
do.call(compile, c(unname(rs_mods), list(output = "test_switch_restart",
                                         cores = test_cores())))
has_cvode <- isTRUE(cvodeConfig$available)
if (has_cvode) {
  rs_cv <- list(
    fwd   = cvode(rs_rhs, modelname = "rs_cv", events = rs_events, deriv = TRUE,
                  compile = FALSE),
    rev   = cvode(rs_rhs, modelname = "rs_cv_rev", events = rs_events,
                  derivMode = "reverse", compile = FALSE),
    slide = cvode(rs_slide, modelname = "rs_cv_slide", deriv = FALSE,
                  compile = FALSE)
  )
  do.call(compile, c(unname(rs_cv), list(output = "test_switch_restart_cv",
                                         cores = test_cores())))
}

tight <- function(m, times, p, ...)
  solveODE(m, times, p, abstol = 1e-10, reltol = 1e-10, roottol = 1e-10, ...)

# The row at an event time shows the state after the event and the switch; the
# t0 row does not depend on the switch time.
rs_times <- c(0, 0.5, 1, 1.5, 2)
rs_closed <- function(ts) {
  on <- rs_times > ts | (rs_times == ts & ts > 0)
  list(x = pmax(rs_times - ts, 0) + 0.5 * (rs_times >= 1),
       dx = cbind(x = 1, ts = -as.numeric(on)))
}
set.seed(11)
rs_W <- array(rnorm(length(rs_times)), c(length(rs_times), 1, 1))

rs_ts <- c(0, 1, 0.75)

expect_closed_form <- function(out, ts, info) {
  ref <- rs_closed(ts)
  expect_equal(unname(out$variable[, "x"]), ref$x, tolerance = 1e-9, info = info)
  expect_equal(unname(out$tangent[, "x", c("x", "ts")]), unname(ref$dx),
               tolerance = 1e-8, info = info)
}
expect_closed_gradient <- function(out, ts, info) {
  ref <- colSums(rs_closed(ts)$dx * rs_W[, 1, 1])
  expect_equal(unname(out$cotangent[c("x", "ts"), 1]), unname(ref),
               tolerance = 1e-8, info = info)
}

test_that("a switch at t0 or at an event time takes the jump of f once", {
  for (ts in rs_ts) for (nm in c("d2", "rb4_grid")) {
    out <- tight(rs_mods[[nm]], rs_times, c(x = 0, ts = ts))
    expect_closed_form(out, ts, paste(nm, ts))
  }
})

test_that("its second-order tangents vanish as those of the closed form", {
  for (ts in rs_ts) {
    out <- tight(rs_mods$d2, rs_times, c(x = 0, ts = ts))
    expect_equal(max(abs(out$hessian)), 0, tolerance = 1e-9, info = ts)
  }
})

test_that("the reverse gradient takes the jump once, the t0 row before it", {
  for (ts in rs_ts) for (nm in c("rev", "fr")) {
    out <- tight(rs_mods[[nm]], rs_times, c(x = 0, ts = ts), cotangent = rs_W)
    expect_closed_gradient(out, ts, paste(nm, ts))
  }
  out <- tight(rs_mods$fr, rs_times, c(x = 0, ts = 0), cotangent = rs_W)
  expect_equal(max(abs(out$curvature)), 0, tolerance = 1e-9)
})

test_that("cvode() switches at t0 and at an event time, forward and reverse", {
  skip_if_not(has_cvode, "CVODE not available")
  for (ts in rs_ts) {
    out <- tight(rs_cv$fwd, rs_times, c(x = 0, ts = ts))
    expect_closed_form(out, ts, ts)
    out <- tight(rs_cv$rev, rs_times, c(x = 0, ts = ts), cotangent = rs_W)
    expect_closed_gradient(out, ts, ts)
  }
})

test_that("cvode() holds a solution that slides along the surface", {
  skip_if_not(has_cvode, "CVODE not available")
  out <- tight(rs_cv$slide, rs_times, c(x = 0.5))
  expect_equal(unname(out$variable[, "x"]), pmin(0.5 + rs_times, 1), tolerance = 1e-9)
})
