# Forward over reverse on the features it runs on: forcings, jumps at fixed,
# parameter and root times, a clock-reading right-hand side and a sparse
# Jacobian. The oracle is forward-forward on bdf, one per model. It keeps the
# sensitivity error control and therefore its own grid, so these are tolerance
# comparisons; the exact claims on the grid are in test-ode-reverse2.R.

skip_on_cran()

eqns <- c(A = "-k1 * A + k2 * B",
          B = "k1 * A - k2 * B - k3 * B * B")
pars  <- c(A = 1.2, B = 0.4, k1 = 0.7, k2 = 0.35, k3 = 1.1)
times <- seq(0, 5, length.out = 21)
tol   <- list(abstol = 1e-10, reltol = 1e-10)

seed_for <- function(n_t, n_x = 2L, n_seed = 1L, seed = 1L) {
  set.seed(seed)
  array(rnorm(n_t * n_x * n_seed), c(n_t, n_x, n_seed))
}

steppers <- c("bdf", "adams", "rb4", "tsit5")

# One model per method, named <prefix><method>.
per_method <- function(rhs, prefix, derivMode, meths = steppers, ...)
  lapply(setNames(nm = meths), function(meth)
    cppODE(rhs, method = meth, derivMode = derivMode,
           modelname = paste0(prefix, meth), compile = FALSE, ...))

# One forward-forward oracle on bdf and the forward-reverse models it checks.
oracle_pair <- function(rhs, tag, meths = steppers, ...)
  list(ff = cppODE(rhs, derivMode = "forward-forward",
                   modelname = paste0("g2_", tag, "_ff"), compile = FALSE, ...),
       fr = per_method(rhs, paste0("g2_", tag, "_fr_"), "forward-reverse", meths, ...))

# Models nest in lists; compile() takes them flat.
flat <- function(x) if (is.list(x)) do.call(c, lapply(unname(x), flat)) else list(x)

# A forcing, and a jump whose height is a parameter. Off the output grid, so
# the jump shows up as its own pair of rows.
eqns_u <- c(A = "-k1 * A + k2 * B + u",
            B = "k1 * A - k2 * B - k3 * B * B")
ev_fixed <- data.frame(var = "A", time = 1.1, value = "d_amt", method = "add",
                       stringsAsFactors = FALSE)
# A jump whose time is a parameter.
ev_ptime <- data.frame(var = "A", time = "t_ev", value = 0.4, method = "add",
                       stringsAsFactors = FALSE)
# A fixed and a root event, every slot an expression.
ev_general <- data.frame(
  var    = c("A", "B"),
  time   = c("0.5 * t_ev + 0.4 * t_ev * t_ev", NA),
  value  = c("k1 * A + d_amt * time", "0.3 * B + 0.2 * d_amt * A * time"),
  root   = c(NA, "B * B + 0.1 * A - 0.40 - 0.02 * time"),
  method = c("add", "add"),
  stringsAsFactors = FALSE)
# A right-hand side that reads the clock.
eqns_t <- c(A = "-k1 * time * A + k2 * B", B = "k1 * A - k2 * B - k3 * B * B")
# A root event whose height reads the clock.
ev_clock <- data.frame(var = "A", time = NA, value = "d_amt * time",
                       root = "B - 0.55", method = "add", stringsAsFactors = FALSE)

pr_forcing <- oracle_pair(eqns_u, "ev", c("bdf", "tsit5"),
                          events = ev_fixed, forcings = "u")
pr_ptime   <- oracle_pair(eqns,   "pt",   "bdf", events = ev_ptime)
pr_general <- oracle_pair(eqns_t, "gen",  events = ev_general)
pr_rtol    <- oracle_pair(eqns,   "rtol", "bdf", events = ev_clock)
do.call(compile, c(flat(list(pr_forcing, pr_ptime, pr_general, pr_rtol)),
                   output = "test_ode_reverse2_events", cores = test_cores()))

# Sparse models need KLU; their test skips without it.
if (isTRUE(cppDE:::cvodeConfig$klu_available)) {
  pr_sparse <- oracle_pair(eqns, "sp", "bdf", sparse = TRUE)
  do.call(compile, c(flat(pr_sparse), output = "test_ode_reverse2_sparse",
                     cores = test_cores()))
}

# w' S contracted over times and states: the gradient the sweep returns.
contract <- function(tangent, W) {
  vapply(seq_len(dim(W)[3]),
         function(k) apply(tangent * as.vector(W[, , k]), 3, sum),
         numeric(dim(tangent)[3]))
}

# The Hessian of w' x from the forward second derivatives.
hess_forward <- function(res, W) {
  ns <- dim(res$hessian)[3]
  outer(seq_len(ns), seq_len(ns),
        Vectorize(function(a, b) sum(as.vector(W) * res$hessian[, , a, b])))
}

expect_second_order <- function(ff, fr, W, info, tol_g = 1e-5, tol_h = 1e-5) {
  ref <- contract(ff$tangent, W)[, 1]
  expect_equal(unname(fr$cotangent[names(ref), 1]), unname(ref),
               tolerance = tol_g, info = info)
  expect_equal(unname(fr$curvature[, , 1]), hess_forward(ff, W),
               tolerance = tol_h, info = info)
}

test_that("forcings and a jump go backwards at second order", {
  # What works: a forcing, and a jump whose height is a parameter. What does
  # not: a jump whose TIME is a parameter, the case below.
  p  <- c(A = 1.2, B = 0.4, k1 = 0.7, k2 = 0.35, k3 = 1.1, d_amt = 0.4)
  fc <- list(u = data.frame(time = c(0, 2, 5), value = c(0.1, 0.25, 0.05)))

  ff <- do.call(solveODE, c(list(pr_forcing$ff, times, p, forcings = fc), tol))
  # The jump has to be in the run, or the test proves nothing.
  expect_gt(nrow(ff$variable), length(times))
  W  <- seed_for(nrow(ff$variable))
  for (m in c("bdf", "tsit5")) {
    fr <- do.call(solveODE, c(list(pr_forcing$fr[[m]], times, p, forcings = fc,
                                   cotangent = W), tol))
    expect_second_order(ff, fr, W, m, tol_h = 1e-4)
  }
})

test_that("a jump whose time is a parameter goes backwards at second order", {
  # The output grid has a row at t*, and that row's TIME moves with the
  # parameter. A cotangent on it makes w.x a different functional, and then no
  # derivative agrees with a difference quotient. The comparison is therefore
  # on the user times alone.
  p  <- c(A = 1.2, B = 0.4, k1 = 0.7, k2 = 0.35, k3 = 1.1, t_ev = 1.1)
  tl <- list(abstol = 1e-12, reltol = 1e-12)

  mf <- pr_ptime$ff
  mr <- pr_ptime$fr$bdf
  ff <- do.call(solveODE, c(list(mf, times, p), tl))
  expect_gt(nrow(ff$variable), length(times))

  W <- seed_for(nrow(ff$variable))
  moving <- which(vapply(ff$time, function(x) min(abs(x - times)) > 1e-9, TRUE))
  expect_length(moving, 1L)
  W[moving, , ] <- 0

  fr <- do.call(solveODE, c(list(mr, times, p, cotangent = W), tl))
  expect_second_order(ff, fr, W, "parameter event time", tol_h = 1e-6)
})

test_that("an event's root, time and value may be any expression", {
  # Every slot at once, because each leaves different terms at zero: a root
  # linear in x and blind to the clock zeroes two thirds of grad g_dot. Here g
  # is quadratic in B, reads A and t, the event time is nonlinear in a
  # parameter, both heights read the state, a parameter and the clock, and so
  # does f, so df/dt reaches the shift of t*. A clock-reading height rides on
  # roottol rather than reltol.
  p  <- c(A = 1.2, B = 0.4, k1 = 0.7, k2 = 0.35, k3 = 1.1,
          d_amt = 0.4, t_ev = 1.0)
  tl <- list(abstol = 1e-12, reltol = 1e-12, roottol = 1e-12)

  ff <- do.call(solveODE, c(list(pr_general$ff, times, p), tl))
  W <- seed_for(nrow(ff$variable))
  moving <- which(vapply(ff$time, function(x) min(abs(x - times)) > 1e-9, TRUE))
  # One row for the fixed jump and a pair for each of the two root crossings,
  # or the test proves nothing: both events have to be in the run.
  expect_length(moving, 5L)
  W[moving, , ] <- 0

  for (m in c("bdf", "adams", "rb4", "tsit5")) {
    fr <- do.call(solveODE, c(list(pr_general$fr[[m]], times, p, cotangent = W), tl))
    expect_second_order(ff, fr, W, m)
  }
})

test_that("a clock-reading jump height rides on roottol", {
  # It reads t* itself, where a height blind to the clock only reads the state
  # there, so the localisation error reaches it undamped. Not a missing term:
  # the gap falls with roottol and floors at reltol.
  p  <- c(A = 1.2, B = 0.4, k1 = 0.7, k2 = 0.35, k3 = 1.1, d_amt = 0.4)
  mf <- pr_rtol$ff
  mr <- pr_rtol$fr$bdf

  gap <- vapply(c(1e-6, 1e-9), function(rt) {
    tl <- list(abstol = 1e-12, reltol = 1e-12, roottol = rt)
    ff <- do.call(solveODE, c(list(mf, times, p), tl))
    W  <- seed_for(nrow(ff$variable))
    W[vapply(ff$time, function(x) min(abs(x - times)) > 1e-9, TRUE), , ] <- 0
    fr <- do.call(solveODE, c(list(mr, times, p, cotangent = W), tl))
    max(abs(unname(fr$curvature[, , 1]) - hess_forward(ff, W)))
  }, numeric(1))

  expect_lt(gap[2], gap[1] * 1e-2)
})

test_that("a sparse Jacobian goes backwards at second order", {
  skip_if_not(isTRUE(cppDE:::cvodeConfig$klu_available), "KLU not available")
  mf <- pr_sparse$ff
  mr <- pr_sparse$fr$bdf
  ff <- do.call(solveODE, c(list(mf, times, pars), tol))
  W  <- seed_for(nrow(ff$variable))
  fr <- do.call(solveODE, c(list(mr, times, pars, cotangent = W), tol))
  expect_second_order(ff, fr, W, "sparse")
})
