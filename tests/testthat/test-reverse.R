# The reverse mode through the R seam: stage 6 of dev/adjoint-plan.md.
#
# The oracle is the forward mode. It is not sharp here and cannot be: the two
# solves adapt independently, the forward one under sensitivities and the
# reverse one in plain double, so they integrate two discretisations that differ
# by O(tol). What the tests below assert is that the gap falls with the
# tolerance the way the value gap does. Where the adjoint itself is checked on
# one shared step sequence, at rounding level, is dev/cxx/test_reverse_*.cpp.

skip_on_cran()

eqns <- c(A = "-k1 * A + k2 * B",
          B = "k1 * A - k2 * B - k3 * B * B")
pars  <- c(A = 1.2, B = 0.4, k1 = 0.7, k2 = 0.35, k3 = 1.1)
times <- c(0, 0.3, 1, 2.5, 5)
tol   <- list(abstol = 1e-10, reltol = 1e-10)

# w' S contracted over times and states, the quantity the sweep returns.
contract <- function(sens1, W) {
  vapply(seq_len(dim(W)[3]),
         function(k) apply(sens1 * as.vector(W[, , k]), 3, sum),
         numeric(dim(sens1)[3]))
}

seed_for <- function(res, n_seed = 1L, seed = 1L) {
  set.seed(seed)
  n <- nrow(res$variable); p <- ncol(res$variable)
  array(rnorm(n * p * n_seed), c(n, p, n_seed))
}

test_that("a reverse model answers what the forward sensitivities answer", {
  mf <- cppODE(eqns, modelname = "rev_plain_f", deriv = TRUE)
  mr <- cppODE(eqns, modelname = "rev_plain_r", sweep = "reverse")

  expect_identical(attr(mr, "sweep"), "reverse")
  expect_false(attr(mr, "deriv"))

  fwd <- do.call(solveODE, c(list(mf, times, pars), tol))
  W   <- seed_for(fwd, n_seed = 2L)
  rev <- do.call(solveODE, c(list(mr, times, pars, seed = W), tol))

  expect_null(rev$sens1)
  expect_equal(dim(rev$adjoint), c(length(pars), 2L))
  expect_equal(rownames(rev$adjoint), names(pars))
  expect_equal(rev$variable, fwd$variable, tolerance = 1e-7)

  ref <- contract(fwd$sens1, W)
  expect_equal(unname(rev$adjoint[rownames(ref), ]), unname(ref), tolerance = 1e-6)
})

test_that("the gap to the forward mode falls with the tolerance", {
  mf <- cppODE(eqns, modelname = "rev_scale_f", deriv = TRUE)
  mr <- cppODE(eqns, modelname = "rev_scale_r", sweep = "reverse")

  rel <- vapply(10^-c(6, 12), function(tt) {
    o   <- list(abstol = tt, reltol = tt)
    fwd <- do.call(solveODE, c(list(mf, times, pars), o))
    W   <- seed_for(fwd)
    rv  <- do.call(solveODE, c(list(mr, times, pars, seed = W), o))
    ref <- contract(fwd$sens1, W)[, 1]
    max(abs(ref - rv$adjoint[names(ref), 1])) / max(abs(ref))
  }, numeric(1))

  # Six decades of tolerance have to buy something close to six decades of
  # agreement. A missing channel would leave a floor instead.
  expect_lt(rel[2], rel[1] * 1e-4)
})

test_that("the reverse mode carries events, roots and forcings", {
  ev <- data.frame(var    = c("A", "B"),
                   time   = c("t_dose", NA),
                   value  = c("d_amt", 1.5),
                   root   = c(NA, "A - 0.6"),
                   method = c("add", "multiply"),
                   stringsAsFactors = FALSE)
  eq <- c(A = "-k1 * A + k2 * B + u",
          B = "k1 * A - k2 * B - k3 * B * B")
  p  <- c(A = 1.2, B = 0.4, k1 = 0.7, k2 = 0.35, k3 = 1.1,
          d_amt = 0.4, t_dose = 1.0)
  fc <- list(u = data.frame(time = c(0, 2, 5), value = c(0.1, 0.25, 0.05)))

  mf <- cppODE(eq, events = ev, forcings = "u", modelname = "rev_ev_f", deriv = TRUE)
  mr <- cppODE(eq, events = ev, forcings = "u", modelname = "rev_ev_r", sweep = "reverse")

  fwd <- do.call(solveODE, c(list(mf, times, p, forcings = fc), tol))
  # Both jumps have to be in the run, or the test proves nothing.
  expect_gt(nrow(fwd$variable), length(times))

  W   <- seed_for(fwd)
  rev <- do.call(solveODE, c(list(mr, times, p, forcings = fc, seed = W), tol))

  expect_equal(rev$variable, fwd$variable, tolerance = 1e-6)
  ref <- contract(fwd$sens1, W)[, 1]
  expect_equal(unname(rev$adjoint[names(ref), 1]), unname(ref), tolerance = 1e-5)
})

test_that("the reverse mode goes through the batch entry", {
  mf <- cppODE(eqns, modelname = "rev_batch_f", deriv = TRUE)
  mr <- cppODE(eqns, modelname = "rev_batch_r", sweep = "reverse")

  p2 <- pars; p2["k1"] <- 1.3
  conds <- list(one = list(parms = pars), two = list(parms = p2))

  fwd <- do.call(solveODEBatch, c(list(mf, conds, times = times), tol))
  W   <- lapply(fwd, seed_for)
  rev <- do.call(solveODEBatch,
                 c(list(mr, mapply(function(cc, w) c(cc, list(seed = w)),
                                   conds, W, SIMPLIFY = FALSE),
                        times = times), tol))

  for (k in seq_along(conds)) {
    ref <- contract(fwd[[k]]$sens1, W[[k]])[, 1]
    expect_equal(unname(rev[[k]]$adjoint[names(ref), 1]), unname(ref),
                 tolerance = 1e-6)
  }
})

test_that("the seed and the mode have to agree", {
  mf <- cppODE(eqns, modelname = "rev_guard_f", deriv = TRUE)
  mr <- cppODE(eqns, modelname = "rev_guard_r", sweep = "reverse")

  fwd <- solveODE(mf, times, pars)
  W   <- seed_for(fwd)

  expect_error(solveODE(mf, times, pars, seed = W), "sweep")
  expect_error(solveODE(mr, times, pars), "needs a 'seed'")
  expect_error(solveODE(mr, times, pars, seed = W[, 1, , drop = FALSE]),
               "state columns")
  expect_error(cppODE(eqns, modelname = "rev_no2nd", sweep = "reverse",
                      deriv2 = TRUE),
               "second order")
})

test_that("every method carries the reverse mode", {
  for (m in c("bdf", "adams", "rb4", "tsit5")) {
    mf <- cppODE(eqns, method = m, modelname = paste0("rev_m_f_", m), deriv = TRUE)
    mr <- cppODE(eqns, method = m, modelname = paste0("rev_m_r_", m),
                 sweep = "reverse")
    fwd <- do.call(solveODE, c(list(mf, times, pars), tol))
    W   <- seed_for(fwd)
    rv  <- do.call(solveODE, c(list(mr, times, pars, seed = W), tol))
    ref <- contract(fwd$sens1, W)[, 1]
    expect_equal(unname(rv$adjoint[names(ref), 1]), unname(ref),
                 tolerance = 1e-5, info = m)
  }
})
