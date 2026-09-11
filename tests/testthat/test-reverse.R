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
  mr <- cppODE(eqns, modelname = "rev_plain_r", derivMode = "reverse")

  expect_identical(attr(mr, "derivMode"), "reverse")
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
  mr <- cppODE(eqns, modelname = "rev_scale_r", derivMode = "reverse")

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
  mr <- cppODE(eq, events = ev, forcings = "u", modelname = "rev_ev_r", derivMode = "reverse")

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
  mr <- cppODE(eqns, modelname = "rev_batch_r", derivMode = "reverse")

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
  mr <- cppODE(eqns, modelname = "rev_guard_r", derivMode = "reverse")

  fwd <- solveODE(mf, times, pars)
  W   <- seed_for(fwd)

  expect_error(solveODE(mf, times, pars, seed = W), "derivMode")
  expect_error(solveODE(mr, times, pars), "needs a 'seed'")
  expect_error(solveODE(mr, times, pars, seed = W[, 1, , drop = FALSE]),
               "state columns")
  expect_error(cppODE(eqns, modelname = "rev_no2nd", derivMode = "reverse",
                      deriv2 = TRUE),
               "second order")
  # The Rosenbrock replay takes a dense Jacobian only. The written adjoint does
  # not replay, so what is refused here is the combination that still does.
  expect_error(cppODE(eqns, modelname = "rev_rb4_sparse_ev", derivMode = "reverse",
                      method = "rb4", sparse = TRUE,
                      events = data.frame(var = "A", time = 1, value = 0.1,
                                          method = "add",
                                          stringsAsFactors = FALSE)),
               "no reverse mode")
})

test_that("a written Rosenbrock adjoint carries a multiplicative forcing", {
  # Four of the six stages add a multiple of df/dt, so its derivative in the
  # state and in the parameters is part of the adjoint. A forcing reaches df/dt
  # through a chain term the Jacobian emitter appends, and multiplicatively is
  # the way that term keeps a state in it.
  eq <- c(A = "u * A - k1 * A", B = "k1 * A - k2 * B * u")
  pf <- c(A = 1.1, B = 0.3, k1 = 0.8, k2 = 0.45)
  fc <- list(u = data.frame(time = c(0, 1, 3, 5), value = c(0.2, 0.5, 0.1, 0.3)))
  for (m in c("rb4", "bdf")) {
    mf <- cppODE(eq, forcings = "u", method = m,
                 modelname = paste0("rev_fc_f_", m), deriv = TRUE)
    mr <- cppODE(eq, forcings = "u", method = m,
                 modelname = paste0("rev_fc_r_", m), derivMode = "reverse")
    fwd <- do.call(solveODE, c(list(mf, times, pf, forcings = fc), tol))
    W   <- seed_for(fwd)
    rv  <- do.call(solveODE, c(list(mr, times, pf, forcings = fc, seed = W), tol))
    ref <- contract(fwd$sens1, W)[, 1]
    expect_equal(unname(rv$adjoint[names(ref), 1]), unname(ref),
                 tolerance = 1e-5, info = m)
  }
})

test_that("rb4 goes backwards on a sparse Jacobian", {
  mf <- cppODE(eqns, modelname = "rev_rb4_sp_f", method = "rb4", sparse = TRUE,
               deriv = TRUE)
  mr <- cppODE(eqns, modelname = "rev_rb4_sp_r", method = "rb4", sparse = TRUE,
               derivMode = "reverse")
  fwd <- do.call(solveODE, c(list(mf, times, pars), tol))
  W   <- seed_for(fwd)
  rv  <- do.call(solveODE, c(list(mr, times, pars, seed = W), tol))
  ref <- contract(fwd$sens1, W)[, 1]
  expect_equal(unname(rv$adjoint[names(ref), 1]), unname(ref), tolerance = 1e-5)
})

test_that("every method carries the reverse mode", {
  for (m in c("bdf", "adams", "rb4", "tsit5")) {
    mf <- cppODE(eqns, method = m, modelname = paste0("rev_m_f_", m), deriv = TRUE)
    mr <- cppODE(eqns, method = m, modelname = paste0("rev_m_r_", m),
                 derivMode = "reverse")
    fwd <- do.call(solveODE, c(list(mf, times, pars), tol))
    W   <- seed_for(fwd)
    rv  <- do.call(solveODE, c(list(mr, times, pars, seed = W), tol))
    ref <- contract(fwd$sens1, W)[, 1]
    expect_equal(unname(rv$adjoint[names(ref), 1]), unname(ref),
                 tolerance = 1e-5, info = m)
  }
})
# A model wider than the Nordsieck history is deep. Every other model in this
# file has fewer states than a step has slots, so a buffer sized by one and used
# for the other fits, and a contraction writing per state stays inside it.
wide  <- local({
  n <- 10L
  v <- paste0("A", seq_len(n))
  k <- paste0("k", seq_len(n))
  eq <- setNames(character(n), v)
  eq[1] <- paste0("-", k[1], "*", v[1])
  for (i in 2:n)
    eq[i] <- paste0(k[i - 1], "*", v[i - 1], " - ", k[i], "*", v[i], "*", v[i])
  eq
})
wpars <- c(setNames(seq(1.5, 0.6, length.out = 10), paste0("A", 1:10)),
           setNames(seq(0.9, 0.2, length.out = 10), paste0("k", 1:10)))

test_that("the reverse mode carries a model wider than its history is deep", {
  for (m in c("bdf", "adams")) {
    mf <- cppODE(wide, method = m, modelname = paste0("rev_wide_f_", m),
                 deriv = TRUE)
    mr <- cppODE(wide, method = m, modelname = paste0("rev_wide_r_", m),
                 derivMode = "reverse")
    fwd <- do.call(solveODE, c(list(mf, times, wpars), tol))
    W   <- seed_for(fwd)
    rv  <- do.call(solveODE, c(list(mr, times, wpars, seed = W), tol))
    ref <- contract(fwd$sens1, W)[, 1]
    expect_equal(unname(rv$adjoint[names(ref), 1]), unname(ref),
                 tolerance = 1e-5, info = m)
  }
})

test_that("adjointGrid reports the grid the sweep ran on", {
  mr <- cppODE(eqns, modelname = "rev_grid_r", derivMode = "reverse")
  mv <- cppODE(eqns, modelname = "rev_grid_v", deriv = FALSE)

  val <- do.call(solveODE, c(list(mv, times, pars), tol))
  W   <- seed_for(val, n_seed = 2L)

  plain <- do.call(solveODE, c(list(mr, times, pars, seed = W), tol))
  expect_null(plain$adjointGrid)

  rv <- do.call(solveODE,
                c(list(mr, times, pars, seed = W, adjointGrid = TRUE), tol))
  G  <- rv$adjointGrid
  expect_named(G, c("time", "h", "eta", "lambda"))

  n <- length(G$h)
  expect_gt(n, 4L)
  expect_equal(dim(G$eta),    c(n, 2L))
  expect_equal(dim(G$lambda), c(n, length(attr(mr, "variables")), 2L))

  # The reported grid is the run's own: it starts where the run started and is
  # gapless. It is not required to end on the last output time, because a
  # multistep method steps past it and interpolates back, so the span is a
  # lower bound.
  expect_equal(G$time[1], times[1])
  expect_equal(G$time[-1], head(G$time, -1) + head(G$h, -1), tolerance = 1e-8)
  expect_gte(sum(G$h), times[length(times)] - times[1])
  expect_lt(sum(G$h), 1.5 * (times[length(times)] - times[1]))

  # Turning the trace on must not move the answer, to the last bit: it is the
  # same sweep either way, with two more vectors written down.
  expect_identical(rv$adjoint, plain$adjoint)
  expect_equal(diagnostics(rv)$accepted, diagnostics(plain)$accepted)
})

test_that("the refinement indicator is alive on every method", {
  # eta = lambda' e_k is what a lambda-weighted step-size controller would read.
  # It reached the R seam as roundoff on bdf and adams until 2026-09-09: a
  # corrector method's step end is the value its equation solved for, and the
  # sweep leaves that adjoint at zero on its way past. Anything at 1e-15 here
  # means wout() or error_scale() has come undone again.
  for (m in c("bdf", "adams", "rb4", "tsit5")) {
    mr <- cppODE(eqns, method = m, modelname = paste0("rev_eta_", m),
                 derivMode = "reverse")
    mv <- cppODE(eqns, method = m, modelname = paste0("rev_eta_v_", m),
                 deriv = FALSE)
    val <- solveODE(mv, times, pars, abstol = 1e-8, reltol = 1e-6)
    W   <- seed_for(val)
    rv  <- solveODE(mr, times, pars, abstol = 1e-8, reltol = 1e-6,
                    seed = W, adjointGrid = TRUE)
    G <- rv$adjointGrid

    expect_gt(max(abs(G$eta)), 1e-12, label = paste("max |eta| on", m))
    expect_gt(max(abs(G$lambda)), 1e-3, label = paste("max |lambda| on", m))
    # Loosening the tolerance a hundredfold has to raise the indicator: it is
    # an error estimate, not a property of the trajectory.
    rv2 <- solveODE(mr, times, pars, abstol = 1e-6, reltol = 1e-4,
                    seed = W, adjointGrid = TRUE)
    expect_gt(sum(abs(rv2$adjointGrid$eta)), sum(abs(G$eta)), label = m)
  }
})
test_that("lambda weights can only refine the grid, never coarsen it", {
  # Stage 9 of dev/adjoint-plan.md. err = max(err_state, err_lambda), so the
  # grid stays at least as fine as abstol/rtol ask. That is the property the
  # whole scheme rests on: it makes a wrong weight cost time and never
  # accuracy, which is what lets a weight from an earlier run be used at all.
  mr <- cppODE(eqns, modelname = "rev_ew_r", derivMode = "reverse")
  mv <- cppODE(eqns, modelname = "rev_ew_v", deriv = FALSE)

  val <- solveODE(mv, times, pars, abstol = 1e-8, reltol = 1e-6)
  W   <- seed_for(val)

  base <- solveODE(mr, times, pars, abstol = 1e-8, reltol = 1e-6,
                   seed = W, adjointGrid = TRUE)
  G <- base$adjointGrid
  wts <- list(time = G$time, lambda = G$lambda[, , 1])

  # Loose enough that the weighted term cannot bind: the grid has to come back
  # bit for bit, or the term is doing something beyond the maximum.
  slack <- solveODE(mr, times, pars, abstol = 1e-8, reltol = 1e-6, seed = W,
                    adjointGrid = TRUE,
                    errWeights = c(wts, list(gradtol = 1e6)))
  expect_identical(slack$adjointGrid$h, G$h)
  expect_identical(slack$adjoint, base$adjoint)

  # Tight enough that it must bind, and it may only add steps.
  tight <- solveODE(mr, times, pars, abstol = 1e-8, reltol = 1e-6, seed = W,
                    adjointGrid = TRUE,
                    errWeights = c(wts, list(gradtol = 1e-11)))
  expect_gt(length(tight$adjointGrid$h), length(G$h))

  # A finer grid must not move the answer beyond the tolerance it was asked for.
  expect_equal(unname(tight$adjoint[, 1]), unname(base$adjoint[, 1]),
               tolerance = 1e-4)
})

test_that("every method takes lambda weights", {
  for (m in c("bdf", "adams", "rb4", "tsit5")) {
    mr <- cppODE(eqns, method = m, modelname = paste0("rev_ewm_", m),
                 derivMode = "reverse")
    mv <- cppODE(eqns, method = m, modelname = paste0("rev_ewv_", m),
                 deriv = FALSE)
    val <- solveODE(mv, times, pars, abstol = 1e-8, reltol = 1e-6)
    W   <- seed_for(val)
    b   <- solveODE(mr, times, pars, abstol = 1e-8, reltol = 1e-6,
                    seed = W, adjointGrid = TRUE)
    G   <- b$adjointGrid
    wts <- list(time = G$time, lambda = G$lambda[, , 1])

    slack <- solveODE(mr, times, pars, abstol = 1e-8, reltol = 1e-6, seed = W,
                      adjointGrid = TRUE,
                      errWeights = c(wts, list(gradtol = 1e6)))
    expect_identical(slack$adjointGrid$h, G$h, info = m)

    tight <- solveODE(mr, times, pars, abstol = 1e-8, reltol = 1e-6, seed = W,
                      adjointGrid = TRUE,
                      errWeights = c(wts, list(gradtol = 1e-11)))
    expect_gt(length(tight$adjointGrid$h), length(G$h))
  }
})

test_that("malformed lambda weights are an error, not a silent no-op", {
  mr <- cppODE(eqns, modelname = "rev_ewbad", derivMode = "reverse")
  mv <- cppODE(eqns, modelname = "rev_ewbadv", deriv = FALSE)
  val <- solveODE(mv, times, pars)
  W   <- seed_for(val)
  ok  <- list(time = c(0, 1, 2), lambda = matrix(1, 3, length(pars) - 3L))
  nx  <- length(attr(mr, "variables"))
  ok$lambda <- matrix(1, 3, nx)

  run <- function(w) solveODE(mr, times, pars, seed = W, errWeights = w)

  expect_error(run(list(time = c(0, 1))), "missing: lambda")
  expect_error(run(c(ok["lambda"], list(time = c(2, 1, 0)))), "ascending")
  expect_error(run(c(ok["time"], list(lambda = matrix(1, 2, nx)))),
               "rows but 'time' has")
  expect_error(run(c(ok["time"], list(lambda = matrix(1, 3, nx + 1L)))),
               "columns but the model has")
  expect_error(run(c(ok, list(gradtol = 0))), "must be positive")
  expect_error(run(c(ok, list(floor = 1))), "must be in")
  expect_error(run(c(ok, list(breaks = 99L))), "must index")
  # Weights without a seed weight nothing, so saying so beats ignoring them.
  expect_error(solveODE(mv, times, pars, errWeights = ok), "needs a 'seed'")
})
# -- CVODES adjoint sensitivity analysis, stage 8 -------------------------------
#
# ASA is not an oracle. It solves the adjoint as its own ODE over checkpointed
# forward states rather than differentiating the steps the forward run took, so
# it is a third discretisation. What it is good for is a second implementation
# of the same mathematics, by a different group, against which a systematic
# error in ours would show.

test_that("the CVODE backend takes derivatives backwards", {
  skip_if_not(isTRUE(cvodeConfig$available), "CVODE backend not available")

  mf <- cppODE(eqns, modelname = "asa_fwd", deriv = TRUE)
  ma <- cvode(eqns, modelname = "asa_rev", derivMode = "reverse")

  expect_identical(attr(ma, "derivMode"), "reverse")
  expect_false(attr(ma, "deriv"))

  fwd <- do.call(solveODE, c(list(mf, times, pars), tol))
  W   <- seed_for(fwd, n_seed = 2L)
  asa <- do.call(solveODE, c(list(ma, times, pars, seed = W), tol))

  expect_null(asa$sens1)
  expect_equal(dim(asa$adjoint), c(length(pars), 2L))
  expect_equal(asa$variable, fwd$variable, tolerance = 1e-7)

  ref <- contract(fwd$sens1, W)
  expect_equal(unname(asa$adjoint[rownames(ref), ]), unname(ref),
               tolerance = 1e-6)
})

test_that("ASA goes through the batch entry too", {
  skip_if_not(isTRUE(cvodeConfig$available), "CVODE backend not available")

  ma <- cvode(eqns, modelname = "asa_batch", derivMode = "reverse")
  mv <- cppODE(eqns, modelname = "asa_batch_v", deriv = FALSE)

  p2 <- pars; p2["k1"] <- 1.3
  conds <- list(one = list(parms = pars), two = list(parms = p2))
  W <- lapply(conds, function(cc)
    seed_for(do.call(solveODE, c(list(mv, times, cc$parms), tol))))

  one <- lapply(seq_along(conds), function(k)
    do.call(solveODE, c(list(ma, times, conds[[k]]$parms, seed = W[[k]]), tol)))
  bat <- do.call(solveODEBatch,
                 c(list(ma, mapply(function(cc, w) c(cc, list(seed = w)),
                                   conds, W, SIMPLIFY = FALSE),
                        times = times), tol))

  # The batch used to size its results before the solve, which left no slot for
  # the adjoint and reported success without it.
  for (k in seq_along(conds)) {
    expect_false(is.null(bat[[k]]$adjoint))
    expect_identical(bat[[k]]$adjoint, one[[k]]$adjoint)
  }
})

test_that("ASA and the native adjoint answer the same question", {
  skip_if_not(isTRUE(cvodeConfig$available), "CVODE backend not available")

  mv <- cppODE(eqns, modelname = "asa_cmp_v", deriv = FALSE)
  mr <- cppODE(eqns, modelname = "asa_cmp_r", derivMode = "reverse")
  ma <- cvode(eqns,  modelname = "asa_cmp_a", derivMode = "reverse")

  val <- do.call(solveODE, c(list(mv, times, pars), tol))
  W   <- seed_for(val)
  rev <- do.call(solveODE, c(list(mr, times, pars, seed = W), tol))
  asa <- do.call(solveODE, c(list(ma, times, pars, seed = W), tol))

  # Two implementations with nothing in common but the mathematics. They differ
  # by the discretisation, so this is a loose tolerance on purpose; a systematic
  # error in either would be orders wider.
  expect_equal(unname(asa$adjoint[, 1]),
               unname(rev$adjoint[rownames(asa$adjoint), 1]),
               tolerance = 1e-6)
})

test_that("the two CVODE directions refuse each other's arguments", {
  skip_if_not(isTRUE(cvodeConfig$available), "CVODE backend not available")

  # Second-order and forward sensitivities have no reverse counterpart here.
  expect_error(cvode(eqns, modelname = "asa_guard_d", derivMode = "reverse",
                     deriv = TRUE),
               "no forward sensitivities")

  # CVODES integrates the adjoint over checkpointed states, so a jump in the
  # state is a jump it cannot be told about. Saying so beats a wrong number.
  ev <- data.frame(var = "A", time = 1.0, value = 0.3, method = "add",
                   stringsAsFactors = FALSE)
  expect_error(cvode(eqns, events = ev, modelname = "asa_guard_e",
                     derivMode = "reverse"),
               "does not support events")

  mf <- cvode(eqns, modelname = "asa_guard_f")
  ma <- cvode(eqns, modelname = "asa_guard_r", derivMode = "reverse")
  val <- solveODE(mf, times, pars)
  W   <- seed_for(val)
  expect_error(solveODE(mf, times, pars, seed = W), "derivMode")
  expect_error(solveODE(ma, times, pars), "needs a .seed.")
})
test_that("the ASA backward problem gets the caller's step budget", {
  skip_if_not(isTRUE(cvodeConfig$available), "CVODE backend not available")

  # It did not, and kept the CVODES default of 500 steps while the forward pass
  # ran with the caller's 1e6. A stiff model reaches 500 long before t0, and
  # CVODES reports that as an integration failure rather than as a limit, so
  # the cause reads as the adjoint and is the budget. Robertson needs far more
  # than 500 backward steps, which is what makes it the test.
  rob <- c(y1 = "-k1*y1 + k2*y2*y3",
           y2 = "k1*y1 - k2*y2*y3 - k3*y2*y2",
           y3 = "k3*y2*y2")
  pr  <- c(y1 = 1, y2 = 0, y3 = 0, k1 = 0.04, k2 = 1e4, k3 = 3e7)
  tr  <- c(0, 10^seq(-5, 4, length.out = 20))
  otol <- list(abstol = 1e-10, reltol = 1e-8)

  mf <- cppODE(rob, modelname = "asa_ms_f", deriv = TRUE)
  ma <- cvode(rob,  modelname = "asa_ms_a", derivMode = "reverse")

  fwd <- do.call(solveODE, c(list(mf, tr, pr), otol))
  # A weight that no linear invariant annihilates: Robertson conserves
  # y1 + y2 + y3, so a constant seed would make the gradient exactly zero and
  # the test would pass on a solver that did nothing.
  W <- array(rep(1 / seq_len(3), each = nrow(fwd$variable)),
             c(nrow(fwd$variable), 3L, 1L))

  asa <- do.call(solveODE, c(list(ma, tr, pr, seed = W), otol))
  expect_equal(diagnostics(asa)$return_code, 0)

  ref <- contract(fwd$sens1, W)[, 1]
  expect_equal(unname(asa$adjoint[names(ref), 1]), unname(ref), tolerance = 1e-4)
})
test_that("a reverse solve can hand its checkpoints to the next one", {
  # A gradient needs two solves at the same theta: one for the values the seed
  # is built from, one for the sweep. Without this the states are integrated
  # twice. What has to hold is that the pair answers exactly what one seeded
  # call answers: the store is a saving, never a different number.
  mr <- cppODE(eqns, modelname = "rev_store_r", derivMode = "reverse")

  one <- do.call(solveODE, c(list(mr, times, pars, seed = NULL,
                                  keepStore = TRUE), tol))
  expect_identical(typeof(one$store), "externalptr")
  expect_null(one$adjoint)

  W <- seed_for(one, n_seed = 2L)
  plain <- do.call(solveODE, c(list(mr, times, pars, seed = W), tol))
  reuse <- do.call(solveODE, c(list(mr, times, pars, seed = W,
                                    store = one$store), tol))

  expect_equal(one$variable, plain$variable)
  expect_equal(reuse$variable, plain$variable)
  # Bit for bit: the same checkpoints replayed the same way, so anything but
  # exact equality means the reuse changed what was differentiated.
  expect_identical(reuse$adjoint, plain$adjoint)

  # The store may be spent more than once.
  again <- do.call(solveODE, c(list(mr, times, pars, seed = W,
                                    store = one$store), tol))
  expect_identical(again$adjoint, plain$adjoint)
})

test_that("a store from another point is refused, not quietly used", {
  # This is the whole safety of the scheme. A store carries the run it was made
  # from; handed back at a different theta it would give a gradient at one point
  # reported at another, and nothing downstream could tell.
  mr <- cppODE(eqns, modelname = "rev_store_g", derivMode = "reverse")
  one <- do.call(solveODE, c(list(mr, times, pars, keepStore = TRUE), tol))
  W   <- seed_for(one)

  p2 <- pars; p2["k1"] <- pars[["k1"]] * 1.1
  expect_error(do.call(solveODE, c(list(mr, times, p2, seed = W,
                                        store = one$store), tol)),
               "different times or parameters")

  t2 <- c(times, max(times) + 1)
  expect_error(do.call(solveODE, c(list(mr, t2, pars,
                                        seed = seed_for(one), store = one$store),
                                   tol)),
               "rows|different times or parameters")

  # And the two arguments belong to the reverse mode alone.
  mf <- cppODE(eqns, modelname = "rev_store_f", deriv = TRUE)
  expect_error(solveODE(mf, times, pars, keepStore = TRUE), "derivMode")
  expect_error(solveODE(mr, times, pars, store = "not a pointer", seed = W),
               "element of an earlier solve")
  expect_error(solveODE(mr, times, pars), "needs a 'seed'")
})
