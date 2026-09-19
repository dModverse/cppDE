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
methods <- c("bdf", "adams", "rb4", "tsit5")

# A dose and a root-triggered jump on a forced variant of eqns.
ev_eqns <- c(A = "-k1 * A + k2 * B + u",
             B = "k1 * A - k2 * B - k3 * B * B")
ev_def  <- data.frame(var    = c("A", "B"),
                      time   = c("t_dose", NA),
                      value  = c("d_amt", 1.5),
                      root   = c(NA, "A - 0.6"),
                      method = c("add", "multiply"),
                      stringsAsFactors = FALSE)
# A forcing that multiplies the states it acts on.
fc_eqns <- c(A = "u * A - k1 * A", B = "k1 * A - k2 * B * u")
# A cascade that settles, for rootfunc = "equilibrate".
eq_eqns <- c(R  = "k_act - k_deact * R",
             A  = "-k1 * A * R + k2 * pA",
             pA = "k1 * A * R - k2 * pA")

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

# One uncompiled cppODE model per method, named by method.
per_method <- function(rhs, prefix, ms = methods, ...)
  lapply(setNames(nm = ms), function(m)
    cppODE(rhs, method = m, modelname = paste0(prefix, m), compile = FALSE, ...))

# Links a list of lists of models into one shared object.
compile_all <- function(ms, output)
  do.call(compile, c(unname(unlist(ms, recursive = FALSE)),
                     list(output = output, cores = 1)))

# Every native model the tests solve without asking for KLU, one build per
# configuration. wide still turns sparse by auto-detection when KLU is present.
models <- list(
  fwd      = per_method(eqns, "rev_m_f_", deriv = TRUE),
  rev      = per_method(eqns, "rev_m_r_", derivMode = "reverse"),
  val      = per_method(eqns, "rev_m_v_", deriv = FALSE),
  ff       = per_method(eqns, "rev2_ff_", derivMode = "forward-forward"),
  fr       = per_method(eqns, "rev2_fr_", derivMode = "forward-reverse"),
  ev_fwd   = per_method(ev_eqns, "rev_ev_f_", events = ev_def,
                        forcings = "u", deriv = TRUE),
  ev_rev   = per_method(ev_eqns, "rev_ev_r_", events = ev_def,
                        forcings = "u", derivMode = "reverse"),
  fc_fwd   = per_method(fc_eqns, "rev_fc_f_", c("rb4", "bdf"),
                        forcings = "u", deriv = TRUE),
  fc_rev   = per_method(fc_eqns, "rev_fc_r_", c("rb4", "bdf"),
                        forcings = "u", derivMode = "reverse"),
  eq_val   = per_method(eq_eqns, "rev_eq_v_", c("bdf", "rb4"),
                        rootfunc = "equilibrate", deriv = FALSE),
  eq_rev   = per_method(eq_eqns, "rev_eq_r_", c("bdf", "rb4"),
                        rootfunc = "equilibrate", derivMode = "reverse"),
  wide_fwd = per_method(wide, "rev_wide_f_", c("bdf", "adams"), deriv = TRUE),
  wide_rev = per_method(wide, "rev_wide_r_", c("bdf", "adams"),
                        derivMode = "reverse"))
compile_all(models, "test_reverse")

# The sparse models need KLU. They all have the Jacobian of eqns, so the KLU
# defines each one records coincide and one shared object serves them all.
if (isTRUE(cvodeConfig$klu_available)) {
  sp_ev <- data.frame(var = "A", time = "t_dose", value = "d_amt",
                      method = "add", stringsAsFactors = FALSE)
  sparse_models <- list(
    ev_fwd = per_method(eqns, "rev_sp_ev_f_", c("bdf", "rb4"), events = sp_ev,
                        sparse = TRUE, deriv = TRUE),
    ev_rev = per_method(eqns, "rev_sp_ev_r_", c("bdf", "rb4"), events = sp_ev,
                        sparse = TRUE, derivMode = "reverse"),
    fwd    = per_method(eqns, "rev_sp_f_", "rb4", sparse = TRUE, deriv = TRUE),
    rev    = per_method(eqns, "rev_sp_r_", "rb4", sparse = TRUE,
                        derivMode = "reverse"))
  compile_all(sparse_models, "test_reverse_klu")
}

# w' S contracted over times and states, the quantity the sweep returns.
contract <- function(tangent, W) {
  vapply(seq_len(dim(W)[3]),
         function(k) apply(tangent * as.vector(W[, , k]), 3, sum),
         numeric(dim(tangent)[3]))
}

seed_for <- function(res, n_seed = 1L, seed = 1L) {
  set.seed(seed)
  n <- nrow(res$variable); p <- ncol(res$variable)
  array(rnorm(n * p * n_seed), c(n, p, n_seed))
}

test_that("a reverse model answers what the forward sensitivities answer", {
  mf <- models$fwd$bdf
  mr <- models$rev$bdf

  expect_identical(attr(mr, "derivMode"), "reverse")
  expect_false(attr(mr, "deriv"))

  fwd <- do.call(solveODE, c(list(mf, times, pars), tol))
  W   <- seed_for(fwd, n_seed = 2L)
  rev <- do.call(solveODE, c(list(mr, times, pars, cotangent = W), tol))

  expect_null(rev$tangent)
  expect_equal(dim(rev$cotangent), c(length(pars), 2L))
  expect_equal(rownames(rev$cotangent), names(pars))
  expect_equal(rev$variable, fwd$variable, tolerance = 1e-7)

  ref <- contract(fwd$tangent, W)
  expect_equal(unname(rev$cotangent[rownames(ref), ]), unname(ref), tolerance = 1e-6)
})

test_that("the gap to the forward mode falls with the tolerance", {
  mf <- models$fwd$bdf
  mr <- models$rev$bdf

  rel <- vapply(10^-c(6, 12), function(tt) {
    o   <- list(abstol = tt, reltol = tt)
    fwd <- do.call(solveODE, c(list(mf, times, pars), o))
    W   <- seed_for(fwd)
    rv  <- do.call(solveODE, c(list(mr, times, pars, cotangent = W), o))
    ref <- contract(fwd$tangent, W)[, 1]
    max(abs(ref - rv$cotangent[names(ref), 1])) / max(abs(ref))
  }, numeric(1))

  # Six decades of tolerance have to buy something close to six decades of
  # agreement. A missing channel would leave a floor instead.
  expect_lt(rel[2], rel[1] * 1e-4)
})

test_that("the reverse mode carries events, roots and forcings", {
  p  <- c(A = 1.2, B = 0.4, k1 = 0.7, k2 = 0.35, k3 = 1.1,
          d_amt = 0.4, t_dose = 1.0)
  fc <- list(u = data.frame(time = c(0, 2, 5), value = c(0.1, 0.25, 0.05)))

  for (m in methods) {
    mf <- models$ev_fwd[[m]]
    mr <- models$ev_rev[[m]]

    fwd <- do.call(solveODE, c(list(mf, times, p, forcings = fc), tol))
    # Both jumps have to be in the run, or the test proves nothing.
    expect_gt(nrow(fwd$variable), length(times))

    W   <- seed_for(fwd)
    rev <- do.call(solveODE, c(list(mr, times, p, forcings = fc, cotangent = W), tol))

    expect_equal(rev$variable, fwd$variable, tolerance = 1e-6, info = m)
    ref <- contract(fwd$tangent, W)[, 1]
    expect_equal(unname(rev$cotangent[names(ref), 1]), unname(ref),
                 tolerance = 1e-5, info = m)
  }
})

test_that("a sparse model jumps backwards too", {
  skip_if_not(isTRUE(cvodeConfig$klu_available), "KLU not available")
  pe <- c(pars, t_dose = 1.0, d_amt = 0.4)
  for (m in c("bdf", "rb4")) {
    mf <- sparse_models$ev_fwd[[m]]
    mr <- sparse_models$ev_rev[[m]]
    fwd <- do.call(solveODE, c(list(mf, times, pe), tol))
    W   <- seed_for(fwd)
    rv  <- do.call(solveODE, c(list(mr, times, pe, cotangent = W), tol))
    ref <- contract(fwd$tangent, W)[, 1]
    expect_equal(unname(rv$cotangent[names(ref), 1]), unname(ref),
                 tolerance = 1e-5, info = m)
  }
})

test_that("the reverse mode goes through the batch entry", {
  mf <- models$fwd$bdf
  mr <- models$rev$bdf

  p2 <- pars; p2["k1"] <- 1.3
  conds <- list(one = list(parms = pars), two = list(parms = p2))

  fwd <- do.call(solveODEBatch, c(list(mf, conds, times = times), tol))
  W   <- lapply(fwd, seed_for)
  rev <- do.call(solveODEBatch,
                 c(list(mr, mapply(function(cc, w) c(cc, list(cotangent = w)),
                                   conds, W, SIMPLIFY = FALSE),
                        times = times), tol))

  for (k in seq_along(conds)) {
    ref <- contract(fwd[[k]]$tangent, W[[k]])[, 1]
    expect_equal(unname(rev[[k]]$cotangent[names(ref), 1]), unname(ref),
                 tolerance = 1e-6)
  }
})

test_that("the cotangent and the mode have to agree", {
  mf <- models$fwd$bdf
  mr <- models$rev$bdf

  fwd <- solveODE(mf, times, pars)
  W   <- seed_for(fwd)

  expect_error(solveODE(mf, times, pars, cotangent = W), "derivMode")
  expect_error(solveODE(mr, times, pars), "needs a 'cotangent'")
  expect_error(solveODE(mr, times, pars, cotangent = W[, 1, , drop = FALSE]),
               "state columns")
  expect_error(cppODE(eqns, modelname = "rev_no2nd", derivMode = "reverse",
                      deriv2 = TRUE),
               "forward-reverse")
})

test_that("a written Rosenbrock adjoint carries a multiplicative forcing", {
  # Four of the six stages add a multiple of df/dt, so its derivative in the
  # state and in the parameters is part of the adjoint. A forcing reaches df/dt
  # through a chain term the Jacobian emitter appends, and multiplicatively is
  # the way that term keeps a state in it.
  pf <- c(A = 1.1, B = 0.3, k1 = 0.8, k2 = 0.45)
  fc <- list(u = data.frame(time = c(0, 1, 3, 5), value = c(0.2, 0.5, 0.1, 0.3)))
  for (m in c("rb4", "bdf")) {
    mf <- models$fc_fwd[[m]]
    mr <- models$fc_rev[[m]]
    fwd <- do.call(solveODE, c(list(mf, times, pf, forcings = fc), tol))
    W   <- seed_for(fwd)
    rv  <- do.call(solveODE, c(list(mr, times, pf, forcings = fc, cotangent = W), tol))
    ref <- contract(fwd$tangent, W)[, 1]
    expect_equal(unname(rv$cotangent[names(ref), 1]), unname(ref),
                 tolerance = 1e-5, info = m)
  }
})

test_that("an equilibrated run goes backwards", {
  # rootfunc stops the run when the state stops moving, so the output grid is
  # the run's own. A sensitivity run does not produce the same one: it adapts
  # on the tangents as well and reaches the root elsewhere, so there is no
  # forward gradient on this grid to compare against. What is asserted is that
  # the backward pass answers on the grid its own forward pass produced, and
  # that this pass is the value run to the last bit.
  pe <- c(R = 1, A = 1, pA = 0, k_act = 0.1, k_deact = 0.7, k1 = 0.1, k2 = 0.05)
  tt <- seq(0, 1e3, length.out = 50)
  for (m in c("bdf", "rb4")) {
    mv <- models$eq_val[[m]]
    mr <- models$eq_rev[[m]]
    val <- do.call(solveODE, c(list(mv, tt, pe), tol))
    expect_lt(nrow(val$variable), length(tt))   # it really did stop early
    W   <- seed_for(val)
    rv  <- do.call(solveODE, c(list(mr, tt, pe, cotangent = W), tol))
    expect_identical(rv$variable, val$variable)
    expect_equal(dim(rv$cotangent), c(length(pe), 1L), info = m)
    expect_true(all(is.finite(rv$cotangent)), info = m)
    expect_gt(max(abs(rv$cotangent)), 1e-6)
  }
})

test_that("rb4 goes backwards on a sparse Jacobian", {
  skip_if_not(isTRUE(cvodeConfig$klu_available), "KLU not available")
  mf <- sparse_models$fwd$rb4
  mr <- sparse_models$rev$rb4
  fwd <- do.call(solveODE, c(list(mf, times, pars), tol))
  W   <- seed_for(fwd)
  rv  <- do.call(solveODE, c(list(mr, times, pars, cotangent = W), tol))
  ref <- contract(fwd$tangent, W)[, 1]
  expect_equal(unname(rv$cotangent[names(ref), 1]), unname(ref), tolerance = 1e-5)
})

test_that("every method carries the reverse mode", {
  for (m in methods) {
    mf <- models$fwd[[m]]
    mr <- models$rev[[m]]
    fwd <- do.call(solveODE, c(list(mf, times, pars), tol))
    W   <- seed_for(fwd)
    rv  <- do.call(solveODE, c(list(mr, times, pars, cotangent = W), tol))
    ref <- contract(fwd$tangent, W)[, 1]
    expect_equal(unname(rv$cotangent[names(ref), 1]), unname(ref),
                 tolerance = 1e-5, info = m)
  }
})

test_that("the reverse mode carries a model wider than its history is deep", {
  for (m in c("bdf", "adams")) {
    mf <- models$wide_fwd[[m]]
    mr <- models$wide_rev[[m]]
    fwd <- do.call(solveODE, c(list(mf, times, wpars), tol))
    W   <- seed_for(fwd)
    rv  <- do.call(solveODE, c(list(mr, times, wpars, cotangent = W), tol))
    ref <- contract(fwd$tangent, W)[, 1]
    expect_equal(unname(rv$cotangent[names(ref), 1]), unname(ref),
                 tolerance = 1e-5, info = m)
  }
})

test_that("adjointGrid reports the grid the sweep ran on", {
  mr <- models$rev$bdf
  mv <- models$val$bdf

  val <- do.call(solveODE, c(list(mv, times, pars), tol))
  W   <- seed_for(val, n_seed = 2L)

  plain <- do.call(solveODE, c(list(mr, times, pars, cotangent = W), tol))
  expect_null(plain$adjointGrid)

  rv <- do.call(solveODE,
                c(list(mr, times, pars, cotangent = W, adjointGrid = TRUE), tol))
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
  expect_identical(rv$cotangent, plain$cotangent)
  expect_equal(diagnostics(rv)$accepted, diagnostics(plain)$accepted)
})

test_that("the refinement indicator is alive on every method", {
  # eta = lambda' e_k is what a lambda-weighted step-size controller would read.
  # It reached the R seam as roundoff on bdf and adams until 2026-09-09: a
  # corrector method's step end is the value its equation solved for, and the
  # sweep leaves that adjoint at zero on its way past. Anything at 1e-15 here
  # means wout() or error_scale() has come undone again.
  for (m in methods) {
    mr <- models$rev[[m]]
    mv <- models$val[[m]]
    val <- solveODE(mv, times, pars, abstol = 1e-8, reltol = 1e-6)
    W   <- seed_for(val)
    rv  <- solveODE(mr, times, pars, abstol = 1e-8, reltol = 1e-6,
                    cotangent = W, adjointGrid = TRUE)
    G <- rv$adjointGrid

    expect_gt(max(abs(G$eta)), 1e-12, label = paste("max |eta| on", m))
    expect_gt(max(abs(G$lambda)), 1e-3, label = paste("max |lambda| on", m))
    # Loosening the tolerance a hundredfold has to raise the indicator: it is
    # an error estimate, not a property of the trajectory.
    rv2 <- solveODE(mr, times, pars, abstol = 1e-6, reltol = 1e-4,
                    cotangent = W, adjointGrid = TRUE)
    expect_gt(sum(abs(rv2$adjointGrid$eta)), sum(abs(G$eta)), label = m)
  }
})
test_that("lambda weights can only refine the grid, never coarsen it", {
  # Stage 9 of dev/adjoint-plan.md. err = max(err_state, err_lambda), so the
  # grid stays at least as fine as abstol/rtol ask. That is the property the
  # whole scheme rests on: it makes a wrong weight cost time and never
  # accuracy, which is what lets a weight from an earlier run be used at all.
  mr <- models$rev$bdf
  mv <- models$val$bdf

  val <- solveODE(mv, times, pars, abstol = 1e-8, reltol = 1e-6)
  W   <- seed_for(val)

  base <- solveODE(mr, times, pars, abstol = 1e-8, reltol = 1e-6,
                   cotangent = W, adjointGrid = TRUE)
  G <- base$adjointGrid
  wts <- list(time = G$time, lambda = G$lambda[, , 1])

  # Loose enough that the weighted term cannot bind: the grid has to come back
  # bit for bit, or the term is doing something beyond the maximum.
  slack <- solveODE(mr, times, pars, abstol = 1e-8, reltol = 1e-6, cotangent = W,
                    adjointGrid = TRUE,
                    errWeights = c(wts, list(gradtol = 1e6)))
  expect_identical(slack$adjointGrid$h, G$h)
  expect_identical(slack$cotangent, base$cotangent)

  # Tight enough that it must bind, and it may only add steps.
  tight <- solveODE(mr, times, pars, abstol = 1e-8, reltol = 1e-6, cotangent = W,
                    adjointGrid = TRUE,
                    errWeights = c(wts, list(gradtol = 1e-11)))
  expect_gt(length(tight$adjointGrid$h), length(G$h))

  # A finer grid must not move the answer beyond the tolerance it was asked for.
  expect_equal(unname(tight$cotangent[, 1]), unname(base$cotangent[, 1]),
               tolerance = 1e-4)
})

test_that("every method takes lambda weights", {
  for (m in methods) {
    mr <- models$rev[[m]]
    mv <- models$val[[m]]
    val <- solveODE(mv, times, pars, abstol = 1e-8, reltol = 1e-6)
    W   <- seed_for(val)
    b   <- solveODE(mr, times, pars, abstol = 1e-8, reltol = 1e-6,
                    cotangent = W, adjointGrid = TRUE)
    G   <- b$adjointGrid
    wts <- list(time = G$time, lambda = G$lambda[, , 1])

    slack <- solveODE(mr, times, pars, abstol = 1e-8, reltol = 1e-6, cotangent = W,
                      adjointGrid = TRUE,
                      errWeights = c(wts, list(gradtol = 1e6)))
    expect_identical(slack$adjointGrid$h, G$h, info = m)

    tight <- solveODE(mr, times, pars, abstol = 1e-8, reltol = 1e-6, cotangent = W,
                      adjointGrid = TRUE,
                      errWeights = c(wts, list(gradtol = 1e-11)))
    expect_gt(length(tight$adjointGrid$h), length(G$h))
  }
})

test_that("malformed lambda weights are an error, not a silent no-op", {
  mr <- models$rev$bdf
  mv <- models$val$bdf
  val <- solveODE(mv, times, pars)
  W   <- seed_for(val)
  ok  <- list(time = c(0, 1, 2), lambda = matrix(1, 3, length(pars) - 3L))
  nx  <- length(attr(mr, "variables"))
  ok$lambda <- matrix(1, 3, nx)

  run <- function(w) solveODE(mr, times, pars, cotangent = W, errWeights = w)

  expect_error(run(list(time = c(0, 1))), "missing: lambda")
  expect_error(run(c(ok["lambda"], list(time = c(2, 1, 0)))), "ascending")
  expect_error(run(c(ok["time"], list(lambda = matrix(1, 2, nx)))),
               "rows but 'time' has")
  expect_error(run(c(ok["time"], list(lambda = matrix(1, 3, nx + 1L)))),
               "columns but the model has")
  expect_error(run(c(ok, list(gradtol = 0))), "must be positive")
  expect_error(run(c(ok, list(floor = 1))), "must be in")
  expect_error(run(c(ok, list(breaks = 99L))), "must index")
  # Weights without a cotangent weight nothing: saying so beats ignoring them.
  expect_error(solveODE(mv, times, pars, errWeights = ok), "needs a 'cotangent'")
})
# -- CVODES adjoint sensitivity analysis, stage 8 -------------------------------
#
# ASA is not an oracle. It solves the adjoint as its own ODE over checkpointed
# forward states rather than differentiating the steps the forward run took, so
# it is a third discretisation. What it is good for is a second implementation
# of the same mathematics, by a different group, against which a systematic
# error in ours would show.

# The CVODES models and the native forward model of Robertson that only the
# ASA tests use, linked against SUNDIALS into one shared object.
if (isTRUE(cvodeConfig$available)) {
  rob <- c(y1 = "-k1*y1 + k2*y2*y3",
           y2 = "k1*y1 - k2*y2*y3 - k3*y2*y2",
           y3 = "k3*y2*y2")
  asa_models <- list(
    rev   = cvode(eqns, modelname = "asa_rev", derivMode = "reverse",
                  compile = FALSE),
    plain = cvode(eqns, modelname = "asa_plain", compile = FALSE),
    rob_f = cppODE(rob, modelname = "asa_ms_f", deriv = TRUE, compile = FALSE),
    rob_a = cvode(rob, modelname = "asa_ms_a", derivMode = "reverse",
                  compile = FALSE))
  compile_all(list(asa_models), "test_reverse_cvode")
}

test_that("the CVODE backend takes derivatives backwards", {
  skip_if_not(isTRUE(cvodeConfig$available), "CVODE backend not available")

  mf <- models$fwd$bdf
  ma <- asa_models$rev

  expect_identical(attr(ma, "derivMode"), "reverse")
  expect_false(attr(ma, "deriv"))

  fwd <- do.call(solveODE, c(list(mf, times, pars), tol))
  W   <- seed_for(fwd, n_seed = 2L)
  asa <- do.call(solveODE, c(list(ma, times, pars, cotangent = W), tol))

  expect_null(asa$tangent)
  expect_equal(dim(asa$cotangent), c(length(pars), 2L))
  expect_equal(asa$variable, fwd$variable, tolerance = 1e-7)

  ref <- contract(fwd$tangent, W)
  expect_equal(unname(asa$cotangent[rownames(ref), ]), unname(ref),
               tolerance = 1e-6)
})

test_that("ASA goes through the batch entry too", {
  skip_if_not(isTRUE(cvodeConfig$available), "CVODE backend not available")

  ma <- asa_models$rev
  mv <- models$val$bdf

  p2 <- pars; p2["k1"] <- 1.3
  conds <- list(one = list(parms = pars), two = list(parms = p2))
  W <- lapply(conds, function(cc)
    seed_for(do.call(solveODE, c(list(mv, times, cc$parms), tol))))

  one <- lapply(seq_along(conds), function(k)
    do.call(solveODE, c(list(ma, times, conds[[k]]$parms, cotangent = W[[k]]), tol)))
  bat <- do.call(solveODEBatch,
                 c(list(ma, mapply(function(cc, w) c(cc, list(cotangent = w)),
                                   conds, W, SIMPLIFY = FALSE),
                        times = times), tol))

  # The batch used to size its results before the solve, which left no slot for
  # the cotangent and reported success without it.
  for (k in seq_along(conds)) {
    expect_false(is.null(bat[[k]]$cotangent))
    expect_identical(bat[[k]]$cotangent, one[[k]]$cotangent)
  }
})

test_that("ASA and the native adjoint answer the same question", {
  skip_if_not(isTRUE(cvodeConfig$available), "CVODE backend not available")

  mv <- models$val$bdf
  mr <- models$rev$bdf
  ma <- asa_models$rev

  val <- do.call(solveODE, c(list(mv, times, pars), tol))
  W   <- seed_for(val)
  rev <- do.call(solveODE, c(list(mr, times, pars, cotangent = W), tol))
  asa <- do.call(solveODE, c(list(ma, times, pars, cotangent = W), tol))

  # Two implementations with nothing in common but the mathematics. They differ
  # by the discretisation, so this is a loose tolerance on purpose; a systematic
  # error in either would be orders wider.
  expect_equal(unname(asa$cotangent[, 1]),
               unname(rev$cotangent[rownames(asa$cotangent), 1]),
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

  mf <- asa_models$plain
  ma <- asa_models$rev
  val <- solveODE(mf, times, pars)
  W   <- seed_for(val)
  expect_error(solveODE(mf, times, pars, cotangent = W), "derivMode")
  expect_error(solveODE(ma, times, pars), "needs a .cotangent.")
})
test_that("the ASA backward problem gets the caller's step budget", {
  skip_if_not(isTRUE(cvodeConfig$available), "CVODE backend not available")

  # It did not, and kept the CVODES default of 500 steps while the forward pass
  # ran with the caller's 1e6. A stiff model reaches 500 long before t0, and
  # CVODES reports that as an integration failure rather than as a limit, so
  # the cause reads as the adjoint and is the budget. Robertson needs far more
  # than 500 backward steps, which is what makes it the test.
  pr  <- c(y1 = 1, y2 = 0, y3 = 0, k1 = 0.04, k2 = 1e4, k3 = 3e7)
  tr  <- c(0, 10^seq(-5, 4, length.out = 20))
  otol <- list(abstol = 1e-10, reltol = 1e-8)

  mf <- asa_models$rob_f
  ma <- asa_models$rob_a

  fwd <- do.call(solveODE, c(list(mf, tr, pr), otol))
  # A weight that no linear invariant annihilates: Robertson conserves
  # y1 + y2 + y3, so a constant cotangent would make the gradient exactly zero
  # and the test would pass on a solver that did nothing.
  W <- array(rep(1 / seq_len(3), each = nrow(fwd$variable)),
             c(nrow(fwd$variable), 3L, 1L))

  asa <- do.call(solveODE, c(list(ma, tr, pr, cotangent = W), otol))
  expect_equal(diagnostics(asa)$return_code, 0)

  ref <- contract(fwd$tangent, W)[, 1]
  expect_equal(unname(asa$cotangent[names(ref), 1]), unname(ref), tolerance = 1e-4)
})
test_that("a reverse solve can hand its checkpoints to the next one", {
  # A gradient needs two solves at the same theta: one for the values the
  # cotangent is built from, one for the sweep. Without this the states are
  # integrated twice. What has to hold is that the pair answers exactly what one
  # seeded call answers: the store is a saving, never a different number.
  mr <- models$rev$bdf

  one <- do.call(solveODE, c(list(mr, times, pars, cotangent = NULL,
                                  keepStore = TRUE), tol))
  expect_identical(typeof(one$store), "externalptr")
  expect_null(one$cotangent)

  W <- seed_for(one, n_seed = 2L)
  plain <- do.call(solveODE, c(list(mr, times, pars, cotangent = W), tol))
  reuse <- do.call(solveODE, c(list(mr, times, pars, cotangent = W,
                                    store = one$store), tol))

  expect_equal(one$variable, plain$variable)
  expect_equal(reuse$variable, plain$variable)
  # Bit for bit: the same checkpoints replayed the same way, so anything but
  # exact equality means the reuse changed what was differentiated.
  expect_identical(reuse$cotangent, plain$cotangent)

  # The store may be spent more than once.
  again <- do.call(solveODE, c(list(mr, times, pars, cotangent = W,
                                    store = one$store), tol))
  expect_identical(again$cotangent, plain$cotangent)
})

test_that("a store from another point is refused, not quietly used", {
  # This is the whole safety of the scheme. A store carries the run it was made
  # from; handed back at a different theta it would give a gradient at one point
  # reported at another, and nothing downstream could tell.
  mr <- models$rev$bdf
  one <- do.call(solveODE, c(list(mr, times, pars, keepStore = TRUE), tol))
  W   <- seed_for(one)

  p2 <- pars; p2["k1"] <- pars[["k1"]] * 1.1
  expect_error(do.call(solveODE, c(list(mr, times, p2, cotangent = W,
                                        store = one$store), tol)),
               "different times or parameters")

  t2 <- c(times, max(times) + 1)
  expect_error(do.call(solveODE, c(list(mr, t2, pars,
                                        cotangent = seed_for(one), store = one$store),
                                   tol)),
               "rows|different times or parameters")

  # And the two arguments belong to the reverse mode alone.
  mf <- models$fwd$bdf
  expect_error(solveODE(mf, times, pars, keepStore = TRUE), "derivMode")
  expect_error(solveODE(mr, times, pars, store = "not a pointer", cotangent = W),
               "element of an earlier solve")
  expect_error(solveODE(mr, times, pars), "needs a 'cotangent'")
})

# ---------------------------------------------------------------------------
#  Second order: the same Hessian from both directions
#
#  forward-forward propagates a nested dual through the states and answers with
#  the hessian; forward-reverse runs the backward sweep itself over a dual and
#  answers with the curvature. Contracted against the same cotangent the two are
#  the same matrix, so each is the other's oracle, and the gap is the
#  discretisation gap the file header describes.
# ---------------------------------------------------------------------------

# The Hessian of w' x, from the forward second derivatives.
hess_forward <- function(res, W) {
  n_sens <- dim(res$hessian)[3]
  outer(seq_len(n_sens), seq_len(n_sens),
        Vectorize(function(a, b) sum(as.vector(W) * res$hessian[, , a, b])))
}

test_that("forward-reverse answers the Hessian forward-forward answers", {
  for (m in c("bdf", "rb4", "tsit5")) {
    mf <- models$ff[[m]]
    mr <- models$fr[[m]]

    expect_identical(attr(mf, "derivMode"), "forward-forward")
    expect_true(attr(mf, "deriv2"))
    expect_identical(attr(mr, "derivMode"), "forward-reverse")
    expect_true(attr(mr, "deriv"))
    expect_false(attr(mr, "deriv2"))

    ff <- do.call(solveODE, c(list(mf, times, pars), tol))
    W  <- seed_for(ff)
    fr <- do.call(solveODE, c(list(mr, times, pars, cotangent = W), tol))

    # The gradient first: forward-reverse keeps the first-order answer.
    expect_equal(unname(fr$cotangent[, 1]),
                 unname(contract(ff$tangent, W)[, 1]),
                 tolerance = 1e-6, info = m)
    expect_equal(unname(fr$curvature[, , 1]), unname(hess_forward(ff, W)),
                 tolerance = 1e-6, info = m)
  }
})

test_that("every method answers the same Hessian backwards", {
  # Four discretisations of one Hessian. Sharper than it looks: the four adapt
  # their own grids, so agreement at 1e-5 says the sweep differentiates what
  # each of them actually did.
  W <- NULL; ref <- NULL
  for (m in methods) {
    mr <- models$fr[[m]]
    if (is.null(W)) {
      W <- seed_for(do.call(solveODE, c(list(mr, times, pars,
                                             cotangent = array(0, c(length(times), 2L, 1L))), tol)))
    }
    fr <- do.call(solveODE, c(list(mr, times, pars, cotangent = W), tol))
    if (is.null(ref)) ref <- unname(fr$curvature[, , 1])
    else expect_equal(unname(fr$curvature[, , 1]), ref, tolerance = 1e-5, info = m)
  }
})

test_that("forward-reverse names its answer and refuses the older spelling", {
  mr <- models$fr$bdf
  ff <- do.call(solveODE, c(list(mr, times, pars,
                                 cotangent = array(1, c(length(times), 2L, 1L))), tol))
  expect_identical(dim(ff$curvature), c(5L, 5L, 1L))
  expect_identical(dimnames(ff$curvature)[[1]], names(pars))
  expect_identical(dimnames(ff$curvature)[[2]], names(pars))

  expect_error(cppODE(eqns, modelname = "rev2_refused", derivMode = "reverse",
                      deriv2 = TRUE),
               "forward-reverse")
})

test_that("the second-order forward mode is repeatable on every method", {
  # A solve must not depend on what ran before it. The dense-output wrapper owns
  # the two state buffers a single-step method reads and writes, and while it
  # went unprimed those buffers took their tangents from the arena: values and
  # first derivatives were bit-identical between repeats, second derivatives
  # were not. Found through cppODE's CPPDE_POISON_ARENA switch.
  for (m in methods) {
    mm <- models$ff[[m]]
    a <- do.call(solveODE, c(list(mm, times, pars), tol))
    b <- do.call(solveODE, c(list(mm, times, pars), tol))
    expect_identical(a$hessian, b$hessian, info = m)
  }
})
