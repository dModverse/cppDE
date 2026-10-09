# The CVODES adjoint (ASA) through the R seam.
#
# ASA is not an oracle. It solves the adjoint as its own ODE over checkpointed
# forward states rather than differentiating the steps the forward run took, so
# it is a third discretisation. What it is good for is a second implementation
# of the same mathematics, by a different group, against which a systematic
# error in ours would show.

skip_on_cran()

eqns <- c(A = "-k1 * A + k2 * B",
          B = "k1 * A - k2 * B - k3 * B * B")
pars  <- c(A = 1.2, B = 0.4, k1 = 0.7, k2 = 0.35, k3 = 1.1)
times <- c(0, 0.3, 1, 2.5, 5)
tol   <- list(abstol = 1e-10, reltol = 1e-10)

# The CVODES models, and the native models they are checked against, linked
# against SUNDIALS into one shared object.
if (isTRUE(cvodeConfig$available)) {
  rob <- c(y1 = "-k1*y1 + k2*y2*y3",
           y2 = "k1*y1 - k2*y2*y3 - k3*y2*y2",
           y3 = "k3*y2*y2")
  asa_models <- list(
    rev   = cvode(eqns, modelname = "asa_rev", derivMode = "reverse",
                  compile = FALSE),
    plain = cvode(eqns, modelname = "asa_plain", compile = FALSE),
    fwd   = cppODE(eqns, modelname = "asa_nat_f", deriv = TRUE, compile = FALSE),
    nat   = cppODE(eqns, modelname = "asa_nat_r", derivMode = "reverse",
                   compile = FALSE),
    rob_f = cppODE(rob, modelname = "asa_ms_f", deriv = TRUE, compile = FALSE),
    rob_a = cvode(rob, modelname = "asa_ms_a", derivMode = "reverse",
                  compile = FALSE))
  if (isTRUE(cvodeConfig$klu_available))
    asa_models$rob_s <- cvode(rob, modelname = "asa_ms_s", derivMode = "reverse",
                              sparse = TRUE, compile = FALSE)
  # Events for ASA: each with its forward-sensitivity counterpart.
  asa_events <- list(
    add  = data.frame(var = "A", time = 1, value = "d", method = "add"),
    tpar = data.frame(var = "A", time = "td", value = "d", method = "add"),
    root = data.frame(var = "B", time = NA, value = "r", root = "A - 0.6",
                      method = "replace"),
    repl = data.frame(var = c("B", "A"), time = c("2", "td"),
                      value = c("r*A", "d*B"), method = c("replace", "multiply")))
  for (n in names(asa_events)) {
    asa_models[[paste0("ev_", n, "_a")]] <-
      cvode(eqns, events = asa_events[[n]], derivMode = "reverse",
            modelname = paste0("asa_ev_", n, "_a"), compile = FALSE)
    asa_models[[paste0("ev_", n, "_f")]] <-
      cvode(eqns, events = asa_events[[n]], deriv = TRUE,
            modelname = paste0("asa_ev_", n, "_f"), compile = FALSE)
  }
  if (isTRUE(cvodeConfig$klu_available)) {
    asa_models$ev_klu_a <- cvode(eqns, events = asa_events$root,
                                 derivMode = "reverse", sparse = TRUE,
                                 modelname = "asa_ev_klu_a", compile = FALSE)
    asa_models$ev_klu_f <- cvode(eqns, events = asa_events$root, deriv = TRUE,
                                 sparse = TRUE, modelname = "asa_ev_klu_f",
                                 compile = FALSE)
  }
  do.call(compile, c(unname(asa_models),
                     list(output = "test_reverse_cvode", cores = test_cores())))
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

test_that("the CVODE backend takes derivatives backwards", {
  skip_if_not(isTRUE(cvodeConfig$available), "CVODE backend not available")

  mf <- asa_models$fwd
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
  mv <- asa_models$fwd

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

test_that("ASA on KLU replays its checkpoints and matches the dense solve", {
  skip_if_not(isTRUE(cvodeConfig$klu_available), "KLU not available")
  # The replay from a checkpoint must take the run's steps; KLU's refactorisation
  # depends on its history, so every setup factorises afresh.
  src <- readLines(attr(asa_models$rob_s, "srcfile"))
  expect_true(any(grepl("LS->ops->setup = klu_setup_fresh", src, fixed = TRUE)))
  expect_true(any(grepl("SUNLinSol_KLUReInit", src, fixed = TRUE)))
  rt <- c(0, 10^seq(-2, 3, length.out = 30))
  rp <- c(y1 = 1, y2 = 0, y3 = 0, k1 = 0.04, k2 = 1e4, k3 = 3e7)
  W  <- cbind(y1 = rep(1, length(rt)), y2 = 1e4, y3 = 1)
  ctl <- adjointControl(gradtol = 1e-12)
  ad <- solveODE(asa_models$rob_a, rt, rp, cotangent = W, abstol = 1e-12, reltol = 1e-10, adjoint = ctl)
  as <- solveODE(asa_models$rob_s, rt, rp, cotangent = W, abstol = 1e-12, reltol = 1e-10, adjoint = ctl)
  expect_equal(as$cotangent, ad$cotangent, tolerance = 1e-6)
})

test_that("ASA passes the adjoint through events", {
  skip_if_not(isTRUE(cvodeConfig$available), "CVODE backend not available")
  # A dose, a dose at a parameter time, a root that replaces, and a replace and
  # a multiply at two times; ASA against the forward sensitivities of CVODES.
  ep <- c(pars, d = 0.5, td = 1.7, r = 0.3)
  et <- c(0, 0.3, 1, 2, 2.5, 5)
  ctl <- adjointControl(gradtol = 1e-12)
  cases <- c(names(asa_events), if (!is.null(asa_models$ev_klu_a)) "klu")
  for (n in cases) {
    ma <- asa_models[[paste0("ev_", n, "_a")]]
    mf <- asa_models[[paste0("ev_", n, "_f")]]
    p  <- ep[unique(c("A", "B", attr(ma, "parameters")))]
    f  <- solveODE(mf, et, p, abstol = 1e-12, reltol = 1e-10)
    set.seed(3)
    W <- array(rnorm(2 * length(f$variable)), c(dim(f$variable), 2),
               dimnames = list(NULL, colnames(f$variable), NULL))
    a <- solveODE(ma, et, p, cotangent = W, abstol = 1e-12, reltol = 1e-10,
                  adjoint = ctl)
    for (j in 1:2) {
      ref <- apply(f$tangent * as.vector(W[, , j]), 3, sum)
      expect_equal(unname(a$cotangent[names(ref), j]), unname(ref),
                   tolerance = 1e-7, label = paste(n, "column", j))
    }
  }
})

test_that("ASA and the native adjoint answer the same question", {
  skip_if_not(isTRUE(cvodeConfig$available), "CVODE backend not available")

  mv <- asa_models$fwd
  mr <- asa_models$nat
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

test_that("ASA holds the gradient to gradtol and refuses the rest", {
  skip_if_not(isTRUE(cvodeConfig$available), "CVODE backend not available")
  ma <- asa_models$rev
  W  <- seed_for(solveODE(asa_models$fwd, times, pars))
  ref <- solveODE(asa_models$nat, times, pars, cotangent = W, abstol = 1e-13,
                  reltol = 1e-12)$cotangent
  err <- function(r) max(abs(r$cotangent[, 1] - ref[rownames(r$cotangent), 1]))
  loose <- solveODE(ma, times, pars, cotangent = W, abstol = 1e-6, reltol = 1e-4)
  held  <- solveODE(ma, times, pars, cotangent = W, abstol = 1e-6, reltol = 1e-4,
                    adjoint = adjointControl(gradtol = 1e-9))
  expect_lt(err(held), err(loose))
  for (a in list(adjointControl(trace = TRUE), adjointControl(refine = TRUE)))
    expect_error(solveODE(ma, times, pars, cotangent = W, adjoint = a),
                 "'gradtol' alone")
})

test_that("the two CVODE directions refuse each other's arguments", {
  skip_if_not(isTRUE(cvodeConfig$available), "CVODE backend not available")

  # Second-order and forward sensitivities have no reverse counterpart here.
  expect_error(cvode(eqns, modelname = "asa_guard_d", derivMode = "reverse",
                     deriv = TRUE),
               "no forward sensitivities")

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
