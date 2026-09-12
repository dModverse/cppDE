# Forward over reverse: the grid it differentiates, and the features it had
# never been run on. Stages 2 to 4 of dev/closed-adjoint-plan.md's successor.
#
# The first block is exact and not a tolerance test. A forward-reverse solve
# takes its step sequence from value arithmetic alone, so it lands on the grid
# a plain value run takes, whatever the tangent count is and whatever the
# tangents contain. That is what makes a Hessian assembled from blocks of
# directions one matrix rather than several.
#
# The feature block's oracle is forward-forward, which keeps the sensitivity
# error control and therefore its own grid, so those comparisons are O(tol).

skip_on_cran()

eqns <- c(A = "-k1 * A + k2 * B",
          B = "k1 * A - k2 * B - k3 * B * B")
pars  <- c(A = 1.2, B = 0.4, k1 = 0.7, k2 = 0.35, k3 = 1.1)
times <- seq(0, 5, length.out = 21)
tol   <- list(abstol = 1e-10, reltol = 1e-10)
n_phi <- length(pars)

seed_for <- function(n_t, n_x = 2L, n_seed = 1L, seed = 1L) {
  set.seed(seed)
  array(rnorm(n_t * n_x * n_seed), c(n_t, n_x, n_seed))
}

# B columns of the identity: the directions a block of a chunked Hessian takes.
block_dirs <- function(B) { S <- matrix(0, n_phi, B); diag(S) <- 1; S }

test_that("a forward-reverse solve lands on the value run's grid at any width", {
  mr <- do.call(cppODE, list(eqns, modelname = "g2_r", derivMode = "reverse"))
  W  <- seed_for(length(times))
  rr <- do.call(solveODE, c(list(mr, times, pars, seed = W), tol))

  for (B in c(1L, 3L, 5L)) {
    m <- cppODE(eqns, modelname = paste0("g2_fr", B),
                derivMode = "forward-reverse", nStack = B)
    fr <- do.call(solveODE,
                  c(list(m, times, pars, sens1ini = block_dirs(B), seed = W), tol))
    expect_identical(fr$time, rr$time, info = paste("B =", B))
    expect_identical(unname(fr$variable), unname(rr$variable), info = paste("B =", B))
    expect_identical(unname(fr$adjoint), unname(rr$adjoint), info = paste("B =", B))
  }
})

test_that("the grid does not depend on what the tangents contain", {
  m <- cppODE(eqns, modelname = "g2_content",
              derivMode = "forward-reverse", nStack = 3L)
  W <- seed_for(length(times))
  set.seed(7)
  dirs <- list(identity = block_dirs(3L),
               random   = matrix(rnorm(n_phi * 3L), n_phi, 3L),
               zero     = matrix(0, n_phi, 3L))
  ref <- NULL
  for (nm in names(dirs)) {
    fr <- do.call(solveODE,
                  c(list(m, times, pars, sens1ini = dirs[[nm]], seed = W), tol))
    if (is.null(ref)) ref <- fr
    expect_identical(fr$time, ref$time, info = nm)
    expect_identical(unname(fr$variable), unname(ref$variable), info = nm)
    expect_identical(unname(fr$adjoint), unname(ref$adjoint), info = nm)
  }
})

test_that("every method takes the value grid backwards", {
  # The grid is exactly the value run's on all four. The numbers on it are not
  # bit-identical, and cannot be: a corrector runs a fused expression over
  # double and a copy plus two axpys over the AD type, which is a different
  # summation order. What matters for a blocked Hessian is the grid, and the
  # exactness across widths is asserted above.
  W <- seed_for(length(times))
  for (meth in c("bdf", "adams", "rb4", "tsit5")) {
    mr <- cppODE(eqns, modelname = paste0("g2m_r_", meth),
                 derivMode = "reverse", method = meth)
    rr <- do.call(solveODE, c(list(mr, times, pars, seed = W), tol))
    m  <- cppODE(eqns, modelname = paste0("g2m_fr_", meth),
                 derivMode = "forward-reverse", method = meth, nStack = 5L)
    fr <- do.call(solveODE,
                  c(list(m, times, pars, sens1ini = block_dirs(5L), seed = W), tol))
    expect_identical(fr$time, rr$time, info = meth)
    expect_identical(fr$diagnostics$accepted, rr$diagnostics$accepted, info = meth)
    expect_equal(unname(fr$variable), unname(rr$variable),
                 tolerance = 1e-9, info = meth)
    expect_equal(unname(fr$adjoint), unname(rr$adjoint),
                 tolerance = 1e-8, info = meth)
  }
})

test_that("a Hessian assembled from blocks is the one a single pass gives", {
  m5 <- cppODE(eqns, modelname = "g2_h5",
               derivMode = "forward-reverse", nStack = 5L)
  m2 <- cppODE(eqns, modelname = "g2_h2",
               derivMode = "forward-reverse", nStack = 2L)
  W  <- seed_for(length(times))
  one <- do.call(solveODE,
                 c(list(m5, times, pars, sens1ini = block_dirs(5L), seed = W), tol))

  # Three blocks of two, the last one short, tiled into the full matrix.
  H <- matrix(NA_real_, n_phi, n_phi)
  for (start in c(1L, 3L, 5L)) {
    idx <- seq(start, min(start + 1L, n_phi))
    S <- matrix(0, n_phi, 2L)
    for (j in seq_along(idx)) S[idx[j], j] <- 1
    blk <- do.call(solveODE,
                   c(list(m2, times, pars, sens1ini = S, seed = W), tol))
    expect_identical(unname(blk$adjoint), unname(one$adjoint),
                     info = as.character(start))
    H[, idx] <- blk$adjoint2[, seq_along(idx), 1L]
  }
  expect_identical(H, unname(one$adjoint2[, , 1L]))
  # Symmetric by construction on one grid, though not bit for bit: each column
  # is a different sequence of the same arithmetic.
  expect_equal(H, t(H), tolerance = 1e-9)
})

# ---------------------------------------------------------------------------
#  Stage 3: the features forward-reverse had never been run on.
#
#  Oracle is forward-forward on the same model. It keeps the sensitivity error
#  control and therefore its own grid, so these are tolerance comparisons, not
#  the exact ones above.
# ---------------------------------------------------------------------------

# w' S contracted over times and states: the gradient the sweep returns.
contract <- function(sens1, W) {
  vapply(seq_len(dim(W)[3]),
         function(k) apply(sens1 * as.vector(W[, , k]), 3, sum),
         numeric(dim(sens1)[3]))
}

# The Hessian of w' x from the forward second derivatives.
hess_forward <- function(res, W) {
  ns <- dim(res$sens2)[3]
  outer(seq_len(ns), seq_len(ns),
        Vectorize(function(a, b) sum(as.vector(W) * res$sens2[, , a, b])))
}

expect_second_order <- function(ff, fr, W, info, tol_g = 1e-5, tol_h = 1e-5) {
  ref <- contract(ff$sens1, W)[, 1]
  expect_equal(unname(fr$adjoint[names(ref), 1]), unname(ref),
               tolerance = tol_g, info = info)
  expect_equal(unname(fr$adjoint2[, , 1]), hess_forward(ff, W),
               tolerance = tol_h, info = info)
}

test_that("forcings and a jump go backwards at second order", {
  # What works: a forcing, and a jump whose height is a parameter. What does
  # not: a jump whose TIME is a parameter, the case below.
  eq <- c(A = "-k1 * A + k2 * B + u",
          B = "k1 * A - k2 * B - k3 * B * B")
  # Off the output grid, so the jump shows up as its own pair of rows.
  ev <- data.frame(var = "A", time = 1.1, value = "d_amt", method = "add",
                   stringsAsFactors = FALSE)
  p  <- c(A = 1.2, B = 0.4, k1 = 0.7, k2 = 0.35, k3 = 1.1, d_amt = 0.4)
  fc <- list(u = data.frame(time = c(0, 2, 5), value = c(0.1, 0.25, 0.05)))
  N  <- length(p)

  for (m in c("bdf", "tsit5")) {
    mf <- cppODE(eq, events = ev, forcings = "u", method = m,
                 modelname = paste0("g2_ev_ff_", m),
                 derivMode = "forward-forward", nStack = N)
    mr <- cppODE(eq, events = ev, forcings = "u", method = m,
                 modelname = paste0("g2_ev_fr_", m),
                 derivMode = "forward-reverse", nStack = N)
    ff <- do.call(solveODE, c(list(mf, times, p, forcings = fc), tol))
    # The jump has to be in the run, or the test proves nothing.
    expect_gt(nrow(ff$variable), length(times))
    W  <- seed_for(nrow(ff$variable))
    fr <- do.call(solveODE, c(list(mr, times, p, forcings = fc, seed = W), tol))
    expect_second_order(ff, fr, W, m, tol_h = 1e-4)
  }
})

test_that("a jump whose time is a parameter is first order only", {
  skip(paste("known gap: the saltation adjoint carries the event time at first",
             "order but not in its tangents. The gradient falls with the",
             "tolerance, the Hessian sits at a floor. Measured: a constant or",
             "parameter-valued jump height at a fixed time 1e-8, a root event",
             "3e-2, a parameter event time 8e-2."))
})

test_that("a sparse Jacobian goes backwards at second order", {
  skip_if_not(isTRUE(cppDE:::cvodeConfig$klu_available), "KLU not available")
  mf <- cppODE(eqns, modelname = "g2_sp_ff", sparse = TRUE,
               derivMode = "forward-forward", nStack = n_phi)
  mr <- cppODE(eqns, modelname = "g2_sp_fr", sparse = TRUE,
               derivMode = "forward-reverse", nStack = n_phi)
  ff <- do.call(solveODE, c(list(mf, times, pars), tol))
  W  <- seed_for(nrow(ff$variable))
  fr <- do.call(solveODE, c(list(mr, times, pars, seed = W), tol))
  expect_second_order(ff, fr, W, "sparse")
})

test_that("several seed columns each carry their own second order", {
  mf <- cppODE(eqns, modelname = "g2_ns_ff",
               derivMode = "forward-forward", nStack = n_phi)
  mr <- cppODE(eqns, modelname = "g2_ns_fr",
               derivMode = "forward-reverse", nStack = n_phi)
  ff <- do.call(solveODE, c(list(mf, times, pars), tol))
  W  <- seed_for(nrow(ff$variable), n_seed = 3L)
  fr <- do.call(solveODE, c(list(mr, times, pars, seed = W), tol))

  expect_identical(dim(fr$adjoint2), c(n_phi, n_phi, 3L))
  for (k in seq_len(3L)) {
    Wk  <- W[, , k, drop = FALSE]
    ref <- contract(ff$sens1, Wk)[, 1]
    expect_equal(unname(fr$adjoint[names(ref), k]), unname(ref),
                 tolerance = 1e-5, info = as.character(k))
    expect_equal(unname(fr$adjoint2[, , k]), hess_forward(ff, Wk),
                 tolerance = 1e-5, info = as.character(k))
  }
})

test_that("a non-identity sens1ini reads the Hessian along its own directions", {
  # adjoint2[, k] is the gradient's derivative along direction k, so with mixed
  # directions it is the full Hessian times S, not S' H S.
  set.seed(11)
  S  <- matrix(rnorm(n_phi * 3L), n_phi, 3L)
  mf <- cppODE(eqns, modelname = "g2_rp_ff",
               derivMode = "forward-forward", nStack = n_phi)
  mr <- cppODE(eqns, modelname = "g2_rp_fr",
               derivMode = "forward-reverse", nStack = n_phi)
  ff <- do.call(solveODE, c(list(mf, times, pars), tol))
  W  <- seed_for(nrow(ff$variable))
  fr <- do.call(solveODE, c(list(mr, times, pars, sens1ini = S, seed = W), tol))

  expect_identical(dim(fr$adjoint2), c(n_phi, 3L, 1L))
  expect_equal(unname(fr$adjoint2[, , 1]), hess_forward(ff, W) %*% S,
               tolerance = 1e-5)
})

test_that("a store is refused under second order rather than answered wrongly", {
  # A checkpoint's tangents live in the arena of the solve that filled it, so a
  # store handed to a later solve would give exact values and a wrong Hessian.
  m <- cppODE(eqns, modelname = "g2_store",
              derivMode = "forward-reverse", nStack = n_phi)
  S <- block_dirs(n_phi)
  expect_error(
    do.call(solveODE, c(list(m, times, pars, sens1ini = S, keepStore = TRUE), tol)),
    "forward-reverse")
})

test_that("the batch entry carries the second order per condition", {
  m <- cppODE(eqns, modelname = "g2_batch",
              derivMode = "forward-reverse", nStack = n_phi)
  S <- block_dirs(n_phi)
  p2 <- pars; p2["k1"] <- 0.9
  one <- do.call(solveODE, c(list(m, times, pars, sens1ini = S,
                                  seed = seed_for(length(times))), tol))
  two <- do.call(solveODE, c(list(m, times, p2, sens1ini = S,
                                  seed = seed_for(length(times))), tol))

  bt <- do.call(solveODEBatch,
                c(list(m, conditions = list(
                    list(times = times, parms = pars, sens1ini = S,
                         seed = seed_for(length(times))),
                    list(times = times, parms = p2, sens1ini = S,
                         seed = seed_for(length(times))))), tol))
  expect_length(bt, 2L)
  expect_identical(unname(bt[[1]]$adjoint2), unname(one$adjoint2))
  expect_identical(unname(bt[[2]]$adjoint2), unname(two$adjoint2))
})
