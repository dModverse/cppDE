# Forward over reverse: the grid it differentiates. The features it runs on are
# in test-ode-reverse2-events.R.
#
# The grid claims are exact and not tolerance tests. A forward-reverse solve
# takes its step sequence from value arithmetic alone, so its grid does not
# depend on how many tangents it follows or what they contain. That is what
# makes a Hessian assembled from blocks of directions one matrix rather than
# several.

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

# ---------------------------------------------------------------------------
#  Models, generated uncompiled and linked into one shared object. Tests that
#  need the same model share it.
# ---------------------------------------------------------------------------

steppers <- c("bdf", "adams", "rb4", "tsit5")

# One model per method, named <prefix><method>.
per_method <- function(rhs, prefix, derivMode, meths = steppers, ...)
  lapply(setNames(nm = meths), function(meth)
    cppODE(rhs, method = meth, derivMode = derivMode,
           modelname = paste0(prefix, meth), compile = FALSE, ...))

# Models nest in lists; compile() takes them flat.
flat <- function(x) if (is.list(x)) do.call(c, lapply(unname(x), flat)) else list(x)

m_r  <- per_method(eqns, "g2_r_",  "reverse")
m_fr <- per_method(eqns, "g2_fr_", "forward-reverse")
m_ff <- cppODE(eqns, modelname = "g2_ff", derivMode = "forward-forward",
               compile = FALSE)
do.call(compile, c(flat(list(m_r, m_fr, m_ff)),
                   output = "test_ode_reverse2", cores = test_cores()))

test_that("a forward-reverse solve lands on the reverse run's grid at any width", {
  mr <- m_r$bdf
  m  <- m_fr$bdf
  W  <- seed_for(length(times))
  rr <- do.call(solveODE, c(list(mr, times, pars, cotangent = W), tol))

  for (B in c(1L, 3L, 5L)) {
    fr <- do.call(solveODE,
                  c(list(m, times, pars, tangent = block_dirs(B), cotangent = W), tol))
    # The grid is the exact claim. The numbers on it are not bit-identical:
    # a corrector sums in a different order over the AD type than over double.
    expect_identical(fr$time, rr$time, info = paste("B =", B))
    expect_equal(unname(fr$variable), unname(rr$variable),
                 tolerance = 1e-9, info = paste("B =", B))
    expect_equal(unname(fr$cotangent), unname(rr$cotangent),
                 tolerance = 1e-8, info = paste("B =", B))
  }
})

test_that("the grid does not depend on what the tangents contain", {
  m <- m_fr$bdf
  W <- seed_for(length(times))
  set.seed(7)
  dirs <- list(identity = block_dirs(3L),
               random   = matrix(rnorm(n_phi * 3L), n_phi, 3L),
               zero     = matrix(0, n_phi, 3L))
  ref <- NULL
  for (nm in names(dirs)) {
    fr <- do.call(solveODE,
                  c(list(m, times, pars, tangent = dirs[[nm]], cotangent = W), tol))
    if (is.null(ref)) ref <- fr
    expect_identical(fr$time, ref$time, info = nm)
    expect_identical(unname(fr$variable), unname(ref$variable), info = nm)
    expect_identical(unname(fr$cotangent), unname(ref$cotangent), info = nm)
  }
})

test_that("every method takes the reverse run's grid backwards", {
  # Both runs control their steps on the state alone, so their grids differ by a
  # few steps at most; that blocks share one grid is asserted exactly below.
  W <- seed_for(length(times))
  for (meth in c("bdf", "adams", "rb4", "tsit5")) {
    mr <- m_r[[meth]]
    rr <- do.call(solveODE, c(list(mr, times, pars, cotangent = W), tol))
    m  <- m_fr[[meth]]
    fr <- do.call(solveODE,
                  c(list(m, times, pars, tangent = block_dirs(5L), cotangent = W), tol))
    expect_identical(fr$time, rr$time, info = meth)
    expect_lt(abs(fr$diagnostics$accepted - rr$diagnostics$accepted) /
                rr$diagnostics$accepted, 0.15, label = meth)
    expect_equal(unname(fr$variable), unname(rr$variable),
                 tolerance = 1e-5, info = meth)
    expect_equal(unname(fr$cotangent), unname(rr$cotangent),
                 tolerance = 1e-4, info = meth)
  }
})

test_that("blocks ride one grid on every method", {
  # The property a chunked second-order objective stands on, across the four
  # steppers and with the runtime width differing from block to block.
  W <- seed_for(length(times))
  for (meth in c("bdf", "adams", "rb4", "tsit5")) {
    m <- m_fr[[meth]]
    one <- do.call(solveODE, c(list(m, times, pars, cotangent = W), tol))
    H <- matrix(0, n_phi, n_phi)
    for (start in c(1L, 3L, 5L)) {
      idx <- seq(start, min(start + 1L, n_phi))
      S <- matrix(0, n_phi, length(idx))
      for (j in seq_along(idx)) S[idx[j], j] <- 1
      blk <- do.call(solveODE, c(list(m, times, pars, tangent = S, cotangent = W), tol))
      # One grid is the exact claim; the block width changes the order the same
      # arithmetic runs in, so the numbers agree to rounding.
      expect_identical(blk$time, one$time, info = meth)
      expect_equal(unname(blk$cotangent), unname(one$cotangent),
                   tolerance = 1e-12, info = meth)
      H[, idx] <- blk$curvature[, seq_along(idx), 1L]
    }
    expect_equal(H, unname(one$curvature[, , 1L]), tolerance = 1e-12, info = meth)
  }
})

test_that("a Hessian assembled from blocks is the one a single pass gives", {
  m  <- m_fr$bdf
  W  <- seed_for(length(times))
  one <- do.call(solveODE,
                 c(list(m, times, pars, tangent = block_dirs(5L), cotangent = W), tol))

  # Three blocks of two, the last one short, tiled into the full matrix.
  H <- matrix(NA_real_, n_phi, n_phi)
  for (start in c(1L, 3L, 5L)) {
    idx <- seq(start, min(start + 1L, n_phi))
    S <- matrix(0, n_phi, 2L)
    for (j in seq_along(idx)) S[idx[j], j] <- 1
    blk <- do.call(solveODE,
                   c(list(m, times, pars, tangent = S, cotangent = W), tol))
    expect_identical(blk$time, one$time, info = as.character(start))
    expect_equal(unname(blk$cotangent), unname(one$cotangent),
                 tolerance = 1e-12, info = as.character(start))
    H[, idx] <- blk$curvature[, seq_along(idx), 1L]
  }
  expect_equal(H, unname(one$curvature[, , 1L]), tolerance = 1e-12)
  # Symmetric by construction on one grid, though not bit for bit: each column
  # is a different sequence of the same arithmetic.
  expect_equal(H, t(H), tolerance = 1e-9)
})

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

test_that("every direction runs on the heap and blocks ride one grid", {
  # Tangent storage is heap-allocated at the width ncol(tangent) gives. What
  # has to hold is that a Hessian assembled from blocks of directions, each a
  # different runtime width, is the one a single pass gives.
  #
  # The trap is silent: a heap dual with no width cannot arm, and a tangent read
  # off an unarmed one returns the out-of-bounds zero, which leaves the gradient
  # bit-identical and the Hessian wrong. Measure the curvature, not just the
  # cotangent.
  mh <- m_fr$bdf
  W  <- seed_for(length(times))
  one <- do.call(solveODE, c(list(mh, times, pars, cotangent = W), tol))

  H <- matrix(0, n_phi, n_phi)
  for (start in c(1L, 3L, 5L)) {
    idx <- seq(start, min(start + 1L, n_phi))
    S <- matrix(0, n_phi, length(idx))
    for (j in seq_along(idx)) S[idx[j], j] <- 1
    blk <- do.call(solveODE, c(list(mh, times, pars, tangent = S, cotangent = W), tol))
    expect_identical(blk$time, one$time, info = as.character(start))
    expect_equal(unname(blk$cotangent), unname(one$cotangent),
                 tolerance = 1e-12, info = as.character(start))
    H[, idx] <- blk$curvature[, seq_along(idx), 1L]
  }
  expect_equal(H, unname(one$curvature[, , 1L]), tolerance = 1e-12)
})

test_that("a cotangent that moves with theta has its own curvature", {
  # A cotangent handed down from above the ODE depends on theta too, and
  # `curvature` is where that enters. Without it the Hessian loses the cross
  # term sum_r (dw_r/dtheta_b)(dx_r/dtheta_a), which is not small.
  mf <- m_ff
  mr <- m_fr$bdf
  ff <- do.call(solveODE, c(list(mf, times, pars), tol))
  nt <- nrow(ff$variable)

  # A cotangent that is itself a function of the state: w_ri = c_i x_1(t_r), so
  # dw_ri/dtheta_j = c_i S_1j(t_r) and the cross term cannot vanish.
  set.seed(7); cvec <- rnorm(2)
  W <- array(0, c(nt, 2L, 1L))
  W[, , 1] <- outer(rep(1, nt), cvec) * as.vector(ff$variable[, 1])
  Wdot <- array(0, c(nt, 2L, 1L, n_phi))
  for (j in seq_len(n_phi)) Wdot[, , 1, j] <- outer(ff$tangent[, 1, j], cvec)

  fr <- do.call(solveODE,
                c(list(mr, times, pars, cotangent = W, curvature = Wdot), tol))

  both <- outer(seq_len(n_phi), seq_len(n_phi), Vectorize(function(a, b)
    sum(Wdot[, , 1, b] * ff$tangent[, , a]) + sum(W[, , 1] * ff$hessian[, , a, b])))
  expect_equal(unname(fr$curvature[, , 1]), both, tolerance = 1e-6)

  # And the term it adds is worth having: drop the curvature and the answer is
  # the one that ignores the cotangent's own motion.
  fr0 <- do.call(solveODE, c(list(mr, times, pars, cotangent = W), tol))
  only <- outer(seq_len(n_phi), seq_len(n_phi), Vectorize(function(a, b)
    sum(W[, , 1] * ff$hessian[, , a, b])))
  expect_equal(unname(fr0$curvature[, , 1]), only, tolerance = 1e-6)
  expect_gt(max(abs(both - only)), 1)
})

test_that("several seed columns each get their own second order", {
  mf <- m_ff
  mr <- m_fr$bdf
  ff <- do.call(solveODE, c(list(mf, times, pars), tol))
  W  <- seed_for(nrow(ff$variable), n_seed = 3L)
  fr <- do.call(solveODE, c(list(mr, times, pars, cotangent = W), tol))

  expect_identical(dim(fr$curvature), c(n_phi, n_phi, 3L))
  for (k in seq_len(3L)) {
    Wk  <- W[, , k, drop = FALSE]
    ref <- contract(ff$tangent, Wk)[, 1]
    expect_equal(unname(fr$cotangent[names(ref), k]), unname(ref),
                 tolerance = 1e-5, info = as.character(k))
    expect_equal(unname(fr$curvature[, , k]), hess_forward(ff, Wk),
                 tolerance = 1e-5, info = as.character(k))
  }
})

test_that("a non-identity tangent reads the Hessian along its own directions", {
  # curvature[, k] is the gradient's derivative along direction k, so with mixed
  # directions it is the full Hessian times S, not S' H S.
  set.seed(11)
  S  <- matrix(rnorm(n_phi * 3L), n_phi, 3L)
  mf <- m_ff
  mr <- m_fr$bdf
  ff <- do.call(solveODE, c(list(mf, times, pars), tol))
  W  <- seed_for(nrow(ff$variable))
  fr <- do.call(solveODE, c(list(mr, times, pars, tangent = S, cotangent = W), tol))

  expect_identical(dim(fr$curvature), c(n_phi, 3L, 1L))
  expect_equal(unname(fr$curvature[, , 1]), hess_forward(ff, W) %*% S,
               tolerance = 1e-5)
})

test_that("a store is refused under second order rather than answered wrongly", {
  # A checkpoint's tangents live in the arena of the solve that filled it, so a
  # store handed to a later solve would give exact values and a wrong Hessian.
  m <- m_fr$bdf
  S <- block_dirs(n_phi)
  expect_error(
    do.call(solveODE, c(list(m, times, pars, tangent = S, keepStore = TRUE), tol)),
    "forward-reverse")
})

test_that("the batch entry returns the second order per condition", {
  m <- m_fr$bdf
  S <- block_dirs(n_phi)
  p2 <- pars; p2["k1"] <- 0.9
  one <- do.call(solveODE, c(list(m, times, pars, tangent = S,
                                  cotangent = seed_for(length(times))), tol))
  two <- do.call(solveODE, c(list(m, times, p2, tangent = S,
                                  cotangent = seed_for(length(times))), tol))

  bt <- do.call(solveODEBatch,
                c(list(m, conditions = list(
                    list(times = times, parms = pars, tangent = S,
                         cotangent = seed_for(length(times))),
                    list(times = times, parms = p2, tangent = S,
                         cotangent = seed_for(length(times))))), tol))
  expect_length(bt, 2L)
  expect_identical(unname(bt[[1]]$curvature), unname(one$curvature))
  expect_identical(unname(bt[[2]]$curvature), unname(two$curvature))
})
