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

# A forward-forward oracle and the forward-reverse model it checks, per method.
oracle_pair <- function(rhs, tag, meths = steppers, ...)
  list(ff = per_method(rhs, paste0("g2_", tag, "_ff_"), "forward-forward", meths, ...),
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
# A root event whose height is a parameter.
ev_root <- data.frame(var = "A", time = NA, value = "d_amt", root = "B - 0.55",
                      method = "add", stringsAsFactors = FALSE)
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

m_r  <- per_method(eqns, "g2_r_",  "reverse")
m_fr <- per_method(eqns, "g2_fr_", "forward-reverse")
m_ff <- cppODE(eqns, modelname = "g2_ff", derivMode = "forward-forward",
               compile = FALSE)
pr_forcing <- oracle_pair(eqns_u, "ev", c("bdf", "tsit5"),
                          events = ev_fixed, forcings = "u")
pr_ptime   <- oracle_pair(eqns,   "pt",   "bdf", events = ev_ptime)
pr_root    <- oracle_pair(eqns,   "rt",   events = ev_root)
pr_general <- oracle_pair(eqns,   "gen",  events = ev_general)
pr_clock   <- oracle_pair(eqns_t, "td",   events = ev_root)
pr_rtol    <- oracle_pair(eqns,   "rtol", "bdf", events = ev_clock)
do.call(compile, c(flat(list(m_r, m_fr, m_ff, pr_forcing, pr_ptime, pr_root,
                             pr_general, pr_clock, pr_rtol)),
                   output = "test_ode_reverse2", cores = 1))

# Sparse models need KLU; their test skips without it.
if (isTRUE(cppDE:::cvodeConfig$klu_available)) {
  pr_sparse <- oracle_pair(eqns, "sp", "bdf", sparse = TRUE)
  do.call(compile, c(flat(pr_sparse), output = "test_ode_reverse2_sparse",
                     cores = 1))
}

test_that("a forward-reverse solve lands on the value run's grid at any width", {
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

test_that("every method takes the value grid backwards", {
  # A forward-reverse solve chooses its steps from value arithmetic alone, so it
  # takes a value run's grid. The numbers on it are not bit-identical: a
  # corrector sums in a different order over the AD type than over double. The
  # step count agrees wherever those last bits do not straddle an acceptance
  # threshold, which is three of the four methods; adams carries twelve orders
  # of history and comes out a few steps apart. What a blocked Hessian stands on
  # is not this but that blocks ride one grid, asserted below at tolerance
  # zero.
  W <- seed_for(length(times))
  for (meth in c("bdf", "adams", "rb4", "tsit5")) {
    mr <- m_r[[meth]]
    rr <- do.call(solveODE, c(list(mr, times, pars, cotangent = W), tol))
    m  <- m_fr[[meth]]
    fr <- do.call(solveODE,
                  c(list(m, times, pars, tangent = block_dirs(5L), cotangent = W), tol))
    expect_identical(fr$time, rr$time, info = meth)
    if (identical(meth, "adams")) {
      expect_lt(abs(fr$diagnostics$accepted - rr$diagnostics$accepted) /
                  rr$diagnostics$accepted, 0.05)
    } else {
      expect_identical(fr$diagnostics$accepted, rr$diagnostics$accepted, info = meth)
    }
    expect_equal(unname(fr$variable), unname(rr$variable),
                 tolerance = 1e-9, info = meth)
    expect_equal(unname(fr$cotangent), unname(rr$cotangent),
                 tolerance = 1e-8, info = meth)
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

# ---------------------------------------------------------------------------
#  Stage 3: the features forward-reverse had never been run on.
#
#  Oracle is forward-forward on the same model. It keeps the sensitivity error
#  control and therefore its own grid, so these are tolerance comparisons, not
#  the exact ones above.
# ---------------------------------------------------------------------------

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
  N  <- length(p)

  for (m in c("bdf", "tsit5")) {
    mf <- pr_forcing$ff[[m]]
    mr <- pr_forcing$fr[[m]]
    ff <- do.call(solveODE, c(list(mf, times, p, forcings = fc), tol))
    # The jump has to be in the run, or the test proves nothing.
    expect_gt(nrow(ff$variable), length(times))
    W  <- seed_for(nrow(ff$variable))
    fr <- do.call(solveODE, c(list(mr, times, p, forcings = fc, cotangent = W), tol))
    expect_second_order(ff, fr, W, m, tol_h = 1e-4)
  }
})

test_that("a jump whose time is a parameter goes backwards at second order", {
  # The output grid carries a row at t*, and that row's TIME moves with the
  # parameter. A cotangent on it makes w.x a different functional, and then no
  # derivative agrees with a difference quotient. The comparison is therefore
  # on the user times alone.
  p  <- c(A = 1.2, B = 0.4, k1 = 0.7, k2 = 0.35, k3 = 1.1, t_ev = 1.1)
  N  <- length(p)
  tl <- list(abstol = 1e-12, reltol = 1e-12)

  mf <- pr_ptime$ff$bdf
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

test_that("a root event goes backwards at second order", {
  # A root's t* moves with theta, and the grid carries the state either side of
  # the jump at that time. A cotangent there is not the same functional at two
  # parameter values, so both rows are zeroed.
  p  <- c(A = 1.2, B = 0.4, k1 = 0.7, k2 = 0.35, k3 = 1.1, d_amt = 0.4)
  tl <- list(abstol = 1e-12, reltol = 1e-12)

  for (m in c("bdf", "adams", "rb4", "tsit5")) {
    mf <- pr_root$ff[[m]]
    mr <- pr_root$fr[[m]]
    ff <- do.call(solveODE, c(list(mf, times, p), tl))
    # The jump has to be in the run, or the test proves nothing.
    expect_gt(nrow(ff$variable), length(times))

    W <- seed_for(nrow(ff$variable))
    moving <- which(vapply(ff$time, function(x) min(abs(x - times)) > 1e-9, TRUE))
    expect_length(moving, 2L)
    W[moving, , ] <- 0

    fr <- do.call(solveODE, c(list(mr, times, p, cotangent = W), tl))
    expect_second_order(ff, fr, W, m)
  }
})

test_that("an event's root, time and value may be any expression", {
  # Every slot at once, because each leaves different terms at zero: a root
  # linear in x and blind to the clock zeroes two thirds of grad g_dot. Here g
  # is quadratic in B, reads A and carries t, the event time is nonlinear in a
  # parameter, and both heights read the state, a parameter and the clock.
  # A clock-reading height rides on roottol rather than reltol.
  p  <- c(A = 1.2, B = 0.4, k1 = 0.7, k2 = 0.35, k3 = 1.1,
          d_amt = 0.4, t_ev = 1.0)
  tl <- list(abstol = 1e-12, reltol = 1e-12, roottol = 1e-12)

  for (m in c("bdf", "adams", "rb4", "tsit5")) {
    mf <- pr_general$ff[[m]]
    mr <- pr_general$fr[[m]]
    ff <- do.call(solveODE, c(list(mf, times, p), tl))

    W <- seed_for(nrow(ff$variable))
    moving <- which(vapply(ff$time, function(x) min(abs(x - times)) > 1e-9, TRUE))
    # One row for the fixed jump and the pair a root emits, or the test proves
    # nothing: both events have to be in the run.
    expect_length(moving, 3L)
    W[moving, , ] <- 0

    fr <- do.call(solveODE, c(list(mr, times, p, cotangent = W), tl))
    expect_second_order(ff, fr, W, m)
  }
})

test_that("a right-hand side that reads the clock goes backwards too", {
  # f_e and f_a are evaluated at t*, which moves, so df/dt reaches the shift.
  # Every other model here is autonomous and leaves that term at zero.
  p  <- c(A = 1.2, B = 0.4, k1 = 0.7, k2 = 0.35, k3 = 1.1, d_amt = 0.4)
  tl <- list(abstol = 1e-12, reltol = 1e-12, roottol = 1e-12)

  for (m in c("bdf", "adams", "rb4", "tsit5")) {
    mf <- pr_clock$ff[[m]]
    mr <- pr_clock$fr[[m]]
    ff <- do.call(solveODE, c(list(mf, times, p), tl))
    W  <- seed_for(nrow(ff$variable))
    moving <- which(vapply(ff$time, function(x) min(abs(x - times)) > 1e-9, TRUE))
    expect_length(moving, 2L)
    W[moving, , ] <- 0
    fr <- do.call(solveODE, c(list(mr, times, p, cotangent = W), tl))
    expect_second_order(ff, fr, W, m)
  }
})

test_that("a clock-reading jump height rides on roottol", {
  # It reads t* itself, where a height blind to the clock only reads the state
  # there, so the localisation error reaches it undamped. Not a missing term:
  # the gap falls with roottol and floors at reltol.
  p  <- c(A = 1.2, B = 0.4, k1 = 0.7, k2 = 0.35, k3 = 1.1, d_amt = 0.4)
  mf <- pr_rtol$ff$bdf
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

test_that("a cotangent that moves with theta carries its own curvature", {
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

test_that("a sparse Jacobian goes backwards at second order", {
  skip_if_not(isTRUE(cppDE:::cvodeConfig$klu_available), "KLU not available")
  mf <- pr_sparse$ff$bdf
  mr <- pr_sparse$fr$bdf
  ff <- do.call(solveODE, c(list(mf, times, pars), tol))
  W  <- seed_for(nrow(ff$variable))
  fr <- do.call(solveODE, c(list(mr, times, pars, cotangent = W), tol))
  expect_second_order(ff, fr, W, "sparse")
})

test_that("several seed columns each carry their own second order", {
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

test_that("the batch entry carries the second order per condition", {
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
