# The reverse mode through the R seam, first and second order, on one model.
# Events, forcings and sparse Jacobians are in test-reverse-events.R, the CVODES
# adjoint in test-reverse-cvode.R.
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

# Models nest in lists; compile() takes them flat.
flat <- function(x) if (is.list(x)) do.call(c, lapply(unname(x), flat)) else list(x)

# One uncompiled cppODE model per method, named by method.
per_method <- function(rhs, prefix, ms = methods, ...)
  lapply(setNames(nm = ms), function(m)
    cppODE(rhs, method = m, modelname = paste0(prefix, m), compile = FALSE, ...))

# The forward oracle on bdf, and per method the reverse, forward-forward and
# forward-reverse models, linked into one shared object.
models <- list(
  fwd = cppODE(eqns, modelname = "rev_m_f", deriv = TRUE, compile = FALSE),
  rev = per_method(eqns, "rev_m_r_", derivMode = "reverse"),
  ff  = per_method(eqns, "rev2_ff_", derivMode = "forward-forward"),
  fr  = per_method(eqns, "rev2_fr_", derivMode = "forward-reverse"))
do.call(compile, c(flat(models),
                   list(output = "test_reverse", cores = test_cores())))

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
  mf <- models$fwd
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
  mf <- models$fwd
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

test_that("the reverse mode goes through the batch entry", {
  mf <- models$fwd
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
  mf <- models$fwd
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

test_that("every method supports the reverse mode", {
  fwd <- do.call(solveODE, c(list(models$fwd, times, pars), tol))
  W   <- seed_for(fwd)
  ref <- contract(fwd$tangent, W)[, 1]
  for (m in methods) {
    rv  <- do.call(solveODE, c(list(models$rev[[m]], times, pars, cotangent = W), tol))
    expect_equal(unname(rv$cotangent[names(ref), 1]), unname(ref),
                 tolerance = 1e-5, info = m)
  }
})

test_that("$adjoint reports the grid the sweep ran on", {
  mr <- models$rev$bdf
  mv <- models$fwd

  val <- do.call(solveODE, c(list(mv, times, pars), tol))
  W   <- seed_for(val, n_seed = 2L)

  plain <- do.call(solveODE, c(list(mr, times, pars, cotangent = W), tol))
  expect_null(plain$adjoint)

  rv <- do.call(solveODE,
                c(list(mr, times, pars, cotangent = W, adjoint = adjointControl(trace = TRUE)), tol))
  G  <- rv$adjoint
  expect_s3_class(G, "cppDEadjoint")
  expect_true(all(c("time", "h", "eta", "lambda") %in% names(G)))

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
  W <- seed_for(solveODE(models$fwd, times, pars))
  for (m in methods) {
    mr <- models$rev[[m]]
    rv  <- solveODE(mr, times, pars, abstol = 1e-8, reltol = 1e-6,
                    cotangent = W, adjoint = adjointControl(trace = TRUE))
    G <- rv$adjoint

    expect_gt(max(abs(G$eta)), 1e-12, label = paste("max |eta| on", m))
    expect_gt(max(abs(G$lambda)), 1e-3, label = paste("max |lambda| on", m))
    # Loosening the tolerance a hundredfold has to raise the indicator: it is
    # an error estimate, not a property of the trajectory.
    rv2 <- solveODE(mr, times, pars, abstol = 1e-6, reltol = 1e-4,
                    cotangent = W, adjoint = adjointControl(trace = TRUE))
    expect_gt(sum(abs(rv2$adjoint$eta)), sum(abs(G$eta)), label = m)
  }
})

test_that("malformed adjoint controls are an error, not a silent no-op", {
  mr <- models$rev$bdf
  mv <- models$fwd
  W  <- seed_for(solveODE(mv, times, pars))

  expect_error(adjointControl(gradtol = 0), "positive")
  expect_error(adjointControl(trace = NA), "TRUE or FALSE")
  expect_error(solveODE(mr, times, pars, cotangent = W, adjoint = list(trace = TRUE)),
               "adjointControl")
  expect_error(solveODE(mv, times, pars, adjoint = adjointControl(trace = TRUE)),
               "backward pass")
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
  # This is the whole safety of the scheme. A store holds the run it was made
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
  mf <- models$fwd
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

test_that("a prepared batch sweeps the stores it is handed", {
  mr <- models$rev$bdf
  p2 <- pars * c(1, 1, 1.1, 0.9, 1.05)
  mk <- function(p) list(times = times, parms = p, keepStore = TRUE)
  fw <- do.call(solveODEBatch, c(list(mr, conditions = list(mk(pars), mk(p2))), tol))
  W  <- lapply(fw, seed_for)
  conds <- lapply(1:2, function(i)
    list(times = times, parms = list(pars, p2)[[i]], cotangent = W[[i]],
         store = fw[[i]]$store))
  ref <- do.call(solveODEBatch, c(list(mr, conditions = conds), tol))
  h <- do.call(prepareBatch, c(list(mr, conditions = conds), tol))

  # New numbers, cotangents and stores, the prepared ones swapped round.
  got <- solveBatch(h, parms = list(p2, pars), cotangent = rev(W),
                    store = list(fw[[2]]$store, fw[[1]]$store))
  expect_equal(got[[1]]$cotangent, ref[[2]]$cotangent, tolerance = 0)
  expect_equal(got[[2]]$cotangent, ref[[1]]$cotangent, tolerance = 0)

  # A store from the other point is refused, not replayed.
  expect_error(solveBatch(h, parms = list(pars, p2),
                          store = list(fw[[2]]$store, fw[[1]]$store)))

  # Without a store the sweep integrates, on the grid of its own forward pass.
  free <- solveBatch(h, store = list(NULL, NULL))
  expect_equal(free[[1]]$cotangent, ref[[1]]$cotangent, tolerance = 1e-10)
})
