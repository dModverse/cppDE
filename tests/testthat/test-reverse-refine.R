# The sweep checked under adjointControl(refine = TRUE), and the
# multistep sweep's operators on a grid whose tails differ step by step.

skip_on_cran()

methods <- c("bdf", "adams", "rb4", "tsit5")

rest_eqns  <- c(x = "-k * x", y = "k * x - y")
rest_pars  <- c(x = 0, y = 0, k = 0.5)
rest_times <- 0:20
q_eqns  <- c(x = "-k * x + sin(100 * time) * q")
q_pars  <- c(x = 1, k = 1, q = 0)
q_times <- seq(0, 1, by = 0.05)
ab_eqns  <- c(A = "-k1 * A + k2 * B", B = "k1 * A - k2 * B - k3 * B * B")
ab_pars  <- c(A = 1.2, B = 0.4, k1 = 0.7, k2 = 0.35, k3 = 1.1)
ab_times <- seq(0, 5, by = 0.25)
ev <- data.frame(var = "A", time = 2.1, value = "d", method = "add")
ev_pars <- c(ab_pars, d = 0.3)
# At rest behind a switch: x stays zero over 0:50 since ts = 60 lies beyond.
sw_eqns <- c(x = "-k*(x - c*piecewise(0, time - ts < 0, 1))", y = "k*x - y")
sw_pars <- c(k = 1, c = 1000, ts = 60, x = 0, y = 0)

per_method <- function(eq, tag, ...)
  lapply(setNames(nm = methods), function(m)
    cppODE(eq, method = m, derivMode = "reverse", modelname = paste0("rf_", tag, "_", m),
           compile = FALSE, ...))
mods <- list(rest = per_method(rest_eqns, "rest"), q = per_method(q_eqns, "q"),
             ab = per_method(ab_eqns, "ab"),
             ev = per_method(ab_eqns, "ev", events = ev),
             sw = per_method(sw_eqns, "sw"))
fwd <- list(rest = cppODE(rest_eqns, modelname = "rf_fwd_rest", compile = FALSE),
            q = cppODE(q_eqns, modelname = "rf_fwd_q", compile = FALSE),
            ab = cppODE(ab_eqns, modelname = "rf_fwd_ab", compile = FALSE),
            ev = cppODE(ab_eqns, events = ev, modelname = "rf_fwd_ev", compile = FALSE))
fr_bdf <- cppODE(ab_eqns, method = "bdf", derivMode = "forward-reverse",
                 modelname = "rf_fr_bdf", compile = FALSE)
do.call(compile, c(unname(c(unlist(mods, recursive = FALSE), fwd, list(fr_bdf))),
                   list(output = "test_reverse_refine", cores = test_cores())))

cases <- list(rest = list(p = rest_pars, t = rest_times),
              q = list(p = q_pars, t = q_times),
              ab = list(p = ab_pars, t = ab_times),
              ev = list(p = ev_pars, t = ab_times))

seed_for <- function(nm) {
  rows <- nrow(solveODE(fwd[[nm]], cases[[nm]]$t, cases[[nm]]$p)$variable)
  set.seed(1)
  nx <- length(attr(fwd[[nm]], "variables"))
  array(rnorm(rows * nx), c(rows, nx, 1))
}

# The gradient a forward solve at a tight tolerance gives, and the error of g
# against it per component relative to rtol |g| + rtol 1e-3 max |g|.
reference <- function(nm, W) {
  r <- solveODE(fwd[[nm]], cases[[nm]]$t, cases[[nm]]$p, abstol = 1e-13, reltol = 1e-11)
  apply(r$tangent * as.vector(W), 3, sum)
}
err_tau <- function(g, ref, rtol) {
  ref <- ref[names(g)]
  max(abs(g - ref) / (rtol * abs(ref) + rtol * 1e-3 * max(abs(ref))))
}

test_that("refine = FALSE leaves the sweep as it was", {
  W <- seed_for("ab")
  for (m in methods) {
    a <- solveODE(mods$ab[[m]], ab_times, ab_pars, cotangent = W, abstol = 1e-8, reltol = 1e-6)
    b <- solveODE(mods$ab[[m]], ab_times, ab_pars, cotangent = W, abstol = 1e-8, reltol = 1e-6,
                  adjoint = adjointControl(refine = FALSE))
    expect_identical(a$cotangent, b$cotangent, info = m)
  }
})

test_that("the checked sweep is accurate where the value grid is not", {
  # The value grid steps over a trajectory at rest and over a fast forcing the
  # state does not see; under gradtol the checked sweep refines the steps whose
  # share of the gradient misses the tolerance and comes within a few tau.
  o <- list(abstol = 1e-6, reltol = 1e-4)
  for (nm in c("rest", "q", "ab")) {
    W <- seed_for(nm); ref <- reference(nm, W)
    for (m in methods) {
      plain <- do.call(solveODE, c(list(mods[[nm]][[m]], cases[[nm]]$t, cases[[nm]]$p,
                                        cotangent = W, adjoint = adjointControl()), o))
      ga <- 1e-4 * 1e-3 * max(abs(ref))
      ref_r <- do.call(solveODE, c(list(mods[[nm]][[m]], cases[[nm]]$t, cases[[nm]]$p,
                                        cotangent = W,
                                        adjoint = adjointControl(refine = TRUE, gradtol = ga)), o))
      e_plain <- err_tau(plain$cotangent[, 1], ref, 1e-4)
      e_ref <- err_tau(ref_r$cotangent[, 1], ref, 1e-4)
      # The fast forcing sums local errors over many periods, as any local
      # control does, the more so on tsit5, which has no embedded estimate.
      expect_lt(e_ref, if (nm == "q") 100 else 10, label = paste(nm, m))
      if (e_plain > 20) expect_lt(3 * e_ref, e_plain, label = paste(nm, m))
    }
  }
})

test_that("the checked sweep subdivides the steps the value grid takes too long", {
  # At rest the value grid takes few long steps, whose share of the gradient
  # misses the tolerance on every method.
  W <- seed_for("rest")
  for (m in methods) {
    r <- solveODE(mods$rest[[m]], rest_times, rest_pars, cotangent = W, abstol = 1e-6,
                  reltol = 1e-4, adjoint = adjointControl(refine = TRUE, trace = TRUE))
    expect_gt(max(r$adjoint$substeps), 1L, label = m)
  }
})

test_that("the checked sweep goes through events and reports the substeps", {
  W <- seed_for("ev"); ref <- reference("ev", W)
  for (m in methods) {
    r <- solveODE(mods$ev[[m]], ab_times, ev_pars, cotangent = W, abstol = 1e-8, reltol = 1e-6,
                  adjoint = adjointControl(refine = TRUE, gradtol = 1e-9 * max(abs(ref)),
                                           trace = TRUE))
    expect_lt(err_tau(r$cotangent[, 1], ref, 1e-6), 10, label = m)
    expect_identical(length(r$adjoint$substeps), length(r$adjoint$time))
    expect_true(all(r$adjoint$substeps >= 1L))
  }
})

test_that("refine is refused where it does not apply", {
  W <- seed_for("ab")
  expect_error(solveODE(fr_bdf, ab_times, ab_pars, cotangent = W,
                        adjoint = adjointControl(refine = TRUE)), "derivMode")
  expect_error(adjointControl(refine = NA), "TRUE or FALSE")
})

# The value correction of the refined sweep, per seed column.
value_of <- function(res, W) apply(W * as.vector(res$variable), 3, sum)
seeds2 <- function(rows, nx = 2L) {
  set.seed(3)
  array(rnorm(rows * nx * 2), c(rows, nx, 2), dimnames = list(NULL, NULL, c("a", "b")))
}

test_that("on a linear model the corrected value is the cotangent times the start", {
  # Every accepted map is linear in the state, so value plus correction is the
  # functional of the differentiated trajectory, x0' w_x0, to rounding. Coarse
  # steps put observations inside them; gradtol makes the flow halve.
  p <- c(x = 1, y = 0.5, k = 0.5)
  o <- list(abstol = 1e-5, reltol = 1e-3)
  for (m in methods) for (gt in list(NULL, 1e-9)) {
    W <- seeds2(length(rest_times))
    r <- do.call(solveODE, c(list(mods$rest[[m]], rest_times, p, cotangent = W,
                                  adjoint = adjointControl(refine = TRUE, gradtol = gt)), o))
    expect_identical(names(r$correction), c("a", "b"))
    expect_gt(max(abs(r$correction)), 0)
    lhs <- value_of(r, W) + r$correction
    rhs <- colSums(c(1, 0.5) * r$cotangent[c("x", "y"), ])
    expect_equal(lhs, rhs, tolerance = 1e-10, info = paste(m, is.null(gt)))
  }
})

test_that("the correction moves a coarse value toward the converged one", {
  # Where every step is substepped, the corrected value is the value of the
  # finer trajectory the gradient belongs to.
  o <- list(abstol = 1e-5, reltol = 1e-3)
  for (nm in c("ab", "ev")) {
    W <- seed_for(nm)
    ref <- value_of(solveODE(fwd[[nm]], cases[[nm]]$t, cases[[nm]]$p,
                             abstol = 1e-13, reltol = 1e-11), W)
    for (m in c("bdf", "tsit5")) {
      r <- do.call(solveODE, c(list(mods[[nm]][[m]], cases[[nm]]$t, cases[[nm]]$p,
                                    cotangent = W, adjoint = adjointControl(refine = TRUE)), o))
      v <- value_of(r, W)
      expect_lt(abs(v + r$correction - ref), 0.2 * abs(v - ref), label = paste(nm, m))
    }
  }
})

test_that("the correction is the same from a store, a batch and a prepared batch", {
  p <- c(x = 1, y = 0.5, k = 0.5)
  o <- list(abstol = 1e-5, reltol = 1e-3)
  ctl <- adjointControl(refine = TRUE)
  W <- seeds2(length(rest_times))
  for (m in c("bdf", "rb4")) {
    mm <- mods$rest[[m]]
    one <- do.call(solveODE, c(list(mm, rest_times, p, cotangent = W, adjoint = ctl), o))
    plain <- do.call(solveODE, c(list(mm, rest_times, p, cotangent = W), o))
    expect_null(plain$correction)

    st <- do.call(solveODE, c(list(mm, rest_times, p, keepStore = TRUE), o))
    expect_null(st$correction)
    fs <- do.call(solveODE, c(list(mm, rest_times, p, cotangent = W, store = st$store,
                                   adjoint = ctl), o))
    expect_identical(fs$correction, one$correction, info = m)

    cond <- list(times = rest_times, parms = p, cotangent = W, adjoint = ctl)
    bt <- do.call(solveODEBatch, c(list(mm, conditions = list(cond, cond)), o))
    expect_identical(bt[[2]]$correction, one$correction, info = m)

    # A handle prepared with one seed column takes two.
    h <- do.call(prepareBatch, c(list(mm, conditions = list(
      modifyList(cond, list(cotangent = W[, , 1, drop = FALSE])))), o))
    got <- solveBatch(h, cotangent = list(W))[[1]]
    expect_identical(got$correction, one$correction, info = m)
    expect_identical(got$cotangent, one$cotangent, info = m)
  }
})

test_that("a trajectory behind a switch gets a correction within the tolerance", {
  # x = e^-t and y = (y0 + t) e^-t while the switch lies beyond the horizon. The
  # corrected value is that of the substepped flow, accurate to the refinement test.
  tt <- 0:50
  p <- c(k = 1, c = 1000, ts = 60, x = 1, y = 0.5)
  set.seed(1)
  W <- array(rnorm(length(tt) * 2), c(length(tt), 2, 1))
  ref <- sum(W[, 1, 1] * exp(-tt) + W[, 2, 1] * (0.5 + tt) * exp(-tt))
  for (m in methods) {
    r <- solveODE(mods$sw[[m]], tt, p, cotangent = W, abstol = 1e-8, reltol = 1e-5,
                  adjoint = adjointControl(refine = TRUE, gradtol = 1e-10))
    v <- value_of(r, W)
    expect_true(is.finite(r$correction), info = m)
    expect_lt(abs(v + r$correction - ref), 10 * 1e-5 * abs(ref), label = m)
  }
})

test_that("the multistep sweep takes each step's own tail", {
  # Operators kept for one step must not serve a later step of the same size
  # whose history rescales by another factor. attr(times, "hmax"), internal
  # (start times, step bounds), forces alternating bounds.
  W <- seed_for("ab")
  tt <- ab_times
  attr(tt, "hmax") <- list(c(0, 1, 2, 3, 4), c(0.02, 0.03, 0.02, 0.03, 0.02))
  r <- solveODE(mods$ab$bdf, tt, ab_pars, cotangent = W, abstol = 1e-6, reltol = 1e-4)
  f <- solveODE(fr_bdf, tt, ab_pars, cotangent = W, abstol = 1e-6, reltol = 1e-4)
  gf <- apply(f$tangent * as.vector(W), 3, sum)
  g <- r$cotangent[names(gf), 1]
  expect_lt(max(abs(g - gf)) / max(abs(gf)), 1e-3)
})

test_that("a trajectory at rest behind a switch gets the exact gradient", {
  # The state's error estimate is zero, so the grid of a value-only run grows
  # without bound and the adjoint on it is off by more than the adjoint size. The
  # checked sweep, the gradient under gradtol, takes those steps in substeps.
  # With k = 1 the tangents are dx/dx0 = e^-t, dy/dx0 = t e^-t, dy/dy0 = e^-t,
  # and zero in k, c and ts while the switch lies beyond the horizon.
  tt <- 0:50
  set.seed(1)
  W <- array(rnorm(length(tt) * 2 * 2), c(length(tt), 2, 2))
  ref <- rbind(x = colSums(W[, 1, ] * exp(-tt) + W[, 2, ] * tt * exp(-tt)),
               y = colSums(W[, 2, ] * exp(-tt)), k = 0, c = 0, ts = 0)
  o <- list(abstol = 1e-12, reltol = 1e-10)
  for (m in methods) {
    rv <- do.call(solveODE, c(list(mods$sw[[m]], tt, sw_pars, cotangent = W,
                                   adjoint = adjointControl(refine = TRUE,
                                                            gradtol = 1e-13)), o))
    expect_equal(unname(rv$cotangent[rownames(ref), ]), unname(ref),
                 tolerance = 1e-6, info = m)
  }
})
