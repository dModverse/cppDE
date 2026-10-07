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

per_method <- function(eq, tag, ...)
  lapply(setNames(nm = methods), function(m)
    cppODE(eq, method = m, derivMode = "reverse", modelname = paste0("rf_", tag, "_", m),
           compile = FALSE, ...))
mods <- list(rest = per_method(rest_eqns, "rest"), q = per_method(q_eqns, "q"),
             ab = per_method(ab_eqns, "ab"),
             ev = per_method(ab_eqns, "ev", events = ev))
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
      if (e_plain > 10) expect_lt(3 * e_ref, e_plain, label = paste(nm, m))
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

test_that("the checked sweep goes through events and reports its substeps", {
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
