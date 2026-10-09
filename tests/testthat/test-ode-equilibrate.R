# Steady-state equilibration (rootfunc = "equilibrate") across methods, and
# warm starts.

skip_on_cran()

# -- Shared setup --------------------------------------------------------------
# The system has a closed-form steady state, so the equilibrated result can be
# asserted against it rather than against a long transient solve.

rhs <- c(
  R  = "k_act - k_deact * R",
  A  = "-k1 * A * R + k2 * pA",
  pA = " k1 * A * R - k2 * pA"
)

pars <- c(R = 1, A = 1, pA = 0, k_act = 0.1, k_deact = 0.7,
          k1 = 0.1, k2 = 0.05)
times <- seq(0, 1e3, length.out = 500)

stiff_methods <- c("bdf", "rb4")
all_methods   <- c("bdf", "adams", "rb4", "tsit5")

# Sensitivities on the stiff methods only, which the tangent tests use; all
# models are generated uncompiled and linked into one shared object.
eq_mod <- lapply(setNames(nm = all_methods), function(m)
  cppODE(rhs, rootfunc = "equilibrate", method = m, deriv = m %in% stiff_methods,
         modelname = paste0("eq_", m), compile = FALSE))
# Two states whose steady states lie 15 orders apart.
rhs_scale <- c(small = "ks - d * small", big = "kb - d * big")
eq_scale  <- cppODE(rhs_scale, rootfunc = "equilibrate", deriv = FALSE,
                    modelname = "eq_scale", compile = FALSE)
do.call(compile, c(unname(eq_mod), list(eq_scale),
                   output = "test_ode_equilibrate", cores = test_cores()))

# The closed-form steady state and its derivatives in every parameter and
# initial value, the latter entering through the conserved total.
ss_expr <- expression(
  R  = k_act / k_deact,
  A  = (A + pA) * k2 / (k2 + k1 * k_act / k_deact),
  pA = (A + pA) * k1 * k_act / k_deact / (k2 + k1 * k_act / k_deact)
)
ss_jac <- function(p) {
  out <- t(vapply(ss_expr, function(e)
    vapply(names(p), function(v) eval(D(e, v), as.list(p)), 0), numeric(length(p))))
  dimnames(out) <- list(variable = names(ss_expr), sens = names(p))
  out
}
ss_of <- function(p) vapply(ss_expr, eval, 0, as.list(p))

# -- Steady state and termination ----------------------------------------------

test_that("every method stops early at the closed-form steady state", {
  for (m in all_methods) {
    res <- solveODE(eq_mod[[m]], times, pars, roottol = 1e-06)
    last <- length(res$time)
    expect_equal(res$variable[last, ], ss_of(pars), tolerance = 1e-4, info = m)
    expect_lt(res$time[last], 500)
    expect_gt(res$time[last], 1)
  }
})

test_that("the tangent at the steady state is the derivative of the closed form", {
  for (m in stiff_methods) {
    res <- solveODE(eq_mod[[m]], times, pars, roottol = 1e-08)
    tan <- res$tangent[length(res$time), , ]
    expect_equal(tan[, names(pars)], ss_jac(pars), tolerance = 1e-4, info = m)
  }
})

test_that("a tighter roottol integrates longer", {
  res_loose <- solveODE(eq_mod$bdf, times, pars, roottol = 1e-02)
  res_tight <- solveODE(eq_mod$bdf, times, pars, roottol = 1e-08)
  expect_gt(max(res_tight$time), max(res_loose$time))
})

# -- Warm start ----------------------------------------------------------------

test_that("a warm start reaches the new steady state in fewer steps than a cold one", {
  for (m in stiff_methods) {
    res1 <- solveODE(eq_mod[[m]], times, pars, roottol = 1e-06)
    last <- length(res1$time)
    pars2 <- replace(pars, colnames(res1$variable), res1$variable[last, ])
    pars2[c("k1", "k2")] <- pars[c("k1", "k2")] * c(1.01, 0.99)

    warm <- solveODE(eq_mod[[m]], times, pars2, tangent = res1$tangent[last, , ],
                     roottol = 1e-06)
    cold <- solveODE(eq_mod[[m]], times, pars2, roottol = 1e-06)
    expect_lt(diagnostics(warm)$accepted, diagnostics(cold)$accepted, label = m)
    expect_equal(warm$variable[length(warm$time), ], ss_of(pars2), tolerance = 1e-4,
                 info = m)
  }
})

test_that("a start at the steady state ends at once, unless the tangent must settle", {
  pars_ss <- replace(pars, names(ss_expr), ss_of(pars))
  states_only <- solveODE(eq_mod$adams, times, pars_ss, roottol = 1e-04)
  expect_lt(diagnostics(states_only)$accepted, 5)

  with_tangent <- solveODE(eq_mod$bdf, times, pars_ss, roottol = 1e-04)
  expect_gt(diagnostics(with_tangent)$accepted, 5)
  expect_equal(with_tangent$variable[length(with_tangent$time), ], ss_of(pars),
               tolerance = 1e-4)
})

# -- Relative criterion ---------------------------------------------------------

test_that("equilibrate judges each state relative to its size", {
  # CVODE checks at output times, so both run on a grid
  tt <- c(0, 10^seq(0, 6, by = 0.25))
  p  <- c(small = 0, big = 0, ks = 1e-10, kb = 1e5, d = 0.1)
  check <- function(res) {
    last <- length(res$time)
    expect_lt(res$time[last], 1e6)
    expect_equal(unname(res$variable[last, "small"]), 1e-9, tolerance = 1e-6)
    expect_equal(unname(res$variable[last, "big"]), 1e6, tolerance = 1e-6)
  }
  check(solveODE(eq_scale, tt, p, roottol = 1e-8, abstol = 1e-16, reltol = 1e-10))

  skip_if_not(isTRUE(cvodeConfig$available), "CVODE backend not available")
  cv <- cvode(rhs_scale, rootfunc = "equilibrate", modelname = "eq_scale_cv")
  check(solveODE(cv, tt, p, roottol = 1e-8, abstol = 1e-16, reltol = 1e-10))
})
