# Test all integration methods on problems with known analytical solutions.

skip_on_cran()

# -- Shared setup --------------------------------------------------------------

eqns_decay <- c(A = "-k1 * A", B = "k1 * A - k2 * B")
times <- seq(0, 50, length.out = 200)
pars  <- c(A = 1, B = 0, k1 = 0.1, k2 = 0.2)
#
exact_A <- function(t, A0, k1) A0 * exp(-k1 * t)
exact_B <- function(t, A0, k1, k2) A0 * k1 / (k2 - k1) * (exp(-k1 * t) - exp(-k2 * t))

methods_all <- c("bdf", "adams", "rb4", "tsit5")

# -- Basic solver output structure --------------------------------------------

test_that("solveODE returns correct structure for all methods", {
  for (m in methods_all) {
    mod <- cppODE(eqns_decay, method = m, modelname = paste0("struct_", m))
    res <- solveODE(mod, times, pars)

    expect_type(res$time, "double")
    expect_true(is.matrix(res$variable))
    expect_equal(ncol(res$variable), 2)
    expect_equal(colnames(res$variable), c("A", "B"))
    expect_equal(length(res$time), nrow(res$variable))
  }
})

# -- Accuracy against analytical solution --------------------------------------

test_that("all methods match analytical solution (decay system)", {
  for (m in methods_all) {
    mod <- cppODE(eqns_decay, method = m, modelname = paste0("accuracy_", m))
    res <- solveODE(mod, times, pars, abstol = 1e-10, reltol = 1e-10)

    t_out <- res$time
    A_num <- res$variable[, "A"]
    B_num <- res$variable[, "B"]

    A_exact <- exact_A(t_out, pars["A"], pars["k1"])
    B_exact <- exact_B(t_out, pars["A"], pars["k1"], pars["k2"])

    expect_equal(A_num, A_exact, tolerance = 1e-6, label = paste(m, "A"))
    expect_equal(B_num, B_exact, tolerance = 1e-6, label = paste(m, "B"))
  }
})

# -- First-order sensitivities via AD vs finite differences --------------------

test_that("first-order sensitivities are correct for all methods", {
  eps <- 1e-5

  for (m in methods_all) {
    mod    <- cppODE(eqns_decay, method = m, deriv = TRUE,
                     modelname = paste0("sens1_", m))
    res    <- solveODE(mod, times, pars, abstol = 1e-10, reltol = 1e-10)

    # Check AD sensitivity of A w.r.t. k1 against finite difference
    p_hi      <- pars; p_hi["k1"] <- pars["k1"] + eps
    p_lo      <- pars; p_lo["k1"] <- pars["k1"] - eps
    mod_noad  <- cppODE(eqns_decay, method = m, deriv = FALSE,
                        modelname = paste0("sens1_fd_", m))
    res_hi    <- solveODE(mod_noad, times, p_hi, abstol = 1e-12, reltol = 1e-12)
    res_lo    <- solveODE(mod_noad, times, p_lo, abstol = 1e-12, reltol = 1e-12)

    # Match time grids (both should be identical since same output times)
    fd_dA_dk1 <- (res_hi$variable[, "A"] - res_lo$variable[, "A"]) / (2 * eps)

    sens_names <- attr(mod, "dimNames")$sens
    k1_idx <- which(sens_names == "k1")
    ad_dA_dk1 <- res$sens1[, 1, k1_idx]  # state 1 (A), param k1_idx

    expect_equal(ad_dA_dk1, fd_dA_dk1, tolerance = 1e-3,
                 label = paste(m, "dA/dk1"))
  }
})

# -- Second-order sensitivities ------------------------------------------------

test_that("second-order sensitivities are finite for stiff methods", {
  stiff_methods <- c("bdf", "rb4")

  for (m in stiff_methods) {
    mod <- cppODE(eqns_decay, method = m, deriv = TRUE, deriv2 = TRUE,
                  modelname = paste0("sens2_", m))
    res <- solveODE(mod, times, pars, abstol = 1e-10, reltol = 1e-10)

    expect_true(!is.null(res$sens2), label = paste(m, "sens2 exists"))
    expect_true(all(is.finite(res$sens2)), label = paste(m, "sens2 finite"))
  }
})

# -- Time-triggered events -----------------------------------------------------

test_that("time-triggered dose event works", {
  eqns <- c(A = "-k1 * A")
  evt  <- data.frame(var = "A", time = "t_e", value = "dose",
                     method = "add", root = NA, stringsAsFactors = FALSE)
  pars_ev <- c(A = 1, k1 = 0.1, t_e = 25, dose = 0.5)

  mod <- cppODE(eqns, events = evt, modelname = "event_time")
  res <- solveODE(mod, seq(0, 50, length.out = 200), pars_ev)

  # After event at t=25, A should jump up
  idx_before <- max(which(res$time < 25))
  idx_after  <- min(which(res$time >= 25))
  expect_gt(res$variable[idx_after, "A"], res$variable[idx_before, "A"])
})

# -- Root-triggered events -----------------------------------------------------

test_that("root-triggered event fires correctly", {
  eqns <- c(x = "-k * x")
  evt  <- data.frame(var = "x", time = NA, value = "dose",
                     method = "add", root = "xc - x",
                     stringsAsFactors = FALSE)
  pars_root <- c(x = 1, k = 0.1, xc = 0.5, dose = 0.5)

  mod <- cppODE(eqns, events = evt, modelname = "event_root")
  res <- solveODE(mod, seq(0, 100, length.out = 500), pars_root)

  # x decays below xc, then gets dose added -> should see multiple oscillations
  # Check that x goes back up after crossing xc at least once
  x_vals <- res$variable[, "x"]
  crossings <- sum(diff(x_vals > 0.5) != 0)
  expect_gt(crossings, 0, label = "root event triggers at least once")
})

test_that("a root landing exactly on an output time still fires", {
  # S' = 1 puts the crossing of S - 14 onto the requested grid, where the sign
  # product of g at the sampled points is zero rather than negative. The same
  # jump written as a fixed event at that time is the reference.
  eqns  <- c(S = "1", C = "0")
  evt_r <- data.frame(var = "C", time = NA, value = "C * 2", method = "replace",
                      root = "S - 14", stringsAsFactors = FALSE)
  evt_f <- data.frame(var = "C", time = 12, value = "C * 2", method = "replace",
                      root = NA, stringsAsFactors = FALSE)
  tt    <- seq(0, 20, by = 1)
  pars  <- c(S = 2, C = 4)

  # The last row at a requested time carries the post-event state, whether or
  # not the localised root inserted its own rows next to it.
  atTimes <- function(res, var)
    res$variable[vapply(tt, function(s) max(which(res$time == s)), 1L), var]

  res <- solveODE(cppODE(eqns, events = evt_r, modelname = "ongrid_root",
                         deriv = FALSE), tt, pars)
  ref <- solveODE(cppODE(eqns, events = evt_f, modelname = "ongrid_fixed",
                         deriv = FALSE), tt, pars)

  expect_equal(atTimes(res, "C"), ifelse(tt < 12, 4, 8))
  expect_equal(atTimes(res, "C"), atTimes(ref, "C"))
  expect_equal(atTimes(res, "S"), atTimes(ref, "S"))
})

test_that("both backends localise a root at the same time", {
  skip_if_not(isTRUE(cvodeConfig$available), "CVODE backend not available")
  # The root is off the grid here, so both backends have to place the firing
  # time themselves rather than inherit it from a requested time.
  eqns  <- c(S = "1", C = "0")
  evt_r <- data.frame(var = "C", time = NA, value = "C * 2", method = "replace",
                      root = "S - 14", stringsAsFactors = FALSE)
  tt    <- seq(0, 20, by = 1)
  pars  <- c(S = 2.5, C = 4)

  nat <- solveODE(cppODE(eqns, events = evt_r, modelname = "offgrid_nat",
                         deriv = FALSE), tt, pars)
  cvd <- solveODE(cvode(eqns, events = evt_r, modelname = "offgrid_cv",
                        deriv = FALSE), tt, pars)

  fired <- function(res) res$time[which(diff(res$variable[, "C"]) != 0) + 1L]
  expect_equal(fired(nat), 11.5, tolerance = 1e-6)
  expect_equal(fired(cvd), 11.5, tolerance = 1e-6)
})

test_that("a root event on the grid carries the firing time into the sensitivities", {
  # S' = a fires the event at t* = (c - S0)/a = 12, again exactly on the grid.
  # Every sensitivity of C after the event picks up dt*/dtheta through the
  # saltation term, so a missed or misplaced root shows up here as well.
  evt <- data.frame(var = "C", time = NA, value = "d", method = "add",
                    root = "S - c", stringsAsFactors = FALSE)
  mod <- cppODE(c(S = "a", C = "-b * C"), events = evt, deriv = TRUE,
                deriv2 = TRUE, modelname = "ongrid_sens")
  pars <- c(S = 2, C = 4, a = 1, b = 0.1, c = 14, d = 3)
  res  <- solveODE(mod, seq(0, 20, by = 1), pars, abstol = 1e-10, reltol = 1e-10)

  # C(T) = C0 exp(-b T) + d exp(-b (T - t*)) for T > t*
  i  <- max(which(res$time == 16))
  E1 <- exp(-0.1 * 16)
  E2 <- exp(-0.1 * (16 - 12))

  expect_equal(unname(res$variable[i, "C"]), 4 * E1 + 3 * E2, tolerance = 1e-7)
  expect_equal(res$sens1[i, "C", "C"], E1,                tolerance = 1e-6)
  expect_equal(res$sens1[i, "C", "d"], E2,                tolerance = 1e-6)
  expect_equal(res$sens1[i, "C", "c"], 3 * E2 * 0.1,      tolerance = 1e-6)
  expect_equal(res$sens1[i, "C", "S"], -3 * E2 * 0.1,     tolerance = 1e-6)
  expect_equal(res$sens2[i, "C", "d", "c"], E2 * 0.1,     tolerance = 1e-6)
  expect_equal(res$sens2[i, "C", "c", "d"], E2 * 0.1,     tolerance = 1e-6)
})

test_that("a fixed event switches on a root condition it steps over", {
  # S jumps past the threshold instead of crossing it, so the sign-change search
  # over the continuous solution never sees it. The condition is read on both
  # sides of the jump, fires there, and does not fire again while it stays true.
  evt <- data.frame(var = c("S", "S", "C"), time = c(5, 7, NA),
                    value = c("20", "30", "C * 2"),
                    method = c("replace", "replace", "replace"),
                    root = c(NA, NA, "S - 14"), stringsAsFactors = FALSE)
  mod <- cppODE(c(S = "0", C = "0"), events = evt, deriv = FALSE,
                modelname = "jump_switches_root")
  res <- solveODE(mod, c(0, 4, 5, 6, 7, 8), c(S = 2, C = 4), maxroot = 2L)

  at <- function(s) unname(res$variable[max(which(res$time == s)), "C"])
  expect_equal(at(4), 4)
  expect_equal(at(5), 8)
  expect_equal(at(8), 8)
})

test_that("a reset switched on by a jump transports like a fixed one", {
  # The reset rides on the surface of the jump, so it has to carry the
  # sensitivities exactly like the same reset written as a fixed event at that
  # time, the parameter dependence of the event time included.
  eqns    <- c(S = "0 * S", C = "-b * C")
  by_root <- data.frame(var = c("S", "C"), time = c("te", NA), value = c("20", "d"),
                        method = c("replace", "add"), root = c(NA, "S - 14"),
                        stringsAsFactors = FALSE)
  by_time <- data.frame(var = c("S", "C"), time = c("te", "te"), value = c("20", "d"),
                        method = c("replace", "add"), root = c(NA, NA),
                        stringsAsFactors = FALSE)
  pars <- c(S = 2, C = 4, b = 0.15, d = 3, te = 4)
  tt   <- seq(0, 10, by = 0.5)
  solved <- function(evt, name)
    solveODE(cppODE(eqns, events = evt, deriv = TRUE, deriv2 = TRUE,
                    modelname = name), tt, pars,
             abstol = 1e-12, reltol = 1e-12, roottol = 1e-12)

  a <- solved(by_root, "jump_sens_root")
  b <- solved(by_time, "jump_sens_time")
  expect_identical(a$time, b$time)
  expect_equal(a$variable, b$variable)
  expect_equal(a$sens1, b$sens1)
  expect_equal(a$sens2, b$sens2)
  # the saltation term of the event time is what makes this more than an identity
  expect_gt(abs(a$sens1[max(which(a$time == 8)), "C", "te"]), 0.1)
})

test_that("a root event does not fire twice on the crossing it just handled", {
  # Two elastic walls, each a root event that turns the velocity around. The
  # step restarts on the surface of the wall that just fired, where the root is
  # zero up to round-off; reading that residue as a crossing lets the mass out.
  evt <- data.frame(var = c("v", "v"), time = c(NA, NA), value = c("-1", "-1"),
                    method = c("multiply", "multiply"),
                    root = c("x - L", "x + L"), stringsAsFactors = FALSE)
  pars <- c(x = 0.2, v = 1.2, w = 1, L = 0.6)
  tt   <- seq(0, 4.4, by = 0.1)

  # amplitude, first wall contact and the flight from one wall to the other
  amp    <- sqrt(pars[["x"]]^2 + (pars[["v"]] / pars[["w"]])^2)
  speed  <- sqrt(pars[["v"]]^2 + (pars[["w"]] * pars[["x"]])^2 -
                 (pars[["w"]] * pars[["L"]])^2)
  first  <- 2 * atan((pars[["v"]] - speed) /
                     (pars[["w"]] * (pars[["L"]] + pars[["x"]]))) / pars[["w"]]
  flight <- 2 * asin(pars[["L"]] / amp) / pars[["w"]]

  for (m in methods_all) {
    mod <- cppODE(c(x = "v", v = "-w^2 * x"), events = evt, method = m,
                  deriv = FALSE, modelname = paste0("walls_", m))
    res <- solveODE(mod, tt, pars, maxroot = 2L,
                    abstol = 1e-12, reltol = 1e-12, roottol = 1e-12)

    expect_lte(max(abs(res$variable[, "x"])), pars[["L"]] + 1e-9,
               label = paste(m, "stays inside the walls"))
    bounces <- res$time[which(diff(sign(res$variable[, "v"])) != 0) + 1L]
    expect_equal(bounces, first + (0:3) * flight, tolerance = 1e-6,
                 label = paste(m, "bounce times"))
    energy <- res$variable[, "v"]^2 + (pars[["w"]] * res$variable[, "x"])^2
    expect_lt(diff(range(energy)), 1e-6, label = paste(m, "energy spread"))
  }
})

# -- Diagnostics ---------------------------------------------------------------

test_that("diagnostics() returns solver statistics", {
  mod <- cppODE(eqns_decay, modelname = "diag_test")
  res <- solveODE(mod, times, pars)
  d   <- diagnostics(res)

  expect_true(is.list(d))
  expect_true(d$accepted > 0)
  expect_true(d$fevals > 0)
  expect_equal(d$return_code, 0)  # success
})

# -- Fixed parameters ----------------------------------------------------------

test_that("fixed parameters are excluded from sensitivities", {
  mod <- cppODE(eqns_decay, deriv = TRUE, fixed = "k2",
                modelname = "fixed_test")
  res <- solveODE(mod, times, pars)

  sens_names <- attr(mod, "dimNames")$sens
  expect_false("k2" %in% sens_names)
  expect_true("k1" %in% sens_names)
})

# -- Constant-only math calls in the generated Jacobian ------------------------

# A 10^x term differentiates to a math call whose arguments are all literals.
# Codegen emits every math call as cppde::<fn>, so this compiles only because
# cppde_dual_math.hpp carries arithmetic-type overloads next to the AD ones.
test_that("a 10^x term compiles and differentiates correctly", {
  t10  <- seq(0, 1, length.out = 25)
  p10  <- c(x = 0.3, k = 0.7)
  ln10 <- log(10)

  # closed form: 10^(-x(t)) = 10^(-x0) + k*ln(10)*t
  u <- 10^(-p10[["x"]]) + p10[["k"]] * ln10 * t10

  mod <- cppODE(c(x = "-k * 10^x"), deriv = TRUE, deriv2 = TRUE, nStack = 2,
                modelname = "pow10")
  res <- solveODE(mod, t10, p10, abstol = 1e-12, reltol = 1e-12)

  expect_equal(as.numeric(res$variable[, "x"]), -log10(u), tolerance = 1e-8)
  expect_equal(as.numeric(res$sens1[, "x", "x"]),
               10^(-p10[["x"]]) / u, tolerance = 1e-6)
  expect_equal(as.numeric(res$sens1[, "x", "k"]), -t10 / u, tolerance = 1e-6)
  expect_equal(as.numeric(res$sens2[, "x", "k", "k"]),
               ln10 * t10^2 / u^2, tolerance = 1e-6)
})

# The generated right-hand side indexes x[] and params[] and calls std::pow, so
# a state or parameter carrying one of those names has to be substituted before
# it can be read as part of the surrounding code.
test_that("state and parameter names that are C++ tokens compile and solve", {
  tt <- seq(0, 2, 0.5)

  mod <- cppODE(c(default = "-std * default + operator",
                  int     = "std * default - int * 10^0.5"),
                modelname = "cxx_tokens_ode")
  ref <- cppODE(c(a = "-k * a + b0",
                  b = "k * a - b * 10^0.5"),
                modelname = "cxx_tokens_ref")

  res <- solveODE(mod, tt, c(default = 1, int = 0, std = 0.7, operator = 0.2))
  exp <- solveODE(ref, tt, c(a = 1, b = 0, k = 0.7, b0 = 0.2))

  expect_equal(unname(res$variable), unname(exp$variable), tolerance = 1e-10)
  expect_equal(unname(res$sens1), unname(exp$sens1), tolerance = 1e-10)
})

test_that("a Python keyword as a symbol name is rejected", {
  expect_error(cppODE(c(y = "-class * y"), modelname = "py_kw_ode"),
               "Python keyword used as a symbol name: 'class'")
  expect_error(cppODE(c(lambda = "-k * lambda"), modelname = "py_kw_state"),
               "'lambda'")
  expect_error(cppODE(c(y = "-k * y"), forcings = "global",
                      modelname = "py_kw_forcing"),
               "'global'")
})
