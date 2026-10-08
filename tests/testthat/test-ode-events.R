# Time and root events on every method, against closed forms.

skip_on_cran()

# -- Models --------------------------------------------------------------------

# A dose at a parameter time on A, and a refill on x each time it falls to xc.
evt_dose <- data.frame(var = c("A", "x"), time = c("t_e", NA), value = c("dose", "dose"),
                       method = "add", root = c(NA, "xc - x"), stringsAsFactors = FALSE)
dosing <- cppODE(c(A = "-k1 * A", x = "-k * x"), events = evt_dose,
                 modelname = "evt_dosing", compile = FALSE)

# The same doubling of C, once on the root S - 14 and once at a fixed time.
eqns_ramp  <- c(S = "1", C = "0")
evt_ramp_r <- data.frame(var = "C", time = NA, value = "C * 2", method = "replace",
                         root = "S - 14", stringsAsFactors = FALSE)
evt_ramp_f <- data.frame(var = "C", time = 12, value = "C * 2", method = "replace",
                         root = NA, stringsAsFactors = FALSE)
ramp_root  <- cppODE(eqns_ramp, events = evt_ramp_r, modelname = "evt_ramp_root",
                     deriv = FALSE, compile = FALSE)
ramp_fixed <- cppODE(eqns_ramp, events = evt_ramp_f, modelname = "evt_ramp_fixed",
                     deriv = FALSE, compile = FALSE)

evt_ramp_sens <- data.frame(var = "C", time = NA, value = "d", method = "add",
                            root = "S - c", stringsAsFactors = FALSE)
ramp_sens <- cppODE(c(S = "a", C = "-b * C"), events = evt_ramp_sens, deriv = TRUE,
                    deriv2 = TRUE, modelname = "evt_ramp_sens", compile = FALSE)

evt_jump <- data.frame(var = c("S", "S", "C"), time = c(5, 7, NA),
                       value = c("20", "30", "C * 2"),
                       method = c("replace", "replace", "replace"),
                       root = c(NA, NA, "S - 14"), stringsAsFactors = FALSE)
jump_switch <- cppODE(c(S = "0", C = "0"), events = evt_jump, deriv = FALSE,
                      modelname = "evt_jump_switches_root", compile = FALSE)

# The same reset, once switched on through a root and once at the jump time.
eqns_reset <- c(S = "0 * S", C = "-b * C")
by_root <- data.frame(var = c("S", "C"), time = c("te", NA), value = c("20", "d"),
                      method = c("replace", "add"), root = c(NA, "S - 14"),
                      stringsAsFactors = FALSE)
by_time <- data.frame(var = c("S", "C"), time = c("te", "te"), value = c("20", "d"),
                      method = c("replace", "add"), root = c(NA, NA),
                      stringsAsFactors = FALSE)
reset_root <- cppODE(eqns_reset, events = by_root, deriv = TRUE, deriv2 = TRUE,
                     modelname = "evt_jump_sens_root", compile = FALSE)
reset_time <- cppODE(eqns_reset, events = by_time, deriv = TRUE, deriv2 = TRUE,
                     modelname = "evt_jump_sens_time", compile = FALSE)

methods_all <- c("bdf", "adams", "rb4", "tsit5")
evt_walls <- data.frame(var = c("v", "v"), time = c(NA, NA), value = c("-1", "-1"),
                        method = c("multiply", "multiply"),
                        root = c("x - L", "x + L"), stringsAsFactors = FALSE)
walls <- lapply(setNames(nm = methods_all), function(m)
  cppODE(c(x = "v", v = "-w^2 * x"), events = evt_walls, deriv = FALSE, method = m,
         modelname = paste0("evt_walls_", m), compile = FALSE))

native <- c(list(dosing, ramp_root, ramp_fixed, ramp_sens, jump_switch, reset_root,
                 reset_time), walls)
do.call(compile, c(unname(native), list(output = "test_ode_events", cores = test_cores())))

if (isTRUE(cvodeConfig$available)) {
  ramp_root_cv <- cvode(eqns_ramp, events = evt_ramp_r, modelname = "evt_ramp_root_cv",
                        deriv = FALSE, compile = FALSE)
  compile(ramp_root_cv, output = "test_ode_events_cvode", cores = test_cores())
}

# -- Time and root events -----------------------------------------------------

test_that("a dose at a time and a refill on a root follow their closed forms", {
  p  <- c(A = 1, x = 1, k1 = 0.1, k = 0.1, t_e = 25, xc = 0.5, dose = 0.5)
  tt <- seq(0, 50, by = 2.5)
  res <- solveODE(dosing, tt, p, abstol = 1e-10, reltol = 1e-10, roottol = 1e-10)

  # The rows at the requested times hold the state after any event there.
  rows <- vapply(tt, function(s) max(which(res$time == s)), 1L)
  A <- exp(-0.1 * tt) + ifelse(tt >= 25, 0.5 * exp(-0.1 * (tt - 25)), 0)
  expect_equal(unname(res$variable[rows, "A"]), A, tolerance = 1e-8)

  # One refill (maxroot = 1) where x falls to xc, at log(2) / k.
  t1 <- log(2) / 0.1
  x <- ifelse(tt < t1, exp(-0.1 * tt), exp(-0.1 * (tt - t1)))
  expect_equal(unname(res$variable[rows, "x"]), x, tolerance = 1e-8)
  expect_true(any(abs(res$time - t1) < 1e-8))
  # The dose time moves A afterwards by the slope of the dose's decay.
  expect_equal(unname(res$tangent[rows[tt > 25], "A", "t_e"]),
               0.1 * 0.5 * exp(-0.1 * (tt[tt > 25] - 25)), tolerance = 1e-6)
})

test_that("a root landing exactly on an output time still fires", {
  # S' = 1 puts the crossing of S - 14 onto the requested grid, where the sign
  # product of g at the sampled points is zero rather than negative. The same
  # jump written as a fixed event at that time is the reference.
  tt    <- seq(0, 20, by = 1)
  pars  <- c(S = 2, C = 4)

  # The last row at a requested time holds the post-event state, whether or
  # not the localised root inserted its own rows next to it.
  atTimes <- function(res, var)
    res$variable[vapply(tt, function(s) max(which(res$time == s)), 1L), var]

  res <- solveODE(ramp_root, tt, pars)
  ref <- solveODE(ramp_fixed, tt, pars)

  expect_equal(atTimes(res, "C"), ifelse(tt < 12, 4, 8))
  expect_equal(atTimes(res, "C"), atTimes(ref, "C"))
  expect_equal(atTimes(res, "S"), atTimes(ref, "S"))
})

test_that("both backends localise a root at the same time", {
  skip_if_not(isTRUE(cvodeConfig$available), "CVODE backend not available")
  # The root is off the grid here, so both backends have to place the firing
  # time themselves rather than inherit it from a requested time.
  tt    <- seq(0, 20, by = 1)
  pars  <- c(S = 2.5, C = 4)

  nat <- solveODE(ramp_root, tt, pars)
  cvd <- solveODE(ramp_root_cv, tt, pars)

  fired <- function(res) res$time[which(diff(res$variable[, "C"]) != 0) + 1L]
  expect_equal(fired(nat), 11.5, tolerance = 1e-6)
  expect_equal(fired(cvd), 11.5, tolerance = 1e-6)
})

test_that("a root event on the grid feeds the firing time into the sensitivities", {
  # The event fires exactly on the grid. Every sensitivity after the event
  # picks up the derivative of the firing time through the saltation term, so
  # a missed or misplaced root shows up here as well.
  pars <- c(S = 2, C = 4, a = 1, b = 0.1, c = 14, d = 3)
  res  <- solveODE(ramp_sens, seq(0, 20, by = 1), pars, abstol = 1e-10, reltol = 1e-10)

  # C(T) = C0 exp(-b T) + d exp(-b (T - t*)) for T > t*
  i  <- max(which(res$time == 16))
  E1 <- exp(-0.1 * 16)
  E2 <- exp(-0.1 * (16 - 12))

  expect_equal(unname(res$variable[i, "C"]), 4 * E1 + 3 * E2, tolerance = 1e-7)
  expect_equal(res$tangent[i, "C", "C"], E1,                tolerance = 1e-6)
  expect_equal(res$tangent[i, "C", "d"], E2,                tolerance = 1e-6)
  expect_equal(res$tangent[i, "C", "c"], 3 * E2 * 0.1,      tolerance = 1e-6)
  expect_equal(res$tangent[i, "C", "S"], -3 * E2 * 0.1,     tolerance = 1e-6)
  expect_equal(res$hessian[i, "C", "d", "c"], E2 * 0.1,     tolerance = 1e-6)
  expect_equal(res$hessian[i, "C", "c", "d"], E2 * 0.1,     tolerance = 1e-6)
})

test_that("a fixed event switches on a root condition it steps over", {
  # S jumps past the threshold instead of crossing it, so the sign-change search
  # over the continuous solution never sees it. The condition is read on both
  # sides of the jump, fires there, and does not fire again while it stays true.
  res <- solveODE(jump_switch, c(0, 4, 5, 6, 7, 8), c(S = 2, C = 4), maxroot = 2L)

  at <- function(s) unname(res$variable[max(which(res$time == s)), "C"])
  expect_equal(at(4), 4)
  expect_equal(at(5), 8)
  expect_equal(at(8), 8)
})

test_that("a reset switched on by a jump transports like a fixed one", {
  # The reset rides on the surface of the jump, so it has to transport the
  # sensitivities exactly like the same reset written as a fixed event at that
  # time, the parameter dependence of the event time included.
  pars <- c(S = 2, C = 4, b = 0.15, d = 3, te = 4)
  tt   <- seq(0, 10, by = 0.5)
  solved <- function(mod)
    solveODE(mod, tt, pars, abstol = 1e-12, reltol = 1e-12, roottol = 1e-12)

  a <- solved(reset_root)
  b <- solved(reset_time)
  expect_identical(a$time, b$time)
  expect_equal(a$variable, b$variable)
  expect_equal(a$tangent, b$tangent)
  expect_equal(a$hessian, b$hessian)
  # the saltation term of the event time is what makes this more than an identity
  expect_gt(abs(a$tangent[max(which(a$time == 8)), "C", "te"]), 0.1)
})

test_that("a root event does not fire twice on the crossing it just handled", {
  # Two elastic walls, each a root event that turns the velocity around. The
  # step restarts on the surface of the wall that just fired, where the root is
  # zero up to round-off; reading that residue as a crossing lets the mass out.
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
    res <- solveODE(walls[[m]], tt, pars, maxroot = 2L,
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
