# A fixed-time event fires on the window of its solve, from the first time of
# the grid up to but not including the last. Outside it the event adds no row
# and moves neither end of the integration; at the first time it is applied
# before the solve starts, at the last time not at all. Two consecutive solves
# over [0, 5] and [5, 8] then reproduce one over [0, 8], which is what multiple
# shooting relies on. The rule lives in inst/include/cppde/cppde_event_window.hpp.

skip_on_cran()

# FitzHugh-Nagumo with "add 1 to V at t = 5", once at a fixed time and once at a
# time and with a dose that are parameters.
fhn  <- c(V = "c * (V - V^3/3 + R)", R = "-(V - a + b*R)/c")
pars <- c(V = -1, R = 1, a = 0.2, b = 0.2, c = 3)
ppar <- c(pars, t_e = 5, dose = 1)
tol  <- list(abstol = 1e-10, reltol = 1e-10)
methods <- c("bdf", "adams", "rb4", "tsit5")

ev_fixed <- data.frame(var = "V", time = 5, value = 1, method = "add",
                       root = NA, stringsAsFactors = FALSE)
ev_param <- data.frame(var = "V", time = "t_e", value = "dose", method = "add",
                       root = NA, stringsAsFactors = FALSE)

# One uncompiled model per method, named by method. Without includeTimeZero
# the grid starts where `times` does, which is what the window is about.
per_method <- function(prefix, ms = methods, ...)
  lapply(setNames(nm = ms), function(m)
    cppODE(fhn, method = m, modelname = paste0(prefix, m), compile = FALSE, ...))

models <- list(
  ev   = per_method("ewin_ev_", events = ev_fixed, deriv = TRUE,
                    includeTimeZero = FALSE),
  none = per_method("ewin_none_", deriv = TRUE, includeTimeZero = FALSE),
  par  = per_method("ewin_par_", events = ev_param, deriv = TRUE,
                    includeTimeZero = FALSE),
  rev  = per_method("ewin_rev_", c("bdf", "rb4"), events = ev_param,
                    derivMode = "reverse", includeTimeZero = FALSE),
  rev0 = per_method("ewin_rev0_", c("bdf", "rb4"), derivMode = "reverse",
                    includeTimeZero = FALSE),
  ff   = per_method("ewin_ff_", c("bdf", "rb4"), events = ev_param,
                    derivMode = "forward-forward", includeTimeZero = FALSE),
  fr   = per_method("ewin_fr_", c("bdf", "rb4"), events = ev_param,
                    derivMode = "forward-reverse", includeTimeZero = FALSE),
  tz   = per_method("ewin_tz_", "bdf", events = ev_fixed, deriv = TRUE),
  ctl  = per_method("ewin_ctl_", "rb4", events = ev_fixed, deriv = TRUE,
                    useDenseOutput = FALSE, includeTimeZero = FALSE))
do.call(compile, c(unname(unlist(models, recursive = FALSE)),
                   list(output = "test_event_window", cores = test_cores())))

if (isTRUE(cvodeConfig$available)) {
  cv_par <- cvode(fhn, events = ev_param, modelname = "ewin_cv_par",
                  deriv = TRUE, includeTimeZero = FALSE, compile = FALSE)
  compile(cv_par, output = "test_event_window_cvode", cores = test_cores())
}

solve <- function(model, times, parms, ...)
  do.call(solveODE, c(list(model, times, parms, ...), tol))

# The tangent of a solve over [t1, t2] started from the state `first` ends on,
# chained through the tangent of `first` at its last row: the first-order
# sensitivities of one solve over both windows.
chain_tangent <- function(first, second, row) {
  T1 <- first$tangent[nrow(first$variable), , ]
  T2 <- second$tangent[row, , ]
  st <- colnames(first$variable)
  out <- T2[, st, drop = FALSE] %*% T1[, colnames(T2), drop = FALSE]
  par <- setdiff(colnames(T2), st)
  out[, par] <- out[, par] + T2[, par]
  dimnames(out) <- dimnames(T2)
  out
}

# w' S contracted over times and states, the quantity a reverse sweep returns.
contract <- function(tangent, W) apply(tangent * as.vector(W), 3, sum)

# w' H contracted the same way, the curvature a forward-reverse sweep returns.
hess_forward <- function(res, W) {
  ns <- dim(res$hessian)[3]
  outer(seq_len(ns), seq_len(ns),
        Vectorize(function(a, b) sum(as.vector(W) * res$hessian[, , a, b])))
}

test_that("an event before the window is ignored", {
  for (m in methods) {
    r   <- solve(models$ev[[m]], 6:8, pars)
    ref <- solve(models$none[[m]], 6:8, pars)
    expect_identical(r$time, c(6, 7, 8), info = m)
    expect_identical(r$variable[1, ], pars[c("V", "R")], info = m)
    expect_identical(r$variable, ref$variable, info = m)
    expect_identical(r$tangent, ref$tangent, info = m)
  }
  # The controlled loop, which lands on every requested time.
  r <- solve(models$ctl$rb4, 6:8, pars)
  expect_identical(r$time, c(6, 7, 8))
  expect_identical(r$variable[1, ], pars[c("V", "R")])
})

test_that("an event after the window adds no row", {
  tt <- seq(0, 4.5, by = 0.5)
  for (m in methods) {
    r   <- solve(models$ev[[m]], tt, pars)
    ref <- solve(models$none[[m]], tt, pars)
    expect_identical(r$time, tt, info = m)
    expect_identical(r$variable, ref$variable, info = m)
    expect_identical(r$tangent, ref$tangent, info = m)
  }
  expect_identical(solve(models$ctl$rb4, tt, pars)$time, tt)
})

test_that("an event at the last time is not applied and the solve completes", {
  # The multistep methods used to restart past the end of the grid here and
  # fail with too many failed steps.
  for (m in methods) {
    r   <- solve(models$ev[[m]], 0:5, pars)
    ref <- solve(models$none[[m]], 0:5, pars)
    expect_identical(r$diagnostics$return_code, 0L, info = m)
    expect_identical(r$time, as.numeric(0:5), info = m)
    # The last row holds the state before the event.
    expect_identical(r$variable, ref$variable, info = m)
    expect_identical(r$tangent, ref$tangent, info = m)
  }
  r <- solve(models$ctl$rb4, 0:5, pars)
  expect_equal(r$variable, solve(models$none$rb4, 0:5, pars)$variable,
               tolerance = 1e-7)
})

test_that("an event at the first time is applied before the solve starts", {
  shifted <- pars; shifted["V"] <- pars[["V"]] + 1
  for (m in methods) {
    r   <- solve(models$ev[[m]], 5:8, pars)
    ref <- solve(models$none[[m]], 5:8, shifted)
    expect_identical(r$time, as.numeric(5:8), info = m)
    # The first row holds the state after the event, as every row at an
    # event time does.
    expect_identical(r$variable[1, ], shifted[c("V", "R")], info = m)
    expect_equal(r$variable, ref$variable, tolerance = 1e-8, info = m)
    expect_equal(r$tangent, ref$tangent, tolerance = 1e-6, info = m)
  }
  # A dose at time 0 on a grid that starts there.
  p0 <- ppar; p0[c("t_e", "dose")] <- c(0, 0.5)
  r <- solve(models$par$bdf, 0:3, p0)
  expect_identical(r$time, as.numeric(0:3))
  expect_identical(unname(r$variable[1, "V"]), -0.5)
  # includeTimeZero puts 0 in front of the grid, so the window starts there
  # and an event at the first requested time lies inside it.
  r <- solve(models$tz$bdf, 5:8, pars)
  expect_identical(r$time, c(0, 5, 6, 7, 8))
  pre <- solve(models$ev$bdf, 0:5, pars)
  expect_equal(unname(r$variable[2, ]), unname(pre$variable[6, ]) + c(1, 0),
               tolerance = 1e-7)
})

test_that("two consecutive windows reproduce one solve over both", {
  for (case in list(list(models$ev, pars), list(models$par, ppar))) {
    for (m in methods) {
      model <- case[[1]][[m]]
      p     <- case[[2]]
      full  <- solve(model, 0:8, p)
      one   <- solve(model, 0:5, p)
      p2    <- p; p2[c("V", "R")] <- one$variable[nrow(one$variable), ]
      two   <- solve(model, 5:8, p2)
      expect_identical(two$time, as.numeric(5:8), info = m)
      rows_full <- match(6:8, full$time)
      expect_equal(two$variable[2:4, ], full$variable[rows_full, ],
                   tolerance = 1e-7, info = m)
      for (i in 1:3)
        expect_equal(chain_tangent(one, two, i + 1),
                     full$tangent[rows_full[i], , ], tolerance = 1e-6,
                     info = paste(m, "t =", 5 + i))
    }
  }
})

test_that("the batch agrees with solveODE on every side of the window", {
  grids <- list(before = 6:8, after = seq(0, 4.5, by = 0.5), last = 0:5,
                first = 5:8)
  conds <- list(list(parms = pars),
                list(parms = replace(pars, "V", -0.5)))
  for (m in c("bdf", "rb4")) {
    for (g in names(grids)) {
      tt  <- grids[[g]]
      bat <- solveODEBatch(models$ev[[m]], conds, times = tt, cores = 2L)
      for (i in seq_along(conds)) {
        ser <- solveODE(models$ev[[m]], tt, conds[[i]]$parms)
        expect_identical(bat[[i]]$time, ser$time, info = paste(m, g))
        expect_identical(bat[[i]]$variable, ser$variable, info = paste(m, g))
        expect_identical(bat[[i]]$tangent, ser$tangent, info = paste(m, g))
      }
    }
  }

  # A parameter event time is evaluated per condition to size each result; one
  # batch holds an event on every side of the window [5, 8).
  t_es  <- c(before = 2, first = 5, inside = 6.5, on_grid = 7, last = 8, after = 9)
  conds <- lapply(t_es, function(te) list(parms = replace(ppar, "t_e", te)))
  for (m in c("bdf", "rb4")) {
    bat <- solveODEBatch(models$par[[m]], conds, times = 5:8, cores = 2L)
    for (i in seq_along(conds)) {
      ser <- solveODE(models$par[[m]], 5:8, conds[[i]]$parms)
      expect_identical(bat[[i]]$time, ser$time, info = names(t_es)[i])
      expect_identical(bat[[i]]$variable, ser$variable, info = names(t_es)[i])
      expect_identical(bat[[i]]$tangent, ser$tangent, info = names(t_es)[i])
    }
    n_rows <- vapply(bat, function(b) length(b$time), 0L)
    expect_identical(unname(n_rows), c(4L, 4L, 5L, 4L, 4L, 4L))
  }
})

test_that("the reverse mode applies the same window, with and without a store", {
  grids <- list(before = 6:8, after = seq(0, 4.5, by = 0.5), last = 0:5,
                first = 5:8, inside = 0:8)
  for (m in c("bdf", "rb4")) {
    for (g in names(grids)) {
      tt  <- grids[[g]]
      fwd <- solve(models$par[[m]], tt, ppar)
      set.seed(1)
      W   <- array(rnorm(length(fwd$variable)), c(dim(fwd$variable), 1L))
      rev <- solve(models$rev[[m]], tt, ppar, cotangent = W)
      expect_identical(length(rev$time), length(fwd$time), info = paste(m, g))
      ref <- contract(fwd$tangent, W)
      expect_equal(unname(rev$cotangent[names(ref), 1]), unname(ref),
                   tolerance = 1e-6, info = paste(m, g))

      one   <- solve(models$rev[[m]], tt, ppar, keepStore = TRUE)
      reuse <- solve(models$rev[[m]], tt, ppar, cotangent = W, store = one$store)
      expect_identical(reuse$cotangent, rev$cotangent, info = paste(m, g))
    }
  }
})

test_that("a reverse solve over a single time answers what the forward one does", {
  # One time integrates nothing: the cotangent of the initial states is the
  # seed, and the parameters get only what an event at that time contributes.
  # Such a run leaves no checkpoint, and the sweep used to read the state count
  # off one and crash.
  W <- array(c(0.3, -0.2), c(1L, 2L, 1L))
  for (m in c("bdf", "rb4")) {
    rv <- solve(models$rev0[[m]], 5, pars, cotangent = W)
    expect_identical(unname(rv$variable[1, ]), unname(pars[c("V", "R")]))
    expect_identical(unname(rv$cotangent[c("V", "R", "a", "b", "c"), 1]),
                     c(0.3, -0.2, 0, 0, 0), info = m)
    ref <- contract(solve(models$none[[m]], 5, pars)$tangent, W)
    expect_equal(unname(rv$cotangent[names(ref), 1]), unname(ref), info = m)

    # The event at that time fires, as it would at the start of a longer grid;
    # one at another time does not.
    for (te in c(5, 3)) {
      info <- paste(m, "t_e =", te)
      p   <- replace(ppar, "t_e", te)
      fwd <- solve(models$par[[m]], 5, p)
      expect_identical(unname(fwd$variable[1, "V"]), if (te == 5) 0 else -1,
                       info = info)
      ref <- contract(fwd$tangent, W)
      rev <- solve(models$rev[[m]], 5, p, cotangent = W)
      expect_identical(rev$variable, fwd$variable, info = info)
      expect_equal(unname(rev$cotangent[names(ref), 1]), unname(ref),
                   tolerance = 1e-12, info = info)
      if (te == 5) expect_true(ref[["t_e"]] != 0, info = info)

      one   <- solve(models$rev[[m]], 5, p, keepStore = TRUE)
      reuse <- solve(models$rev[[m]], 5, p, cotangent = W, store = one$store)
      expect_identical(reuse$cotangent, rev$cotangent, info = info)

      bat <- solveODEBatch(models$rev[[m]],
                           list(list(parms = p, cotangent = W),
                                list(parms = p, cotangent = 2 * W)),
                           times = 5, cores = 2L, abstol = tol$abstol,
                           reltol = tol$reltol)
      expect_identical(bat[[1]]$cotangent, rev$cotangent, info = info)
      expect_equal(bat[[2]]$cotangent, 2 * rev$cotangent, info = info)

      ff <- solve(models$ff[[m]], 5, p)
      fr <- solve(models$fr[[m]], 5, p, cotangent = W)
      expect_equal(unname(fr$cotangent[names(ref), 1]), unname(ref),
                   tolerance = 1e-12, info = info)
      expect_equal(unname(fr$curvature[, , 1]), hess_forward(ff, W),
                   tolerance = 1e-10, info = info)
    }
  }
})

test_that("the CVODE backend applies the same window", {
  skip_if_not(isTRUE(cvodeConfig$available), "CVODE backend not available")
  # The tangent's rows follow the model's flat vector, states then parameters.
  flat <- c(names(fhn), attr(cv_par, "parameters"))
  si <- diag(length(flat)); dimnames(si) <- list(NULL, flat)
  grids <- list(before = 6:8, after = seq(0, 4.5, by = 0.5), last = 0:5,
                first = 5:8, inside = 0:8)
  for (g in names(grids)) {
    tt <- grids[[g]]
    rc <- solve(cv_par, tt, ppar, tangent = si)
    rn <- solve(models$par$bdf, tt, ppar)
    expect_identical(rc$time, rn$time, info = g)
    expect_equal(rc$variable, rn$variable, tolerance = 1e-6, info = g)
    expect_equal(rc$tangent[, , colnames(si)], rn$tangent[, , colnames(si)],
                 tolerance = 1e-5, info = g)
  }
  # An event at the first time: the first row holds the state after it.
  r <- solve(cv_par, 5:8, ppar, tangent = si)
  expect_identical(unname(r$variable[1, ]), c(0, 1))

  t_es  <- c(before = 2, first = 5, inside = 6.5, last = 8, after = 9)
  conds <- lapply(t_es, function(te)
    list(parms = replace(ppar, "t_e", te), tangent = si))
  bat <- solveODEBatch(cv_par, conds, times = 5:8, cores = 2L)
  for (i in seq_along(conds)) {
    ser <- solveODE(cv_par, 5:8, conds[[i]]$parms, tangent = si)
    expect_identical(bat[[i]]$time, ser$time, info = names(t_es)[i])
    expect_identical(bat[[i]]$variable, ser$variable, info = names(t_es)[i])
    expect_identical(bat[[i]]$tangent, ser$tangent, info = names(t_es)[i])
  }
})
