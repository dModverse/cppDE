# Every integration method on problems with closed-form solutions and
# sensitivities, and the solver options of solveODE().

skip_on_cran()

# -- Models --------------------------------------------------------------------

# A decay chain with a closed form, and a state driven explicitly by the clock.
eqns_decay <- c(A = "-k1 * A", B = "k1 * A - k2 * B", x = "-x + cos(5 * time)")
times <- seq(0, 50, length.out = 200)
pars  <- c(A = 1, B = 0, x = 1, k1 = 0.1, k2 = 0.2)
methods_all <- c("bdf", "adams", "rb4", "tsit5")

# A, B and their sensitivities in closed form, [length(t), state, A/B/k1/k2].
exact_decay <- function(t, A0 = 1, B0 = 0, k1 = 0.1, k2 = 0.2) {
  e1 <- exp(-k1 * t); e2 <- exp(-k2 * t); r <- k1 / (k2 - k1)
  A <- A0 * e1
  B <- B0 * e2 + A0 * r * (e1 - e2)
  dr1 <- k2 / (k2 - k1)^2; dr2 <- -k1 / (k2 - k1)^2
  sens <- array(0, c(length(t), 2, 4),
                list(NULL, c("A", "B"), c("A", "B", "k1", "k2")))
  sens[, "A", "A"]  <- e1
  sens[, "A", "k1"] <- -t * A
  sens[, "B", "A"]  <- r * (e1 - e2)
  sens[, "B", "B"]  <- e2
  sens[, "B", "k1"] <- A0 * (dr1 * (e1 - e2) - r * t * e1)
  sens[, "B", "k2"] <- -B0 * t * e2 + A0 * (dr2 * (e1 - e2) + r * t * e2)
  list(A = A, B = B, sens = sens)
}

per_method <- function(ms, prefix, ...)
  lapply(setNames(nm = ms), function(m)
    cppODE(..., method = m, modelname = paste0(prefix, m), compile = FALSE))

decay <- per_method(methods_all, "mth_decay_", eqns_decay)

# Second order on the decay chain and on a 10^x term, whose derivatives are
# math calls on literals.
eqns_d2 <- c(A = "-k1 * A", B = "k1 * A - k2 * B", z = "-kz * 10^z")
decay_d2 <- per_method(c("bdf", "rb4"), "mth_d2_", eqns_d2, deriv = TRUE, deriv2 = TRUE)

robertson <- cppODE(c(y1 = "-k1*y1 + k3*y2*y3", y2 = "k1*y1 - k2*y2^2 - k3*y2*y3",
                      y3 = "k2*y2^2"),
                    deriv2 = TRUE, fixed = c("y1", "y2", "y3"),
                    modelname = "mth_sens2_sym", compile = FALSE)

# The BDF corrector without NDF coefficients, and a compile-time fixed rate.
ndf_fixed <- cppODE(eqns_decay[c("A", "B")], useNDF = FALSE, fixed = "k2",
                    modelname = "mth_bdf_fixed", compile = FALSE)

cxx_tokens <- cppODE(c(default = "-std * default + operator",
                       int     = "std * default - int * 10^0.5"),
                     modelname = "mth_cxx_tokens_ode", compile = FALSE)

forced <- cppODE(c(A = "-k1 * A * u + k2 * B",
                   B = "k1 * A * u - k2 * B"),
                 forcings = "u", modelname = "mth_forcing_in_jac", compile = FALSE)

# Integration stops where the root expression first crosses zero.
stopper <- cppODE(c(A = "-k1 * A"), rootfunc = "A - 0.25", modelname = "mth_rootfunc",
                  compile = FALSE)

native <- c(decay, decay_d2, list(robertson, ndf_fixed, cxx_tokens, forced, stopper))
do.call(compile, c(unname(native), list(output = "test_ode_methods", cores = test_cores())))

has_cvode <- isTRUE(cvodeConfig$available)
if (has_cvode) {
  # Adams, a forcing in the Jacobian and a compile-time fixed rate on CVODES.
  cv_forced <- cvode(c(A = "-k1 * A * u + k2 * B", B = "k1 * A * u - k2 * B"),
                     forcings = "u", method = "adams", fixed = "k2", deriv = TRUE,
                     modelname = "mth_cv_forced", compile = FALSE)
  cv_stopper <- cvode(c(A = "-k1 * A"), rootfunc = "A - 0.25", modelname = "mth_cv_rootfunc",
                      compile = FALSE)
  compile(cv_forced, cv_stopper, output = "test_ode_methods_cvode", cores = test_cores())
}

u_data <- data.frame(time = c(0, 0.5, 1, 2), value = c(0.4, 1.1, 0.7, 1.5))

# -- Accuracy against closed forms --------------------------------------------

test_that("every method returns the closed-form states and sensitivities", {
  ex <- exact_decay(times)
  for (m in methods_all) {
    res <- solveODE(decay[[m]], times, pars, abstol = 1e-10, reltol = 1e-10)
    expect_identical(res$time, times, info = m)
    expect_identical(colnames(res$variable), c("A", "B", "x"), info = m)
    expect_identical(dimnames(res$tangent)$sens, names(pars), info = m)
    expect_equal(unname(res$variable[, "A"]), ex$A, tolerance = 1e-7, info = m)
    expect_equal(unname(res$variable[, "B"]), ex$B, tolerance = 1e-7, info = m)
    expect_equal(unname(res$tangent[, c("A", "B"), c("A", "B", "k1", "k2")]),
                 unname(ex$sens), tolerance = 1e-6, info = m)
  }
})

test_that("all methods keep their order on an explicitly time-dependent system", {
  # A wrong df/dt term costs a Rosenbrock method its order and inflates its
  # step count.
  exact <- function(t) (1 - 1/26) * exp(-t) + (cos(5 * t) + 5 * sin(5 * t)) / 26
  for (m in methods_all) {
    lo <- solveODE(decay[[m]], c(0, 10), pars, abstol = 1e-6, reltol = 1e-6)
    hi <- solveODE(decay[[m]], c(0, 10), pars, abstol = 1e-9, reltol = 1e-9)
    expect_lt(abs(hi$variable[2, "x"] - exact(10)), 1e-7, label = m)
    expect_lt(hi$diagnostics$accepted, 15 * lo$diagnostics$accepted, label = m)
  }
})

# -- Second order --------------------------------------------------------------

test_that("second-order sensitivities match their closed forms", {
  tt <- seq(0, 1, length.out = 25)
  p  <- c(A = 1, B = 0, z = 0.3, k1 = 0.1, k2 = 0.2, kz = 0.7)
  ln10 <- log(10)
  # 10^(-z(t)) = 10^(-z0) + kz ln(10) t
  u <- 10^(-p[["z"]]) + p[["kz"]] * ln10 * tt
  A <- exp(-p[["k1"]] * tt)
  for (m in names(decay_d2)) {
    res <- solveODE(decay_d2[[m]], tt, p, abstol = 1e-12, reltol = 1e-12)
    expect_equal(unname(res$hessian[, "A", "k1", "k1"]), tt^2 * A, tolerance = 1e-6, info = m)
    expect_equal(unname(res$hessian[, "A", "A", "k1"]), -tt * A, tolerance = 1e-6, info = m)
    expect_equal(unname(res$hessian[, "A", "k2", "k2"]), rep(0, length(tt)),
                 tolerance = 1e-10, info = m)
    expect_equal(unname(res$variable[, "z"]), -log10(u), tolerance = 1e-8, info = m)
    expect_equal(unname(res$tangent[, "z", "z"]), 10^(-p[["z"]]) / u, tolerance = 1e-6,
                 info = m)
    expect_equal(unname(res$tangent[, "z", "kz"]), -tt / u, tolerance = 1e-6, info = m)
    expect_equal(unname(res$hessian[, "z", "kz", "kz"]), ln10 * tt^2 / u^2,
                 tolerance = 1e-6, info = m)
    h <- res$hessian
    expect_identical(unname(h), unname(aperm(h, c(1, 2, 4, 3))), info = m)
  }
})

test_that("second-order sensitivities are symmetric over a long stiff run", {
  res <- solveODE(robertson, c(0, 10^seq(-2, 3, length.out = 30)),
                  c(y1 = 1, y2 = 0, y3 = 0, k1 = 0.04, k2 = 3e7, k3 = 1e4),
                  abstol = 1e-12, reltol = 1e-10)
  s <- res$hessian[31, , , ]
  expect_equal(unname(s), unname(aperm(s, c(1, 3, 2))), tolerance = 1e-12)
  expect_lt(max(abs(s)), 1e3)
})

# -- Solver options ------------------------------------------------------------

test_that("BDF without NDF coefficients is accurate and says so", {
  p <- pars[c("A", "B", "k1", "k2")]
  res <- solveODE(ndf_fixed, times, p, abstol = 1e-10, reltol = 1e-10)
  ex <- exact_decay(times)
  expect_equal(unname(res$variable[, "A"]), ex$A, tolerance = 1e-7)
  expect_equal(unname(res$variable[, "B"]), ex$B, tolerance = 1e-7)
  expect_false(res$diagnostics$useNDF)
  expect_true(solveODE(decay$bdf, times, pars)$diagnostics$useNDF)
  # A compile-time fixed parameter has no sensitivity column at all.
  expect_identical(attr(ndf_fixed, "dimNames")$sens, c("A", "B", "k1"))
  expect_identical(dimnames(res$tangent)$sens, c("A", "B", "k1"))
  expect_equal(unname(res$tangent[, , "k1"]), unname(ex$sens[, , "k1"]), tolerance = 1e-6)
})

test_that("sensErrCon = FALSE drops the sensitivities from the error test", {
  # Fewer components in the error norm change the steps, in either direction.
  ex <- exact_decay(times)
  on  <- solveODE(decay$bdf, times, pars, abstol = 1e-8, reltol = 1e-8)
  off <- solveODE(decay$bdf, times, pars, abstol = 1e-8, reltol = 1e-8,
                  sensErrCon = FALSE)
  expect_false(identical(off$tangent, on$tangent))
  expect_equal(unname(off$variable[, "A"]), ex$A, tolerance = 1e-6)
  expect_equal(unname(off$tangent[, "A", "k1"]), unname(ex$sens[, "A", "k1"]),
               tolerance = 1e-4)
})

test_that("a given first step size is taken and leaves the answer alone", {
  ex <- exact_decay(times)
  res <- solveODE(decay$rb4, times, pars, abstol = 1e-10, reltol = 1e-10, hini = 1e-4)
  expect_equal(unname(res$variable[, "A"]), ex$A, tolerance = 1e-7)
})

test_that("onFailure decides between an error, a warning and silence", {
  tt <- c(0, 1, 2)
  p  <- pars
  expect_error(solveODE(decay$bdf, tt, p, maxsteps = 2L), "did not complete")
  expect_warning(solveODE(decay$bdf, tt, p, maxsteps = 2L, onFailure = "warn"),
                 "did not complete")
  partial <- suppressWarnings(solveODE(decay$bdf, tt, p, maxsteps = 2L, onFailure = "warn"))
  expect_lt(length(partial$time), length(tt))
  expect_lt(partial$diagnostics$return_code, 0L)
  quiet <- expect_silent(solveODE(decay$bdf, tt, p, maxsteps = 2L, onFailure = "silent"))
  expect_identical(quiet$variable, partial$variable)
})

test_that("diagnostics() reports the statistics of a solve", {
  res <- solveODE(decay$bdf, times, pars)
  expect_output(d <- diagnostics(res), "NDF solver statistics")
  expect_identical(d$return_code, 0L)
  expect_gt(d$accepted, 0)
  expect_gt(d$fevals, d$accepted)
  expect_equal(d$t_reached, max(times))
  expect_identical(d$method, "bdf")
})

test_that("a forgotten symbol pairing is found again", {
  ref <- solveODE(decay$tsit5, c(0, 1), pars)
  expect_null(clearNativeSymbols())
  expect_identical(solveODE(decay$tsit5, c(0, 1), pars), ref)
})

test_that("rootfunc stops the integration at the first zero crossing", {
  # The crossing of A = 0.25 under A' = -k1 A lies at log(4) / k1.
  res <- solveODE(stopper, 0:20, c(A = 1, k1 = 0.2), abstol = 1e-10, reltol = 1e-10,
                  roottol = 1e-10)
  last <- length(res$time)
  expect_equal(res$time[last], log(4) / 0.2, tolerance = 1e-8)
  expect_equal(unname(res$variable[last, "A"]), 0.25, tolerance = 1e-8)
  expect_identical(res$time[-last], as.numeric(0:6))

  skip_if_not(has_cvode, "CVODE backend not available")
  cv <- solveODE(cv_stopper, 0:20, c(A = 1, k1 = 0.2), abstol = 1e-10, reltol = 1e-10,
                 roottol = 1e-10)
  last <- length(cv$time)
  expect_equal(cv$time[last], log(4) / 0.2, tolerance = 1e-6)
  expect_identical(cv$time[-last], as.numeric(0:6))
})

# -- Symbol names -------------------------------------------------------------

# The generated right-hand side indexes x[] and params[] and calls std::pow, so
# a state or parameter with one of those names has to be substituted before
# it can be read as part of the surrounding code.
test_that("state and parameter names that are C++ tokens compile and solve", {
  tt <- seq(0, 2, 0.5)
  p  <- c(default = 1, int = 0.3, std = 0.7, operator = 0.2)
  res <- solveODE(cxx_tokens, tt, p, abstol = 1e-11, reltol = 1e-11)

  s <- p[["std"]]; cc <- 10^0.5; dinf <- p[["operator"]] / s
  d <- dinf + (p[["default"]] - dinf) * exp(-s * tt)
  int <- p[["int"]] * exp(-cc * tt) + s * dinf * (1 - exp(-cc * tt)) / cc +
    s * (p[["default"]] - dinf) * (exp(-s * tt) - exp(-cc * tt)) / (cc - s)
  expect_equal(unname(res$variable[, "default"]), d, tolerance = 1e-9)
  expect_equal(unname(res$variable[, "int"]), int, tolerance = 1e-9)
  expect_equal(unname(res$tangent[, "default", "default"]), exp(-s * tt), tolerance = 1e-8)
  expect_equal(unname(res$tangent[, "default", "operator"]), (1 - exp(-s * tt)) / s,
               tolerance = 1e-8)
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

# -- Forcings -----------------------------------------------------------------

test_that("a forcing that multiplies a state reaches the Jacobian", {
  # Only additive forcings vanish from df/dx. The derivative of the solution has
  # no closed form under an interpolated forcing, so central differences of two
  # tight solves are the reference.
  tt  <- seq(0, 2, 0.25)
  p <- c(A = 1, B = 0, k1 = 0.8, k2 = 0.3)
  res <- solveODE(forced, tt, p, forcings = list(u = u_data),
                  abstol = 1e-10, reltol = 1e-10)
  fd <- vapply(names(p), function(nm) {
    h <- 1e-6 * max(abs(p[[nm]]), 1)
    pp <- pm <- p; pp[nm] <- pp[nm] + h; pm[nm] <- pm[nm] - h
    a <- solveODE(forced, tt, pp, forcings = list(u = u_data), abstol = 1e-12, reltol = 1e-12)
    b <- solveODE(forced, tt, pm, forcings = list(u = u_data), abstol = 1e-12, reltol = 1e-12)
    (a$variable - b$variable) / (2 * h)
  }, matrix(0, length(tt), 2L))
  # Loose because each adaptive solve takes its own grid; a missing forcing term
  # is orders of magnitude larger.
  for (k in seq_along(p))
    expect_equal(unname(res$tangent[, , k]), unname(fd[, , k]), tolerance = 1e-3,
                 info = names(p)[k])
})

test_that("a forcing of one point is the constant of two", {
  tt <- seq(0, 3, 0.5)
  p <- c(A = 1, B = 0, k1 = 0.8, k2 = 0.3)
  one <- solveODE(forced, tt, p, forcings = list(u = data.frame(time = 1, value = 0.7)),
                  abstol = 1e-10, reltol = 1e-10)
  two <- solveODE(forced, tt, p,
                  forcings = list(u = cbind(c(0, 3), 0.7)),
                  abstol = 1e-10, reltol = 1e-10)
  expect_equal(one$variable, two$variable, tolerance = 1e-8)
  expect_equal(one$tangent, two$tangent, tolerance = 1e-8)
})

test_that("CVODES Adams takes the forcing and the fixed rate as cppODE does", {
  skip_if_not(has_cvode, "CVODE backend not available")
  tt <- seq(0, 2, 0.25)
  p <- c(A = 1, B = 0, k1 = 0.8, k2 = 0.3)
  nat <- solveODE(forced, tt, p, forcings = list(u = u_data), abstol = 1e-10, reltol = 1e-10)
  cv  <- solveODE(cv_forced, tt, p, forcings = list(u = u_data),
                  abstol = 1e-10, reltol = 1e-10)
  expect_identical(cv$diagnostics$method, "adams")
  expect_identical(dimnames(cv$tangent)$sens, c("A", "B", "k1"))
  expect_equal(cv$variable, nat$variable, tolerance = 1e-7)
  expect_equal(unname(cv$tangent), unname(nat$tangent[, , c("A", "B", "k1")]),
               tolerance = 1e-6)
})
