# Argument checks and guards that run without generating or compiling a model.

# A directory holding executables of the given names, each a no-op script.
fake_path <- function(tools) {
  d <- tempfile("fakebin_")
  dir.create(d)
  for (t in tools) {
    f <- file.path(d, t)
    writeLines("#!/bin/sh\nexit 0", f)
    Sys.chmod(f, "0755")
  }
  d
}

# Runs `code` with PATH set to `path`, restoring it afterwards.
with_path <- function(path, code) {
  old <- Sys.getenv("PATH")
  on.exit(Sys.setenv(PATH = old), add = TRUE)
  Sys.setenv(PATH = path)
  force(code)
}

# -- codegenAvailable() ---------------------------------------------------------

test_that("codegenAvailable() answers TRUE or FALSE and checks its backend", {
  expect_true(isTRUE(codegenAvailable()) || isFALSE(codegenAvailable()))
  expect_error(codegenAvailable("python"), "should be one of")
  # The cvode backend needs the native requirements and SUNDIALS on top.
  if (!isTRUE(cppDE:::cvodeConfig$available)) expect_false(codegenAvailable("cvode"))
  if (isTRUE(cppDE:::cvodeConfig$available))
    expect_identical(codegenAvailable("cvode"), codegenAvailable("native"))
})

test_that("codegenAvailable() is FALSE without a compiler on the PATH", {
  skip_on_os("windows")
  expect_false(with_path(fake_path(character(0)), codegenAvailable()))
})

# -- installLibs() --------------------------------------------------------------

test_that("installLibs() checks its tools before it touches anything", {
  skip_on_os("windows")
  dir <- tempfile("libs_")
  expect_error(installLibs("all", dir = dir), "should be one of")
  expect_error(with_path(fake_path(character(0)), installLibs(dir = dir, ask = FALSE)),
               "cmake was not found")
  expect_error(with_path(fake_path("cmake"), installLibs(dir = dir, ask = FALSE)),
               "need one of curl, wget or git")
  expect_false(dir.exists(dir))
})

test_that("a declined installLibs() downloads and writes nothing", {
  skip_on_os("windows")
  local_mocked_bindings(askYesNo = function(...) FALSE, .package = "utils")
  dir <- tempfile("libs_")
  expect_message(
    out <- with_path(fake_path(c("cmake", "curl")), installLibs("suitesparse", dir = dir,
                                                               ask = TRUE)),
    "nothing was downloaded")
  expect_null(out)
  expect_false(dir.exists(dir))
})

test_that("installLibs() points Windows to Rtools", {
  skip_if(.Platform$OS.type != "windows", "not Windows")
  expect_error(installLibs(ask = FALSE), "pacman")
})

test_that("install_libs() is a deprecated alias of installLibs()", {
  expect_warning(expect_error(install_libs(which = "none"), "should be one of"),
                 "deprecated")
})

# -- adjointControl() -----------------------------------------------------------

test_that("adjointControl() returns its three settings and checks each", {
  ctl <- adjointControl()
  expect_s3_class(ctl, "cppDEadjointControl")
  expect_identical(unclass(ctl), list(refine = FALSE, gradtol = NULL, trace = FALSE))
  expect_identical(adjointControl(TRUE, 1e-8, TRUE)$gradtol, 1e-8)

  for (bad in list(NA, "yes", c(TRUE, FALSE), 1))
    expect_error(adjointControl(refine = bad), "'refine' must be TRUE or FALSE")
  for (bad in list(NA, "yes", c(TRUE, FALSE)))
    expect_error(adjointControl(trace = bad), "'trace' must be TRUE or FALSE")
  for (bad in list(0, -1, c(1, 2), "1e-8", NA_real_))
    expect_error(adjointControl(gradtol = bad), "'gradtol' must be a positive number")
})

test_that("an adjoint trace prints its size", {
  g <- structure(list(time = c(0, 1, 2), h = c(1, 1, 1),
                      lambda = array(0, c(3L, 2L, 4L)), eta = matrix(0, 3L, 4L)),
                 class = "cppDEadjoint")
  expect_output(out <- print(g), "<cppDEadjoint> 3 steps, 2 states, 4 seed\\(s\\)")
  expect_identical(out, g)
})

# -- diagnostics() --------------------------------------------------------------

test_that("diagnostics() prints the statistics and returns them", {
  d <- list(return_code = -1L, accepted = 3L, rejected = 1L, fevals = 10L,
            jevals = 2L, setups = 1L, last_dt = 0.1, last_order = 2L,
            t_reached = 1.5, method = "bdf", useNDF = FALSE)
  expect_output(out <- diagnostics(list(diagnostics = d)), "BDF solver statistics")
  expect_identical(out, d)
  expect_output(diagnostics(list(diagnostics = d)), "Too much work")
  expect_output(diagnostics(list(diagnostics = modifyList(d, list(useNDF = TRUE)))),
                "NDF solver statistics")
  expect_output(diagnostics(list(diagnostics = modifyList(d, list(backend = "cvode",
                                                                   method = "adams")))),
                "CVODE: ADAMS solver statistics")
  expect_output(diagnostics(list(diagnostics = modifyList(d, list(return_code = -77L)))),
                "Unknown return code: -77")
  expect_message(expect_null(diagnostics(list(time = 0))), "No solver diagnostics")
})

# -- forcingValues() ------------------------------------------------------------

test_that("forcingValues() interpolates through the points and holds the ends", {
  u <- data.frame(time = c(0, 1, 3), value = c(1, 3, 7))
  # Points on one line have equal secants, so the PCHIP interpolant is that line.
  tt <- c(-1, 0, 0.25, 1, 2.5, 3, 4)
  v <- forcingValues(tt, list(u = u))
  expect_identical(dim(v), c(length(tt), 1L))
  expect_identical(colnames(v), "u")
  expect_equal(unname(v[, "u"]), c(1, 1, 1.5, 3, 6, 7, 7), tolerance = 1e-14)

  # A monotone forcing stays within the range of each interval, unsorted points
  # are sorted, a matrix gives what a data.frame gives, and columns follow names.
  w <- data.frame(time = c(2, 0, 1, 3), value = c(1, 0, 0.9, 5))
  tt <- seq(0, 3, by = 0.05)
  ww <- forcingValues(tt, list(a = w, b = cbind(w$time, w$value)))
  expect_identical(ww[, "a"], ww[, "b"])
  expect_true(all(diff(ww[, "a"]) >= 0))
  expect_equal(unname(forcingValues(c(0, 1, 2, 3), list(w = w))[, 1]), c(0, 0.9, 1, 5))
  expect_equal(unname(forcingValues(c(-5, 5), list(c = data.frame(time = 1, value = 2)))[, 1]),
               c(2, 2))
})

test_that("forcingValues() rejects malformed input by name", {
  ok <- data.frame(time = 0:1, value = 1:2)
  expect_error(forcingValues(c(1, NA), list(u = ok)), "'times' must be finite")
  expect_error(forcingValues(c(1, Inf), list(u = ok)), "'times' must be finite")
  expect_error(forcingValues(1, ok), "'forcings' must be a named list")
  expect_error(forcingValues(1, list(ok)), "'forcings' must be a named list")
  expect_error(forcingValues(1, list(u = data.frame(t = 0, value = 1))),
               "Forcing 'u' needs columns 'time' and 'value'")
  expect_error(forcingValues(1, list(u = ok[0, ])), "Forcing 'u' needs a time point")
  expect_error(forcingValues(1, list(u = data.frame(time = c(0, 0), value = 1:2))),
               "Forcing 'u': duplicate times")
  expect_error(forcingValues(1, list(u = data.frame(time = c(0, Inf), value = 1:2))),
               "Forcing 'u': non-finite time")
  expect_error(forcingValues(1, list(u = data.frame(time = 0:1, value = c(1, NA)))),
               "Forcing 'u': non-finite value")
})

# -- Model builders -------------------------------------------------------------

test_that("the builders reject contradictory arguments before generating code", {
  ev_both <- data.frame(var = "A", time = 1, root = "A - 1", value = 1, method = "add")
  expect_error(cppODE(c(A = "-k*A"), derivMode = "reverse", deriv2 = TRUE),
               "forward-reverse")
  expect_error(cppODE(c(A = "-k*A"), events = ev_both),
               "exactly one of 'time' or 'root'")
  expect_error(cppODE(c(A = "-k*A"), forcings = "u"), "Unknown forcing symbols: u")
  expect_error(cppODE(c(A = "-k*A*u"), forcings = "A"),
               "Forcing names cannot be state variables: A")
  expect_error(cppODE(c(A = "-k*A"), method = "rk4"), "should be one of")
  expect_error(cppFUN(c(y = "a"), parameters = "a", deriv2 = TRUE, derivMode = "reverse"),
               "no reverse counterpart")
  expect_error(cppFUN(c(y = "a"), outdir = file.path(tempfile(), "none")),
               "outdir does not exist")

  skip_if_not(isTRUE(cppDE:::cvodeConfig$available), "CVODE backend not available")
  expect_error(cvode(c(A = "-k*A"), derivMode = "reverse", deriv = TRUE),
               "no forward sensitivities")
  expect_error(cvode(c(A = "-k*A"), asaCheckpoints = 0), "'asaCheckpoints' must be")
  expect_error(cvode(c("-k*A")), "'rhs' must be a named character vector")
  expect_error(cvode(c(A = "-k*A"), events = ev_both), "exactly one of 'time' or 'root'")
  expect_error(cvode(c(A = "-k*A"), forcings = "u"), "Unknown forcing symbols: u")
  expect_error(cvode(c(A = "-k*A*u"), forcings = "A"),
               "Forcing names cannot be state variables: A")
  expect_error(cvode(c(A = "-k*A"), method = "rb4"), "should be one of")
})
