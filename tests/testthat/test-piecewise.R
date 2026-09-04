# Piecewise support in generated models. A piecewise is emitted as
# cppde::select(cond, a, b): the branches of a ternary cannot share a type once
# one of them is an expression-template node, and the condition reads the value.

skip_on_cran()

library(cppDE)

# -- Emitted form -------------------------------------------------------------

test_that("a piecewise is emitted as cppde::select, not as a ternary", {
  outdir <- tempfile("pw_emit_")
  dir.create(outdir)
  cppODE(c(A = "-k*A*piecewise(0, A - thr > 0, 1)"),
         modelname = "pw_emit", outdir = outdir, deriv = TRUE, compile = FALSE)
  body <- grep("dxdt\\[", readLines(file.path(outdir, "pw_emit.cpp")), value = TRUE)
  expect_true(any(grepl("cppde::select(", body, fixed = TRUE)))
  expect_false(any(grepl("?", body, fixed = TRUE)))
})

# -- Solve and sensitivities --------------------------------------------------

test_that("a time switch integrates and differentiates like its closed form", {
  # Switching on time keeps the crossing independent of the parameters, so the
  # forward tangents of the switched system are the exact derivatives.
  f <- cppODE(c(A = "-piecewise(kf*A, time - ts < 0, ks*A)",
                B = " piecewise(kf*A, time - ts < 0, ks*A)"),
              modelname = "pw_time_switch", deriv = TRUE)

  times <- seq(0, 6, by = 0.5)
  p <- c(kf = 0.5, ks = 0.05, ts = 2, A = 1, B = 0)
  out <- solveODE(f, times = times, parms = p, abstol = 1e-10, reltol = 1e-10)

  tk  <- pmin(times, p[["ts"]])          # time spent on the fast branch
  tl  <- pmax(times - p[["ts"]], 0)      # time spent on the slow branch
  A   <- p[["A"]] * exp(-p[["kf"]] * tk - p[["ks"]] * tl)
  expect_equal(unname(out$variable[, "A"]), A, tolerance = 1e-6)
  expect_equal(unname(out$variable[, "B"]), p[["A"]] - A, tolerance = 1e-6)

  expect_equal(unname(out$sens1[, "B", "kf"]), tk * A, tolerance = 1e-6)
  expect_equal(unname(out$sens1[, "B", "ks"]), tl * A, tolerance = 1e-6)
  expect_equal(unname(out$sens1[, "B", "A"]), 1 - A / p[["A"]], tolerance = 1e-6)
})

test_that("both branches of a state switch are taken", {
  f <- cppODE(c(A = "-piecewise(kf*A, A - thr > 0, ks*A)"),
              modelname = "pw_state_switch", deriv = TRUE)

  p <- c(kf = 0.5, ks = 0.05, thr = 0.6, A = 1)
  ts <- -log(p[["thr"]] / p[["A"]]) / p[["kf"]]     # crossing of the threshold
  times <- c(0, ts / 2, ts, ts + 1, ts + 3)
  out <- solveODE(f, times = times, parms = p, abstol = 1e-10, reltol = 1e-10)

  tk <- pmin(times, ts)
  tl <- pmax(times - ts, 0)
  expect_equal(unname(out$variable[, "A"]),
               p[["A"]] * exp(-p[["kf"]] * tk - p[["ks"]] * tl),
               tolerance = 1e-6)
})

# -- Second order -------------------------------------------------------------

test_that("select carries value, gradient and Hessian on both branches", {
  # Same oracle as test-dual2nd-primitives.R: the dual2nd path and the
  # SymPy-derived closed forms are independent implementations.
  expr <- c(y = "piecewise(a^2*b, a - 1 > 0, b*a + a^3)")
  nms <- c("a", "b")
  dP  <- diag(2); dimnames(dP) <- list(nms, nms)
  dP2 <- array(0, c(2, 2, 2), dimnames = list(nms, nms, nms))

  out <- lapply(c(dual = "dual", symbolic = "symbolic"), function(mode) {
    # The R fallback cannot parse a piecewise; only the compiled path is used.
    f <- suppressWarnings(
      funCpp(expr, parameters = nms, deriv = TRUE, deriv2 = TRUE,
             derivMode = mode, compile = TRUE,
             modelname = paste0("pw_d2_", mode)))
    lapply(c(0.5, 2), function(a)
      f$evaluate(a = a, b = 3, dP = dP, dP2 = dP2, deriv2 = TRUE))
  })

  for (i in seq_along(out$dual)) {
    expect_equal(unname(out$dual[[i]]$y),   unname(out$symbolic[[i]]$y))
    expect_equal(unname(out$dual[[i]]$dy),  unname(out$symbolic[[i]]$dy),
                 tolerance = 1e-10)
    expect_equal(unname(out$dual[[i]]$d2y), unname(out$symbolic[[i]]$d2y),
                 tolerance = 1e-10)
  }
})

test_that("a branch that is a literal clears the tangents it replaces", {
  # The second-order materialiser leaves the target's tangent buffers alone
  # when a tree carries no dependence, and dxdt is reused across calls, so a
  # select landing on a literal has to stay on the writing path.
  f <- cppODE(c(A = "-piecewise(ks*A, time - ts < 0, 0)"),
              modelname = "pw_freeze", deriv = TRUE, deriv2 = TRUE)

  times <- seq(0, 4, by = 0.5)
  p <- c(ks = 0.4, ts = 2, A = 1)
  out <- solveODE(f, times = times, parms = p, abstol = 1e-11, reltol = 1e-11)

  tk <- pmin(times, p[["ts"]])          # the decay freezes at ts
  A  <- p[["A"]] * exp(-p[["ks"]] * tk)
  expect_equal(unname(out$variable[, "A"]), A, tolerance = 1e-7)
  expect_equal(unname(out$sens1[, "A", "ks"]), -tk * A, tolerance = 1e-7)
  expect_equal(unname(out$sens1[, "A", "A"]), A / p[["A"]], tolerance = 1e-7)
  expect_equal(unname(out$sens2[, "A", "ks", "ks"]), tk^2 * A, tolerance = 1e-7)
  expect_equal(unname(out$sens2[, "A", "ks", "A"]), -tk * A / p[["A"]],
               tolerance = 1e-7)
})

# -- Conditions --------------------------------------------------------------

test_that("the R and C spellings of the logical operators parse", {
  # && and || bind below the comparisons while & and | bind above them, so the
  # grouping has to survive the rewrite even without the parentheses.
  eqn <- function(cond) c(A = paste0("-k*A*piecewise(1, ", cond, ", 0)"))
  times <- seq(0, 6, by = 0.5)
  p <- c(k = 0.4, t1 = 1, t2 = 3, A = 1)
  run <- function(cond, nm)
    solveODE(cppODE(eqn(cond), modelname = nm, deriv = FALSE),
             times = times, parms = p, abstol = 1e-11, reltol = 1e-11)$variable[, "A"]

  bare  <- run("time > t1 && time <= t2",     "pw_and_bare")
  wrapped <- run("(time > t1) && (time <= t2)", "pw_and_wrapped")
  expect_equal(unname(bare), unname(wrapped))

  # Decay runs between t1 and t2 only, so A is flat on either side.
  tk <- pmin(pmax(times - p[["t1"]], 0), p[["t2"]] - p[["t1"]])
  expect_equal(unname(bare), p[["A"]] * exp(-p[["k"]] * tk), tolerance = 1e-7)

  # The complement, spelled with || and !, has to give the mirror image.
  outside <- run("!(time > t1) || time > t2", "pw_or_not")
  expect_equal(unname(outside),
               p[["A"]] * exp(-p[["k"]] * (times - tk)), tolerance = 1e-7)
})

test_that("an expression that does not parse names itself", {
  # The reason has to survive the trip through reticulate, which truncates a
  # long message and then indexes it with an offset from the untruncated one.
  long <- paste0("-k*A ** * A", strrep(" + 0*A", 200))
  expect_error(cppODE(c(A = long), modelname = "pw_unparseable", compile = FALSE),
               "cannot parse expression")
})
