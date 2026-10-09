# Piecewise support in generated models. A piecewise is emitted as
# cppde::select(cond, a, b): the branches of a ternary cannot share a type once
# one of them is an expression-template node, and the condition reads the value.

skip_on_cran()

library(cppDE)

# Every model the solve tests below run, built once and linked into one shared
# object. The emitted-form and parse-error tests build their own.
pw_cond <- function(cond, nm)
  cppODE(c(A = paste0("-k*A*piecewise(1, ", cond, ", 0)")), modelname = nm,
         deriv = FALSE, compile = FALSE)
# A step input switched on at ts, in both spellings of the condition: x sees it
# from ts on, y from the next double. a and b decay smoothly and keep the order
# of a multistep method up.
pw_step <- function(nm, ...)
  cppODE(c(a = "-0.1*a", b = "0.1*a - 0.05*b",
           x = "-k*(x - c*piecewise(0, time - ts < 0, 1))",
           y = "-k*y + k*c*piecewise(1, time > ts, 0)"),
         modelname = nm, compile = FALSE, ...)
# A pulse between ts and t2 on a state at rest, where f and with it the error
# estimate are zero until ts.
pw_pulse <- function(nm, ...)
  cppODE(c(x = "-k*x + c*piecewise(1, time > ts && time <= t2, 0)",
           y = "k*x - y"),
         modelname = nm, compile = FALSE, ...)
# Every operator, each in a state of its own. x runs with time, so a condition on
# x is a switch on a state and one on time a switch in time.
pw_ops <- c(
  x = "1",
  arith = "a*x^2 - x**3/b + (x - a)/(b + x)",
  cmp_time = "(time > t1) + 2*(time <= t2)",
  cmp_state = "(x >= t1)*(x < t2)",
  logic = "a*(x > t1 && x < t2) + (time > t2 || x < h)",
  not_mul = "3*!(x > t1)",
  not_prec = "!x > t1",
  eq_param = "(mode == 1)*x + (mode != 1)",
  pw_and = "piecewise(a, time >= t1 && x < t2, b)",
  heaviside = "a*Heaviside(x - t1)"
)
has_cvode <- isTRUE(cvodeConfig$available)
pw_mod <- list(
  time_switch = cppODE(c(A = "-piecewise(kf*A, time - ts < 0, ks*A)",
                         B = " piecewise(kf*A, time - ts < 0, ks*A)"),
                       modelname = "pw_time_switch", deriv = TRUE,
                       compile = FALSE),
  d2 = cppFUN(c(y = "piecewise(a^2*b, a - 1 > 0, b*a + a^3)"),
              parameters = c("a", "b"), deriv = TRUE, deriv2 = TRUE,
              derivMode = "forward", modelname = "pw_d2"),
  freeze = cppODE(c(A = "-piecewise(ks*A, time - ts < 0, 0)"),
                  modelname = "pw_freeze", deriv = TRUE, deriv2 = TRUE,
                  compile = FALSE),
  and_bare = pw_cond("time > t1 && time <= t2", "pw_and_bare"),
  and_wrapped = pw_cond("(time > t1) && (time <= t2)", "pw_and_wrapped"),
  or_not = pw_cond("!(time > t1) || time > t2", "pw_or_not"),
  heaviside = cppODE(c(A = "-k * A * Heaviside(A - 0.5)", B = "k * A"),
                     modelname = "heaviside_jac", compile = FALSE),
  step_bdf = pw_step("pw_step_bdf", method = "bdf", deriv = FALSE),
  step_adams = pw_step("pw_step_adams", method = "adams", deriv = FALSE),
  step_fwd = pw_step("pw_step_fwd", method = "bdf", deriv = TRUE),
  step_rev = pw_step("pw_step_rev", method = "bdf", derivMode = "reverse"),
  pulse_bdf = pw_pulse("pw_pulse_bdf", method = "bdf", deriv = FALSE),
  pulse_bdf_s = pw_pulse("pw_pulse_bdf_s", method = "bdf", deriv = TRUE),
  pulse_adams = pw_pulse("pw_pulse_adams", method = "adams", deriv = FALSE),
  pulse_adams_s = pw_pulse("pw_pulse_adams_s", method = "adams", deriv = TRUE),
  pulse_rb4 = pw_pulse("pw_pulse_rb4", method = "rb4", deriv = FALSE),
  pulse_rb4_s = pw_pulse("pw_pulse_rb4_s", method = "rb4", deriv = TRUE),
  pulse_rb4_grid = pw_pulse("pw_pulse_rb4_grid", method = "rb4", deriv = FALSE,
                            useDenseOutput = FALSE),
  ops = cppODE(pw_ops, modelname = "pw_ops", deriv = TRUE, compile = FALSE)
)
if (has_cvode) {
  pulse_eqns <- c(x = "-k*x + c*piecewise(1, time > ts && time <= t2, 0)",
                  y = "k*x - y")
  pw_mod$pulse_cvode <- cvode(pulse_eqns, modelname = "pw_pulse_cvode", compile = FALSE)
  pw_mod$pulse_cvode_s <- cvode(pulse_eqns, modelname = "pw_pulse_cvode_s",
                                deriv = TRUE, compile = FALSE)
  pw_mod$ops_cvode <- cvode(pw_ops, modelname = "pw_ops_cvode", deriv = TRUE,
                            compile = FALSE)
}
do.call(compile, c(unname(pw_mod), list(output = "test_piecewise", cores = test_cores())))

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
  f <- pw_mod$time_switch

  times <- seq(0, 6, by = 0.5)
  p <- c(kf = 0.5, ks = 0.05, ts = 2, A = 1, B = 0)
  out <- solveODE(f, times = times, parms = p, abstol = 1e-10, reltol = 1e-10)

  tk  <- pmin(times, p[["ts"]])          # time spent on the fast branch
  tl  <- pmax(times - p[["ts"]], 0)      # time spent on the slow branch
  A   <- p[["A"]] * exp(-p[["kf"]] * tk - p[["ks"]] * tl)
  expect_equal(unname(out$variable[, "A"]), A, tolerance = 1e-6)
  expect_equal(unname(out$variable[, "B"]), p[["A"]] - A, tolerance = 1e-6)

  expect_equal(unname(out$tangent[, "B", "kf"]), tk * A, tolerance = 1e-6)
  expect_equal(unname(out$tangent[, "B", "ks"]), tl * A, tolerance = 1e-6)
  expect_equal(unname(out$tangent[, "B", "A"]), 1 - A / p[["A"]], tolerance = 1e-6)
})

test_that("a step in time is crossed on a dense output grid", {
  # x and y rest until ts. No step containing the switch passes the error
  # test, however few doubles it spans, and a step that crosses it leaves the
  # right-hand side from before it in the history. Either takes a multistep
  # method without sensitivities to the floor of the step size.
  times <- 0:180
  p <- c(k = 10, c = 1000, ts = 60, a = 1, b = 0, x = 0, y = 0)
  on <- p[["c"]] * pmax(1 - exp(-p[["k"]] * (times - p[["ts"]])), 0)
  exact <- cbind(a = exp(-0.1 * times),
                 b = 2 * (exp(-0.05 * times) - exp(-0.1 * times)),
                 x = on, y = on)
  ref <- solveODE(pw_mod$step_fwd, times, p, abstol = 1e-12, reltol = 1e-10)

  for (m in c("bdf", "adams")) {
    out <- solveODE(pw_mod[[paste0("step_", m)]], times, p,
                    abstol = 1e-12, reltol = 1e-10)
    expect_equal(nrow(out$variable), length(times), info = m)
    expect_equal(unname(out$variable[, colnames(exact)]), unname(exact),
                 tolerance = 1e-8, info = m)
    expect_equal(out$variable, ref$variable, tolerance = 1e-8, info = m)
  }
})

test_that("the reverse mode differentiates across a step in time", {
  times <- 0:180
  p <- c(k = 10, c = 1000, ts = 60, a = 1, b = 0, x = 0.5, y = 0.5)
  fwd <- solveODE(pw_mod$step_fwd, times, p, abstol = 1e-12, reltol = 1e-10)
  set.seed(1)
  W <- array(rnorm(length(fwd$variable)), c(dim(fwd$variable), 1))
  rev <- solveODE(pw_mod$step_rev, times, p, cotangent = W,
                  abstol = 1e-12, reltol = 1e-10)

  expect_equal(rev$variable, fwd$variable, tolerance = 1e-8)
  ref <- apply(fwd$tangent * as.vector(W[, , 1]), 3, sum)
  expect_equal(unname(rev$cotangent[names(ref), 1]), unname(ref),
               tolerance = 1e-6)
})

test_that("a pulse on a state at rest is not stepped over", {
  # Nothing but the switching times at ts and t2 keeps a step from growing
  # across the pulse: without sensitivities from 60 to 90, and for a short one
  # late in the grid with them too. Closed form for k = 1.
  cases <- list(list(ts = 60, t2 = 90, times = 0:180),
                list(ts = 600, t2 = 602, times = c(0, seq(590, 620, 2))))
  for (cs in cases) {
    times <- cs$times
    p <- c(k = 1, c = 1000, ts = cs$ts, t2 = cs$t2, x = 0, y = 0)
    u <- pmin(pmax(times - p[["ts"]], 0), p[["t2"]] - p[["ts"]])
    v <- pmax(times - p[["t2"]], 0)
    x2 <- p[["c"]] * (1 - exp(-u))
    y2 <- p[["c"]] * (1 - exp(-u) - u * exp(-u))
    exact <- cbind(x = x2 * exp(-v), y = (y2 + x2 * v) * exp(-v))
    run <- function(nm)
      solveODE(pw_mod[[nm]], times, p, abstol = 1e-12, reltol = 1e-10)$variable

    for (m in c("bdf", "adams", "rb4")) {
      info <- paste(m, cs$ts)
      plain <- run(paste0("pulse_", m))
      sens <- run(paste0("pulse_", m, "_s"))
      expect_equal(nrow(plain), length(times), info = info)
      expect_equal(unname(plain), unname(exact), tolerance = 1e-8, info = info)
      expect_equal(unname(sens), unname(exact), tolerance = 1e-8, info = info)
      expect_equal(plain, sens, tolerance = 1e-8, info = info)
    }
    # The controlled loop lands on every output time and crosses the switches
    # between them.
    expect_equal(unname(run("pulse_rb4_grid")), unname(exact),
                 tolerance = 1e-8, info = cs$ts)
    # cvode() locates the switching times as roots.
    if (!has_cvode) next
    expect_equal(unname(run("pulse_cvode")), unname(exact), tolerance = 1e-7,
                 info = paste("cvode", cs$ts))
    expect_equal(unname(run("pulse_cvode_s")), unname(exact), tolerance = 1e-7,
                 info = paste("cvode sens", cs$ts))
  }
})

# -- Second order -------------------------------------------------------------

test_that("select propagates value, gradient and Hessian on both branches", {
  # Reference: stats::D() of the branch the condition selects, an independent
  # derivation of the same closed forms.
  nms <- c("a", "b")
  dP  <- diag(2); dimnames(dP) <- list(nms, nms)
  dP2 <- array(0, c(2, 2, 2), dimnames = list(nms, nms, nms))

  f <- pw_mod$d2
  branch <- list(quote(b*a + a^3), quote(a^2*b))
  for (a in c(0.5, 2)) {
    out <- f$evaluate(a = a, b = 3, tangentP = dP, hessianP = dP2, deriv2 = TRUE)
    e <- branch[[1L + (a - 1 > 0)]]
    env <- list(a = a, b = 3)
    H <- matrix(0, 2, 2)
    for (k in 1:2) for (l in 1:2) H[k, l] <- eval(D(D(e, nms[k]), nms[l]), env)
    expect_equal(unname(out$y[1, 1]), eval(e, env))
    expect_equal(unname(out$tangent[1, 1, ]),
                 vapply(nms, function(v) eval(D(e, v), env), 0),
                 tolerance = 1e-12, ignore_attr = TRUE)
    expect_equal(unname(out$hessian[1, 1, , ]), H, tolerance = 1e-12)
  }
})

test_that("a branch that is a literal clears the tangents it replaces", {
  # The second-order materialiser leaves the target's tangent buffers alone
  # when a tree has no dependence, and dxdt is reused across calls, so a
  # select landing on a literal has to stay on the writing path.
  f <- pw_mod$freeze

  times <- seq(0, 4, by = 0.5)
  p <- c(ks = 0.4, ts = 2, A = 1)
  out <- solveODE(f, times = times, parms = p, abstol = 1e-11, reltol = 1e-11)

  tk <- pmin(times, p[["ts"]])          # the decay freezes at ts
  A  <- p[["A"]] * exp(-p[["ks"]] * tk)
  expect_equal(unname(out$variable[, "A"]), A, tolerance = 1e-7)
  expect_equal(unname(out$tangent[, "A", "ks"]), -tk * A, tolerance = 1e-7)
  expect_equal(unname(out$tangent[, "A", "A"]), A / p[["A"]], tolerance = 1e-7)
  expect_equal(unname(out$hessian[, "A", "ks", "ks"]), tk^2 * A, tolerance = 1e-7)
  expect_equal(unname(out$hessian[, "A", "ks", "A"]), -tk * A / p[["A"]],
               tolerance = 1e-7)
})

# -- Conditions --------------------------------------------------------------

test_that("the R and C spellings of the logical operators parse", {
  # && and || bind below the comparisons while & and | bind above them, so the
  # grouping has to survive the rewrite even without the parentheses.
  times <- seq(0, 6, by = 0.5)
  p <- c(k = 0.4, t1 = 1, t2 = 3, A = 1)
  run <- function(mod)
    solveODE(mod, times = times, parms = p, abstol = 1e-11,
             reltol = 1e-11)$variable[, "A"]

  bare    <- run(pw_mod$and_bare)
  wrapped <- run(pw_mod$and_wrapped)
  expect_equal(unname(bare), unname(wrapped))

  # Decay runs between t1 and t2 only, so A is flat on either side.
  tk <- pmin(pmax(times - p[["t1"]], 0), p[["t2"]] - p[["t1"]])
  expect_equal(unname(bare), p[["A"]] * exp(-p[["k"]] * tk), tolerance = 1e-7)

  # The complement, spelled with || and !, has to give the mirror image.
  outside <- run(pw_mod$or_not)
  expect_equal(unname(outside),
               p[["A"]] * exp(-p[["k"]] * (times - tk)), tolerance = 1e-7)
})

test_that("every operator reads as in R, on both backends", {
  # The reference evaluates the same strings in R and integrates them piece by
  # piece between the switching points; x equals time. Its finite differences
  # in t1 and t2 take the jumps at the switching times the tangents carry.
  times <- seq(0, 4, by = 0.25)
  pars <- c(a = 2, b = 3, t1 = 1.3, t2 = 2.6, h = 0.5, mode = 1)
  p <- c(setNames(rep(0, length(pw_ops)), names(pw_ops)), pars)
  env <- list(piecewise = function(v, cond, otherwise) if (cond) v else otherwise,
              Heaviside = function(z) if (z > 0) 1 else if (z == 0) 0.5 else 0)
  reference <- function(pars) {
    cuts <- sort(unique(c(times, pars[c("h", "t1", "t2")])))
    sapply(setdiff(names(pw_ops), "x"), function(nm) {
      f <- function(s) vapply(s, function(si) eval(str2lang(pw_ops[[nm]]),
        c(list(x = si, time = si), as.list(pars), env)), 0)
      piece <- vapply(seq_len(length(cuts) - 1), function(i)
        integrate(f, cuts[i], cuts[i + 1], rel.tol = 1e-12)$value, 0)
      cumsum(c(0, piece))[match(times, cuts)]
    })
  }
  ref <- reference(pars)
  h <- 1e-5
  dref <- lapply(c(t1 = "t1", t2 = "t2"), function(n)
    (reference(replace(pars, n, pars[[n]] + h)) -
       reference(replace(pars, n, pars[[n]] - h))) / (2 * h))

  for (m in c("ops", if (has_cvode) "ops_cvode")) {
    out <- solveODE(pw_mod[[m]], times, p, abstol = 1e-11, reltol = 1e-11)
    expect_equal(out$variable[, colnames(ref)], ref, tolerance = 1e-7, info = m)
    for (n in names(dref))
      expect_equal(unname(out$tangent[, colnames(ref), n]), unname(dref[[n]]),
                   tolerance = 1e-5, info = paste(m, n))
  }
})

test_that("an expression that does not parse names itself", {
  # The reason has to survive the trip through reticulate, which truncates a
  # long message and then indexes it with an offset from the untruncated one.
  long <- paste0("-k*A ** * A", strrep(" + 0*A", 200))
  expect_error(cppODE(c(A = long), modelname = "pw_unparseable", compile = FALSE),
               "cannot parse expression")
})

test_that("Heaviside survives differentiation", {
  # d/dx Heaviside(x) is DiracDelta, which no printer knows. A discrete model
  # cannot mean an impulse of infinite height, so it is emitted as one at the
  # switching point and zero either side, the way cppFUN already takes it.
  mod <- pw_mod$heaviside
  tt <- seq(0, 2, 0.5)
  res <- solveODE(mod, tt, c(A = 1, B = 0, k = 0.7),
                  abstol = 1e-10, reltol = 1e-10)

  # The decay switches off at 0.5 and the state holds there.
  expect_equal(unname(res$variable[nrow(res$variable), 1]), 0.5, tolerance = 1e-3)
  expect_false(anyNA(res$tangent))
})
