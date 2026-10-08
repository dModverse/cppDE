# What solveODE(), solveODEBatch(), prepareBatch() and solveBatch() reject,
# and with which message. Every check runs before the compiled entry point is
# looked up, so the models are generated and never compiled.

skip_on_cran()
skip_if_not(codegenAvailable())

fwd <- cppODE(c(A = "-k * A * u"), forcings = "u", deriv2 = TRUE,
              modelname = "chk_fwd", compile = FALSE)
plain <- cppODE(c(A = "-k * A"), deriv = FALSE, modelname = "chk_plain", compile = FALSE)
rev <- cppODE(c(A = "-k * A"), derivMode = "reverse", modelname = "chk_rev", compile = FALSE)
fr  <- cppODE(c(A = "-k * A"), derivMode = "forward-reverse", modelname = "chk_fr",
              compile = FALSE)

tt <- c(0, 1)
p  <- c(A = 1, k = 0.5)
u  <- list(u = data.frame(time = 0, value = 1))
W  <- matrix(1, 2, 1)

test_that("solveODE() checks the times, parameters and forcings it is given", {
  expect_error(solveODE(fwd, numeric(0), p, forcings = u), "'times' must be a non-empty")
  expect_error(solveODE(fwd, c(0, NA), p, forcings = u), "'times' must be a non-empty")
  expect_error(solveODE(fwd, "1", p, forcings = u), "'times' must be a non-empty")
  expect_error(solveODE(fwd, tt, unname(p), forcings = u), "'parms' must be a named numeric")
  expect_error(solveODE(fwd, tt, p["A"], forcings = u), "'parms' missing: k")
  expect_error(solveODE(fwd, tt, c(A = 1, k = Inf), forcings = u), "'parms' must be finite")
  expect_error(solveODE(fwd, tt, p), "Model requires forcings: u")
  expect_error(solveODE(fwd, tt, p, forcings = unname(u)), "'forcings' must be a named list")
  expect_error(solveODE(fwd, tt, p, forcings = list(v = u[[1]])), "Missing forcings: u")
  expect_error(solveODE(fwd, tt, p, forcings = list(u = data.frame(time = c(1, 1), value = 1))),
               "Forcing 'u': duplicate times")
  expect_error(solveODE(structure("x", variables = "A"), tt, p),
               "'model' is missing attributes")
  expect_error(solveODE(fwd, tt, p, forcings = u), "Model not loaded. Run compile\\(\\) first.")
})

test_that("solveODE() checks its solver options", {
  chk <- function(...) solveODE(fwd, tt, p, forcings = u, ...)
  expect_error(chk(abstol = 0), "'abstol' must be positive")
  expect_error(chk(reltol = -1), "'reltol' must be positive")
  expect_error(chk(hini = -1), "'hini' must be non-negative")
  expect_error(chk(roottol = 0), "'roottol' must be positive")
  expect_error(chk(maxattempts = 0), "'maxattempts' must be positive")
  expect_error(chk(maxsteps = 0), "'maxsteps' must be positive")
  expect_error(chk(maxroot = 0), "'maxroot' must be positive")
  expect_error(chk(sensErrCon = NA), "'sensErrCon' must be TRUE or FALSE")
  expect_error(chk(onFailure = "ignore"), "should be one of")
  expect_error(solveODE(plain, tt, p, sensErrCon = FALSE), "computes no sensitivities")
})

test_that("solveODE() checks tangent, hessian and fixed against the model", {
  chk <- function(...) solveODE(fwd, tt, p, forcings = u, ...)
  expect_error(solveODE(plain, tt, p, tangent = diag(2)), "model has deriv = FALSE")
  expect_error(solveODE(rev, tt, p, cotangent = W, hessian = array(0, c(2, 1, 1))),
               "model has deriv2 = FALSE")
  expect_error(chk(tangent = "a"), "'tangent' must be numeric")
  expect_error(chk(tangent = 1:3), "'tangent' must have length")
  expect_error(chk(tangent = matrix(0, 1, 2, dimnames = list(NULL, c("a", "b")))),
               "column names must match: A, k")
  expect_error(chk(tangent = matrix(0, 2, 2, dimnames = list(c("A", "x"), NULL))),
               "row names must be c\\(variables, parameters\\)")
  expect_error(chk(tangent = matrix(0, 1, 2, dimnames = list("x", NULL))),
               "unknown row names: x")
  expect_error(chk(tangent = matrix(0, 3, 1, dimnames = list(c("k", "k", "A"), NULL))),
               "duplicate row names")
  expect_error(chk(fixed = "z"), "Unknown 'fixed' names: z")
  expect_error(chk(fixed = 1), "'fixed' must be a character vector")
  expect_warning(expect_error(solveODE(plain, tt, p, fixed = "k"), "Model not loaded"),
                 "'fixed' ignored")
  expect_error(chk(hessian = "a"), "'hessian' must be numeric")
  expect_error(chk(hessian = array(0, c(2, 1, 1))), "dim 2 and 3 equal to 2")
  expect_error(chk(hessian = 1:3), "'hessian' must have length 8")
  expect_error(chk(hessian = array(0, c(1, 2, 2), list("x", NULL, NULL))),
               "unknown dim-1 names: x")
})

test_that("solveODE() keeps the reverse arguments to the reverse modes", {
  expect_error(solveODE(plain, tt, p, cotangent = W), "not compiled with derivMode")
  expect_error(solveODE(rev, tt, p), "needs a 'cotangent'")
  expect_error(solveODE(rev, tt, p, cotangent = "w"), "'cotangent' must be numeric")
  expect_error(solveODE(rev, tt, p, cotangent = 1:2), "\\[n_out, n_states\\] matrix")
  expect_error(solveODE(rev, tt, p, cotangent = matrix(1, 2, 3)), "3 state columns")
  expect_error(solveODE(rev, tt, p, cotangent = W, curvature = array(0, c(2, 1, 1, 2))),
               "needs derivMode = \"forward-reverse\"")
  expect_error(solveODE(plain, tt, p, curvature = array(0, c(2, 1, 1, 2))),
               "derivative of a cotangent and needs one")
  expect_error(solveODE(fr, tt, p, cotangent = W, curvature = array(0, c(2, 1, 2, 2))),
               "cotangent's own first three dimensions")
  expect_error(solveODE(fr, tt, p, cotangent = W, curvature = "c"),
               "'curvature' must be numeric")
  expect_error(solveODE(plain, tt, p, keepStore = TRUE), "belongs to a model compiled")
  expect_error(solveODE(rev, tt, p, cotangent = W, store = "s"),
               "'store' must be the `store` element")
  expect_error(solveODE(fr, tt, p, keepStore = TRUE), "forward-reverse")
  expect_error(solveODE(rev, tt, p, cotangent = W, adjoint = list(trace = TRUE)),
               "must come from adjointControl")
  expect_error(solveODE(fwd, tt, p, forcings = u, adjoint = adjointControl()),
               "backward pass of a model compiled")
  expect_error(solveODE(fr, tt, p, cotangent = W, adjoint = adjointControl(refine = TRUE)),
               "backward pass of a model compiled")
})

test_that("the CVODE backend refuses what it cannot do", {
  skip_if_not(isTRUE(cvodeConfig$available), "CVODE backend not available")
  cf <- cvode(c(A = "-k * A"), deriv = TRUE, modelname = "chk_cv", compile = FALSE)
  cr <- cvode(c(A = "-k * A"), derivMode = "reverse", modelname = "chk_cv_rev",
              compile = FALSE)
  expect_error(solveODE(cf, tt, p, sensErrCon = FALSE), "not available on the CVODE")
  expect_error(solveODE(cr, tt, p, keepStore = TRUE), "not available on the CVODE")
  expect_error(solveODE(cr, tt, p, cotangent = W, adjoint = adjointControl(trace = TRUE)),
               "'gradtol' alone")
  expect_error(solveODE(cf, tt, p), "Model not loaded")
})

test_that("solveODEBatch() and prepareBatch() check the conditions", {
  cs <- list(a = list(parms = p), b = list(parms = p * 2))
  for (f in list(solveODEBatch, prepareBatch)) {
    expect_error(f(plain, list(), times = tt), "'conditions' must be a non-empty list")
    expect_error(f(plain, "a", times = tt), "'conditions' must be a non-empty list")
    expect_error(f(plain, list(p), times = tt), "every element of 'conditions' must be a list")
    expect_error(f(plain, list(list(parms = p, tol = 1)), times = tt),
                 "unknown per-condition argument\\(s\\): tol")
    expect_error(f(plain, list(list(parms = p))), "condition 1 has no 'times'")
    expect_error(f(plain, list(list(times = tt))), "condition 1 has no 'times' or no 'parms'")
    expect_error(f(plain, list(a = list(parms = p), b = list(parms = p["A"])), times = tt),
                 "'parms' missing: k")
  }
  expect_error(solveODEBatch(plain, cs, times = tt), "Model not loaded")
})

test_that("solveBatch() checks what replaces the prepared inputs", {
  h <- prepareBatch(plain, list(a = list(parms = p), b = list(parms = p)), times = tt)
  expect_s3_class(h, "cppDEbatch")
  expect_identical(h$names, c("a", "b"))
  expect_identical(h$required, c("A", "k"))
  expect_error(solveBatch(list()), "'handle' must come from prepareBatch")
  expect_error(solveBatch(h, parms = list(p)), "'parms' must be a list with one element")
  expect_error(solveBatch(h, parms = list(p, unname(p))), "'parms' elements must be named")
  expect_error(solveBatch(h, parms = list(p, p["A"])), "condition 2: 'parms' missing k")
  expect_error(solveBatch(h, parms = list(p, c(A = NaN, k = 1))),
               "condition 2: 'parms' must be finite")
  expect_error(solveBatch(h, store = list(NULL, NULL)), "belongs to a model compiled")
  expect_error(solveBatch(h, cotangent = list(W, W)), "prepared without a 'cotangent'")
  expect_error(solveBatch(h, onFailure = "x"), "should be one of")
  expect_error(solveBatch(h), "Model not loaded")

  hf <- prepareBatch(fwd, list(list(parms = p, forcings = u)), times = tt)
  expect_error(solveBatch(hf, tangent = list(diag(3))), "tangent shape changed")

  hr <- prepareBatch(rev, list(list(parms = p, cotangent = W)), times = tt)
  expect_error(solveBatch(hr, cotangent = list(matrix(1, 3, 1))), "cotangent shape changed")
  expect_error(solveBatch(hr, store = list("s")), "condition 1: 'store' must be the")
  expect_error(solveBatch(hr, curvature = list(array(0, c(2, 1, 1, 2)))),
               "prepared without a 'curvature'")

  hc <- prepareBatch(fr, list(list(parms = p, cotangent = W,
                                   curvature = array(0, c(2, 1, 1, 2)))), times = tt)
  expect_error(solveBatch(hc, curvature = list(array(0, c(2, 1, 1, 3)))),
               "curvature shape changed")
})

test_that("a trace file needs a path per condition or one template", {
  expect_error(cppDE:::.batchTraceFiles("", 2L, NULL),
               "'traceFile' must be a non-empty character vector")
  expect_identical(cppDE:::.batchTraceFiles(file.path("d", "t.csv"), 2L, c("a b", "c")),
                   list(file.path("d", "t_a.b.csv"), file.path("d", "t_c.csv")))
  expect_identical(cppDE:::.batchTraceFiles("t", 2L, NULL), list("t_1", "t_2"))
})
