# Unit tests for cppde::dual2nd math primitives. One primitive per model,
# asserted equal between derivMode "dual" and "symbolic": independent code
# paths, with the symbolic one as the trusted oracle.

skip_on_cran()

library(cppDE)

# Helper: turn an expression into a model name. Operators are spelled out, so
# a + b, a - b, a * b, a / b and a^b do not collapse onto one name and make
# unique_modelname() warn its way through a_b_2, a_b_3.
name_of <- function(expr) {
  s <- paste(expr, collapse = "_")
  ops <- c("+" = "_add_", "-" = "_sub_", "*" = "_mul_",
           "/" = "_div_", "^" = "_pow_")
  for (op in names(ops)) s <- gsub(op, ops[[op]], s, fixed = TRUE)
  gsub("^_|_$", "", gsub("[^A-Za-z0-9]+", "_", s))
}

# Helper: run both modes and return (y, dy, d2y) arrays.
run_modes <- function(expr, parameters, x_vals, dP, dP2 = NULL) {
  out <- list()
  for (mode in c("dual", "symbolic")) {
    f <- funCpp(expr, parameters = parameters,
                deriv = TRUE, deriv2 = TRUE, derivMode = mode,
                compile = TRUE,
                modelname = paste0("d2prim_", mode, "_", name_of(expr)))
    args <- as.list(x_vals)
    args$dP <- dP
    if (!is.null(dP2)) args$dP2 <- dP2
    args$deriv2 <- TRUE
    out[[mode]] <- do.call(f$evaluate, args)
  }
  out
}

# Helper: assert dual and symbolic agree on all three derivative levels.
expect_modes_agree <- function(out, tol_y = 1e-12, tol_dy = 1e-10, tol_d2y = 1e-10) {
  expect_equal(unname(out$dual$y),   unname(out$symbolic$y),   tolerance = tol_y)
  expect_equal(unname(out$dual$dy),  unname(out$symbolic$dy),  tolerance = tol_dy)
  expect_equal(unname(out$dual$d2y), unname(out$symbolic$d2y), tolerance = tol_d2y)
}

# Identity Phi(theta) = theta seed: dP = I, dP2 = 0. Each parameter is its
# own theta direction so the raw Jacobian / Hessian come through.
identity_seeds <- function(par_names) {
  k <- length(par_names)
  dP  <- diag(k); dimnames(dP)  <- list(par_names, par_names)
  dP2 <- array(0, c(k, k, k), dimnames = list(par_names, par_names, par_names))
  list(dP = dP, dP2 = dP2)
}

# -- Binary arithmetic --------------------------------------------------------

test_that("dual2nd matches symbolic on a + b", {
  s <- identity_seeds(c("a", "b"))
  out <- run_modes(c(y = "a + b"), parameters = c("a", "b"),
                   x_vals = list(a = 1.5, b = 2.5), dP = s$dP, dP2 = s$dP2)
  expect_modes_agree(out)
})

test_that("dual2nd matches symbolic on a * b (cross-Hessian)", {
  s <- identity_seeds(c("a", "b"))
  out <- run_modes(c(y = "a * b"), parameters = c("a", "b"),
                   x_vals = list(a = 1.5, b = 2.5), dP = s$dP, dP2 = s$dP2)
  expect_modes_agree(out)
})

test_that("dual2nd matches symbolic on a / b", {
  s <- identity_seeds(c("a", "b"))
  out <- run_modes(c(y = "a / b"), parameters = c("a", "b"),
                   x_vals = list(a = 1.5, b = 2.5), dP = s$dP, dP2 = s$dP2)
  expect_modes_agree(out)
})

test_that("dual2nd matches symbolic on a - b", {
  s <- identity_seeds(c("a", "b"))
  out <- run_modes(c(y = "a - b"), parameters = c("a", "b"),
                   x_vals = list(a = 1.5, b = 2.5), dP = s$dP, dP2 = s$dP2)
  expect_modes_agree(out)
})

# -- Transcendentals ----------------------------------------------------------

test_that("dual2nd matches symbolic on sin / cos / tan", {
  s <- identity_seeds(c("x"))
  for (fn in c("sin", "cos", "tan")) {
    expr <- setNames(sprintf("%s(x)", fn), "y")
    out  <- run_modes(expr, parameters = "x",
                      x_vals = list(x = 0.7), dP = s$dP, dP2 = s$dP2)
    expect_modes_agree(out)
  }
})

test_that("dual2nd matches symbolic on exp / log / sqrt", {
  s <- identity_seeds(c("x"))
  for (fn in c("exp", "log", "sqrt")) {
    expr <- setNames(sprintf("%s(x)", fn), "y")
    out  <- run_modes(expr, parameters = "x",
                      x_vals = list(x = 1.7), dP = s$dP, dP2 = s$dP2)
    expect_modes_agree(out)
  }
})

test_that("dual2nd matches symbolic on hyperbolic trig", {
  s <- identity_seeds(c("x"))
  for (fn in c("sinh", "cosh", "tanh")) {
    expr <- setNames(sprintf("%s(x)", fn), "y")
    out  <- run_modes(expr, parameters = "x",
                      x_vals = list(x = 0.5), dP = s$dP, dP2 = s$dP2)
    expect_modes_agree(out)
  }
})

# -- pow ----------------------------------------------------------------------

test_that("dual2nd matches symbolic on a^b (both AD)", {
  s <- identity_seeds(c("a", "b"))
  out <- run_modes(c(y = "a^b"), parameters = c("a", "b"),
                   x_vals = list(a = 1.7, b = 2.3), dP = s$dP, dP2 = s$dP2)
  expect_modes_agree(out)
})

test_that("dual2nd matches symbolic on a^2 (scalar exponent)", {
  s <- identity_seeds(c("a"))
  out <- run_modes(c(y = "a^2"), parameters = c("a"),
                   x_vals = list(a = 1.7), dP = s$dP, dP2 = s$dP2)
  expect_modes_agree(out)
})

# -- Composite ----------------------------------------------------------------

test_that("dual2nd matches symbolic on a*sin(b) + exp(a)", {
  s <- identity_seeds(c("a", "b"))
  out <- run_modes(c(y = "a*sin(b) + exp(a)"), parameters = c("a", "b"),
                   x_vals = list(a = 0.4, b = 1.1), dP = s$dP, dP2 = s$dP2)
  expect_modes_agree(out)
})

# -- Hessian symmetry exposed in output ---------------------------------------

test_that("dual2nd output is Hessian-symmetric (mirror via dd_raw)", {
  s <- identity_seeds(c("a", "b", "c"))
  out <- run_modes(c(y = "a*b + b*c + a*c"), parameters = c("a", "b", "c"),
                   x_vals = list(a = 1, b = 2, c = 3), dP = s$dP, dP2 = s$dP2)
  d2y <- out$dual$d2y
  for (i in seq_len(dim(d2y)[3]))
    for (j in seq_len(dim(d2y)[4]))
      expect_equal(d2y[, , i, j], d2y[, , j, i])
})
