# Test cppFUN() algebraic function compilation with sensitivities.

skip_on_cran()

# Closed forms of y1 = a*x^2 + b, y2 = sin(c*x) over (a, b, c, x).
.xs_jac <- function(a, b, c, x) {
  rbind(y1 = c(a = x^2, b = 1, c = 0, x = 2 * a * x),
        y2 = c(a = 0, b = 0, c = x * cos(c * x), x = c * cos(c * x)))
}
.xs_hess <- function(a, b, c, x) {
  s <- c("a", "b", "c", "x")
  H <- array(0, c(2, 4, 4), list(c("y1", "y2"), s, s))
  H["y1", "a", "x"] <- H["y1", "x", "a"] <- 2 * x
  H["y1", "x", "x"] <- 2 * a
  H["y2", "c", "c"] <- -x^2 * sin(c * x)
  H["y2", "c", "x"] <- H["y2", "x", "c"] <- cos(c * x) - c * x * sin(c * x)
  H["y2", "x", "x"] <- -c^2 * sin(c * x)
  H
}

# Every object the tests below evaluate, built once and linked into one shared
# object. The tests on uncompiled objects and on rejected names build their own.
fx_trafo  <- c(y1 = "a * exp(-b * x)", y2 = "a + b * x")
vjp_trafo <- c(y1 = "a * exp(-k * t) + b", y2 = "log(a + k * k) * t")
collide_orders <- list(c("p", "x", "y", "k", "x_obs", "y_local"),
                       c("y_local", "x_obs", "k", "y", "x", "p"),
                       c("k", "p", "y_local", "x", "y", "x_obs"))

fx <- list(
  basic  = cppFUN(fx_trafo, parameters = c("a", "b", "x"), deriv = TRUE,
                  modelname = "fun_basic", convenient = TRUE),
  tokens = cppFUN(c(o1 = "std * exp(ini * log(10)) + default",
                    o2 = "sqrt(int) + std^2"),
                  variables = "int", parameters = c("std", "ini", "default"),
                  deriv = TRUE, derivMode = c("forward", "reverse"),
                  modelname = "cxx_tokens", convenient = TRUE),
  jac    = cppFUN(fx_trafo, parameters = c("a", "b", "x"),
                  deriv = TRUE, derivMode = "forward",
                  modelname = "fun_jac", convenient = TRUE),
  hess   = cppFUN(c(y = "a * b * x^2"), parameters = c("a", "b", "x"),
                  deriv = TRUE, deriv2 = TRUE, derivMode = "forward",
                  modelname = "fun_hess", convenient = TRUE),
  fixed  = cppFUN(c(y = "a * b + c"), parameters = c("a", "b", "c"),
                  fixed = "c", deriv = TRUE, derivMode = "forward",
                  modelname = "fun_fixed", convenient = TRUE),
  xs     = cppFUN(c(y1 = "a*x^2 + b", y2 = "sin(c*x)"),
                  parameters = c("a", "b", "c", "x"),
                  deriv = TRUE, deriv2 = TRUE, derivMode = "forward",
                  modelname = "xs_d2"),
  passthru = cppFUN(c(la = "la", y2 = "la^2 + b", zero = "0"),
                    parameters = c("la", "b"),
                    deriv = TRUE, deriv2 = TRUE, derivMode = "forward",
                    modelname = "xs_passthru"),
  vjp_rev = cppFUN(vjp_trafo, variables = "t", parameters = c("a", "b", "k"),
                   deriv = TRUE, derivMode = "reverse", modelname = "vjp_rev",
                   convenient = FALSE),
  vjp_fwd = cppFUN(vjp_trafo, variables = "t", parameters = c("a", "b", "k"),
                   deriv = TRUE, derivMode = "forward", modelname = "vjp_fwd",
                   convenient = FALSE),
  vjp_seeds = cppFUN(c(y1 = "a * b", y2 = "sin(a) + b * b"), variables = NULL,
                     parameters = c("a", "b"), deriv = TRUE,
                     derivMode = "reverse", modelname = "vjp_seeds",
                     convenient = FALSE),
  dm_fwd  = cppFUN(c(y = "a * a"), parameters = "a", deriv = TRUE,
                   derivMode = "forward", modelname = "dm_fwd",
                   convenient = FALSE),
  dm_rev  = cppFUN(c(y = "a * a"), parameters = "a", deriv = TRUE,
                   derivMode = "reverse", modelname = "dm_rev",
                   convenient = FALSE),
  dm_both = cppFUN(c(y = "a * a"), parameters = "a", deriv = TRUE,
                   derivMode = c("forward", "reverse"), modelname = "dm_both",
                   convenient = FALSE),
  vjp_fr = cppFUN(c(y1 = "a * x1^2 + b * x1 * x2", y2 = "sin(a * x2) + b^2 * x1"),
                  variables = c("x1", "x2"), parameters = c("a", "b"),
                  deriv2 = TRUE, derivMode = c("forward", "forward-reverse"),
                  modelname = "cf_vjpfr")
)
fx_collide <- lapply(seq_along(collide_orders), function(i)
  cppFUN(c(o1 = "p + 2 * x", o2 = "y * k",
           o3 = "x_obs + y_local", o4 = "p * y_local - k"),
         parameters = collide_orders[[i]], deriv = TRUE, derivMode = "forward",
         modelname = paste0("collide_", i), convenient = TRUE))
do.call(compile, c(unname(fx), fx_collide,
                   list(output = "test_cppFUN", cores = 1)))

# -- Basic cppFUN output structure ---------------------------------------------

test_that("cppFUN returns correct output", {
  f <- fx$basic
  res <- f$func(a = 2, b = 0.5, x = 1)

  expect_true(is.matrix(res))
  expect_equal(colnames(res), c("y1", "y2"))
  expect_equal(unname(res[1, "y1"]), 2 * exp(-0.5), tolerance = 1e-10)
  expect_equal(unname(res[1, "y2"]), 2.5, tolerance = 1e-10)
})

test_that("an object that was not compiled says so", {
  # Every entry needs compiled code.
  f <- cppFUN(c(y = "a * x"), variables = "x", parameters = "a",
              derivMode = c("forward", "forward-reverse"), modelname = "fun_uncompiled",
              convenient = FALSE)
  M <- matrix(2, 1, 1, dimnames = list(NULL, "x"))
  msg <- "'fun_uncompiled' is not compiled; call compile\\(\\)"
  expect_error(f$func(M, c(a = 3)), msg)
  expect_error(f$jac(M, c(a = 3)), msg)
  expect_error(f$evaluate(M, c(a = 3)), msg)
  expect_error(f$vjp(M, c(a = 3), matrix(1, 1, 1)), msg)
  expect_error(f$vjp(M, c(a = 3), matrix(1, 1, 1), tangentP = matrix(1, 1, 1)),
               msg)
})

test_that("parameters named like the generated arrays do not collide", {
  # The emitted code indexes p[] and x_obs[]. A parameter of that name is
  # substituted for its own slot, in whatever order the parameters are listed.
  pars <- c(p = 2, x = 3, y = 5, k = 7, x_obs = 11, y_local = 13)

  for (i in seq_along(collide_orders)) {
    f <- fx_collide[[i]]
    res <- do.call(f$func, as.list(pars[collide_orders[[i]]]))
    jac <- do.call(f$jac, as.list(pars[collide_orders[[i]]]))[1, , ]

    expect_equal(unname(res[1, ]), c(8, 35, 24, 19), tolerance = 1e-10,
                 label = paste("values, order", i))
    expect_equal(unname(jac["o1", "x"]), 2, tolerance = 1e-10)
    expect_equal(unname(jac["o4", "p"]), 13, tolerance = 1e-10)
    expect_equal(unname(jac["o4", "y_local"]), 2, tolerance = 1e-10)
  }
})

test_that("symbols named after C++ tokens do not reach the generated source", {
  # The printer writes std::pow and spells a reserved word default_. A symbol
  # of either name is substituted for its slot, so neither reaches the source.
  f <- fx$tokens
  res <- f$func(int = 4, std = 2, ini = 0.5, default = 3)
  jac <- f$jac(int = 4, std = 2, ini = 0.5, default = 3)[1, , ]

  expect_equal(unname(res[1, ]), c(2 * 10^0.5 + 3, 6), tolerance = 1e-10)
  expect_equal(unname(jac["o1", "std"]), 10^0.5, tolerance = 1e-8)
  expect_equal(unname(jac["o1", "default"]), 1, tolerance = 1e-10)
  expect_equal(unname(jac["o2", "int"]), 0.25, tolerance = 1e-8)

  r <- f$vjp(matrix(4, 1, 1, dimnames = list(NULL, "int")),
             c(std = 2, ini = 0.5, default = 3), matrix(c(0, 1), 1, 2))
  expect_equal(unname(r$cotangentX[1, 1, 1]), 0.25, tolerance = 1e-10)
  expect_equal(unname(r$cotangentP["std", 1]), 4, tolerance = 1e-10)
})

test_that("a Python keyword as a symbol name is rejected", {
  # The generator parses through Python, where the name is a syntax error, and
  # True, False and None read as constants and drop the symbol from the model.
  expect_error(cppFUN(c(y = "class * 2"), modelname = "py_kw_eqn"),
               "Python keyword used as a symbol name: 'class'")
  expect_error(cppFUN(c(y = "a * 2"), parameters = c("a", "lambda"),
                      modelname = "py_kw_par"),
               "'lambda'")
  expect_error(cppFUN(c(lambda = "a * 2"), parameters = "a",
                      modelname = "py_kw_out"),
               "'lambda'")
  expect_error(cppFUN(c(y = "a * True"), parameters = c("a", "True"),
                      modelname = "py_kw_const"),
               "'True'")
})

# -- Jacobian correctness -----------------------------------------------------

test_that("cppFUN Jacobian matches analytical derivatives", {
  f <- fx$jac
  jac <- f$jac(a = 2, b = 0.5, x = 1)

  # jac is [obs, outputs, params] array
  expect_equal(length(dim(jac)), 3)
  j <- jac[1, , ]  # first (only) observation

  # dy1/da = exp(-b*x) = exp(-0.5)
  expect_equal(j["y1", "a"], exp(-0.5), tolerance = 1e-8)
  # dy1/db = -a*x*exp(-b*x) = -2*1*exp(-0.5)
  expect_equal(j["y1", "b"], -2 * exp(-0.5), tolerance = 1e-8)
  # dy1/dx = -a*b*exp(-b*x) = -2*0.5*exp(-0.5)
  expect_equal(j["y1", "x"], -1 * exp(-0.5), tolerance = 1e-8)
  # dy2/da = 1
  expect_equal(j["y2", "a"], 1, tolerance = 1e-10)
  # dy2/db = x = 1
  expect_equal(j["y2", "b"], 1, tolerance = 1e-10)
  # dy2/dx = b = 0.5
  expect_equal(j["y2", "x"], 0.5, tolerance = 1e-10)
})

# -- Hessian structure ---------------------------------------------------------

test_that("cppFUN Hessian has correct dimensions and is symmetric", {
  f <- fx$hess
  hess_arr <- f$hess(a = 2, b = 3, x = 4)

  expect_true(!is.null(hess_arr))
  # [obs, outputs, params, params]
  expect_equal(dim(hess_arr)[2], 1)   # 1 output
  expect_equal(dim(hess_arr)[3], 3)   # 3 params
  expect_equal(dim(hess_arr)[4], 3)   # 3 params

  hess <- hess_arr[1, 1, , ]
  # Hessian should be symmetric
  expect_equal(hess, t(hess), tolerance = 1e-10)

  # d2y/da db = x^2 = 16
  expect_equal(hess["a", "b"], 16, tolerance = 1e-8)
  # d2y/da dx = 2*b*x = 24
  expect_equal(hess["a", "x"], 24, tolerance = 1e-8)
  # d2y/db dx = 2*a*x = 16
  expect_equal(hess["b", "x"], 16, tolerance = 1e-8)
})

# -- Fixed parameters in cppFUN ------------------------------------------------

test_that("cppFUN fixed parameters are excluded from derivatives", {
  f <- fx$fixed
  jac <- f$jac(a = 2, b = 3, c = 1)

  # Only 2 params in Jacobian (a, b), not c
  expect_equal(dim(jac)[3], 2)
  jac_names <- dimnames(jac)[[3]]
  expect_false("c" %in% jac_names)
  expect_equal(unname(jac[1, "y", ]), c(3, 2), tolerance = 1e-12)
})

# -- raw derivatives against closed forms ---------------------------------------

test_that("cppFUN raw jac/hess match the closed forms", {
  pars  <- list(a = 2, b = 1, c = 0.3, x = 3)
  f <- fx$xs
  y    <- do.call(f$func, pars)
  jac  <- do.call(f$jac,  pars)
  hess <- do.call(f$hess, pars)
  expect_equal(unname(y[1, ]), c(2 * 9 + 1, sin(0.9)), tolerance = 1e-12)
  expect_equal(unname(jac[1, , ]), unname(do.call(.xs_jac, pars)),
               tolerance = 1e-12)
  expect_equal(unname(hess[1, , , ]), unname(do.call(.xs_hess, pars)),
               tolerance = 1e-12)
})

# -- second-order chain rule, including hessianX/hessianP -----------------------

test_that("cppFUN evaluates the second-order chain rule", {
  th    <- c("th1", "th2", "th3")
  s     <- c("a", "b", "c", "x")
  # Linear part of Phi: theta -> (a, b, c) = (2*th1, 1*th2, 1*th3 + th1)
  dP <- matrix(0, 4, length(th), dimnames = list(s, th))
  dP["a", "th1"] <- 2
  dP["b", "th2"] <- 1
  dP["c", "th3"] <- 1
  dP["c", "th1"] <- 1
  # Nonlinear quadratic part: a depends on th1*th2 with coefficient 0.5.
  dP2 <- array(0, c(4, length(th), length(th)), dimnames = list(s, th, th))
  dP2["a", "th1", "th2"] <- 0.5
  dP2["a", "th2", "th1"] <- 0.5
  pars <- list(a = 2, b = 1, c = 0.3, x = 3)
  f <- fx$xs
  out <- do.call(f$evaluate,
                 c(pars, list(tangentP = dP, hessianP = dP2, deriv2 = TRUE)))

  J <- do.call(.xs_jac, pars)
  H <- do.call(.xs_hess, pars)
  d2 <- array(0, c(2, 3, 3))
  for (o in 1:2) {
    d2[o, , ] <- t(dP) %*% H[o, , ] %*% dP
    for (i in 1:4) d2[o, , ] <- d2[o, , ] + J[o, i] * dP2[i, , ]
  }
  expect_equal(unname(out$tangent[1, , ]), unname(J %*% dP), tolerance = 1e-12)
  expect_equal(unname(out$hessian[1, , , ]), d2, tolerance = 1e-12)
})

# -- forward + deriv2 + identity pass-through (regression) ------------------------
# An output that is just a bare param or var arms no inner Hessian seed, so the
# dual2nd writeback has to reach it through the bounds-safe const overload.

test_that("cppFUN forward deriv2 handles identity pass-through", {
  pars  <- list(la = 1.5, b = 0.7)
  f <- fx$passthru
  d <- do.call(f$evaluate, c(pars, list(deriv2 = TRUE)))
  H <- array(0, c(3, 2, 2))
  H[2, 1, 1] <- 2
  expect_equal(unname(d$y[1, ]), c(1.5, 1.5^2 + 0.7, 0), tolerance = 1e-12)
  expect_equal(unname(d$tangent[1, , ]), rbind(c(1, 0), c(3, 1), c(0, 0)),
               tolerance = 1e-12)
  expect_equal(unname(d$hessian[1, , , ]), H, tolerance = 1e-12)
})


# -- Reverse mode ---------------------------------------------------------------

test_that("vjp contracts the Jacobian the forward path returns", {
  pars  <- c(a = 2, b = -0.5, k = 0.7)
  M     <- matrix(c(0.3, 1.1, 2.7), ncol = 1, dimnames = list(NULL, "t"))
  w     <- matrix(c(0.4, -1.3, 2.2, 0.9, -0.6, 1.7), nrow = 3, ncol = 2)

  fr <- fx$vjp_rev
  ff <- fx$vjp_fwd

  r <- fr$vjp(M, pars, w)
  J <- ff$jac(M, pars)

  expect_equal(unname(r$y), unname(ff$func(M, pars)), tolerance = 1e-12)

  # A variable is per observation, a parameter is shared, so the parameter
  # cotangent sums over observations and the variable one does not.
  wt <- vapply(seq_len(nrow(M)),
               function(o) sum(w[o, ] * J[o, , "t"]), numeric(1))
  expect_equal(unname(r$cotangentX[, 1, 1]), wt, tolerance = 1e-12)

  for (nm in names(pars))
    expect_equal(unname(r$cotangentP[nm, 1]), sum(w * J[, , nm]),
                 tolerance = 1e-12, label = paste("cotangentP", nm))
})

test_that("vjp sweeps several seeds against one recording", {
  pars  <- c(a = 0.6, b = 1.4)
  f     <- fx$vjp_seeds

  # Seeding the identity over the outputs recovers the full Jacobian row by row.
  w <- array(0, c(1, 2, 2))
  w[1, 1, 1] <- 1
  w[1, 2, 2] <- 1
  r <- f$vjp(NULL, pars, w)

  expect_equal(dim(r$cotangentP), c(2L, 2L))
  expect_equal(unname(r$cotangentP[, 1]), c(pars[["b"]], pars[["a"]]),
               tolerance = 1e-12)
  expect_equal(unname(r$cotangentP[, 2]), c(cos(pars[["a"]]), 2 * pars[["b"]]),
               tolerance = 1e-12)
})

test_that("derivMode builds exactly the directions it names", {
  eq <- c(y = "a * a")

  expect_error(cppFUN(eq, parameters = "a", derivMode = "symbolic",
                      modelname = "dm_symb"),
               'derivMode = "symbolic" is gone')
  # Second order is forward-only, so asking for it with the reverse direction
  # alone would otherwise return no Hessian without saying so.
  expect_error(cppFUN(eq, parameters = "a", deriv2 = TRUE, derivMode = "reverse",
                      modelname = "dm_bad2"),
               "no reverse counterpart")

  # Each direction is its own build product, and naming one omits the other.
  ff <- fx$dm_fwd
  expect_null(ff$vjp)
  expect_false(is.null(ff$jac))

  fr <- fx$dm_rev
  expect_null(fr$jac)
  expect_false(is.null(fr$vjp))
  expect_false(is.null(fr$func))

  fb <- fx$dm_both
  expect_false(is.null(fb$jac))
  expect_false(is.null(fb$vjp))

  # The two objects agree where they overlap.
  p <- c(a = 1.3)
  w <- matrix(1, 1, 1)
  expect_equal(unname(fr$vjp(NULL, p, w)$cotangentP[1, 1]),
               unname(fb$jac(NULL, p)[1, 1, "a"]), tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
#  vjp over a dual: forward over reverse on an observation function.
#
#  Oracle is the forward Hessian of the same object, contracted with the
#  cotangent and read along the tangents the inputs carry. Both are exact
#  derivatives of the same expressions, so the gap is rounding.
# ---------------------------------------------------------------------------

test_that("the dual vjp keeps the first order it already answered", {
  f <- fx$vjp_fr
  set.seed(3)
  X <- matrix(rnorm(8), 4L, 2L, dimnames = list(NULL, c("x1", "x2")))
  P <- c(a = 0.7, b = -0.4)
  W <- matrix(rnorm(8), 4L, 2L)

  r1 <- f$vjp(X, P, W)
  r2 <- f$vjp(X, P, W, tangentX = array(rnorm(24), c(4L, 2L, 3L)),
              tangentP = matrix(rnorm(6), 2L, 3L))
  # The dual instantiation sums the same expressions in its own order, so the
  # first order comes back to rounding rather than to the last bit.
  expect_equal(r2$cotangentX, r1$cotangentX, tolerance = 1e-12)
  expect_equal(r2$cotangentP, r1$cotangentP, tolerance = 1e-12)
  expect_identical(dim(r2$curvatureX), c(4L, 2L, 1L, 3L))
  expect_identical(dim(r2$curvatureP), c(2L, 1L, 3L))
})

test_that("the dual vjp answers the curvature the forward Hessian carries", {
  f <- fx$vjp_fr
  set.seed(3)
  n <- 4L; nd <- 3L
  X <- matrix(rnorm(n * 2), n, 2L, dimnames = list(NULL, c("x1", "x2")))
  P <- c(a = 0.7, b = -0.4)
  W <- matrix(rnorm(n * 2), n, 2L)
  VX <- array(rnorm(n * 2 * nd), c(n, 2L, nd))
  VP <- matrix(rnorm(2 * nd), 2L, nd)

  r <- f$vjp(X, P, W, tangentX = VX, tangentP = VP)
  H <- f$hess(x1 = X[, 1], x2 = X[, 2], a = P[["a"]], b = P[["b"]])
  nsym <- dim(H)[3L]

  V <- array(0, c(n, nsym, nd))
  V[, 1:2, ] <- VX
  for (k in seq_len(nd)) for (q in 1:2) V[, 2L + q, k] <- VP[q, k]

  ref_x <- array(0, c(n, 2L, 1L, nd))
  ref_p <- array(0, c(2L, 1L, nd))
  for (k in seq_len(nd)) for (o in seq_len(n)) {
    M <- W[o, 1L] * H[o, 1L, , ] + W[o, 2L] * H[o, 2L, , ]
    contrib <- as.numeric(M %*% V[o, , k])
    ref_x[o, , 1L, k] <- contrib[1:2]
    ref_p[, 1L, k] <- ref_p[, 1L, k] + contrib[3:4]
  }
  expect_equal(r$curvatureX, ref_x, tolerance = 1e-12)
  expect_equal(r$curvatureP, ref_p, tolerance = 1e-12)
})

test_that("a curvature alone goes through linearly", {
  # The vjp is linear in the cotangent, so a curvature alone has to reproduce a
  # first-order vjp taken with that direction as the cotangent.
  f <- fx$vjp_fr
  set.seed(5)
  n <- 4L; nd <- 2L
  X <- matrix(rnorm(n * 2), n, 2L, dimnames = list(NULL, c("x1", "x2")))
  P <- c(a = 0.7, b = -0.4)
  W <- matrix(rnorm(n * 2), n, 2L)
  DW <- array(rnorm(n * 2 * nd), c(n, 2L, 1L, nd))

  r <- f$vjp(X, P, W, curvature = DW)
  for (k in seq_len(nd)) {
    rk <- f$vjp(X, P, matrix(DW[, , 1L, k], n, 2L))
    expect_equal(unname(r$curvatureX[, , 1L, k]), unname(rk$cotangentX[, , 1L]),
                 tolerance = 1e-12, info = as.character(k))
    expect_equal(unname(r$curvatureP[, 1L, k]), unname(rk$cotangentP[, 1L]),
                 tolerance = 1e-12, info = as.character(k))
  }
})
