# Test cppFUN() algebraic function compilation with sensitivities.

skip_on_cran()

# -- Basic cppFUN output structure ---------------------------------------------

test_that("cppFUN returns correct output", {
  trafo <- c(
    y1 = "a * exp(-b * x)",
    y2 = "a + b * x"
  )

  f <- cppFUN(trafo, parameters = c("a", "b", "x"),
              deriv = TRUE, modelname = "fun_basic", convenient = TRUE)

  res <- f$func(a = 2, b = 0.5, x = 1)

  expect_true(is.matrix(res))
  expect_equal(colnames(res), c("y1", "y2"))
  expect_equal(unname(res[1, "y1"]), 2 * exp(-0.5), tolerance = 1e-10)
  expect_equal(unname(res[1, "y2"]), 2.5, tolerance = 1e-10)
})

test_that("parameters named like the generated arrays do not collide", {
  # The emitted code indexes p[] and x_obs[]. A parameter of that name is
  # substituted for its own slot, in whatever order the parameters are listed.
  trafo <- c(o1 = "p + 2 * x", o2 = "y * k",
             o3 = "x_obs + y_local", o4 = "p * y_local - k")
  pars   <- c(p = 2, x = 3, y = 5, k = 7, x_obs = 11, y_local = 13)
  orders <- list(names(pars), rev(names(pars)),
                 c("k", "p", "y_local", "x", "y", "x_obs"))

  for (i in seq_along(orders)) {
    f <- cppFUN(trafo, parameters = orders[[i]], deriv = TRUE,
                derivMode = "symbolic", modelname = paste0("collide_", i),
                convenient = TRUE)
    res <- do.call(f$func, as.list(pars[orders[[i]]]))
    jac <- do.call(f$jac, as.list(pars[orders[[i]]]))[1, , ]

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
  trafo <- c(o1 = "std * exp(ini * log(10)) + default",
             o2 = "sqrt(int) + std^2")

  for (mode in c("forward", "symbolic")) {
    f <- cppFUN(trafo, variables = "int", parameters = c("std", "ini", "default"),
                deriv = TRUE, derivMode = mode, compile = TRUE,
                modelname = paste0("cxx_tokens_", mode), convenient = TRUE)
    res <- f$func(int = 4, std = 2, ini = 0.5, default = 3)
    jac <- f$jac(int = 4, std = 2, ini = 0.5, default = 3)[1, , ]

    expect_equal(unname(res[1, ]), c(2 * 10^0.5 + 3, 6), tolerance = 1e-10,
                 label = paste("values,", mode))
    expect_equal(unname(jac["o1", "std"]), 10^0.5, tolerance = 1e-8)
    expect_equal(unname(jac["o1", "default"]), 1, tolerance = 1e-10)
    expect_equal(unname(jac["o2", "int"]), 0.25, tolerance = 1e-8)
  }
})

test_that("a Python keyword as a symbol name is rejected", {
  # SymPy parses through Python, where the name is a syntax error, and True,
  # False and None read as constants and drop the symbol from the model.
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
  trafo <- c(
    y1 = "a * exp(-b * x)",
    y2 = "a + b * x"
  )

  f <- cppFUN(trafo, parameters = c("a", "b", "x"),
              deriv = TRUE, derivMode = "symbolic",
              modelname = "fun_jac", convenient = TRUE)

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
  trafo <- c(y = "a * b * x^2")

  f <- cppFUN(trafo, parameters = c("a", "b", "x"),
              deriv = TRUE, deriv2 = TRUE, derivMode = "symbolic",
              modelname = "fun_hess", convenient = TRUE)

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
  trafo <- c(y = "a * b + c")

  f <- cppFUN(trafo, parameters = c("a", "b", "c"),
              fixed = "c", deriv = TRUE, derivMode = "symbolic",
              modelname = "fun_fixed", convenient = TRUE)

  jac <- f$jac(a = 2, b = 3, c = 1)

  # Only 2 params in Jacobian (a, b), not c
  expect_equal(dim(jac)[3], 2)
  jac_names <- dimnames(jac)[[3]]
  expect_false("c" %in% jac_names)
})

# -- forward vs symbolic agree on raw derivatives --------------------------------

test_that("cppFUN forward and symbolic give identical raw jac/hess", {
  trafo <- c(y1 = "a*x^2 + b", y2 = "sin(c*x)")
  pars  <- list(a = 2, b = 1, c = 0.3, x = 3)
  out_d <- list(); out_s <- list()
  for (mode in c("forward", "symbolic")) {
    f <- cppFUN(trafo, parameters = c("a", "b", "c", "x"),
                deriv = TRUE, deriv2 = TRUE, derivMode = mode,
                compile = TRUE, modelname = paste0("xs_raw_", mode))
    res <- list()
    res$y    <- do.call(f$func, pars)
    res$jac  <- do.call(f$jac,  pars)
    res$hess <- do.call(f$hess, pars)
    if (mode == "forward") out_d <- res else out_s <- res
  }
  expect_equal(unname(out_d$y),    unname(out_s$y),    tolerance = 1e-12)
  expect_equal(unname(out_d$jac),  unname(out_s$jac),  tolerance = 1e-10)
  expect_equal(unname(out_d$hess), unname(out_s$hess), tolerance = 1e-10)
})

# -- forward vs symbolic agree under chain rule, including dX2/dP2 ---------------

test_that("cppFUN forward and symbolic agree under second-order chain rule", {
  trafo <- c(y1 = "a*x^2 + b", y2 = "sin(c*x)")
  th    <- c("th1", "th2", "th3")
  # Linear part of Phi: theta -> (a, b, c) = (2*th1, 1*th2, 1*th3 + th1)
  dP <- matrix(0, 4, length(th),
               dimnames = list(c("a", "b", "c", "x"), th))
  dP["a", "th1"] <- 2
  dP["b", "th2"] <- 1
  dP["c", "th3"] <- 1
  dP["c", "th1"] <- 1
  # Nonlinear quadratic part: a depends on th1*th2 with coefficient 0.5.
  dP2 <- array(0, c(4, length(th), length(th)),
               dimnames = list(c("a", "b", "c", "x"), th, th))
  dP2["a", "th1", "th2"] <- 0.5
  dP2["a", "th2", "th1"] <- 0.5
  pars <- list(a = 2, b = 1, c = 0.3, x = 3)
  out  <- list()
  for (mode in c("forward", "symbolic")) {
    f <- cppFUN(trafo, parameters = c("a", "b", "c", "x"),
                deriv = TRUE, deriv2 = TRUE, derivMode = mode,
                compile = TRUE, modelname = paste0("xs_chain_", mode))
    out[[mode]] <- do.call(f$evaluate, c(pars, list(dP = dP, dP2 = dP2,
                                                    deriv2 = TRUE)))
  }
  expect_equal(out$forward$y,   out$symbolic$y,   tolerance = 1e-12)
  expect_equal(out$forward$dy,  out$symbolic$dy,  tolerance = 1e-10)
  expect_equal(out$forward$d2y, out$symbolic$d2y, tolerance = 1e-9)
})

# -- forward + deriv2 + identity pass-through (regression) ------------------------
# An output that is just a bare param or var arms no inner Hessian seed, so the
# dual2nd writeback has to reach it through the bounds-safe const overload.

test_that("cppFUN forward deriv2 handles identity pass-through", {
  trafo <- c(la = "la", y2 = "la^2 + b", zero = "0")
  pars  <- list(la = 1.5, b = 0.7)
  for (mode in c("forward", "symbolic")) {
    f <- cppFUN(trafo, parameters = c("la", "b"),
                deriv = TRUE, deriv2 = TRUE, derivMode = mode,
                compile = TRUE, modelname = paste0("xs_passthru_", mode))
    out <- do.call(f$evaluate, c(pars, list(deriv2 = TRUE)))
    if (mode == "forward") d <- out else s <- out
  }
  expect_equal(unname(d$y),   unname(s$y),   tolerance = 1e-12)
  expect_equal(unname(d$dy),  unname(s$dy),  tolerance = 1e-10)
  expect_equal(unname(d$d2y), unname(s$d2y), tolerance = 1e-10)
})


# -- Reverse mode ---------------------------------------------------------------

test_that("vjp contracts the Jacobian the symbolic path returns", {
  trafo <- c(y1 = "a * exp(-k * t) + b", y2 = "log(a + k * k) * t")
  pars  <- c(a = 2, b = -0.5, k = 0.7)
  M     <- matrix(c(0.3, 1.1, 2.7), ncol = 1, dimnames = list(NULL, "t"))
  w     <- matrix(c(0.4, -1.3, 2.2, 0.9, -0.6, 1.7), nrow = 3, ncol = 2)

  fr <- cppFUN(trafo, variables = "t", parameters = names(pars), deriv = TRUE,
               derivMode = "reverse", modelname = "vjp_rev", compile = TRUE,
               convenient = FALSE)
  fs <- cppFUN(trafo, variables = "t", parameters = names(pars), deriv = TRUE,
               derivMode = "symbolic", modelname = "vjp_symb", compile = TRUE,
               convenient = FALSE)

  r <- fr$vjp(M, pars, w)
  J <- fs$jac(M, pars)

  expect_equal(unname(r$y), unname(fs$func(M, pars)), tolerance = 1e-12)

  # A variable is per observation, a parameter is shared, so the parameter
  # cotangent sums over observations and the variable one does not.
  wt <- vapply(seq_len(nrow(M)),
               function(o) sum(w[o, ] * J[o, , "t"]), numeric(1))
  expect_equal(unname(r$wx[, 1, 1]), wt, tolerance = 1e-12)

  for (nm in names(pars))
    expect_equal(unname(r$wp[nm, 1]), sum(w * J[, , nm]), tolerance = 1e-12,
                 label = paste("wp", nm))
})

test_that("vjp sweeps several seeds against one recording", {
  trafo <- c(y1 = "a * b", y2 = "sin(a) + b * b")
  pars  <- c(a = 0.6, b = 1.4)

  f <- cppFUN(trafo, variables = NULL, parameters = names(pars), deriv = TRUE,
              derivMode = "reverse", modelname = "vjp_seeds", compile = TRUE,
              convenient = FALSE)

  # Seeding the identity over the outputs recovers the full Jacobian row by row.
  w <- array(0, c(1, 2, 2))
  w[1, 1, 1] <- 1
  w[1, 2, 2] <- 1
  r <- f$vjp(NULL, pars, w)

  expect_equal(dim(r$wp), c(2L, 2L))
  expect_equal(unname(r$wp[, 1]), c(pars[["b"]], pars[["a"]]), tolerance = 1e-12)
  expect_equal(unname(r$wp[, 2]), c(cos(pars[["a"]]), 2 * pars[["b"]]),
               tolerance = 1e-12)
})

test_that("derivMode builds exactly the directions it names", {
  eq <- c(y = "a * a")

  # symbolic is a backend for the forward Jacobian, not a direction.
  fs <- cppFUN(eq, parameters = "a", deriv = TRUE, derivMode = "symbolic",
               modelname = "dm_symb", convenient = FALSE)
  expect_null(fs$vjp)
  expect_false(is.null(fs$jac))
  expect_error(cppFUN(eq, parameters = "a", derivMode = c("symbolic", "reverse"),
                      modelname = "dm_bad"),
               "cannot be combined")
  # Second order is forward-only, so asking for it with the reverse direction
  # alone would otherwise return no Hessian without saying so.
  expect_error(cppFUN(eq, parameters = "a", deriv2 = TRUE, derivMode = "reverse",
                      modelname = "dm_bad2"),
               "no reverse counterpart")

  # Each direction is its own build product, and naming one omits the other.
  ff <- cppFUN(eq, parameters = "a", deriv = TRUE, derivMode = "forward",
               modelname = "dm_fwd", compile = TRUE, convenient = FALSE)
  expect_null(ff$vjp)
  expect_false(is.null(ff$jac))

  fr <- cppFUN(eq, parameters = "a", deriv = TRUE, derivMode = "reverse",
               modelname = "dm_rev", compile = TRUE, convenient = FALSE)
  expect_null(fr$jac)
  expect_false(is.null(fr$vjp))
  expect_false(is.null(fr$func))

  fb <- cppFUN(eq, parameters = "a", deriv = TRUE,
               derivMode = c("forward", "reverse"),
               modelname = "dm_both", compile = TRUE, convenient = FALSE)
  expect_false(is.null(fb$jac))
  expect_false(is.null(fb$vjp))

  # The two objects agree where they overlap.
  p <- c(a = 1.3)
  w <- matrix(1, 1, 1)
  expect_equal(unname(fr$vjp(NULL, p, w)$wp[1, 1]),
               unname(fb$jac(NULL, p)[1, 1, "a"]), tolerance = 1e-12)
})
