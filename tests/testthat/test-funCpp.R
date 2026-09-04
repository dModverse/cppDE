# Test funCpp() algebraic function compilation with sensitivities.

skip_on_cran()

# -- Basic funCpp output structure ---------------------------------------------

test_that("funCpp returns correct output", {
  trafo <- c(
    y1 = "a * exp(-b * x)",
    y2 = "a + b * x"
  )

  f <- funCpp(trafo, parameters = c("a", "b", "x"),
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
    f <- funCpp(trafo, parameters = orders[[i]], deriv = TRUE,
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

  for (mode in c("dual", "symbolic")) {
    f <- funCpp(trafo, variables = "int", parameters = c("std", "ini", "default"),
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
  expect_error(funCpp(c(y = "class * 2"), modelname = "py_kw_eqn"),
               "Python keyword used as a symbol name: 'class'")
  expect_error(funCpp(c(y = "a * 2"), parameters = c("a", "lambda"),
                      modelname = "py_kw_par"),
               "'lambda'")
  expect_error(funCpp(c(lambda = "a * 2"), parameters = "a",
                      modelname = "py_kw_out"),
               "'lambda'")
  expect_error(funCpp(c(y = "a * True"), parameters = c("a", "True"),
                      modelname = "py_kw_const"),
               "'True'")
})

# -- Jacobian correctness -----------------------------------------------------

test_that("funCpp Jacobian matches analytical derivatives", {
  trafo <- c(
    y1 = "a * exp(-b * x)",
    y2 = "a + b * x"
  )

  f <- funCpp(trafo, parameters = c("a", "b", "x"),
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

test_that("funCpp Hessian has correct dimensions and is symmetric", {
  trafo <- c(y = "a * b * x^2")

  f <- funCpp(trafo, parameters = c("a", "b", "x"),
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

# -- Fixed parameters in funCpp ------------------------------------------------

test_that("funCpp fixed parameters are excluded from derivatives", {
  trafo <- c(y = "a * b + c")

  f <- funCpp(trafo, parameters = c("a", "b", "c"),
              fixed = "c", deriv = TRUE, derivMode = "symbolic",
              modelname = "fun_fixed", convenient = TRUE)

  jac <- f$jac(a = 2, b = 3, c = 1)

  # Only 2 params in Jacobian (a, b), not c
  expect_equal(dim(jac)[3], 2)
  jac_names <- dimnames(jac)[[3]]
  expect_false("c" %in% jac_names)
})

# -- dual vs symbolic agree on raw derivatives --------------------------------

test_that("funCpp dual and symbolic give identical raw jac/hess", {
  trafo <- c(y1 = "a*x^2 + b", y2 = "sin(c*x)")
  pars  <- list(a = 2, b = 1, c = 0.3, x = 3)
  out_d <- list(); out_s <- list()
  for (mode in c("dual", "symbolic")) {
    f <- funCpp(trafo, parameters = c("a", "b", "c", "x"),
                deriv = TRUE, deriv2 = TRUE, derivMode = mode,
                compile = TRUE, modelname = paste0("xs_raw_", mode))
    res <- list()
    res$y    <- do.call(f$func, pars)
    res$jac  <- do.call(f$jac,  pars)
    res$hess <- do.call(f$hess, pars)
    if (mode == "dual") out_d <- res else out_s <- res
  }
  expect_equal(unname(out_d$y),    unname(out_s$y),    tolerance = 1e-12)
  expect_equal(unname(out_d$jac),  unname(out_s$jac),  tolerance = 1e-10)
  expect_equal(unname(out_d$hess), unname(out_s$hess), tolerance = 1e-10)
})

# -- dual vs symbolic agree under chain rule, including dX2/dP2 ---------------

test_that("funCpp dual and symbolic agree under second-order chain rule", {
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
  for (mode in c("dual", "symbolic")) {
    f <- funCpp(trafo, parameters = c("a", "b", "c", "x"),
                deriv = TRUE, deriv2 = TRUE, derivMode = mode,
                compile = TRUE, modelname = paste0("xs_chain_", mode))
    out[[mode]] <- do.call(f$evaluate, c(pars, list(dP = dP, dP2 = dP2,
                                                    deriv2 = TRUE)))
  }
  expect_equal(out$dual$y,   out$symbolic$y,   tolerance = 1e-12)
  expect_equal(out$dual$dy,  out$symbolic$dy,  tolerance = 1e-10)
  expect_equal(out$dual$d2y, out$symbolic$d2y, tolerance = 1e-9)
})

# -- dual + deriv2 + identity pass-through (regression) ------------------------
# An output that is just a bare param or var arms no inner Hessian seed, so the
# dual2nd writeback has to reach it through the bounds-safe const overload.

test_that("funCpp dual deriv2 handles identity pass-through", {
  trafo <- c(la = "la", y2 = "la^2 + b", zero = "0")
  pars  <- list(la = 1.5, b = 0.7)
  for (mode in c("dual", "symbolic")) {
    f <- funCpp(trafo, parameters = c("la", "b"),
                deriv = TRUE, deriv2 = TRUE, derivMode = mode,
                compile = TRUE, modelname = paste0("xs_passthru_", mode))
    out <- do.call(f$evaluate, c(pars, list(deriv2 = TRUE)))
    if (mode == "dual") d <- out else s <- out
  }
  expect_equal(unname(d$y),   unname(s$y),   tolerance = 1e-12)
  expect_equal(unname(d$dy),  unname(s$dy),  tolerance = 1e-10)
  expect_equal(unname(d$d2y), unname(s$d2y), tolerance = 1e-10)
})
