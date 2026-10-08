# The cppde::dual2nd math primitives against stats::D(), an independent
# derivation in R. One output per primitive, all in one model.

skip_on_cran()

prims <- c(add = "a + b", mul = "a * b", div = "a / b", sub = "a - b",
           pow = "a^b", sq = "a^2", comp = "a*sin(b) + exp(a)",
           sym = "a*b + b*c + a*c",
           sin = "sin(x)", cos = "cos(x)", tan = "tan(x)", exp = "exp(x)",
           log = "log(x)", sqrt = "sqrt(x)", sinh = "sinh(x)", cosh = "cosh(x)",
           tanh = "tanh(x)")
pars <- c("a", "b", "c", "x")
d2prim <- cppFUN(prims, parameters = pars, deriv = TRUE, deriv2 = TRUE,
                 derivMode = "forward", modelname = "d2prim")
compile(d2prim, output = "test_dual2nd_primitives", cores = test_cores())

# Points away from every singularity of the primitives above, one per row.
pts <- rbind(c(a = 1.5, b = 2.5, c = 3, x = 0.7),
             c(a = 1.7, b = 2.3, c = 1, x = 1.7),
             c(a = 0.4, b = 1.1, c = 2, x = 0.5))

test_that("value, gradient and Hessian of every primitive match D()", {
  # Identity seeds: each parameter is its own direction, no curvature in Phi.
  dP  <- diag(4); dimnames(dP) <- list(pars, pars)
  dP2 <- array(0, c(4, 4, 4), dimnames = list(pars, pars, pars))
  for (i in seq_len(nrow(pts))) {
    p <- pts[i, ]
    out <- do.call(d2prim$evaluate,
                   c(as.list(p), list(tangentP = dP, hessianP = dP2, deriv2 = TRUE)))
    for (nm in names(prims)) {
      e <- str2lang(prims[[nm]])
      at <- function(z) eval(z, as.list(p))
      dy <- vapply(pars, function(v) at(D(e, v)), 0)
      d2y <- outer(pars, pars, Vectorize(function(u, v) at(D(D(e, u), v))))
      info <- paste(nm, "at point", i)
      expect_equal(unname(out$y[1, nm]), at(e), tolerance = 1e-12, info = info)
      expect_equal(unname(out$tangent[1, nm, ]), unname(dy), tolerance = 1e-10, info = info)
      expect_equal(unname(out$hessian[1, nm, , ]), unname(d2y), tolerance = 1e-10,
                   info = info)
    }
  }
})

test_that("the Hessian comes out symmetric entry by entry", {
  out <- do.call(d2prim$hess, as.list(pts[1, ]))
  for (nm in names(prims))
    expect_identical(out[1, nm, , ], t(out[1, nm, , ]), info = nm)
})
