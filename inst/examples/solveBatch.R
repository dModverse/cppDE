\donttest{
if (codegenAvailable()) {
  ## Marshal the conditions once, then re-solve with new numbers, as inside an
  ## optimiser.
  f <- cppODE(c(x = "-k*x"), deriv = FALSE, compile = FALSE)
  r <- cppODE(c(x = "-k*x"), derivMode = "reverse", compile = FALSE)
  compile(f, r, output = "solveBatch_example")
  conditions <- list(a = list(parms = c(x = 1, k = 0.1)),
                     b = list(parms = c(x = 2, k = 0.1)))
  h   <- prepareBatch(f, conditions, times = 0:5)
  res <- solveBatch(h, parms = list(c(x = 1, k = 0.5), c(x = 2, k = 0.5)))
  res$b$variable

  ## Reverse mode: the forward passes keep their checkpoints, and a prepared
  ## sweep reuses them for a new cotangent.
  fw <- solveODEBatch(r, lapply(conditions, c, keepStore = TRUE), times = 0:5)
  W  <- matrix(1, 6, 1, dimnames = list(NULL, "x"))
  hr <- prepareBatch(r, conditions, times = 0:5, cotangent = W)
  solveBatch(hr, cotangent = list(W, 2 * W),
             store = lapply(fw, `[[`, "store"))$b$cotangent
}
}
