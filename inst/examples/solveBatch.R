\donttest{
if (codegenAvailable()) {
  ## Marshal the conditions once, then re-solve with new numbers, as inside an
  ## optimiser.
  f <- cppODE(c(x = "-k*x"), deriv = FALSE)
  conditions <- list(a = list(parms = c(x = 1, k = 0.1)),
                     b = list(parms = c(x = 2, k = 0.1)))
  h   <- prepareBatch(f, conditions, times = 0:5)
  res <- solveBatch(h, parms = list(c(x = 1, k = 0.5), c(x = 2, k = 0.5)))
  res$b$variable
}
}
