\donttest{
if (codegenAvailable()) {
  ## One model, several conditions, solved on threads in one call. A condition's
  ## list overrides the batch-wide arguments.
  f <- cppODE(c(x = "-k*x"), deriv = FALSE)
  conditions <- list(slow = list(parms = c(x = 1, k = 0.1)),
                     fast = list(parms = c(x = 1, k = 1)))
  res <- solveODEBatch(f, conditions, times = 0:5)
  sapply(res, function(r) r$variable[, "x"])
}
}
