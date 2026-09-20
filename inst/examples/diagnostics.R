\donttest{
if (codegenAvailable()) {
  f <- cppODE(c(x = "-k*x"), deriv = FALSE)
  diagnostics(solveODE(f, times = 0:5, parms = c(x = 1, k = 0.3)))
}
}
