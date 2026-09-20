\donttest{
if (codegenAvailable()) {
  f <- cppODE(c(x = "-k*x"), deriv = FALSE)
  batchAvailable(f)
}
}
