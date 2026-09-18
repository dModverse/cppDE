\donttest{
f <- cppODE(c(x = "-k*x"), deriv = FALSE)
batchAvailable(f)
}
