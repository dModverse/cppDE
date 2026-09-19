\donttest{
if (codegenAvailable()) {
  ## Two observables of two states and a scale, compiled with their Jacobian.
  ## The functions take the inputs as named arguments.
  g <- cppFUN(c(y = "s*(x1 + x2)", z = "log(x1)"), variables = c("x1", "x2"),
              parameters = "s", compile = TRUE)
  g$func(x1 = 1, x2 = 2, s = 0.5)
  g$jac(x1 = 1, x2 = 2, s = 0.5)[1, , ]
}
}
