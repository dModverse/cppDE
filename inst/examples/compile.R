\donttest{
## Generate an ODE and an observable without compiling, then link both into
## one shared object next to their sources in tempdir().
f <- cppODE(c(x = "-k*x"), deriv = FALSE, compile = FALSE)
g <- cppFUN(c(y = "2*x"), variables = "x", compile = FALSE)
compile(f, g, output = "decay_and_double")
solveODE(f, times = 0:3, parms = c(x = 1, k = 1))$variable
g$func(x = 3)
}
