\dontrun{
## Needs SUNDIALS at install time, see install_libs(). The models are
## generated first and compiled together, into one shared object in tempdir().
eqns <- c(A = "-k1*A", B = "k1*A - k2*B")

plain <- cvode(eqns, compile = FALSE)
sens  <- cvode(eqns, deriv = TRUE, compile = FALSE)
adj   <- cvode(eqns, derivMode = "reverse", compile = FALSE)
compile(plain, sens, adj, output = "cvode_examples")

times <- 0:10
parms <- c(A = 1, B = 0, k1 = 0.5, k2 = 0.2)
solveODE(plain, times, parms)$variable[11, ]
solveODE(sens, times, parms)$tangent[11, "B", ] # forward sensitivities

## CVODES adjoint: the gradient of sum(B) over the grid
onB <- cbind(A = 0, B = rep(1, length(times)))
solveODE(adj, times, parms, cotangent = onB)$cotangent
}
