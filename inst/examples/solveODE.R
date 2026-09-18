\donttest{
## A forward and a reverse build of one model, compiled together into one
## shared object in tempdir().
eqns <- c(A = "-k1*A", B = "k1*A - k2*B")
fwd  <- cppODE(eqns, compile = FALSE)
rev  <- cppODE(eqns, derivMode = "reverse", compile = FALSE)
compile(fwd, rev, output = "solveODE_examples")

times <- 0:10
parms <- c(A = 1, B = 0, k1 = 0.5, k2 = 0.2)

## Values and sensitivities at tighter tolerances than the default 1e-6
out <- solveODE(fwd, times, parms, abstol = 1e-10, reltol = 1e-10)
out$variable[11, ]
out$sens1[11, , ]                     # [state, parameter] at t = 10

## Only two directions: sens1ini seeds the rate constants, one column each
S <- matrix(0, 2, 2, dimnames = list(c("k1", "k2"), c("k1", "k2")))
diag(S) <- 1
solveODE(fwd, times, parms, sens1ini = S)$sens1[11, , ]

## Reverse mode: a seed weights the states at every time, and $adjoint is the
## gradient of sum(seed * x) in the initial values and parameters. keepStore
## keeps the forward pass, so each further seed costs one backward sweep.
first <- solveODE(rev, times, parms, keepStore = TRUE)
onB   <- cbind(A = 0, B = rep(1, length(times)))
onA   <- cbind(A = rep(1, length(times)), B = 0)
solveODE(rev, times, parms, seed = onB, store = first$store)$adjoint
solveODE(rev, times, parms, seed = onA, store = first$store)$adjoint
}
