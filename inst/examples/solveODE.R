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
out$tangent[11, , ]                   # [state, parameter] at t = 10

## Only two directions: the tangent seeds the rate constants, one column each
S <- matrix(0, 2, 2, dimnames = list(c("k1", "k2"), c("k1", "k2")))
diag(S) <- 1
solveODE(fwd, times, parms, tangent = S)$tangent[11, , ]

## Reverse mode: a cotangent weights the states at every time, and $cotangent is
## the gradient of sum(cotangent * x) in the initial values and parameters.
## keepStore keeps the forward pass, so each further cotangent costs one
## backward sweep.
first <- solveODE(rev, times, parms, keepStore = TRUE)
onB   <- cbind(A = 0, B = rep(1, length(times)))
onA   <- cbind(A = rep(1, length(times)), B = 0)
solveODE(rev, times, parms, cotangent = onB, store = first$store)$cotangent
solveODE(rev, times, parms, cotangent = onA, store = first$store)$cotangent
}
