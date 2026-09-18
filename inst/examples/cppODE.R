\donttest{
## Six variants of one decay chain A -> B. Each is generated first and all
## are compiled together afterwards, into one shared object in tempdir().
eqns <- c(A = "-k1*A", B = "k1*A - k2*B")

## First-order sensitivities, the default
fwd <- cppODE(eqns, compile = FALSE)

## A dose of A at the parameter time t_e
events <- data.frame(var = "A", time = "t_e", value = 1, method = "add")
dosed  <- cppODE(eqns, events = events, compile = FALSE)

## A constant inflow k0, integrated until the steady state
steady <- cppODE(c(A = "k0 - k1*A", B = "k1*A - k2*B"),
                 rootfunc = "equilibrate", compile = FALSE)

## A time-dependent input u, given as data when solving
forced <- cppODE(c(A = "u - k1*A", B = "k1*A - k2*B"), forcings = "u",
                 compile = FALSE)

## Second-order sensitivities
fwd2 <- cppODE(eqns, deriv2 = TRUE, compile = FALSE)

## Reverse mode: the states in plain double, the derivatives from one
## backward sweep whose cost does not grow with the number of parameters
rev <- cppODE(eqns, derivMode = "reverse", compile = FALSE)

compile(fwd, dosed, steady, forced, fwd2, rev, output = "cppODE_examples")

times <- seq(0, 10, by = 0.5)
parms <- c(A = 1, B = 0, k1 = 0.5, k2 = 0.2)

out <- solveODE(fwd, times, parms)
head(out$variable)
out$sens1[21, "B", ]                  # dB(10) / d(A, B, k1, k2)

out <- solveODE(dosed, times, c(parms, t_e = 4))
out$variable[times %in% c(3.5, 4, 4.5), ]

out <- solveODE(steady, times = c(0, 1000), c(parms, k0 = 1))
tail(out$time, 1)                     # stopped before t = 1000
tail(out$variable, 1)                 # the steady state

u   <- list(u = data.frame(time = times, value = 1 + sin(times)))
out <- solveODE(forced, times, parms, forcings = u)
out$variable[21, ]

out <- solveODE(fwd2, times, parms)
out$sens2[21, "B", , ]                # Hessian of B(10)

## The gradient of sum(B) over the grid: a seed of ones on B
seed <- cbind(A = 0, B = rep(1, length(times)))
solveODE(rev, times, parms, seed = seed)$adjoint
}
