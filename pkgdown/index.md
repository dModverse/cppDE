# cppDE

cppDE writes, compiles and runs a C++ solver for a system of ordinary
differential equations given as equations in R, with the derivatives of the
solution in the parameters, forward or reverse, to first and second order.

## Features

- Integrators: variable-order BDF/NDF and Adams in Nordsieck form, a
  Rosenbrock method of order 4 and the Tsitouras 5(4) Runge-Kutta pair, with
  dense or sparse (KLU) linear algebra
- Forward sensitivities on dual numbers, gradients by a discrete adjoint at a
  cost independent of the parameter count, and exact Hessians forward over
  forward or forward over reverse
- Time- and root-triggered events, switches and sliding, forcings, with
  sensitivities across every jump
- Jacobians and derivative contractions by automatic differentiation on an
  expression graph, emitted as straight-line C++
- The same models on SUNDIALS CVODES with `cvode()`, algebraic functions with
  `cppFUN()`, and batches of conditions solved in parallel

## Example

A decay chain, compiled once and solved with the sensitivities of every state
in the initial values and parameters:

``` r
library(cppDE)

model <- cppODE(c(A = "-k1*A", B = "k1*A - k2*B"))
out <- solveODE(model, times = 0:10, parms = c(A = 1, B = 0, k1 = 0.5, k2 = 0.2))

out$variable[11, ]       # states at t = 10
#>           A           B 
#> 0.006738166 0.214328533
out$tangent[11, , ]      # [state, parameter] at t = 10
#>         sens
#> variable           A         B          k1        k2
#>        A 0.006738166 0.0000000 -0.06737862  0.000000
#>        B 0.214328533 0.1353353 -0.17347368 -1.541159
```

## Where to go next

- [Get started](articles/cppDE.html): a decay chain from equations to a fit
- [Guide](articles/guide-model.html): every entry point, chapter by chapter
- [Installation and optional libraries](articles/installation.html): the
  requirements, SUNDIALS, KLU and OpenMP, and how to diagnose a failure
- [Methods](articles/methods-overview.html): what the generated code computes
- Examples: [sensitivities across events, checked
  symbolically](articles/saltation.html), and [fitting magnetization
  dynamics](articles/fitCoPt.html)
- [Benchmarks](articles/benchmarks.html): the native integrators against
  SUNDIALS CVODES
- [dMod2](https://dmodverse.github.io/dMod2/) builds a modeling framework on
  cppDE: observation functions, parameter transformations, objectives,
  optimizers and profile likelihoods.
