# cppDE

cppDE writes, compiles and runs a C++ solver for a system of ordinary
differential equations given as equations in R, with the derivatives of the
solution in the initial values and parameters, forward or reverse, to first
and second order.

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
  cppDE: observation functions, parameter transformations, objective
  functions, optimizers and profile likelihoods.
