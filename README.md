
<!-- README.md is generated from README.Rmd. Please edit that file -->

# cppDE

<!-- badges: start -->

[![R-CMD-check](https://github.com/dModverse/cppDE/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/dModverse/cppDE/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

**C++ code generation for ODEs with sensitivity propagation.**

cppDE writes, compiles and runs a C++ solver for a system of ordinary
differential equations, or C++ code for algebraic functions, given as
equations in R. Derivatives with respect to the parameters are
propagated forward on dual numbers or backward by a discrete adjoint, to
first and second order.

[dMod2](https://dmodverse.github.io/dMod2/) builds a modeling framework
on cppDE.

## Features

- **Solvers**: variable-order BDF/NDF and Adams multistep in Nordsieck
  form, adapted from
  [SUNDIALS/CVODE](https://computing.llnl.gov/projects/sundials);
  Rosenbrock4 extended from
  [Boost.Odeint](https://www.boost.org/doc/libs/release/libs/numeric/odeint/);
  the Tsitouras 5(4) pair, an explicit Runge-Kutta method with dense
  output
- **Sensitivities in both directions**: forward mode integrates one
  tangent per parameter through the solve; reverse mode returns the
  gradient of a functional in one backward sweep, at a cost independent
  of the parameter count
- **Second order**: exact Hessians, forward over forward or forward over
  reverse
- **Events and switches**: time- and root-triggered events, and
  `piecewise()` or comparisons on states and time located as roots; a
  switch whose branches both push onto the switching surface slides
  along it (Filippov). First- and second-order sensitivities take every
  jump by the saltation correction
- **Code generation on an expression graph**: Jacobians and derivative
  contractions by automatic differentiation, emitted as straight-line
  C++, with automatic sparsity detection and KLU factorization for large
  systems
- **`cvode()`**: the same models on SUNDIALS CVODES, with forward
  sensitivities and adjoint
- **`cppFUN()`**: the same code generation for algebraic functions, with
  Jacobian, Hessian and vector-Jacobian products
- **Batched solving**: `solveODEBatch()` integrates many conditions in
  one `.Call` over OpenMP threads, without the fork and serialization of
  `parallel::mclapply()`

## Installation

``` r
# install.packages("devtools")
devtools::install_github("dModverse/cppDE")
```

Model code is generated and compiled at run time, so a **C++17
compiler** and **Python with SymPy** are required. Windows users need
[Rtools](https://cran.r-project.org/bin/windows/Rtools/). SUNDIALS (for
`cvode()`), SuiteSparse/KLU (for sparse Jacobians) and OpenMP (for
parallel batches) are optional; [Installation and optional
libraries](https://dmodverse.github.io/cppDE/articles/installation.html)
lists how to install them and how to diagnose a failing build or solve.

## Example

A two-state model, compiled once and solved with the sensitivities of
every state in the initial values and parameters:

``` r
library(cppDE)

model <- cppODE(c(A = "-k1*A", B = "k1*A - k2*B"))
out <- solveODE(model, times = 0:10, parms = c(A = 1, B = 0, k1 = 0.5, k2 = 0.2))

out$variable[11, ]       # states at t = 10
#>           A           B 
#> 0.006738521 0.214327979
out$tangent[11, , ]      # [state, parameter] at t = 10
#>         sens
#> variable           A         B          k1        k2
#>        A 0.006738521 0.0000000 -0.06737644  0.000000
#>        B 0.214327979 0.1353353 -0.17347569 -1.541163
```

`derivMode = "reverse"` builds the same model for gradients by the
adjoint, `cvode()` for SUNDIALS CVODES; `?expressions` lists what an
equation can contain.

## Documentation

The site, <https://dmodverse.github.io/cppDE/>, has the reference and
these articles:

- [Get started](https://dmodverse.github.io/cppDE/articles/cppDE.html)
- [Guide](https://dmodverse.github.io/cppDE/articles/guide-model.html):
  every entry point, chapter by chapter
- [Installation and optional
  libraries](https://dmodverse.github.io/cppDE/articles/installation.html)
- [Methods](https://dmodverse.github.io/cppDE/articles/methods-overview.html):
  what the generated code computes
- Examples: [Sensitivities across events, checked
  symbolically](https://dmodverse.github.io/cppDE/articles/saltation.html)
  and [Fitting magnetization
  dynamics](https://dmodverse.github.io/cppDE/articles/fitCoPt.html)
- [Benchmarks](https://dmodverse.github.io/cppDE/articles/benchmarks.html)
  against SUNDIALS CVODES

## License

cppDE is distributed under the MIT License (`LICENSE`, `LICENSE.md`).
Portions of the C++ solver core in `inst/include/cppde/` are derived
from third-party projects and remain under their upstream licenses: the
Nordsieck multistep stepper from
[SUNDIALS/CVODE(S)](https://computing.llnl.gov/projects/sundials)
(BSD-3-Clause, Lawrence Livermore National Security and Southern
Methodist University), the Rosenbrock4 stepper from
[Boost.Numeric.Odeint](https://www.boost.org/doc/libs/release/libs/numeric/odeint/)
(Boost Software License 1.0, Karsten Ahnert, Mario Mulansky and
Christoph Koke). Both licenses are permissive and compatible with MIT.
Using or redistributing cppDE requires preserving their notices in
[`inst/COPYRIGHTS`](https://github.com/dModverse/cppDE/blob/master/inst/COPYRIGHTS),
which lists the files concerned and reproduces the license texts, and
the attribution blocks in the headers.

The optional libraries linked at run time are not vendored: SUNDIALS
(BSD-3-Clause) and SuiteSparse/KLU (LGPL-2.1+) are governed by their own
licenses at the site where they are installed.
