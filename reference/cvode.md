# Generate and Compile an ODE Solver Linked Against SUNDIALS CVODE(S)

Generates a C++ ODE solver linked against SUNDIALS CVODES, from the
system or from
[`installLibs()`](https://dmodverse.github.io/cppDE/reference/installLibs.md),
compiles it via `R CMD SHLIB`, and returns a handle for use with
[`solveODE()`](https://dmodverse.github.io/cppDE/reference/solveODE.md).
Sensitivities are computed by the CVODES forward sensitivity solver
(`deriv = TRUE`) or by the CVODES adjoint (`derivMode = "reverse"`). The
compiled model exposes the same R interface as a model from
[`cppODE()`](https://dmodverse.github.io/cppDE/reference/cppODE.md); the
differences between the two backends are described in [cvode() and the
CVODES
adjoint](https://dmodverse.github.io/cppDE/articles/methods-cvode.html).

## Usage

``` r
cvode(
  rhs,
  events = NULL,
  rootfunc = NULL,
  fixed = NULL,
  forcings = NULL,
  compile = TRUE,
  modelname = NULL,
  outdir = tempdir(),
  deriv = FALSE,
  derivMode = c("forward", "reverse"),
  asaCheckpoints = 200L,
  sparse = NULL,
  method = c("bdf", "adams"),
  includeTimeZero = TRUE,
  stepTrace = FALSE,
  verbose = FALSE
)
```

## Arguments

- rhs:

  Named character vector of ODE right-hand sides, see
  [expressions](https://dmodverse.github.io/cppDE/reference/expressions.md).
  Names are the state variables.

- events:

  Optional event `data.frame`, see Details. Default `NULL`.

- rootfunc:

  Optional integration-termination root: `"equilibrate"` or a character
  vector of expressions. Default `NULL`.

- fixed:

  Optional character vector of state or parameter names excluded from
  sensitivities at compile time. Default `NULL`.

- forcings:

  Optional character vector, the names of the forcings in `rhs`. Their
  data go to `forcings` of
  [`solveODE()`](https://dmodverse.github.io/cppDE/reference/solveODE.md).
  Default `NULL`.

- compile:

  Logical. Compile and load the generated C++ code. Default `TRUE`; with
  `FALSE`, compile several models together with
  [`compile()`](https://dmodverse.github.io/cppDE/reference/compile.md).

- modelname:

  Optional base name for the generated source file and exported C
  symbols. A random identifier is used when `NULL`.

- outdir:

  Directory for the generated C++ source. Default
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html).

- deriv:

  Logical. Compute first-order forward sensitivities. Default `FALSE`.

- derivMode:

  Direction the derivatives are taken in. `"forward"` (default) is the
  CVODES forward sensitivity solver, driven by `deriv`. `"reverse"` is
  the CVODES adjoint, one backward solve per cotangent column; it needs
  `deriv = FALSE`.

- asaCheckpoints:

  Number of forward steps between two checkpoints of the adjoint, under
  `derivMode = "reverse"`. Fewer steps take more memory and interpolate
  the forward state more accurately. Default `200`.

- sparse:

  Logical or `NULL`. `NULL` (default) selects sparse or dense LU from
  the Jacobian sparsity; `TRUE` forces sparse, `FALSE` forces dense.
  Sparse LU requires KLU at install time, see
  [`installLibs()`](https://dmodverse.github.io/cppDE/reference/installLibs.md).

- method:

  One of `"bdf"` (default) or `"adams"`.

- includeTimeZero:

  Logical. Ensure that `0` is part of the integration times, as
  [`cppODE()`](https://dmodverse.github.io/cppDE/reference/cppODE.md)
  does. Default `TRUE`.

- stepTrace:

  Logical, default `FALSE`. Record per-step diagnostics, returned as
  `$trace` by
  [`solveODE()`](https://dmodverse.github.io/cppDE/reference/solveODE.md);
  with `events` or `rootfunc` one row per output time and event.

- verbose:

  Logical. Print progress messages. Default `FALSE`.

## Value

The compiled model name (character) with the attributes of a model from
[`cppODE()`](https://dmodverse.github.io/cppDE/reference/cppODE.md),
except that `deriv2` is `FALSE` and `useNDF` is `NA`, plus
`backend = "cvode"` and `lapackDense`, whether dense Jacobians use
`SUNLinSol_LapackDense` rather than `SUNLinSol_Dense`.

## Details

Available methods are `"bdf"` (default) and `"adams"`. Forward
sensitivities are of first order only; `deriv2` is not supported.
Events, forcings, `rootfunc`, `fixed` and a switch on a state or the
time behave as in
[`cppODE()`](https://dmodverse.github.io/cppDE/reference/cppODE.md); the
switch is located by CVODES' root finding.

Needs SUNDIALS (\>= 6.0) at install time, see
[`installLibs()`](https://dmodverse.github.io/cppDE/reference/installLibs.md);
otherwise `cvode()` stops with platform-specific install instructions.
KLU, required when `sparse = TRUE`, is detected the same way.

## References

Hindmarsh, A. C., Brown, P. N., Grant, K. E., Lee, S. L., Serban, R.,
Shumaker, D. E., and Woodward, C. S. (2005). SUNDIALS: Suite of
Nonlinear and Differential/Algebraic Equation Solvers. *ACM Transactions
on Mathematical Software* **31**(3), 363-396.

## See also

[expressions](https://dmodverse.github.io/cppDE/reference/expressions.md),
[`cppODE()`](https://dmodverse.github.io/cppDE/reference/cppODE.md),
[`solveODE()`](https://dmodverse.github.io/cppDE/reference/solveODE.md),
[`cppFUN()`](https://dmodverse.github.io/cppDE/reference/cppFUN.md);
[cvode() and the CVODES
adjoint](https://dmodverse.github.io/cppDE/articles/methods-cvode.html)
on the cppDE site.

## Examples

``` r
# \donttest{
if (codegenAvailable("cvode")) {
  ## Needs SUNDIALS at install time, see installLibs(). The models are
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
#> using C++ compiler: g++ [-O2 -DNDEBUG -w -fPIC -I/usr/lib/x86_64-linux-gnu/openmpi/include -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi -fopenmp]
#>          [,1]
#> A    3.957127
#> B    4.905396
#> k1   1.036174
#> k2 -12.529511
# }
```
