# Prepare a Batch for Repeated Solving

Validates and converts a set of conditions once, so that repeated solves
with
[`solveBatch()`](https://dmodverse.github.io/cppDE/reference/solveBatch.md),
an optimizer evaluating the same model at new parameters, only pass the
numbers that changed.
[`solveODEBatch()`](https://dmodverse.github.io/cppDE/reference/solveODEBatch.md)
repeats that work on every call.

## Usage

``` r
prepareBatch(
  model,
  conditions,
  times = NULL,
  parms = NULL,
  tangent = NULL,
  hessian = NULL,
  fixed = NULL,
  forcings = NULL,
  abstol = 1e-06,
  reltol = 1e-06,
  maxattempts = 50L,
  maxsteps = 1000000L,
  hini = 0,
  roottol = 1e-06,
  maxroot = 1L,
  cotangent = NULL,
  curvature = NULL,
  keepStore = FALSE,
  store = NULL,
  sensErrCon = TRUE,
  adjoint = NULL,
  maxOrder = NULL,
  stald = FALSE
)
```

## Arguments

- model:

  A model handle from
  [`cppODE()`](https://dmodverse.github.io/cppDE/reference/cppODE.md) or
  [`cvode()`](https://dmodverse.github.io/cppDE/reference/cvode.md).

- conditions:

  A named or unnamed list of per-condition argument lists. Recognized
  names are `times`, `parms`, `tangent`, `hessian`, `cotangent`,
  `curvature`, `fixed`, `forcings`, `abstol`, `reltol`, `maxattempts`,
  `maxsteps`, `hini`, `roottol`, `maxroot`, `keepStore`, `store`,
  `sensErrCon`, `adjoint`, `maxOrder` and `stald`; each overrides the
  batch-wide argument of the same name.

- times, parms:

  Batch-wide `times` and `parms` as in
  [`solveODE()`](https://dmodverse.github.io/cppDE/reference/solveODE.md),
  for the conditions that do not give their own. Default `NULL`.

- tangent:

  Optional numeric matrix, the tangent of the inputs; shapes in Details.
  Default `NULL`: the identity on the active sensitivities.

- hessian:

  Optional numeric array, the second derivative of the inputs; shapes in
  Details. Needs a model with `attr(model, "deriv2")`. Default `NULL`:
  the inputs are affine in \\\theta\\.

- fixed:

  Optional character vector, sensitivity names held fixed in this call,
  a subset of `attr(model, "dimNames")$sens`. Unlike `fixed` of
  [`cppODE()`](https://dmodverse.github.io/cppDE/reference/cppODE.md) it
  changes without recompiling. Not with a full or partial `tangent`.
  Default `NULL`.

- forcings:

  Named list of forcing data, one element per name in
  `attr(model, "forcings")`: a `data.frame` with columns `time` and
  `value`, or a two-column matrix.
  [`cppODE()`](https://dmodverse.github.io/cppDE/reference/cppODE.md)
  and [`cvode()`](https://dmodverse.github.io/cppDE/reference/cvode.md)
  take the forcing names, this argument the data. A forcing is the
  monotone cubic Hermite interpolant (PCHIP) of the data points and
  constant beyond the first and the last point, see
  [`forcingValues()`](https://dmodverse.github.io/cppDE/reference/forcingValues.md).
  Default `NULL`.

- abstol:

  Absolute error tolerance. Default `1e-6`.

- reltol:

  Relative error tolerance. Default `1e-6`.

- maxattempts:

  Maximum number of consecutive step attempts without progress in time
  before the solve stops with `return_code = -4`. Default `50`.

- maxsteps:

  Maximum total number of integration steps. Default `1e6`.

- hini:

  Initial step size; `0` (default) triggers automatic estimation.

- roottol:

  Tolerance for locating root events. Default `1e-6`. Under
  `rootfunc = "equilibrate"` the solve stops once every \\\|\dot x_i\|
  \le \mathrm{roottol}\\ \|x_i\| + \mathrm{abstol}\\; on a
  [`cppODE()`](https://dmodverse.github.io/cppDE/reference/cppODE.md)
  model the test covers the sensitivities as well.

- maxroot:

  Maximum number of triggers per root event. Default `1`.

- cotangent:

  The cotangent of the outputs, `[n_t, n_x]` or `[n_t, n_x, n_seed]`.
  Required by a model compiled with `derivMode = "reverse"` or
  `"forward-reverse"` unless `keepStore = TRUE`; an error on any other
  model. `result$cotangent` is then \\\sum\_{t,i} w\_{t,i}\\ \partial
  x_i(t) / \partial \theta\\, one column per seed. Default `NULL`.

- curvature:

  Optional, `"forward-reverse"` only: the derivative of `cotangent`
  along the tangent, `[n_t, n_x, n_seed, n_s]`. For a cotangent that is
  the gradient of a functional of the outputs, this is that functional's
  Hessian applied to the output tangent. Default `NULL`: the cotangent
  is constant in \\\theta\\.

- keepStore:

  Logical, default `FALSE`. Return the checkpoints of a reverse solve as
  `$store`, for a later solve to reuse through `store`. The `cotangent`
  may then be omitted, and the solve returns the values alone. Only on a
  [`cppODE()`](https://dmodverse.github.io/cppDE/reference/cppODE.md)
  model compiled with `derivMode = "reverse"`.

- store:

  The `$store` of an earlier solve of the same model at the same `times`
  and `parms`. The solve then integrates nothing and runs the backward
  sweep alone. A store from a different point is an error. A store may
  be reused any number of times. Default `NULL`.

- sensErrCon:

  Logical. Include the sensitivities, under
  `derivMode = "forward-reverse"` the tangents, in the error test and
  the first step. Default `TRUE`. `FALSE` takes the steps and corrector
  iterations of a value run, up to rounding: cheaper and less accurate,
  the sensitivities are the derivative of the computed solution, as the
  reverse mode of `"adams"` returns; that of `"bdf"` solves the
  corrector equations instead. Under `"forward-reverse"` the step grid
  does not depend on the directions of the tangent. Needs a model with
  sensitivities, has no effect under `derivMode = "reverse"`; a
  [`cvode()`](https://dmodverse.github.io/cppDE/reference/cvode.md)
  model takes `TRUE` only.

- adjoint:

  Optional
  [`adjointControl()`](https://dmodverse.github.io/cppDE/reference/adjointControl.md)
  object for a model compiled with `derivMode = "reverse"`. Default
  `NULL`, the same as
  [`adjointControl()`](https://dmodverse.github.io/cppDE/reference/adjointControl.md).

- maxOrder:

  Optional upper bound on the order of a multistep method, 1 to 5 for
  `"bdf"` and 1 to 12 for `"adams"`. Default `NULL`, the method's own
  maximum. An error on `"rb4"` and `"tsit5"`.

- stald:

  Logical, default `FALSE`. Stability limit detection for method
  `"bdf"`: the order drops when the step size leaves the stability
  region of the current order. See [Stepping and error
  control](https://dmodverse.github.io/cppDE/articles/methods-stepping.html#adaptive-step-size-control).

## Value

An object of class `"cppDEbatch"` for
[`solveBatch()`](https://dmodverse.github.io/cppDE/reference/solveBatch.md).

## See also

[`solveBatch()`](https://dmodverse.github.io/cppDE/reference/solveBatch.md),
[`solveODEBatch()`](https://dmodverse.github.io/cppDE/reference/solveODEBatch.md)

## Examples

``` r
# \donttest{
if (codegenAvailable()) {
  ## Marshal the conditions once, then re-solve with new numbers, as inside an
  ## optimizer.
  f <- cppODE(c(x = "-k*x"), deriv = FALSE)
  conditions <- list(a = list(parms = c(x = 1, k = 0.1)),
                     b = list(parms = c(x = 2, k = 0.1)))
  h   <- prepareBatch(f, conditions, times = 0:5)
  res <- solveBatch(h, parms = list(c(x = 1, k = 0.5), c(x = 2, k = 0.5)))
  res$b$variable
}
#> using C++ compiler: g++ [-O2 -DNDEBUG -w -fPIC -fopenmp]
#>              x
#> [1,] 2.0000000
#> [2,] 1.2130658
#> [3,] 0.7357696
#> [4,] 0.4462715
#> [5,] 0.2706804
#> [6,] 0.1641778
# }
```
