# Run a Compiled ODE Model

Integrates a model compiled by
[`cppODE()`](https://dmodverse.github.io/cppDE/reference/cppODE.md) or
[`cvode()`](https://dmodverse.github.io/cppDE/reference/cvode.md) over
the time points `times`. Returns the state trajectory and, when the
model was compiled with derivatives, their tangent and Hessian (forward)
or the cotangent and curvature of a seeded functional (reverse).

## Usage

``` r
solveODE(
  model,
  times,
  parms,
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
  onFailure = c("stop", "warn", "silent"),
  traceFile = NULL,
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

  A compiled ODE model returned by
  [`cppODE()`](https://dmodverse.github.io/cppDE/reference/cppODE.md) or
  [`cvode()`](https://dmodverse.github.io/cppDE/reference/cvode.md).

- times:

  Numeric vector of output times, non-empty and finite.

- parms:

  Named numeric vector of initial conditions and parameters. Names must
  include all of
  `c(attr(model, "variables"), attr(model, "parameters"))`.

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

- onFailure:

  How to react when the solver returns a non-zero return code. One of
  `"stop"` (default; raise an error with the solver message and no
  partial results), `"warn"` (emit a warning and return partial results
  up to `t_reached`), or `"silent"` (return the partial result without
  any signal).

- traceFile:

  Optional path of a CSV file the step trace is written to. Needs a
  model compiled with `stepTrace = TRUE`, ignored otherwise. Default
  `NULL`.

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
  the first step. Default `TRUE`; `FALSE` is cheaper and less accurate,
  gives the steps and corrector iterations of a value run, so the
  sensitivities are the derivative of the computed solution and equal
  those of the reverse mode to rounding, and under `"forward-reverse"`
  gives a step grid that does not depend on the directions of the
  tangent. Needs a model with sensitivities, has no effect under
  `derivMode = "reverse"`, and a
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

A named list:

- `time`: the output times, length `n_t`.

- `variable`: the states, `[n_t, n_x]`.

- `tangent`: when `attr(model, "deriv")` is `TRUE`, `[n_t, n_x, n_s]`.
  The third dimension is named by the active sensitivities or the
  columns of `tangent`.

- `hessian`: when `attr(model, "deriv2")` is `TRUE`,
  `[n_t, n_x, n_s, n_s]`.

- `cotangent`: under `derivMode = "reverse"` or `"forward-reverse"`, the
  cotangent of the inputs, `[n_x + n_p, n_seed]`, rows as in a full
  `tangent`.

- `curvature`: under `"forward-reverse"`, `[n_x + n_p, n_s, n_seed]`,
  the derivative of each `cotangent` entry along the tangent. Under the
  identity tangent the columns are the Hessian of the seeded functional.
  The order of the last two dimensions is that of a Hessian, not that of
  the argument `curvature`.

- `correction`: under `derivMode = "reverse"` with
  `adjointControl(refine = TRUE)`, length `n_seed`, named as the seed
  columns: the first-order change of each seeded functional from the
  computed trajectory to the one the refined sweep differentiated. Value
  plus `correction` is the value whose gradient `cotangent` is. See [The
  value of the refined
  sweep](https://dmodverse.github.io/cppDE/articles/methods-reverse-checks.html#value-correction).

- `diagnostics`: solver statistics, printed by
  [`diagnostics()`](https://dmodverse.github.io/cppDE/reference/diagnostics.md):
  `return_code` (0 on success, negative codes as in SUNDIALS CVODE),
  `message`, `accepted` and `rejected` (step counts), `fevals`, `jevals`
  and `setups` (right-hand side evaluations, Jacobian evaluations,
  matrix factorizations), `last_dt` and `last_order` (size and order of
  the last successful step), `t_reached`, and the model's `method`,
  `useNDF` and, on a
  [`cvode()`](https://dmodverse.github.io/cppDE/reference/cvode.md)
  model, `backend`.

- `trace`: with `stepTrace = TRUE`, a `data.frame` with one row per step
  attempt: `nst` (accepted steps before the attempt), `t` (time at the
  end of the step), `h` (signed step size), `q` (order), `dsm` (error
  test value, the step passes below 1), `acnrm` and `acnrm_state`
  (weighted RMS norm of the local error over all components and over the
  states alone), `tq2` (error constant), `gamma` and `gamrat` (step size
  times the leading coefficient, and the ratio of that product to the
  value at the last factorization), `newton_conv` (1 when the corrector
  converged), `mode` (`"BDF"`, `"NDF"`, `"ADAMS"`, `"ONESTEP"`, or
  `"CVODE"` and `"CVODE_event"` on a
  [`cvode()`](https://dmodverse.github.io/cppDE/reference/cvode.md)
  model), `nfe`, `njev` and `nsetups` (cumulative counters) and
  `setup_reason` (what triggered the last factorization). A
  [`cvode()`](https://dmodverse.github.io/cppDE/reference/cvode.md)
  model reports `NaN` where CVODES exposes no value.

- `adjoint`: under `adjointControl(trace = TRUE)`, see
  [`adjointControl()`](https://dmodverse.github.io/cppDE/reference/adjointControl.md).

- `store`: with `keepStore = TRUE`, an external pointer to the
  checkpoints.

## Details

### Dimensions

`n_t` is the number of output rows, `nrow(result$variable)`:
`length(times)`, plus one row for each event inside the window and for
the time 0 that `includeTimeZero` adds. `n_x` is the number of states,
`n_p` the number of parameters, `n_s` the number of tangent directions
and `n_seed` the number of cotangent columns.

### Derivative arguments and results

A derivative argument and the corresponding result have the same name.
`tangent` and `hessian` are the first and second derivative in
\\\theta\\ of the inputs going in, and of the states coming out.
`cotangent` is the gradient of a functional with respect to the outputs
going in, and with respect to the inputs coming out. `curvature` is the
derivative of the cotangent along the tangent: that functional's Hessian
applied to the tangent, on the outputs going in and on the inputs coming
out.

### Shapes of `tangent` and `hessian`

`tangent` is the Jacobian \\\Phi'(\theta)\\ and `hessian` the second
derivative \\\Phi''(\theta)\\ of a reparametrization \\p =
\Phi(\theta)\\ of the initial states and parameters; the results are
then derivatives in \\\theta\\. `tangent` takes one of three shapes,
told apart by the row count and row names. `hessian` has the rows of
`tangent` and two dimensions of length `n_s`.

- State-only, `[n_x, n_s]`, rows unnamed or named by the states. The
  parameter rows are the identity on the active sensitivities, the
  model's sensitivity names minus `fixed`. `result$tangent[t, , ]` has
  this shape and can seed a following solve.

- Full, `[n_x + n_p, n_s]`, states first: \\\Phi'(\theta)\\ itself.
  `n_s` may change from call to call.

- Partial, `[k, n_s]` with `k < n_x + n_p`. The row names, a subset of
  `c(variables, parameters)`, place the rows; missing rows are zero, so
  those inputs are fixed.

Run-time `fixed` takes the state-only shape only. Column names, when
given, are the active sensitivity names (state-only) or name the
directions \\\theta\\ (full and partial); unnamed full or partial
columns are named `theta1`, `theta2`, and so on.

## See also

[`cppODE()`](https://dmodverse.github.io/cppDE/reference/cppODE.md) and
[`cvode()`](https://dmodverse.github.io/cppDE/reference/cvode.md) for
model compilation;
[`solveODEBatch()`](https://dmodverse.github.io/cppDE/reference/solveODEBatch.md)
for many conditions in one call;
[`adjointControl()`](https://dmodverse.github.io/cppDE/reference/adjointControl.md);
[`forcingValues()`](https://dmodverse.github.io/cppDE/reference/forcingValues.md);
[`diagnostics()`](https://dmodverse.github.io/cppDE/reference/diagnostics.md)
for printing solver statistics. [Stepping and error
control](https://dmodverse.github.io/cppDE/articles/methods-stepping.html),
[Forward
sensitivities](https://dmodverse.github.io/cppDE/articles/methods-forward.html),
[Reverse mode: the discrete
adjoint](https://dmodverse.github.io/cppDE/articles/methods-reverse.html)
and [Second derivatives: forward over
reverse](https://dmodverse.github.io/cppDE/articles/methods-forward-reverse.html)
on the cppDE site; [Installation and optional
libraries](https://dmodverse.github.io/cppDE/articles/installation.html)
for diagnosing a failing solve.

## Examples

``` r
# \donttest{
if (codegenAvailable()) {
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
# }
```
