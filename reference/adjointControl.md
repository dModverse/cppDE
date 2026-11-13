# Controls of a Reverse Solve

Options for the backward sweep of a model compiled with
`derivMode = "reverse"`: check the sweep against the tolerances and
report the grid it ran on. Handed to
[`solveODE()`](https://dmodverse.github.io/cppDE/reference/solveODE.md),
[`solveODEBatch()`](https://dmodverse.github.io/cppDE/reference/solveODEBatch.md),
[`prepareBatch()`](https://dmodverse.github.io/cppDE/reference/prepareBatch.md)
and
[`solveBatch()`](https://dmodverse.github.io/cppDE/reference/solveBatch.md)
as `adjoint`.

## Usage

``` r
adjointControl(refine = FALSE, gradtol = NULL, trace = FALSE)
```

## Arguments

- refine:

  Logical, default `FALSE`. Check each step of the sweep against
  `reltol` and `abstol`, and with `gradtol` the contribution of the step
  to the gradient as well; a step that fails is swept again in substeps.
  The solve then also returns `$correction`, see
  [`solveODE()`](https://dmodverse.github.io/cppDE/reference/solveODE.md).

- gradtol:

  Positive number, the absolute tolerance on the gradient. Default
  `NULL`: no test on the gradient. On a
  [`cppODE()`](https://dmodverse.github.io/cppDE/reference/cppODE.md)
  model used only with `refine = TRUE`; on a
  [`cvode()`](https://dmodverse.github.io/cppDE/reference/cvode.md)
  model the absolute tolerance of the gradient quadrature in the
  backward problem's error test.

- trace:

  Logical, default `FALSE`. Return `$adjoint`.

## Value

A list of class `cppDEadjointControl`.

A solve with `trace` returns `$adjoint`, of class `cppDEadjoint`: `time`
and `h`, the start and length of each step the sweep ran on; `lambda`,
`[n_steps, n_states, n_seed]`, the adjoint state at each step's start;
`eta`, `[n_steps, n_seed]`, lambda times each step's local error
estimate, `NA` for `bdf` and `adams` under `refine`. Under `refine` also
`substeps`, the substeps each step was swept in, and `failures`, the
steps over all seed columns that failed the test at the largest substep
count.

## Details

`refine` and `trace` are an error on a
[`cvode()`](https://dmodverse.github.io/cppDE/reference/cvode.md) model.

## See also

[`solveODE()`](https://dmodverse.github.io/cppDE/reference/solveODE.md);
[Reverse mode: the discrete
adjoint](https://dmodverse.github.io/cppDE/articles/methods-reverse.html)
and [Grid and refinement of the
sweep](https://dmodverse.github.io/cppDE/articles/methods-reverse-checks.html)
on the cppDE site.

## Examples

``` r
# \donttest{
if (codegenAvailable()) {
  ## A fast forcing that the state does not see, but the derivative of the state in q does.
  rev <- cppODE(c(x = "-k*x + sin(100*time)*q"), derivMode = "reverse",
                modelname = "adjointControl_example")
  times <- seq(0, 2, 0.1)
  parms <- c(x = 1, k = 1, q = 0)
  W <- cbind(x = rep(1, length(times)))
  ## The reference is a refined sweep at tight tolerances: without refine, the
  ## sweep stays on the grid the state sets, however tight the tolerances.
  ref <- solveODE(rev, times, parms, cotangent = W, abstol = 1e-12, reltol = 1e-12,
                  adjoint = adjointControl(refine = TRUE))$cotangent[, 1]

  ## The sweep on the grid of the value run, and checked against the
  ## tolerances, with the gradient held to gradtol as well.
  val <- solveODE(rev, times, parms, cotangent = W, reltol = 1e-6)
  chk <- solveODE(rev, times, parms, cotangent = W, reltol = 1e-6,
                  adjoint = adjointControl(refine = TRUE, gradtol = 1e-9,
                                           trace = TRUE))
  print(rbind(value_grid = val$cotangent[, 1] - ref, checked = chk$cotangent[, 1] - ref))
  sum(chk$adjoint$substeps)
}
# }
```
