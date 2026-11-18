# Solve a Prepared Batch

Re-solves the conditions of a
[`prepareBatch()`](https://dmodverse.github.io/cppDE/reference/prepareBatch.md)
handle with new numbers. Only `parms`, `tangent`, `hessian`,
`cotangent`, `curvature`, `adjoint` and `store` may change; anything
else needs a fresh handle.

Each of these arguments is a list with one element per condition,
matched to the prepared conditions by position; names are ignored. A
`NULL` element keeps that condition's prepared value, and `NULL`
(default) for the whole argument keeps all of them.

## Usage

``` r
solveBatch(
  handle,
  parms = NULL,
  tangent = NULL,
  hessian = NULL,
  cotangent = NULL,
  curvature = NULL,
  adjoint = NULL,
  store = NULL,
  cores = NULL,
  traceFile = NULL,
  onFailure = c("stop", "warn", "silent")
)
```

## Arguments

- handle:

  A `"cppDEbatch"` object from
  [`prepareBatch()`](https://dmodverse.github.io/cppDE/reference/prepareBatch.md).

- parms:

  List of named numeric vectors.

- tangent, hessian, cotangent, curvature:

  Lists of arrays, each as the argument of the same name in
  [`solveODE()`](https://dmodverse.github.io/cppDE/reference/solveODE.md)
  and of the prepared shape. A `cotangent` may have another number of
  seed columns than the prepared one; rows and state columns stay.

- adjoint:

  List of
  [`adjointControl()`](https://dmodverse.github.io/cppDE/reference/adjointControl.md)
  objects.

- store:

  List of checkpoint stores, each the `$store` of an earlier reverse
  solve. Here a `NULL` element integrates that condition afresh.

- cores, traceFile, onFailure:

  As in
  [`solveODEBatch()`](https://dmodverse.github.io/cppDE/reference/solveODEBatch.md).

## Value

A list of
[`solveODE()`](https://dmodverse.github.io/cppDE/reference/solveODE.md)
results, named as the prepared conditions.

## See also

[`prepareBatch()`](https://dmodverse.github.io/cppDE/reference/prepareBatch.md),
[`solveODEBatch()`](https://dmodverse.github.io/cppDE/reference/solveODEBatch.md)

## Examples

``` r
# \donttest{
if (codegenAvailable()) {
  ## Marshal the conditions once, then re-solve with new numbers, as inside an
  ## optimizer.
  f <- cppODE(c(x = "-k*x"), deriv = FALSE, compile = FALSE)
  r <- cppODE(c(x = "-k*x"), derivMode = "reverse", compile = FALSE)
  compile(f, r, output = "solveBatch_example")
  conditions <- list(a = list(parms = c(x = 1, k = 0.1)),
                     b = list(parms = c(x = 2, k = 0.1)))
  h   <- prepareBatch(f, conditions, times = 0:5)
  res <- solveBatch(h, parms = list(c(x = 1, k = 0.5), c(x = 2, k = 0.5)))
  res$b$variable

  ## Reverse mode: the forward passes keep their checkpoints, and a prepared
  ## sweep reuses them for a new cotangent.
  fw <- solveODEBatch(r, lapply(conditions, c, keepStore = TRUE), times = 0:5)
  W  <- matrix(1, 6, 1, dimnames = list(NULL, "x"))
  hr <- prepareBatch(r, conditions, times = 0:5, cotangent = W)
  solveBatch(hr, cotangent = list(W, 2 * W),
             store = lapply(fw, `[[`, "store"))$b$cotangent
}
#> using C++ compiler: g++ [-O2 -DNDEBUG -w -fPIC -fopenmp]
#>         [,1]
#> x   9.482463
#> k -41.915020
# }
```
