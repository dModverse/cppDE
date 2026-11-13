# Will a Batch Run in Parallel?

[`solveODEBatch()`](https://dmodverse.github.io/cppDE/reference/solveODEBatch.md)
degrades to a serial loop without saying so. This reports the conditions
it depends on, so a caller can warn instead of silently running
serially.

## Usage

``` r
batchAvailable(model)
```

## Arguments

- model:

  A model handle from
  [`cppODE()`](https://dmodverse.github.io/cppDE/reference/cppODE.md) or
  [`cvode()`](https://dmodverse.github.io/cppDE/reference/cvode.md).

## Value

List with `symbol` (the batch entry resolved), `openmp` (cppDE was built
with OpenMP), `modelOpenmp` (this model's shared object was), and
`parallel` (all three). Two further fallbacks are decided at run time
and are not visible here: inside a forked child, and inside an enclosing
OpenMP region.

## See also

[`solveODEBatch()`](https://dmodverse.github.io/cppDE/reference/solveODEBatch.md)

## Examples

``` r
# \donttest{
if (codegenAvailable()) {
  f <- cppODE(c(x = "-k*x"), deriv = FALSE)
  batchAvailable(f)
}
# }
```
