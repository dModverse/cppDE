# Will a Batch Run in Parallel?

Reports whether
[`solveODEBatch()`](https://dmodverse.github.io/cppDE/reference/solveODEBatch.md)
solves the conditions of a model in parallel. Without a batch entry
point or OpenMP it runs a serial loop and gives no warning.

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
#> using C++ compiler: g++ [-O2 -DNDEBUG -w -fPIC -fopenmp]
#> $symbol
#> [1] TRUE
#> 
#> $openmp
#> [1] TRUE
#> 
#> $modelOpenmp
#> [1] TRUE
#> 
#> $parallel
#> [1] TRUE
#> 
# }
```
