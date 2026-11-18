# Print Solver Diagnostics

Prints a summary of the solver diagnostics returned by
[`solveODE()`](https://dmodverse.github.io/cppDE/reference/solveODE.md).

## Usage

``` r
diagnostics(result)
```

## Arguments

- result:

  A list returned by
  [`solveODE()`](https://dmodverse.github.io/cppDE/reference/solveODE.md),
  containing a `diagnostics` element.

## Value

The `diagnostics` list, invisibly; the fields are listed under Value in
[`solveODE()`](https://dmodverse.github.io/cppDE/reference/solveODE.md).
Without diagnostics, a message and `NULL`.

## Examples

``` r
# \donttest{
if (codegenAvailable()) {
  f <- cppODE(c(x = "-k*x"), deriv = FALSE)
  diagnostics(solveODE(f, times = 0:5, parms = c(x = 1, k = 0.3)))
}
#> using C++ compiler: g++ [-O2 -DNDEBUG -w -fPIC -fopenmp]
#> ----------------------  NDF solver statistics  ----------------------
#>   Return code                  : 0
#>   Message                      : Integration was successful.
#> ---------------------------------------------------------------------
#>   Accepted steps               : 32
#>   Rejected steps               : 0
#>   Function evaluations         : 43
#>   Jacobian evaluations         : 1
#>   LU factorizations            : 8
#>   Last step size (successful)  : 0.501632
#>   Last method order            : 5
#>   Time reached                 : 5
#> ---------------------------------------------------------------------
# }
```
