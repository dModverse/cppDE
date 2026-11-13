# Evaluate Several Compiled Functions at Once

Evaluates argument sets, one per
[`cppFUN()`](https://dmodverse.github.io/cppDE/reference/cppFUN.md)
object, all sets in one call from R and one parallel region.

## Usage

``` r
evaluateMany(funs, sets, cores = 1L, deriv2 = FALSE)
```

## Arguments

- funs:

  List of the `evaluateBatch` entries of
  [`cppFUN()`](https://dmodverse.github.io/cppDE/reference/cppFUN.md)
  objects built with forward derivatives, one per set.

- sets:

  List of argument sets, as `evaluateBatch` takes them.

- cores:

  Number of threads over the sets.

- deriv2:

  Logical. Second derivatives as well. Default `FALSE`.

## Value

A list with one element per set, as `evaluateBatch` returns them.

## See also

[`cppFUN()`](https://dmodverse.github.io/cppDE/reference/cppFUN.md)

## Examples

``` r
# \donttest{
if (codegenAvailable()) {
  ## Two functions of one variable, compiled together and evaluated on separate
  ## argument sets; the parameter seeds give the derivatives in a and b.
  f1 <- cppFUN(c(y = "a * x"), variables = "x", parameters = "a",
               modelname = "many_a", outdir = tempdir())
  f2 <- cppFUN(c(y = "exp(-b * x)"), variables = "x", parameters = "b",
               modelname = "many_b", outdir = tempdir())
  compile(f1, f2, output = "many_example")
  x <- matrix(c(0, 1, 2), 3, 1, dimnames = list(NULL, "x"))
  seed <- function(p) matrix(1, 1, 1, dimnames = list(p, p))
  sets <- list(list(vars = x, params = c(a = 2), tangentP = seed("a")),
               list(vars = x, params = c(b = 0.5), tangentP = seed("b")))
  out <- evaluateMany(list(f1$evaluateBatch, f2$evaluateBatch), sets)
  out[[2]]$y
  out[[2]]$tangent[, , "b"]
}
# }
```
