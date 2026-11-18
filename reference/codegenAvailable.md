# Can Code Be Generated and Compiled Here?

Reports whether
[`cppODE()`](https://dmodverse.github.io/cppDE/reference/cppODE.md),
[`cvode()`](https://dmodverse.github.io/cppDE/reference/cvode.md) and
[`cppFUN()`](https://dmodverse.github.io/cppDE/reference/cppFUN.md) can
generate and compile a model without installing anything first. Code
generation needs Python with SymPy, which reticulate provisions on first
use and may download; this check never triggers that.

## Usage

``` r
codegenAvailable(backend = c("native", "cvode"))
```

## Arguments

- backend:

  `"native"` (default) for
  [`cppODE()`](https://dmodverse.github.io/cppDE/reference/cppODE.md)
  and
  [`cppFUN()`](https://dmodverse.github.io/cppDE/reference/cppFUN.md),
  `"cvode"` to also require the SUNDIALS libraries of the
  [`cvode()`](https://dmodverse.github.io/cppDE/reference/cvode.md)
  backend.

## Value

`TRUE` when a C++ compiler and SymPy are found without installing
anything, `FALSE` otherwise.

## Examples

``` r
codegenAvailable()
#> [1] TRUE
codegenAvailable("cvode")
#> [1] TRUE
```
