# Forget Which Shared Object Exports Which Entry Point

[`cppFUN()`](https://dmodverse.github.io/cppDE/reference/cppFUN.md)
objects and
[`solveODE()`](https://dmodverse.github.io/cppDE/reference/solveODE.md)
remember the shared object that exports each generated entry point.
Recompiling a model into a shared object of another name makes that
pairing stale; call this to drop it.
[`compile()`](https://dmodverse.github.io/cppDE/reference/compile.md)
calls it itself, and loading or unloading a shared object does not need
it.

## Usage

``` r
clearNativeSymbols()
```

## Value

`NULL`, invisibly.

## Examples

``` r
## After a model was recompiled into a shared object of another name
clearNativeSymbols()
```
