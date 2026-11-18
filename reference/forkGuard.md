# Report the BLAS Fork Guard

Reports the guard cppDE installs so that a threaded BLAS does not
deadlock in a forked child, as under
[`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html). The
guard sets BLAS to one thread for the width of every `fork()` and
restores the previous count in the parent.

## Usage

``` r
forkGuard()
```

## Value

List with `api` (the BLAS whose thread count cppDE steers, `NA` when no
entry point resolved), `threads` (the current thread count, `NA` without
a getter) and `guard` (whether the `fork()` handler is installed).

## Details

BLAS also runs on one thread during each solve. Outside a solve and a
`fork()` the thread count is the caller's own; a forked child keeps the
single thread it inherited.

`api` is `NA` on a BLAS with no runtime thread-control entry point, and
the deadlock is then still reachable; the fallback is
`OMP_NUM_THREADS=1` in the environment that starts R. `guard` is `FALSE`
on Windows, which has no `fork()`.

Attaching the package prints one line when BLAS runs more than one
thread or no entry point resolved. Set `options(cppDE.quiet = TRUE)`
beforehand to suppress the line, or attach with
[`suppressPackageStartupMessages()`](https://rdrr.io/r/base/message.html).

## See also

[`batchAvailable()`](https://dmodverse.github.io/cppDE/reference/batchAvailable.md)

## Examples

``` r
forkGuard()
#> $api
#> [1] "OpenBLAS"
#> 
#> $threads
#> [1] 4
#> 
#> $guard
#> [1] TRUE
#> 
```
