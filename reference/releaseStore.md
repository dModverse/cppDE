# Release the Checkpoints of a Reverse Solve

Frees the checkpoints a solve kept with `keepStore = TRUE` at once.
Otherwise they are freed at the next garbage collection, which does not
count their memory. A solve given a released store integrates again.

## Usage

``` r
releaseStore(store)
```

## Arguments

- store:

  The `$store` of a solve, or a list of them; `NULL` entries are
  skipped.

## Value

Invisibly, the number of stores freed.

## See also

[`solveODE()`](https://dmodverse.github.io/cppDE/reference/solveODE.md)

## Examples

``` r
# \donttest{
if (codegenAvailable()) {
  rev <- cppODE(c(A = "-k1*A", B = "k1*A - k2*B"), derivMode = "reverse",
                modelname = "releaseStore_example")
  parms <- c(A = 1, B = 0, k1 = 0.5, k2 = 0.2)
  first <- solveODE(rev, 0:10, parms, keepStore = TRUE)
  onB <- cbind(A = 0, B = rep(1, 11))
  solveODE(rev, 0:10, parms, cotangent = onB, store = first$store)$cotangent
  ## releaseStore() frees the checkpoints; a solve given the store integrates again
  releaseStore(first$store)
  solveODE(rev, 0:10, parms, cotangent = onB, store = first$store)$cotangent
}
# }
```
