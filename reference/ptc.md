# Steady State by Pseudo-Transient Continuation

Solves \\f(x, p) = 0\\ for variables of a
[`cppFUN()`](https://dmodverse.github.io/cppDE/reference/cppFUN.md)
model, one output per solved variable, optionally under linear
constraints \\C x = total\\. `flow = TRUE` follows the ODE flow to a
stable steady state; `flow = FALSE` returns any regular root. The result
depends only on the start, the parameters and the controls.

## Usage

``` r
ptc(
  fun,
  x,
  parms,
  solve = names(x),
  rows = solve,
  C = NULL,
  total = NULL,
  flow = TRUE,
  positive = TRUE,
  controls = list()
)
```

## Arguments

- fun:

  A compiled
  [`cppFUN()`](https://dmodverse.github.io/cppDE/reference/cppFUN.md)
  object with forward derivatives.

- x:

  Named numeric, values of all variables of `fun`. The solved variables
  start there, the others keep their values.

- parms:

  Named numeric, values of the parameters of `fun`.

- solve:

  Names of the variables to solve for. Default all.

- rows:

  Names of the outputs of `fun` paired with `solve`, in the same order.
  Default `solve`.

- C, total:

  Constraint matrix with columns named by `solve`, and the right-hand
  side: `C x = total`. Default none.

- flow:

  Logical, default `TRUE`: the outputs are the right-hand side of an
  ODE. `FALSE` for plain equations.

- positive:

  Logical, default `TRUE`. Keep the solved variables positive; they must
  start positive.

- controls:

  Named list:

  - `reltol` (default `1e-10`) and `abstol` (`1e-14`): each residual is
    measured against `reltol` times the turnover \\\sum_k \|\partial f_i
    / \partial x_k\|\\ \|x_k\|\\ plus `abstol` times the largest
    Jacobian entry of the residual. Newton steps from there, while they
    contract, bring each variable to within `reltol * (abs(x) + abstol)`
    of the root.

  - `flowTol` (`1`): relative local error of a pseudo-time step.

  - `maxit` (`400`): maximum number of iterations of one descent; a
    continuation from an unstable root starts another.

  - `dtInit` (`1e-2`): the first pseudo-time step, relative to the
    fastest rate if `flow = TRUE`.

  - `zmax` (`5`): the largest change of the logarithm of a variable per
    step.

  - `stability` (`TRUE`): with `flow = TRUE`, a root must be stable on
    the kernel of `C`. From an unstable root the search continues along
    the unstable direction, first to the side the iterates came from.

  - `stabilityTol` (`0`): with `stability`, a root counts as stable
    while the real part of the leading Jacobian eigenvalue on the kernel
    of `C` is at most `stabilityTol`. A small positive value accepts a
    state that leaves the root more slowly than a time scale of
    interest, such as near a Hopf bifurcation.

## Value

A list with the named solution `x` (all variables), `converged`,
`iterations` and `message`.

## References

Kelley CT, Keyes DE (1998). Convergence analysis of pseudo-transient
continuation. SIAM J Numer Anal 35(2):508-523.

## See also

[`cppFUN()`](https://dmodverse.github.io/cppDE/reference/cppFUN.md)

## Examples

``` r
# \donttest{
if (codegenAvailable()) {
  ## Enzyme E converts S to P, P returns to S. E + ES and S + ES + P are
  ## conserved, so the steady state is fixed by the two totals.
  f <- c(E  = "-kon * E * S + (koff + kcat) * ES",
         S  = "-kon * E * S + koff * ES + kback * P",
         ES = " kon * E * S - (koff + kcat) * ES",
         P  = " kcat * ES - kback * P")
  enz <- cppFUN(f, variables = names(f), parameters = c("kon", "koff", "kcat", "kback"),
                modelname = "ptc_enzyme", outdir = tempdir(), compile = TRUE)
  C <- rbind(Etot = c(E = 1, S = 0, ES = 1, P = 0),
             Stot = c(E = 0, S = 1, ES = 1, P = 1))
  ss <- ptc(enz, x = c(E = 1, S = 1, ES = 1, P = 1),
            parms = c(kon = 1, koff = 0.5, kcat = 0.3, kback = 0.1),
            C = C, total = c(1, 10))
  ss$x
  enz$func(E = ss$x[["E"]], S = ss$x[["S"]], ES = ss$x[["ES"]], P = ss$x[["P"]],
           kon = 1, koff = 0.5, kcat = 0.3, kback = 0.1)
}
#> using C++ compiler: g++ [-O2 -DNDEBUG -w -fPIC -fopenmp]
#>                  E            S           ES             P
#> [1,] -3.379727e-16 1.287859e-14 3.379727e-16 -1.329308e-14
# }
```
