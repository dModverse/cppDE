# Generate and Compile a C++ ODE Solver

Generates C++ source code for an ODE system \\\dot{x} = f(x, p)\\,
compiles it via `R CMD SHLIB`, and returns a handle for use with
[`solveODE()`](https://dmodverse.github.io/cppDE/reference/solveODE.md).
The four available methods (`"bdf"`, `"adams"`, `"rb4"`, `"tsit5"`) and
their selection criteria, the dual-number sensitivity framework, and the
event handling are described in the articles [Stepping and error
control](https://dmodverse.github.io/cppDE/articles/methods-stepping.html),
[Forward
sensitivities](https://dmodverse.github.io/cppDE/articles/methods-forward.html)
and [Events, switches and
saltation](https://dmodverse.github.io/cppDE/articles/methods-events.html).

## Usage

``` r
cppODE(
  rhs,
  events = NULL,
  rootfunc = NULL,
  fixed = NULL,
  forcings = NULL,
  compile = TRUE,
  modelname = NULL,
  outdir = tempdir(),
  deriv = TRUE,
  deriv2 = FALSE,
  includeTimeZero = TRUE,
  useDenseOutput = TRUE,
  sparse = NULL,
  method = c("bdf", "adams", "rb4", "tsit5"),
  useNDF = TRUE,
  derivMode = c("forward", "reverse", "forward-forward", "forward-reverse"),
  profile = FALSE,
  stepTrace = FALSE,
  verbose = FALSE
)
```

## Arguments

- rhs:

  Named character vector of ODE right-hand sides, see
  [expressions](https://dmodverse.github.io/cppDE/reference/expressions.md).
  Names are the state variables.

- events:

  Optional event `data.frame`, see Details. Default `NULL`.

- rootfunc:

  Optional integration-termination root: `"equilibrate"` or a character
  vector of expressions. Default `NULL`.

- fixed:

  Optional character vector of state or parameter names excluded from
  sensitivities at compile time. Default `NULL`.

- forcings:

  Optional character vector, the names of the forcings in `rhs`. Their
  data go to `forcings` of
  [`solveODE()`](https://dmodverse.github.io/cppDE/reference/solveODE.md).
  Default `NULL`.

- compile:

  Logical. Compile and load the generated C++ code. Default `TRUE`; with
  `FALSE`, compile several models together with
  [`compile()`](https://dmodverse.github.io/cppDE/reference/compile.md).

- modelname:

  Optional base name for the generated source file and exported C
  symbols. A random identifier is used when `NULL`.

- outdir:

  Directory for the generated C++ source. Default
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html).

- deriv:

  Logical. Compute first-order parameter sensitivities. Default `TRUE`.
  Ignored under `derivMode = "reverse"` and `"forward-reverse"`.

- deriv2:

  Logical. Compute second-order parameter sensitivities; implies
  `deriv = TRUE`. Default `FALSE`. `derivMode = "forward-forward"` is
  the preferred spelling.

- includeTimeZero:

  Logical. Ensure that `0` is part of the integration times. Default
  `TRUE`.

- useDenseOutput:

  Logical. Use the dense output of the method for user-requested time
  points. Default `TRUE`. Applies to `method = "rb4"` and `"tsit5"`,
  whose interpolant is of lower order than the method; `FALSE` makes the
  solver land on each requested time exactly. `FALSE` is ignored with a
  warning for `"bdf"` and `"adams"` and under the reverse modes.

- sparse:

  Logical or `NULL`. `NULL` (default) selects sparse or dense LU from
  the Jacobian sparsity; `TRUE` forces sparse, `FALSE` forces dense.
  Sparse LU requires KLU at install time, see
  [`installLibs()`](https://dmodverse.github.io/cppDE/reference/installLibs.md).

- method:

  Integration method: `"bdf"` (default), `"adams"`, `"rb4"`, or
  `"tsit5"`.

- useNDF:

  Logical. Use the NDF coefficients of Shampine and Reichelt in the BDF
  corrector. Default `TRUE`. Applies to `method = "bdf"`; ignored
  otherwise.

- derivMode:

  Direction the derivatives are taken in.

  - `"forward"` (default): forward sensitivities as built by `deriv`,
    returned as `$tangent`.

  - `"reverse"`: no forward sensitivities; the model takes a `cotangent`
    in
    [`solveODE()`](https://dmodverse.github.io/cppDE/reference/solveODE.md)
    and returns `$cotangent`, one row per state and parameter.

  - `"forward-forward"`: the same as `deriv2 = TRUE`, returning
    `$hessian`.

  - `"forward-reverse"`: returns `$tangent`, `$cotangent` and
    `$curvature`, the derivatives of `$cotangent` along each tangent
    direction.

- profile:

  Logical, default `FALSE`. Compile with timers that print a summary per
  category to the standard error stream after each solve.

- stepTrace:

  Logical, default `FALSE`. Compile to record per-step diagnostics,
  returned as `$trace` by
  [`solveODE()`](https://dmodverse.github.io/cppDE/reference/solveODE.md).

- verbose:

  Logical. Print progress messages. Default `FALSE`.

## Value

The compiled model name (character) with the attributes required by
[`solveODE()`](https://dmodverse.github.io/cppDE/reference/solveODE.md):
`equations`, `srcfile`, `variables`, `parameters`, `forcings`, `events`,
`rootfunc`, `fixed` (the compile-time set), `deriv`, `deriv2`,
`derivMode`, `sparse`, `method`, `useNDF`, `dimNames` (with element
`sens` only when `deriv` is `TRUE`), `compileArgs`, and `linkArgs` when
OpenMP is available.

## Details

Events are specified as a `data.frame` with columns `var` (affected
variable), `value` (expression to apply), `method` (`"replace"`,
`"add"`, or `"multiply"`), and exactly one of `time` (event time) or
`root` (root expression). `rootfunc` terminates integration:
`"equilibrate"` stops at steady state (see `roottol` in
[`solveODE()`](https://dmodverse.github.io/cppDE/reference/solveODE.md));
a character vector of expressions stops at the first zero crossing.

A time event fires when the event time lies between the first and the
last integration time, the first included and the last not; with
`includeTimeZero` the integration times contain 0. An event at the first
time is applied before the solve starts, and one inside the window adds
the event time to the output, with the state after the event in that
row. An event before the first or at or after the last time is ignored
and changes neither the output nor the span of the integration. A grid
of a single time applies the events at that time.

A comparison `<`, `<=`, `>`, `>=` whose sides differ by an amount that
reads a state or the time, inside `piecewise()`, a logical operator or a
product, and `Heaviside()` or
[`sign()`](https://rdrr.io/r/base/sign.html) of such an amount, switch
the right-hand side where the amount crosses zero, for instance
`"piecewise(1, time > ts && time <= t2, 0)"` or `"kin - v*(x > 0)"`. The
solver locates each crossing as a root, so a short pulse is not stepped
over, and the sensitivities take the jump. A switch adds no output row
and is not limited by `maxroot`. Each branch is evaluated a step beyond
the switch, so it has to be defined there. Where both branches push the
solution onto the switching surface, it slides along it. [Events,
switches and
saltation](https://dmodverse.github.io/cppDE/articles/methods-events.html#sliding)
describes how crossings and slides are solved. `==` and `!=` are
evaluated as written.

## See also

[expressions](https://dmodverse.github.io/cppDE/reference/expressions.md)
for the syntax of equations;
[`solveODE()`](https://dmodverse.github.io/cppDE/reference/solveODE.md)
for integration;
[`cvode()`](https://dmodverse.github.io/cppDE/reference/cvode.md) for
the SUNDIALS-backed alternative;
[`cppFUN()`](https://dmodverse.github.io/cppDE/reference/cppFUN.md) for
algebraic functions; [Code
generation](https://dmodverse.github.io/cppDE/articles/methods-codegen.html),
[Stepping and error
control](https://dmodverse.github.io/cppDE/articles/methods-stepping.html),
[Forward
sensitivities](https://dmodverse.github.io/cppDE/articles/methods-forward.html),
[Events, switches and
saltation](https://dmodverse.github.io/cppDE/articles/methods-events.html)
[Reverse mode: the discrete
adjoint](https://dmodverse.github.io/cppDE/articles/methods-reverse.html)
and [Events in reverse
mode](https://dmodverse.github.io/cppDE/articles/methods-reverse-events.html)
on the cppDE site.

## Examples

``` r
# \donttest{
if (codegenAvailable()) {
  ## Six variants of one decay chain A -> B. Each is generated first and all
  ## are compiled together afterwards, into one shared object in tempdir().
  eqns <- c(A = "-k1*A", B = "k1*A - k2*B")

  ## First-order sensitivities, the default
  fwd <- cppODE(eqns, compile = FALSE)

  ## A dose of A at the parameter time t_e
  events <- data.frame(var = "A", time = "t_e", value = 1, method = "add")
  dosed  <- cppODE(eqns, events = events, compile = FALSE)

  ## A constant inflow k0, integrated until the steady state
  steady <- cppODE(c(A = "k0 - k1*A", B = "k1*A - k2*B"),
                   rootfunc = "equilibrate", compile = FALSE)

  ## A time-dependent input u, given as data when solving
  forced <- cppODE(c(A = "u - k1*A", B = "k1*A - k2*B"), forcings = "u",
                   compile = FALSE)

  ## Second-order sensitivities
  fwd2 <- cppODE(eqns, deriv2 = TRUE, compile = FALSE)

  ## Reverse mode: the states in plain double, the derivatives from one
  ## backward sweep whose cost does not grow with the number of parameters
  rev <- cppODE(eqns, derivMode = "reverse", compile = FALSE)

  compile(fwd, dosed, steady, forced, fwd2, rev, output = "cppODE_examples")

  times <- seq(0, 10, by = 0.5)
  parms <- c(A = 1, B = 0, k1 = 0.5, k2 = 0.2)

  out <- solveODE(fwd, times, parms)
  head(out$variable)
  out$tangent[21, "B", ]                # dB(10) / d(A, B, k1, k2)

  out <- solveODE(dosed, times, c(parms, t_e = 4))
  out$variable[times %in% c(3.5, 4, 4.5), ]

  out <- solveODE(steady, times = c(0, 1000), c(parms, k0 = 1))
  tail(out$time, 1)                     # stopped before t = 1000
  tail(out$variable, 1)                 # the steady state

  u   <- list(u = data.frame(time = times, value = 1 + sin(times)))
  out <- solveODE(forced, times, parms, forcings = u)
  out$variable[21, ]

  out <- solveODE(fwd2, times, parms)
  out$hessian[21, "B", , ]              # Hessian of B(10)

  ## The gradient of sum(B) over the grid: a cotangent of ones on B
  onB <- cbind(A = 0, B = rep(1, length(times)))
  solveODE(rev, times, parms, cotangent = onB)$cotangent
}
# }
```
