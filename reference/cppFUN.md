# Compile Algebraic Functions with Optional Derivatives

Generates and compiles C++ code that evaluates a system of algebraic
expressions \\y = g(x, p)\\ on one or more rows of input, with optional
first- and second-order derivatives. There is no time integration; the
principal use cases are observation maps for likelihood-based inference
and reparametrizations for
[`solveODE()`](https://dmodverse.github.io/cppDE/reference/solveODE.md).
`derivMode` selects which derivative products are built. Every returned
function runs compiled code: an object built with `compile = FALSE` is
evaluable only after
[`compile()`](https://dmodverse.github.io/cppDE/reference/compile.md).
See [Derivative directions of
cppFUN()](https://dmodverse.github.io/cppDE/articles/methods-cppfun.html)
for the computational paths and the pass-through of unmodeled inputs.

## Usage

``` r
cppFUN(
  eqns,
  variables = getSymbols(eqns, omit = parameters),
  parameters = NULL,
  fixed = NULL,
  modelname = NULL,
  outdir = tempdir(),
  compile = FALSE,
  verbose = FALSE,
  convenient = TRUE,
  deriv = TRUE,
  deriv2 = FALSE,
  derivMode = "forward"
)
```

## Arguments

- eqns:

  Named character vector or list of algebraic expressions, see
  [expressions](https://dmodverse.github.io/cppDE/reference/expressions.md).
  Names define the output variables; defaults to `f1`, `f2`, ... when
  unnamed.

- variables:

  Character vector of variable names supplied per observation. Defaults
  to all symbols in `eqns` not in `parameters`.

- parameters:

  Character vector of parameter names (constant across observations).
  Default `NULL`.

- fixed:

  Optional character vector of symbols excluded from derivative
  computation. Default `NULL`.

- modelname:

  Optional base name for generated C++ symbols and files.

- outdir:

  Directory for generated C++ source files. Default
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html).

- compile:

  Logical. Compile and load the generated C++ code. Default `FALSE`;
  compile several objects together with
  [`compile()`](https://dmodverse.github.io/cppDE/reference/compile.md).

- verbose:

  Logical. Print progress messages. Default `FALSE`.

- convenient:

  Logical, default `TRUE`. Return `func`, `jac`, `hess` and `evaluate`
  as functions of named arguments rather than of `(vars, params)`.
  `evaluateBatch` and `vjp` always take `vars` and `params`.

- deriv:

  Logical. Generate first-order derivative entry points. Default `TRUE`;
  `FALSE` generates `func` alone, whatever `derivMode` says.

- deriv2:

  Logical. Generate Hessian entry points; implies `deriv = TRUE`.
  Default `FALSE`.

- derivMode:

  Which derivative products to build, any of `"forward"` (default),
  `"reverse"` and `"forward-reverse"`.

  - `"forward"`: forward-mode AD on `cppde::dual`, delivering `jac`,
    `hess`, `evaluate` and `evaluateBatch`.

  - `"reverse"`: the vector-Jacobian product `vjp`.

  - `"forward-reverse"`: `vjp` and the derivative of `vjp` along a
    tangent.

  A mode not named is not generated and costs no compile time.

## Value

A list of functions, each `NULL` when not generated. `vars`, `params`,
the seeds, `attach.input` and `fixed` are described in Details.

- `func(..., attach.input = FALSE, fixed = NULL)` returns the outputs,
  `[n_obs, n_out]`.

- `jac(..., tangentX = NULL, tangentP = NULL, attach.input = FALSE, fixed = NULL)`
  returns their tangent, `[n_obs, n_out, n_theta]`.

- `hess(..., tangentX = NULL, tangentP = NULL, hessianX = NULL, hessianP = NULL, attach.input = FALSE, fixed = NULL)`
  returns their second derivative, `[n_obs, n_out, n_theta, n_theta]`.
  Needs `deriv2 = TRUE`.

- `evaluate(..., tangentX = NULL, tangentP = NULL, hessianX = NULL, hessianP = NULL, deriv2 = FALSE, attach.input = FALSE, fixed = NULL)`
  returns a list with `y`, `tangent` and, with `deriv2 = TRUE`,
  `hessian`, shaped as the results of `func`, `jac` and `hess`.

- `evaluateBatch(sets, cores = 1, deriv2 = FALSE)` runs `evaluate` over
  `sets`, a list of argument lists with elements `vars`, `params` and
  optionally the seeds, `fixed` and `attach.input`, in one call and on
  `cores` threads. It returns the list of results.

- `vjp(vars, params, cotangent, tangentX = NULL, tangentP = NULL, curvature = NULL)`
  contracts the Jacobian with `cotangent`, the cotangent of the outputs,
  `[n_obs, n_out]` or `[n_obs, n_out, n_seed]`. It returns `y`,
  `cotangentX`, `[n_obs, n_vars, n_seed]`, and `cotangentP`,
  `[n_params, n_seed]`, summed over the rows because the parameters are
  shared. Given `tangentX`, `[n_obs, n_vars, n_dir]`, `tangentP`,
  `[n_params, n_dir]`, or `curvature`, `[n_obs, n_out, n_seed, n_dir]`,
  the derivative of the cotangent along the tangent, it also returns
  `curvatureX`, `[n_obs, n_vars, n_seed, n_dir]`, and `curvatureP`,
  `[n_params, n_seed, n_dir]`; this needs `"forward-reverse"`. These
  seeds are positional, without name matching.

`jac`, `hess`, `evaluate` and `evaluateBatch` need `"forward"`, `vjp`
needs `"reverse"` or `"forward-reverse"`.

The list has attributes `equations`, `variables`, `parameters`, `fixed`,
`modelname`, `srcfile` and `derivMode`.

## Details

### Inputs

`n_obs` is the number of rows of input, `n_vars` the number of
variables, `n_params` the number of parameters, `n_out` the number of
outputs and `n_theta` the number of directions \\\theta\\. The variables
are a matrix `vars`, `[n_obs, n_vars]` with columns named by the
variables, the parameters a named numeric vector `params`. With
`convenient = TRUE` the returned `func`, `jac`, `hess` and `evaluate`
take them as named arguments instead, each variable a vector of length
`n_obs` and each parameter a number.

The returned `jac`, `hess` and `evaluate` take the chain rule through
the seeds `tangentX` and `tangentP`, the tangent in \\\theta\\ of the
variables and parameters, and `hessianX` and `hessianP`, their second
derivative:

- `tangentX`: `[n_obs, n_vars, n_theta]`, dimension names
  `list(NULL, variables, theta)`.

- `tangentP`: `[n_params, n_theta]`, row names the parameters.

- `hessianX`: `[n_obs, n_vars, n_theta, n_theta]`.

- `hessianP`: `[n_params, n_theta, n_theta]`.

Rows missing from a seed are zero. A seed not given is taken from
attribute `"deriv"` (tangent) or `"deriv2"` (Hessian) of `vars` or
`params`. Without any seed, the derivatives are taken in the variables
and in the parameters not in `fixed`.

[`solveODE()`](https://dmodverse.github.io/cppDE/reference/solveODE.md)
takes the inputs as one vector `parms` and the derivatives as `tangent`,
`hessian`, `cotangent` and `curvature`. Here the inputs are `vars` and
`params`, and each derivative is split into a variable part (suffix `X`)
and a parameter part (suffix `P`).

`attach.input = TRUE` passes inputs that are neither variables nor
parameters through to the result. `fixed` at run time excludes further
parameters from differentiation.

## See also

[expressions](https://dmodverse.github.io/cppDE/reference/expressions.md)
for the syntax of equations;
[`compile()`](https://dmodverse.github.io/cppDE/reference/compile.md)
for compilation;
[`cppODE()`](https://dmodverse.github.io/cppDE/reference/cppODE.md) and
[`cvode()`](https://dmodverse.github.io/cppDE/reference/cvode.md) for
ODE integration;
[`ptc()`](https://dmodverse.github.io/cppDE/reference/ptc.md) for steady
states; [Derivative directions of
cppFUN()](https://dmodverse.github.io/cppDE/articles/methods-cppfun.html)
on the cppDE site.

## Examples

``` r
# \donttest{
if (codegenAvailable()) {
  ## Two observables of two states and a scale, compiled with their Jacobian.
  ## The functions take the inputs as named arguments.
  g <- cppFUN(c(y = "s*(x1 + x2)", z = "log(x1)"), variables = c("x1", "x2"),
              parameters = "s", compile = TRUE)
  g$func(x1 = 1, x2 = 2, s = 0.5)
  g$jac(x1 = 1, x2 = 2, s = 0.5)[1, , ]
}
# }
```
