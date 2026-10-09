#' Expressions in Model Equations
#'
#' Right-hand sides, root functions and event values of [cppODE()] and
#' [cvode()] models and the equations of [cppFUN()] are character strings in R
#' syntax. They are read with SymPy and translated to C++.
#'
#' @section Operators:
#' `+`, `-`, `*`, `/`, powers as `^` or `**`, the comparisons `<`, `<=`, `>`,
#' `>=`, `==` and `!=`, and the logical operators `&&`, `||` and `!`.
#'
#' @section Functions:
#' | Group | Functions |
#' |---|---|
#' | Exponential and logarithm | `exp`, `exp2`, `exp10`, `log`, `ln`, `log2`, `log10`, `log(x, b)` |
#' | Roots and powers | `sqrt`, `cbrt`, `root(x, n)`, `pow(x, y)` |
#' | Trigonometric | `sin`, `cos`, `tan`, `cot`, `sec`, `csc`, `asin`, `acos`, `atan`, `acot`, `asec`, `acsc`, `atan2(y, x)` |
#' | Hyperbolic | `sinh`, `cosh`, `tanh`, `coth`, `sech`, `csch`, `asinh`, `acosh`, `atanh`, `acoth`, `asech`, `acsch` |
#' | Error function | `erf`, `erfc` |
#' | Piecewise | `abs`, `sign`, `Heaviside`, `min`, `max`, `floor`, `ceiling`, `round` |
#' | Switches | `piecewise(v1, c1, v2, c2, ..., otherwise)` |
#'
#' Every function works with every backend and derivative mode. `piecewise()`
#' takes its arguments in SBML order: the first value whose condition holds,
#' else `otherwise`.
#'
#' @section Switches:
#' A comparison `<`, `<=`, `>`, `>=`, also inside `piecewise()` or a logical
#' operator, and `Heaviside()` and `sign()` switch the right-hand side. On time
#' and parameters, affine in time, the switching time is known and the solver
#' stops in front of it; on a state, the switch is located as a root during
#' the solve. Sensitivities take the jump at both kinds, see Details of
#' [cppODE()]. `==` and `!=` are evaluated as written.
#'
#' @section Symbols:
#' `time` is the independent variable. Every other name is a state, a
#' parameter, a forcing or a variable of [cppFUN()], `pi` and `E` included:
#' write `exp(1)` for Euler's number. Python keywords such as `lambda`, `in`
#' or `if`, and `True`, `False` and `None`, cannot be names.
#'
#' @section Further functions:
#' A function not listed here can be requested from the maintainer,
#' `maintainer("cppDE")`, or as an issue at
#' <https://github.com/dModverse/cppDE/issues>.
#'
#' @name expressions
#' @seealso [cppODE()], [cvode()], [cppFUN()]
NULL
