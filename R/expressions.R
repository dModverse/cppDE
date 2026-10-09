#' Expressions in Model Equations
#'
#' Right-hand sides, root functions and event values of [cppODE()] and
#' [cvode()] models and the equations of [cppFUN()] are character strings in R
#' syntax. Both ODE backends read them alike.
#'
#' @section Arithmetic:
#' `+`, `-`, `*` and `/`, and powers written as `^` or `**`, as in
#' `"k1*A^2/(K + A)"`.
#'
#' @section Comparisons and logic:
#' The comparisons `<`, `<=`, `>`, `>=`, `==` and `!=` and the logical
#' operators `&&`, `||` and `!` give 1 where they hold and 0 elsewhere. They
#' can stand anywhere a number can:
#'
#' * as a factor, `"k*(A > A0)"` or `"k*(time >= t1 && time < t2)"`;
#' * as the condition of `piecewise()`, `"piecewise(kon, A > A0, koff)"`.
#'
#' As in R, `!` negates everything up to the next `&&`, `||`, comma or closing
#' bracket: `"!A > A0"` reads `!(A > A0)`. Write the negation in brackets where
#' it is a factor, `"k*(!(A > A0))"`.
#'
#' @section Functions:
#' | Group | Functions |
#' |---|---|
#' | Exponential and logarithm | `exp`, `exp2`, `exp10`, `log`, `ln`, `log2`, `log10`, `log(x, b)` |
#' | Roots and powers | `sqrt`, `cbrt`, `root(x, n)`, `pow(x, y)` |
#' | Trigonometric | `sin`, `cos`, `tan`, `cot`, `sec`, `csc`, `asin`, `acos`, `atan`, `acot`, `asec`, `acsc`, `atan2(y, x)` |
#' | Hyperbolic | `sinh`, `cosh`, `tanh`, `coth`, `sech`, `csch`, `asinh`, `acosh`, `atanh`, `acoth`, `asech`, `acsch` |
#' | Error function | `erf`, `erfc` |
#' | Kinks and steps | `abs`, `min`, `max`, `floor`, `ceiling`, `round` |
#' | Switches | `piecewise(v1, c1, v2, c2, ..., otherwise)`, `Heaviside(x)`, `sign(x)` |
#'
#' Every function works with every backend and derivative mode. `floor()`,
#' `ceiling()` and `round()` have the derivative 0.
#'
#' @section Switches:
#' `piecewise()` takes its arguments in SBML order: the first value whose
#' condition holds, else `otherwise`. A comparison `<`, `<=`, `>` or `>=` on a
#' state or on the time switches the right-hand side, inside `piecewise()` or
#' as a factor, and so do `Heaviside()` and `sign()` of such an amount. The
#' solver locates every such switch as a root and continues on the other
#' branch, so a pulse in time is not stepped over, and the sensitivities take
#' the jump, also where a parameter moves the switch. `==` and `!=` are
#' evaluated as written. Details are in [cppODE()].
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
#' @example inst/examples/expressions.R
NULL
