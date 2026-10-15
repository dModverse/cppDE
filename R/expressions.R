#' Expressions in Model Equations
#'
#' Right-hand sides, root functions and event values of [cppODE()] and
#' [cvode()] models and the equations of [cppFUN()] are character strings in R
#' syntax. Both ODE backends read them alike.
#'
#' @section Arithmetic:
#' `+`, `-`, `*` and `/`, and powers written as `^` or `**`, as in
#' `"k1*A^2/(K + A)"`. Operators R does not have or reads otherwise are an
#' error: `%%` (write `a - floor(a/b)*b`), `%/%` (write `floor(a/b)`), `%`,
#' `//`, `~`, `@`, and `<-` without a space, which R reads as an assignment.
#'
#' @section Comparisons and logic:
#' The comparisons `<`, `<=`, `>`, `>=`, `==` and `!=` and the logical
#' operators `&&`, `||` and `!` give 1 where they hold and 0 elsewhere. `&`
#' and `|` read as `&&` and `||`, which they equal on numbers. They can stand
#' anywhere a number can:
#'
#' * as a factor, `"k*(A > A0)"` or `"k*(time >= t1 && time < t2)"`;
#' * as the condition of `piecewise()`, `"piecewise(kon, A > A0, koff)"`, or
#'   of `ifelse()`, `"ifelse(A > A0, kon, koff)"`.
#'
#' A comparison does not chain: `"t1 < time < t2"` is an error, write
#' `"t1 < time && time < t2"`.
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
#' | Kinks and steps | `abs`, `min`, `max`, `floor`, `ceiling`, `round`, `trunc` |
#' | Switches | `piecewise(v1, c1, v2, c2, ..., otherwise)`, `ifelse(c, v1, v2)`, `Heaviside(x)`, `sign(x)` |
#'
#' Every function works with every backend and derivative mode. `floor()`,
#' `ceiling()`, `round()` and `trunc()` have the derivative 0; `round()`
#' rounds half to even, as R does. At a kink or step the derivative is that of
#' the branch the value takes, forward and reverse alike: `sign()` and
#' `Heaviside()` have the derivative 0, and a tie of `max()` takes its first
#' argument, of `min()` its second. In the right-hand side of an ODE model they
#' take parameters and forcings only: of a state or of the time they jump
#' where no root is located, which is an error. Write such a step with
#' comparisons, `"piecewise(0, x < 1, 1, x < 2, 2)"`.
#'
#' @section Switches:
#' `piecewise()` takes its arguments in SBML order: the first value whose
#' condition holds, else `otherwise`. A comparison `<`, `<=`, `>` or `>=` on a
#' state or on the time switches the right-hand side, inside `piecewise()` or
#' as a factor, and so do `Heaviside()` and `sign()` of such an amount. The
#' solver locates every such switch as a root and continues on the other
#' branch, so a pulse in time is not stepped over, and the sensitivities take
#' the jump, also where a parameter moves the switch. Where both branches push
#' the solution onto the surface it slides along it, so `"v*(x > 0)"` is the
#' Km to 0 limit of `"v*x/(Km + x)"`, also with an input. `==` and `!=` are
#' evaluated as written. Details are in [cppODE()].
#'
#' @section Symbols:
#' `time` is the independent variable, and `pi` is the number, as in R; `pi`
#' cannot name a state, a parameter, a forcing or a variable. Every other name
#' is a state, a parameter, a forcing or a variable of [cppFUN()], `E`
#' included: write `exp(1)` for Euler's number. Python keywords such as
#' `lambda`, `in` or `if`, and `True`, `False` and `None`, cannot be names.
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
