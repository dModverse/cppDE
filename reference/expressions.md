# Expressions in Model Equations

Right-hand sides, root functions and event values of
[`cppODE()`](https://dmodverse.github.io/cppDE/reference/cppODE.md) and
[`cvode()`](https://dmodverse.github.io/cppDE/reference/cvode.md) models
and the equations of
[`cppFUN()`](https://dmodverse.github.io/cppDE/reference/cppFUN.md) are
character strings in R syntax. Both ODE backends read them alike.

## Arithmetic

`+`, `-`, `*` and `/`, and powers written as `^` or `**`, as in
`"k1*A^2/(K + A)"`. Operators R does not have or reads otherwise are an
error: `%%` (write `a - floor(a/b)*b`), `%/%` (write `floor(a/b)`), `%`,
`//`, `~`, `@`, and `<-` without a space, which R reads as an
assignment.

## Comparisons and logic

The comparisons `<`, `<=`, `>`, `>=`, `==` and `!=` and the logical
operators `&&`, `||` and `!` give 1 where they hold and 0 elsewhere. `&`
and `|` read as `&&` and `||`, which they equal on numbers. They can
stand anywhere a number can:

- as a factor, `"k*(A > A0)"` or `"k*(time >= t1 && time < t2)"`;

- as the condition of `piecewise()`, `"piecewise(kon, A > A0, koff)"`,
  or of [`ifelse()`](https://rdrr.io/r/base/ifelse.html),
  `"ifelse(A > A0, kon, koff)"`.

A comparison does not chain: `"t1 < time < t2"` is an error, write
`"t1 < time && time < t2"`.

As in R, `!` negates everything up to the next `&&`, `||`, comma or
closing bracket: `"!A > A0"` reads `!(A > A0)`. Write the negation in
brackets where it is a factor, `"k*(!(A > A0))"`.

## Functions

|  |  |
|----|----|
| Group | Functions |
| Exponential and logarithm | `exp`, `exp2`, `exp10`, `log`, `ln`, `log2`, `log10`, `log(x, b)` |
| Roots and powers | `sqrt`, `cbrt`, `root(x, n)`, `pow(x, y)` |
| Trigonometric | `sin`, `cos`, `tan`, `cot`, `sec`, `csc`, `asin`, `acos`, `atan`, `acot`, `asec`, `acsc`, `atan2(y, x)` |
| Hyperbolic | `sinh`, `cosh`, `tanh`, `coth`, `sech`, `csch`, `asinh`, `acosh`, `atanh`, `acoth`, `asech`, `acsch` |
| Error function | `erf`, `erfc` |
| Kinks and steps | `abs`, `min`, `max`, `floor`, `ceiling`, `round`, `trunc` |
| Switches | `piecewise(v1, c1, v2, c2, ..., otherwise)`, `ifelse(c, v1, v2)`, `Heaviside(x)`, `sign(x)` |

Every function works with every backend and derivative mode.
[`floor()`](https://rdrr.io/r/base/Round.html),
[`ceiling()`](https://rdrr.io/r/base/Round.html),
[`round()`](https://rdrr.io/r/base/Round.html) and
[`trunc()`](https://rdrr.io/r/base/Round.html) have the derivative 0;
[`round()`](https://rdrr.io/r/base/Round.html) rounds half to even, as R
does. At a kink or step the derivative is that of the branch the value
takes, forward and reverse alike:
[`sign()`](https://rdrr.io/r/base/sign.html) and `Heaviside()` have the
derivative 0, and a tie of
[`max()`](https://rdrr.io/r/base/Extremes.html) takes the derivative of
the first argument, of [`min()`](https://rdrr.io/r/base/Extremes.html)
that of the second. In the right-hand side of an ODE model they take
parameters and forcings only: of a state or of the time they jump where
no root is located, which is an error. Write such a step with
comparisons, `"piecewise(0, x < 1, 1, x < 2, 2)"`.

## Switches

`piecewise()` takes the arguments in SBML order: the first value whose
condition holds, else `otherwise`. A comparison `<`, `<=`, `>` or `>=`
on a state or on the time switches the right-hand side, inside
`piecewise()` or as a factor, and so do `Heaviside()` and
[`sign()`](https://rdrr.io/r/base/sign.html) of such an amount. The
solver locates every such switch as a root and continues on the other
branch, so a pulse in time is not stepped over, and the sensitivities
take the jump, also where a parameter moves the switch. Where both
branches push the solution onto the surface it slides along it, so
`"v*(x > 0)"` is the Km to 0 limit of `"v*x/(Km + x)"`, also with an
input. `==` and `!=` are evaluated as written. Details are in
[`cppODE()`](https://dmodverse.github.io/cppDE/reference/cppODE.md).

## Symbols

`time` is the independent variable, and `pi` is the number, as in R;
`pi` cannot name a state, a parameter, a forcing or a variable. Every
other name is a state, a parameter, a forcing or a variable of
[`cppFUN()`](https://dmodverse.github.io/cppDE/reference/cppFUN.md), `E`
included: write `exp(1)` for Euler's number. Python keywords such as
`lambda`, `in` or `if`, and `True`, `False` and `None`, cannot be names.

## Further functions

A function not listed here can be requested from the maintainer,
`maintainer("cppDE")`, or as an issue at
<https://github.com/dModverse/cppDE/issues>.

## See also

[`cppODE()`](https://dmodverse.github.io/cppDE/reference/cppODE.md),
[`cvode()`](https://dmodverse.github.io/cppDE/reference/cvode.md),
[`cppFUN()`](https://dmodverse.github.io/cppDE/reference/cppFUN.md)

## Examples

``` r
# \donttest{
if (codegenAvailable()) {
  ## Comparisons and logic give 1 or 0: as a factor or as a condition.
  ##   u: a pulse between t1 and t2, `&&` joins two comparisons on time
  ##   x: decays fast below xc and slowly above, unless mode is 0
  ##   y: produced once x exceeds xc, with powers as `^` and `**`
  eqns <- c(
    u = "a*(time >= t1 && time < t2) - u/tau",
    x = "k*u - piecewise(kfast, mode == 0 || !(x > xc), kslow)*x",
    y = "kp*Heaviside(x - xc)*(1 - y^2) - y**2/(K + y)"
  )
  parms <- c(u = 0, x = 0, y = 0, a = 2, t1 = 1, t2 = 3, tau = 1, k = 1,
             kslow = 0.1, kfast = 1, xc = 0.5, mode = 1, kp = 1, K = 0.5)
  times <- 0:8

  ## The same equations on both backends, compiled together
  models <- list(cppODE = cppODE(eqns, deriv = FALSE, compile = FALSE))
  if (codegenAvailable("cvode"))
    models$cvode <- cvode(eqns, compile = FALSE)
  do.call(compile, c(unname(models), list(output = "expressions_example")))

  out <- lapply(models, function(m) solveODE(m, times, parms)$variable)
  cbind(time = times, out$cppODE)
  if (!is.null(out$cvode)) max(abs(out$cppODE - out$cvode))
}
# }
```
