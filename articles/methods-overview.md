# Overview and notation

## Entry points

| entry point | problem | numerical core | article |
|:---|:---|:---|:---|
| [`cppODE()`](https://dmodverse.github.io/cppDE/reference/cppODE.md) | ODE \\\dot x = f(x, p, t)\\ | in-tree integrator, `inst/include/cppde/` | [Stepping](https://dmodverse.github.io/cppDE/articles/methods-stepping.md) to [Second derivatives](https://dmodverse.github.io/cppDE/articles/methods-forward-reverse.md) |
| [`cvode()`](https://dmodverse.github.io/cppDE/reference/cvode.md) | ODE \\\dot x = f(x, p, t)\\ | SUNDIALS CVODE/CVODES | [cvode() and the CVODES adjoint](https://dmodverse.github.io/cppDE/articles/methods-cvode.md) |
| [`cppFUN()`](https://dmodverse.github.io/cppDE/reference/cppFUN.md) | algebraic \\y = g(x, p)\\ | none | [Derivative directions of cppFUN()](https://dmodverse.github.io/cppDE/articles/methods-cppfun.md) |
| [`ptc()`](https://dmodverse.github.io/cppDE/reference/ptc.md) | steady state \\f(x, p) = 0\\, \\C x = b\\ | pseudo-transient continuation on a [`cppFUN()`](https://dmodverse.github.io/cppDE/reference/cppFUN.md) object | [Steady states by pseudo-transient continuation](https://dmodverse.github.io/cppDE/articles/methods-ptc.md) |

Entry points. The first three share the [code
generator](https://dmodverse.github.io/cppDE/articles/methods-codegen.md).
{.table}

[`cppODE()`](https://dmodverse.github.io/cppDE/reference/cppODE.md) and
[`cvode()`](https://dmodverse.github.io/cppDE/reference/cvode.md) take
the same model arguments (`rhs`, `events`, `rootfunc`, `fixed`,
`forcings`) and return the same object for
[`solveODE()`](https://dmodverse.github.io/cppDE/reference/solveODE.md)
([comparison](https://dmodverse.github.io/cppDE/articles/methods-cvode.html#comparison)).
[`cppFUN()`](https://dmodverse.github.io/cppDE/reference/cppFUN.md)
derivatives are `"forward"` (second order with `deriv2 = TRUE`),
`"reverse"` or `"forward-reverse"`.

## Order of the articles

1.  [Code
    generation](https://dmodverse.github.io/cppDE/articles/methods-codegen.md):
    expression graph and emitted C++ of all entry points.
2.  [Stepping and error
    control](https://dmodverse.github.io/cppDE/articles/methods-stepping.md):
    the four methods of
    [`cppODE()`](https://dmodverse.github.io/cppDE/reference/cppODE.md),
    error control, dense output.
3.  [Forward
    sensitivities](https://dmodverse.github.io/cppDE/articles/methods-forward.md):
    dual numbers.
4.  [Events, switches and
    saltation](https://dmodverse.github.io/cppDE/articles/methods-events.md):
    sensitivities across resets, root events and switches.
5.  [Reverse mode: the discrete
    adjoint](https://dmodverse.github.io/cppDE/articles/methods-reverse.md):
    gradients at a cost independent of the parameter count.
6.  [Events in reverse
    mode](https://dmodverse.github.io/cppDE/articles/methods-reverse-events.md):
    adjoint of a jump, to second order.
7.  [Grid and refinement of the
    sweep](https://dmodverse.github.io/cppDE/articles/methods-reverse-checks.md):
    the adjoint grid and `refine`.
8.  [Second derivatives: forward over
    reverse](https://dmodverse.github.io/cppDE/articles/methods-forward-reverse.md):
    Hessian-vector products.
9.  [cvode() and the CVODES
    adjoint](https://dmodverse.github.io/cppDE/articles/methods-cvode.md):
    the SUNDIALS backend.
10. [Derivative directions of
    cppFUN()](https://dmodverse.github.io/cppDE/articles/methods-cppfun.md):
    observation maps and reparameterizations.
11. [Steady states by pseudo-transient
    continuation](https://dmodverse.github.io/cppDE/articles/methods-ptc.md):
    [`ptc()`](https://dmodverse.github.io/cppDE/reference/ptc.md) under
    conservation laws, stability test.
12. [Method
    coefficients](https://dmodverse.github.io/cppDE/articles/methods-coefficients.md):
    appendix.

Batches of conditions: [Batches of
conditions](https://dmodverse.github.io/cppDE/articles/guide-batches.md).
[Benchmarks](https://dmodverse.github.io/cppDE/articles/benchmarks.md)
compares the integrators and the derivative modes.

## Notation

| symbol | meaning |
|:---|:---|
| \\t\\, \\t_0\\, \\T\\ | time, initial time, last output time |
| \\x \in \mathbb{R}^{n_x}\\ | state; \\x_0 = x(t_0)\\ |
| \\p \in \mathbb{R}^{n_p}\\ | parameters of the right-hand side |
| \\f(x, p, t)\\ | right-hand side, \\\dot x = f\\ |
| \\J = \partial f/\partial x\\ | state Jacobian |
| \\y = g(x, p)\\ | observables, the prediction compared with data; the algebraic function of [`cppFUN()`](https://dmodverse.github.io/cppDE/reference/cppFUN.md) |
| \\\theta \in \mathbb{R}^{M}\\ | parameters of the derivatives, with \\x_0(\theta)\\, \\p(\theta)\\; \\\theta = (x_0, p)\\, \\M = n_x + n_p\\ under the identity tangent |
| \\\partial(x_0, p)/\partial\theta\\ | tangent of the inputs, argument `tangent` of [`solveODE()`](https://dmodverse.github.io/cppDE/reference/solveODE.md) |
| \\S = \partial x/\partial\theta \in \mathbb{R}^{n_x \times M}\\ | sensitivity matrix |
| \\t_o\\, \\o = 1, \dots, n\_{\mathrm{obs}}\\ | output times, rows of the returned trajectory |
| \\w \in \mathbb{R}^{n\_{\mathrm{obs}} \times n_x}\\ | cotangent of the outputs, row \\w_o\\ at \\t_o\\ |
| \\\lambda \in \mathbb{R}^{n_x}\\ | adjoint state |
| \\\bar\theta\\, \\\bar x_0\\, \\\bar p\\ | cotangent of the inputs, the gradient of \\\sum_o w_o^{\mathsf T} x(t_o)\\ |
| \\\mathcal{J}(\theta)\\ | objective function of a fit, a scalar functional of the outputs; \\\bar\theta = \nabla\_\theta \mathcal{J}\\ for \\w_o = \partial\mathcal{J}/\partial x(t_o)\\ |
| \\t_n\\, \\x_n\\ | step points and computed states |
| \\h_n = t\_{n+1} - t_n\\ | step size |
| \\q\\, \\q\_{\max}\\ | order, maximum order |
| \\e\_{n+1}\\ | local error estimate of the step from \\t_n\\ |
| \\\lVert\cdot\rVert\_{\mathrm{WRMS}}\\ | weighted root-mean-square norm, [Stepping (8)](https://dmodverse.github.io/cppDE/articles/methods-stepping.html#eq-wrms) |
| \\W\\ | iteration matrix of an implicit step, a shift of \\-J\\ |
| \\h_e(x, t, \theta)\\ | root function of an event or switch |
| \\t_e(\theta)\\ | firing time of a time-triggered event, an expression |
| \\t^{\*}\\, \\t\_\*^{(0)}\\ | firing time; the value of \\t^{\*}\\ without derivatives, the grid time |
| \\\delta t^{\*} = t^{\*} - t\_\*^{(0)}\\ | firing-time residual: zero in value, holds the derivatives of \\t^{\*}\\ |
| \\\dot h_e = \partial h_e/\partial t + \nabla_x h_e \cdot f\\ | rate of the root function along the trajectory |
| \\R(x, t)\\ | reset map of an event |
| \\x^{-}\\, \\x^{+}\\; \\f^{-}\\, \\f^{+}\\ | state and right-hand side before and after an event |

Symbols shared by the methods articles; others are defined where they
first appear. {.table}
