# Forward sensitivities

[`cppODE()`](https://dmodverse.github.io/cppDE/reference/cppODE.md)
computes the sensitivities \\S = \partial x/\partial\theta\\ of the
trajectory by dual numbers, under error control. Symbols are those of
the [notation
table](https://dmodverse.github.io/cppDE/articles/methods-overview.html#notation).
R interface: `deriv` and `derivMode` of
[`cppODE()`](https://dmodverse.github.io/cppDE/reference/cppODE.md) and
`tangent` and `sensErrCon` of
[`solveODE()`](https://dmodverse.github.io/cppDE/reference/solveODE.md).

## Dual-number arithmetic

Sensitivities come from forward-mode automatic differentiation
[\[1\]](#ref-griewank2008evaluating), not finite differences. Arithmetic
and elementary functions are overloaded for the *dual number* \\ \tilde
x = x + x'\\\varepsilon, \qquad \varepsilon^2 = 0,\quad \varepsilon \ne
0, \\ with primal value \\x \in \mathbb{R}\\ and directional derivative
\\x' \in \mathbb{R}\\, so the C++ source of \\f\\ compiles unchanged for
primal and sensitivities. As a first-order Taylor expansion in
\\\varepsilon\\, with \\\tilde y = y + y'\varepsilon\\, \\ \tilde x +
\tilde y = (x+y) + (x'+y')\varepsilon, \qquad \tilde x \cdot \tilde y =
xy + (xy' + x'y)\varepsilon, \\ and for any smooth scalar function
\\g\\, \\ g(\tilde x) = g(x) + g'(x)\\x'\\\varepsilon. \\ By
composition, \\f\\ called with dual inputs returns \\f\\ with the
directional derivative.

## Tangent storage

A dual holds a value and a run-time number of tangent slots. A
second-order dual nests a dual in a dual and propagates only the lower
triangle of the Hessian. The state vectors kept across steps are
allocated once per solve, so memory does not grow with the number of
steps.

## Sensitivity equations by dual evaluation

For a right-hand side depending on parameters \\p \in
\mathbb{R}^{n_p}\\, \\ \dot x = f(x, p, t), \qquad x(t_0) = x_0(p), \\
the *sensitivity matrix* is \\ S_p = \frac{\partial x}{\partial p} \in
\mathbb{R}^{n_x \times n_p}, \qquad (S_p)\_{ij} = \frac{\partial
x_i}{\partial p_j}. \\ Seeding state \\i\\ with row \\i\\ of \\S_p\\
gives \\ \tilde x_i = x_i + (S_p)\_{ij}\\\varepsilon_j, \\ with dual
unit \\\varepsilon_j\\ for direction \\p_j\\ and summation over repeated
indices. The \\\varepsilon_j\\-coefficient of the \\i\\-th output of
\\f\\ evaluated over the dual algebra is \\ \bigl\[f(\tilde x, p,
t)\bigr\]\_{i,\\\varepsilon_j} = \frac{\partial f_i}{\partial
x_k}\\(S_p)\_{kj} + \frac{\partial f_i}{\partial p_j}. \\
Differentiating the primal part \\\dot x_i = f_i(x,p,t)\\ in \\p_j\\
identifies this coefficient with \\\dot{(S_p)}\_{ij}\\: \\ \dot S_p =
\frac{\partial f}{\partial x}\\S_p + \frac{\partial f}{\partial p},
\qquad S_p(t_0) = \frac{\partial x_0}{\partial p}. \tag{1} \\ One dual
evaluation of \\f\\ gives [(1)](#eq-sensitivity-p) for all columns, at a
cost proportional to \\n_p\\, without forming \\\partial f/\partial x\\.

## Seeding of the tangent

For a reparametrization \\(x_0, p) = \Phi(\theta)\\, \\\theta \in
\mathbb{R}^M\\, \\ S\_\theta = \frac{\partial x}{\partial \theta} \in
\mathbb{R}^{n_x \times M}, \\ and by the chain rule
[(1)](#eq-sensitivity-p) becomes \\ \dot S\_\theta = \frac{\partial
f}{\partial x}\\S\_\theta + \frac{\partial f}{\partial p}
\\\frac{\partial p}{\partial \theta}, \qquad S\_\theta(t_0) =
\frac{\partial x_0}{\partial \theta}. \tag{2} \\ The *tangent*
\\\partial (x_0, p) / \partial\theta\\ is the `tangent` argument of
[`solveODE()`](https://dmodverse.github.io/cppDE/reference/solveODE.md),
fixed during integration. The identity tangent gives \\\theta = (x_0,
p)\\. Each direction costs one state vector.

## Error norm with sensitivities

Every step produces \\1+M\\ error vectors of length \\n_x\\: \\e^{(x)}\\
for the state, with norm [Stepping
(8)](https://dmodverse.github.io/cppDE/articles/methods-stepping.html#eq-wrms),
and \\e_S^{(k)}\\ for each column \\\partial x/\partial \theta_k\\, with
\\S\_\theta\\ in place of \\x\\ in error vector and weight, \\ \\
e_S^{(k)} \\\_{\mathrm{WRMS}} = \sqrt{\frac{1}{n_x} \sum\_{i=1}^{n_x}
\left( \frac{(e_S)\_{ik}}{w^{(k)}\_{i}} \right)^{\\2}}, \qquad
w^{(k)}\_{i} = \mathtt{abstol} +
\mathtt{reltol}\\\|(S\_\theta)\_{ik}^{\mathrm{ref}}\|, \tag{3} \\ with
\\(S\_\theta)\_{ik}^{\mathrm{ref}}\\ taken like \\x^{\mathrm{ref}}\\:
the accepted value for multistep methods, the maximum over both step
ends for single-step methods. As in CVODES, the controller uses the
largest norm, \\ \\ e\_{n+1} \\\_{\mathrm{WRMS}} = \max\left\\ \\
e^{(x)} \\\_{\mathrm{WRMS}},\\ \max\_{1 \le k \le M} \\ e_S^{(k)}
\\\_{\mathrm{WRMS}} \right\\. \tag{4} \\ A flat WRMS over all
\\n_x(1+M)\\ components would let small errors in some columns mask a
large one in another; the maximum holds each direction to the tolerance.
The BDF/NDF corrector uses the same maximum, so an unconverged
sensitivity direction forces further iterations or a smaller step. With
`sensErrCon = FALSE`, error test, order selection, corrector test and
initial-step estimate use the state only, and steps and iteration counts
are those of a value run.

## Corrector iterations

In `bdf` and `adams` the iterations [Stepping
(3)](https://dmodverse.github.io/cppDE/articles/methods-stepping.html#eq-newton)
run on dual numbers: each tangent takes the iterations of the values,
against the same \\W\\ and without the derivative of \\W\\, which like
the step size is a solver decision. With `sensErrCon = FALSE` the
sensitivities are the derivative of the computed solution on the grid of
a value run, as the reverse mode returns ([Reverse
mode](https://dmodverse.github.io/cppDE/articles/methods-reverse.html#discrete-and-continuous-adjoints)).

## References

\[1\]

A. Griewank and A. Walther, *Evaluating derivatives: Principles and
techniques of algorithmic differentiation*, 2nd ed. in Other titles in
applied mathematics, no. 105. SIAM, 2008. doi:
[10.1137/1.9780898717761](https://doi.org/10.1137/1.9780898717761).
