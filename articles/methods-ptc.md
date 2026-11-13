# Steady states by pseudo-transient continuation

[`ptc()`](https://dmodverse.github.io/cppDE/reference/ptc.md) solves for
a steady state by pseudo-transient continuation. dMod2 calls the same
solver from C++ for `Pimpl()` through `cppde_ptc_v1`.
[Notation](https://dmodverse.github.io/cppDE/articles/methods-overview.html#notation);
\\x\\ are the solved variables of a
[`cppFUN()`](https://dmodverse.github.io/cppDE/reference/cppFUN.md)
object, \\f(x, p)\\ the paired outputs.

## Steady-state problem

[`ptc()`](https://dmodverse.github.io/cppDE/reference/ptc.md) solves

\\G(x) = \begin{pmatrix} f(x, p) \\ C x - b \end{pmatrix} = 0 ,
\tag{1}\\

with \\C\\ the constraint matrix and \\b\\ the argument `total`. For a
reaction network the rows of \\C\\ are conserved quantities, and (1)
selects one root per value of \\b\\. With `flow = TRUE` (default), \\f\\
is the right-hand side of \\\dot x = f(x, p)\\ and the solver seeks a
stable steady state, not any root.

## Continuation step

With \\M\\ the identity on the rows of \\f\\ and zero on the rows of
\\C\\, every step solves

\\\left(\frac{M}{\Delta t} - \frac{\partial G}{\partial x}\right) d =
G(x) \tag{2}\\

for the step \\d\\ [\[1\]](#ref-kelley1998convergence),
[\[2\]](#ref-coffey2003pseudotransient). For small \\\Delta t\\, (2) is
an implicit Euler step of \\\dot x = f\\ on \\C x = b\\; for \\\Delta t
\to \infty\\ it is Newton’s method on (1). Far from the root, where
Newton’s method could jump to an unstable root or leave the positive
orthant, the iteration follows the flow; near it, convergence is
quadratic. QR solves (2) in the least-squares sense, rows scaled by the
residual scale of (4), columns by \\\|x_i\| + \mathrm{abstol}\\.

- `positive = TRUE` (default): the update \\x_i \exp(d_i / x_i)\\, \\x +
  d\\ to first order, exponent capped at \\\pm\\`zmax`, projected back
  onto \\C x = b\\ by \\x \exp(C^\mathsf{T} \mu)\\, \\\mu\\ from a
  damped Newton iteration.
- `flow = FALSE`: the step follows the Newton flow, \\d = -(\partial
  G/\partial x)^{+} G \\ \Delta t / (1 + \Delta t)\\, to any regular
  root.

## Pseudo-time step control

The first step is \\\Delta t_0 = \mathrm{dtInit} / \max_i \|J\_{ii}\|\\.
The local error of an implicit Euler step is estimated as

\\e = \left(\frac{M}{\Delta t} - \frac{\partial G}{\partial
x}\right)^{-1} \begin{pmatrix} \tfrac12 \bigl(f(x + d, p) - f(x,
p)\bigr) \\ 0 \end{pmatrix} , \tag{3}\\

in the root mean square of \\e_i / (\mathrm{flowTol}\\\|x_i\| +
\mathrm{abstol})\\. A step with error above 1 is repeated at smaller
\\\Delta t\\. After an accepted step, \\\Delta t\\ changes by the
residual ratio (2 to 1000) if the residual has halved, otherwise by the
error-estimate factor (0.2 to 10). In the Newton range, \\\Delta t \ge
10^8 / \max_i \|J\_{ii}\|\\, a step halving the residual is accepted
regardless of the error. In the flow range the error decides alone,
since such a step can still leave the basin of the stable root.

## Convergence test

Row \\i\\ of (1) has converged when

\\\|G_i\| \le \mathrm{reltol} \sum_k \left\|\frac{\partial G_i}{\partial
x_k}\right\| \|x_k\| + \mathrm{abstol} \max_k \left\|\frac{\partial
G_i}{\partial x_k}\right\| , \tag{4}\\

the residual against the turnover of the row, the sum of the rates that
cancel in it, so the test scales with the rates of the row. In the
Newton range, a step below \\\mathrm{reltol}\\\|x_i\| +
\mathrm{abstol}\\ in every variable also ends the iteration. Contracting
Newton steps then bring every variable to within
\\\mathrm{reltol}\\(\|x_i\| + \mathrm{abstol})\\ of the root. Equal
start, parameters and controls return the same root to the last digit.

## Stability test

With \\N\\ an orthonormal kernel basis of \\C\\ from the singular value
decomposition, the flow on \\C x = b\\ near the root has the Jacobian

\\B = N^\mathsf{T} J N , \tag{5}\\

and the root is stable if every eigenvalue of \\B\\ has negative real
part. With `stability = TRUE` (default, `flow = TRUE` only),
[`ptc()`](https://dmodverse.github.io/cppDE/reference/ptc.md) accepts
the root while the leading eigenvalue satisfies \\\mathrm{Re}\\\lambda
\le \max(10^{-8} \max_k \|\lambda_k\|, \mathrm{stabilityTol})\\. A small
positive `stabilityTol` accepts a root that perturbations leave slower
than a time scale of interest, as near a Hopf bifurcation.

From an unstable root \\x^\*\\ the search restarts at \\x^\* \pm \eta
v\\, \\v\\ the real part of the leading eigenvector mapped back by
\\N\\, first on the side the iterates came from while they followed the
flow (\\\Delta t\\ \mathrm{Re}\\\lambda \< 1\\). \\\eta\\ moves every
variable by at most a tenth of the variable value. Near \\x^\*\\,
\\\Delta t \le 1 / (2\\\mathrm{Re}\\\lambda)\\, so the implicit Euler
step leaves along the unstable direction, where Newton’s method would
return. The new root is tested again; if neither side leads to a stable
root, [`ptc()`](https://dmodverse.github.io/cppDE/reference/ptc.md)
reports the unstable one, leading eigenvalue in `message`.

## References

\[1\]

C. T. Kelley and D. E. Keyes, “Convergence analysis of pseudo-transient
continuation,” *SIAM Journal on Numerical Analysis*, vol. 35, no. 2, pp.
508–523, 1998.

\[2\]

T. S. Coffey, C. T. Kelley, and D. E. Keyes, “Pseudotransient
continuation and differential-algebraic equations,” *SIAM Journal on
Scientific Computing*, vol. 25, no. 2, pp. 553–569, 2003.
