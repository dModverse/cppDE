# Stepping and error control

Arguments: `method`, `useNDF` and `sparse` of
[`cppODE()`](https://dmodverse.github.io/cppDE/reference/cppODE.md);
`abstol`, `reltol`, `maxOrder`, `stald` and `sensErrCon` of
[`solveODE()`](https://dmodverse.github.io/cppDE/reference/solveODE.md).
Symbols follow the [notation
table](https://dmodverse.github.io/cppDE/articles/methods-overview.html#notation).

## Initial-value problem

[`cppODE()`](https://dmodverse.github.io/cppDE/reference/cppODE.md)
solves \\ \dot{x} = f(x, p, t), \qquad x(t_0) = x_0, \tag{1} \\ with
\\f\\ explicit in time for forcings. Step \\n\\ advances \\x_n\\ at
\\t_n\\ to \\x\_{n+1}\\ at \\t\_{n+1} = t_n + h_n\\
([coefficients](https://dmodverse.github.io/cppDE/articles/methods-coefficients.md)).

## Stiffness and the choice of method

On \\\dot y = \lambda y\\, \\\lambda \in \mathbb{C}^-\\, a one-step
method multiplies \\y\\ per step by \\R(z)\\, \\z = h \lambda\\. The
method is

- *A-stable* if \\\|R(z)\| \le 1\\ on the left half-plane, so the
  solution decays with the exact one at any step size;
- *A(\\\varphi\\)-stable* if the bound holds on the wedge \\\\ z \in
  \mathbb{C} : \|\arg(-z)\| \le \varphi \\\\ for some \\\varphi \<
  90^{\circ}\\;
- *L-stable* if A-stable with \\R(z) \to 0\\ as \\\operatorname{Re}(z)
  \to -\infty\\, which damps very stiff transients within one step.

Problem [(1)](#eq-ivp) is *stiff* when stability, not accuracy, limits
an explicit step: \\h \lambda_k\\, \\\lambda_k\\ the eigenvalues of
\\J\\ along the trajectory, must stay in the stability domain. Stiffness
depends on spectrum, horizon and tolerance; kinetic networks with rate
constants over orders of magnitude and many discretized parabolic PDEs
are stiff. BDF/NDF and Rosenbrock4, stable on the left half-plane or a
wide wedge of it, suit stiff problems; Adams-Moulton and Tsit5 do not,
at any tolerance. An explicit method fails or takes uneconomic steps on
a stiff problem, an implicit one costs a constant overhead per step on a
non-stiff problem; BDF/NDF is the default, with `sparse = TRUE` for
large sparse Jacobians.

## Multistep methods

A linear multistep method of order \\q\\ combines \\\\x\_{n+1-j},
f\_{n+1-j}\\\_{j=0}^{q}\\ linearly, with few right-hand-side evaluations
per step. cppDE stores BDF and Adams in the Nordsieck form of SUNDIALS
[\[1\]](#ref-hindmarsh2005sundials), the *scaled derivative array* at
\\t_n\\ in place of the past \\q+1\\ states, \\ z =\left(x_n,\\ h_n
\dot{x}\_n,\\ \dfrac{h_n^2}{2!} \ddot{x}\_n, \\\dots,\\
\dfrac{h_n^{q}}{q!} x_n^{(q)} \right). \\ A change of step size or order
is a linear transformation of \\z\\, not a refit.

### Backward Differentiation Formula (`bdf`)

The BDF of order \\q\\ defines \\x\_{n+1}\\ implicitly by \\
\sum\_{j=0}^{q} \alpha_j\\ x\_{n+1-j} \\=\\ h_n\\ f(x\_{n+1}, t\_{n+1}),
\tag{2} \\ exact for polynomials of degree \\\le q\\ on the past \\q+1\\
nodes; A-stable for \\q \le 2\\, A(\\\varphi\\)-stable for \\q \in
\\3,4,5\\\\ with \\\varphi\\ from \\86^\circ\\ (\\q=3\\) to \\51^\circ\\
(\\q=5\\).

The corrector [(2)](#eq-bdf) is solved by simplified Newton iteration,
\\ W \Delta x^{(\nu)} = - s\\ \tilde r \left(x\_{n+1}^{(\nu)}\right),
\qquad x\_{n+1}^{(\nu+1)} = x\_{n+1}^{(\nu)} + \Delta x^{(\nu)}, \tag{3}
\\ with \\\tilde r\\ the residual of [(2)](#eq-bdf) divided by \\h_n
\beta_0\\, \\\beta_0 = 1/\alpha_0\\, and iteration matrix \\ W =
\frac{1}{\gamma_p}I - J\bigl(x_s, t_s\bigr), \tag{4} \\ in the SUNDIALS
form [\[1\]](#ref-hindmarsh2005sundials), \\(I - \gamma_p J)/\gamma_p\\.
\\J\\ is taken at the predicted \\x_s\\, \\t_s\\ of the step that last
evaluated it; \\\gamma_p\\ is \\h\beta_0\\ at the last factorization.
With \\\gamma = h_n\beta_0\\ of the current step, the scale \\s =
(\gamma/\gamma_p)\\ 2/(1 + \gamma/\gamma_p)\\ (BDF) or \\s =
\gamma/\gamma_p\\ (Adams) restores the step of a factorization at
\\\gamma\\, as in CVODE. Factorizing \\W\\ dominates the step cost and
is reused across steps.

The default Numerical Differentiation Formula (NDF)
[\[2\]](#ref-MATLABODESuite) adds to the BDF corrector a multiple
\\\kappa_q H_q\\(x\_{n+1} - x\_{n+1}^{\mathrm{pred}})\\ of the Nordsieck
predictor error, \\H_q = \sum\_{m=1}^{q} 1/m\\ ([BDF and NDF
coefficients](https://dmodverse.github.io/cppDE/articles/methods-coefficients.html#bdf-coeffs)).
The term scales \\1/(h_n\beta_0)\\ in [(4)](#eq-bdf-W) by \\(1 -
\kappa_q)\\ and the leading truncation-error coefficient \\1/(q+1)\\ by
\\1 + \kappa_q (q+1) H_q\\. With \\\kappa_q \< 0\\ the step at equal
accuracy grows by \\26\\\\\\ over BDF for \\q \le 3\\ and \\12\\\\\\ for
\\q = 4\\ [\[2\]](#ref-MATLABODESuite), at a small loss of stability
angle. `useNDF = FALSE` gives BDF (\\\kappa_q \equiv 0\\).

The corrector takes at most three iterations, each one right-hand-side
evaluation and one linear solve, and converges at \\\delta\\\min(1,
\rho)/\tau_q \le 1\\, with \\\delta\\ the WRMS norm of the latest
correction, \\\rho\\ the estimated contraction rate and \\\tau_q\\ an
order-dependent constant. As in CVODE:

- \\W\\ is refactorized when \\h_n\beta_0\\ has moved by more than
  \\30\\\\\\ since the last factorization, after \\20\\ steps, or after
  a failed corrector or error test;
- \\J\\ is re-evaluated at a refactorization once it is \\51\\ steps
  old, and after a failed corrector;
- with the old Jacobian, a sparse \\W\\ changes only on the diagonal.

Like the step size, \\W\\ is a solver decision: forward sensitivities
and the reverse sweep iterate against \\W\\ without differentiating it,
and so differentiate the computed solution ([Forward
sensitivities](https://dmodverse.github.io/cppDE/articles/methods-forward.html#corrector),
[Reverse
mode](https://dmodverse.github.io/cppDE/articles/methods-reverse.html#adjoint-of-the-step-map)).

### Adams-Moulton (`adams`)

Adams-Moulton integrates the interpolant of \\f\\ through
\\\\(t\_{n+1-j}, f\_{n+1-j})\\\_{j=0}^{q}\\ over \\\[t_n, t\_{n+1}\]\\,
\\ x\_{n+1} \\=\\ x_n \\+\\ h_n \sum\_{j=0}^{q} \beta_j\\ f(x\_{n+1-j},
t\_{n+1-j}), \tag{5} \\ implicit through the \\j=0\\ term. As in CVODES
with a linear solver, cppDE applies the BDF Newton iteration on
[(4)](#eq-bdf-W) with the Adams weight \\\beta_0\\: same setup triggers,
Jacobian reuse, convergence test and controller, other coefficients and
\\q\_{\max} = 12\\. The stability domain shrinks with order and for \\q
\ge 3\\ contains no wedge of \\\mathbb{C}^-\\, so Adams suits non-stiff,
smooth problems over long horizons. The Newton corrector keeps a mildly
stiff component from stalling the iteration, but the stability domain
still bounds the step.

## Single-step methods

Single-step methods keep only the current state, reach order \\q \> 1\\
through internal stages, and restart at full order after an event, where
multistep methods rebuild the Nordsieck history from order \\1\\.

### Rosenbrock (`rb4`)

A Rosenbrock method linearizes the implicit Runge-Kutta condition about
\\(x_n, t_n)\\; the \\s\\ stages \\k_i\\ solve \\ W k_i = f\left(x_n +
\sum\_{j\<i} \alpha\_{ij}\\ k_j,\\ t_n + c_i\\ h_n \right) +
\dfrac{1}{h_n} \sum\_{j\<i} \gamma\_{ij}\\ k_j + h_n\\ d_i\\
\dfrac{\partial f}{\partial t}(x_n, t_n), \tag{6} \\ with iteration
matrix \\ W = \frac{1}{\gamma\\ h_n}I - J(x_n, t_n). \tag{7} \\

| symbol | meaning |
|:---|:---|
| \\\alpha\_{ij}\\ | weights of earlier stages \\k_j\\, \\j\<i\\, in stage \\i\\ |
| \\c_i\\ | time abscissa of stage \\i\\ |
| \\\gamma\_{ij}\\ | weights of earlier stages in the linear-solve correction of stage \\i\\ |
| \\\gamma\\ | diagonal coefficient in \\W\\ |
| \\d_i\\ | weights of the explicit time dependence \\\partial f/\partial t\\ |
| \\m_i\\ | weights of the advanced solution \\x\_{n+1}\\ |
| \\\hat{m}\_i\\ | weights of the embedded order-3 estimator for step control |

Rosenbrock method coefficients. {.table}

The advanced solution and embedded estimate are \\ x\_{n+1} = x_n +
\sum\_{i=1}^{s} m_i\\ k_i, \qquad \hat{x}\_{n+1} \\=\\ x_n +
\sum\_{i=1}^{s} \hat{m}\_i\\ k_i. \\ Each stage
[(6)](#eq-rosenbrock-stage) is one linear solve, no Newton iteration; an
accepted step costs one \\J\\, one factorization of \\W\\, six
right-hand-side evaluations and six back-substitutions. The [order-\\4\\
coefficients](https://dmodverse.github.io/cppDE/articles/methods-coefficients.html#ros-coeffs)
of Hairer and Wanner [\[3, Sec. IV.7\]](#ref-hairer2010solving), also
used by the `rosenbrock4` stepper of Boost.Numeric.Odeint, form an
L-stable pair with an embedded order-\\3\\ estimator. Rosenbrock4
competes with BDF on stiff problems with a cheap Jacobian and suits
frequent state discontinuities at moderate dimension.

### Tsitouras (`tsit5`)

The explicit \\5(4)\\ pair of Tsitouras [\[4\]](#ref-TSITOURAS2011770)
assumes only the first column simplifying condition, \\ \sum\_{i=1}^{s}
b_i\\ a\_{ij} \\=\\ b_j\\(1 - c_j), \qquad j = 1, \dots, s, \\ not the
row condition \\\sum\_{j} a\_{ij}\\ c_j = c_i^2/2\\, \\i \ge 3\\, of the
Dormand-Prince pair. The freed parameters lower the leading-error norm;
the
[coefficients](https://dmodverse.github.io/cppDE/articles/methods-coefficients.html#tsit5-coeffs)
meet the column condition to rounding and violate the row condition. The
\\s = 7\\ stages are First-Same-As-Last (FSAL), \\ k_7 \\=\\
f(x\_{n+1},\\ t_n + h_n) \\=\\ k_1 \text{ at step } n+1, \\ so an
accepted step costs six right-hand-side evaluations.

| symbol | meaning |
|:---|:---|
| \\a\_{ij}\\ | weights of earlier stages in stage \\i\\ (strictly lower triangular) |
| \\c_i\\ | time abscissa of stage \\i\\, \\c_i = \sum\_{j\<i} a\_{ij}\\ |
| \\b_i\\ | weights of the order-5 solution; \\b_7 = 0\\ |
| \\\hat{b}\_i\\ | weights of the embedded order-4 estimator; \\\hat{b}\_7 \ne 0\\ |
| \\k_i\\ | stage \\i\\: \\k_i = f\\\left(x_n + h_n\sum\_{j\<i}a\_{ij}k_j,\\ t_n + c_i h_n\right)\\ |

Symbols of the Tsit5 tableau. {.table}

The advanced solution, embedded estimate and local error estimate are \\
x\_{n+1} = x_n + h_n \sum\_{i=1}^{7} b_i\\ k_i, \qquad \hat{x}\_{n+1} =
x_n + h_n \sum\_{i=1}^{7} \hat{b}\_i\\ k_i, \\ \\ e\_{n+1} = h_n
\sum\_{i=1}^{7} (\hat{b}\_i - b_i)\\ k_i. \\ Without Jacobian or linear
solve, Tsit5 advances furthest per right-hand-side evaluation of the
four methods and suits non-stiff problems with an expensive right-hand
side.

## Adaptive step-size control

The error estimate is reduced to a scalar by the weighted
root-mean-square norm \\ \\ e\_{n+1} \\\_{\mathrm{WRMS}} =
\sqrt{\frac{1}{n_x} \sum\_{j=1}^{n_x} \left( \frac{e\_{n+1,j}}{w\_{j}}
\right)^{\\2}}, \qquad w\_{j} = \mathtt{abstol} +
\mathtt{reltol}\\\|x^{\mathrm{ref}}\_{j}\|, \tag{8} \\ where
\\x^{\mathrm{ref}}\\ is \\x_n\\ in BDF/NDF and Adams, as in CVODE, and
\\\max(\|x_n\|, \|x\_{n+1}\|)\\ componentwise in the single-step
methods. A step is accepted when \\\\ e\_{n+1} \\\_{\mathrm{WRMS}} \le
1\\.

The single-step methods use the proportional-integral (PI) controller of
[\[3\]](#ref-hairer2010solving). With \\\varepsilon\_{n+1} = \\ e\_{n+1}
\\\_{\mathrm{WRMS}}\\ and \\\hat q\\ the embedded order (\\3\\
Rosenbrock4, \\4\\ Tsit5), an accepted step sets \\ h\_{n+1} = h_n\\
\min\Bigl\\\rho\_{\max},\\ \max\Bigl(\rho\_{\min},\\ S\\
\varepsilon\_{n+1}^{-(k_P+k_I)}\\ \varepsilon\_{n}^{\\k_I}\Bigr)\Bigr\\,
\qquad k_P = \frac{0.7}{\hat q+1}, \quad k_I = \frac{0.4}{\hat q+1}, \\
with \\S = 0.9\\, \\\rho\_{\min} = 0.2\\, \\\rho\_{\max} = 5\\, and
\\\varepsilon\_{n}\\ the norm of the previous accepted step, floored at
\\0.01\\. The first step and the first accepted step after a rejection
use the integral factor \\S\\\varepsilon\_{n+1}^{-1/(\hat q+1)}\\,
capped at \\1\\ after a rejection. A rejected step is retried with
\\S\\\varepsilon\_{n+1}^{-1/(\hat q+1)}\\ clamped to \\\[\rho\_{\min},
0.9\]\\.

The multistep methods select step size and order as CVODE. With
\\E\_{q'}\\ the error norm at order \\q'\\, the candidate ratios are \\
\eta\_{q'} = \frac{1}{\bigl(B\_{q'}\\ E\_{q'}\bigr)^{1/(q'+1)} +
10^{-6}}, \qquad q' \in \\q-1,\\ q,\\ q+1\\, \\ with \\B\_{q-1} = B_q =
6\\ and \\B\_{q+1} = 10\\. Neighboring orders are compared \\q+1\\ steps
after a change to order \\q\\ and every second step at unchanged order;
the largest ratio fixes the next order. A ratio below \\1.5\\ leaves
\\h\\ unchanged; the cap is \\10^4\\ on the first step, \\10\\
afterwards.

The order is at most \\q\_{\max}\\ (\\5\\ for BDF/NDF, \\12\\ for Adams)
or a smaller `maxOrder`. For \\q \ge 3\\ the stability domain excludes a
wedge around the imaginary axis; at the edge of this domain the error
test rejects steps and shrinks \\h\\ where a lower order would be
stable. `stald = TRUE` adds CVODE stability limit detection to the BDF
family [\[1\]](#ref-hindmarsh2005sundials):

- from order \\3\\ on, each step records the squared norms of the scaled
  derivatives of degree \\q-1\\, \\q\\ (Nordsieck columns) and \\q+1\\
  (corrector correction);
- after \\q+5\\ steps at order \\q\\, the last five records give the
  magnitude \\r\\ of the dominant root of step-to-step error
  propagation;
- \\r \> 0.98\\ lowers the order by one at step size \\\eta\_{q-1}
  h_n\\, overriding the order selection unless that already lowers the
  order.

Off by default as in CVODE, it saves steps on weakly damped stiff
oscillations and costs the estimate on every step. On a
[`cvode()`](https://dmodverse.github.io/cppDE/reference/cvode.md) model,
`maxOrder` and `stald` set `CVodeSetMaxOrd` and `CVodeSetStabLimDet`.

### The first step

The single-step methods follow [\[5\]](#ref-hairer2008solving),
Algorithm II.4.4. A rough \\h\_{\mathrm{s}} = 0.01\\\\x_0\\ / \\f_0\\\\,
\\f_0 = f(x_0, p, t_0)\\, fixes a scale, the second derivative \\\ddot
x\\ along the trajectory fixes the curvature, and \\ h\_{\mathrm{c}} =
\left( \frac{0.01}{\max\bigl(\\f_0\\, \\\ddot x\\\bigr)} \right)^{\\1/2}
\tag{9} \\ combines the two. The exponent is that of a first-order
method: \\\\\ddot x\\\\ stands in for the derivative of order \\p+1\\,
exact only at \\p = 1\\; higher orders raise the step within a few
steps. Generated code gives \\\ddot x = \partial f/\partial t + J f\\,
the derivative of \\f\\ along \\(f, 1)\\. Norms are weighted maximum
norms, weights \\\mathtt{abstol} + \mathtt{reltol}\\\|x_0\|\\, including
the sensitivities with `sensErrCon = TRUE`. The step is
\\\min(100\\h\_{\mathrm{s}},\\ h\_{\mathrm{c}},\\ h\_{\mathrm{ub}})\\.
The upper bound \\h\_{\mathrm{ub}}\\, at most a tenth of the interval,
keeps one Euler step from changing any component by more than a fraction
of the component scale; a lower bound keeps \\t_0 + h \ne t_0\\ in
double precision.

The multistep methods use the bounds of CVODES `cvHin`:
\\h\_{\mathrm{lb}}\\ above the rounding of the time variable, the same
\\h\_{\mathrm{ub}}\\, and the step \\\sqrt{2/\\\ddot x\\}\\ in the
weighted RMS norm if below \\h\_{\mathrm{ub}}\\, else
\\h\_{\mathrm{lb}}^{1/4} h\_{\mathrm{ub}}^{3/4}\\, halved and clamped to
the bounds. `cvHin` estimates \\\ddot x\\ by differences of \\f\\ over
trial steps until two estimates agree, one round with the exact \\\ddot
x\\. After every event each method estimates a new first step from the
post-event state.

## Dense output

Each accepted step interpolates output times inside the step without
right-hand-side evaluations. With \\\tau = (t - t_n)/h_n \in \[0, 1\]\\:

- Tsit5: the cubic Hermite interpolant through \\x_n\\, \\x\_{n+1}\\ and
  the FSAL derivatives \\k_1\\, \\k_7\\, \\ x(t) = (1+2\tau)(1-\tau)^2\\
  x_n + \tau^2(3-2\tau)\\ x\_{n+1} + \tau(1-\tau)^2\\ h_n k_1 +
  \tau^2(\tau-1)\\ h_n k_7 . \\
- Rosenbrock4: the stage-based cubic of Hairer and Wanner
  [\[3\]](#ref-hairer2010solving), \\ x(t) = (1-\tau)\\ x_n +
  \tau\bigl(x\_{n+1} + (1-\tau)(u_3 + \tau\\ u_4)\bigr), \qquad u_3 =
  \sum\_{i=1}^{5} d_i^{(2)} k_i, \quad u_4 = \sum\_{i=1}^{5} d_i^{(3)}
  k_i, \\ with the [dense-output
  coefficients](https://dmodverse.github.io/cppDE/articles/methods-coefficients.html#ros-coeffs).
- BDF/NDF and Adams: the Nordsieck polynomial of degree \\q\\ after the
  step, slots \\z_j\\ at \\t\_{n+1}\\ scaled with \\h_n\\, \\ x(t) =
  \sum\_{j=0}^{q} z\_{j}\\ \Bigl(\frac{t - t\_{n+1}}{h_n}\Bigr)^{j}. \\

The interpolant order is \\3\\ for the single-step methods, below the
step order (\\5\\ Tsit5, \\4\\ Rosenbrock4), and \\q\\ for the multistep
methods.

## References

\[1\]

A. C. Hindmarsh *et al.*, “SUNDIALS: Suite of nonlinear and
differential/algebraic equation solvers,” *ACM Transactions on
Mathematical Software (TOMS)*, vol. 31, no. 3, pp. 363–396, 2005, doi:
[10.1145/1089014.1089020](https://doi.org/10.1145/1089014.1089020).

\[2\]

L. F. Shampine and M. W. Reichelt, “The MATLAB ODE suite,” *SIAM Journal
on Scientific Computing*, vol. 18, no. 1, pp. 1–22, 1997, doi:
[10.1137/S1064827594276424](https://doi.org/10.1137/S1064827594276424).

\[3\]

E. Hairer and G. Wanner, *Solving ordinary differential equations II:
Stiff and differential-algebraic problems*, 2nd ed. in Springer series
in computational mathematics, no. 14. Springer, 1996. doi:
[10.1007/978-3-642-05221-7](https://doi.org/10.1007/978-3-642-05221-7).

\[4\]

Ch. Tsitouras, “Runge–Kutta pairs of order 5(4) satisfying only the
first column simplifying assumption,” *Computers & Mathematics with
Applications*, vol. 62, no. 2, pp. 770–775, 2011, doi:
[10.1016/j.camwa.2011.06.002](https://doi.org/10.1016/j.camwa.2011.06.002).

\[5\]

E. Hairer, S. P. Nørsett, and G. Wanner, *Solving ordinary differential
equations I: Nonstiff problems*, 2nd ed. in Springer series in
computational mathematics, no. 8. Springer, 1993. doi:
[10.1007/978-3-540-78862-1](https://doi.org/10.1007/978-3-540-78862-1).
