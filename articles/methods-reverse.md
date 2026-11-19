# Reverse mode: the discrete adjoint

The reverse mode of
[`cppODE()`](https://dmodverse.github.io/cppDE/reference/cppODE.md)
returns a gradient at a cost independent of the number of parameters, as
the discrete adjoint of the computed trajectory;
[`cvode()`](https://dmodverse.github.io/cppDE/reference/cvode.md)
computes the continuous adjoint. Symbols follow the [notation
table](https://dmodverse.github.io/cppDE/articles/methods-overview.html#notation).
Used from R through `derivMode = "reverse"` of
[`cppODE()`](https://dmodverse.github.io/cppDE/reference/cppODE.md), the
arguments `cotangent`, `keepStore`, `store` and `adjoint` of
[`solveODE()`](https://dmodverse.github.io/cppDE/reference/solveODE.md),
and
[`adjointControl()`](https://dmodverse.github.io/cppDE/reference/adjointControl.md)
and
[`releaseStore()`](https://dmodverse.github.io/cppDE/reference/releaseStore.md).

## Gradient of one objective

Parameter estimation needs the gradient of one number, the objective
function \\\mathcal{J}(\theta)\\, a likelihood or a residual sum of
squares, in all \\M\\ parameters \\\theta\\. With [forward
sensitivities](https://dmodverse.github.io/cppDE/articles/methods-forward.md),
one tangent per parameter, it costs about \\M\\ solves; the reverse mode
needs one forward solve and one backward sweep, whatever \\M\\
[\[1\]](#ref-griewank2008evaluating), [\[2\]](#ref-naumann2012art). With
the solver steps \\\theta \to x_1 \to \dots \to x_N\\, the gradient is
the product of the Jacobians, \\
\frac{\mathrm{d}\mathcal{J}}{\mathrm{d}\theta} = \frac{\partial
\mathcal{J}}{\partial x_N}\\ \frac{\partial x_N}{\partial x\_{N-1}}
\cdots \frac{\partial x_2}{\partial x_1}\\ \frac{\partial x_1}{\partial
\theta} . \\ The forward mode multiplies from the right, through the
\\n_x \times M\\ sensitivities. The reverse mode multiplies from the
left, through rows of length \\n_x\\, the *adjoint* \\\lambda\\. The
factors depend on the states, so the reverse mode solves forward, stores
what the factors need, then walks back.

A reverse solve returns a contraction: for a cotangent \\w \in
\mathbb{R}^{n\_{\mathrm{obs}} \times n_x}\\ on the observed trajectory,
one sweep returns \\ \bar\theta_j = \sum\_{o=1}^{n\_{\mathrm{obs}}}
\sum\_{i=1}^{n_x} w\_{oi} \\ \frac{\partial x_i(t_o)}{\partial
\theta_j}, \qquad j = 1,\dots,M, \tag{1} \\ with \\\theta = (x_0, p)\\,
\\M = n_x + n_p\\, and \\w_o \in \mathbb{R}^{n_x}\\ the row of \\w\\ at
\\t_o\\. For \\w\_{oi} = \partial \mathcal{J}/\partial x_i(t_o)\\ (for
least squares, the weighted residual), \\\bar\theta = \nabla\_\theta
\mathcal{J}\\. Several cotangents, as columns of a three-dimensional
array, take one sweep each.

### Example: two Euler steps

For \\\dot x = -k x\\, two explicit Euler steps of size \\h\\ and
\\\mathcal{J} = x_2\\, each step multiplies by \\a = 1 - hk\\: \\ x_1 =
a\\x_0, \qquad x_2 = a\\x_1, \qquad \mathcal{J} = x_2 = a^2 x_0 . \\ A
step has the derivatives \\\partial x\_{n+1}/\partial k = -h\\x_n\\ and
\\\partial x\_{n+1}/\partial x_n = a\\.

*Forward.* The tangent \\s_n = \mathrm{d}x_n/\mathrm{d}k\\ starts at
\\s_0 = 0\\: \\ s_1 = a\\s_0 - h\\x_0 = -h\\x_0, \qquad s_2 = a\\s_1 -
h\\x_1 = -2h\\a\\x_0 . \\ The derivative in \\x_0\\ needs a second
tangent, started at \\1\\.

*Reverse.* The adjoint \\\lambda_n = \partial \mathcal{J}/\partial x_n\\
starts at \\\lambda_2 = 1\\. Each transposed step passes
\\a\lambda\_{n+1}\\ to the previous state and adds
\\-h\\x_n\lambda\_{n+1}\\ to the parameter cotangent \\\bar k\\: \\
\lambda_1 = a\\\lambda_2 = a, \qquad \bar k = -h\\x_1\lambda_2 =
-h\\a\\x_0, \\ \\ \lambda_0 = a\\\lambda_1 = a^2, \qquad \bar k =
-h\\a\\x_0 - h\\x_0\lambda_1 = -2h\\a\\x_0 . \\ One backward pass
returns both derivatives and reads \\x_1\\ and \\x_0\\, which the
forward run stores. A further parameter adds a tangent in the forward
mode and one entry of \\\bar\theta\\ in the reverse mode.

\\-2h(1-hk)\\x_0\\ is the exact derivative of the computed
\\\mathcal{J}\_h = (1-hk)^2 x_0\\, not of the exact solution \\x(2h) =
e^{-2hk}x_0\\, whose derivative is \\-2h\\e^{-2hk}x_0\\. At \\h = 0.5\\,
\\k = 1\\, \\x_0 = 1\\: \\\mathcal{J}\_h = 0.25\\ with derivative
\\-0.5\\, exact value \\0.368\\ with derivative \\-0.368\\ ([Discrete
and continuous adjoints](#discrete-and-continuous-adjoints)).

## Adjoint equation

The time-continuous form evaluates [(1)](#eq-adjoint-contraction)
without the sensitivity \\S_p\\ of [Forward
(1)](https://dmodverse.github.io/cppDE/articles/methods-forward.html#eq-sensitivity-p).
The adjoint \\\lambda(t) \in \mathbb{R}^{n_x}\\ runs backwards from the
end of the observed interval, \\ \dot\lambda = -J(x, p, t)^{\mathsf
T}\lambda, \qquad \lambda(T^+) = 0, \qquad \lambda(t_o^-) =
\lambda(t_o^+) + w_o, \tag{2} \\ the counterpart of \\\lambda_n =
a\\\lambda\_{n+1}\\. Along the solution the Jacobian terms cancel, \\
\frac{\mathrm{d}}{\mathrm{d}t} \left(\lambda^{\mathsf T} S_p\right) =
-\lambda^{\mathsf T} J S_p + \lambda^{\mathsf T}\\ \left(J S_p +
\frac{\partial f}{\partial p}\right) = \lambda^{\mathsf T}
\frac{\partial f}{\partial p}, \tag{3} \\ so \\\lambda\\ does not depend
on the choice of \\\theta\\. Integrating [(3)](#eq-adjoint-identity)
over \\\[t_0, T\]\\, with the jumps giving the terms of
[(1)](#eq-adjoint-contraction), \\ \bar\theta = \left(\frac{\partial
x_0}{\partial \theta}\right)^{\\\mathsf T}\\ \lambda(t_0) +
\left(\frac{\partial p}{\partial \theta}\right)^{\\\mathsf T}\\
\int\_{t_0}^{T} \left(\frac{\partial f}{\partial p}\right)^{\\\mathsf
T}\\\lambda(t)\\ \mathrm{d}t . \tag{4} \\ For \\\theta = (x_0, p)\\ the
outer factors are identities, and a reverse solve reports \\ \bar x_0 =
\lambda(t_0), \qquad \bar p = \int\_{t_0}^{T} \left(\frac{\partial
f}{\partial p}\right)^{\\\mathsf T}\\\lambda(t)\\ \mathrm{d}t, \tag{5}
\\ the counterparts of \\\lambda_0\\ and \\\bar k\\, stacked in the row
order of the `tangent` argument, states first. A reparameterization
\\(x_0, p) = \Phi(\theta)\\ enters afterwards through the outer factors
of [(4)](#eq-adjoint-result), the chain rule of [Forward
(2)](https://dmodverse.github.io/cppDE/articles/methods-forward.html#eq-sensitivity-theta)
applied on the left. The cost is one backward solve of dimension \\n_x\\
and a quadrature of dimension \\n_p\\, whatever \\M\\; \\\partial
f/\partial p\\ enters only contracted with \\\lambda\\.

## Discrete adjoint in cppODE()

`cppODE(..., derivMode = "reverse")` does not integrate
[(2)](#eq-adjoint-ode). As in the Euler example, the forward run stores
one checkpoint per accepted step. The backward sweep applies the
transposed accepted steps, last to first, and accumulates the parameter
part of [(5)](#eq-adjoint-quadrature). Step sizes are constants of the
differentiated map ([Internal numerical
differentiation](https://dmodverse.github.io/cppDE/articles/methods-reverse-checks.html#ind)).

### Adjoint of the step map

An adjoint sums over all later uses of a quantity, so unlike a dual
number it is not computed alongside the value. cppDE keeps no tape: the
step map is one accepted step of a known method on a fixed grid, and the
transpose is written in closed form. A multistep step consists of the
Nordsieck rescale and Pascal shift, the corrector, and the controller
tail: \\ z^{\text{pred}} = A\\ z^{\text{in}}, \tag{6} \\ \\ \Delta^{(0)}
= 0, \qquad \Delta^{(m+1)} = \Delta^{(m)} + s\\ W^{-1}
r\bigl(\Delta^{(m)}, z^{\text{pred}}, \theta\bigr), \quad m = 0, \dots,
k-1, \tag{7} \\ \\ z^{\text{out}} = B\\ z^{\text{pred}} + c\\ \Delta,
\qquad \Delta = \Delta^{(k)}, \tag{8} \\ with \\z\\ the Nordsieck array,
\\\Delta\\ the correction, residual \\r(\Delta, z, \theta) = f(z_0 +
\Delta, t\_{n+1}, \theta) - (\ell z_1 + \Delta)/\gamma\\ with step
coefficients \\\gamma = h_n\beta_0\\ and \\\ell\\, iteration matrix
\\W\\ and scale \\s\\ of [Stepping
(3)-(4)](https://dmodverse.github.io/cppDE/articles/methods-stepping.html#eq-newton),
\\k\\ iterations, and corrected state \\y = z_0^{\text{pred}} +
\Delta\\. \\A\\, \\B\\ and \\c\\ act on the slot index only, equally for
every state component, with coefficients from the *Nordsieck carry*
(order, step-size history and counters). They are read off the stepper:
a probe whose states are the Nordsieck slots, seeded with the identity
in the \\\[\text{slot} \times \text{component}\]\\ array, runs through
the forward-step routines and returns every column in one pass, so the
adjoint uses exactly the forward coefficients.

Under `bdf` the sweep differentiates the corrector equation \\r = 0\\ at
the accepted state \\y\\: the cotangent \\\mu\\ of the residual solves
\\ \bigl(I/\gamma - J(y)\bigr)^{\mathsf T} \mu = \bar\Delta , \\ and
\\J(y)^{\mathsf T}\mu\\, \\-(\ell/\gamma)\\\mu\\ and \\(\partial
f/\partial\theta)(y)^{\mathsf T}\mu\\ pass to \\z_0^{\text{pred}}\\,
\\z_1^{\text{pred}}\\ and \\\theta\\. \\\mu\\ is summed from the
transpose of the iteration [(7)](#eq-step-corr) at \\y\\, \\ \mu_m = s\\
W^{-\mathsf T} \bar\Delta^{(m)}, \qquad \bar\Delta^{(m+1)} =
\bar\Delta^{(m)} + J(y)^{\mathsf T}\mu_m - \mu_m/\gamma, \qquad
\bar\Delta^{(0)} = \bar\Delta , \\ one transposed solve and one product
with \\J^{\mathsf T}\\ per iteration, until the WRMS norm of
\\\bar\Delta^{(m)}\\ against \\\mathtt{reltol}\\\|\lambda\| +
\mathtt{abstol}\\ is at most \\0.1\\, with \\\lambda\\ the cotangent the
step hands down: the share of the local error test of the backward
problem that CVODES gives the Newton iteration
[\[3\]](#ref-hindmarsh2005sundials). Where the frozen \\W\\ cannot reach
the test within ten iterations, the iteration continues against \\W\\
built at \\y\\ and the current \\\gamma\\. The \\k\\ iterations the
forward run took leave the tangent of \\\Delta\\ unconverged even where
\\\Delta\\ has converged, and on a long stiff trajectory their
derivative can be off by percents.

Under `adams` the corrector is differentiated as the \\k\\ forward
iterations [(7)](#eq-step-corr), with \\W\\, the step size, the order
and \\k\\ held fixed. From \\\bar\Delta^{(k)} = \bar\Delta\\, each
backward iteration takes one transposed solve and one product with
\\J^{\mathsf T}\\ at \\y^{(m)} = z_0^{\text{pred}} + \Delta^{(m)}\\: \\
\mu_m = s\\ W^{-\mathsf T} \bar\Delta^{(m+1)}, \qquad \bar\Delta^{(m)} =
\bar\Delta^{(m+1)} + J(y^{(m)})^{\mathsf T}\mu_m - \mu_m/\gamma, \\ with
the three terms of each \\\mu_m\\ at \\y^{(m)}\\ passed on as above. The
forward sensitivities with `sensErrCon = FALSE` take the same iterations
on dual numbers ([Forward
sensitivities](https://dmodverse.github.io/cppDE/articles/methods-forward.html#corrector)),
so both differentiate the same map. An exactly solved Adams-Moulton step
can be unstable on a stiff model at the step sizes the value run chose
while the iterations taken are not: at order 3 a parasitic root leaves
the unit disk for \\h\lambda \< -6\\, and on Robertson’s kinetics at
\\h\lambda \approx -30\\ the solved sweep grows by the modulus of that
root, about \\1.5\\, per step.

A Rosenbrock step transposes to a stage recursion of six solves against
one factorization, an explicit step to the backward recursion over the
tableau.

The code generator emits \\J^{\mathsf T}\lambda\\ and \\(\partial
f/\partial p)^{\mathsf T}\lambda\\ beside the right-hand side. For a
Rosenbrock method it also emits both applied to \\J v\\, a second
derivative, since the stages solve against a matrix built from the
Jacobian.

### Checkpoints and memory

A checkpoint holds the step-start state, the Nordsieck carry, time, step
size, the controller operations applied to the Nordsieck carry
afterwards and, for a multistep method, the corrector iterates. Memory
is linear in the number of steps and holds no derivatives. Each
iteration matrix is recorded once, as the point, time and \\\gamma_p\\
of the factorization, and refactorized in the sweep once per run of
steps that used it; a `bdf` step whose corrector adjoint does not meet
the test against it within ten iterations adds one, built at the
accepted state. Slot operators are cached by Nordsieck carry and step
size, so the probe runs only where the controller changed the control
state.

An observation at an interpolated time goes through the continuous
extension of the step containing it. After the tail the probe holds a
valid interpolant, whose value at the observation time is the row
\\\partial x\_{\text{interp}} / \partial (z^{\text{pred}}, \Delta)\\.

### Reuse of checkpoints

The sweep reads only the checkpoints, so value run and sweep can be
separate calls, as a gradient needs the values to form the cotangent.

- `keepStore = TRUE` without a cotangent runs the values and returns the
  checkpoints as `$store`.
- `store` passes them to a later solve at the same times and parameters,
  which goes directly to the sweep; a store from other parameter values
  is rejected.
- A three-dimensional cotangent sweeps one store per column:
  \\n\_{\text{seed}}\\ sweeps, one integration.
- The checkpoints live until the garbage collector frees the store or
  [`releaseStore()`](https://dmodverse.github.io/cppDE/reference/releaseStore.md)
  frees them; a released store makes the solve integrate again.
- `derivMode = "forward-reverse"` refuses both arguments, since
  checkpoints of that mode hold tangents that do not outlive the
  recording solve.

The reverse mode supports all four methods, time- and root-triggered
events ([Events in reverse
mode](https://dmodverse.github.io/cppDE/articles/methods-reverse-events.md)),
forcings, the batch entry point and the sparse linear solver. The code
is in `cppde_adjoint_step.hpp`, `cppde_reverse_step.hpp` and
`cppde_reverse_trajectory.hpp`.

## Discrete and continuous adjoints

The discrete adjoint of
[`cppODE()`](https://dmodverse.github.io/cppDE/reference/cppODE.md)
discretizes, then differentiates: it is the derivative of the steps the
forward solver took. The continuous adjoint of
[`cvode()`](https://dmodverse.github.io/cppDE/reference/cvode.md), the
adjoint sensitivity analysis (ASA) of CVODES, differentiates, then
discretizes: a second integrator with separate step-size control solves
[(2)](#eq-adjoint-ode) and [(5)](#eq-adjoint-quadrature) backwards over
forward states interpolated between checkpoints
[\[4\]](#ref-cao2003adjoint) ([cvode() and the CVODES
adjoint](https://dmodverse.github.io/cppDE/articles/methods-cvode.html#asa)).

With \\\mathcal{J}\\ the objective of the exact solution and
\\\mathcal{J}\_h\\ that of the computed one, the discrete adjoint is
\\\nabla\_\theta \mathcal{J}\_h\\ on the value-run grid. Under `adams`,
`rb4` and `tsit5`, \\\mathcal{J}\_h\\ is the value the solver returns,
with corrector iterations and matrices as taken, and the adjoint equals
the forward sensitivities on the same grid (`sensErrCon = FALSE`) up to
rounding. Under `bdf`, \\\mathcal{J}\_h\\ is the value of the method
with the corrector equations solved, which the returned value meets to
the corrector test; the forward sensitivities on the same grid differ
from the adjoint by what their corrector iterations leave. The
continuous adjoint approximates \\\nabla\_\theta \mathcal{J}\\ to the
discretization error of both solves, \\O(\mathrm{tol})\\ away from
\\\nabla\_\theta \mathcal{J}\_h\\; in the Euler example it would
approximate \\-0.368\\, not \\-0.5\\.

|  | discrete adjoint | continuous adjoint |
|:---|:---|:---|
| call | `cppODE(..., derivMode = "reverse")` | `cvode(..., derivMode = "reverse")` |
| order | discretize, then differentiate | differentiate, then discretize |
| backward pass | transposed steps, forward grid | second integrator for [(2)](#eq-adjoint-ode) |
| forward states | checkpoint per step | checkpoint every `asaCheckpoints` steps, interpolated |
| implicit step | `bdf`: corrector equation solved; `adams`: \\k\\ transposed iterations | stiff linear ODE under error control |
| event | transposed jump arithmetic | restart with transposed saltation |
| result | \\\nabla\_\theta \mathcal{J}\_h\\ on the forward grid | \\\nabla\_\theta \mathcal{J}\\ to \\O(\mathrm{tol})\\ |
| consistent with value | up to rounding, `bdf` to the corrector test | to \\O(\mathrm{tol})\\ |

The two adjoints. Both converge to \\\nabla\_\theta \mathcal{J}\\ as the
tolerance is lowered. {.table}

**Consistency with the value.** Line searches, trust-region ratios and
gradient stopping tests assume the gradient is the derivative of the
value. Near the optimum \\\nabla\_\theta \mathcal{J}\\ shrinks to the
size of an \\O(\mathrm{tol})\\ error, and an inconsistent gradient
points where the computed value does not decrease: the trust-region
ratio degrades, the line search fails, and a gradient tolerance below
the error is never reached. The discrete adjoint is the gradient of the
value the optimizer sees, under `bdf` of that value with the corrector
equations solved.

**Refined sweep.** `adjointControl(refine = TRUE)` holds the sweep to
the tolerances and differentiates the flow, the solution map of the ODE,
instead of the step map. Like the continuous adjoint, the result is
closer to \\\nabla\_\theta \mathcal{J}\\ and differs from
\\\nabla\_\theta \mathcal{J}\_h\\ by the order of the discretization
error. It pays off where the value grid is too coarse for the
derivative, such as a trajectory at rest crossed in a few steps ([Grid
and refinement of the
sweep](https://dmodverse.github.io/cppDE/articles/methods-reverse-checks.html#state-error)).
A refined solve also returns `$correction`, per seed column the
first-order change of the seeded functional from the computed to the
differentiated trajectory; value plus correction and the refined
gradient belong to one function ([Value correction of the refined
sweep](https://dmodverse.github.io/cppDE/articles/methods-reverse-checks.html#value-correction)).

## Example: Michaelis-Menten conversion

\\x\\ converts into \\y\\ by Michaelis-Menten kinetics, \\y\\ degrades
at first order and is observed at eleven time points, \\\mathcal{J} =
\sum_o \bigl(y(t_o) - d_o\bigr)^2\\. The cotangent \\w_o\\ holds
\\2\bigl(y(t_o) - d_o\bigr)\\ in the column of \\y\\. The
[`cppODE()`](https://dmodverse.github.io/cppDE/reference/cppODE.md)
models use `bdf` and `rb4`.

``` r

library(cppDE)
eqns  <- c(x = "-k1*x/(Km + x)", y = "k1*x/(Km + x) - k2*y")
parms <- c(x = 2, y = 0, k1 = 1, Km = 0.5, k2 = 0.3)
times <- 0:10
data  <- c(0, 0.7, 1.1, 1.2, 1.15, 1.05, 0.9, 0.8, 0.7, 0.6, 0.5)

methods <- c("bdf", "rb4")
fwd <- lapply(setNames(nm = methods), function(m)
  cppODE(eqns, method = m, modelname = paste0("mmFwd_", m), compile = FALSE))
rev <- lapply(setNames(nm = methods), function(m)
  cppODE(eqns, method = m, derivMode = "reverse",
         modelname = paste0("mmRev_", m), compile = FALSE))
do.call(compile, c(unname(c(fwd, rev)), list(output = "mmNative")))

hasCvode <- codegenAvailable("cvode")
if (hasCvode) {
  asa <- cvode(eqns, derivMode = "reverse", modelname = "mmAsa", compile = FALSE)
  compile(asa, output = "mmCvode")
}
```

The CVODES adjoint keeps no store, so the value run takes a zero
cotangent.

``` r

cotangent <- function(y) {
  w <- matrix(0, length(times), 2, dimnames = list(NULL, c("x", "y")))
  w[, "y"] <- 2 * (y - data)
  w
}
gradForward <- function(tol, method = "bdf", sensErrCon = FALSE) {
  out <- solveODE(fwd[[method]], times, parms, abstol = tol, reltol = tol,
                  sensErrCon = sensErrCon)
  2 * colSums((out$variable[, "y"] - data) * out$tangent[, "y", ])
}
gradReverse <- function(tol, method = "bdf", adjoint = NULL) {
  val <- solveODE(rev[[method]], times, parms, abstol = tol, reltol = tol,
                  keepStore = TRUE)
  solveODE(rev[[method]], times, parms, abstol = tol, reltol = tol,
           store = val$store, cotangent = cotangent(val$variable[, "y"]),
           adjoint = adjoint)$cotangent[, 1]
}
gradAsa <- function(tol) {
  if (!hasCvode) return(rep(NA_real_, length(parms)))
  val <- solveODE(asa, times, parms, abstol = tol, reltol = tol,
                  cotangent = cotangent(0 * data))
  solveODE(asa, times, parms, abstol = tol, reltol = tol,
           cotangent = cotangent(val$variable[, "y"]))$cotangent[, 1]
}

exact <- gradForward(1e-12, "rb4", sensErrCon = TRUE)
exact
#>          x          y         k1         Km         k2 
#> -1.4861036 -0.8366524  0.7868458 -0.5163659  9.5726795
```

A forward solve at `reltol = 1e-12` stands in for \\\nabla\_\theta
\mathcal{J}\\. The tables give the largest difference over the five
entries, relative to the largest entry.

``` r

gap <- function(a, b) max(abs(a - b)) / max(abs(b))
tols <- 10^-(4:10)
gapTable <- function(method) {
  t(sapply(tols, function(tol) {
    forward <- gradForward(tol, method)
    discrete <- gradReverse(tol, method)
    refined <- gradReverse(tol, method, adjointControl(refine = TRUE))
    c(reltol = tol,
      "discrete vs forward" = gap(discrete, forward),
      "discrete vs exact" = gap(discrete, exact),
      "refine vs forward" = gap(refined, forward),
      "refine vs exact" = gap(refined, exact),
      "continuous vs exact" = gap(gradAsa(tol), exact))
  }))
}
knitr::kable(formatC(gapTable("bdf"), format = "e", digits = 1), align = "r",
             caption = "bdf")
```

| reltol | discrete vs forward | discrete vs exact | refine vs forward | refine vs exact | continuous vs exact |
|---:|---:|---:|---:|---:|---:|
| 1.0e-04 | 1.3e-04 | 3.1e-03 | 1.3e-03 | 1.8e-03 | 2.0e-03 |
| 1.0e-05 | 1.0e-05 | 5.0e-04 | 2.4e-04 | 2.7e-04 | 9.0e-05 |
| 1.0e-06 | 1.4e-06 | 6.5e-05 | 2.6e-05 | 4.0e-05 | 5.0e-06 |
| 1.0e-07 | 8.7e-08 | 3.2e-06 | 2.6e-06 | 3.1e-06 | 1.6e-06 |
| 1.0e-08 | 3.0e-08 | 5.4e-08 | 6.1e-07 | 6.4e-07 | 6.0e-08 |
| 1.0e-09 | 6.0e-10 | 5.2e-09 | 1.0e-07 | 1.1e-07 | 4.5e-09 |
| 1.0e-10 | 9.6e-11 | 1.9e-09 | 9.0e-09 | 7.0e-09 | 8.4e-10 |

bdf {.table}

``` r

knitr::kable(formatC(gapTable("rb4"), format = "e", digits = 1), align = "r",
             caption = "rb4")
```

| reltol | discrete vs forward | discrete vs exact | refine vs forward | refine vs exact | continuous vs exact |
|---:|---:|---:|---:|---:|---:|
| 1.0e-04 | 1.7e-16 | 1.1e-05 | 1.3e-05 | 1.2e-05 | 2.0e-03 |
| 1.0e-05 | 5.6e-16 | 2.2e-06 | 7.5e-07 | 2.4e-06 | 9.0e-05 |
| 1.0e-06 | 3.7e-16 | 8.2e-08 | 1.9e-07 | 1.7e-07 | 5.0e-06 |
| 1.0e-07 | 1.7e-15 | 2.1e-08 | 2.5e-08 | 8.3e-09 | 1.6e-06 |
| 1.0e-08 | 4.4e-16 | 2.9e-09 | 3.0e-09 | 6.4e-10 | 6.0e-08 |
| 1.0e-09 | 4.6e-15 | 1.1e-10 | 3.0e-10 | 1.8e-10 | 4.5e-09 |
| 1.0e-10 | 8.2e-15 | 4.2e-11 | 4.0e-11 | 8.8e-12 | 8.4e-10 |

rb4 {.table}

Under `rb4` the discrete adjoint matches the forward sensitivities to
rounding, both being \\\nabla\_\theta \mathcal{J}\_h\\. Under `bdf` the
two differ by what the corrector iterations of the forward run leave, a
gap that falls with the tolerance. The refined sweep differs from them
by the order of the tolerance, being the derivative of the value plus
`$correction`.

## Choice of derivative mode

| need | choice |
|:---|:---|
| full sensitivity matrix \\S\\, as for Gauss-Newton; cost grows with \\M\\ | forward mode |
| one contraction per sweep, cost independent of \\M\\, gradient for an optimizer | discrete adjoint |
| SUNDIALS backend required | continuous adjoint |
| value grid does not resolve the derivative | `refine = TRUE` |
| second derivatives | [Second derivatives](https://dmodverse.github.io/cppDE/articles/methods-forward-reverse.md) |

The dMod2 article on trust-region optimization compares [Gauss-Newton
with forward derivatives against SR1 with reverse
derivatives](https://dmodverse.github.io/dMod2/articles/methods-trust.html#forward-and-reverse-derivatives-on-the-oppelt-fit)
on a real objective.

## References

\[1\]

A. Griewank and A. Walther, *Evaluating derivatives: Principles and
techniques of algorithmic differentiation*, 2nd ed. in Other titles in
applied mathematics, no. 105. SIAM, 2008. doi:
[10.1137/1.9780898717761](https://doi.org/10.1137/1.9780898717761).

\[2\]

U. Naumann, *The art of differentiating computer programs: An
introduction to algorithmic differentiation*. in Software, environments,
and tools, no. 24. SIAM, 2012. doi:
[10.1137/1.9781611972078](https://doi.org/10.1137/1.9781611972078).

\[3\]

A. C. Hindmarsh *et al.*, “SUNDIALS: Suite of nonlinear and
differential/algebraic equation solvers,” *ACM Transactions on
Mathematical Software (TOMS)*, vol. 31, no. 3, pp. 363–396, 2005, doi:
[10.1145/1089014.1089020](https://doi.org/10.1145/1089014.1089020).

\[4\]

Y. Cao, S. Li, L. Petzold, and R. Serban, “Adjoint sensitivity analysis
for differential-algebraic equations: The adjoint DAE system and its
numerical solution,” *SIAM Journal on Scientific Computing*, vol. 24,
no. 3, pp. 1076–1089, 2003, doi:
[10.1137/S1064827501380630](https://doi.org/10.1137/S1064827501380630).
