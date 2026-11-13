# Events, switches and saltation

Symbols are those of the [notation
table](https://dmodverse.github.io/cppDE/articles/methods-overview.html#notation).
R interface: `events` and the comparisons, `Heaviside()` and
[`sign()`](https://rdrr.io/r/base/sign.html) in the `rhs` of
[`cppODE()`](https://dmodverse.github.io/cppDE/reference/cppODE.md);
`roottol` and `maxroot` of
[`solveODE()`](https://dmodverse.github.io/cppDE/reference/solveODE.md).

A *time-triggered event* fires at \\t^{\*} = t_e(\theta)\\, a
closed-form expression. A *root-triggered event* fires where a root
function \\g_e(x, t, \theta)\\ crosses zero, at \\t^{\*}(\theta)\\
defined by \\g_e(x(t^{\*}), t^{\*}, \theta) = 0\\. Either applies a
reset \\ x^{+} = R(x^{-}, t), \\ such as a dose, a clip to non-negative
values or a change of regime.

## Restart after an event

Single-step methods (Rosenbrock4, Tsit5) resume at full order. Multistep
methods (BDF/NDF, Adams) restart at order \\1\\ and rebuild the
Nordsieck array over several steps, the dominant cost under frequent
events. Every method re-estimates the step size from the post-event
state with the port of the CVODES routine `cvHin` ([the first
step](https://dmodverse.github.io/cppDE/articles/methods-stepping.html#the-first-step)),
bounded by the remaining interval; for a single-step method this is the
whole restart.

## Sensitivity transport across a discontinuity

The update \\S\_\theta^{+} = (\partial R/\partial x)\\S\_\theta^{-}\\
misses the \\\theta\\-dependence of the firing time. With \\f^{-}\\ the
right-hand side at \\x^{-}\\ and \\f^{+}\\ that at \\R(x^{-}, t^{\*})\\,
it adds two drifts of order \\f^{\pm}\\(\partial t^{\*}/\partial
\theta)\\:

1.  *Pre-event drift:* \\x^{-}\\ flows along \\f^{-}\\ from the grid
    time \\t\_\*^{(0)} := \mathrm{scalar}(t^{\*})\\ to \\t^{\*}\\, over
    \\\delta t^{\*} = t^{\*} - t\_\*^{(0)}\\, zero in value but not in
    the \\\theta\\-derivatives.
2.  *Post-event drift:* \\R(x^{-}, t^{\*})\\ flows along \\f^{+}\\ from
    \\t^{\*}\\ back to the grid time, over \\-\delta t^{\*}\\.

To first order the two give the [saltation matrix](#saltation). cppDE
applies forward drift, reset and backward drift in dual arithmetic, each
drift one step of Heun’s method (explicit trapezoidal rule). The state
receives only the \\\theta\\-derivatives of the displacement, and the
first and second derivatives of \\R\\. Without sensitivities, \\R\\ is a
plain map. \\\delta t^{\*}\\ comes from \\t_e(\theta)\\ or, for a root,
from the implicit function theorem.

## Time-triggered events

The integrator runs to \\t\_\*^{(0)} := \mathrm{scalar}(t_e)\\ and
treats \\ \delta t^{\*} = t_e(\theta) - t\_\*^{(0)}, \qquad
\mathrm{scalar}(\delta t^{\*}) = 0, \tag{1} \\ as a dual number with
tangents \\\partial t_e/\partial \theta_a\\ and \\\partial^2
t_e/\partial \theta_a \partial \theta_b\\. The state \\x^{-}\\ at
\\t\_\*^{(0)}\\ has the pre-event sensitivities \\S\_\theta\\ as
tangents. The transport reads \\ \begin{aligned} &\text{forward:}\quad
f_1 = f(x^{-}, t\_\*^{(0)}), \quad x_E = x^{-} + f_1\\\delta t^{\*},
\quad f_2 = f(x_E, t^{\*}), \\ &\hphantom{\text{forward:}\quad} x^{\*} =
x^{-} + \tfrac{1}{2}(f_1 + f_2)\\\delta t^{\*}, \\ &\text{reset:}\quad
x^{\*}\_{+} = R(x^{\*}, t^{\*}), \\ &\text{backward:}\quad \tilde f_1 =
f(x^{\*}\_{+}, t^{\*}), \quad x_B = x^{\*}\_{+} - \tilde f_1\\\delta
t^{\*}, \quad \tilde f_2 = f(x_B, t\_\*^{(0)}), \\
&\hphantom{\text{backward:}\quad} x^{+} = x^{\*}\_{+} -
\tfrac{1}{2}(\tilde f_1 + \tilde f_2)\\\delta t^{\*}. \end{aligned}
\tag{2} \\ An Euler step suffices at first order, since the powers
\\(\delta t^{\*})^{k}\\, \\k \ge 2\\, have no first-order dual content.
At second order the curvature term \\ \tfrac{1}{2}\bigl(J\\f + \partial
f/\partial t\bigr)\\(\delta t^{\*})^{2} \\ remains, because \\(\delta
t^{\*})^{2}\\ has the mixed derivative \\ \frac{\partial^{2}(\delta
t^{\*})^{2}}{\partial \theta_a\\\partial \theta_b} = 2\\ \frac{\partial
t_e}{\partial \theta_a}\\ \frac{\partial t_e}{\partial \theta_b}. \\ The
second Heun evaluation, at shifted state and time \\t^{\*}\\, supplies
it with an error of \\O((\delta t^{\*})^{3})\\. Evaluating \\f_2\\ and
\\\tilde f_1\\ at \\t\_\*^{(0)}\\ would drop \\\tfrac{1}{2}(\partial
f/\partial t)\\(\delta t^{\*})^{2}\\ for an explicitly time-dependent
right-hand side.

## Root-triggered events

A root event fires at \\t^{\*}(\theta)\\ defined by \\ G(t^{\*}(\theta),
\theta) := g_e\bigl(x(t^{\*}(\theta), \theta), t^{\*}(\theta),
\theta\bigr) = 0, \tag{3} \\ in three steps: localization on the scalar
trajectory, \\\delta t^{\*}\\ by the implicit function theorem,
transport [(2)](#eq-fixed-heun).

### Localization of the root

\\g_e\\ is evaluated at the end of every accepted step and at every
output time. A crossing on \\\[t_n, t\_{n+1}\]\\ is accepted on a sign
change or on \\g_e(x\_{n+1}, t\_{n+1}) = 0\\ exactly. A zero at the
lower end belongs to the previous sample, so an event does not fire
again on the surface it was just applied on, where \\g_e\\ is zero up to
round-off.

The crossing is refined on the dense output of the step, without
right-hand-side evaluations, by false position with the interpolation
factor clamped to \\\[0.1, 0.9\]\\. The clamp bounds the stagnation of
the classical rule and keeps the speed on a nearly linear \\g_e\\. It
stops at a bracket narrower than the root tolerance, after fifty
iterations, or at \\g_e = 0\\. It returns as \\t\_\*^{(0)} :=
\mathrm{scalar}(t^{\*})\\ the end of the final bracket with the
pre-crossing sign of \\g_e\\.

### Dual residual of the event time

For smooth \\G\\ with \\\partial G/\partial t \neq 0\\, the implicit
function theorem applied to [(3)](#eq-root-condition) gives \\
\frac{\partial t^{\*}}{\partial \theta} = -\frac{G\_\theta}{G_t} =
-\frac{\nabla_x g_e \cdot S\_\theta + \partial g_e/\partial \theta}
{\nabla_x g_e \cdot f + \partial g_e/\partial t}\bigg\|\_{\*}, \tag{4}
\\ with the rate \\\dot g_e\\ at the crossing as denominator. A
transversal crossing has \\\dot g_e \neq 0\\; at a grazing one
\\t^{\*}\\ is not differentiable in \\\theta\\ and the fallback of
simultaneous events applies. In dual arithmetic [(4)](#eq-dt-star-first)
is the quotient \\ \delta t^{\*} = -\frac{g_e(x^{-}, t\_\*^{(0)},
\theta)} {\dot g_e(x^{-}, t\_\*^{(0)}, \theta)}, \tag{5} \\ on the dual
state and parameters at the scalar grid time. The generated partials
\\\nabla_x g_e\\ and \\\partial g_e/\partial t\\ serve only \\\dot
g_e\\; the dual arithmetic forms the numerator. The scalar value of the
numerator, of the order of the root tolerance, is set to zero before the
division, so \\\delta t^{\*}\\ has only derivatives.

At second order the quotient lacks the \\G\_{tt}\\ contribution, since
\\g_e\\ is evaluated at \\t\_\*^{(0)}\\, not at the dual time
\\t^{\*}\\. Expanding [(3)](#eq-root-condition) with \\\delta t^{\*} =
\delta t^{\*}\_{(1)} + \delta t^{\*}\_{(2)} + \cdots\\ and
\\G(t\_\*^{(0)}, \theta^{(0)}) = 0\\ gives \\ 0 = G_t\\\delta
t^{\*}\_{(2)} + \tfrac{1}{2} G\_{tt}\\(\delta t^{\*}\_{(1)})^{2} +
G\_{t\theta}\\\delta t^{\*}\_{(1)}\\\delta\theta + \tfrac{1}{2}
G\_{\theta\theta}\\\delta\theta^{2}, \\ of which the quotient supplies
the \\G\_{t\theta}\\ and \\G\_{\theta\theta}\\ terms. The corrected
update, which changes only second-order components, is \\ \delta t^{\*}
\leftarrow \delta t^{\*} - \tfrac{1}{2}\\\frac{G\_{tt}}{G_t}\\(\delta
t^{\*})^{2}, \tag{6} \\ with \\ G\_{tt} =
\frac{\mathrm{d}^{2}}{\mathrm{d} t^{2}} g_e\bigl(x(t), t,
\theta\bigr)\bigg\|\_{\*} \\ from generated code, including second time
derivatives of forcings. Where a slide ends, \\f\\ is continuous and
\\G\_{tt}\\ enters scaled by a zero jump of \\f\\. The transport
[(2)](#eq-fixed-heun) then uses \\t^{\*} = t\_\*^{(0)} + \delta
t^{\*}\\.

### Simultaneous events

Crossings in one step are localized separately; those sharing
\\t\_\*^{(0)}\\ within the localization tolerance form a cluster. One
forward shift, every reset of the cluster (the resets are assumed to
commute) and one backward shift cost four right-hand-side evaluations
for any number of events. \\\delta t^{\*}\\ comes from the first event
of the cluster with non-zero \\\dot g_e\\ and is shared by the others,
which sit at the same grid time with the same parameter dependence
through \\x^{-}\\.

If every event of the cluster has \\\dot g_e = 0\\, detected by the size
of \\\|\dot g_e\|\\, there is no \\\partial t^{\*}/\partial \theta\\ and
the reset is a plain map, exact in the limit where both
\\f^{\pm}\\(\partial t^{\*}/\partial \theta)\\ vanish.

### Conditions switched by a reset

A jump is not a crossing of the continuous solution, so every condition
is evaluated on both sides of a jump:

- a condition whose sign changes across the jump fires there, one per
  sweep, so that a reset can switch on the next condition;
- a condition that held before the jump does not fire again;
- terminal conditions stop the integration instead.

Switched resets are applied between the forward and backward shift of
the causing event. One transport then applies the composition of both
resets, with the sensitivities of a time-triggered event at that time.

## Switches on states and time

\\f\\ switches where \\g_e(x, t, \theta)\\, reading a state or the time,
changes sign: in a comparison whose sides differ by \\g_e\\, and in
\\\mathrm{Heaviside}(g_e)\\ or \\\mathrm{sign}(g_e)\\. Read literally,
\\f\\ is discontinuous in \\x\\: a step across the switch fails the
error test until short enough, and sensitivities miss the dependence of
the switching time on \\\theta\\.

The code generator keeps each condition as a *mode* \\m \in \\0, 1\\\\,
\\1\\ where \\g_e \> 0\\ (or \\g_e \ge 0\\ for a non-strict comparison);
\\a \< b\\ becomes the mode of \\g_e = b - a\\. All generated code reads
the mode in place of the comparison, so between switches \\f(x, t; m)\\
is as smooth as the branches. \\\mathrm{Heaviside}(g_e)\\ and
\\\mathrm{sign}(g_e)\\ take two modes, on \\g_e\\ and \\-g_e\\, to keep
the value at \\g_e = 0\\. Equality and inequality stay as written: they
switch on a set of measure zero.

A switch is a root event with condition \\g_e\\ and the identity as
reset. The crossing step stays on the branch it started on. It is
localized to a few doubles in \\t\\ rather than to the root tolerance,
since the switching time enters the state through the jump of \\f\\. The
mode is set between forward and backward shift, so [(2)](#eq-fixed-heun)
uses \\f^{-} = f(\\\cdot\\; m^{-})\\ and \\f^{+} = f(\\\cdot\\;
m^{+})\\. To first order this is the [saltation matrix](#saltation) with
\\\partial R/\partial x = I\\, \\ S\_\theta^{+} = S\_\theta^{-} +
(f^{-} - f^{+})\\ \frac{\partial t^{\*}}{\partial \theta}, \qquad
\frac{\partial t^{\*}}{\partial \theta} = -\frac{\nabla_x g_e \cdot
S\_\theta^{-} + \partial g_e/\partial \theta} {\nabla_x g_e \cdot
f^{-} + \partial g_e/\partial t}, \\ and second order follows from the
same dual arithmetic, with the \\G\_{tt}\\ correction
[(6)](#eq-dt-star-second) along \\f^{-}\\. At a kink the first-order
jump vanishes, and the transport supplies the second-order term from the
jump of \\\partial f/\partial x\\.

The modes are read off the state at the start of the solve and after
every reset, which can move the state across a switch; after a crossing
the mode takes the side the root crossed to. A switch adds no output row
and does not count towards `maxroot`. A state leaving a root zero at the
start of a step, to the side the mode does not describe, shows no sign
change; the departure counts as a crossing once \\\|g_e\|\\ exceeds root
residual and rounding.

### Sliding along a switching surface

On the surface a switch takes the side both branches move to, by the
sign of \\\dot g_e = \partial g_e/\partial t + \nabla_x g_e \cdot f\\
under each mode. For \\\dot g\_{e,\mathrm{on}} \< 0 \< \dot
g\_{e,\mathrm{off}}\\ the solution slides
[\[1\]](#ref-filippov1988differential) on \\ f = a\\ f\_{\mathrm{on}} +
(1 - a)\\ f\_{\mathrm{off}}, \qquad a = \frac{\dot
g\_{e,\mathrm{off}}}{\dot g\_{e,\mathrm{off}} - \dot
g\_{e,\mathrm{on}}}, \\ with \\\dot g_e = 0\\. A slide ends where \\a\\
reaches \\0\\ or \\1\\, another root, and the solution leaves to that
side with continuous \\f\\. Switches sliding at once take \\\sum_c
w_c(a)\\ f_c\\ over the mode combinations \\c\\, with \\w_c\\ the
product of \\a_k\\ or \\1 - a_k\\ over the sliding switches \\k\\ and
\\a\\ chosen to zero every sliding rate. The solve stops with an error
without such a mix, or when a switch leaves the surface twice in a row
without sliding. The weights are formed in the scalar type of the model,
so the dual sensitivities keep the dependence on \\x\\ and \\\theta\\.

### Switches in the reverse mode and in `cvode()`

The reverse mode stores the modes of each step with the checkpoint, sets
them before the adjoint of the step, and transposes the switch as a root
event ([the adjoint of an
event](https://dmodverse.github.io/cppDE/articles/methods-reverse-events.html#event-adjoint)),
with \\f^{-}\\ and the derivatives of \\f^{-}\\ in the old mode and
\\f^{+}\\ in the new one. While sliding, the weights depend on the state
through the rates \\D_k(a) = \sum_c w_c(a)\\ \dot g\_{k,c}\\ held at
zero, \\\nabla a = -D_a^{-1} D_x\\, and the pullback of \\\lambda\\ is
\\ J^{\mathsf T}\lambda = \sum_c w_c\\ J_c^{\mathsf T}\lambda - \sum_k
u_k \sum_c w_c\\ \nabla \dot g\_{k,c}, \qquad u = D_a^{-\mathsf T} b,
\quad b_j = \frac{\partial f}{\partial a_j}\cdot\lambda, \\ with \\\dot
g\_{k,c}\\ the rate of switch \\k\\ at combination \\c\\ and \\\nabla
\dot g\_{k,c}\\ the generated gradient of that rate, likewise in the
parameters.

[`cvode()`](https://dmodverse.github.io/cppDE/reference/cvode.md) hands
\\g_e\\ to CVODES root finding and at the returned root

- sets the mode,
- adds \\(f^{-} - f^{+})\\\partial t^{\*}/\partial \theta\\ to the
  sensitivities before `CVodeSensReInit()`,
- under adjoint sensitivity analysis starts a new stretch, whose
  backward problem reads the modes of that stretch.

Slides follow the same rule, with two further roots per switch for the
slide ends.

## Saltation matrix

The classical *saltation matrix* adds to \\\partial R/\partial x\\ a
term for the \\\theta\\-dependence of the firing time: \\ S\_\theta^{+}
= \frac{\partial R}{\partial x}\\S\_\theta^{-} + \left( \frac{\partial
R}{\partial t} + \frac{\partial R}{\partial x}\\f^{-} - f^{+} \right)
\frac{\partial t^{\*}}{\partial \theta}. \\ The first-order part of the
Heun transport equals it, with \\\partial t^{\*}/\partial \theta =
\partial t_e/\partial \theta\\ for time-triggered events and
[(4)](#eq-dt-star-first) for root events; the transport also gives
second order. Implementation: `cppde_saltation.hpp`, event loop in
`cppde_event_engine.hpp`.

## References

\[1\]

A. F. Filippov, *Differential equations with discontinuous righthand
sides*. in Mathematics and its applications, no. 18. Kluwer Academic
Publishers, 1988. doi:
[10.1007/978-94-015-7793-9](https://doi.org/10.1007/978-94-015-7793-9).
