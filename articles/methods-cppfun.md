# Derivative directions of cppFUN()

Here \\x\\ denotes the variables of the function, not ODE states; other
symbols are those of the [notation
table](https://dmodverse.github.io/cppDE/articles/methods-overview.html#notation).
R interface:
[`cppFUN()`](https://dmodverse.github.io/cppDE/reference/cppFUN.md) and
[`evaluateMany()`](https://dmodverse.github.io/cppDE/reference/evaluateMany.md).

[`cppFUN()`](https://dmodverse.github.io/cppDE/reference/cppFUN.md)
compiles an algebraic system \\ y_i = g_i(x_1, \dots, x\_{n_x}, p_1,
\dots, p\_{n_p}), \qquad i = 1, \dots, n_y, \tag{1} \\ with the first
and second derivatives of \\g\\, named as the derivatives of
[`solveODE()`](https://dmodverse.github.io/cppDE/reference/solveODE.md).
The variables \\x\\ are the rows of an \\n\_{\mathrm{obs}} \times n_x\\
matrix, one per observation; the parameters \\p\\ are shared by all
rows. Typical uses: observation maps and reparameterizations \\(x_0, p)
= \Phi(\theta)\\, the Jacobian being the tangent of [Forward
(2)](https://dmodverse.github.io/cppDE/articles/methods-forward.html#eq-sensitivity-theta).

From the [expression
graph](https://dmodverse.github.io/cppDE/articles/methods-codegen.html#codegen)
the generator emits two templates over the scalar type:

- \\g\\ for one observation, instantiated on `double`,
  `cppde::dual<double>` and the nested `cppde::dual2nd<double>`;
- the vector-Jacobian product as a reverse sweep on the same graph,
  instantiated on `double` and the dual.

With the default `compile = FALSE` the constructor only writes the
source, and the functions stop with an error until
[`compile()`](https://dmodverse.github.io/cppDE/reference/compile.md).

## Entry points and derivative directions

`derivMode` names the products built, any of `"forward"` (default),
`"reverse"` and `"forward-reverse"` (implies `"reverse"`); a product
left out costs no entries and no compile time. `deriv` and `deriv2`
switch the first- and second-order entries; `deriv2 = TRUE` requires
`"forward"`.

| R function | direction | compiled entry |
|:---|:---|:---|
| `func(x, p)` | any | `_eval_c` |
| `jac(x, p[, tangentX, tangentP])` | forward | `_eval_ad_c` |
| `hess(x, p[, tangentX, tangentP, hessianX, hessianP])` | forward | `_eval_ad2_c` |
| `evaluate(..., deriv2)` | forward | `_eval_ad_c` or `_eval_ad2_c` |
| `evaluateBatch(sets, cores, deriv2)` | forward | `_eval_ad_batch` or `_eval_ad2_batch` |
| `vjp(x, p, cotangent)` | reverse | `_vjp_c` |
| `vjp(x, p, cotangent, tangentX, tangentP, curvature)` | forward-reverse | `_vjp_ad_c` |

R functions and compiled entries, prefixed by the model name; `x`, `p`
stand for `vars`, `params`. `_c` entries are `.Call` interfaces taking
arguments by reference; `.C` versions of `_eval`, `_eval_ad`,
`_eval_ad2` and `_vjp` serve a shared object without the `.Call` symbol.
Under `"forward-reverse"`, `_vjp_ad` is also a thread-safe plain C
entry. {.table}

**Forward direction.** `jac` evaluates \\g\\ on duals, `hess` on nested
duals; they return arrays \\\[n\_{\mathrm{obs}}, n_y, M\]\\ and
\\\[n\_{\mathrm{obs}}, n_y, M, M\]\\. Tangents (next section) set the
directions \\\theta \in \mathbb{R}^M\\, by default the identity over the
variables and differentiated parameters. `evaluate` returns `y`,
`tangent` and, with `deriv2 = TRUE`, `hessian` from one pass instead of
up to three. `evaluateBatch` runs `evaluate` over a list of argument
sets on OpenMP threads with dynamic scheduling, as a loop for fewer than
two sets.
[`evaluateMany()`](https://dmodverse.github.io/cppDE/reference/evaluateMany.md)
does the same for one argument set each of several
[`cppFUN()`](https://dmodverse.github.io/cppDE/reference/cppFUN.md)
objects, in one parallel region.

**Reverse direction.** For a cotangent \\w\\ of the outputs,
\\\[n\_{\mathrm{obs}}, n_y\]\\ or \\\[n\_{\mathrm{obs}}, n_y,
n\_{\mathrm{seed}}\]\\, `vjp` returns the value \\y\\ and \\ (w_x)\_{o j
s} = \sum\_{i} w\_{o i s}\\ \frac{\partial g_i}{\partial x_j}(x_o, p),
\qquad (w_p)\_{k s} = \sum\_{o} \sum\_{i} w\_{o i s}\\ \frac{\partial
g_i}{\partial p_k}(x_o, p), \tag{2} \\ as `cotangentX` and `cotangentP`.
The contraction is one reverse sweep per observation and cotangent
column, at a cost independent of the number of upstream parameters,
without forming Jacobian or Hessian. Statements independent of the
cotangent run once per observation.

With input tangents \\v_x\\, \\v_p\\ in \\n\_{\mathrm{dir}}\\ directions
(`tangentX`, `tangentP`) or the derivative \\\dot w\\ of the cotangent
along them (`curvature`), `vjp` runs the same template on duals and also
returns \\ \frac{\mathrm{d}}{\mathrm{d}v}\bigl(w^{\mathsf T} J\bigr) =
w^{\mathsf T}\\\frac{\partial J}{\partial (x, p)}\\v + \dot w^{\mathsf
T} J, \qquad J = \frac{\partial g}{\partial (x, p)}, \tag{3} \\ as
`curvatureX` and `curvatureP`, forward over reverse ([second
derivatives](https://dmodverse.github.io/cppDE/articles/methods-forward-reverse.md))
without a Hessian.

## Tangent, hessian, cotangent and curvature

The upstream map \\\Phi: \theta \mapsto (x, p)\\ has the `tangent`
\\\partial\Phi(\theta)\\ and the `hessian` \\\partial^2\Phi(\theta)\\. A
scalar functional \\L\\ of the outputs has the `cotangent` \\\partial
L(y)\\ and the `curvature` \\\partial^2 L(y)\\\partial y\\, the
derivative of the cotangent along the tangent. Forward arguments live on
the inputs, reverse arguments on the outputs, results on the other side.

| [`solveODE()`](https://dmodverse.github.io/cppDE/reference/solveODE.md) | [`cppFUN()`](https://dmodverse.github.io/cppDE/reference/cppFUN.md) argument | [`cppFUN()`](https://dmodverse.github.io/cppDE/reference/cppFUN.md) result |
|:---|:---|:---|
| `tangent` | `tangentX`, `tangentP` | `tangent` |
| `hessian` | `hessianX`, `hessianP` | `hessian` |
| `cotangent` | `cotangent` | `cotangentX`, `cotangentP` |
| `curvature` | `curvature` | `curvatureX`, `curvatureP` |

The four derivative objects in R.
[`solveODE()`](https://dmodverse.github.io/cppDE/reference/solveODE.md)
uses one name for argument and result;
[`cppFUN()`](https://dmodverse.github.io/cppDE/reference/cppFUN.md)
results are elements of the lists from `evaluate` and `vjp`. {.table}

Forward over forward maps input to output `hessian` in every pair of
tangent directions. Forward over reverse maps output to input
`curvature`, the Hessian of the seeded functional applied to the
tangent; the identity tangent gives the Hessian columns.

**Chain rule.** Tangent and hessian of the inputs come in two blocks
each:

- `tangentX[obs, name, theta]`: \\\partial x_i/\partial\theta_a\\, per
  observation row;
- `tangentP[name, theta]`: \\\partial p_j/\partial\theta_a\\, shared by
  all observations;
- `hessianX[obs, name, theta1, theta2]`: \\\partial^2
  x_i/\partial\theta_a\partial\theta_b\\, optional, zero by default;
- `hessianP[name, theta1, theta2]`: \\\partial^2
  p_j/\partial\theta_a\partial\theta_b\\, optional, zero by default.

Attributes `deriv` and `deriv2` of variables and parameters are an
alternative. Parameters in `fixed` and symbols absent from the tangent
have zero tangent rows.

With the tangent \\\partial\Phi \in \mathbb{R}^{n\_{\mathrm{obs}} \times
n\_{\mathrm{diff}} \times M}\\ over the \\n\_{\mathrm{diff}} = n_x +
n_p\\ inputs and the hessian \\\partial^2\Phi\\, the chain rule gives
the Jacobian \\ \frac{\partial y}{\partial\theta}\bigg\|\_{\mathrm{obs}}
= J\bigl(x\_{\mathrm{obs}}, p\bigr) \cdot \partial\Phi\_{\mathrm{obs}},
\tag{4} \\ and the Hessian of output \\i\\ \\ \frac{\partial^2
y_i}{\partial\theta\\\partial\theta^{\mathsf T}}\bigg\|\_{\mathrm{obs}}
= \partial\Phi\_{\mathrm{obs}}^{\mathsf T} \\
H_i\bigl(x\_{\mathrm{obs}}, p\bigr) \\ \partial\Phi\_{\mathrm{obs}} +
\sum\_{j} J\_{i, j}\bigl(x\_{\mathrm{obs}}, p\bigr) \\
\partial^2\Phi\_{\mathrm{obs}, j, \cdot, \cdot}, \tag{5} \\ whose second
term vanishes for affine \\\Phi\\. Neither \\J\\ nor \\H\\ is formed;
seeding the dual tangent slots gives both contractions. Each input
becomes \\ \tilde x_i = (x_i +
\partial\Phi\_{i,\cdot}\\\varepsilon\_{\mathrm{in}}) +
\bigl(\partial\Phi\_{i,\cdot} +
\partial^2\Phi\_{i,\cdot,\cdot}\\\varepsilon\_{\mathrm{in}}\bigr)\\\varepsilon\_{\mathrm{out}},
\\ a nested dual: inner tangents \\\partial x_i/\partial\theta\\, outer
slot \\k\\ an inner dual with value \\\partial x_i/\partial\theta_k\\
and inner tangents \\\partial^2 x_i/\partial\theta_k\partial\theta\\.
One evaluation propagates both terms of (5). The inner tangents of
\\\tilde y\\ are \\\partial y/\partial\theta\\, the cross slots \\\tilde
y\[k\]\[m\]\\ are \\\partial^2 y/\partial\theta_k\partial\theta_m\\.

## Pass-through of unmodeled inputs

With `attach.input = TRUE`, inputs that are not variables or parameters
of [(1)](#eq-funcpp) pass to the output unchanged: values copied,
Jacobian rows zero for a variable and the identity for a parameter,
Hessian rows zero. An observation map then need not repeat every
integrator column as an identity, as with the partial-shape `tangent` of
[`solveODE()`](https://dmodverse.github.io/cppDE/reference/solveODE.md).
