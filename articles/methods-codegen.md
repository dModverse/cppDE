# Code generation

[`cppODE()`](https://dmodverse.github.io/cppDE/reference/cppODE.md),
[`cvode()`](https://dmodverse.github.io/cppDE/reference/cvode.md) and
[`cppFUN()`](https://dmodverse.github.io/cppDE/reference/cppFUN.md)
share one code generator for equations in the syntax of
[`expressions`](https://dmodverse.github.io/cppDE/reference/expressions.md).
Symbols follow the [notation
table](https://dmodverse.github.io/cppDE/articles/methods-overview.html#notation).

## Expression graph and emitted code

The generator emits straight-line C++ templated over the scalar type,
one source for value, dual and nested-dual instantiations. Automatic
differentiation [\[1\]](#ref-griewank2008evaluating) on an expression
graph forms the Jacobian, the discrete-adjoint contractions, and the
time derivatives of right-hand side and event roots at code-generation
time.

**Expression graph.** Python’s parser reads each expression into a
directed acyclic graph. States, parameters, forcings, time and constants
are leaves; sums, products, powers and function calls are inner nodes.
Nodes are hash-consed [\[2\]](#ref-ershov1958programming),
[\[3\]](#ref-filliatre2006hashconsing): identified by operation,
operands and constant attributes, each is created once, so a repeated
subexpression is one node without an elimination pass. Sums and products
are flat with numeric coefficients and exponents; equal constants are
one node. SymPy reads syntax the parser does not and differentiates
functions without a graph rule.

![Expression graph of \$f_A = f_B = -k_1 A B + k_2 C\$ and \$f_C = k_1 A
B - k_2 C\$; edge labels are
coefficients.](methods-codegen_files/figure-html/expression-graph-1.png)

Expression graph of \\f_A = f_B = -k_1 A B + k_2 C\\ and \\f_C = k_1 A
B - k_2 C\\; edge labels are coefficients.

**Example.** Nodes used twice or more become temporaries:

``` cpp
const double _t0 = params[4]*x[2];       // k2*C
const double _t1 = params[3]*x[0]*x[1];  // k1*A*B
const double _t2 = _t0 - _t1;
dxdt[0] = _t2;
dxdt[1] = _t2;
dxdt[2] = -_t0 + _t1;
```

**Derivatives.** Forward and reverse sweeps over the graph form
Jacobian-vector and vector-Jacobian products in time linear in the graph
size [\[1\]](#ref-griewank2008evaluating), as new nodes of the same
graph, so derivative and value share common parts. The reverse
contractions \\J^{\mathsf T}\lambda\\, \\(\partial f/\partial
p)^{\mathsf T}\lambda\\ are vector-Jacobian products; Rosenbrock stages
and event localization use Jacobian-vector products along \\\dot x =
f\\.

\\J = \partial f/\partial x\\ is assembled by entries or by coloring:

- *Entries.* Each state-dependent node holds a sparse state gradient;
  those of the graph roots are the nonzero entries of \\J\\. The cost
  \\W\_{\mathrm{ent}} = \sum\_{v} \|\mathrm{supp}\\\nabla v\|\\ over the
  \\E_s\\ state-dependent nodes \\v\\ of \\f\\ grows with the density of
  \\J\\.
- *Coloring*, when \\W\_{\mathrm{ent}} \> \max(4E_s, 50\\000)\\. A
  greedy distance-2 coloring of the column intersection graph groups
  structurally orthogonal columns [\[4\]](#ref-curtis1974estimation),
  [\[5\]](#ref-coleman1983estimation),
  [\[6\]](#ref-gebremedhin2005color). With at most \\n_x/4\\ colors, one
  Jacobian-vector product per color, seeded with ones on the columns of
  that color, recovers those columns; with more colors, entries are
  kept.

**Linear tables.** A sum of twelve or more terms \\c_j\\x_j\\\phi\\,
numeric \\c_j\\, common state-free \\\phi\\, becomes a row of a constant
matrix \\C\\, \\f(x) = \psi(x, Cx)\\. The code evaluates \\z = Cx\\ once
per call from a sparse or dense table and assembles \\ J =
\frac{\partial \psi}{\partial x} + \frac{\partial \psi}{\partial z}\\C,
\tag{1} \\ while the reverse contractions apply \\C^{\mathsf T}\\. A
dense coupling, such as a dipolar interaction of all states, costs one
table, not a sum of expressions per entry.

**Loops.** Statements whose graphs agree up to state indices and
constants, such as method-of-lines grid points, form a class. Sixteen or
more instances become one loop over index and constant tables, a sum of
many like terms an inner loop. Differentiation stays per instance, so
generator work grows with the instance count; the emitted code of a
regular model does not.

## References

\[1\]

A. Griewank and A. Walther, *Evaluating derivatives: Principles and
techniques of algorithmic differentiation*, 2nd ed. in Other titles in
applied mathematics, no. 105. SIAM, 2008. doi:
[10.1137/1.9780898717761](https://doi.org/10.1137/1.9780898717761).

\[2\]

A. P. Ershov, “On programming of arithmetic operations,” *Communications
of the ACM*, vol. 1, no. 8, pp. 3–6, 1958, doi:
[10.1145/368892.368907](https://doi.org/10.1145/368892.368907).

\[3\]

J.-C. Filliâtre and S. Conchon, “Type-safe modular hash-consing,” in
*Proceedings of the 2006 workshop on ML*, ACM, 2006, pp. 12–19. doi:
[10.1145/1159876.1159880](https://doi.org/10.1145/1159876.1159880).

\[4\]

A. R. Curtis, M. J. D. Powell, and J. K. Reid, “On the estimation of
sparse Jacobian matrices,” *IMA Journal of Applied Mathematics*, vol.
13, no. 1, pp. 117–119, 1974, doi:
[10.1093/imamat/13.1.117](https://doi.org/10.1093/imamat/13.1.117).

\[5\]

T. F. Coleman and J. J. Moré, “Estimation of sparse Jacobian matrices
and graph coloring problems,” *SIAM Journal on Numerical Analysis*, vol.
20, no. 1, pp. 187–209, 1983, doi:
[10.1137/0720013](https://doi.org/10.1137/0720013).

\[6\]

A. H. Gebremedhin, F. Manne, and A. Pothen, “What color is your
Jacobian? Graph coloring for computing derivatives,” *SIAM Review*, vol.
47, no. 4, pp. 629–705, 2005, doi:
[10.1137/S0036144504444711](https://doi.org/10.1137/S0036144504444711).
