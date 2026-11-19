# Method coefficients

Coefficients and sources of the four methods of [Stepping and error
control](https://dmodverse.github.io/cppDE/articles/methods-stepping.md),
selected by `method` and `useNDF` of
[`cppODE()`](https://dmodverse.github.io/cppDE/reference/cppODE.md).
Symbols follow the [notation
table](https://dmodverse.github.io/cppDE/articles/methods-overview.html#notation).

## BDF and NDF coefficients

*Source.* Shampine and Reichelt [\[1\]](#ref-MATLABODESuite), Table 1;
implemented in `cppde_multistepper.hpp`.

The NDF of order \\q\\ adds to the BDF\\q\\ corrector a term
proportional to the Nordsieck predictor error, \\ \sum\_{m=1}^{q}
\frac{1}{m}\\ \nabla^{m} x\_{n+1} \\-\\ h_n\\ f(x\_{n+1}, t\_{n+1})
\\-\\ \kappa_q\\ H_q\\ \bigl(x\_{n+1} - x\_{n+1}^{\mathrm{pred}}\bigr)
\\=\\ 0, \qquad H_q \\=\\ \sum\_{m=1}^{q} \tfrac{1}{m}, \tag{1} \\ with
\\\nabla\\ the backward-difference operator; \\\kappa_q = 0\\ gives
BDF\\q\\. `useNDF` (default `TRUE`) selects NDF or BDF at compile time.

| \\q\\ | \\\kappa_q\\ | step-size gain over BDF\\q\\ | stability angle           |
|:------|-------------:|:-----------------------------|:--------------------------|
| 1     |      -0.1850 | \\+26\\\\\\                  | A-stable                  |
| 2     |      -0.1111 | \\+26\\\\\\                  | A-stable                  |
| 3     |      -0.0823 | \\+26\\\\\\                  | \\80^\circ\\              |
| 4     |      -0.0415 | \\+12\\\\\\                  | \\66^\circ\\              |
| 5     |       0.0000 | \\0\\\\\\                    | \\51^\circ\\ (\\=\\ BDF5) |

NDF coefficients \\\kappa_q\\ of cppDE, empirical step-size gain over
BDF\\q\\ at matched local error, and stability angle. NDF5 is BDF5,
since the BDF5 angle of \\51^\circ\\ leaves no room for reduction.
{.table}

## Adams-Moulton coefficients

*Source.* Hindmarsh et al. [\[2\]](#ref-hindmarsh2005sundials), the
\\\ell\\-vector formulation of CVODE; implemented in
`cppde_multistepper.hpp` by `adams_set_coefficients()`.

The coefficients depend on the past step sizes and are recomputed every
step. With \\\xi_j = (t\_{n+1} - t\_{n+1-j})/h_n\\, the correction
vector \\(\ell_0, \dots, \ell_q)\\ of the Nordsieck array holds the
coefficients of the polynomial \\\Lambda(u) = \sum\_{j=0}^{q} \ell_j
u^j\\ with \\ \Lambda(-1) \\=\\ 0, \qquad \Lambda(0) \\=\\ 1, \qquad
\Lambda'(u) \\\propto\\ \prod\_{j=1}^{q-1} \bigl(1 + u/\xi_j\bigr), \\
and the weight of the implicit term is \\\beta_0 = 1/\ell_1\\. At
constant step size the corrector weights \\\beta_0, \dots, \beta_q\\ are
the fixed-step values of [\[3\]](#ref-hairer2008solving), Table III.1.1.
The controller selects \\q \in \\1, \dots, 12\\\\; \\q\_{\max} = 12\\ is
the CVODE default.

## Rosenbrock4 coefficients

*Source.* Hairer and Wanner [\[4, Sec. IV.7\]](#ref-hairer2010solving),
the default coefficients of the `rosenbrock4` stepper of
Boost.Numeric.Odeint except the sign of \\d_4\\, which is that of
`rodas.f`; implemented in `default_rosenbrock_coefficients` in
`cppde_rosenbrock4.hpp`.

With \\W = (\gamma h_n)^{-1} I - J(x_n, t_n)\\, the stage equations are
\\ \begin{aligned} W\\ k_i &\\=\\ f\\\left(x_n + \sum\_{j\<i}
\alpha\_{ij}\\ k_j,\\ t_n + c_i\\ h_n \right) \\+\\ \frac{1}{h_n}
\sum\_{j\<i} \gamma\_{ij}\\ k_j \\+\\ h_n\\ d_i\\ \frac{\partial
f}{\partial t}, \\ x\_{n+1} &\\=\\ x_n + \sum\_{i=1}^{s} m_i\\ k_i,
\qquad \hat x\_{n+1} \\=\\ x_n + \sum\_{i=1}^{s} \hat m_i\\ k_i.
\end{aligned} \\ cppDE uses \\s = 6\\ and \\\gamma = 1/4\\, indexed as
in the source.

|            |             |         |        |         |       |
|:-----------|------------:|--------:|-------:|--------:|------:|
| \\\gamma\\ |        0.25 | \\d_1\\ |   0.25 | \\c_2\\ | 0.386 |
| \\d_2\\    | \\-0.1043\\ | \\d_3\\ | 0.1035 | \\c_3\\ |  0.21 |
| \\d_4\\    | \\-0.0362\\ |         |        | \\c_4\\ |  0.63 |

Rosenbrock4 diagonal coefficient \\\gamma\\, time-partial coefficients
\\d_i\\, nodes \\c_i\\. {.table}

|                   |      \\j=1\\ |      \\j=2\\ |       \\j=3\\ |       \\j=4\\ |
|:------------------|-------------:|-------------:|--------------:|--------------:|
| \\\alpha\_{2,j}\\ | 1.5440000000 |              |               |               |
| \\\alpha\_{3,j}\\ | 0.9466785281 | 0.2557011699 |               |               |
| \\\alpha\_{4,j}\\ | 3.3148251871 | 2.8961240160 |  0.9986419140 |               |
| \\\alpha\_{5,j}\\ | 1.2212245092 | 6.0191344813 | 12.5370833293 | -0.6878860361 |

Rosenbrock4 stage couplings \\\alpha\_{ij}\\: row \\i\\ the stage
formed, column \\j\\ an earlier stage. {.table}

|                   |       \\j=1\\ |        \\j=2\\ |        \\j=3\\ |       \\j=4\\ |
|:------------------|--------------:|---------------:|---------------:|--------------:|
| \\\gamma\_{2,j}\\ | -5.6688000000 |                |                |               |
| \\\gamma\_{3,j}\\ | -2.4300933568 |  -0.2063599157 |                |               |
| \\\gamma\_{4,j}\\ | -0.1073529058 |  -9.5945622510 | -20.4702861481 |               |
| \\\gamma\_{5,j}\\ |  7.4964433140 | -10.2468043146 | -33.9999035282 | 11.7089089321 |

Rosenbrock4 stage couplings \\\gamma\_{ij}\\ of the term
\\(1/h_n)\sum\_{j\<i}\gamma\_{ij}\\k_j\\, indexed as \\\alpha\_{ij}\\.
{.table}

The sixth stage, at the embedded order-\\3\\ solution \\X_6\\, is the
error estimate \\e\_{n+1}\\, which the step solution adds to \\X_6\\: \\
X_6 \\=\\ x_n + \sum\_{i=1}^{5} \alpha\_{6,i}\\ k_i, \qquad W\\ e\_{n+1}
\\=\\ f(X_6,\\ t_n + h_n) + \frac{1}{h_n}\sum\_{i=1}^{5} \gamma\_{6,i}\\
k_i, \qquad x\_{n+1} \\=\\ X_6 + e\_{n+1}, \\ with the sixth-row
coefficients

|  | \\i=1\\ | \\i=2\\ | \\i=3\\ | \\i=4\\ | \\i=5\\ |
|:---|---:|---:|---:|---:|---:|
| \\\alpha\_{6,i}\\ | 1.2212245092 | 6.0191344813 | 12.5370833293 | -0.6878860361 | 1.0 |
| \\\gamma\_{6,i}\\ | 8.0832467959 | -7.9811329881 | -31.5215943287 | 16.3193054312 | -6.0588182388 |

Rosenbrock4 sixth-row coefficients: \\\alpha\_{6,i}\\ form \\X_6\\,
\\\gamma\_{6,i}\\ enter the error estimate \\e\_{n+1}\\. {.table
style="width:100%;"}

The order-\\3\\ interpolant on \\\[t_n, t\_{n+1}\]\\ combines the five
stages and the step solution.

|  | \\i=1\\ | \\i=2\\ | \\i=3\\ | \\i=4\\ | \\i=5\\ |
|:---|---:|---:|---:|---:|---:|
| \\d_i^{(2)}\\ | 10.1262350834 | -7.4879958776 | -34.8009186156 | -7.9927717076 | 1.0251377233 |
| \\d_i^{(3)}\\ | -0.6762803393 | 6.0877146517 | 16.4308432089 | 24.7672251142 | -6.5943891257 |

Rosenbrock4 dense-output coefficients \\d_i^{(2)}\\ and \\d_i^{(3)}\\.
{.table style="width:100%;"}

## Tsitouras 5(4) coefficients

*Source.* Tsitouras [\[5\]](#ref-TSITOURAS2011770), Table 1 and eqs.
(7), (8); implemented in `cppde_tsit5.hpp`.

The pair has \\s = 7\\ stages and satisfies the first column simplifying
assumption \\\sum_i b_i a\_{ij} = b_j(1 - c_j)\\ for every \\j\\, but
not the row condition \\\sum_j a\_{ij} c_j = c_i^2/2\\. The seventh
stage is the FSAL stage: row \\a\_{7,j}\\ equals \\b_j\\ for \\j = 1,
\dots, 6\\.

| \\i\\ |      \\c_i\\ |   \\a\_{i1}\\ |    \\a\_{i2}\\ |  \\a\_{i3}\\ |
|:------|-------------:|--------------:|---------------:|-------------:|
| 1     |            0 |               |                |              |
| 2     |        0.161 |         0.161 |                |              |
| 3     |        0.327 | -0.0084806555 |   0.3354806555 |              |
| 4     |          0.9 |  2.8971530571 |  -6.3594484900 | 4.3622954329 |
| 5     | 0.9800255409 |  5.3258648284 | -11.7488835641 | 7.4955393429 |
| 6     |          1.0 |  5.8614554429 | -12.9209693178 | 8.1593678986 |
| 7     |          1.0 |  0.0964607668 |           0.01 | 0.4798896504 |

Tsitouras 5(4) nodes \\c_i\\ and coefficients \\a\_{i1}\\ to
\\a\_{i3}\\. {.table}

| \\i\\ |   \\a\_{i4}\\ |   \\a\_{i5}\\ |  \\a\_{i6}\\ |
|:------|--------------:|--------------:|-------------:|
| 5     | -0.0924950664 |               |              |
| 6     | -0.0715849733 | -0.0282690504 |              |
| 7     |  1.3790085741 | -3.2900695154 | 2.3247105241 |

Tsitouras 5(4) coefficients \\a\_{i4}\\ to \\a\_{i6}\\, nonzero from row
\\5\\ on. Row \\7\\ of both tables is \\b_j\\ (FSAL). {.table}

\\k_7 = f(x\_{n+1},\\ t_n + h_n)\\ is \\k_1\\ of step \\n+1\\. With
\\b_7 = 0\\ and \\\hat b_7 \neq 0\\, it enters the order-\\4\\ estimate
only. The controller uses \\ e\_{n+1} \\=\\ h_n \sum\_{i=1}^{7} (\hat
b_i - b_i)\\ k_i. \\

| \\i\\ |       \\b_i\\ |  \\\hat b_i\\ | \\\hat b_i - b_i\\ |
|:------|--------------:|--------------:|-------------------:|
| 1     |  0.0964607668 |  0.0946807557 |      -0.0017800111 |
| 2     |          0.01 |  0.0091835655 |      -0.0008164345 |
| 3     |  0.4798896504 |  0.4877705284 |       0.0078808780 |
| 4     |  1.3790085741 |  1.2342975669 |      -0.1447110072 |
| 5     | -3.2900695154 | -2.7077123499 |       0.5823571655 |
| 6     |  2.3247105241 |  1.8666284182 |      -0.4580821059 |
| 7     |             0 |  0.0151515152 |       0.0151515152 |

Tsitouras 5(4) weights of order \\5\\ (\\b_i\\) and \\4\\ (\\\hat
b_i\\). The differences, `e1`, …, `e7` in `cppde_tsit5.hpp`, sum to
zero. {.table}

## References

\[1\]

L. F. Shampine and M. W. Reichelt, “The MATLAB ODE suite,” *SIAM Journal
on Scientific Computing*, vol. 18, no. 1, pp. 1–22, 1997, doi:
[10.1137/S1064827594276424](https://doi.org/10.1137/S1064827594276424).

\[2\]

A. C. Hindmarsh *et al.*, “SUNDIALS: Suite of nonlinear and
differential/algebraic equation solvers,” *ACM Transactions on
Mathematical Software (TOMS)*, vol. 31, no. 3, pp. 363–396, 2005, doi:
[10.1145/1089014.1089020](https://doi.org/10.1145/1089014.1089020).

\[3\]

E. Hairer, S. P. Nørsett, and G. Wanner, *Solving ordinary differential
equations I: Nonstiff problems*, 2nd ed. in Springer series in
computational mathematics, no. 8. Springer, 1993. doi:
[10.1007/978-3-540-78862-1](https://doi.org/10.1007/978-3-540-78862-1).

\[4\]

E. Hairer and G. Wanner, *Solving ordinary differential equations II:
Stiff and differential-algebraic problems*, 2nd ed. in Springer series
in computational mathematics, no. 14. Springer, 1996. doi:
[10.1007/978-3-642-05221-7](https://doi.org/10.1007/978-3-642-05221-7).

\[5\]

Ch. Tsitouras, “Runge–Kutta pairs of order 5(4) satisfying only the
first column simplifying assumption,” *Computers & Mathematics with
Applications*, vol. 62, no. 2, pp. 770–775, 2011, doi:
[10.1016/j.camwa.2011.06.002](https://doi.org/10.1016/j.camwa.2011.06.002).
