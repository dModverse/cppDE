# cppDE 0.13.16

* A product with a constant and a quotient by a constant have finite first and second derivatives
  where the other operand is infinite: `exp10(x)` at `x = -Inf` gives 0, not NaN.

# cppDE 0.13.15

* The reverse and forward-reverse modes of `bdf` solve the corrector equation of every step in the
  backward sweep, to a tenth of the local error test of the backward problem against `reltol` and
  `abstol`, tested on the values: against the iteration matrix of the forward run and, where that
  cannot reach the test within ten iterations, against the iteration matrix at the accepted state.
  The gradient is the derivative of the method with the corrector equations solved, and so is the
  Hessian of forward-reverse. The forward sensitivities with `sensErrCon = FALSE` take the corrector
  iterations of the value run and differ from it by what those iterations leave unconverged. Under
  `adams` the sweep transposes the iterations taken: the exactly solved Adams-Moulton step can be
  unstable on a stiff model at the step sizes of the value run.
* A forward solve with `sensErrCon = FALSE` takes the steps of a value run up to rounding; on a long
  stiff run the two step sequences can part.
* Under `adjointControl(refine = TRUE)` a `bdf` or `adams` sweep reuses the factorizations of the
  extrapolated flow over intervals of equal length and sweeps all seed columns at once.
* A root event at the last output time ends the solve on the state the event leaves.
* The event table of `cppODE()` and `cvode()` has the columns `var`, `value`, `method` and `time`
  or `root`; any other column is an error that names it. The columns `terminal` and `direction`
  are removed: a root event fires on crossings in both directions and the solve continues, and
  `rootfunc` ends a solve.
* `compile()` builds the objects of an `output` with the makefiles of `R CMD SHLIB`, so a user
  Makevars (`R_MAKEVARS_USER` or `~/.R/Makevars`) reaches the compile of every object as well as
  the link, and compiles every object anew. A source that fails to compile is an error that names
  it and shows the compiler output, also under `cores > 1`.

# cppDE 0.13.14

* The tests that call `cvode()` skip without SUNDIALS also where KLU is installed, and the states
  of a forward and a reverse solve on one grid are compared to rounding, not bit for bit.

# cppDE 0.13.13

* `compile()` applies the `compileArgs` of a model to the source of that model only, also when
  `output` links several models into one library; the step trace of one model and the KLU settings
  of a sparse one reach no other model.
* `solveBatch()` keeps the prepared parameters of a condition whose element of `parms` is `NULL`.
* `cppFUN()` functions with `attach.input = TRUE` append an extra input with one value per row as
  a column.

# cppDE 0.13.12

* The reverse mode of `bdf` and `adams` transposes the corrector iterations the forward run took,
  against the iteration matrix it factorized: the gradient equals the forward sensitivities on the
  same grid (`sensErrCon = FALSE`) to rounding, under `adams` on stiff models as well.
* The forward sensitivities of `bdf` and `adams` pass the corrector against the iteration matrix
  without the derivative of that matrix, at first and second order.
* Under `adjointControl(refine = TRUE)` a reverse solve also returns `$correction`, per seed column
  the first-order change of the seeded functional from the computed trajectory to the one the
  refined sweep differentiated, so value plus correction and the gradient belong to one function.
  A model compiled against earlier headers is refused under `refine` until it is compiled again.
* `solveBatch()` takes a `cotangent` with another number of seed columns than the prepared one.

# cppDE 0.13.11

* The backward sweep of a second-order reverse solve with a multistep method sweeps all
  cotangent columns at once, step by step: the Jacobian and factorization of a step serve every
  column.
* A reverse solve that integrates for a sweep in the same call keeps only the part of the multistep history
  that the sweep reads: the initial state; under `adjointControl(refine = TRUE)` the state at the start of
  every step; under `adjointControl(trace = TRUE)` the whole history. A store kept with
  `keepStore = TRUE` for a later solve keeps the whole history.
* An arena of dual numbers goes back to the pool with at most 128 MB of storage, so that a solve
  with many tangents returns the rest of the memory to the system.

# cppDE 0.13.10

* Every derivative the solvers use is exact; no difference quotient is left:
  * The first step of every method comes from the second time derivative of the solution,
    `J f + df/dt`, generated with the model. A multistep method that restarts at order 1 after
    repeated error test failures takes the step this curvature allows.
  * The second time derivative of a root condition along the solution, which the second order
    of an event time reads, is generated for models with forcings too, from the second derivative
    of the interpolant of the forcing.
  * Where switches slide, the derivatives of the Filippov weights come from the implicit function
    theorem on the generated gradients of the switches' rates: the Jacobian of the slide and the
    time derivative of that Jacobian, the derivatives of the exit roots, the terms of the second-order adjoint and
    the sensitivities and adjoint of `cvode()`.

* The dual types have a run-time tangent width only. Sources generated by an earlier cppDE are
  generated again before they compile against these headers.

# cppDE 0.13.9

* New `releaseStore()` frees the checkpoints of a reverse solve kept with `keepStore = TRUE` at
  once instead of at the next garbage collection; a solve given a released store integrates
  again.

* Under `derivMode = "forward-reverse"` the tangents enter the error test and the first step as
  in a forward solve, by `sensErrCon` (default `TRUE`); `sensErrCon = FALSE` keeps the step grid
  independent of the directions of the tangent. A model generated by an earlier cppDE ignores
  `sensErrCon` under forward-reverse until it is generated and compiled again.

* `install_libs()` is removed; the function is `installLibs()`.

* cppDE has a site, <https://dmodverse.github.io/cppDE/>: the reference, Get started, a Guide
  in eight chapters, the methods in ten articles, two worked examples and the benchmarks against
  SUNDIALS CVODES. The methods articles replace the PDF vignette `Methods`, and `R.rsp` is not
  suggested. The demos are articles of the site: `basics` is Get started, `saltation` and
  `fitCoPt` are the worked examples; `demo()` lists none.

# cppDE 0.13.8

* The code generators read no environment variables; `CPPDE_PY_DIR`, `CPPDE_JAC`, `CPPDE_LINEAR`
  and `CPPDE_VECTORISE` are removed.

# cppDE 0.13.7

* The dense LU factorization of the step path, the solves with it and the products of the corrector and
  Nordsieck updates run in plain loops for systems up to 64 states, outside the lock under which
  OpenBLAS hands out work buffers.

* Batches, `evaluateBatch()` and `evaluateMany()` run on the cores asked for, at most the
  processors, also with fewer conditions or sets than cores, so that OpenMP keeps the thread
  team from one call to the next; `attr(, "threads")` of a batch reports that count.

* New `evaluateMany()` evaluates argument sets, one per `cppFUN()` object, in one
  parallel region.

* `cppFUN()` with `derivMode = "forward-reverse"` also emits `<model>_vjp_ad`, a plain C entry
  of forward over reverse that native code can call from several threads at once.

# cppDE 0.13.6

* The arena that holds the tangents of temporaries and the scratch buffers of the steppers goes
  back to a pool when the outermost scope of a solve closes, with at most 1 GB of storage,
  and the next solve on any thread takes it from there. The memory of a process stays bounded
  when OpenMP ends and starts threads.

* Solves with second-order sensitivities release the temporaries of every step, so that their
  memory does not grow with the steps. An allocation that fails in the arena is an R error.

* The first step of a forward solve is estimated from the tangents the right-hand side returns;
  the step grid of forward solves with sensitivities changes within the tolerance.

* `cppde::dual2nd` has compound assignments with a `dual2nd` operand that keep the Hessian.
  Copies, moves and the results of LU solves of a `dual2nd` write the value layer only where it is
  bound.

# cppDE 0.13.5

* `ptc()` ends with Newton steps until every solved variable is within
  `reltol * (abs(x) + abstol)` of the root, so that a variable many orders below the others is
  resolved relative to itself and the root within one basin does not depend on the start.

* With `flow = TRUE`, `ptc()` returns a stable steady state: at an unstable root it continues
  along the unstable direction; `controls = list(stability = FALSE)` accepts any root of the
  flow. `controls = list(stabilityTol = )` counts a root as stable while the real part of the
  leading Jacobian eigenvalue on the kernel of `C` is at most `stabilityTol`, default 0.

* The solver of `ptc()` is available to the compiled code of other packages through
  `R_GetCCallable("cppDE", "cppde_ptc_v1")`, safe in threads; the signature is in `src/ptc.cpp`.

* `ptc()` takes the tolerances as `reltol` and `abstol` only.

# cppDE 0.13.4

* `solveODE()` with `deriv2 = TRUE` and a `tangent` but no `hessian` takes the
  second derivative of the inputs as zero.

# cppDE 0.13.3

* `compile()` builds in parallel on Windows too: separate libraries side by side on
  `cores` workers, a combined `output` with `make -j`.

* Models of the Sundials backend with events compile on Windows.

# cppDE 0.13.2

* Breaking: the generated solve entries take one named argument list. Models
  generated by an earlier cppDE are refused with a message; generate and
  compile them again.

* `solveODE()`, `solveODEBatch()` and `prepareBatch()` take `maxOrder`, an
  upper bound on the order of `"bdf"` and `"adams"`, and `stald`, stability
  limit detection for `"bdf"` (off by default), on both backends. Both are
  set per solve; the compiled model does not change.

# cppDE 0.13.1

* The demos use ASCII only, and `scales`, which `demo("fitCoPt")` calls, is
  in `Suggests`.

# cppDE 0.13.0

* Breaking: `pi` is the number in every expression, as in R; a state,
  parameter, forcing or variable named `pi` is an error. `E` is a name;
  Euler's number is `exp(1)`.

* Breaking: expressions read as R reads them. `&` and `|` are `&&` and `||`,
  `round()` rounds half to even, and `ifelse()` and `trunc()` are new. A
  chained comparison such as `t1 < time < t2`, `%%`, `%/%`, `%`, `//`, `~`,
  `@`, `<-` without a space and an unknown function are errors that name the
  construct.

* A state switch whose branches both push the solution onto the switching
  surface slides along it (Filippov) instead of stopping with an error, and several
  switches slide at once. Values and first- and second-order sensitivities on
  the native backend, first-order sensitivities on `cvode()`, reverse
  gradients on both, and the forward-reverse curvature.

* Switches that cross at the same time take one saltation each. A root
  crossed at a small rate takes the saltation: the crossing counts as
  tangential only where the rate is negligible against the terms it is summed
  from, not below a fixed 1e-15.

* A switch whose root is zero at `t0`, at a fixed event or after a crossing of
  that root takes the side the solution leaves to, on both backends, and the
  sensitivities take the jump once.

* `abs()`, `min()` and `max()` of a state or of the time are switches: the
  solver locates the kink as a root, as for a comparison. A comparison of
  switches follows their modes. The sensitivities of every order take the kink
  or jump. `floor()`, `ceiling()`,
  `round()` and `trunc()` of a state or of the time in an ODE right-hand side
  are an error.

* At a kink or step the derivative is that of the branch the value takes,
  forward and reverse alike, in `cppFUN()` and the symbolic Jacobians.

* `solveODE()`, `solveODEBatch()` and `solveBatch()` reject a parameter
  named twice in `parms`.

* The demos `basics`, `saltation` and `fitCoPt` replace the scripts in
  `inst/examples/`. `saltation` checks the sensitivities across events,
  switches and sliding against symbolic solutions.

# cppDE 0.12.1

* A switch on the time, such as the pulse
  `piecewise(1, time > ts && time <= t2, 0)`, is located as a root like a
  switch on a state, on both backends. `cvode()` resolves such a pulse also
  where the right-hand side is zero before it, and the sensitivities of the
  native backend take the jump where a parameter moves the switching time.
* `!` after an operator, as in `k*!(x > 1)`, parses. It negates as in R, up to
  the next `&&`, `||`, comma or closing bracket.
* `?expressions` explains arithmetic, comparisons and logic, with an example
  that solves the same equations with `cppODE()` and `cvode()`.

# cppDE 0.12.0

* New `adjointControl()`, passed as `adjoint` to `solveODE()`,
  `solveODEBatch()`, `prepareBatch()` and `solveBatch()`, controls what a
  reverse solve does beyond the gradient. With `refine = TRUE` the sweep is
  checked step by step against the local error test of the CVODES backward
  problem: the cotangent of the states to `reltol` and `abstol`, and with
  `gradtol` each step's share of the gradient to `reltol` and `gradtol`. A
  tolerance thereby means the same as on a `cvode()` reverse model. A step
  that fails is swept again in substeps. `bdf` and `adams` sweep the flow
  between their grid points by extrapolated linearly implicit Euler instead
  of the multistep scheme, so their gradient is that of the flow along the
  computed trajectory. With `trace = TRUE` the solve returns `$adjoint`, the
  grid and adjoint of the sweep and, under `refine`, the substeps each step
  took.
* Breaking: `solveODE()`, `solveODEBatch()` and `prepareBatch()` lose
  `adjointGrid` and `errWeights`, and `solveBatch()` loses `errWeights`.
  `adjoint = adjointControl(trace = TRUE)` returns the grid as `$adjoint` in
  place of `$adjointGrid`; the step-size weights have no replacement.
* Breaking: the value pass of a reverse solve takes the grid of a value run,
  as the forward run of CVODES' adjoint analysis does. `sensErrCon` integrates
  no control tangent there and has no effect on a reverse model.
  Where that grid misses the gradient, as on a trajectory at rest or under a
  forcing the state does not see, use `adjointControl(refine = TRUE)`.
* Breaking: the limit of consecutive rejected steps is the argument
  `maxattempts` of `solveODE()`, `solveODEBatch()`, `prepareBatch()` and of the
  per-condition lists.
* `install_libs()` is renamed `installLibs()`; the old name remains as a
  deprecated alias.
* `ptc()` takes the tolerances in `controls` as `reltol` and `abstol`, as the
  solvers do. The old names `rtol` and `atol` still work and warn.
* The help pages state every default, list the fields of `$diagnostics` and
  the columns of `$trace` in `?solveODE`, and give the signatures and shapes
  of the functions `cppFUN()` returns. `?solveODE` marks `cotangent` as
  optional where `keepStore = TRUE` lets a reverse solve run without one.
* `cvode()` reverse models take events and `rootfunc`: the forward run is
  split at the jumps, one CVODES memory per stretch, and the backward sweep
  takes the adjoint and the gradient through each jump, including event times
  that move with the parameters or the state. They also take
  `adjointControl(gradtol =)`, which puts the gradient quadrature into the
  backward problem's error test.
* In `cvode()` reverse models with a sparse Jacobian every KLU setup
  factorizes afresh, so the replay from a checkpoint in the backward pass
  takes the steps the run stored.
* `adams` solves the corrector by Newton's method, as CVODES does with a linear
  solver attached, in place of the PECE fixed point; the step trace loses
  `pece_iters` and `pece_diverged`.
* The Newton step on an iteration matrix factorized at an earlier gamma takes
  the factor gamma / gammap, and the BDF scaling 2 / (1 + gamrat) only under
  BDF.
* The reverse sweep of `bdf` and `adams` reuses the operators of a step only
  for a later step of the same size whose history rescales by the same
  factor.
* `rb4` has fourth order also on a right-hand side that reads the time or a
  forcing: the coefficient of df/dt in the fourth stage has the sign of the
  method.
* Reverse sweeps: the transposed solve refines against a held
  factorization where that is cheaper than factorizing each step, generated
  kernels keep their temporaries on the stack, the pullbacks of the
  right-hand side and of the Jacobian take one pass each, and KLU skips the
  check of the matrices the generated code assembles.
* `solveBatch()` gains `store`, the checkpoint stores of a prepared reverse
  batch.
* `compile(output =)` builds the objects of the combined library `cores` at a
  time.
* `installLibs()` and `./configure` build SUNDIALS 7.9.0 and SuiteSparse
  7.14.1 by default, with the CMake option names of SUNDIALS 7.7 and later;
  older pins keep the old names. `./configure` prefers the cache prefix of the
  pinned versions over the lexically last one.
* `adjointControl(trace = TRUE)` returns `failures` in `$adjoint`, the checked
  steps that met the test at no number of substeps; `eta` is `NA` where `bdf`
  and `adams` sweep the flow.
* The colored Jacobian, which `CPPDE_JAC = "colour"` forces and large sparse
  models choose, compiles: the reverse Jacobian-vector product seeds no
  parameters.
* The checked sweep keeps the seeds before a step's start, and the
  extrapolated flow takes intervals in either direction of time.
* `cvode()` reverse models on KLU factorize afresh through
  `SUNLinSol_KLUReInit()` rather than a SUNDIALS internal.
* A comparison in the right-hand side whose sides differ by an amount that
  reads a state, in `piecewise()` or a logical operator, and `Heaviside()` or
  `sign()` of such an amount switch where that amount crosses zero, on both
  backends. The solve holds each switch as a mode, locates the root and sets
  the mode there, so forward tangents of first and second order, the reverse
  gradient and the curvature, and `cvode()`'s forward and adjoint
  sensitivities take the jump of f. A switch adds no output row and is not
  limited by `maxroot`. A solution that slides along a switching surface stops
  with an error.
* `cvode()` reverse models with root events interpolate the forward state by
  Hermite, also at the start of every stretch after a root.
* Root events that cross zero within one step fire one by one at their
  crossing times, the earliest first.
* A terminal `rootfunc` ends the output with one row at the root.
* `./configure` finds a SUNDIALS built with MPI, as Debian ships it: it adds
  the MPI link flags, so the probe links and `cvode()` is enabled.
* `erf()`, `erfc()` and `atan2()` take forward derivatives of first and second
  order, so `cppODE()` with `derivMode = "forward"` or `"forward-forward"` and
  `cppFUN()` accept them as `cvode()` and the reverse modes do.

# cppDE 0.11.4

* `forcingValues()` builds with clang and libc++ (macOS): `R_NO_REMAP` before
  the R headers, as in `ptc()`.

# cppDE 0.11.3

* A forcing holds the value of the first and last point outside them, and a
  single point gives a constant. Both apply to
  `cppODE()` and `cvode()` models.
* New `forcingValues()` evaluates forcings as the models see them.

# cppDE 0.11.2

* `ptc()` builds with clang and libc++ (macOS): `R_NO_REMAP` before the R
  headers.

# cppDE 0.11.1

* `ptc()`: the absolute floor of the residual scale of a row is `atol` times the largest
  partial derivative in that row, not the diagonal entry.

# cppDE 0.11.0

* New `ptc()`: steady state of a `cppFUN()` model by pseudo-transient
  continuation, with linear constraints `C x = total` and positive variables.
  Implicit Euler steps under local error control follow the flow to a stable
  steady state and become Newton steps near it; `flow = FALSE` solves plain
  equations.
* `rootfunc = "equilibrate"` stops once every `|dx/dt| <= roottol |x| + abstol`,
  in both backends, in place of the absolute `|dx/dt| < roottol`.

# cppDE 0.10.5

* `derivMode = "reverse"` and `"forward-reverse"` return the gradient the
  forward mode returns also where the trajectory is at rest or small against
  `abstol / reltol`. The forward run integrates a control tangent,
  `z' = J(x) z` along a fixed direction of the initial state, on the same steps
  as the state, and the error of that tangent joins the error of the state in
  the error test, as the error of a sensitivity does under `sensErrCon`.
  `sensErrCon = FALSE`, accepted by both reverse modes, keeps the grid of a
  value-only run. Reverse models with `method = "tsit5"` generate their
  Jacobian for it.
* With `sensErrCon = FALSE` the order-decrease candidate of `"bdf"` and
  `"adams"` reads the state alone, as the error test does.

# cppDE 0.10.4

* A pulse such as `piecewise(1, time > ts && time <= t2, 0)` on a state at
  rest is not stepped over, under every method, with and without
  sensitivities. `cppODE()` solves every comparison, `Heaviside()` and `sign()`
  of the right-hand side on time and parameters alone that is affine in time
  for the switching time. The solve stops a few doubles in front of each
  switching time inside the grid, follows `f` over them to the jump and
  restarts past it, as at a jump in time. The switching times add no output
  rows, and the reverse mode transposes the crossing.
* After an event, `"rb4"` and `"tsit5"` with dense output also restart their
  step-size controller, and `"tsit5"` evaluates the first stage anew.

# cppDE 0.10.3

* A `"bdf"` or `"adams"` solve crosses a jump of the right-hand side in time,
  such as the step input `piecewise(0, time - ts < 0, 1)`, also from rest under
  a tight `abstol`. When the step size fails at the resolution of `t`, the
  solve locates the jump of `f(t, x)` in `t` within the failed
  step, follows `f` over the few doubles up to it and restarts on the far side
  as after an event. With no jump ahead, it restarts where it stands, once per
  point of time. The reverse mode transposes both. A stall that persists stops
  at once with "Step size fell below the resolution of t".

# cppDE 0.10.2

* A `"bdf"` or `"adams"` solve with sensitivities that starts late under a
  tight tolerance takes a first step. The first-step estimate takes the
  geometric mean of the bounds where a sensitivity crosses them, as CVODES'
  `cvHin` does, and the first step of every method is at least four ticks of
  `max(|t0|, |t_final|, 1)`.
* A fixed-time event fires only from the first time of the integration grid
  up to, but not including, the last one, and `includeTimeZero` adds 0 to that
  grid. An event before the grid or at the last time is not applied, and none
  adds a row beyond the last requested time. Two consecutive solves over
  `[t0, t1]` and `[t1, t2]` reproduce one solve over `[t0, t2]`, which
  multiple shooting relies on. A grid of a single time applies the events at
  that time, as the start of a longer grid does. The rule holds for event
  times that depend on parameters, in the batch, in the reverse mode and in
  `cvode()`.
* The reverse mode counts the cotangent of the first output row once when a
  fixed event sits at the first time.
* A reverse solve over a single time integrates nothing: the cotangent of the
  initial state is the seed and the parameters receive only what an event at
  that time contributes.
* `cvode()` reports the state after an event at the first time in the first
  row, as the native backend does, and applies no event that lies before the
  first time.

# cppDE 0.10.1

* `cppFUN()` builds `"forward"` only by default. `"forward-reverse"` is a
  separate mode; `"reverse"` does not build it.
* The forward-reverse `vjp` is generated as derivative code instead of running
  `vjp` over a dual.
* Products of four or more factors are differentiated through shared prefix and
  suffix products, and chained adjoints and tangents stay one factor. Reverse
  and forward-reverse code grows linearly, not quadratically, with product length.

# cppDE 0.10.0

* **Reverse mode.** `cppODE(..., derivMode = "reverse")` computes the gradient
  of a seeded functional by a discrete adjoint, at a cost that does not grow
  with the number of parameters. `"forward-reverse"` adds the curvature, the
  Hessian applied to the tangent, and `"forward-forward"` the full Hessian. All
  four methods, events, roots, forcings, the sparse solver and the batch entries
  are covered. `cvode(..., derivMode = "reverse")` wraps the CVODES adjoint
  behind the same interface. `keepStore`/`store` reuse a forward pass,
  `errWeights` weights the step size by the adjoint and `adjointGrid` returns
  the backward grid.
* **Derivative names.** Arguments and results are called `tangent`, `hessian`,
  `cotangent` and `curvature`, replacing `sens1ini`/`sens1`, `sens2ini`/`sens2`,
  `seed`/`adjoint` and `seedTangent`/`adjoint2`. `cppFUN()` takes `tangentX`,
  `tangentP`, `hessianX` and `hessianP`, `evaluate()` returns `y`, `tangent` and
  `hessian`, and `vjp()` covers both orders. `funCpp()` is renamed `cppFUN()`; no
  old name is kept as an alias.
* **Code generation on an expression graph.** The code generator builds a
  hash-consed expression graph of the model and differentiates it by automatic
  differentiation; SymPy is left for unusual syntax. Generation grows linearly
  with the model, long linear sums become tables and regular structure becomes
  loops. `derivMode = "symbolic"` and `derivSymb()` are gone, and a `cppFUN()`
  object runs compiled code only.
* `solveODE(..., sensErrCon = FALSE)` lets forward sensitivities ride the grid
  of a value run. Tangent storage is allocated on the heap, without a
  compile-time limit (`nStack` is gone).
* **Bug fixes.** Second derivatives through fixed and root events, including
  resets and roots that read the clock, agree with forward over forward to
  1e-10. The upper half of `$hessian` does not drift on long stiff runs.
  Unnamed models get distinct names also under `set.seed()`. The CVODES
  adjoint converges on models with many parameters and returns the result
  through the batch. A forcing that multiplies a state compiles, `cvode()`
  reads a `terminal` column given as text, and the `rb4` initial step and
  time-only sparse Jacobian entries are right.
* `codegenAvailable()` reports whether a model can be generated and compiled
  without installing anything; the examples run only then.
* Python 3.9 or newer and R 4.3 or newer are required. The methods vignette
  covers the reverse mode and benchmarks against CVODES, and every export has
  an example.

# cppDE 0.9.5

* OpenMP is detected on Windows: `configure.win` reads
  `SHLIB_OPENMP_CXXFLAGS` from the `Makeconf` under the architecture
  subdirectory, and models build with `-fopenmp`.

# cppDE 0.9.4

* A threaded BLAS does not deadlock a forked worker. BLAS is pinned to one
  thread for the width of every `fork()` and restored in the parent.
* `forkGuard()` reports which BLAS cppDE steers, the BLAS thread count and whether
  the handler is installed. The guard is installed when the DLL loads.
* The BLAS thread-count probe covers FlexiBLAS and BLIS alongside MKL and
  OpenBLAS.
* The Windows branch of that probe finds the BLAS, so the two-runtime guard
  applies there too.
* The package has a `src/`, holding the fork guard's entry point.

# cppDE 0.9.3

* A `piecewise` translates, as `cppde::select(cond, a, b)`. Both branches are
  evaluated, so each has to be safe to evaluate.
* `&&`, `||` and `!` are accepted in equations, grouped by Python's parser.
* An expression that does not parse names itself and gives a reason.
* A model symbol cannot collide with an identifier the generator emits, such
  as `std`.
* A symbol named after a C++ keyword compiles, `default` and `int` included.
* A symbol named after a Python keyword is rejected and named. SymPy parses
  through Python's parser.
* The `double` locals of a root event's `G_tt` lambda are named by position.
* `cppde::value_of(x)` is the value accessor across arithmetic types, both dual
  orders and their expression templates.

# cppDE 0.9.2

* A root event whose crossing falls exactly on an evaluated time fires.
* A root event fires once per crossing.
* A fixed event that makes a root condition true fires it, with the resets
  riding on the jump's surface. Terminal conditions are excluded.
* The step size is re-estimated after every event, not only for the multistep
  methods.
* `funCpp()` substitutes all symbols in one pass.
* `compile()` passes the OpenMP and KLU flags that the constructors recorded
  on the model once.
* `inst/examples/example_saltation.R` checks the sensitivity transport against
  a SymPy solution to first and second order, over five models.
* Generated entry points are dispatched by name and shared object rather than
  by a cached address, so loading and unloading are without consequence. A
  model whose library is gone names what it is missing.
* Scoping the lookup to one shared object also stops two models sharing an
  entry point name from reaching into each other.
* `clearNativeSymbols()` drops the remembered name pairings. It is not needed
  after loading or unloading.

# cppDE 0.9.1

* Test suite over solvers, sensitivities, events, reparameterization and the
  batch path.

# cppDE 0.9.0

* `solveODEBatch()` solves many conditions in one `.Call`, over OpenMP
  threads, with results identical to the serial path.
* `prepareBatch()` and `solveBatch()` reuse a validated handle across solves.
* `batchAvailable()` reports why a batch would run serially.

# cppDE 0.8.4

* `./configure` detects SUNDIALS, SuiteSparse/KLU and OpenMP by linking a
  test binary; a missing library disables only the feature that needs it.
* `install_libs()` builds SUNDIALS and SuiteSparse into a per-user cache.

# cppDE 0.8.3

* `funCpp()` compiles algebraic functions with optional derivatives, through
  forward-mode AD or analytic SymPy derivatives.

# cppDE 0.8.2

* `solveODE()` integrates a model and returns states with first and second
  order parameter sensitivities.
* `diagnostics()` prints the solver statistics.

# cppDE 0.8.1

* `cppODE()` and `cvode()` generate, compile and load a solver for an ODE
  system and return a model handle.

# cppDE 0.8.0

* `compile()` builds generated sources through `R CMD SHLIB`.
* Native symbol lookups are cached; `clearNativeSymbols()` drops the cache.

# cppDE 0.7.10

* C++ generation for ODE, CVODE and funCpp models, with common
  subexpression elimination and Jacobian sparsity detection.

# cppDE 0.7.9

* `derivSymb()` exposes symbolic first and second derivatives through SymPy.

# cppDE 0.7.8

* Event engine with saltation corrections, root finding and PCHIP forcings.

# cppDE 0.7.7

* Rosenbrock4 and Tsit5 single-step steppers with embedded error estimate
  and dense output.

# cppDE 0.7.6

* Variable-order BDF/NDF and Adams multistep steppers in Nordsieck form,
  with a Newton corrector.

# cppDE 0.7.5

* AD aware dense LU through LAPACK and sparse LU through KLU.

# cppDE 0.7.4

* Thread-local bump arena and contiguous tangent slab as tangent storage.

# cppDE 0.7.3

* `dual2nd<T, N>` for second order sensitivities.

# cppDE 0.7.2

* Expression templates collapse right-hand side temporaries into one fused
  chain-rule loop.

# cppDE 0.7.1

* Forward-mode dual numbers `dual<T, N>` for first order sensitivities,
  with a static and a heap allocated specialization.

# cppDE 0.7.0

* Package skeleton and license.
