# cppDE (development version)

* A reverse seed can carry its own tangents. `solveODE(..., seedTangent = )`
  takes `[n_out, n_states, n_seed, n_sens]` beside the seed, and forward over
  reverse fills the sweep's dual seed from it. A cotangent handed down by
  something sitting above the ODE moves with the parameters, and without this
  the second order silently dropped that motion: the gradient was right and the
  Hessian was wrong in the columns of whatever the node above depended on. Rides
  as an attribute of the seed, for the same reason `errWeights` does.
* `derivMode = "reverse"` on `cppFUN()` builds `vjp2` beside `vjp`. The written
  vector-Jacobian product no longer names a scalar type, and instantiated over a
  dual it returns `J' dW` and `W' H V` in one pass, so an observation function
  contributes its curvature to a backward chain without a symbolic Hessian ever
  being formed. `vjp2(vars, params, w, vx, vp, dw)` takes the tangents the
  inputs carry and those of the cotangent, and answers with `dwx` and `dwp`
  beside what `vjp` returns.
* `derivMode = "forward-forward"` compiles with a sparse Jacobian. The sparse
  AD solver keeps its iteration matrix and asks how many derivative directions
  it carries, and `max_deriv_size` had overloads for a vector and for a dense
  matrix but not for a compressed-column one. Only the nested-dual branch needs
  it, and that branch had never been compiled.
* `store` and `keepStore` are refused under `derivMode = "forward-reverse"`. A
  checkpoint keeps its value inline and its tangents behind a pointer into the
  arena, and the arena is reset when the solve that filled it returns. Handing
  such a store to a later solve gave exact values, an exact gradient and a wrong
  Hessian, with nothing to show for it. Until the store owns its tangents, the
  second solve integrates.
* `solveODE(..., sensErrCon = FALSE)` is a cheaper and coarser forward mode.
  The step size, the order and the corrector's convergence test then read value
  arithmetic only, so the sensitivities ride the grid a value run takes instead
  of the finer one their own error would ask for. Worth it where the step count
  matters more than the last digits of a tangent: the maximum is taken over the
  state and every direction, so on a wide sensitivity set the worst-resolved
  direction otherwise sets the step for all of them, and the step count grows
  by about a factor of three before it saturates. This is the convention
  SUNDIALS ships, where the same switch is off by default. It carries through
  `solveODEBatch()` and `prepareBatch()`, batch-wide or per condition, and
  needs a model that has sensitivities.
* A model built with `derivMode = "forward-reverse"` differentiates the grid a
  value run takes. The error test, the order choice, the corrector's
  convergence test and the first step now read value arithmetic only, so the
  step sequence no longer depends on how many tangent directions ride along or
  on what they contain. Three things follow. A Hessian assembled from several
  blocks of directions is one matrix rather than several: `sens1ini` may be
  split into chunks of any width and the result is the same to the last bit.
  Such a model now takes the same grid as `derivMode = "reverse"`, so its
  gradient belongs to the trajectory a value run produces rather than to a
  finer one of its own; the two agree to rounding rather than exactly, because
  a corrector sums in a different order over the AD type than over `double`.
  And the step count stops growing with the direction count. `derivMode = "forward"`
  and `"forward-forward"` are untouched and keep the sensitivity error control.
* The methods vignette covers the reverse mode as it now stands. It derives
  the adjoint equation and its quadrature, separates the discrete adjoint from
  the continuous one, gives the transposed saltation relation and the restart
  collapse for events, and states the second-order adjoint system that forward
  over reverse discretises. The checkpoint store and the lambda-weighted
  controller are documented, and the CVODE chapter gains a section on CVODES
  adjoint sensitivity analysis: the backward problem, the checkpoint
  interpolation, the seeded jumps and the restrictions.
* **Bug fix.** A forcing that multiplies a state produced a source that did not
  compile. The Jacobian entries were printed without the forcing list, so a
  forcing surviving differentiation came out as a bare identifier that nothing
  declares. Only an additive forcing vanishes from `df/dx`, which is why every
  example and every test carried one.
* `Heaviside` can be differentiated. Its derivative is `DiracDelta`, which no
  printer knows, so the ODE generator gave up on any model that used it. A
  discrete model cannot mean an impulse of infinite height, so it is emitted as
  one at the switching point and zero either side, the way `cppFUN()` has always
  taken it.
* A model built with `derivMode = "reverse"` also emits the two contractions a
  written step adjoint asks for, `J' lambda` and `(df/dp)' lambda`, in plain
  `double`. They replace what the tape derives per step. Internal for now: the
  step adjoints that call them are still being written.
* `cppODE(method = "rb4", sparse = TRUE, derivMode = "reverse")` says no rather
  than failing in the compiler. The Rosenbrock replay builds its own dense
  Jacobian and forms the residual over every entry, so a sparse model has no
  backward path there. Every other combination of method, linear algebra and
  direction stands.
* **Bug fix.** A CVODES adjoint solved through `solveODEBatch()` returned no
  `$adjoint` and reported success. The batch sizes its results before the solve
  where the output grid is fixed by `times`, and that skeleton has no slot for
  an adjoint whose width is only known afterwards; the native backend declines
  the same shortcut under `derivMode = "reverse"` and the CVODE emitter did
  not. A caller reading `$adjoint` saw `NULL`, which reads as a zero gradient
  rather than as a failure.
* The batch entry points carry the whole reverse mode. `solveODEBatch()` and
  `prepareBatch()` accept `adjointGrid`, `errWeights`, `keepStore` and `store`
  beside `seed`, batch-wide or per condition, so the checkpoint reuse and the
  weighted controller reach a batch of conditions and not only `solveODE()`.
  `solveBatch()` takes a new `seed` and new `errWeights` on a prepared handle:
  both change with every objective evaluation while their shapes do not, which
  is what makes a prepared batch usable in reverse mode at all.
* A reverse gradient can integrate the states once instead of twice. It needs
  two solves at the same parameter, one for the values the seed is built from
  and one for the sweep. `solveODE(..., keepStore = TRUE)` now returns the
  checkpoints of the first as `$store`, which the second takes through
  `store = `. A store handed back at a different point is an error rather than a
  silent reuse. Worth about a tenth of a reverse gradient on a stiff model; the
  checkpoint capture the first solve now pays takes back part of it.
* The reverse mode's tape records about three times faster. Liveness was tested
  four times per binary operation: twice to decide whether to record, twice
  again inside the recording. Each test was three comparisons. It is now
  one comparison, evaluated once. Recording fell from 9.3 to 3.1 nanoseconds per
  tape node, and a reverse step from 167 to 79 times a plain evaluation of the
  same right-hand side.
* The CVODE backend takes derivatives backwards too. `cvode(..., derivMode =
  "reverse")` compiles CVODES adjoint sensitivity analysis: the forward pass
  stores checkpoints and one backward solve per seed column integrates the
  adjoint equation with the parameter quadrature riding along. It answers the
  same `solveODE(..., seed = W)` interface and returns the same `$adjoint`. It
  refuses events and `rootfunc`, which the native backend carries, because
  CVODES integrates the adjoint over checkpointed states and has no way to be
  told about a jump. Two independent implementations of the same mathematics
  are the point: neither is an oracle for the other, but a systematic error in
  one would show.
* The step-size controller can be weighted by the adjoint. `solveODE(..., seed
  = W, errWeights = list(time = , lambda = ))` adds a term to the error norm
  under the same maximum the sensitivity columns already use, so the grid stays
  at least as fine as `abstol` and `reltol` ask and is finer only where the
  adjoint says a step carries objective error. A weight that is wrong therefore
  costs time and never accuracy, which is what makes an estimate from an
  earlier run at a nearby parameter usable. Both controllers carry it.
* A reverse solve can report the grid it swept. `solveODE(..., seed = W,
  adjointGrid = TRUE)` returns `$adjointGrid` with the step times and sizes, the
  cotangents of both, the adjoint state per step, and `eta`, the dual-weighted
  residual that says how much of the objective's error each step carries. The
  flag rides on the seed rather than on a new argument, so no compiled model
  needs rebuilding.
* Derivatives can be taken backwards. `cppODE(..., derivMode = "reverse")` compiles a
  fourth object beside the value, first- and second-order ones: it integrates the
  states in plain `double`, keeps a checkpoint per accepted step, and replays each
  step under a new reverse-mode scalar to sweep one tape backwards.
  `solveODE(..., seed = W)` hands it a cotangent of the trajectory and gets back
  `$adjoint`, one row per state and parameter, at a cost that does not grow with
  the number of parameters. All four methods carry it, with events, root events,
  forcings, the batch entry and the sparse solver.
* The reverse mode's answer belongs to the trajectory a value-only solve produces,
  so a value and its gradient are consistent with each other. Forward
  sensitivities adapt on a finer step sequence, because the error norm takes the
  maximum over every tangent column, so their gradient belongs to a different
  discretisation than the value they are reported with. The two agree to O(tol).
* `cppFUN(derivMode = )` names derivative *directions* rather than a backend, and
  more than one may be asked for. `"dual"` is gone; it becomes `"forward"` for
  the Jacobian, `"reverse"` for the vector-Jacobian product, or both. The two are
  now built independently: the reverse entry instantiates the expression body a
  second time over `cppde::codual`, so a caller that only multiplies by the
  Jacobian no longer compiles it. `"symbolic"` is unchanged and, being a backend
  for the forward Jacobian rather than a direction, cannot be combined with
  either.
* `funCpp()` is now `cppFUN()`, for symmetry with `cppODE()`. The old name is
  gone rather than deprecated, the package never having been released under it.
* A sparse Jacobian's reused pivot order is checked. `klu_refactor` keeps the
  ordering the first factorisation chose and reports success even where it has
  become numerically hopeless; a Newton corrector converges anyway, so nothing
  forward ever noticed, while a reverse solve uses the result once and directly
  and was four orders of magnitude out on a stiff model. The reciprocal pivot
  growth now decides whether to refactor fully.
* The code generator no longer recognises scalar types by name. What it needs to
  know about the type, how many derivative layers and whether they live in the
  arena, is stated by the caller, and the generated code spells
  `cppde::dual<double, N>` and `cppde::codual<double>` out. A name it did not
  recognise used to pass silently through every AD branch.

# cppDE 0.9.5

* **Bug fix.** OpenMP is detected on Windows. `configure.win` read
  `SHLIB_OPENMP_CXXFLAGS` from `R_HOME/etc/Makeconf`, which does not exist
  there -- R keeps that file under the architecture subdirectory. Detection
  therefore reported "no SHLIB_OPENMP_CXXFLAGS in Makeconf" on every Windows
  install, `solveODEBatch()` ran serially and generated models were built
  without `-fopenmp`. Both configure scripts now look under `etc/$R_ARCH`
  first, and the Windows summary line reports OpenMP alongside CVODE and KLU.

# cppDE 0.9.4

* A threaded BLAS no longer deadlocks a forked worker. Its worker threads do not
  survive `fork()`, and the first call in the child large enough to thread hangs
  on a lock they held. BLAS is now pinned to one thread for the width of every
  `fork()` and the previous count is restored in the parent.
* `forkGuard()` reports which BLAS cppDE steers, its thread count and whether the
  handler is installed. Attaching the package says the same in one line, which
  `options(cppDE.quiet = TRUE)` suppresses. The guard is installed when the DLL
  loads, so it is in place whether or not the package is attached.
* The BLAS thread-count probe covers FlexiBLAS and BLIS alongside MKL and
  OpenBLAS.
* The Windows branch of that probe searched the process image, which never
  exports these entry points, so the guard against two live OpenMP runtimes did
  nothing there.
* The package has a `src/`, holding the fork guard's entry point and nothing
  else. `./configure` writes `src/Makevars` after probing whether `dlsym()`
  needs `-ldl`.

# cppDE 0.9.3

* A `piecewise` translates. Comparisons are defined on the AD nodes of both
  orders and a piecewise is emitted as `cppde::select(cond, a, b)`. Both
  branches are evaluated, so each has to be safe to evaluate.
* `&&`, `||` and `!` are accepted in equations. Python's own parser does the
  grouping, so `a > b && c > d` keeps its meaning without parentheses.
* An expression that does not parse names itself and gives a reason on one
  line. The old message was truncated by reticulate and then indexed out of
  range.
* A model symbol can no longer collide with an identifier the generator emits.
  Symbols are substituted for their slot while the expression is printed, not
  in the finished source, where a parameter named `std` rewrote `std::pow`
  into `p[18]::pow`.
* A symbol named after a C++ keyword compiles, `default` and `int` included.
* A symbol named after a Python keyword is rejected and named in the message,
  in `funCpp()` as well as in `cppODE()` and `cvode()`. It used to be renamed,
  which left the caller holding the old name. SymPy parses through Python's
  parser, where such a name is a syntax error.
* The `double` locals of a root event's `G_tt` lambda are named by position,
  not after the model's own symbols.
* `cppde::value_of(x)` is the value accessor across arithmetic types, both dual
  orders and their expression templates.

# cppDE 0.9.2

* A root event whose crossing falls exactly on an evaluated time now fires.
  Detection and bisection both tested the sign product strictly, so an exact
  zero counted as no crossing and the event was lost or localised past the root.
* A root event no longer fires again on the crossing it just handled. The
  restart sits on the event surface, where the round-off residue of the root
  function carried a sign that read as a second crossing.
* A fixed event that makes a root condition true now fires it. The conditions
  are read on both sides of the jump and the resets ride on its surface, so
  they transport the sensitivities like a fixed event at that time. Terminal
  conditions are excluded.
* The step size is re-estimated after every event, not only for the multistep
  methods.
* `funCpp()` substitutes all symbols in one pass. A parameter carrying the name
  of a generated array, `p` for instance, rewrote the slots already emitted for
  the others, so the result depended on the order the parameters were listed in.
* `compile()` no longer repeats the OpenMP and KLU flags that the constructors
  already recorded on the model.
* `inst/examples/example_saltation.R` checks the sensitivity transport against a
  SymPy solution of the same model, to first and second order, over five models
  covering both event kinds, repeated firings, an explicitly time-dependent
  right-hand side and an oscillator between two elastic walls.
* Generated entry points are dispatched by name and shared object instead of by
  a cached address. `dyn.unload()` nulls an address in place and nothing
  resolves it again, so a reload left every caller that had already resolved a
  symbol pointing at nothing. Loading and unloading are now without
  consequence, and a model whose library is gone names the entry point and the
  library it is missing instead of dying on a null address.
* Scoping the lookup to one shared object also stops two models that export the
  same entry point name from reaching into each other.
* `clearNativeSymbols()` drops the remembered name pairings, which only a
  recompile into a differently named shared object can make stale. It is no
  longer needed after loading or unloading.

# cppDE 0.9.1

* Test suite over solvers, sensitivities, events, reparametrisation and the
  batch path.

# cppDE 0.9.0

* `solveODEBatch()` solves many conditions in one `.Call`, over OpenMP
  threads, with results identical to the serial path.
* `prepareBatch()` and `solveBatch()` reuse a validated handle across solves.
* `batchAvailable()` reports why a batch would run serially.

# cppDE 0.8.4

* `./configure` detects SUNDIALS, SuiteSparse/KLU and OpenMP by linking a
  test binary; a missing library disables only its own feature.
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
  with a static and a heap allocated specialisation.

# cppDE 0.7.0

* Package skeleton and licence.
