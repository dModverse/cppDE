# cppDE 0.11.1

* `ptc()`: the absolute floor of a row's residual scale is `atol` times its largest
  partial derivative, not its diagonal. A row whose own rate is tiny but whose
  partner's is not (Smad2 next to pSmad2) was scaled up to 1e12 times too tight,
  the step matrix became numerically singular and `ptc()` only shrank `dt`.

# cppDE 0.11.0

* New `ptc()`: steady state of a `cppFUN()` model by pseudo-transient
  continuation, with linear constraints `C x = total` and positive variables.
  Implicit Euler steps under local error control follow the flow to a stable
  steady state and become Newton steps near it; `flow = FALSE` solves plain
  equations.
* `rootfunc = "equilibrate"` stops once every `|dx/dt| <= roottol |x| + abstol`,
  in both backends. The absolute `|dx/dt| < roottol` stopped small states too
  early and large ones too late.

# cppDE 0.10.5

* **Bug fix.** `derivMode = "reverse"` and `"forward-reverse"` return the
  gradient the forward mode returns where the trajectory is at rest. Their
  forward run took its steps from the state's error alone. On a trajectory
  that stays exactly zero, such as a state that starts at 0 and is not driven
  within the grid, that error is zero, the step grew without bound, and the
  adjoint on that grid was off by more than its size: +0.80 against -0.43 for a
  derivative with respect to a zero initial state on `0:50`. Where a state is
  small against `abstol / reltol` the grid was too coarse for the same reason,
  which left the gradients of the PEtab benchmark models up to 1e-5 off at a
  tolerance of 1e-8. The forward run now integrates a control tangent,
  `z' = J(x) z` along a fixed direction of the initial state, on the same steps
  as the state, and its error joins the state's in the error test as a
  sensitivity's does under `sensErrCon`. On those models the reverse solve
  now takes about as many steps as the forward mode and gets its accuracy, at
  1.4 to 6 times the time it took before. `sensErrCon = FALSE`, now accepted by
  both reverse modes, keeps the grid of a value-only run. Reverse models with
  `method = "tsit5"` generate their Jacobian for it.
* With `sensErrCon = FALSE` the order-decrease candidate of `"bdf"` and
  `"adams"` reads the state alone, as the error test does. It read the
  sensitivities, so the grid could depend on them.

# cppDE 0.10.4

* **Bug fix.** A pulse such as `piecewise(1, time > ts && time <= t2, 0)` on a
  state at rest is no longer stepped over. Before `ts` the right-hand side and
  with it the error estimate are zero, so the step size grew until one step
  went across the whole pulse: without sensitivities under `"bdf"` and
  `"adams"` for the pulse from 60 to 90 on `0:180`, and for a short pulse late
  in the grid under every method, with sensitivities too. `"rb4"` with
  `useDenseOutput = FALSE` stopped at such a jump with "Maximum number of steps
  exceeded". `cppODE()` now solves every comparison, `Heaviside()` and `sign()`
  of the right-hand side on time and parameters alone that is affine in time
  for its switching time. The solve stops a few doubles in front of each
  switching time inside the grid, follows `f` over them to the jump and
  restarts past it, as at a jump in time since 0.10.3. The switching times add
  no output rows, and the reverse mode transposes the crossing.
* After an event, `"rb4"` and `"tsit5"` with dense output also restart their
  step-size controller, and `"tsit5"` evaluates its first stage anew rather
  than reusing the last stage from before the event.

# cppDE 0.10.3

* **Bug fix.** A `"bdf"` or `"adams"` solve crosses a jump of the right-hand
  side in time, such as the step input `piecewise(0, time - ts < 0, 1)`. From
  rest under a tight `abstol`, no step containing the jump passed the error
  test above the resolution of `t`, and a solve without sensitivities stopped
  in front of it with "Too many failed steps in dense output stepper". A step
  that did cross it left the history with the right-hand side from before the
  jump, and the solve stopped just behind it. When the step size fails at that
  floor, the solve now locates the jump of `f(t, x)` in `t` within the failed
  step, follows `f` over the few doubles up to it and restarts on the far side
  as after an event. With no jump ahead, it restarts where it stands, once per
  point of time. The reverse mode transposes both. A stall that persists stops
  at once with "Step size fell below the resolution of t".

# cppDE 0.10.2

* **Bug fix.** A `"bdf"` or `"adams"` solve with sensitivities that starts late
  under a tight tolerance takes its first step again. A sensitivity starting
  at zero with a rate far above `abstol` crosses the bounds of the first-step
  estimate; the port of CVODES' `cvHin` then took the upper bound instead of
  their geometric mean, and at `t0 = 280` with `abstol = 1e-12` that step lay
  below a tick of `t`. `t + h` rounded back to `t`, and the solve stopped with
  "Too many failed steps in dense output stepper". The estimate now takes the
  geometric mean, as `cvHin` does, and the first step of every method is at
  least four ticks of `max(|t0|, |t_final|, 1)`.
* **Bug fix.** A fixed-time event fires only from the first time of the
  integration grid up to, but not including, the last one, and
  `includeTimeZero` adds 0 to that grid. An event before the grid used to move
  the start of the integration back to its own time, and one after it added a
  row beyond the last requested time. An event at the last time is no longer
  applied, which also ends the failure it caused under `"bdf"` and `"adams"`.
  Two consecutive solves over `[t0, t1]` and `[t1, t2]` now reproduce one solve
  over `[t0, t2]`, which multiple shooting relies on. A grid of a single time
  applies the events at that time, as the start of a longer grid does. The
  rule holds for event times that depend on parameters, in the batch, in the
  reverse mode and in `cvode()`.
* **Bug fix.** The reverse mode counted the cotangent of the first output row
  twice when a fixed event sat at the first time, so the derivatives with
  respect to the initial state were wrong.
* **Bug fix.** A reverse solve over a single time no longer crashes R. It
  integrates nothing, so the cotangent of the initial state is the seed and the
  parameters receive only what an event at that time contributes. The sweep
  read the state count off the first checkpoint, which such a run never
  writes.
* `cvode()` reports the state after an event at the first time in the first
  row, as the native backend does, and no longer applies an event that lies
  before the first time.

# cppDE 0.10.1

* `cppFUN()` builds `"forward"` only by default. `"forward-reverse"` is a mode of
  its own; `"reverse"` no longer builds it.
* The forward-reverse `vjp` is generated as derivative code instead of running
  `vjp` over a dual. It compiles in a fraction of the time and memory and runs
  about twice as fast.
* Products of four or more factors are differentiated through shared prefix and
  suffix products, and chained adjoints and tangents stay one factor. Reverse
  and forward-reverse code grows linearly, not quadratically, with product length.

# cppDE 0.10.0

* **Reverse mode.** `cppODE(..., derivMode = "reverse")` computes the gradient
  of a seeded functional by a discrete adjoint, at a cost that does not grow
  with the number of parameters. `"forward-reverse"` adds its curvature, the
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
  `hessian`, and `vjp()` covers both orders. `funCpp()` is now `cppFUN()`; no
  old name is kept as an alias.
* **Code generation on an expression graph.** The code generator builds a
  hash-consed expression graph of the model and differentiates it by automatic
  differentiation; SymPy is left for unusual syntax. Generation grows linearly
  with the model, long linear sums become tables and regular structure becomes
  loops, so models with thousands of states generate and compile in about a
  minute. `derivMode = "symbolic"` and `derivSymb()` are gone, and a `cppFUN()`
  object runs compiled code only.
* `solveODE(..., sensErrCon = FALSE)` lets forward sensitivities ride the grid
  of a value run. Tangent storage is allocated on the heap, without a
  compile-time limit (`nStack` is gone).
* **Bug fixes.** Second derivatives through fixed and root events, including
  resets and roots that read the clock, agree with forward over forward to
  1e-10. The upper half of `$hessian` no longer drifts on long stiff runs.
  Unnamed models no longer repeat their names under `set.seed()`. The CVODES
  adjoint converges on models with many parameters and returns its result
  through the batch. A forcing that multiplies a state compiles, `cvode()`
  reads a `terminal` column given as text, and the `rb4` initial step and
  time-only sparse Jacobian entries are right.
* `codegenAvailable()` reports whether a model can be generated and compiled
  without installing anything; the examples run only then.
* Python 3.9 or newer and R 4.3 or newer are required. The methods vignette
  covers the reverse mode and benchmarks against CVODES, and every export has
  an example.

# cppDE 0.9.5

* **Bug fix.** OpenMP is detected on Windows. `configure.win` read
  `SHLIB_OPENMP_CXXFLAGS` from `R_HOME/etc/Makeconf`, which on Windows lives
  under the architecture subdirectory, so every install ran serially and built
  models without `-fopenmp`.

# cppDE 0.9.4

* A threaded BLAS no longer deadlocks a forked worker. BLAS is pinned to one
  thread for the width of every `fork()` and restored in the parent.
* `forkGuard()` reports which BLAS cppDE steers, its thread count and whether
  the handler is installed. The guard is installed when the DLL loads.
* The BLAS thread-count probe covers FlexiBLAS and BLIS alongside MKL and
  OpenBLAS.
* The Windows branch of that probe searched the process image, which exports
  none of these, so the two-runtime guard did nothing there.
* The package has a `src/`, holding the fork guard's entry point.

# cppDE 0.9.3

* A `piecewise` translates, as `cppde::select(cond, a, b)`. Both branches are
  evaluated, so each has to be safe to evaluate.
* `&&`, `||` and `!` are accepted in equations, grouped by Python's parser.
* An expression that does not parse names itself and gives a reason.
* A model symbol can no longer collide with an identifier the generator emits.
  A parameter named `std` used to rewrite `std::pow` into `p[18]::pow`.
* A symbol named after a C++ keyword compiles, `default` and `int` included.
* A symbol named after a Python keyword is rejected and named, rather than
  silently renamed. SymPy parses through Python's parser.
* The `double` locals of a root event's `G_tt` lambda are named by position.
* `cppde::value_of(x)` is the value accessor across arithmetic types, both dual
  orders and their expression templates.

# cppDE 0.9.2

* A root event whose crossing falls exactly on an evaluated time now fires.
  Detection and bisection both tested the sign product strictly.
* A root event no longer fires again on the crossing it just handled.
* A fixed event that makes a root condition true now fires it, with the resets
  riding on the jump's surface. Terminal conditions are excluded.
* The step size is re-estimated after every event, not only for the multistep
  methods.
* `funCpp()` substitutes all symbols in one pass. A parameter named after a
  generated array rewrote the slots already emitted.
* `compile()` no longer repeats the OpenMP and KLU flags that the constructors
  already recorded on the model.
* `inst/examples/example_saltation.R` checks the sensitivity transport against
  a SymPy solution to first and second order, over five models.
* Generated entry points are dispatched by name and shared object rather than
  by a cached address, so loading and unloading are without consequence. A
  model whose library is gone names what it is missing.
* Scoping the lookup to one shared object also stops two models sharing an
  entry point name from reaching into each other.
* `clearNativeSymbols()` drops the remembered name pairings. It is no longer
  needed after loading or unloading.

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
