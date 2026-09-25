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
