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
