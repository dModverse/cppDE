# Package index

## Models

Generate and compile the solver of an ODE model, natively or against
SUNDIALS CVODES.

- [`cppODE()`](https://dmodverse.github.io/cppDE/reference/cppODE.md) :
  Generate and Compile a C++ ODE Solver
- [`cvode()`](https://dmodverse.github.io/cppDE/reference/cvode.md) :
  Generate and Compile an ODE Solver Linked Against SUNDIALS CVODE(S)
- [`expressions`](https://dmodverse.github.io/cppDE/reference/expressions.md)
  : Expressions in Model Equations
- [`compile()`](https://dmodverse.github.io/cppDE/reference/compile.md)
  : Compile Generated C++ Model Code

## Solving

Solve a compiled model with sensitivities, read forcings and solver
statistics.

- [`solveODE()`](https://dmodverse.github.io/cppDE/reference/solveODE.md)
  : Run a Compiled ODE Model
- [`forcingValues()`](https://dmodverse.github.io/cppDE/reference/forcingValues.md)
  : Values of Forcings
- [`diagnostics()`](https://dmodverse.github.io/cppDE/reference/diagnostics.md)
  : Print Solver Diagnostics

## Reverse mode

Control the backward sweep and free the checkpoints of a kept forward
pass.

- [`adjointControl()`](https://dmodverse.github.io/cppDE/reference/adjointControl.md)
  : Controls of a Reverse Solve
- [`releaseStore()`](https://dmodverse.github.io/cppDE/reference/releaseStore.md)
  : Release the Checkpoints of a Reverse Solve

## Batches of conditions

Many conditions of one model in one call, in parallel.

- [`solveODEBatch()`](https://dmodverse.github.io/cppDE/reference/solveODEBatch.md)
  : Solve Many Conditions in One Call
- [`prepareBatch()`](https://dmodverse.github.io/cppDE/reference/prepareBatch.md)
  : Prepare a Batch for Repeated Solving
- [`solveBatch()`](https://dmodverse.github.io/cppDE/reference/solveBatch.md)
  : Solve a Prepared Batch

## Algebraic functions and steady states

Compiled algebraic functions with their derivatives, and steady states
by pseudo-transient continuation.

- [`cppFUN()`](https://dmodverse.github.io/cppDE/reference/cppFUN.md) :
  Compile Algebraic Functions with Optional Derivatives
- [`evaluateMany()`](https://dmodverse.github.io/cppDE/reference/evaluateMany.md)
  : Evaluate Several Compiled Functions at Once
- [`ptc()`](https://dmodverse.github.io/cppDE/reference/ptc.md) : Steady
  State by Pseudo-Transient Continuation

## Setup

Check the toolchain and the optional libraries, build them, and manage
threads and native symbols.

- [`codegenAvailable()`](https://dmodverse.github.io/cppDE/reference/codegenAvailable.md)
  : Can Code Be Generated and Compiled Here?
- [`batchAvailable()`](https://dmodverse.github.io/cppDE/reference/batchAvailable.md)
  : Will a Batch Run in Parallel?
- [`forkGuard()`](https://dmodverse.github.io/cppDE/reference/forkGuard.md)
  : Report the BLAS Fork Guard
- [`installLibs()`](https://dmodverse.github.io/cppDE/reference/installLibs.md)
  : Build the Optional Solver Libraries from Source
- [`clearNativeSymbols()`](https://dmodverse.github.io/cppDE/reference/clearNativeSymbols.md)
  : Forget Which Shared Object Exports Which Entry Point
