# Installation and optional libraries

## Package installation

``` r

# install.packages("remotes")
remotes::install_github("dModverse/cppDE")
```

Code generation at run time needs a C++17 compiler (on Windows
[Rtools](https://cran.r-project.org/bin/windows/Rtools/)) and Python
with SymPy, which reticulate provides on first use if none is
configured.
[`codegenAvailable()`](https://dmodverse.github.io/cppDE/reference/codegenAvailable.md)
checks both without installing anything, `codegenAvailable("cvode")`
also SUNDIALS.

``` r

library(cppDE)
codegenAvailable()
#> [1] TRUE
codegenAvailable("cvode")
#> [1] TRUE
```

## Optional system libraries

Three independent features need a system library, detected by
`./configure` at install time:

| Feature | Needs | Without it |
|:---|:---|:---|
| [`cvode()`](https://dmodverse.github.io/cppDE/reference/cvode.md) backend | SUNDIALS (\>= 6.0) | [`cvode()`](https://dmodverse.github.io/cppDE/reference/cvode.md) stops with instructions; native solvers unaffected |
| sparse Jacobians | SuiteSparse / KLU | dense LU; `sparse = TRUE` stops with instructions |
| parallel [`solveODEBatch()`](https://dmodverse.github.io/cppDE/reference/solveODEBatch.md) | OpenMP toolchain | same results, conditions solved in sequence |

[`cppODE()`](https://dmodverse.github.io/cppDE/reference/cppODE.md) uses
KLU without SUNDIALS.

| Platform | SUNDIALS | SuiteSparse / KLU |
|:---|:---|:---|
| Debian/Ubuntu | `sudo apt install libsundials-dev` | `sudo apt install libsuitesparse-dev` |
| Fedora | `sudo dnf install sundials-devel` | `sudo dnf install suitesparse-devel` |
| Arch | `sudo pacman -S sundials` | `sudo pacman -S suitesparse` |
| macOS (Homebrew) | `brew install sundials` | `brew install suite-sparse` |

### Library build by installLibs()

[`installLibs()`](https://dmodverse.github.io/cppDE/reference/installLibs.md)
builds pinned releases into a per-user cache, for missing administrator
rights or outdated system versions:

``` r

cppDE::installLibs("sundials")         # SUNDIALS and SuiteSparse
cppDE::installLibs("suitesparse")      # SuiteSparse/KLU only
```

The build needs `cmake` (\>= 3.18) and network access. cppDE records
library paths at install time and is then reinstalled
([`installLibs()`](https://dmodverse.github.io/cppDE/reference/installLibs.md)
prints the command); `./configure` finds the cache on every later
install.
[`?installLibs`](https://dmodverse.github.io/cppDE/reference/installLibs.md)
covers `CPPDE_BUILD_SUNDIALS`, `CPPDE_BUILD_SUITESPARSE`, CMake options
and removal of the cache.

### OpenMP

`./configure` enables OpenMP with R’s `SHLIB_OPENMP_CXXFLAGS` if a test
program links. Apple’s clang needs `libomp` (`brew install libomp`);
Linux compilers and Rtools have OpenMP. Without it,
[`batchAvailable()`](https://dmodverse.github.io/cppDE/reference/batchAvailable.md)
is `FALSE`.

### Detected configuration

`./configure` writes the result to `cvodeConfig.dcf`:

``` r

writeLines(readLines(system.file("cvodeConfig.dcf", package = "cppDE")))
#> available: TRUE
#> cflags: -I/usr/lib/x86_64-linux-gnu/openmpi/include -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi
#> libs: -lsundials_cvodes -lsundials_nvecserial -lsundials_sunmatrixdense -lsundials_sunlinsoldense -lsundials_core -L/usr/lib/x86_64-linux-gnu/openmpi/lib -lmpi
#> klu_available: TRUE
#> klu_cflags: -I/usr/include/suitesparse
#> klu_libs: -lsundials_sunmatrixsparse -lsundials_sunlinsolklu -lklu
#> cvode_lapack_available: FALSE
#> cvode_lapack_libs: 
#> openmp_available: TRUE
#> openmp_cxxflags: -fopenmp
#> openmp_libs: -fopenmp
```

`available`, `klu_available` and `openmp_available` report the three
features, the other fields the compile and link flags of generated
models.

### Compiler flags

[`compile()`](https://dmodverse.github.io/cppDE/reference/compile.md)
appends `CPPDE_EXTRA_CXXFLAGS` and `CPPDE_EXTRA_LIBS` to the compiler
and linker flags of every model, for instance for `-march` or a
sanitizer:

``` r

Sys.setenv(CPPDE_EXTRA_CXXFLAGS = "-march=native")
```

`PKG_CXXFLAGS` has no effect, since
[`compile()`](https://dmodverse.github.io/cppDE/reference/compile.md)
sets it.

## Failures and diagnosis

### Errors before the solve

| Error | Cause |
|:---|:---|
| `Model not loaded. Run compile() first.` | generated with `compile = FALSE`, not yet compiled |
| `Compilation failed:` and compiler output | no C++17 compiler, or an expression that does not compile; `compile(verbose = TRUE)` shows the commands |
| `'parms' missing: ...` | every name of `attr(model, "variables")` and `attr(model, "parameters")` is needed |
| `The CVODE backend was disabled at install time ...` | SUNDIALS not found at install time |
| `Sparse Jacobian requested but the KLU linear solver was not available ...` | KLU not found at install time |

``` r

eqns  <- c(A = "-k1*A", B = "k1*A - k2*B")
parms <- c(A = 1, B = 0, k1 = 0.5, k2 = 0.2)

pending <- cppODE(eqns, modelname = "pending", compile = FALSE)
solveODE(pending, 0:10, parms)
#> Error:
#> ! Model not loaded. Run compile() first.
```

### Solver failures

The solver stops after `maxsteps` steps (code `-1`) or `maxattempts`
consecutive attempts without progress (`-4`).
[`solveODE()`](https://dmodverse.github.io/cppDE/reference/solveODE.md)
then stops with the solver message; `onFailure = "warn"` or `"silent"`
returns the partial result.
[`diagnostics()`](https://dmodverse.github.io/cppDE/reference/diagnostics.md)
prints the solver statistics:

``` r

traced <- cppODE(eqns, modelname = "traced", stepTrace = TRUE)
solveODE(traced, 0:10, parms, maxsteps = 5)
#> Error:
#> ! Solver did not complete: Maximum number of steps exceeded (5).
#>   Reached t = 0 (1 of 11 time points).

out <- solveODE(traced, 0:10, parms, maxsteps = 5, onFailure = "warn")
#> Warning: Solver did not complete: Maximum number of steps exceeded (5).
#>   Reached t = 0 (1 of 11 time points).
#>   Returning partial results.
diagnostics(out)
#> ----------------------  NDF solver statistics  ----------------------
#>   Return code                  : -1
#>   Message                      : Too much work: maximum number of integration steps exceeded (CV_TOO_MUCH_WORK).
#> ---------------------------------------------------------------------
#>   Accepted steps               : 6
#>   Rejected steps               : 0
#>   Function evaluations         : 7
#>   Jacobian evaluations         : 1
#>   LU factorizations            : 5
#>   Last step size (successful)  : 0.00171081
#>   Last method order            : 2
#>   Time reached                 : 0
#> ---------------------------------------------------------------------
```

### Step trace and timing profile

With `stepTrace = TRUE`, `$trace` has one row per step attempt (columns
in
[`?solveODE`](https://dmodverse.github.io/cppDE/reference/solveODE.md));
`traceFile` writes it to CSV.

``` r

head(out$trace[, c("nst", "t", "h", "q", "dsm", "newton_conv", "setup_reason")])
#>   nst            t            h q          dsm newton_conv setup_reason
#> 1   0 1.085376e-08 1.085376e-08 1 3.458854e-11           1       no_jac
#> 2   1 1.085484e-04 1.085376e-04 1 3.458728e-03           1        DGMAX
#> 3   2 1.193924e-03 1.085376e-03 2 9.747118e-02           1        DGMAX
#> 4   3 2.279300e-03 1.085376e-03 2 4.255788e-02           1             
#> 5   4 3.990112e-03 1.710812e-03 2 2.779881e-02           1        DGMAX
#> 6   5 7.098073e-03 3.107962e-03 2 1.760033e-02           1        DGMAX
```

With `profile = TRUE`, each solve prints the time per category to
standard error. Both options slow every step.
