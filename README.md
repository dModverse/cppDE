
<!-- README.md is generated from README.Rmd. Please edit that file -->

# cppDE

<!-- badges: start -->

[![R-CMD-check](https://github.com/simonbeyer1/cppDE/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/simonbeyer1/cppDE/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

**Automated C++ code generation for ODE integration with sensitivity
calculation via in-tree forward-mode dual numbers.**

cppDE generates, compiles, and runs C++ solvers for ODE systems, with
optional first- and second-order parameter sensitivities computed
automatically via forward-mode automatic differentiation.

## Features

- **BDF/NDF and Adams multistep solvers** (variable-order 1–5,
  variable-step, Nordsieck form) adapted from
  [SUNDIALS/CVODE](https://computing.llnl.gov/projects/sundials)
- **Rosenbrock4 solver** extended from
  [Boost.Odeint](https://www.boost.org/doc/libs/release/libs/numeric/odeint/doc/html/index.html)
- **Tsit5 solver**: 7-stage explicit FSAL Runge–Kutta (Tsitouras
  2011) with embedded error estimate and dense output
- **Automatic differentiation**: in-tree dual numbers for exact first-
  and second-order parameter sensitivities in a single integration pass
- **Event handling**: time-based and root-triggered events with
  saltation matrix corrections for first- and second-order sensitivity
  continuity across events
- **Symbolic Jacobian**: automatically derived and compiled as analytic
  C++ code
- **Sparse support**: auto-detection of Jacobian sparsity with KLU
  factorization for large systems
- **LAPACK-backed dense solver** for `cvode()` when SUNDIALS provides
  it, routing the factorization through the BLAS that R is linked
  against
- **Batched solving**: `solveODEBatch()` integrates many conditions – or
  the subjects of a mixed-effects fit – in one `.Call`, over OpenMP
  threads, avoiding the fork and result serialization that a
  `parallel::mclapply()` loop pays per evaluation

## Installation

``` r
install.packages("devtools")
devtools::install_github("simonbeyer1/cppDE")
```

### Optional system dependencies

The core package installs without any system libraries. Three features
are gated on optional dependencies, detected by `./configure`
(Unix/macOS) or `./configure.win` (Windows) at install time:

- **`cvode()` backend** requires **SUNDIALS** (\>= 6.0).
- **Sparse Jacobian path** of both `cppODE()` and `cvode()` requires
  **SuiteSparse / KLU**.
- **Parallel `solveODEBatch()`** requires a toolchain with **OpenMP**.
  Nothing to install on Linux or Rtools; on macOS, Apple’s clang needs
  `libomp` (`brew install libomp`). `configure` links a test binary
  rather than trusting the flag, because macOS reports a non-empty
  `SHLIB_OPENMP_CXXFLAGS` even when libomp is absent.

They are independent: KLU is not part of SUNDIALS, and `cppODE()` calls
KLU without involving SUNDIALS at all. A missing library disables only
its own feature; the install succeeds either way. Without OpenMP,
`solveODEBatch()` still returns correct results – it just runs the
conditions serially.

Five ways to provide them:

| \# | Route | What you do |
|:---|:---|:---|
| 1 | Nothing | `cvode()` and `sparse = TRUE` report how to enable themselves. |
| 2 | System packages | Install them, then install cppDE. |
| 3 | Build during install | `CPPDE_BUILD_SUNDIALS=1` (or `CPPDE_BUILD_SUITESPARSE=1`), then install as usual. |
| 4 | Build from R | `cppDE::install_libs("sundials")`, then install cppDE. |
| 5 | Existing build | `CPPDE_SUNDIALS_HOME=<prefix>`, or the four `CPPDE_CVODE_*` flags. |

Routes 3 and 4 need no administrator rights and download nothing unless
asked. `./configure` also searches the system directories, `~/.local`,
`SUNDIALS_ROOT` / `SUNDIALS_DIR`, and the cache written by routes 3 and
4, so later re-installs need no variable.

#### System packages

| Platform | SUNDIALS | SuiteSparse / KLU |
|:---|:---|:---|
| Debian/Ubuntu | `sudo apt install libsundials-dev` | `sudo apt install libsuitesparse-dev` |
| Fedora | `sudo dnf install sundials-devel` | `sudo dnf install suitesparse-devel` |
| Arch | `sudo pacman -S sundials` | `sudo pacman -S suitesparse` |
| macOS (Homebrew) | `brew install sundials` | `brew install suite-sparse` |

On **Windows**, call Rtools’ bundled `pacman` from any shell,
substituting your Rtools version for `<ver>`:

``` powershell
C:\rtools<ver>\usr\bin\pacman.exe -Sy --noconfirm `
    mingw-w64-ucrt-x86_64-sundials mingw-w64-ucrt-x86_64-suitesparse
```

#### Building from source

For machines without administrator rights or with versions too old to be
useful. Set the flag and install as usual:

``` r
Sys.setenv(CPPDE_BUILD_SUNDIALS = 1)   # or CPPDE_BUILD_SUITESPARSE
devtools::install()
```

or run the build separately:

``` r
cppDE::install_libs("sundials")        # or "suitesparse"
```

`CPPDE_BUILD_SUNDIALS` covers SuiteSparse too, because SUNDIALS links
its sparse wrappers against an external KLU. Requires `cmake` and one of
`curl`, `wget`, `git`. The libraries land in
`tools::R_user_dir("cppDE", "cache")`; the flag is needed once, later
installs find them by themselves. To remove them:

``` r
unlink(tools::R_user_dir("cppDE", "cache"), recursive = TRUE)
```

`CPPDE_SUITESPARSE_CMAKE_ARGS` and `CPPDE_SUNDIALS_CMAKE_ARGS` are
appended to the respective `cmake` call for different build options.
SUNDIALS is built against the BLAS and LAPACK that R reports, which is
also what enables its LAPACK-backed dense solver.

#### An existing build

``` sh
CPPDE_SUNDIALS_HOME=$HOME/opt/sundials-7.4.0 R CMD INSTALL .
```

Point it at the `CMAKE_INSTALL_PREFIX`, i.e. the directory containing
`include/` and `lib/` (or `lib64/`), **not** a build tree, which has no
usable `include/` layout. Everything else is derived and verified by a
test compile. Not implemented on Windows.

<details>

<summary>

<b>Expert: hand-written compiler and linker flags</b>
</summary>

For layouts the prefix scan cannot express, split include/lib trees,
static linking, extra `-l` entries. Setting any of these skips the
corresponding probe; the flags are still verified by a test compile.

| Variable | Purpose |
|:---|:---|
| `CPPDE_CVODE_CFLAGS` | `-I` flags for the SUNDIALS headers |
| `CPPDE_CVODE_LIBS` | full link line for SUNDIALS |
| `CPPDE_CVODE_KLU_CFLAGS` | `-I` flags for `klu.h` |
| `CPPDE_CVODE_KLU_LIBS` | full link line for KLU + the SUNDIALS sparse wrappers |

``` sh
PREFIX=$HOME/opt/sundials-7.4.0
LIBDIR=$PREFIX/lib                 # $PREFIX/lib64 on Fedora/RHEL/SUSE

export CPPDE_CVODE_CFLAGS="-I$PREFIX/include"
export CPPDE_CVODE_LIBS="-L$LIBDIR -Wl,-rpath,$LIBDIR \
  -lsundials_cvodes -lsundials_nvecserial \
  -lsundials_sunmatrixdense -lsundials_sunlinsoldense -lsundials_core"
export CPPDE_CVODE_KLU_CFLAGS="-I/usr/include/suitesparse"
export CPPDE_CVODE_KLU_LIBS="-L$LIBDIR -Wl,-rpath,$LIBDIR \
  -lsundials_sunmatrixsparse -lsundials_sunlinsolklu -lklu"
```

You own the whole line, including `-Wl,-rpath`: without it linking
succeeds but `dyn.load()` fails later. `-lsundials_core` must come last
and exists only in SUNDIALS \>= 7. The `KLU_*` variables may point
elsewhere: a source-built SUNDIALS against the distribution’s
SuiteSparse is a normal combination.

</details>

<details>

<summary>

<b>Troubleshooting</b>
</summary>

``` sh
Rscript -e 'writeLines(readLines(system.file("cvodeConfig.dcf", package = "cppDE")))'
```

| Symptom | Cause |
|:---|:---|
| `cvodes/cvodes.h: No such file` | wrong prefix, or `CPPDE_CVODE_CFLAGS` unset |
| `sundials/sundials_config.h: No such file` | headers from a git checkout; `cmake --install` was not run |
| `cannot find -lsundials_core` | `lib` vs `lib64` mismatch |
| `undefined symbol: SUNContext_Create` | `-lsundials_core` missing or not last |
| `cannot open shared object file` at load | `-Wl,-rpath` missing |
| `sparse = TRUE` fails to link, dense works | SUNDIALS built without `-DENABLE_KLU=ON` |

`configure` prints the failing compile command when a library is present
but its probe failed. To see the probes regardless:

``` sh
CPPDE_CONFIGURE_VERBOSE=1 R CMD INSTALL .
```

</details>

<details>

<summary>

<b>Compiler flags for generated models</b>
</summary>

The variables above are read by `./configure` at install time. One
further variable is read at model-compile time: `CPPDE_EXTRA_CXXFLAGS`
is appended to the flags of every model `compile()` builds in the
session.

``` r
Sys.setenv(CPPDE_EXTRA_CXXFLAGS = "-march=native")
```

Exporting `PKG_CXXFLAGS` instead has no effect, `compile()` sets that
variable itself, from each model’s `compileArgs` attribute.

</details>

## License

cppDE is distributed under the **MIT License**; see `LICENSE` /
`LICENSE.md` for the full text.

Portions of the C++ solver core in `inst/include/cppde/` are derived
from third-party projects and remain subject to their upstream licenses:

- The Nordsieck multistep stepper (BDF/NDF and Adams-Moulton) in
  `cppde_multistepper.hpp`, `cppde_multistepper_controller.hpp` and
  `cppde_newton.hpp` is a port of the corresponding routines in
  [SUNDIALS/CVODE(S)](https://computing.llnl.gov/projects/sundials)
  (Copyright © 2002–2024 Lawrence Livermore National Security and
  Southern Methodist University), distributed under the **BSD-3-Clause**
  license.
- The Rosenbrock4 stepper architecture in `cppde_rosenbrock4.hpp` (and
  the surrounding stage-matrix / dense-output / stepper-protocol
  infrastructure) is derived from
  [Boost.Numeric.Odeint](https://www.boost.org/doc/libs/release/libs/numeric/odeint/)
  by Karsten Ahnert, Mario Mulansky and Christoph Koke, distributed
  under the **Boost Software License 1.0**.

Full upstream copyright notices and license texts are reproduced in
[`inst/COPYRIGHTS`](inst/COPYRIGHTS). Both upstream licenses are
permissive and compatible with MIT; using or redistributing cppDE
requires preserving the notices in `inst/COPYRIGHTS` and the relevant
header attribution blocks.

The optional system libraries linked at runtime are not vendored:
SUNDIALS (BSD-3-Clause) for the `cvode()` backend and SuiteSparse-KLU
(LGPL-2.1+) for the sparse-Jacobian path. Their own licenses govern the
libraries at the user’s installation site.
