# Build the Optional Solver Libraries from Source

Downloads pinned upstream releases of
[SuiteSparse](https://github.com/DrTimothyAldenDavis/SuiteSparse) and
[SUNDIALS](https://github.com/LLNL/sundials), builds them, and installs
them into a per-user cache directory. This is the alternative to
installing them as system packages, for machines without administrator
rights or with versions too old to be useful.

## Usage

``` r
installLibs(
  which = c("sundials", "suitesparse"),
  dir = NULL,
  sundials_version = "7.9.0",
  suitesparse_version = "7.14.1",
  quiet = FALSE,
  ask = interactive()
)
```

## Arguments

- which:

  Which libraries to build. Default `"sundials"`, which installs
  SUNDIALS and SuiteSparse; `"suitesparse"` installs SuiteSparse/KLU
  only.

- dir:

  Directory to install into. Defaults to
  `tools::R_user_dir("cppDE", "cache")`. A version-specific subdirectory
  is created inside it.

- sundials_version, suitesparse_version:

  Upstream release tags. The defaults are a pair verified to build
  together.

- quiet:

  If `TRUE`, suppress compiler and CMake output. Default `FALSE`.

- ask:

  If `TRUE` (the default in an interactive session), ask for
  confirmation before downloading and writing to `dir`.

## Value

The install prefix, invisibly.

## Details

The two libraries serve different features and can be installed
independently:

- `"suitesparse"` installs SuiteSparse/KLU, which the sparse Jacobian
  path of
  [`cppODE()`](https://dmodverse.github.io/cppDE/reference/cppODE.md)
  and [`cvode()`](https://dmodverse.github.io/cppDE/reference/cvode.md)
  requires. It does not involve SUNDIALS.

- `"sundials"` installs SUNDIALS/CVODES for the
  [`cvode()`](https://dmodverse.github.io/cppDE/reference/cvode.md)
  backend. SUNDIALS links the sparse solver wrappers against an external
  KLU when SUNDIALS is built, so this installs SuiteSparse as well.

The build needs network access and runs for several minutes. It is never
started as part of installing cppDE unless explicitly requested; see the
section below.

cppDE records library paths at install time, so **cppDE must be
re-installed** for a build to take effect. The function prints the
required command when it finishes.

## Building during installation

To avoid installing cppDE twice, set the corresponding environment
variable before installing and repeat the same install command.
`configure` then performs the build and links against the result
immediately:

    Sys.setenv(CPPDE_BUILD_SUNDIALS = 1)      # or CPPDE_BUILD_SUITESPARSE
    devtools::install()

The variable is only needed for the build itself. `configure` scans the
cache on every run, so later re-installs and upgrades find an existing
build without it.

The trade-off is that R holds a lock on the package library for the
duration of the build. Interrupting it leaves a stale `00LOCK-cppDE`
which the next install refuses to overwrite; remove it with
`unlink(file.path(.libPaths()[1], "00LOCK-cppDE"), recursive = TRUE)`.
Calling this function first and installing afterwards avoids that, since
only the short install then holds the lock. Passing `args = "--no-lock"`
to the install avoids it as well, at the cost of the backup R would
otherwise keep of the previous version.

An interrupted build is safe to repeat: each library is marked complete
only after it has been installed and verified, so a half-finished one is
neither used nor kept, and a rerun resumes with what is missing.

## Removing the libraries

Everything lives in one directory; no system paths are written and no
system packages are touched. Deleting it is the complete uninstall:

    unlink(tools::R_user_dir("cppDE", "cache"), recursive = TRUE)

Re-install cppDE afterwards. Until then the installed package still
points at the deleted directory and affected models fail to load.

## What gets built

Only the components cppDE calls: the `klu` module of SuiteSparse with
the `amd`, `colamd` and `btf` dependencies, and CVODES plus the dense
and KLU sparse linear solvers of SUNDIALS. MPI and OpenMP are disabled
and no examples are built. The prefix is keyed by the version pair, so
changing either version produces a separate build.

Requires `cmake` (\>= 3.18) and one of `curl`, `wget` or `git`.

## Custom build options

`CPPDE_SUITESPARSE_CMAKE_ARGS` and `CPPDE_SUNDIALS_CMAKE_ARGS` are
appended to the respective `cmake` configure step and therefore override
the defaults, for instance to select a compiler:

    CPPDE_SUNDIALS_CMAKE_ARGS="-DCMAKE_C_COMPILER=icx" \
      CPPDE_BUILD_SUNDIALS=1 R CMD INSTALL .

## BLAS and LAPACK

SUNDIALS is built against R's BLAS and LAPACK; there is no option to
choose another. When R reports a BLAS, the build enables
`SUNLinSol_LapackDense`, which
[`cvode()`](https://dmodverse.github.io/cppDE/reference/cvode.md) then
uses for dense Jacobians; otherwise
[`cvode()`](https://dmodverse.github.io/cppDE/reference/cvode.md) uses
`SUNLinSol_Dense`. During a solve and across `fork()` cppDE runs BLAS on
one thread, see
[`forkGuard()`](https://dmodverse.github.io/cppDE/reference/forkGuard.md).

## See also

[`cvode()`](https://dmodverse.github.io/cppDE/reference/cvode.md) and
the `sparse` argument of
[`cppODE()`](https://dmodverse.github.io/cppDE/reference/cppODE.md),
which report the relevant option when a library is unavailable.

## Examples

``` r
if (FALSE) { # \dontrun{
## Downloads and builds SuiteSparse/KLU into the per-user cache
installLibs("suitesparse")
} # }
```
