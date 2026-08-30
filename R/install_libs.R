#' Build the optional solver libraries from source
#'
#' Downloads pinned upstream releases of
#' [SuiteSparse](https://github.com/DrTimothyAldenDavis/SuiteSparse) and
#' [SUNDIALS](https://github.com/LLNL/sundials), builds them, and installs
#' them into a per-user cache directory. This is the alternative to
#' installing them as system packages, for machines without administrator
#' rights or with versions too old to be useful.
#'
#' The two libraries serve different features and can be installed
#' independently:
#'
#' * `"suitesparse"` installs SuiteSparse/KLU, which the sparse Jacobian
#'   path of [cppODE()] and [cvode()] requires. It does not involve
#'   SUNDIALS.
#' * `"sundials"` installs SUNDIALS/CVODES for the [cvode()] backend.
#'   SUNDIALS links its sparse solver wrappers against an external KLU at
#'   its own build time, so this installs SuiteSparse as well.
#'
#' The build needs network access and runs for several minutes. It is never
#' started as part of installing cppDE unless explicitly requested; see
#' the section below.
#'
#' cppDE records library paths at install time, so **cppDE must be
#' re-installed** for a build to take effect. The function prints the
#' required command when it finishes.
#'
#' @section Building during installation:
#' To avoid installing cppDE twice, set the corresponding environment
#' variable before installing and repeat the same install command.
#' `configure` then performs the build and links against the result
#' immediately:
#'
#' ```r
#' Sys.setenv(CPPDE_BUILD_SUNDIALS = 1)      # or CPPDE_BUILD_SUITESPARSE
#' devtools::install()
#' ```
#'
#' The variable is only needed for the build itself. `configure` scans the
#' cache on every run, so later re-installs and upgrades find an existing
#' build without it.
#'
#' The trade-off is that R holds a lock on the package library for the
#' duration of the build. Interrupting it leaves a stale `00LOCK-cppDE`
#' which the next install refuses to overwrite; remove it with
#' `unlink(file.path(.libPaths()[1], "00LOCK-cppDE"), recursive = TRUE)`.
#' Calling this function first and installing afterwards avoids that, since
#' only the short install then holds the lock. Passing
#' `args = "--no-lock"` to the install avoids it as well, at the cost of
#' the backup R would otherwise keep of the previous version.
#'
#' An interrupted build is safe to repeat: each library is marked complete
#' only after it has been installed and verified, so a half-finished one is
#' neither used nor kept, and a rerun resumes with what is missing.
#'
#' @section Removing the libraries:
#' Everything lives in one directory; no system paths are written and no
#' system packages are touched. Deleting it is the complete uninstall:
#'
#' ```r
#' unlink(tools::R_user_dir("cppDE", "cache"), recursive = TRUE)
#' ```
#'
#' Re-install cppDE afterwards. Until then the installed package still
#' points at the deleted directory and affected models fail to load.
#'
#' @section What gets built:
#' Only the components cppDE calls: the `klu` module of SuiteSparse with
#' its `amd`, `colamd` and `btf` dependencies, and CVODES plus the dense
#' and KLU sparse linear solvers of SUNDIALS. MPI and OpenMP are disabled
#' and no examples are built. The prefix is keyed by the version pair, so
#' changing either version produces a separate build.
#'
#' Requires `cmake` (>= 3.18) and one of `curl`, `wget` or `git`.
#'
#' @section Custom build options:
#' `CPPDE_SUITESPARSE_CMAKE_ARGS` and `CPPDE_SUNDIALS_CMAKE_ARGS` are
#' appended to the respective `cmake` configure step and therefore override
#' the defaults, for instance to select a compiler:
#'
#' ```sh
#' CPPDE_SUNDIALS_CMAKE_ARGS="-DCMAKE_C_COMPILER=icx" \
#'   CPPDE_BUILD_SUNDIALS=1 R CMD INSTALL .
#' ```
#'
#' @section BLAS and LAPACK:
#' There is no option for selecting a BLAS: everything follows R. The
#' `cppODE()` backend calls `dgetrf`/`dgetrs` through R's own BLAS and
#' LAPACK, and SUNDIALS is built against the same libraries R reports, so a
#' compiled model links one implementation rather than two. Whatever R uses
#' -- reference LAPACK, OpenBLAS, MKL via FlexiBLAS -- is what cppDE uses.
#' cppDE sets that library to a single thread, because its factorisations
#' are small and frequent and thread startup costs more than it saves.
#'
#' This enables SUNDIALS' `SUNLinSol_LapackDense`, which [cvode()] then
#' uses for dense Jacobians. If R reports no BLAS, or the LAPACK-enabled
#' configure step fails, the build continues without that solver and
#' `cvode()` uses `SUNLinSol_Dense`. There is no per-model switch: the
#' choice is fixed when cppDE is installed.
#'
#' The SuiteSparse modules built here (`klu`, `amd`, `colamd`, `btf`) do
#' not use BLAS at all.
#'
#' @param which Which libraries to build. `"sundials"` (the default)
#'   installs SUNDIALS and SuiteSparse; `"suitesparse"` installs
#'   SuiteSparse/KLU only.
#' @param dir Directory to install into. Defaults to
#'   `tools::R_user_dir("cppDE", "cache")`. A version-specific
#'   subdirectory is created inside it.
#' @param sundials_version,suitesparse_version Upstream release tags. The
#'   defaults are a pair verified to build together.
#' @param quiet If `TRUE`, suppress compiler and CMake output.
#' @param ask If `TRUE` (the default in an interactive session), ask for
#'   confirmation before downloading and writing to `dir`.
#'
#' @return The install prefix, invisibly.
#'
#' @seealso [cvode()] and the `sparse` argument of [cppODE()], which report
#'   the relevant option when a library is unavailable.
#'
#' @export
install_libs <- function(which = c("sundials", "suitesparse"),
                         dir = NULL,
                         sundials_version = "7.4.0",
                         suitesparse_version = "7.10.0",
                         quiet = FALSE,
                         ask = interactive()) {

  which <- match.arg(which)

  if (.Platform$OS.type == "windows")
    stop("install_libs() is not supported on Windows.\n",
         "  Use Rtools' package manager instead, substituting your Rtools\n",
         "  version for <ver>:\n",
         "    C:/rtools<ver>/usr/bin/pacman.exe -Sy --noconfirm \\\n",
         "      mingw-w64-ucrt-x86_64-sundials mingw-w64-ucrt-x86_64-suitesparse\n",
         "  Then: R CMD INSTALL <path/to/cppDE>", call. = FALSE)

  script <- system.file("tools", "build-libs.sh", package = "cppDE")
  if (!nzchar(script)) {
    # Source checkout, e.g. under devtools::load_all().
    script <- file.path("inst", "tools", "build-libs.sh")
    if (!file.exists(script))
      stop("cannot locate build-libs.sh; is cppDE installed correctly?",
           call. = FALSE)
  }

  if (is.null(dir)) dir <- tools::R_user_dir("cppDE", "cache")

  if (Sys.which("cmake") == "")
    stop("cmake was not found on the PATH but is required for the build.\n",
         "  Install it (e.g. 'sudo apt install cmake', 'sudo dnf install ",
         "cmake',\n  'brew install cmake') and retry.", call. = FALSE)
  if (all(Sys.which(c("curl", "wget", "git")) == ""))
    stop("need one of curl, wget or git to download the sources, but none ",
         "was found on the PATH.", call. = FALSE)

  if (isTRUE(ask)) {
    what <- if (which == "sundials") {
      paste0("  SuiteSparse ", suitesparse_version, " (klu, amd, colamd, btf)\n",
             "  SUNDIALS    ", sundials_version, " (cvodes, dense + KLU solvers)\n")
    } else {
      paste0("  SuiteSparse ", suitesparse_version, " (klu, amd, colamd, btf)\n")
    }
    msg <- paste0("This downloads and builds from source:\n", what,
                  "into: ", dir, "\n",
                  "It requires network access and takes several minutes.\n",
                  "Proceed?")
    if (!isTRUE(utils::askYesNo(msg, default = TRUE))) {
      message("Aborted; nothing was downloaded or written.")
      return(invisible(NULL))
    }
  }

  dir.create(dir, recursive = TRUE, showWarnings = FALSE)

  # The build itself lives in the shell script so that configure can run it
  # unchanged at install time.
  out <- system2(
    "sh", c(shQuote(normalizePath(script)), which),
    env = c(paste0("CPPDE_LIBS_CACHE=", shQuote(dir)),
            paste0("CPPDE_SUNDIALS_VERSION=", shQuote(sundials_version)),
            paste0("CPPDE_SUITESPARSE_VERSION=", shQuote(suitesparse_version))),
    stdout = TRUE,
    stderr = if (quiet) FALSE else ""
  )

  if (!is.null(attr(out, "status")) && attr(out, "status") != 0L)
    stop("the build failed; see the messages above for the failing step.",
         call. = FALSE)

  prefix <- utils::tail(out[nzchar(out)], 1L)
  if (length(prefix) != 1L || !dir.exists(prefix))
    stop("the build reported success but produced no usable prefix.",
         call. = FALSE)

  built <- if (which == "sundials") "SUNDIALS and SuiteSparse/KLU" else
    "SuiteSparse/KLU"

  # configure scans the default cache by itself; only a prefix outside it
  # has to be named explicitly on every re-install.
  in_default_cache <- identical(
    normalizePath(dirname(prefix), mustWork = FALSE),
    normalizePath(tools::R_user_dir("cppDE", "cache"), mustWork = FALSE))

  message("\n", built, " are ready in:\n  ", prefix,
          "\n\ncppDE records library paths at install time. Re-install it ",
          "to enable\nthe affected features:\n\n",
          if (in_default_cache) {
            paste0("  devtools::install()\n\n",
                   "No environment variable is needed: configure finds this ",
                   "build\nautomatically, now and on later re-installs.")
          } else {
            paste0("  Sys.setenv(CPPDE_SUNDIALS_HOME = \"", prefix, "\")\n",
                   "  devtools::install()\n\n",
                   "This prefix is outside the default cache, so the variable ",
                   "is\nrequired on every re-install. Add it to ~/.Renviron to ",
                   "persist it.")
          },
          "\n\nTo remove it again:\n  unlink(\"", prefix, "\", recursive = TRUE)")

  invisible(prefix)
}
