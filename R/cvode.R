#' Generate and Compile an ODE Solver Linked Against SUNDIALS CVODE(S)
#'
#' Generates a C++ ODE solver linked against SUNDIALS CVODES, from the system
#' or from [installLibs()], compiles it via `R CMD SHLIB`, and returns a handle
#' for use with [solveODE()]. Sensitivities are computed by the CVODES forward
#' sensitivity solver (`deriv = TRUE`) or by its adjoint
#' (`derivMode = "reverse"`). The compiled model exposes the same R interface
#' as a model from [cppODE()]; the differences between the two backends are
#' described in `vignette("Methods", package = "cppDE")`.
#'
#' Available methods are `"bdf"` (default) and `"adams"`. Forward
#' sensitivities are of first order only; `deriv2` is not supported. Events,
#' forcings, `rootfunc`, `fixed` and a switch on a state behave as in
#' [cppODE()]; the switch is located by CVODES' root finding.
#'
#' Needs SUNDIALS (>= 6.0) at install time, see [installLibs()]; otherwise
#' `cvode()` errors at the first call with platform-specific install hints.
#' KLU, required when `sparse = TRUE`, is detected the same way.
#'
#' @inheritParams cppODE
#' @param deriv Logical. Compute first-order forward sensitivities. Default
#'   `FALSE`.
#' @param compile Logical. Compile and load the generated C++ code. Default
#'   `TRUE`; with `FALSE`, compile several models together with [compile()].
#' @param includeTimeZero Logical. Ensure that `0` is part of the integration
#'   times, as [cppODE()] does. Default `TRUE`.
#' @param method One of `"bdf"` (default) or `"adams"`.
#' @param asaCheckpoints Number of forward steps between two checkpoints of the
#'   adjoint, under `derivMode = "reverse"`. Fewer steps take more memory and
#'   interpolate the forward state more accurately. Default `200`.
#' @param derivMode Direction the derivatives are taken in. `"forward"`
#'   (default) is the CVODES forward sensitivity solver, driven by `deriv`.
#'   `"reverse"` is the CVODES adjoint, one backward solve per cotangent
#'   column; it needs `deriv = FALSE`.
#' @param stepTrace Logical, default `FALSE`. Record per-step diagnostics,
#'   returned as `$trace` by [solveODE()]; with `events` or `rootfunc` one row
#'   per output time and event.
#'
#' @return The compiled model name (character) with the attributes of a
#'   model from [cppODE()], except that `deriv2` is `FALSE` and `useNDF` is
#'   `NA`, plus `backend = "cvode"` and `lapackDense`, whether dense Jacobians
#'   use `SUNLinSol_LapackDense` rather than `SUNLinSol_Dense`.
#'
#' @references
#' Hindmarsh, A. C., Brown, P. N., Grant, K. E., Lee, S. L., Serban, R.,
#' Shumaker, D. E., and Woodward, C. S. (2005). SUNDIALS: Suite of
#' Nonlinear and Differential/Algebraic Equation Solvers.
#' \emph{ACM Transactions on Mathematical Software} \strong{31}(3), 363-396.
#'
#' @seealso [cppODE()], [solveODE()], [cppFUN()];
#'   `vignette("Methods", package = "cppDE")`.
#' @example inst/examples/cvode.R
#' @export
cvode <- function(rhs, events = NULL, rootfunc = NULL, fixed = NULL, forcings = NULL,
                  compile = TRUE, modelname = NULL, outdir = tempdir(),
                  deriv = FALSE,
                  derivMode = c("forward", "reverse"),
                  asaCheckpoints = 200L,
                  sparse = NULL,
                  method = c("bdf", "adams"),
                  includeTimeZero = TRUE,
                  stepTrace = FALSE,
                  verbose = FALSE) {

  method <- match.arg(method)
  derivMode <- match.arg(derivMode)
  is_reverse <- identical(derivMode, "reverse")
  # The two directions are separate compilations, as they are on the native
  # backend: the direction decides what the generated code is.
  if (is_reverse && deriv)
    stop("derivMode = \"reverse\" computes no forward sensitivities; use ",
         "deriv = FALSE.", call. = FALSE)
  asaCheckpoints <- as.integer(asaCheckpoints)[1]
  if (is.na(asaCheckpoints) || asaCheckpoints < 1L)
    stop("'asaCheckpoints' must be a positive integer", call. = FALSE)

  # --- Availability check (populated by configure at install time) ---
  if (!isTRUE(cvodeConfig$available)) {
    stop(
      "The CVODE backend was disabled at install time because SUNDIALS ",
      "(>= 6.0) was not found on the build host.\n",
      "  cppODE()'s own solvers are unaffected.\n",
      "  Build SUNDIALS (and SuiteSparse/KLU) from source into a per-user\n",
      "  cache, no administrator rights required:\n",
      "      cppDE::installLibs(\"sundials\")\n",
      "  then run the re-install command it prints. Alternatively install\n",
      "  the SUNDIALS development headers system-wide and re-install:\n",
      "    Debian/Ubuntu : sudo apt install libsundials-dev\n",
      "    Fedora        : sudo dnf install sundials-devel\n",
      "    macOS (brew)  : brew install sundials\n",
      "    Windows       : from any shell (PowerShell / cmd / Git Bash)\n",
      "                    call Rtools' pacman by full path, substitute\n",
      "                    your installed version for <ver> (e.g. 44 or 45):\n",
      "                      C:/rtools<ver>/usr/bin/pacman.exe -Sy --noconfirm mingw-w64-ucrt-x86_64-sundials\n",
      "                    The .pc files land in C:/rtools<ver>/ucrt64/\n",
      "                    where the package's configure.win picks them up\n",
      "                    automatically on re-install.\n",
      "  Then: R CMD INSTALL <path/to/cppDE>",
      call. = FALSE)
  }

  # --- Normalize rhs (same as cppODE) ---
  rhs <- unclass(rhs)
  rhs <- gsub("\n", "", rhs)

  # Every expression the model is built from, as in cppODE().
  checkSymbolNames(rhs, forcings, fixed)
  if (!is.null(events)) {
    for (col in c("var", "value", "time", "root"))
      if (col %in% names(events) && is.character(events[[col]]))
        checkSymbolNames(events[[col]])
  }
  if (!is.null(rootfunc) && !identical(tolower(rootfunc), "equilibrate"))
    checkSymbolNames(rootfunc)

  variables <- names(rhs)
  if (is.null(variables) || any(!nzchar(variables)))
    stop("'rhs' must be a named character vector")

  # --- Identify parameters via getSymbols (same helper as cppODE) ---
  # Collect symbols from rhs and any event/rootfunc expressions too,
  # so params captures everything the generated code will reference.
  all_expressions <- rhs
  if (!is.null(events)) {
    bad <- which(!xor(!is.na(events$time), !is.na(events$root)))
    if (length(bad) > 0) {
      stop(sprintf(
        "Each event must define exactly one of 'time' or 'root'. Invalid event(s): %s",
        paste(bad, collapse = ", ")))
    }
    all_expressions <- c(all_expressions,
                         events$value,
                         if ("time" %in% names(events)) events$time,
                         if ("root" %in% names(events)) events$root)
  }
  if (!is.null(rootfunc) && !identical(tolower(rootfunc), "equilibrate")) {
    all_expressions <- c(all_expressions, rootfunc)
  }
  symbols <- getSymbols(all_expressions)

  if (is.null(forcings)) forcings <- character(0)
  if (length(forcings) > 0) {
    unknown_forcings <- setdiff(forcings, symbols)
    if (length(unknown_forcings) > 0)
      stop("Unknown forcing symbols: ", paste(unknown_forcings, collapse = ", "))
    forcing_states <- intersect(forcings, variables)
    if (length(forcing_states) > 0)
      stop("Forcing names cannot be state variables: ", paste(forcing_states, collapse = ", "))
  }

  params <- setdiff(symbols, c(variables, forcings, "time"))

  # --- Handle fixed ---
  if (is.null(fixed)) fixed <- character(0)
  fixed_initials <- if (deriv) intersect(fixed, variables) else character(0)
  fixed_params   <- if (deriv) intersect(fixed, params)    else character(0)
  sens_initials  <- if (deriv) setdiff(variables, fixed_initials) else character(0)
  sens_params    <- if (deriv) setdiff(params,    fixed_params)   else character(0)
  sens_names     <- c(sens_initials, sens_params)
  n_total_sens   <- length(sens_names)

  # --- Unique model name ---
  if (is.null(modelname)) {
    modelname <- randomModelname("c")
  }
  modelname <- unique_modelname(modelname)

  if (!dir.exists(outdir)) stop("outdir does not exist: ", outdir)

  # --- Early KLU check: explicit sparse = TRUE with no KLU is fatal ---
  if (isTRUE(sparse) && !isTRUE(cvodeConfig$klu_available)) {
    stop(
      "sparse = TRUE requested but the KLU linear solver was not available\n",
      "at install time.\n",
      "  Build SuiteSparse/KLU from source into a per-user cache, no\n",
      "  administrator rights required:\n",
      "      cppDE::installLibs(\"suitesparse\")\n",
      "  then run the re-install command it prints. Alternatively install\n",
      "  the SuiteSparse development headers system-wide and re-install:\n",
      "    Debian/Ubuntu : sudo apt install libsuitesparse-dev\n",
      "    Fedora        : sudo dnf install suitesparse-devel\n",
      "    macOS (brew)  : brew install suite-sparse\n",
      "    Windows       : from any shell (PowerShell / cmd / Git Bash)\n",
      "                    call Rtools' pacman by full path, substitute\n",
      "                    your installed version for <ver> (e.g. 44 or 45):\n",
      "                      C:/rtools<ver>/usr/bin/pacman.exe -Sy --noconfirm mingw-w64-ucrt-x86_64-suitesparse\n",
      "                    then re-run R CMD INSTALL <path/to/cppDE>",
      call. = FALSE)
  }
  # Auto-selected sparse without KLU -> force dense.
  sparse_for_codegen <- sparse
  if (is.null(sparse) && !isTRUE(cvodeConfig$klu_available)) {
    sparse_for_codegen <- FALSE
  }

  # The dense linear solver follows what ./configure found at install
  # time; SUNLinSol_LapackDense when SUNDIALS provides it, otherwise
  # SUNLinSol_Dense.
  lapack_for_codegen <- isTRUE(cvodeConfig$cvode_lapack_available)

  # --- Codegen ---
  codegen <- get_codegen_cvode_py()
  if (verbose) message("Generating CVODE C++ source...")

  res <- codegen$generate_cvode_cpp(
    rhs_dict = as.list(setNames(rhs, variables)),
    params_list = params,
    modelname = modelname,
    outdir = normalizePath(outdir, winslash = "/", mustWork = FALSE),
    deriv = deriv,
    reverse = is_reverse,
    asa_checkpoints = asaCheckpoints,
    fixed_states = fixed_initials,
    fixed_params = fixed_params,
    sparse = sparse_for_codegen,
    lapack = lapack_for_codegen,
    method = method,
    forcings_list = forcings,
    events = events,
    rootfunc = rootfunc,
    include_time_zero = includeTimeZero,
    version = as.character(utils::packageVersion("cppDE"))
  )

  use_sparse <- isTRUE(res$use_sparse)
  use_lapack <- isTRUE(res$use_lapack)
  if (use_sparse && verbose) {
    message(sprintf("  Sparse Jacobian detected (%d states, %d nnz)",
                    length(variables), length(res$jac_nnz_rows)))
  }

  # --- Attributes (mirror cppODE so solveODE works unchanged) ---
  attr(modelname, "equations")   <- rhs
  attr(modelname, "srcfile")     <- normalizePath(res$srcfile, winslash = "/", mustWork = FALSE)
  attr(modelname, "variables")   <- variables
  attr(modelname, "parameters")  <- params
  attr(modelname, "forcings")    <- forcings
  attr(modelname, "events")      <- events
  attr(modelname, "rootfunc")    <- rootfunc
  attr(modelname, "fixed")       <- c(fixed_initials, fixed_params)
  attr(modelname, "deriv")       <- isTRUE(deriv)
  attr(modelname, "deriv2")      <- FALSE
  attr(modelname, "derivMode")   <- derivMode
  attr(modelname, "sparse")      <- use_sparse
  attr(modelname, "lapackDense") <- use_lapack
  attr(modelname, "method")      <- method
  attr(modelname, "useNDF")      <- NA  # not meaningful for CVODE
  attr(modelname, "backend")     <- "cvode"

  # The sens dim defaults to model-parameter names; solveODE() overrides it per
  # call when the tangent has a full Phi' shape.
  attr(modelname, "dimNames") <- if (deriv) {
    list(time = "time", variable = variables, sens = sens_names)
  } else {
    list(time = "time", variable = variables)
  }

  # --- Compile args: codegen preprocessor defs (+ -DCVODE_KLU in sparse mode)
  # plus the SUNDIALS include path discovered at install time. Linker flags
  # come from `cvodeConfig` (populated by ./configure), not from codegen.
  compile_args <- c(unlist(res$compile_defs), cvodeConfig$cflags)
  link_libs    <- cvodeConfig$libs
  if (isTRUE(use_lapack)) {
    link_libs <- paste(cvodeConfig$cvode_lapack_libs, link_libs)
  }
  if (use_sparse) {
    compile_args <- c(compile_args, cvodeConfig$klu_cflags)
    link_libs    <- paste(link_libs, cvodeConfig$klu_libs)
  }
  if (isTRUE(stepTrace)) {
    compile_args <- c(compile_args, "-DCVODE_STEP_TRACE")
  }
  attr(modelname, "compileArgs") <- paste(compile_args[nzchar(compile_args)], collapse = " ")
  if (isTRUE(cvodeConfig$openmp_available)) {
    attr(modelname, "compileArgs") <- paste(attr(modelname, "compileArgs"),
                                            cvodeConfig$openmp_cxxflags)
    link_libs <- paste(link_libs, cvodeConfig$openmp_libs)
  }
  attr(modelname, "linkArgs")    <- link_libs

  if (compile) {
    compile(modelname, verbose = verbose)
  }
  modelname
}
