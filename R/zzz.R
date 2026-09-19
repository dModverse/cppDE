#' Package-internal CVODE backend configuration
#'
#' Populated at load time from `inst/cvodeConfig.dcf` (written by
#' `./configure` or `./configure.win` at install time).  Defaults mark
#' the backend as disabled so the package loads cleanly even if the DCF
#' is missing; in that case `cvode()` errors with a clear hint.
#'
#' @keywords internal
#' @noRd
cvodeConfig <- new.env(parent = emptyenv())

#' Package initialization
#'
#' @keywords internal
#' @importFrom reticulate py_require
#' @useDynLib cppDE, .registration = TRUE, .fixes = "C_"
#' @noRd
.onLoad <- function(libname, pkgname) {
  reticulate::py_require("sympy", python_version = ">=3.9")

  cvodeConfig$available        <- FALSE
  cvodeConfig$cflags           <- ""
  cvodeConfig$libs             <- ""
  cvodeConfig$klu_available    <- FALSE
  cvodeConfig$klu_cflags       <- ""
  cvodeConfig$klu_libs         <- ""
  cvodeConfig$cvode_lapack_available <- FALSE
  cvodeConfig$cvode_lapack_libs      <- ""
  cvodeConfig$runtime_dll_path <- ""
  cvodeConfig$openmp_available <- FALSE
  cvodeConfig$openmp_cxxflags  <- ""
  cvodeConfig$openmp_libs      <- ""

  file <- system.file("cvodeConfig.dcf", package = pkgname)
  if (nzchar(file) && file.exists(file)) {
    d <- tryCatch(read.dcf(file), error = function(e) NULL)
    if (!is.null(d) && nrow(d) >= 1L) {
      get_str <- function(k) {
        if (!(k %in% colnames(d))) return("")
        v <- d[1L, k]
        if (is.na(v)) "" else as.character(v)
      }
      cvodeConfig$available        <- identical(get_str("available"), "TRUE")
      cvodeConfig$cflags           <- get_str("cflags")
      cvodeConfig$libs             <- get_str("libs")
      cvodeConfig$klu_available    <- identical(get_str("klu_available"), "TRUE")
      cvodeConfig$klu_cflags       <- get_str("klu_cflags")
      cvodeConfig$klu_libs         <- get_str("klu_libs")
      cvodeConfig$cvode_lapack_available <- identical(get_str("cvode_lapack_available"), "TRUE")
      cvodeConfig$cvode_lapack_libs      <- get_str("cvode_lapack_libs")
      cvodeConfig$runtime_dll_path <- get_str("runtime_dll_path")
      cvodeConfig$openmp_available <- identical(get_str("openmp_available"), "TRUE")
      cvodeConfig$openmp_cxxflags  <- get_str("openmp_cxxflags")
      cvodeConfig$openmp_libs      <- get_str("openmp_libs")
    }
  }

  # The SUNDIALS / SuiteSparse DLLs from the Rtools ucrt64 sysroot live outside
  # R's default DLL search path, so the recorded bin/ goes on PATH for dyn.load().
  # A no-op off Windows, or when configure.win populated nothing.
  if (.Platform$OS.type == "windows" && nzchar(cvodeConfig$runtime_dll_path)) {
    dll_path  <- gsub("/", "\\\\", cvodeConfig$runtime_dll_path)
    cur_path  <- Sys.getenv("PATH")
    has_entry <- vapply(strsplit(cur_path, ";", fixed = TRUE)[[1]], function(p) {
      identical(tolower(gsub("/", "\\\\", p)), tolower(dll_path))
    }, logical(1))
    if (!any(has_entry)) {
      Sys.setenv(PATH = paste(dll_path, cur_path, sep = ";"))
    }
  }
}

#' Package attach
#'
#' @keywords internal
#' @noRd
.onAttach <- function(libname, pkgname) {
  .announceForkGuard()
}

# One-line report of the BLAS fork guard, or a warning when no entry point
# resolved. Silent when there is nothing to pin. Called from .onAttach, not
# .onLoad, since tooling loads namespaces without attaching them.
.announceForkGuard <- function() {
  if (isTRUE(getOption("cppDE.quiet"))) return(invisible(NULL))
  g <- tryCatch(forkGuard(), error = function(e) NULL)
  if (is.null(g) || !isTRUE(g$guard)) return(invisible(NULL))

  if (is.na(g$api)) {
    packageStartupMessage(
      "cppDE: no BLAS thread-control entry point found. If this BLAS is ",
      "threaded, a forked worker can deadlock; start R with OMP_NUM_THREADS=1. ",
      "See ?forkGuard.")
  } else if (!is.na(g$threads) && g$threads > 1L) {
    packageStartupMessage(sprintf(
      "cppDE: %s runs %d threads; each fork() pins it to 1 and restores it afterwards (?forkGuard).",
      g$api, g$threads))
  }
  invisible(NULL)
}

#' Lazy import of internal Python modules
#'
#' @keywords internal
#' @noRd
.cppde_py_cache <- new.env(parent = emptyenv())

# Directory of the Python generators; CPPDE_PY_DIR overrides the installed one.
.cppde_py_dir <- function() {
  d <- Sys.getenv("CPPDE_PY_DIR")
  if (nzchar(d)) normalizePath(d, "/", mustWork = TRUE)
  else system.file("python", package = "cppDE")
}

#' @keywords internal
#' @importFrom reticulate import_from_path
#' @noRd
get_codegen_cppODE_py <- function() {
  if (!exists("codegen_cppODE", envir = .cppde_py_cache, inherits = FALSE)) {
    .cppde_py_cache$codegen_cppODE <-
      reticulate::import_from_path(
        "codegen_cppODE",
        path = .cppde_py_dir(),
        delay_load = TRUE
      )
  }
  .cppde_py_cache$codegen_cppODE
}

#' @keywords internal
#' @importFrom reticulate import_from_path
#' @noRd
get_codegen_cppFUN_py <- function() {
  if (!exists("codegen_cppFUN", envir = .cppde_py_cache, inherits = FALSE)) {
    .cppde_py_cache$codegen_cppFUN <-
      reticulate::import_from_path(
        "codegen_cppFUN",
        path = .cppde_py_dir(),
        delay_load = TRUE
      )
  }
  .cppde_py_cache$codegen_cppFUN
}

#' @keywords internal
#' @importFrom reticulate import_from_path
#' @noRd
get_codegen_cvode_py <- function() {
  if (!exists("codegen_cvode", envir = .cppde_py_cache, inherits = FALSE)) {
    .cppde_py_cache$codegen_cvode <-
      reticulate::import_from_path(
        "codegen_cvode",
        path = .cppde_py_dir(),
        delay_load = TRUE
      )
  }
  .cppde_py_cache$codegen_cvode
}
