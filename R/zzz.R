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
#' @noRd
.onLoad <- function(libname, pkgname) {
  reticulate::py_require("sympy")

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

#' Lazy import of internal Python modules
#'
#' @keywords internal
#' @noRd
.cppde_py_cache <- new.env(parent = emptyenv())

#' @keywords internal
#' @importFrom reticulate import_from_path
#' @noRd
get_codegen_cppODE_py <- function() {
  if (!exists("codegen_cppODE", envir = .cppde_py_cache, inherits = FALSE)) {
    .cppde_py_cache$codegen_cppODE <-
      reticulate::import_from_path(
        "codegen_cppODE",
        path = system.file("python", package = "cppDE"),
        delay_load = TRUE
      )
  }
  .cppde_py_cache$codegen_cppODE
}

#' @keywords internal
#' @importFrom reticulate import_from_path
#' @noRd
get_codegen_funCpp_py <- function() {
  if (!exists("codegen_funCpp", envir = .cppde_py_cache, inherits = FALSE)) {
    .cppde_py_cache$codegen_funCpp <-
      reticulate::import_from_path(
        "codegen_funCpp",
        path = system.file("python", package = "cppDE"),
        delay_load = TRUE
      )
  }
  .cppde_py_cache$codegen_funCpp
}

#' @keywords internal
#' @importFrom reticulate import_from_path
#' @noRd
get_codegen_cvode_py <- function() {
  if (!exists("codegen_cvode", envir = .cppde_py_cache, inherits = FALSE)) {
    .cppde_py_cache$codegen_cvode <-
      reticulate::import_from_path(
        "codegen_cvode",
        path = system.file("python", package = "cppDE"),
        delay_load = TRUE
      )
  }
  .cppde_py_cache$codegen_cvode
}

#' @keywords internal
#' @importFrom reticulate import_from_path
#' @noRd
get_derivSymb_py <- function() {
  if (!exists("derivSymb", envir = .cppde_py_cache, inherits = FALSE)) {
    .cppde_py_cache$derivSymb <-
      reticulate::import_from_path(
        "derivSymb",
        path = system.file("python", package = "cppDE"),
        delay_load = TRUE
      )
  }
  .cppde_py_cache$derivSymb
}
