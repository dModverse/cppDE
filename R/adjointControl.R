#' Controls of a Reverse Solve
#'
#' @description
#' Options for the backward sweep of a model compiled with
#' `derivMode = "reverse"`: check the sweep against the tolerances and report
#' the grid it ran on. Handed to [solveODE()], [solveODEBatch()],
#' [prepareBatch()] and [solveBatch()] as `adjoint`.
#'
#' @param refine Logical, default `FALSE`. Check each step of the sweep against
#'   `reltol` and `abstol`, and with `gradtol` its share of the gradient as
#'   well; a step that fails is swept again in substeps.
#' @param gradtol Positive number, the absolute tolerance on the gradient.
#'   Default `NULL`: no test on the gradient. On a [cppODE()] model used only
#'   with `refine = TRUE`; on a [cvode()] model the absolute tolerance of the
#'   gradient quadrature in the backward problem's error test.
#' @param trace Logical, default `FALSE`. Return `$adjoint`.
#'
#' @details `refine` and `trace` are an error on a [cvode()] model.
#'
#' @return A list of class `cppDEadjointControl`.
#'
#' A solve with `trace` returns `$adjoint`, of class `cppDEadjoint`: `time`
#' and `h`, the start and length of each step the sweep ran on; `lambda`,
#' `[n_steps, n_states, n_seed]`, the adjoint state at each step's start;
#' `eta`, `[n_steps, n_seed]`, lambda times each step's local error estimate,
#' `NA` for `bdf` and `adams` under `refine`. Under `refine` also `substeps`,
#' the substeps each step was swept in, and `failures`, the steps over all seed
#' columns that failed the test at the largest substep count.
#'
#' @seealso [solveODE()]; `vignette("Methods", package = "cppDE")`, "Checking
#'   the sweep".
#' @example inst/examples/adjointControl.R
#' @export
adjointControl <- function(refine = FALSE, gradtol = NULL, trace = FALSE) {
  if (!is.logical(refine) || length(refine) != 1L || is.na(refine))
    stop("'refine' must be TRUE or FALSE", call. = FALSE)
  if (!is.logical(trace) || length(trace) != 1L || is.na(trace))
    stop("'trace' must be TRUE or FALSE", call. = FALSE)
  if (!is.null(gradtol) &&
      (!is.numeric(gradtol) || length(gradtol) != 1L || !isTRUE(gradtol > 0)))
    stop("'gradtol' must be a positive number", call. = FALSE)
  structure(list(refine = refine, gradtol = gradtol, trace = trace),
            class = "cppDEadjointControl")
}

# Whether the solve returns $adjoint.
.adjointWanted <- function(ctl)
  !is.null(ctl) && ctl$trace

# The control onto the cotangent, the channel the generated sweep reads, so the
# .Call signature stays fixed.
.adjointApply <- function(ctl, model, cotangent, is_cvode) {
  if (is.null(ctl)) return(cotangent)
  if (!inherits(ctl, "cppDEadjointControl"))
    stop("'adjoint' must come from adjointControl()", call. = FALSE)
  mode <- attr(model, "derivMode")
  if (!identical(mode, "reverse"))
    stop("'adjoint' controls the backward pass of a model compiled with ",
         "derivMode = \"reverse\"", call. = FALSE)
  if (is_cvode) {
    if (ctl$trace || ctl$refine)
      stop("on the CVODE backend 'adjoint' takes 'gradtol' alone: the backward ",
           "problem is CVODES' own continuous adjoint, with no discrete steps ",
           "to check or report", call. = FALSE)
    if (!is.null(ctl$gradtol) && !is.null(cotangent))
      attr(cotangent, "gradtol") <- ctl$gradtol
    return(cotangent)
  }
  if (!is.null(cotangent) && .adjointWanted(ctl))
    attr(cotangent, "adjointGrid") <- TRUE
  if (!is.null(cotangent) && ctl$refine) {
    attr(cotangent, "refine") <- 1
    if (!is.null(ctl$gradtol)) attr(cotangent, "gradtol") <- ctl$gradtol
  }
  cotangent
}

# The sweep's grid as $adjoint.
.adjointResult <- function(result, prep) {
  g <- result$adjointGrid
  result$adjointGrid <- NULL
  if (is.null(g) || !isTRUE(prep$adjoint_wanted)) return(result)
  if (is.null(dimnames(g$lambda)))
    dimnames(g$lambda) <- list(step = NULL, variable = prep$variables,
                               seed = prep$seed_names)
  if (is.null(dimnames(g$eta))) dimnames(g$eta) <- list(NULL, prep$seed_names)
  result$adjoint <- structure(g, class = "cppDEadjoint")
  result
}

#' @export
print.cppDEadjoint <- function(x, ...) {
  cat("<cppDEadjoint>", length(x$time), "steps,",
      dim(x$lambda)[2L], "states,", dim(x$lambda)[3L], "seed(s)\n")
  invisible(x)
}
