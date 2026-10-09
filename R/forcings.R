#' Values of Forcings
#'
#' Evaluates forcings the way a model from [cppODE()] or [cvode()] sees them:
#' the monotone cubic Hermite interpolant (PCHIP) through the points of each forcing,
#' holding the value of the first and the last point outside them.
#'
#' @param times Numeric vector of times.
#' @param forcings Named list of forcing data, as `forcings` of [solveODE()].
#'
#' @return Numeric matrix `[length(times), length(forcings)]`, the columns
#'   named by the forcings.
#'
#' @seealso [solveODE()]
#' @example inst/examples/forcingValues.R
#' @export
forcingValues <- function(times, forcings) {
  times <- as.double(times)
  if (anyNA(times) || any(!is.finite(times))) stop("'times' must be finite")
  if (!is.list(forcings) || is.data.frame(forcings) || is.null(names(forcings)) ||
      any(!nzchar(names(forcings))))
    stop("'forcings' must be a named list")
  parsed <- Map(.parseForcing, forcings, names(forcings))
  out <- .Call(C_cppde_forcing_values, times,
               lapply(parsed, `[[`, "times"), lapply(parsed, `[[`, "values"))
  dimnames(out) <- list(NULL, names(forcings))
  out
}
