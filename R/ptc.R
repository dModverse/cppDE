#' Steady State by Pseudo-Transient Continuation
#'
#' @description
#' Solves \eqn{f(x, p) = 0} for variables of a [cppFUN()] model, one output per
#' solved variable, optionally under linear constraints \eqn{C x = total}.
#' `flow = TRUE` follows the ODE flow to a stable steady state; `flow = FALSE`
#' returns any regular root.
#'
#' @param fun A compiled [cppFUN()] object with forward derivatives.
#' @param x Named numeric, values of all variables of `fun`. The solved
#'   variables start there, the others keep their values.
#' @param parms Named numeric, values of the parameters of `fun`.
#' @param solve Names of the variables to solve for. Default all.
#' @param rows Names of the outputs of `fun` paired with `solve`, in the same
#'   order. Default `solve`.
#' @param C,total Constraint matrix with columns named by `solve`, and its
#'   right-hand side: `C x = total`. Default none.
#' @param flow Logical, default `TRUE`: the outputs are the right-hand side of
#'   an ODE. `FALSE` for plain equations.
#' @param positive Logical, default `TRUE`. Keep the solved variables positive;
#'   they must start positive.
#' @param controls Named list:
#'   * `reltol` (default `1e-10`) and `abstol` (`1e-14`): each residual is
#'     measured against `reltol` times its turnover
#'     \eqn{\sum_k |\partial f_i / \partial x_k|\, |x_k|}{
#'     sum_k |df_i/dx_k| |x_k|} plus `abstol` times its largest Jacobian
#'     entry. The old names `rtol` and `atol` are
#'     deprecated.
#'   * `flowTol` (`1`): relative local error of a pseudo-time step.
#'   * `maxit` (`400`): maximum number of iterations.
#'   * `dtInit` (`1e-2`): the first pseudo-time step, relative to the fastest
#'     rate if `flow = TRUE`.
#'   * `zmax` (`5`): the largest change of the logarithm of a variable per step.
#' @return A list with the named solution `x` (all variables), `converged`,
#'   `iterations` and `message`.
#' @references Kelley CT, Keyes DE (1998). Convergence analysis of
#'   pseudo-transient continuation. SIAM J Numer Anal 35(2):508-523.
#' @seealso [cppFUN()]
#' @example inst/examples/ptc.R
#' @export
ptc <- function(fun, x, parms, solve = names(x), rows = solve, C = NULL, total = NULL,
                flow = TRUE, positive = TRUE, controls = list()) {
  vars <- attr(fun, "variables"); pars <- attr(fun, "parameters")
  outs <- names(attr(fun, "equations")); mn <- attr(fun, "modelname")
  if (is.null(vars) || is.null(mn)) stop("'fun' must be a cppFUN() object")
  if (!"forward" %in% attr(fun, "derivMode")) stop("'fun' needs forward derivatives")
  if (!all(vars %in% names(x))) stop("'x' misses ", toString(setdiff(vars, names(x))))
  if (!all(pars %in% names(parms))) stop("'parms' misses ", toString(setdiff(pars, names(parms))))
  if (!all(solve %in% vars)) stop("unknown variables in 'solve': ", toString(setdiff(solve, vars)))
  if (length(rows) != length(solve) || !all(rows %in% outs))
    stop("'rows' must name one output of 'fun' per solved variable")
  if (is.null(C)) {
    C <- matrix(0, 0, length(solve)); total <- numeric(0)
  } else {
    if (!all(solve %in% colnames(C))) stop("'C' needs a column per solved variable")
    C <- C[, solve, drop = FALSE]
    if (length(total) != nrow(C)) stop("'total' needs one value per row of 'C'")
  }
  def <- list(reltol = 1e-10, abstol = 1e-14, flowTol = 1, maxit = 400L, dtInit = 1e-2, zmax = 5)
  controls <- .renamedControls(as.list(controls), c(rtol = "reltol", atol = "abstol"))
  bad <- setdiff(names(controls), names(def))
  if (length(bad)) stop("unknown controls: ", toString(bad))
  ctrl <- utils::modifyList(def, controls)
  # The native side reads the tolerances as rtol and atol.
  names(ctrl)[match(c("reltol", "abstol"), names(ctrl))] <- c("rtol", "atol")
  ctrl$flow <- as.numeric(isTRUE(flow)); ctrl$positive <- as.numeric(isTRUE(positive))
  ctrl <- lapply(ctrl, as.double)
  x0 <- as.double(x[vars])
  if (isTRUE(positive) && any(x0[match(solve, vars)] <= 0))
    stop("with positive = TRUE the solved variables must start positive")

  sym <- .nativeSym(paste0(mn, "_eval_ad"))
  if (is.null(sym)) stop("compile '", mn, "' first")
  addr <- getNativeSymbolInfo(sym$name, PACKAGE = sym$dll)$address
  res <- .Call(C_cppde_ptc, addr, length(outs), x0, as.double(parms[pars]),
               as.integer(match(solve, vars) - 1L), as.integer(match(rows, outs) - 1L),
               matrix(as.double(C), nrow(C)), as.double(total), ctrl)
  xo <- setNames(x0, vars); xo[solve] <- res$x
  list(x = xo, converged = res$converged, iterations = res$iterations, message = res$message)
}


# `controls` with deprecated names replaced by their new ones, `map` giving the
# new name for each old one. Warns once per call; both names at once is an
# error.
.renamedControls <- function(controls, map) {
  old <- intersect(names(controls), names(map))
  if (!length(old)) return(controls)
  clash <- intersect(map[old], names(controls))
  if (length(clash))
    stop("give ", toString(sQuote(clash, FALSE)), " only; ",
         toString(sQuote(old[map[old] %in% clash], FALSE)),
         " is deprecated", call. = FALSE)
  warning("deprecated control ", toString(sQuote(old, FALSE)), "; use ",
          toString(sQuote(map[old], FALSE)), call. = FALSE)
  names(controls)[match(old, names(controls))] <- map[old]
  controls
}
