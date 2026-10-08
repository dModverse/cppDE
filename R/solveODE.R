# The points of one forcing as times and values, from a data.frame with
# columns time and value or a two-column matrix.
.parseForcing <- function(f, nm) {
  if (is.matrix(f)) f <- data.frame(time = f[,1L], value = f[,2L])
  else if (!is.data.frame(f)) f <- as.data.frame(f)
  if (!all(c("time","value") %in% names(f)))
    stop("Forcing '", nm, "' needs columns 'time' and 'value'")
  ft <- as.double(f$time); fv <- as.double(f$value)
  if (length(ft) < 1L) stop("Forcing '", nm, "' needs a time point")
  if (length(ft) != length(fv)) stop("Forcing '", nm, "': length mismatch")
  if (anyNA(ft) || any(!is.finite(ft))) stop("Forcing '", nm, "': non-finite time")
  if (anyNA(fv) || any(!is.finite(fv))) stop("Forcing '", nm, "': non-finite value")
  if (anyDuplicated(ft)) stop("Forcing '", nm, "': duplicate times")
  list(times = ft, values = fv)
}

# Marshal one condition into the 15 positional .Call arguments.  Shared by
# solveODE() and solveODEBatch() so both see identical validation.
.odeCallArgs <- function(model, times, parms,
                         tangent = NULL, hessian = NULL,
                         fixed = NULL, forcings = NULL,
                         abstol = 1e-6, reltol = 1e-6,
                         maxattempts = 50L, maxsteps = 1e6L,
                         hini = 0, roottol = 1e-6, maxroot = 1L,
                         cotangent = NULL, curvature = NULL,
                         keepStore = FALSE, store = NULL,
                         sensErrCon = TRUE, adjoint = NULL) {

  ## --- Unpack model attributes ---
  stopifnot(is.character(model), length(model) == 1L)
  required_attrs <- c("variables", "parameters", "forcings", "deriv", "deriv2", "dimNames")
  missing_attrs <- setdiff(required_attrs, names(attributes(model)))
  if (length(missing_attrs))
    stop("'model' is missing attributes: ", paste(missing_attrs, collapse = ", "))

  variables     <- attr(model, "variables")
  parameters    <- attr(model, "parameters")
  forcing_names <- attr(model, "forcings")
  deriv         <- attr(model, "deriv")
  deriv2        <- attr(model, "deriv2")
  all_sens      <- if (deriv) attr(model, "dimNames")$sens else character(0)
  backend       <- attr(model, "backend")  # "cvode" for CVODE, NULL/other for native
  is_cvode      <- identical(backend, "cvode")

  n_states  <- length(variables)
  n_params  <- length(parameters)
  n_phi_rows <- n_states + n_params

  ## --- Detect tangent shape per call ---
  ## State-only [n_states, n_active] leaves the parameter rows at the identity, full
  ## Phi' supplies them, partial names the rows it supplies and zero-pads the rest.
  if (!is.null(tangent) && !deriv)
    stop("'tangent' supplied but model has deriv = FALSE")
  if (!is.null(hessian) && !deriv2)
    stop("'hessian' supplied but model has deriv2 = FALSE")

  ## --- Reverse mode: which direction the model was built for ---
  ## The mode is stamped on the model, not passed per call, exactly as `deriv`
  ## is: it decides which code was emitted and cannot be chosen afterwards.
  is_reverse <- attr(model, "derivMode") %in% c("reverse", "forward-reverse")
  if (!is.null(cotangent) && !is_reverse)
    stop("'cotangent' supplied but the model was not compiled with derivMode = \"reverse\"")
  ## A reverse call without a cotangent is the value half of the pair: it
  ## integrates, fills the store and sweeps nothing.
  if (is.null(cotangent) && is_reverse && !isTRUE(keepStore))
    stop("a model compiled with derivMode = \"reverse\" needs a 'cotangent', or ",
         "keepStore = TRUE to run it for its values alone")
  if (!is.null(cotangent)) {
    if (!is.numeric(cotangent)) stop("'cotangent' must be numeric")
    d <- dim(cotangent)
    if (is.null(d) || length(d) < 2L || length(d) > 3L)
      stop("'cotangent' must be a [n_out, n_states] matrix or an ",
           "[n_out, n_states, n_seed] array")
    if (d[2] != n_states)
      stop("'cotangent' has ", d[2], " state columns, the model has ", n_states)
    storage.mode(cotangent) <- "double"
    ## The cotangent's derivative along the tangent, [n_out, n_states, n_seed,
    ## n_sens]. It rides on the cotangent so the .Call signature stays fixed.
    if (!is.null(curvature)) {
      if (!identical(attr(model, "derivMode"), "forward-reverse"))
        stop("a curvature needs derivMode = \"forward-reverse\"; the first ",
             "order has no slot for it", call. = FALSE)
      if (!is.numeric(curvature)) stop("'curvature' must be numeric", call. = FALSE)
      sd <- dim(curvature)
      if (length(d) == 2L) d <- c(d, 1L)
      if (length(sd) != 4L || !identical(as.integer(sd[1:3]), as.integer(d[1:3])))
        stop("'curvature' must be [n_out, n_states, n_seed, n_sens] on the ",
             "cotangent's own first three dimensions", call. = FALSE)
      storage.mode(curvature) <- "double"
      attr(cotangent, "curvature") <- curvature
    }

  }
  if (!is.null(curvature) && is.null(cotangent))
    stop("'curvature' is the derivative of a cotangent and needs one",
         call. = FALSE)


  is_2d_tangent <- !is.null(tangent) &&
    (is.matrix(tangent) ||
       (is.array(tangent) && length(dim(tangent)) == 2L))
  tangent_is_legacy <-
    is_2d_tangent && n_phi_rows != n_states && nrow(tangent) == n_states && {
      rn <- rownames(tangent)
      is.null(rn) || all(rn %in% variables)
    }
  ## "is_full" means a full or partial Phi', anything but the state-only shape.
  ## It gates runtime `fixed` and sets n_theta_active and the hessian shape.
  tangent_is_full <- is_2d_tangent && !tangent_is_legacy

  ## --- Runtime fixed: incompatible with full-shape tangent ---
  fixed_indices <- integer(0)
  if (!is.null(fixed)) {
    if (!deriv) { warning("'fixed' ignored when deriv = FALSE") }
    else if (tangent_is_full) {
      stop("'fixed' is not supported together with full-shape tangent; ",
           "express fixedness through zero rows in tangent instead")
    } else {
      if (!is.character(fixed)) stop("'fixed' must be a character vector")
      bad <- setdiff(fixed, all_sens)
      if (length(bad)) stop("Unknown 'fixed' names: ", paste(bad, collapse = ", "),
                            "\nValid: ", paste(all_sens, collapse = ", "))
      fixed_indices <- match(fixed, all_sens) - 1L  # 0-based for C++
    }
  }
  active_sens <- if (deriv && length(fixed_indices))
    all_sens[-( fixed_indices + 1L)] else all_sens
  n_active  <- length(active_sens)

  ## --- Per-call active sens dimension ---
  ## tangent full shape: M = ncol(tangent) (theta count, may differ from n_active).
  ## otherwise: M = n_active (state-only or identity seeding uses the active basis).
  n_theta_active <- if (tangent_is_full) as.integer(ncol(tangent)) else n_active

  ## --- Identity-on-active-params padding for the state-only shape ---
  build_param_identity <- function(col_names) {
    pad <- matrix(0, nrow = n_params, ncol = length(col_names),
                  dimnames = list(parameters, col_names))
    for (j in seq_len(n_params)) {
      col_idx <- match(parameters[j], col_names)
      if (!is.na(col_idx)) pad[j, col_idx] <- 1.0
    }
    pad
  }

  ## --- Output sens column names (per call) ---
  ## Full Phi' uses its own colnames, or theta1..thetaM when it has none;
  ## state-only and NULL use active_sens, the model-parameter basis.
  sens_col_names <- if (tangent_is_full) {
    cn <- colnames(tangent)
    if (!is.null(cn)) cn else sprintf("theta%d", seq_len(n_theta_active))
  } else {
    active_sens
  }

  ## --- Coerce tangent to flat [n_phi_rows, n_theta_active] ---
  ## The reorder is skipped when the order already matches, the common case when
  ## an optimiser reuses one shape. State-only padding fills a preallocated matrix.
  coerce_tangent <- function(x, n_cols, col_names, arg) {
    if (!is.numeric(x)) stop("'", arg, "' must be numeric")

    is_2d <- is.matrix(x) || (is.array(x) && length(dim(x)) == 2L)
    if (!is_2d) {
      ## Vector form: state-only if length == n_states*n_cols, else full.
      ## (Partial-row form requires names and is therefore matrix-only.)
      if (length(x) == n_states * n_cols && n_phi_rows != n_states) {
        out <- matrix(0, n_phi_rows, n_cols)
        out[seq_len(n_states), ] <- x
        out[(n_states + 1L):n_phi_rows, ] <- build_param_identity(col_names)
        return(as.double(out))
      }
      if (length(x) != n_phi_rows * n_cols)
        stop(sprintf("'%s' must have length %d (n_phi_rows * n_cols) or %d (legacy)",
                     arg, n_phi_rows * n_cols, n_states * n_cols))
      return(as.double(x))
    }

    nr <- nrow(x); nc <- ncol(x)
    rn <- rownames(x); cn <- colnames(x)
    expected_rows <- c(variables, parameters)

    ## Column reorder/check (once, regardless of row interpretation).
    if (!is.null(cn)) {
      if (!setequal(cn, col_names))
        stop("'", arg, "' column names must match: ",
             paste(col_names, collapse = ", "))
      if (!identical(cn, col_names))
        x <- x[, col_names, drop = FALSE]
    } else if (nc != n_cols) {
      stop(sprintf("'%s' must have %d columns", arg, n_cols))
    }

    ## Full Phi' shape: [n_phi_rows, n_cols].
    if (nr == n_phi_rows) {
      if (!is.null(rn)) {
        if (!setequal(rn, expected_rows))
          stop("'", arg, "' row names must be c(variables, parameters)")
        if (!identical(rn, expected_rows))
          x <- x[expected_rows, , drop = FALSE]
      }
      return(as.double(x))
    }

    ## State-only shape: [n_states, n_cols], no rownames or all-variable rownames.
    if (nr == n_states && n_phi_rows != n_states &&
        (is.null(rn) || all(rn %in% variables))) {
      if (!is.null(rn)) {
        if (!setequal(rn, variables))
          stop("'", arg, "' row names must match variables (legacy shape)")
        if (!identical(rn, variables))
          x <- x[variables, , drop = FALSE]
      }
      out <- matrix(0, n_phi_rows, n_cols)
      out[seq_len(n_states), ] <- x
      out[(n_states + 1L):n_phi_rows, ] <- build_param_identity(col_names)
      return(as.double(out))
    }

    ## Partial-row shape: rownames required, subset of c(variables, parameters).
    ## Missing rows are padded with zeros (= implicit fixed for those slots).
    if (is.null(rn))
      stop(sprintf(
        "'%s' has shape [%d, %d]; expected [%d, %d] (full Phi'), [%d, %d] (legacy), or a partial-row matrix with rownames identifying a subset of c(variables, parameters)",
        arg, nr, nc, n_phi_rows, n_cols, n_states, n_cols))
    bad <- setdiff(rn, expected_rows)
    if (length(bad))
      stop("'", arg, "' has unknown row names: ", paste(bad, collapse = ", "))
    if (anyDuplicated(rn))
      stop("'", arg, "' has duplicate row names")
    out <- matrix(0, n_phi_rows, n_cols)
    out[match(rn, expected_rows), ] <- x
    as.double(out)
  }

  ## --- Build tangent for the C++ side ---
  ## CVODE always needs a full Phi', the generated cppDE code accepts NULL and
  ## identity-seeds via diff(ai). Runtime `fixed` becomes zero rows in a default.
  if (deriv && is.null(tangent) && is_cvode) {
    default_pp <- matrix(0, nrow = n_phi_rows, ncol = n_active,
                         dimnames = list(c(variables, parameters), active_sens))
    for (j in seq_along(active_sens)) {
      r <- match(active_sens[j], c(variables, parameters))
      if (!is.na(r)) default_pp[r, j] <- 1.0
    }
    tangent <- as.double(default_pp)
    dim(tangent) <- c(n_phi_rows, n_active)
  } else if (!is.null(tangent)) {
    flat <- coerce_tangent(tangent, n_theta_active, sens_col_names, "tangent")
    ## Preserve 2-D shape for C++ Rf_ncols(), distinguishes [phi_rows, M] from a flat vector.
    dim(flat) <- c(n_phi_rows, n_theta_active)
    tangent <- flat
  }

  if (!is.null(adjoint) && !inherits(adjoint, "cppDEadjointControl"))
    stop("'adjoint' must come from adjointControl()", call. = FALSE)

  if (!is.null(hessian)) {
    if (!is.numeric(hessian)) stop("'hessian' must be numeric")
    ## State-only [n_states, M, M] is accepted only alongside a state-only or absent
    ## tangent, where Phi'' vanishes on the parameter block. Full is the unified
    ## shape; partial names the rows it supplies and the C++ side gets it padded.
    legacy_d2_ok <- !tangent_is_full && n_phi_rows != n_states
    expected_rows <- c(variables, parameters)

    if (is.array(hessian) && length(dim(hessian)) == 3) {
      d <- dim(hessian); nr <- d[1L]; nc1 <- d[2L]; nc2 <- d[3L]
      dn <- dimnames(hessian)
      rn  <- if (length(dn) >= 1L) dn[[1L]] else NULL
      cn1 <- if (length(dn) >= 2L) dn[[2L]] else NULL
      cn2 <- if (length(dn) >= 3L) dn[[3L]] else NULL

      if (nc1 != n_theta_active || nc2 != n_theta_active)
        stop(sprintf("'hessian' must have dim 2 and 3 equal to %d", n_theta_active))

      ## Column reorder/check on dims 2 and 3 (once each).
      if (!is.null(cn1)) {
        if (!setequal(cn1, sens_col_names))
          stop("'hessian' dim 2 must match sens columns")
        if (!identical(cn1, sens_col_names))
          hessian <- hessian[, sens_col_names, , drop = FALSE]
      }
      if (!is.null(cn2)) {
        if (!setequal(cn2, sens_col_names))
          stop("'hessian' dim 3 must match sens columns")
        if (!identical(cn2, sens_col_names))
          hessian <- hessian[, , sens_col_names, drop = FALSE]
      }

      ## Full shape: [n_phi_rows, M, M].
      if (nr == n_phi_rows) {
        if (!is.null(rn)) {
          if (!setequal(rn, expected_rows))
            stop("'hessian' dim 1 must be c(variables, parameters)")
          if (!identical(rn, expected_rows))
            hessian <- hessian[expected_rows, , , drop = FALSE]
        }
        hessian <- as.double(hessian)

      ## State-only shape: [n_states, M, M], state rownames or none.
      } else if (nr == n_states && legacy_d2_ok &&
                 (is.null(rn) || all(rn %in% variables))) {
        if (!is.null(rn)) {
          if (!setequal(rn, variables))
            stop("'hessian' dim 1 must match variables (legacy shape)")
          if (!identical(rn, variables))
            hessian <- hessian[variables, , , drop = FALSE]
        }
        full <- array(0, dim = c(n_phi_rows, n_theta_active, n_theta_active))
        full[seq_len(n_states), , ] <- hessian
        hessian <- as.double(full)

      ## Partial-row shape: dim-1 names required, subset of expected_rows.
      } else {
        if (is.null(rn))
          stop(sprintf(
            "'hessian' has shape [%d, %d, %d]; expected [%d, %d, %d] (full Phi''), [%d, %d, %d] (legacy), or a partial-row array with dim-1 names identifying a subset of c(variables, parameters)",
            nr, nc1, nc2, n_phi_rows, n_theta_active, n_theta_active,
            n_states, n_theta_active, n_theta_active))
        bad <- setdiff(rn, expected_rows)
        if (length(bad))
          stop("'hessian' has unknown dim-1 names: ", paste(bad, collapse = ", "))
        if (anyDuplicated(rn))
          stop("'hessian' has duplicate dim-1 names")
        full <- array(0, dim = c(n_phi_rows, n_theta_active, n_theta_active))
        full[match(rn, expected_rows), , ] <- hessian
        hessian <- as.double(full)
      }

    } else {
      ## Vector form (no partial; needs names).
      len <- length(hessian)
      if (len == n_phi_rows * n_theta_active^2) {
        hessian <- as.double(hessian)
      } else if (len == n_states * n_theta_active^2 && legacy_d2_ok) {
        pad <- array(0, dim = c(n_phi_rows, n_theta_active, n_theta_active))
        pad[seq_len(n_states), , ] <- array(as.double(hessian),
                                            dim = c(n_states, n_theta_active, n_theta_active))
        hessian <- as.double(pad)
      } else {
        stop(sprintf("'hessian' must have length %d", n_phi_rows * n_theta_active^2))
      }
    }
  }

  ## --- times ---
  if (!is.numeric(times) || !length(times) || anyNA(times) || any(!is.finite(times)))
    stop("'times' must be a non-empty finite numeric vector")
  # An upper bound on the step size, list(time, h), rides on `times`; only the
  # tests set it.
  hmax <- attr(times, "hmax")
  times <- as.double(times)
  attr(times, "hmax") <- hmax

  ## Store flags travel as attributes of `times`, since the call that makes a
  ## store has no cotangent. Native backend only: CVODES keeps its checkpoints.
  if (isTRUE(keepStore) || !is.null(store)) {
    what <- if (isTRUE(keepStore)) "keepStore" else "store"
    if (!is_reverse)
      stop("'", what, "' belongs to a model compiled with derivMode = \"reverse\"",
           call. = FALSE)
    if (is_cvode)
      stop("'", what, "' is not available on the CVODE backend: CVODES holds ",
           "its checkpoints itself. Use cppODE() for a pair of solves that ",
           "share one integration.", call. = FALSE)
    ## A checkpoint's tangents point into the arena of the solve that took it,
    ## which is reset on return, so forward-reverse cannot reuse a store.
    if (identical(attr(model, "derivMode"), "forward-reverse"))
      stop("'", what, "' is not available under derivMode = ",
           "\"forward-reverse\": a checkpoint's tangents live in the arena of ",
           "the solve that took them and do not outlive it. Let the second ",
           "solve integrate.", call. = FALSE)
  }
  ## sensErrCon: whether error control reads the tangents (see ?solveODE).
  if (!is.logical(sensErrCon) || length(sensErrCon) != 1L || is.na(sensErrCon))
    stop("'sensErrCon' must be TRUE or FALSE", call. = FALSE)
  if (!sensErrCon) {
    if (!isTRUE(attr(model, "deriv")) &&
        !identical(attr(model, "derivMode"), "reverse"))
      stop("'sensErrCon' weighs the sensitivity error against the state error, ",
           "and this model computes no sensitivities", call. = FALSE)
    if (is_cvode)
      stop("'sensErrCon' is not available on the CVODE backend", call. = FALSE)
    attr(times, "sensErrCon") <- FALSE
  }
  if (isTRUE(keepStore)) attr(times, "keepStore") <- TRUE
  if (!is.null(store)) {
    if (!inherits(store, "externalptr"))
      stop("'store' must be the `store` element of an earlier solve",
           call. = FALSE)
    attr(times, "store") <- store
  }

  ## --- parms ---
  if (!is.numeric(parms) || is.null(names(parms)))
    stop("'parms' must be a named numeric vector")
  required_nms <- c(variables, parameters)
  miss <- setdiff(required_nms, names(parms))
  if (length(miss)) stop("'parms' missing: ", paste(miss, collapse = ", "))
  parms_ordered <- as.double(parms[required_nms])
  if (anyNA(parms_ordered) || any(!is.finite(parms_ordered)))
    stop("'parms' must be finite")

  ## --- forcings ---
  n_forcings <- length(forcing_names)
  if (n_forcings && is.null(forcings))
    stop("Model requires forcings: ", paste(forcing_names, collapse = ", "))

  if (!n_forcings) {
    forcing_times_list <- forcing_values_list <- list()
  } else {
    if (!is.list(forcings) || is.null(names(forcings)))
      stop("'forcings' must be a named list")
    miss_f <- setdiff(forcing_names, names(forcings))
    if (length(miss_f)) stop("Missing forcings: ", paste(miss_f, collapse = ", "))
    parsed <- lapply(forcing_names, function(nm) .parseForcing(forcings[[nm]], nm))
    forcing_times_list  <- lapply(parsed, `[[`, "times")
    forcing_values_list <- lapply(parsed, `[[`, "values")
  }

  ## --- solver options ---
  if (!is.numeric(abstol)  || abstol  <= 0) stop("'abstol' must be positive")
  if (!is.numeric(reltol)  || reltol  <= 0) stop("'reltol' must be positive")
  if (!is.numeric(hini)    || hini    <  0) stop("'hini' must be non-negative")
  if (!is.numeric(roottol) || roottol <= 0) stop("'roottol' must be positive")
  maxattempts <- as.integer(maxattempts); maxsteps <- as.integer(maxsteps); maxroot <- as.integer(maxroot)
  if (maxattempts <= 0L) stop("'maxattempts' must be positive")
  if (maxsteps    <= 0L) stop("'maxsteps' must be positive")
  if (maxroot     <= 0L) stop("'maxroot' must be positive")

  cotangent <- .adjointApply(adjoint, model, cotangent, is_cvode)

  list(call_args = list(times, parms_ordered, tangent, hessian, fixed_indices,
                        as.double(abstol), as.double(reltol), maxattempts, maxsteps,
                        as.double(hini), as.double(roottol), maxroot,
                        forcing_times_list, forcing_values_list, cotangent),
       times = times, variables = variables, sens_col_names = sens_col_names,
       adjoint_wanted = .adjointWanted(adjoint),
       theta_names = c(variables, parameters),
       seed_names = if (!is.null(cotangent) && length(dim(cotangent)) == 3L)
                      dimnames(cotangent)[[3]] else NULL)
}


# Attach dimnames, enrich diagnostics, apply onFailure, shape the trace.
.odeFinish <- function(result, model, prep, onFailure, traceFile = NULL) {
  variables <- prep$variables
  out_sens  <- prep$sens_col_names

  ## Only when the generated code did not already set them.  It does so on
  ## the batch pre-allocation path, while the arrays are still unaliased;
  ## repeating it here would see a refcount above one and duplicate each array.
  if (!is.null(result$variable) && is.null(dimnames(result$variable)))
    colnames(result$variable) <- variables
  if (!is.null(result$tangent) && is.null(dimnames(result$tangent)))
    dimnames(result$tangent) <- list(time = NULL, variable = variables, sens = out_sens)
  if (!is.null(result$hessian) && is.null(dimnames(result$hessian)))
    dimnames(result$hessian) <- list(time = NULL, variable = variables,
                                     sens1 = out_sens, sens2 = out_sens)
  ## The cotangent of the inputs, one row per model parameter, states first,
  ## indexed as the argument `tangent`.
  if (!is.null(result$cotangent) && is.null(dimnames(result$cotangent)))
    dimnames(result$cotangent) <- list(prep$theta_names, prep$seed_names)
  ## Forward-reverse: the cotangent's derivatives, one block per tangent
  ## direction. Under the identity tangent that block is a Hessian.
  if (!is.null(result$curvature) && is.null(dimnames(result$curvature)))
    dimnames(result$curvature) <- list(theta = prep$theta_names, sens = out_sens,
                                       seed = prep$seed_names)
  ## The sweep's own grid, as $adjoint.
  result <- .adjointResult(result, prep)

  diag <- result$diagnostics
  if (!is.null(diag)) {
    diag$method  <- attr(model, "method")
    diag$useNDF  <- attr(model, "useNDF")
    diag$backend <- attr(model, "backend")
    result$diagnostics <- diag
  }
  if (!is.null(diag) && diag$return_code != 0L) {
    msg <- paste0(
      "Solver did not complete: ", diag$message,
      "\n  Reached t = ", format(diag$t_reached, digits = 6),
      " (", length(result$time), " of ", length(prep$times), " time points)."
    )
    switch(onFailure,
           stop   = stop(msg, call. = FALSE),
           warn   = warning(paste0(msg, "\n  Returning partial results."),
                            call. = FALSE, immediate. = TRUE),
           silent = invisible(NULL))
  }

  if (!is.null(result$trace) && length(result$trace$nst) > 0L) {
    result$trace <- as.data.frame(result$trace, stringsAsFactors = FALSE)
    if (!is.null(traceFile)) {
      stopifnot(is.character(traceFile), length(traceFile) == 1L, nzchar(traceFile))
      utils::write.csv(result$trace, traceFile, row.names = FALSE)
    }
  } else {
    result$trace <- NULL
  }

  result
}


#' Run a Compiled ODE Model
#'
#' @description
#' Integrates a model compiled by [cppODE()] or [cvode()] over the time points
#' `times`. Returns the state trajectory and, when the model was compiled with
#' derivatives, their tangent and Hessian (forward) or the cotangent and
#' curvature of a seeded functional (reverse).
#'
#' @details
#' ## Dimensions
#'
#' `n_t` is the number of output rows, `nrow(result$variable)`: `length(times)`,
#' plus one row for each event inside the window and for the time 0 that
#' `includeTimeZero` adds. `n_x` is the number of states, `n_p` the number of
#' parameters, `n_s` the number of tangent directions and `n_seed` the number of
#' cotangent columns.
#'
#' ## Derivative arguments and results
#'
#' A derivative argument and its result have the same name. `tangent` and
#' `hessian` are the first and second derivative in \eqn{\theta} of the
#' inputs going in, and of the states coming out. `cotangent` is the gradient
#' of a functional with respect to the outputs going in, and with respect to
#' the inputs coming out. `curvature` is the derivative of the cotangent
#' along the tangent: that functional's Hessian applied to the tangent, on the
#' outputs going in and on the inputs coming out.
#'
#' ## Shapes of `tangent` and `hessian`
#'
#' `tangent` is the Jacobian \eqn{\Phi'(\theta)} and `hessian` the second
#' derivative \eqn{\Phi''(\theta)} of a reparametrisation
#' \eqn{p = \Phi(\theta)} of the initial states and parameters; the results are
#' then derivatives in \eqn{\theta}. `tangent` takes one of three shapes, told
#' apart by its row count and row names. `hessian` has the rows of `tangent` and
#' two dimensions of length `n_s`.
#'
#' - State-only, `[n_x, n_s]`, rows unnamed or named by the states. The
#'   parameter rows are the identity on the active sensitivities, the model's
#'   sensitivity names minus `fixed`. `result$tangent[t, , ]` has this shape and
#'   can seed a following solve.
#' - Full, `[n_x + n_p, n_s]`, states first: \eqn{\Phi'(\theta)} itself. `n_s`
#'   may change from call to call.
#' - Partial, `[k, n_s]` with `k < n_x + n_p`. The row names, a subset of
#'   `c(variables, parameters)`, place the rows; missing rows are zero, so
#'   those inputs are fixed.
#'
#' Run-time `fixed` takes the state-only shape only. Column names, when given,
#' are the active sensitivity names (state-only) or name the directions
#' \eqn{\theta} (full and partial); unnamed full or partial columns are named
#' `theta1`, `theta2`, and so on.
#'
#' @param model A compiled ODE model returned by [cppODE()] or [cvode()].
#' @param times Numeric vector of output times, non-empty and finite.
#' @param parms Named numeric vector of initial conditions and parameters.
#'   Names must include all of
#'   `c(attr(model, "variables"), attr(model, "parameters"))`.
#' @param tangent Optional numeric matrix, the tangent of the inputs; shapes in
#'   Details. Default `NULL`: the identity on the active sensitivities.
#' @param hessian Optional numeric array, the second derivative of the inputs;
#'   shapes in Details. Needs a model with `attr(model, "deriv2")`. Default
#'   `NULL`: the inputs are affine in \eqn{\theta}.
#' @param fixed Optional character vector, sensitivity names held fixed in this
#'   call, a subset of `attr(model, "dimNames")$sens`. Unlike `fixed` of
#'   [cppODE()] it changes without recompiling. Not with a full or partial
#'   `tangent`. Default `NULL`.
#' @param forcings Named list of forcing data, one element per name in
#'   `attr(model, "forcings")`: a `data.frame` with columns `time` and `value`,
#'   or a two-column matrix. [cppODE()] and [cvode()] take the forcing names,
#'   this argument the data. A forcing is the monotone cubic Hermite
#'   interpolant (PCHIP) of its points and constant beyond the first and the
#'   last point, see [forcingValues()]. Default `NULL`.
#' @param abstol Absolute error tolerance. Default `1e-6`.
#' @param reltol Relative error tolerance. Default `1e-6`.
#' @param maxattempts Maximum number of consecutive step attempts without
#'   progress in time before the solve stops with `return_code = -4`. Default
#'   `50`.
#' @param maxsteps Maximum total number of integration steps. Default
#'   `1e6`.
#' @param hini Initial step size; `0` (default) triggers automatic
#'   estimation.
#' @param roottol Tolerance for locating root events. Default `1e-6`. Under
#'   `rootfunc = "equilibrate"` the solve stops once every
#'   \eqn{|\dot x_i| \le \mathrm{roottol}\, |x_i| + \mathrm{abstol}}{
#'   |dx_i/dt| <= roottol |x_i| + abstol}; on a [cppODE()] model the test
#'   covers the sensitivities as well.
#' @param maxroot Maximum number of triggers per root event. Default `1`.
#' @param onFailure How to react when the solver returns a non-zero
#'   return code. One of `"stop"` (default; raise an error with the solver
#'   message and no partial results), `"warn"` (emit a warning and return
#'   partial results up to `t_reached`), or `"silent"` (return the partial
#'   result without any signal).
#' @param traceFile Optional path of a CSV file the step trace is written to.
#'   Needs a model compiled with `stepTrace = TRUE`, ignored otherwise. Default
#'   `NULL`.
#' @param cotangent The cotangent of the outputs, `[n_t, n_x]` or
#'   `[n_t, n_x, n_seed]`. Required by a model compiled with
#'   `derivMode = "reverse"` or `"forward-reverse"` unless `keepStore = TRUE`;
#'   an error on any other model. `result$cotangent` is then
#'   \eqn{\sum_{t,i} w_{t,i}\, \partial x_i(t) / \partial \theta}{
#'   sum_{t,i} w_{t,i} dx_i(t)/dtheta}, one column per seed. Default `NULL`.
#' @param curvature Optional, `"forward-reverse"` only: the derivative of
#'   `cotangent` along the tangent, `[n_t, n_x, n_seed, n_s]`. For a cotangent
#'   that is the gradient of a functional of the outputs, this is that
#'   functional's Hessian applied to the output tangent. Default `NULL`: the
#'   cotangent is constant in \eqn{\theta}.
#' @param keepStore Logical, default `FALSE`. Return the checkpoints of a
#'   reverse solve as `$store`, for a later solve to reuse through `store`. The
#'   `cotangent` may then be omitted, and the solve returns the values alone.
#'   Only on a [cppODE()] model compiled with `derivMode = "reverse"`.
#' @param store The `$store` of an earlier solve of the same model at the same
#'   `times` and `parms`. The solve then integrates nothing and runs the
#'   backward sweep alone. A store from a different point is an error. A store
#'   may be reused any number of times. Default `NULL`.
#' @param sensErrCon Logical. Include the sensitivities in the error test.
#'   Default `TRUE`; `FALSE` is cheaper and less accurate. Needs a model with
#'   sensitivities, has no effect under `derivMode = "reverse"`, and a
#'   [cvode()] model takes `TRUE` only.
#' @param adjoint Optional [adjointControl()] object for a model compiled with
#'   `derivMode = "reverse"`. Default `NULL`, the same as `adjointControl()`.
#' @param maxattemps Deprecated spelling of `maxattempts`.
#'
#' @return A named list:
#' - `time`: the output times, length `n_t`.
#' - `variable`: the states, `[n_t, n_x]`.
#' - `tangent`: when `attr(model, "deriv")` is `TRUE`, `[n_t, n_x, n_s]`. The
#'   third dimension is named by the active sensitivities or the columns of
#'   `tangent`.
#' - `hessian`: when `attr(model, "deriv2")` is `TRUE`, `[n_t, n_x, n_s, n_s]`.
#' - `cotangent`: under `derivMode = "reverse"` or `"forward-reverse"`, the
#'   cotangent of the inputs, `[n_x + n_p, n_seed]`, rows as in a full
#'   `tangent`.
#' - `curvature`: under `"forward-reverse"`, `[n_x + n_p, n_s, n_seed]`, the
#'   derivative of each `cotangent` entry along the tangent. Under the identity
#'   tangent its columns are the Hessian of the seeded functional. The order of
#'   the last two dimensions is that of a Hessian, not that of the argument
#'   `curvature`.
#' - `diagnostics`: solver statistics, printed by [diagnostics()]:
#'   `return_code` (0 on success, negative codes as in SUNDIALS CVODE),
#'   `message`, `accepted` and `rejected` (step counts), `fevals`, `jevals`
#'   and `setups` (right-hand side evaluations, Jacobian evaluations, matrix
#'   factorisations), `last_dt` and `last_order` (size and order of the last
#'   successful step), `t_reached`, and the model's `method`, `useNDF` and,
#'   on a [cvode()] model, `backend`.
#' - `trace`: with `stepTrace = TRUE`, a `data.frame` with one row per step
#'   attempt: `nst` (accepted steps before the attempt), `t` (time at the end
#'   of the step), `h` (signed step size), `q` (order), `dsm` (error test
#'   value, the step passes below 1), `acnrm` and `acnrm_state` (weighted RMS
#'   norm of the local error over all components and over the states alone),
#'   `tq2` (error constant), `gamma` and `gamrat` (step size times the leading
#'   coefficient, and its ratio to the value at the last factorisation),
#'   `newton_conv` (1 when the corrector converged), `mode` (`"BDF"`, `"NDF"`,
#'   `"ADAMS"`, `"ONESTEP"`, or `"CVODE"` and `"CVODE_event"` on a [cvode()]
#'   model), `nfe`, `njev` and `nsetups` (cumulative counters) and
#'   `setup_reason` (what triggered the last factorisation). A [cvode()] model
#'   reports `NaN` where CVODES exposes no value.
#' - `adjoint`: under `adjointControl(trace = TRUE)`, see [adjointControl()].
#' - `store`: with `keepStore = TRUE`, an external pointer to the checkpoints.
#'
#' @seealso [cppODE()] and [cvode()] for model compilation;
#'   [solveODEBatch()] for many conditions in one call; [adjointControl()];
#'   [forcingValues()]; [diagnostics()] for printing solver statistics.
#'
#' @example inst/examples/solveODE.R
#' @export
solveODE <- function(model, times, parms,
                     tangent = NULL, hessian = NULL,
                     fixed = NULL, forcings = NULL,
                     abstol = 1e-6, reltol = 1e-6,
                     maxattempts = 50L, maxsteps = 1e6L,
                     hini = 0, roottol = 1e-6, maxroot = 1L,
                     onFailure = c("stop", "warn", "silent"),
                     traceFile = NULL, cotangent = NULL, curvature = NULL,
                     keepStore = FALSE, store = NULL,
                     sensErrCon = TRUE, adjoint = NULL, maxattemps = NULL) {

  onFailure <- match.arg(onFailure)
  maxattempts <- .renamedArg(maxattempts, maxattemps, !missing(maxattempts),
                             "maxattempts", "maxattemps")

  prep <- .odeCallArgs(model, times, parms, tangent, hessian, fixed, forcings,
                       abstol, reltol, maxattempts, maxsteps, hini, roottol, maxroot,
                       cotangent, curvature, keepStore, store, sensErrCon,
                       adjoint)

  SYM <- .nativeSym(paste0("solve_", as.character(model)))
  if (is.null(SYM)) stop("Model not loaded. Run compile() first.", call. = FALSE)

  ## The dimnames are passed in so the generated code can attach them while the
  ## arrays are unaliased; setting them here would duplicate every array.
  dn <- list(prep$variables, prep$sens_col_names)
  result <- tryCatch(
    do.call(.callSym, c(list(SYM), prep$call_args, list(dn))),
    error = function(e) stop("ODE solver error: ", e$message, call. = FALSE))

  .odeFinish(result, model, prep, onFailure, traceFile)
}


#' Solve Many Conditions in One Call
#'
#' @description
#' Integrates a compiled ODE model over several independent parameter sets
#' (experimental conditions, or the subjects of a mixed-effects fit) inside a
#' single `.Call`, using OpenMP where the toolchain provides it.
#'
#' @details
#' Results are bit-identical to the serial path. The batch falls back to a
#' serial loop in a forked child, inside an OpenMP region, or without OpenMP;
#' see [batchAvailable()]. A model without a batch entry point is solved by
#' [solveODE()] per condition.
#'
#' @param model A model handle from [cppODE()] or [cvode()].
#' @param conditions A named or unnamed list of per-condition argument lists.
#'   Recognised names are `times`, `parms`, `tangent`, `hessian`, `cotangent`,
#'   `curvature`, `fixed`, `forcings`, `abstol`, `reltol`, `maxattempts`,
#'   `maxsteps`, `hini`, `roottol`, `maxroot`, `keepStore`, `store`,
#'   `sensErrCon` and `adjoint`; each overrides the batch-wide argument of the
#'   same name.
#' @param times,parms Batch-wide `times` and `parms` as in [solveODE()], for
#'   the conditions that do not give their own. Default `NULL`.
#' @param traceFile Optional. Either one path per condition, or a single path
#'   used as a template, in which case the condition's name (or its index) is
#'   inserted before the extension. Needs a model built with `stepTrace = TRUE`.
#'   Default `NULL`.
#' @param cores Number of threads, capped by `length(conditions)`. `NULL`
#'   (default) takes the first of `getOption("cppDE.cores")`,
#'   `getOption("Ncpus")` and `parallel::detectCores(logical = FALSE)` that is
#'   set. `1` forces the serial loop.
#' @param onFailure Applied once after all conditions have run. `"stop"`
#'   (default) runs every condition, then raises one error naming the failed
#'   ones; `"warn"` and `"silent"` as in [solveODE()].
#' @inheritParams solveODE
#'
#' @return A list of [solveODE()] results, one per condition, named after
#'   `conditions`, with attribute `threads`, the thread count used.
#'
#' @seealso [solveODE()]; [prepareBatch()] and [solveBatch()] for repeated
#'   solves of the same conditions; [batchAvailable()].
#' @example inst/examples/solveODEBatch.R
#' @export
solveODEBatch <- function(model, conditions,
                          times = NULL, parms = NULL,
                          tangent = NULL, hessian = NULL,
                          fixed = NULL, forcings = NULL,
                          abstol = 1e-6, reltol = 1e-6,
                          maxattempts = 50L, maxsteps = 1e6L,
                          hini = 0, roottol = 1e-6, maxroot = 1L,
                          cores = NULL,
                          traceFile = NULL,
                          onFailure = c("stop", "warn", "silent"),
                          cotangent = NULL, curvature = NULL,
                          keepStore = FALSE, store = NULL,
                          sensErrCon = TRUE, adjoint = NULL,
                          maxattemps = NULL) {

  onFailure <- match.arg(onFailure)
  maxattempts <- .renamedArg(maxattempts, maxattemps, !missing(maxattempts),
                             "maxattempts", "maxattemps")
  preps <- .batchPreps(model, conditions, times, parms, tangent, hessian,
                       fixed, forcings, abstol, reltol, maxattempts, maxsteps,
                       hini, roottol, maxroot, cotangent, curvature,
                       keepStore, store, sensErrCon, adjoint,
                       warned = !is.null(maxattemps))

  SYM <- .nativeSym(paste0("solve_", as.character(model), "_batch"))
  .batchRun(model, preps, SYM, .batchDimnames(preps, SYM), names(conditions),
            .batchCores(cores, length(conditions)), onFailure, traceFile)
}


# Dimnames handed to the C++ side so it can label the freshly allocated arrays.
# Setting them in R afterwards duplicates the sensitivity array, which is the
# largest object in flight. Conditions that disagree get one entry each.
.batchDimnames <- function(preps, sym) {
  if (is.null(sym)) return(NULL)
  sens_nms <- lapply(preps, `[[`, "sens_col_names")
  if (length(unique(sens_nms)) == 1L)
    return(list(preps[[1L]]$variables, sens_nms[[1L]]))
  lapply(preps, function(p) list(p$variables, p$sens_col_names))
}


# Validate and marshal every condition. Serial R work, shared by solveODEBatch()
# and prepareBatch(). `warned` is TRUE when the caller has already warned about
# the old spelling `maxattemps`.
.batchPreps <- function(model, conditions, times, parms, tangent, hessian,
                        fixed, forcings, abstol, reltol, maxattempts, maxsteps,
                        hini, roottol, maxroot, cotangent = NULL, curvature = NULL,
                        keepStore = FALSE, store = NULL, sensErrCon = TRUE,
                        adjoint = NULL, warned = FALSE) {

  if (!is.list(conditions) || !length(conditions))
    stop("'conditions' must be a non-empty list", call. = FALSE)
  if (!all(vapply(conditions, is.list, logical(1))))
    stop("every element of 'conditions' must be a list of arguments", call. = FALSE)
  conditions <- .renameInConditions(conditions, warned)

  known <- c("times", "parms", "tangent", "hessian", "fixed", "forcings",
             "abstol", "reltol", "maxattempts", "maxsteps", "hini", "roottol",
             "maxroot", "cotangent", "curvature", "keepStore", "store",
             "sensErrCon", "adjoint")
  bad <- setdiff(unlist(lapply(conditions, names)), known)
  if (length(bad))
    stop("unknown per-condition argument(s): ", paste(unique(bad), collapse = ", "),
         "\n  Per-condition arguments are: ", paste(known, collapse = ", "),
         call. = FALSE)

  shared <- list(times = times, parms = parms, tangent = tangent,
                 hessian = hessian, fixed = fixed, forcings = forcings,
                 abstol = abstol, reltol = reltol, maxattempts = maxattempts,
                 maxsteps = maxsteps, hini = hini, roottol = roottol,
                 maxroot = maxroot, cotangent = cotangent, curvature = curvature,
                 keepStore = keepStore, store = store, sensErrCon = sensErrCon,
                 adjoint = adjoint)

  lapply(seq_along(conditions), function(i) {
    a <- utils::modifyList(shared, conditions[[i]])
    if (is.null(a$times) || is.null(a$parms))
      stop("condition ", i, " has no 'times' or no 'parms', and none was given ",
           "batch-wide", call. = FALSE)
    .odeCallArgs(model, a$times, a$parms, a$tangent, a$hessian, a$fixed,
                 a$forcings, a$abstol, a$reltol, a$maxattempts, a$maxsteps,
                 a$hini, a$roottol, a$maxroot, a$cotangent, a$curvature,
                 a$keepStore, a$store, a$sensErrCon, a$adjoint)
  })
}


# The value of a renamed argument. The old name warns and is used when given;
# giving both names is an error.
.renamedArg <- function(new, old, newGiven, newName, oldName) {
  if (is.null(old)) return(new)
  if (newGiven)
    stop("give '", newName, "' only; '", oldName, "' is its deprecated ",
         "spelling", call. = FALSE)
  warning("'", oldName, "' is deprecated; use '", newName, "'", call. = FALSE)
  old
}


# Per-condition lists with `maxattemps` renamed to `maxattempts`, warning once
# per call unless the caller has warned already.
.renameInConditions <- function(conditions, warned = FALSE) {
  old <- vapply(conditions, function(a) "maxattemps" %in% names(a), logical(1))
  if (!any(old)) return(conditions)
  both <- vapply(conditions, function(a) "maxattempts" %in% names(a), logical(1))
  if (any(old & both))
    stop("give 'maxattempts' only; 'maxattemps' is its deprecated spelling",
         call. = FALSE)
  if (!warned)
    warning("'maxattemps' is deprecated; use 'maxattempts'", call. = FALSE)
  conditions[old] <- lapply(conditions[old], function(a) {
    names(a)[names(a) == "maxattemps"] <- "maxattempts"
    a
  })
  conditions
}


# Resolve `traceFile` to one path per condition. A single path is a template:
# the condition's name (or index) goes in front of the extension, so a batch
# never has its conditions overwrite one another's trace.
.batchTraceFiles <- function(traceFile, k, nms) {
  if (is.null(traceFile)) return(vector("list", k))
  if (!is.character(traceFile) || !all(nzchar(traceFile)))
    stop("'traceFile' must be a non-empty character vector.", call. = FALSE)
  if (length(traceFile) == k) return(as.list(traceFile))
  if (length(traceFile) != 1L)
    stop("'traceFile' must be a single path or one per condition (",
         k, ").", call. = FALSE)
  tag <- if (!is.null(nms)) make.names(nms) else as.character(seq_len(k))
  ext <- sub("^.*(\\.[^.]+)$", "\\1", traceFile)
  if (identical(ext, traceFile)) ext <- ""
  as.list(paste0(sub("\\.[^.]+$", "", traceFile), "_", tag, ext))
}


# Fire the batch (or loop when the entry point is missing) and finish results.
.batchRun <- function(model, preps, SYM, dn, nms, nt, onFailure,
                      traceFile = NULL) {

  raw <- if (is.null(SYM)) {
    single <- .nativeSym(paste0("solve_", as.character(model)))
    if (is.null(single)) stop("Model not loaded. Run compile() first.", call. = FALSE)
    lapply(preps, function(p)
      do.call(.callSym, c(list(single), p$call_args,
                       list(list(p$variables, p$sens_col_names)))))
  } else {
    tryCatch(.callSym(SYM, lapply(preps, `[[`, "call_args"), as.integer(nt), dn),
             error = function(e) stop("ODE solver error: ", e$message, call. = FALSE))
  }

  # One trace per condition. The workers only fill the buffers; the files are
  # written here, on the R thread, after the parallel region.
  tf <- .batchTraceFiles(traceFile, length(preps), nms)
  out <- Map(function(r, p, f) .odeFinish(r, model, p, "silent", f),
             raw, preps, tf)
  names(out) <- nms

  rc <- vapply(out, function(r) {
    d <- r$diagnostics
    if (is.null(d)) 0L else as.integer(d$return_code)
  }, integer(1))
  if (any(rc != 0L) && onFailure != "silent") {
    who <- if (is.null(names(out))) which(rc != 0L) else names(out)[rc != 0L]
    msg <- paste0(sum(rc != 0L), " of ", length(rc),
                  " conditions did not complete: ", paste(who, collapse = ", "))
    switch(onFailure,
           stop = stop(msg, call. = FALSE),
           warn = warning(paste0(msg, "\n  Returning partial results."),
                          call. = FALSE, immediate. = TRUE))
  }
  # The C++ side reports the thread count it used; fall back to the request when
  # the batch entry is missing.
  if (is.null(attr(raw, "threads"))) attr(out, "threads") <- as.integer(nt)
  else attr(out, "threads") <- attr(raw, "threads")
  out
}


#' Will a Batch Run in Parallel?
#'
#' @description
#' `solveODEBatch()` degrades to a serial loop without saying so. This reports
#' the conditions it depends on, so a caller can warn instead of silently
#' losing its threads.
#'
#' @param model A model handle from [cppODE()] or [cvode()].
#' @return List with `symbol` (the batch entry resolved), `openmp` (cppDE was
#'   built with OpenMP), `modelOpenmp` (this model's shared object was), and
#'   `parallel` (all three). Two further fallbacks are decided at run time and
#'   are not visible here: inside a forked child, and inside an enclosing
#'   OpenMP region.
#' @seealso [solveODEBatch()]
#' @example inst/examples/batchAvailable.R
#' @export
batchAvailable <- function(model) {
  sym  <- .nativeSym(paste0("solve_", as.character(model), "_batch"))
  flags <- paste(attr(model, "compileArgs"), attr(model, "linkArgs"))
  mo <- grepl("fopenmp|openmp", flags, ignore.case = TRUE)
  ok <- !is.null(sym) && isTRUE(cvodeConfig$openmp_available) && mo
  list(symbol = !is.null(sym), openmp = isTRUE(cvodeConfig$openmp_available),
       modelOpenmp = mo, parallel = ok)
}


#' Prepare a Batch for Repeated Solving
#'
#' @description
#' Validates and converts a set of conditions once, so that repeated solves
#' with [solveBatch()], an optimiser evaluating the same model at new
#' parameters, only pass the numbers that changed. [solveODEBatch()] repeats
#' that work on every call.
#'
#' @inheritParams solveODEBatch
#' @return An object of class `"cppDEbatch"` for [solveBatch()].
#' @seealso [solveBatch()], [solveODEBatch()]
#' @example inst/examples/prepareBatch.R
#' @export
prepareBatch <- function(model, conditions,
                         times = NULL, parms = NULL,
                         tangent = NULL, hessian = NULL,
                         fixed = NULL, forcings = NULL,
                         abstol = 1e-6, reltol = 1e-6,
                         maxattempts = 50L, maxsteps = 1e6L,
                         hini = 0, roottol = 1e-6, maxroot = 1L,
                         cotangent = NULL, curvature = NULL,
                         keepStore = FALSE, store = NULL,
                         sensErrCon = TRUE, adjoint = NULL,
                         maxattemps = NULL) {

  maxattempts <- .renamedArg(maxattempts, maxattemps, !missing(maxattempts),
                             "maxattempts", "maxattemps")
  preps <- .batchPreps(model, conditions, times, parms, tangent, hessian,
                       fixed, forcings, abstol, reltol, maxattempts, maxsteps,
                       hini, roottol, maxroot, cotangent, curvature,
                       keepStore, store, sensErrCon, adjoint,
                       warned = !is.null(maxattemps))

  sym <- .nativeSym(paste0("solve_", as.character(model), "_batch"))
  structure(list(
    model     = model,
    names     = names(conditions),
    preps     = preps,
    required  = c(attr(model, "variables"), attr(model, "parameters")),
    sens_dim  = lapply(preps, function(p) dim(p$call_args[[3L]])),
    sym       = sym,
    dn        = .batchDimnames(preps, sym)),
    class = "cppDEbatch")
}


#' Solve a Prepared Batch
#'
#' @description
#' Re-solves the conditions of a [prepareBatch()] handle with new numbers.
#' Only `parms`, `tangent`, `hessian`, `cotangent`, `curvature`, `adjoint`
#' and `store` may change; anything else needs a fresh handle.
#'
#' Each of these arguments is a list with one element per condition, matched
#' to the prepared conditions by position; names are ignored. A `NULL` element
#' keeps that condition's prepared value, and `NULL` (default) for the whole
#' argument keeps all of them.
#'
#' @param handle A `"cppDEbatch"` object from [prepareBatch()].
#' @param parms List of named numeric vectors.
#' @param tangent,hessian,cotangent,curvature Lists of arrays, each as the
#'   argument of the same name in [solveODE()] and of the prepared shape.
#' @param adjoint List of [adjointControl()] objects.
#' @param store List of checkpoint stores, each the `$store` of an earlier
#'   reverse solve. Here a `NULL` element integrates that condition afresh.
#' @param cores,traceFile,onFailure As in [solveODEBatch()].
#' @return A list of [solveODE()] results, named as the prepared conditions.
#' @seealso [prepareBatch()], [solveODEBatch()]
#' @example inst/examples/solveBatch.R
#' @export
solveBatch <- function(handle, parms = NULL, tangent = NULL, hessian = NULL,
                       cotangent = NULL, curvature = NULL, adjoint = NULL,
                       store = NULL, cores = NULL, traceFile = NULL,
                       onFailure = c("stop", "warn", "silent")) {

  onFailure <- match.arg(onFailure)
  if (!inherits(handle, "cppDEbatch"))
    stop("'handle' must come from prepareBatch()", call. = FALSE)

  preps <- handle$preps
  K <- length(preps)
  chk <- function(x, what) {
    if (is.null(x)) return(NULL)
    if (!is.list(x) || length(x) != K)
      stop("'", what, "' must be a list with one element per condition",
           call. = FALSE)
    x
  }
  parms <- chk(parms, "parms"); tangent <- chk(tangent, "tangent")
  hessian <- chk(hessian, "hessian"); cotangent <- chk(cotangent, "cotangent")
  curvature <- chk(curvature, "curvature"); adjoint <- chk(adjoint, "adjoint")
  store <- chk(store, "store")
  n_states <- length(attr(handle$model, "variables"))
  if (!is.null(store)) {
    if (!identical(attr(handle$model, "derivMode"), "reverse"))
      stop("'store' belongs to a model compiled with derivMode = \"reverse\"",
           call. = FALSE)
    for (k in seq_len(K))
      if (!is.null(store[[k]]) && !inherits(store[[k]], "externalptr"))
        stop("condition ", k, ": 'store' must be the `store` element of an ",
             "earlier solve", call. = FALSE)
  }

  for (k in seq_len(K)) {
    if (!is.null(parms)) {
      pk <- parms[[k]]
      if (is.null(names(pk))) stop("'parms' elements must be named", call. = FALSE)
      miss <- setdiff(handle$required, names(pk))
      if (length(miss))
        stop("condition ", k, ": 'parms' missing ", paste(miss, collapse = ", "),
             call. = FALSE)
      v <- as.double(pk[handle$required])
      if (anyNA(v) || any(!is.finite(v)))
        stop("condition ", k, ": 'parms' must be finite", call. = FALSE)
      preps[[k]]$call_args[[2L]] <- v
    }
    # Assigning NULL into a list slot removes it and shifts every later
    # positional argument, so a NULL element means "keep the prepared value".
    if (!is.null(tangent) && !is.null(tangent[[k]])) {
      if (!identical(dim(tangent[[k]]), handle$sens_dim[[k]]))
        stop("condition ", k, ": tangent shape changed; call prepareBatch() again.",
             call. = FALSE)
      preps[[k]]$call_args[[3L]] <- tangent[[k]]
    }
    if (!is.null(hessian) && !is.null(hessian[[k]]))
      preps[[k]]$call_args[[4L]] <- hessian[[k]]
    # The cotangent holds the curvature and the adjoint control as attributes,
    # so a replacement takes the prepared ones over unless new ones are given
    # below.
    if (!is.null(cotangent) && !is.null(cotangent[[k]])) {
      sk  <- cotangent[[k]]
      old <- preps[[k]]$call_args[[15L]]
      if (is.null(old))
        stop("condition ", k, ": the handle was prepared without a 'cotangent'; ",
             "call prepareBatch() again.", call. = FALSE)
      if (!is.numeric(sk) || !identical(dim(sk), dim(old)))
        stop("condition ", k, ": cotangent shape changed; call prepareBatch() again.",
             call. = FALSE)
      storage.mode(sk) <- "double"
      attributes(sk) <- attributes(old)
      preps[[k]]$call_args[[15L]] <- sk
    }
    if (!is.null(curvature) && !is.null(curvature[[k]])) {
      ck  <- curvature[[k]]
      sk  <- preps[[k]]$call_args[[15L]]
      old <- attr(sk, "curvature")
      if (is.null(old))
        stop("condition ", k, ": the handle was prepared without a 'curvature'; ",
             "call prepareBatch() again.", call. = FALSE)
      if (!is.numeric(ck) || !identical(dim(ck), dim(old)))
        stop("condition ", k, ": curvature shape changed; call prepareBatch() again.",
             call. = FALSE)
      storage.mode(ck) <- "double"
      attr(sk, "curvature") <- ck
      preps[[k]]$call_args[[15L]] <- sk
    }
    # A new control replaces every attribute the prepared one wrote.
    if (!is.null(adjoint) && !is.null(adjoint[[k]])) {
      a <- adjoint[[k]]
      ck <- preps[[k]]$call_args[[15L]]
      if (!is.null(ck)) {
        for (nm in c("adjointGrid", "refine", "gradtol")) attr(ck, nm) <- NULL
        preps[[k]]$call_args[[15L]] <- .adjointApply(a, handle$model, ck, FALSE)
      }
      preps[[k]]$adjoint_wanted <- .adjointWanted(a)
    }
    # The store travels on `times`, the first argument; NULL drops it.
    if (!is.null(store)) attr(preps[[k]]$call_args[[1L]], "store") <- store[[k]]
  }

  .batchRun(handle$model, preps, handle$sym, handle$dn, handle$names,
            .batchCores(cores, K), onFailure, traceFile)
}


# Thread count for a batch of K conditions. Deliberately not read from
# OMP_NUM_THREADS: nothing keeps a process-wide variable in step with the
# BLAS thread count cppDE pins, so it is not a trustworthy source.
.batchCores <- function(cores, K) {
  if (!is.null(cores)) {
    cores <- as.integer(cores)
    if (is.na(cores) || cores < 1L) stop("'cores' must be a positive integer", call. = FALSE)
    return(min(cores, K))
  }
  n <- getOption("cppDE.cores", NULL)
  if (is.null(n)) n <- getOption("Ncpus", NULL)
  if (is.null(n)) n <- parallel::detectCores(logical = FALSE)
  if (!is.numeric(n) || !is.finite(n) || n < 1) n <- 1L
  min(as.integer(n), K)
}



#' Print Solver Diagnostics
#'
#' @description
#' Prints a summary of the solver diagnostics returned by [solveODE()].
#'
#' @param result A list returned by [solveODE()], containing a
#'   `diagnostics` element.
#'
#' @return The `diagnostics` list, invisibly; its fields are listed under
#'   Value in [solveODE()]. Without diagnostics, a message and `NULL`.
#'
#' @example inst/examples/diagnostics.R
#'
#' @export
diagnostics <- function(result) {
  UseMethod("diagnostics")
}

#' @export
diagnostics.default <- function(result) {
  diag <- result$diagnostics
  if (is.null(diag)) {
    message("No solver diagnostics available.")
    return(invisible(NULL))
  }

  rc <- diag$return_code
  # Return codes follow the SUNDIALS CVODE flag scheme
  # (see inst/include/cppde/cppde_return_codes.hpp).
  rc_text <- switch(as.character(rc),
                    "0"   = "Integration was successful.",
                    "1"   = "Reached TSTOP (requested stop time).",
                    "2"   = "Root function returned a root.",
                    "-1"  = "Too much work: maximum number of integration steps exceeded (CV_TOO_MUCH_WORK).",
                    "-2"  = "Too much accuracy requested for the requested step (CV_TOO_MUCH_ACC).",
                    "-3"  = "Error test failures too numerous (CV_ERR_FAILURE).",
                    "-4"  = "Convergence test failures too numerous / no progress (CV_CONV_FAILURE).",
                    "-5"  = "Linear solver initialisation failed (CV_LINIT_FAIL).",
                    "-6"  = "Linear solver setup failed unrecoverably (CV_LSETUP_FAIL).",
                    "-7"  = "Linear solver solve failed unrecoverably (CV_LSOLVE_FAIL).",
                    "-8"  = "RHS function failed unrecoverably (CV_RHSFUNC_FAIL).",
                    "-9"  = "RHS function failed at the first call (CV_FIRST_RHSFUNC_ERR).",
                    "-10" = "RHS function had repeated recoverable errors (CV_REPTD_RHSFUNC_ERR).",
                    "-11" = "RHS function had a recoverable error that could not be handled (CV_UNREC_RHSFUNC_ERR).",
                    "-12" = "Root function failed unrecoverably (CV_RTFUNC_FAIL).",
                    "-13" = "Nonlinear solver initialisation failed (CV_NLS_INIT_FAIL).",
                    "-14" = "Nonlinear solver setup failed (CV_NLS_SETUP_FAIL).",
                    "-15" = "Inequality constraint check failed (CV_CONSTR_FAIL).",
                    "-16" = "Nonlinear solver failed (CV_NLS_FAIL).",
                    "-20" = "Memory allocation request failed (CV_MEM_FAIL).",
                    "-21" = "Integrator memory is NULL (CV_MEM_NULL).",
                    "-22" = "Illegal input provided (CV_ILL_INPUT).",
                    "-23" = "Integrator memory was not allocated (CV_NO_MALLOC).",
                    "-24" = "Bad k value in CVodeGetDky (CV_BAD_K).",
                    "-25" = "Bad t value in CVodeGetDky (CV_BAD_T).",
                    "-26" = "Bad dky argument in CVodeGetDky (CV_BAD_DKY).",
                    "-27" = "Output time too close to initial time (CV_TOO_CLOSE).",
                    "-99" = "Unrecognised error (CV_UNRECOGNIZED_ERR).",
                    paste0("Unknown return code: ", rc)
  )

  label <- if (identical(diag$backend, "cvode")) {
    paste0("CVODE: ", toupper(diag$method %||% "bdf"))
  } else {
    switch(diag$method %||% "bdf",
           bdf         = if (isFALSE(diag$useNDF)) "BDF" else "NDF",
           adams       = "ADAMS",
           rosenbrock4 = "ROSENBROCK4",
           tsit5       = "TSIT5",
           toupper(diag$method))
  }
  label <- paste(label, "solver statistics")
  total <- 69L
  inner <- paste0("  ", label, "  ")
  dash_n <- max(total - nchar(inner), 2L)
  left  <- dash_n %/% 2L
  right <- dash_n - left
  cat(strrep("-", left), inner, strrep("-", right), "\n", sep = "")
  cat(sprintf("  Return code                  : %d\n", rc))
  cat(sprintf("  Message                      : %s\n", rc_text))
  cat("---------------------------------------------------------------------\n")
  cat(sprintf("  Accepted steps               : %d\n", diag$accepted))
  cat(sprintf("  Rejected steps               : %d\n", diag$rejected))
  cat(sprintf("  Function evaluations         : %d\n", diag$fevals))
  cat(sprintf("  Jacobian evaluations         : %d\n", diag$jevals))
  cat(sprintf("  LU factorizations            : %d\n", diag$setups))
  cat(sprintf("  Last step size (successful)  : %g\n", diag$last_dt))
  cat(sprintf("  Last method order            : %d\n", diag$last_order))
  cat(sprintf("  Time reached                 : %g\n", diag$t_reached))
  cat("---------------------------------------------------------------------\n")

  invisible(diag)
}
