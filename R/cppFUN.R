#' Compile Algebraic Functions with Optional Derivatives
#'
#' Generates and compiles C++ code that evaluates a system of algebraic
#' expressions \eqn{y = g(x, p)} on one or more rows of input, with
#' optional first- and second-order derivatives. There is no time
#' integration; the principal use cases are observation maps for
#' likelihood-based inference and reparametrisation Jacobians for
#' [solveODE()]. `derivMode` selects which derivative products are built.
#' The chain rule is available through the optional arguments `tangentX`,
#' `tangentP`, `hessianX` and `hessianP`, the tangent and the Hessian in
#' \eqn{\theta} of the variables and parameters. Every entry runs compiled
#' code: an object built
#' with `compile = FALSE` is evaluable only after [compile()]. See
#' `vignette("Methods", package = "cppDE")` for the computational paths and
#' the pass-through convention for unmodelled inputs.
#'
#' @param eqns Named character vector or list of algebraic expressions.
#'   Names define the output variables; defaults to `f1`, `f2`, ... when
#'   unnamed.
#' @param variables Character vector of variable names supplied per
#'   observation. Defaults to all symbols in `eqns` not in `parameters`.
#' @param parameters Character vector of parameter names (constant
#'   across observations).
#' @param fixed Optional character vector of symbols excluded from
#'   derivative computation.
#' @param modelname Optional base name for generated C++ symbols and
#'   files.
#' @param outdir Directory for generated C++ source files. Default
#'   `tempdir()`.
#' @param compile Logical. Compile and load the generated C++ code.
#'   Default `FALSE`.
#' @param verbose Logical. Print progress messages.
#' @param convenient Logical. Return wrappers that accept named
#'   arguments rather than the low-level `(vars, params)` signature.
#' @param deriv Logical. Generate first-order derivative entry points.
#' @param deriv2 Logical. Generate Hessian entry points; implies
#'   `deriv = TRUE`.
#' @param derivMode Which derivative products to build, any of `"forward"`
#'   (default), `"reverse"` and `"forward-reverse"`.
#'   * `"forward"`: forward-mode AD on `cppde::dual`, delivering `jac`,
#'     `hess`, `evaluate` and `evaluateBatch`.
#'   * `"reverse"`: the vector-Jacobian product `vjp`, differentiated at
#'     code-generation time.
#'   * `"forward-reverse"`: `vjp` and its derivative along a tangent, also
#'     differentiated at code-generation time.
#'   A mode not named is not generated and costs no compile time.
#'
#' @return A list with components `func`, `jac`, `hess`, `evaluate`,
#'   `evaluateBatch` and `vjp`, each `NULL` when not generated.
#'   `jac`, `hess`, `evaluate` and `evaluateBatch` need `"forward"`.
#'   `evaluate(..., tangentX, tangentP, hessianX, hessianP, deriv2)` returns
#'   `y`, `tangent` and, with `deriv2 = TRUE`, `hessian`: the outputs' tangent
#'   and Hessian in \eqn{\theta}. `evaluateBatch(sets, cores, deriv2)` runs
#'   `evaluate` over a list of argument lists in one call.
#'
#'   `vjp` needs `"reverse"`. `vjp(vars, params, cotangent)` contracts the
#'   Jacobian with a cotangent of the outputs, `[n_obs, n_out]` or
#'   `[n_obs, n_out, n_seed]`, and returns `y`, `cotangentX` and `cotangentP`;
#'   `cotangentP` sums over observations because the parameters are shared
#'   across them. Given `tangentX`, `tangentP` or `curvature`, the derivative
#'   of the cotangent along the tangent, the same call runs forward-reverse and
#'   adds `curvatureX` and `curvatureP`; this needs `"forward-reverse"`. Carries attributes `equations`,
#'   `variables`, `parameters`, `fixed`, `modelname`, `srcfile` and
#'   `derivMode`.
#'
#' @seealso [compile()] for compilation; [cppODE()] and [cvode()] for ODE
#'   integration; `vignette("Methods", package = "cppDE")`.
#' @example inst/examples/cppFUN.R
#' @export
cppFUN <- function(eqns, variables = getSymbols(eqns, omit = parameters), parameters = NULL,
                   fixed = NULL, modelname = NULL, outdir = tempdir(), compile = FALSE,
                   verbose = FALSE, convenient = TRUE, deriv = TRUE, deriv2 = FALSE,
                   derivMode = "forward") {

  derivMode <- matchDerivMode(derivMode, c("forward", "reverse", "forward-reverse"))
  if (deriv2 && !deriv) { warning("deriv2 requires deriv. Setting deriv = TRUE."); deriv <- TRUE }
  emit_deriv <- deriv || deriv2
  use_ad     <- emit_deriv && "forward" %in% derivMode
  use_vjpfr  <- emit_deriv && "forward-reverse" %in% derivMode
  use_vjp    <- emit_deriv && ("reverse" %in% derivMode || use_vjpfr)
  ## Second order is a forward-mode facility. Asking for it with only the
  ## reverse direction would silently return no Hessian.
  if (deriv2 && !use_ad)
    stop("deriv2 = TRUE has no reverse counterpart; add \"forward\" to ",
         "derivMode.", call. = FALSE)

  # The symbol arguments name the same symbols as the equations and are checked
  # with them.
  checkSymbolNames(eqns, variables, parameters, fixed)

  outnames <- names(eqns) %||% paste0("f", seq_along(eqns))
  if (!is.null(fixed)) { variables <- setdiff(variables, fixed); parameters <- union(parameters, fixed) }
  innames <- variables; diff_params <- setdiff(parameters, fixed); diff_syms <- c(variables, diff_params)
  if (!dir.exists(outdir)) stop("outdir does not exist: ", outdir)
  modelname <- modelname %||% randomModelname("f")
  modelname <- unique_modelname(modelname)

  # --- C++ codegen ---
  codegen <- get_codegen_cppFUN_py()
  cpp_file <- file.path(outdir, paste0(modelname, ".cpp"))
  if (file.exists(cpp_file)) message("Overwriting: ", normalizePath(cpp_file, "/", FALSE))
  codegen$generate_fun_cpp(exprs = setNames(as.list(eqns), outnames), variables = as.list(variables),
                           parameters = as.list(parameters),
                           ad = use_ad, deriv2 = deriv2, vjp = use_vjp,
                           vjp_fr = use_vjpfr,
                           modelname = modelname, outdir = normalizePath(outdir, "/", FALSE), version = as.character(utils::packageVersion("cppDE")))

  # --- Instance state and thin wrappers ---
  ## All the engine needs and nothing else. The returned closures carry only
  ## this environment, so a per-condition instance costs data, not code.
  st <- list2env(list(innames = innames, parameters = parameters,
                      outnames = outnames, modelname = modelname,
                      diff_syms = diff_syms, fixed = intersect(fixed, parameters),
                      use_ad = use_ad, use_vjp = use_vjp, use_vjpfr = use_vjpfr),
                 parent = emptyenv())

  fun_impl      <- function(...) .fun_impl(st, ...)
  ## Entries of the requested directions only.
  fwd           <- emit_deriv && use_ad
  jac_impl      <- if (deriv  && fwd) function(...) .jac_impl(st, ...)
  hess_impl     <- if (deriv2 && fwd) function(...) .hess_impl(st, ...)
  evaluate_impl <- if (fwd)     function(...) .evaluate_impl(st, ...)
  evaluateBatch_impl <- if (fwd) function(...) .evaluateBatch_impl(st, ...)
  vjp_impl      <- if (use_vjp) function(...) .vjp_impl(st, ...)

  # --- Output ---
  ## Installed with keep.source, cppFUN's body carries srcrefs, so the wrappers
  ## built here would hand the caller a copy each.
  .stripSource(environment())
  outfn <- list(
    func     = if (convenient) .makeFunWrapper(st, fun_impl) else fun_impl,
    jac      = if (convenient) .makeDerivWrapper(st, jac_impl, FALSE) else jac_impl,
    hess     = if (convenient) .makeDerivWrapper(st, hess_impl, TRUE) else hess_impl,
    evaluate = if (convenient) .makeEvalWrapper(st, evaluate_impl) else evaluate_impl,
    evaluateBatch = evaluateBatch_impl,
    vjp      = vjp_impl
  )
  attr(outfn, "equations") <- eqns; attr(outfn, "variables") <- variables; attr(outfn, "parameters") <- parameters
  attr(outfn, "fixed") <- fixed; attr(outfn, "modelname") <- modelname; attr(outfn, "srcfile") <- normalizePath(cpp_file, "/", FALSE)
  attr(outfn, "derivMode") <- derivMode
  for (nm in c("func", "jac", "hess", "evaluate", "evaluateBatch", "vjp")) {
    if (!is.null(outfn[[nm]])) {
      attr(outfn[[nm]], "modelname") <- modelname
      attr(outfn[[nm]], "srcfile")   <- attr(outfn, "srcfile")
    }
  }
  if (compile) compile(outfn, verbose = verbose)
  outfn
}


# ============================================================================
# Runtime engine
# ============================================================================

# Package level, so a model with one instance per condition carries this code
# once rather than once per condition; `st` holds the per-instance data.


# .Call needs plain doubles; a matrix already is one, so this is a no-op there.
.asdbl <- function(x) if (is.double(x)) x else as.double(x)

# Stops for an entry point that is not loaded.
.notCompiled <- function(st) {
  stop("cppFUN object '", st$modelname, "' is not compiled; call compile().",
       call. = FALSE)
}

# Parameters excluded from differentiation: construction- and call-time `fixed`.
.fixedAt <- function(st, fixed) {
  if (is.null(fixed)) st$fixed else union(st$fixed, intersect(fixed, st$parameters))
}

# --- Input validation ---

.checkInputs <- function(st, vars, params, attach = FALSE) {
  n_obs <- if (is.matrix(vars) || is.data.frame(vars)) nrow(vars)
  else if (is.vector(vars) && !is.list(vars)) length(vars) / max(length(st$innames), 1L) else 1L
  n_obs <- max(as.integer(n_obs), 1L); extra_vars <- extra_params <- NULL

  if (!length(st$innames)) {
    M <- matrix(0, n_obs, 0)
    if (attach && !is.null(vars) && (is.matrix(vars) || is.data.frame(vars)) && ncol(vars) > 0) {
      extra_vars <- as.matrix(vars); n_obs <- nrow(extra_vars)
    }
  } else {
    if (is.null(vars)) stop("Variables defined but 'vars' is NULL.")
    if (is.vector(vars) && !is.list(vars)) vars <- matrix(vars, ncol = length(st$innames), dimnames = list(NULL, st$innames))
    cn <- colnames(vars)
    if (is.null(cn)) { colnames(vars) <- st$innames; cn <- st$innames }
    # match() does the same test as setdiff() without the unique()/coerce pass,
    # and an already-matching column set needs no subset copy at all.
    if (is.matrix(vars) && identical(cn, st$innames)) {
      M <- vars
    } else {
      hit <- match(st$innames, cn, 0L)
      if (any(hit == 0L))
        stop("Missing variables: ", paste(st$innames[hit == 0L], collapse = ", "))
      M <- vars[, hit, drop = FALSE]
    }
    n_obs <- nrow(M)
    if (attach) { ex <- setdiff(cn, st$innames); if (length(ex)) extra_vars <- vars[, ex, drop = FALSE] }
  }

  if (!length(st$parameters)) {
    p <- numeric(0); if (attach && length(params)) extra_params <- params
  } else {
    pn <- names(params)
    if (is.null(pn)) stop("params must be named.")
    hit <- match(st$parameters, pn, 0L)
    if (any(hit == 0L))
      stop("Missing parameters: ", paste(st$parameters[hit == 0L], collapse = ", "))
    p <- params[hit]
    if (attach) { ex <- setdiff(pn, st$parameters); if (length(ex)) extra_params <- params[ex] }
  }
  # M stays as R holds it, [n_obs, n_vars]: the generated entries index it as
  # obs + n_obs * j, so no per-condition matrix has to be transposed.
  list(M = M, p = p, n_obs = n_obs, extra_vars = extra_vars, extra_params = extra_params)
}

# --- Attach helpers (pass-through of unmodelled inputs) ---
.abind3 <- function(a, b) { d <- dim(a); db <- dim(b); r <- array(0, c(d[1], d[2]+db[2], d[3]), list(dimnames(a)[[1]], c(dimnames(a)[[2]], dimnames(b)[[2]]), dimnames(a)[[3]])); r[,1:d[2],] <- a; r[,d[2]+1:db[2],] <- b; r }
.abind4 <- function(a, b) { d <- dim(a); db <- dim(b); r <- array(0, c(d[1], d[2]+db[2], d[3], d[4]), list(dimnames(a)[[1]], c(dimnames(a)[[2]], dimnames(b)[[2]]), dimnames(a)[[3]], dimnames(a)[[4]])); r[,1:d[2],,] <- a; r[,d[2]+1:db[2],,] <- b; r }

.attachExtras <- function(res, n_obs, ev, ep, type) {
  if (is.null(ev) && is.null(ep)) return(res)
  if (type == "fun") {
    if (!is.null(ev)) res <- cbind(res, ev)
    if (!is.null(ep)) res <- cbind(res, matrix(rep(ep, each = n_obs), n_obs, length(ep), dimnames = list(NULL, names(ep))))
  } else if (type == "jac") {
    cs <- dimnames(res)[[3]]; ncs <- length(cs)
    if (!is.null(ev)) res <- .abind3(res, array(0, c(n_obs, ncol(ev), ncs), list(NULL, colnames(ev), cs)))
    if (!is.null(ep)) { np <- length(ep); pn <- names(ep); d <- dim(res); new <- array(0, c(d[1], d[2]+np, d[3]+np), list(dimnames(res)[[1]], c(dimnames(res)[[2]], pn), c(dimnames(res)[[3]], pn))); new[,1:d[2],1:d[3]] <- res; for (k in seq_len(np)) new[,d[2]+k,d[3]+k] <- 1; res <- new }
  } else {
    cs <- dimnames(res)[[3]]; ncs <- length(cs)
    if (!is.null(ev)) res <- .abind4(res, array(0, c(n_obs, ncol(ev), ncs, ncs), list(NULL, colnames(ev), cs, cs)))
    if (!is.null(ep)) { np <- length(ep); pn <- names(ep); d <- dim(res); new <- array(0, c(d[1], d[2]+np, d[3]+np, d[4]+np), list(dimnames(res)[[1]], c(dimnames(res)[[2]], pn), c(dimnames(res)[[3]], pn), c(dimnames(res)[[4]], pn))); new[,1:d[2],1:d[3],1:d[4]] <- res; res <- new }
  }
  res
}

# --- Chain-rule helpers ---

# Pull theta names from any seed; raise if seeds disagree.
.resolveTheta <- function(tangentX, tangentP, hessianX = NULL, hessianP = NULL) {
  cands <- list()
  if (!is.null(tangentX))  cands[[length(cands) + 1L]] <- dimnames(tangentX)[[3]]
  if (!is.null(tangentP))  cands[[length(cands) + 1L]] <- colnames(tangentP)
  if (!is.null(hessianX)) cands[[length(cands) + 1L]] <- dimnames(hessianX)[[3]]
  if (!is.null(hessianP)) cands[[length(cands) + 1L]] <- dimnames(hessianP)[[2]]
  cands <- cands[lengths(cands) > 0]
  if (!length(cands)) return(character(0))
  theta <- cands[[1]]
  for (k in seq_along(cands)[-1])
    if (!identical(theta, cands[[k]]) && !setequal(theta, cands[[k]]))
      stop("Seed theta names disagree across tangentX/tangentP/hessianX/hessianP")
  theta
}

# Align seeds onto the function's internal (vars, params) order, flat and ready
# for the dual-mode .C() entries. An absent hessianX or hessianP goes in as a
# length-zero vector, which the C side guards with `has_dX2` / `has_dP2`.
.alignSeedsDual <- function(st, tangentX, tangentP, hessianX, hessianP, n_obs,
                            theta, fixed_rt) {
  n_vars <- length(st$innames); n_params <- length(st$parameters); n_theta <- length(theta)
  dX_arr <- array(0, c(n_obs, n_vars, n_theta))
  dP_mat <- matrix(0, n_params, n_theta)
  if (n_theta > 0) {
    if (!is.null(tangentX) && n_vars > 0) {
      idx <- match(st$innames, dimnames(tangentX)[[2]]); pres <- !is.na(idx)
      # Character selection over n_theta names costs a lookup per element.
      tX <- if (identical(dimnames(tangentX)[[3]], theta)) TRUE else theta
      if (any(pres)) dX_arr[, pres, ] <- tangentX[, idx[pres], tX, drop = FALSE]
    }
    if (!is.null(tangentP) && n_params > 0) {
      idx <- match(st$parameters, rownames(tangentP)); pres <- !is.na(idx)
      tP <- if (identical(colnames(tangentP), theta)) TRUE else theta
      if (any(pres)) dP_mat[pres, ] <- tangentP[idx[pres], tP, drop = FALSE]
    }
    if (length(fixed_rt) && n_params > 0) {
      fp <- match(fixed_rt, st$parameters); fp <- fp[!is.na(fp)]
      if (length(fp)) dP_mat[fp, ] <- 0
    }
  }
  has_dX2 <- !is.null(hessianX) && n_theta > 0 && n_vars > 0
  has_dP2 <- !is.null(hessianP) && n_theta > 0 && n_params > 0
  dX2_arr <- if (has_dX2) {
    r <- array(0, c(n_obs, n_vars, n_theta, n_theta))
    idx <- match(st$innames, dimnames(hessianX)[[2]]); pres <- !is.na(idx)
    dn <- dimnames(hessianX)
    t3 <- if (identical(dn[[3]], theta)) TRUE else theta
    t4 <- if (identical(dn[[4]], theta)) TRUE else theta
    if (any(pres)) r[, pres, , ] <- hessianX[, idx[pres], t3, t4, drop = FALSE]
    r
  } else double(0)
  dP2_arr <- if (has_dP2) {
    r <- array(0, c(n_params, n_theta, n_theta))
    idx <- match(st$parameters, dimnames(hessianP)[[1]]); pres <- !is.na(idx)
    dn <- dimnames(hessianP)
    t2 <- if (identical(dn[[2]], theta)) TRUE else theta
    t3 <- if (identical(dn[[3]], theta)) TRUE else theta
    if (any(pres)) r[pres, , ] <- hessianP[idx[pres], t2, t3, drop = FALSE]
    if (length(fixed_rt)) {
      fp <- match(fixed_rt, st$parameters); fp <- fp[!is.na(fp)]
      if (length(fp)) r[fp, , ] <- 0
    }
    r
  } else double(0)
  # Kept as arrays: .C() and .Call() both read the REALSXP directly, and
  # as.double() on an array is a full copy of the largest object in play.
  list(tangentX = dX_arr, tangentP = dP_mat,
       hessianX = dX2_arr, hessianP = dP2_arr,
       has_dX2 = as.integer(has_dX2), has_dP2 = as.integer(has_dP2))
}

# Identity seeds for raw J/H: tangentX = I on vars, tangentP = I on params,
# over the combined basis c(st$innames, st$parameters). A `fixed` parameter
# seeds zero but keeps its column.
.identitySeedsRaw <- function(st, n_obs, fixed_rt) {
  n_vars <- length(st$innames); n_params <- length(st$parameters)
  theta_full <- c(st$innames, st$parameters)
  n_theta <- length(theta_full)
  tangentX <- array(0, c(n_obs, n_vars, n_theta), dimnames = list(NULL, st$innames, theta_full))
  if (n_vars > 0) for (i in seq_along(st$innames)) tangentX[, st$innames[i], st$innames[i]] <- 1
  tangentP <- matrix(0, n_params, n_theta, dimnames = list(st$parameters, theta_full))
  for (i in seq_along(st$parameters))
    if (!(st$parameters[i] %in% fixed_rt))
      tangentP[st$parameters[i], st$parameters[i]] <- 1
  list(tangentX = tangentX, tangentP = tangentP, theta = theta_full)
}

# --- Core implementations (outputs time-first) ---

.fun_impl <- function(st, vars, params = numeric(0), attach.input = FALSE, fixed = NULL) {
  chk <- .checkInputs(st, vars, params, attach.input); M <- chk$M; p <- chk$p; n_obs <- chk$n_obs
  # The _c entry takes its arguments by reference; the .C() entry, used when
  # no _c symbol exists, copies every argument in and every result out.
  symc <- .nativeSym(paste0(st$modelname, "_eval_c"))
  sym <- if (is.null(symc)) .nativeSym(paste0(st$modelname, "_eval")) else NULL
  if (!is.null(symc)) {
    res <- .callSym(symc, .asdbl(M), .asdbl(p), as.integer(n_obs))
    dimnames(res) <- list(NULL, st$outnames)
  } else if (!is.null(sym)) {
    out <- .cSym(sym, x = as.double(M), y = double(length(st$outnames) * n_obs), p = as.double(p), n = as.integer(n_obs), k = as.integer(length(st$innames)), l = as.integer(length(st$outnames)))
    res <- matrix(out$y, n_obs, length(st$outnames), dimnames = list(NULL, st$outnames))
  } else .notCompiled(st)
  .attachExtras(res, n_obs, chk$extra_vars, chk$extra_params, "fun")
}

# --- Reverse path ---

# Vector-Jacobian product: cotangents of the variables and parameters from a
# cotangent of the outputs, [n_obs, n_out] or [n_obs, n_out, n_seed].
# cotangentP sums over observations because the parameters are shared.
# With a tangent or a curvature the same contraction runs over a dual, forward
# over reverse, and adds the curvature of the inputs.
.vjp_impl <- function(st, vars, params = numeric(0), cotangent,
                      tangentX = NULL, tangentP = NULL, curvature = NULL) {
  chk <- .checkInputs(st, vars, params); M <- chk$M; p <- chk$p; n_obs <- chk$n_obs
  n_vars <- length(st$innames); n_params <- length(st$parameters)
  n_out  <- length(st$outnames)

  w <- cotangent
  if (is.null(dim(w))) w <- matrix(w, n_obs, n_out)
  n_seed <- if (length(dim(w)) == 3L) dim(w)[3L] else 1L
  if (dim(w)[1L] != n_obs || dim(w)[2L] != n_out)
    stop("cotangent must be [n_obs, n_out] or [n_obs, n_out, n_seed].")

  dn_x <- list(NULL, st$innames, NULL)
  dn_p <- list(st$parameters, NULL)
  second <- !is.null(tangentX) || !is.null(tangentP) || !is.null(curvature)
  if (second) return(.vjpDual(st, M, p, w, tangentX, tangentP, curvature,
                              n_obs, n_seed, dn_x, dn_p))

  funsym <- paste0(st$modelname, "_vjp")
  symc <- .nativeSym(paste0(funsym, "_c"))
  if (!is.null(symc)) {
    r <- .callSym(symc, .asdbl(M), .asdbl(p), .asdbl(w),
                  as.integer(n_obs), as.integer(n_seed))
    y <- r[[1L]]; dimnames(y) <- list(NULL, st$outnames)
    cx <- r[[2L]]; dimnames(cx) <- dn_x
    cp <- r[[3L]]; dimnames(cp) <- dn_p
    return(list(y = y, cotangentX = cx, cotangentP = cp))
  }

  sym <- .nativeSym(funsym)
  if (is.null(sym)) .notCompiled(st)
  out <- .cSym(sym,
            x        = as.double(M),
            p        = as.double(p),
            w        = as.double(w),
            y        = double(n_out * n_obs),
            wx       = double(n_obs * n_vars * n_seed),
            wp       = double(n_params * n_seed),
            n_obs    = as.integer(n_obs),
            n_vars   = as.integer(n_vars),
            n_params = as.integer(n_params),
            n_out    = as.integer(n_out),
            n_seed   = as.integer(n_seed))
  list(y          = matrix(out$y, n_obs, n_out, dimnames = list(NULL, st$outnames)),
       cotangentX = array(out$wx, c(n_obs, n_vars, n_seed), dimnames = dn_x),
       cotangentP = matrix(out$wp, n_params, n_seed, dimnames = dn_p))
}

# Forward over reverse: the tangents of the inputs and the curvature of the
# cotangent go in, and both terms of the derivative of the cotangent contraction
# come back in one pass. Compiled entry only.
.vjpDual <- function(st, M, p, w, tangentX, tangentP, curvature,
                     n_obs, n_seed, dn_x, dn_p) {
  n_dir <- if (!is.null(tangentX)) utils::tail(dim(tangentX), 1L)
           else if (!is.null(tangentP)) utils::tail(dim(tangentP), 1L)
           else utils::tail(dim(curvature), 1L)

  if (!isTRUE(st$use_vjpfr))
    stop("tangents need derivMode \"forward-reverse\".", call. = FALSE)
  symc <- .nativeSym(paste0(st$modelname, "_vjp_ad_c"))
  if (is.null(symc)) .notCompiled(st)
  r <- .callSym(symc, .asdbl(M), .asdbl(p), .asdbl(w),
                if (is.null(tangentX)) NULL else .asdbl(tangentX),
                if (is.null(tangentP)) NULL else .asdbl(tangentP),
                if (is.null(curvature)) NULL else .asdbl(curvature),
                as.integer(n_obs), as.integer(n_seed), as.integer(n_dir))
  y <- r[[1L]]; dimnames(y) <- list(NULL, st$outnames)
  cx <- r[[2L]]; dimnames(cx) <- dn_x
  cp <- r[[3L]]; dimnames(cp) <- dn_p
  list(y = y, cotangentX = cx, cotangentP = cp,
       curvatureX = r[[4L]], curvatureP = r[[5L]])
}

# --- Dual-path .C() helpers ---

# Single-dual AD pass; returns list(y, dy) with dy of shape [n_obs, n_out, n_theta]
# (or [..., 0] if n_theta == 0). Used by the dual path for jac and for
# evaluate() at deriv2 = FALSE.
.call_eval_ad <- function(st, M, p, dX_seed, dP_seed, n_obs, theta) {
  funsym <- paste0(st$modelname, "_eval_ad"); n_out <- length(st$outnames); n_theta <- length(theta)
  symc <- .nativeSym(paste0(funsym, "_c"))
  if (!is.null(symc)) {
    r <- .callSym(symc, .asdbl(M), .asdbl(p), .asdbl(dX_seed), .asdbl(dP_seed),
               as.integer(n_obs), as.integer(n_theta))
    y <- r[[1L]]; dimnames(y) <- list(NULL, st$outnames)
    dy <- if (n_theta > 0) {
      d <- r[[2L]]; dimnames(d) <- list(NULL, st$outnames, theta); d
    } else array(0, c(n_obs, n_out, 0L), list(NULL, st$outnames, NULL))
    return(list(y = y, dy = dy))
  }
  sym <- .nativeSym(funsym)
  if (is.null(sym)) .notCompiled(st)
  out <- .cSym(sym,
            x        = as.double(M),
            p        = as.double(p),
            dX       = dX_seed,
            dP       = dP_seed,
            y        = double(n_out * n_obs),
            dy       = double(n_out * max(n_theta, 1L) * n_obs),
            n_obs    = as.integer(n_obs),
            n_vars   = as.integer(length(st$innames)),
            n_params = as.integer(length(st$parameters)),
            n_out    = as.integer(n_out),
            n_theta  = as.integer(n_theta))
  y <- matrix(out$y, n_obs, n_out, dimnames = list(NULL, st$outnames))
  dy <- if (n_theta > 0)
    array(out$dy[seq_len(n_out * n_theta * n_obs)], c(n_obs, n_out, n_theta),
          list(NULL, st$outnames, theta))
  else array(0, c(n_obs, n_out, 0L), list(NULL, st$outnames, NULL))
  list(y = y, dy = dy)
}

# Nested-dual AD pass; returns list(y, dy, d2y).
.call_eval_ad2 <- function(st, M, p, aligned, n_obs, theta) {
  funsym <- paste0(st$modelname, "_eval_ad2"); n_out <- length(st$outnames); n_theta <- length(theta)
  symc <- .nativeSym(paste0(funsym, "_c"))
  if (!is.null(symc)) {
    r <- .callSym(symc, .asdbl(M), .asdbl(p), .asdbl(aligned$tangentX), .asdbl(aligned$tangentP),
               if (aligned$has_dX2 != 0L) .asdbl(aligned$hessianX) else NULL,
               if (aligned$has_dP2 != 0L) .asdbl(aligned$hessianP) else NULL,
               as.integer(n_obs), as.integer(n_theta))
    y <- r[[1L]]; dimnames(y) <- list(NULL, st$outnames)
    if (n_theta > 0) {
      dy <- r[[2L]];  dimnames(dy)  <- list(NULL, st$outnames, theta)
      d2y <- r[[3L]]; dimnames(d2y) <- list(NULL, st$outnames, theta, theta)
    } else {
      dy  <- array(0, c(n_obs, n_out, 0L), list(NULL, st$outnames, NULL))
      d2y <- array(0, c(n_obs, n_out, 0L, 0L), list(NULL, st$outnames, NULL, NULL))
    }
    return(list(y = y, dy = dy, d2y = d2y))
  }
  sym <- .nativeSym(funsym)
  if (is.null(sym)) .notCompiled(st)
  out <- .cSym(sym,
            x        = as.double(M),
            p        = as.double(p),
            dX       = as.double(aligned$tangentX),
            dP       = as.double(aligned$tangentP),
            dX2_in   = as.double(aligned$hessianX),
            dP2_in   = as.double(aligned$hessianP),
            has_dX2  = aligned$has_dX2,
            has_dP2  = aligned$has_dP2,
            y        = double(n_out * n_obs),
            dy       = double(n_out * max(n_theta, 1L) * n_obs),
            d2y      = double(n_out * max(n_theta, 1L)^2 * n_obs),
            n_obs    = as.integer(n_obs),
            n_vars   = as.integer(length(st$innames)),
            n_params = as.integer(length(st$parameters)),
            n_out    = as.integer(n_out),
            n_theta  = as.integer(n_theta))
  y <- matrix(out$y, n_obs, n_out, dimnames = list(NULL, st$outnames))
  if (n_theta > 0) {
    dy  <- array(out$dy[seq_len(n_out * n_theta * n_obs)], c(n_obs, n_out, n_theta),
                 list(NULL, st$outnames, theta))
    d2y <- array(out$d2y[seq_len(n_out * n_theta^2 * n_obs)], c(n_obs, n_out, n_theta, n_theta),
                 list(NULL, st$outnames, theta, theta))
  } else {
    dy  <- array(0, c(n_obs, n_out, 0L), list(NULL, st$outnames, NULL))
    d2y <- array(0, c(n_obs, n_out, 0L, 0L), list(NULL, st$outnames, NULL, NULL))
  }
  list(y = y, dy = dy, d2y = d2y)
}

# --- Public derivative implementations ---

.jac_impl <- function(st, vars, params = numeric(0), tangentX = NULL, tangentP = NULL,
                                attach.input = FALSE, fixed = NULL) {
  if (is.null(tangentX)) tangentX <- attr(vars, "deriv")
  if (is.null(tangentP)) tangentP <- attr(params, "deriv")
  has_seeds <- !is.null(tangentX) || !is.null(tangentP)
  chk <- .checkInputs(st, vars, params, attach.input); M <- chk$M; p <- chk$p; n_obs <- chk$n_obs
  fixed_rt <- .fixedAt(st, fixed)

  # Identity seed for the raw case.
  if (!has_seeds) {
    seeds <- .identitySeedsRaw(st, n_obs, fixed_rt)
    tangentX <- seeds$tangentX; tangentP <- seeds$tangentP
  }
  theta <- .resolveTheta(tangentX, tangentP)
  aligned <- .alignSeedsDual(st, tangentX, tangentP, NULL, NULL, n_obs, theta, fixed_rt)
  res <- .call_eval_ad(st, M, p, aligned$tangentX, aligned$tangentP, n_obs, theta)
  arr <- res$dy
  if (!has_seeds) {
    # Drop runtime-fixed columns from the canonical-basis output.
    dsyms <- setdiff(theta, fixed_rt)
    arr <- arr[, , dsyms, drop = FALSE]
  }
  .attachExtras(arr, n_obs, chk$extra_vars, chk$extra_params, "jac")
}

.hess_impl <- function(st, vars, params = numeric(0),
                                  tangentX = NULL, tangentP = NULL,
                                  hessianX = NULL, hessianP = NULL,
                                  attach.input = FALSE, fixed = NULL) {
  if (is.null(tangentX)) tangentX <- attr(vars,   "deriv")
  if (is.null(tangentP)) tangentP <- attr(params, "deriv")
  if (is.null(hessianX)) hessianX <- attr(vars,   "deriv2")
  if (is.null(hessianP)) hessianP <- attr(params, "deriv2")
  has_seeds <- !is.null(tangentX) || !is.null(tangentP) ||
               !is.null(hessianX) || !is.null(hessianP)
  chk <- .checkInputs(st, vars, params, attach.input); M <- chk$M; p <- chk$p; n_obs <- chk$n_obs
  fixed_rt <- .fixedAt(st, fixed)

  if (!has_seeds) {
    seeds <- .identitySeedsRaw(st, n_obs, fixed_rt)
    tangentX <- seeds$tangentX; tangentP <- seeds$tangentP
  }
  theta <- .resolveTheta(tangentX, tangentP, hessianX, hessianP)
  aligned <- .alignSeedsDual(st, tangentX, tangentP, hessianX, hessianP, n_obs, theta, fixed_rt)
  res <- .call_eval_ad2(st, M, p, aligned, n_obs, theta)
  arr <- res$d2y
  if (!has_seeds) {
    dsyms <- setdiff(theta, fixed_rt)
    arr <- arr[, , dsyms, dsyms, drop = FALSE]
  }
  .attachExtras(arr, n_obs, chk$extra_vars, chk$extra_params, "hess")
}

# Many argument sets in one .Call. Only the first-order dual path is batched,
# everything else loops, so the caller never has to branch. `sets` is a list of
# evaluate() argument lists: vars, params and optionally tangentX, tangentP, fixed.
.evaluateBatch_impl <- function(st, sets, cores = 1L, deriv2 = FALSE) {

  one <- function(a) do.call(.evaluate_impl,
                             c(list(st), a, list(deriv2 = deriv2)))
  sym <- .nativeSym(paste0(st$modelname,
                           if (deriv2) "_eval_ad2_batch" else "_eval_ad_batch"))
  if (is.null(sym) || length(sets) < 2L) return(lapply(sets, one))

  n_out <- length(st$outnames)
  # The requests of one batch are the same function at different numbers: the
  # `fixed` set and the theta basis are almost always shared. Carry the last
  # result forward instead of redoing the name algebra per request.
  fx_in <- NULL; fx_out <- character(0)
  th_key <- NULL; th_val <- NULL
  prep <- lapply(sets, function(a) {
    vars <- a$vars; params <- if (is.null(a$params)) numeric(0) else a$params
    tangentX <- if (is.null(a$tangentX)) attr(vars, "deriv")   else a$tangentX
    tangentP <- if (is.null(a$tangentP)) attr(params, "deriv") else a$tangentP
    hessianX <- if (deriv2) (if (is.null(a$hessianX)) attr(vars,   "deriv2") else a$hessianX)
    hessianP <- if (deriv2) (if (is.null(a$hessianP)) attr(params, "deriv2") else a$hessianP)
    has_seeds <- !is.null(tangentX) || !is.null(tangentP) ||
                 !is.null(hessianX) || !is.null(hessianP)
    att <- isTRUE(a$attach.input)
    chk <- .checkInputs(st, vars, params, att)
    if (is.null(a$fixed)) {
      fixed_rt <- .fixedAt(st, NULL)
    } else if (!is.null(fx_in) && identical(fx_in, a$fixed)) {
      fixed_rt <- fx_out
    } else {
      fixed_rt <- .fixedAt(st, a$fixed)
      fx_in <<- a$fixed; fx_out <<- fixed_rt
    }
    if (!has_seeds) {
      sd <- .identitySeedsRaw(st, chk$n_obs, fixed_rt)
      tangentX <- sd$tangentX; tangentP <- sd$tangentP
    }
    key <- list(dimnames(tangentX)[[3]], colnames(tangentP),
                if (deriv2) dimnames(hessianX)[[3]], if (deriv2) dimnames(hessianP)[[2]])
    if (!is.null(th_key) && identical(th_key, key)) {
      theta <- th_val
    } else {
      theta <- .resolveTheta(tangentX, tangentP, hessianX, hessianP)
      th_key <<- key; th_val <<- theta
    }
    al <- .alignSeedsDual(st, tangentX, tangentP, hessianX, hessianP, chk$n_obs, theta, fixed_rt)
    list(M = chk$M, p = chk$p, n_obs = chk$n_obs, theta = theta,
         aligned = al, has_seeds = has_seeds, fixed_rt = fixed_rt,
         extra_vars = chk$extra_vars, extra_params = chk$extra_params)
  })

  nullIf <- function(x)
    if (is.null(x) || !length(x)) NULL else if (is.double(x)) x else as.double(x)
  call_sets <- lapply(prep, function(q) {
    head <- list(as.double(q$M), as.double(q$p),
                 nullIf(q$aligned$tangentX), nullIf(q$aligned$tangentP))
    if (deriv2)
      head <- c(head, list(if (identical(q$aligned$has_dX2, 1L))
                             nullIf(q$aligned$hessianX) else NULL,
                           if (identical(q$aligned$has_dP2, 1L))
                             nullIf(q$aligned$hessianP) else NULL))
    c(head, list(as.integer(q$n_obs), as.integer(length(st$innames)),
                 as.integer(length(st$parameters)), as.integer(n_out),
                 as.integer(length(q$theta))))
  })

  raw <- .callSym(sym, call_sets, as.integer(cores))

  lapply(seq_along(prep), function(i) {
    q <- prep[[i]]; nt <- length(q$theta)
    keep <- if (q$has_seeds) q$theta else setdiff(q$theta, q$fixed_rt)
    y <- matrix(raw[[i]][[1L]], q$n_obs, n_out,
                dimnames = list(NULL, st$outnames))
    dy <- if (nt > 0)
      array(raw[[i]][[2L]][seq_len(n_out * nt * q$n_obs)],
            c(q$n_obs, n_out, nt), list(NULL, st$outnames, q$theta))
    else array(0, c(q$n_obs, n_out, 0L), list(NULL, st$outnames, NULL))
    if (!q$has_seeds) dy <- dy[, , keep, drop = FALSE]
    res <- list(y       = .attachExtras(y,  q$n_obs, q$extra_vars, q$extra_params, "fun"),
                tangent = .attachExtras(dy, q$n_obs, q$extra_vars, q$extra_params, "jac"))
    if (deriv2) {
      d2y <- if (nt > 0)
        array(raw[[i]][[3L]][seq_len(n_out * nt * nt * q$n_obs)],
              c(q$n_obs, n_out, nt, nt),
              list(NULL, st$outnames, q$theta, q$theta))
      else array(0, c(q$n_obs, n_out, 0L, 0L),
                 list(NULL, st$outnames, NULL, NULL))
      res$hessian <- if (q$has_seeds) d2y else d2y[, , keep, keep, drop = FALSE]
    }
    res
  })
}


.evaluate_impl <- function(st, vars, params = numeric(0),
                                          tangentX = NULL, tangentP = NULL,
                                          hessianX = NULL, hessianP = NULL,
                                          deriv2 = FALSE,
                                          attach.input = FALSE, fixed = NULL) {
  if (is.null(tangentX)) tangentX <- attr(vars,   "deriv")
  if (is.null(tangentP)) tangentP <- attr(params, "deriv")
  if (is.null(hessianX)) hessianX <- attr(vars,   "deriv2")
  if (is.null(hessianP)) hessianP <- attr(params, "deriv2")
  has_seeds <- !is.null(tangentX) || !is.null(tangentP) ||
               !is.null(hessianX) || !is.null(hessianP)
  chk <- .checkInputs(st, vars, params, attach.input); M <- chk$M; p <- chk$p; n_obs <- chk$n_obs
  fixed_rt <- .fixedAt(st, fixed)
  n_out <- length(st$outnames)

  if (!has_seeds) {
    seeds <- .identitySeedsRaw(st, n_obs, fixed_rt)
    tangentX <- seeds$tangentX; tangentP <- seeds$tangentP
  }
  theta <- .resolveTheta(tangentX, tangentP, hessianX, hessianP)
  aligned <- .alignSeedsDual(st, tangentX, tangentP, hessianX, hessianP, n_obs, theta, fixed_rt)
  if (deriv2) {
    res <- .call_eval_ad2(st, M, p, aligned, n_obs, theta)
    y <- res$y; dy <- res$dy; d2y <- res$d2y
    if (!has_seeds) {
      dsyms <- setdiff(theta, fixed_rt)
      dy  <- dy [, , dsyms, drop = FALSE]
      d2y <- d2y[, , dsyms, dsyms, drop = FALSE]
    }
  } else {
    res <- .call_eval_ad(st, M, p, aligned$tangentX, aligned$tangentP, n_obs, theta)
    y <- res$y; dy <- res$dy; d2y <- NULL
    if (!has_seeds) {
      dsyms <- setdiff(theta, fixed_rt)
      dy <- dy[, , dsyms, drop = FALSE]
    }
  }

  y   <- .attachExtras(y,   n_obs, chk$extra_vars, chk$extra_params, "fun")
  dy  <- .attachExtras(dy,  n_obs, chk$extra_vars, chk$extra_params, "jac")
  out <- list(y = y, tangent = dy)
  if (deriv2) out$hessian <- .attachExtras(d2y, n_obs, chk$extra_vars, chk$extra_params, "hess")
  out
}

# --- Convenient wrappers ---

.makeFunWrapper <- function(st, impl) {
  if (is.null(impl)) return(NULL)
  function(..., attach.input = FALSE, fixed = NULL) {
    args <- list(...); M <- if (length(st$innames)) do.call(cbind, args[st$innames]); p <- if (length(st$parameters)) do.call(c, args[st$parameters]) else numeric(0)
    if (attach.input) { extra <- setdiff(names(args), c(st$innames, st$parameters)); n_obs <- if (!is.null(M)) nrow(M) else 1L
    for (nm in extra) { v <- args[[nm]]; if (length(v) == n_obs) { M <- if (is.null(M)) matrix(v, ncol=1, dimnames=list(NULL,nm)) else cbind(M, setNames(data.frame(v), nm)) } else if (length(v) == 1) p <- c(p, setNames(v, nm)) else warning("Extra '", nm, "' ignored") } }
    impl(M, p, attach.input, fixed)
  }
}

.makeDerivWrapper <- function(st, impl, has_d2 = FALSE) {
  if (is.null(impl)) return(NULL)
  if (has_d2) {
    function(..., tangentX = NULL, tangentP = NULL, hessianX = NULL,
             hessianP = NULL, attach.input = FALSE, fixed = NULL) {
      args <- list(...); M <- if (length(st$innames)) do.call(cbind, args[st$innames]); p <- if (length(st$parameters)) do.call(c, args[st$parameters]) else numeric(0)
      if (attach.input) { extra <- setdiff(names(args), c(st$innames, st$parameters)); n_obs <- if (!is.null(M)) nrow(M) else 1L
      for (nm in extra) { v <- args[[nm]]; if (length(v) == n_obs) { M <- if (is.null(M)) matrix(v, ncol=1, dimnames=list(NULL,nm)) else cbind(M, setNames(data.frame(v), nm)) } else if (length(v) == 1) p <- c(p, setNames(v, nm)) else warning("Extra '", nm, "' ignored") } }
      impl(M, p, tangentX, tangentP, hessianX, hessianP, attach.input, fixed)
    }
  } else {
    function(..., tangentX = NULL, tangentP = NULL,
             attach.input = FALSE, fixed = NULL) {
      args <- list(...); M <- if (length(st$innames)) do.call(cbind, args[st$innames]); p <- if (length(st$parameters)) do.call(c, args[st$parameters]) else numeric(0)
      if (attach.input) { extra <- setdiff(names(args), c(st$innames, st$parameters)); n_obs <- if (!is.null(M)) nrow(M) else 1L
      for (nm in extra) { v <- args[[nm]]; if (length(v) == n_obs) { M <- if (is.null(M)) matrix(v, ncol=1, dimnames=list(NULL,nm)) else cbind(M, setNames(data.frame(v), nm)) } else if (length(v) == 1) p <- c(p, setNames(v, nm)) else warning("Extra '", nm, "' ignored") } }
      impl(M, p, tangentX, tangentP, attach.input, fixed)
    }
  }
}

.makeEvalWrapper <- function(st, impl) {
  if (is.null(impl)) return(NULL)
  function(..., tangentX = NULL, tangentP = NULL, hessianX = NULL,
           hessianP = NULL, deriv2 = FALSE, attach.input = FALSE, fixed = NULL) {
    args <- list(...); M <- if (length(st$innames)) do.call(cbind, args[st$innames]); p <- if (length(st$parameters)) do.call(c, args[st$parameters]) else numeric(0)
    if (attach.input) { extra <- setdiff(names(args), c(st$innames, st$parameters)); n_obs <- if (!is.null(M)) nrow(M) else 1L
    for (nm in extra) { v <- args[[nm]]; if (length(v) == n_obs) { M <- if (is.null(M)) matrix(v, ncol=1, dimnames=list(NULL,nm)) else cbind(M, setNames(data.frame(v), nm)) } else if (length(v) == 1) p <- c(p, setNames(v, nm)) else warning("Extra '", nm, "' ignored") } }
    impl(M, p, tangentX, tangentP, hessianX, hessianP, deriv2, attach.input, fixed)
  }
}
