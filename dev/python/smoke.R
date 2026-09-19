## Compile checks of the graph-AD code generators, one PSOCK worker per case.
##
##   Rscript dev/python/smoke.R <phase> [workers]
##
## The generators are read from inst/python (CPPDE_PY_DIR); cppDE's R code and
## headers come from the installed package. A `[legacy]` twin reads the frozen
## generators of dev/python/legacy (e55c599) instead.

args <- commandArgs(trailingOnly = TRUE)
phase <- if (length(args)) args[1] else "2"
workers <- if (length(args) > 1) as.integer(args[2]) else 4L

root <- normalizePath(local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  file.path(dirname(f[1]), "..", "..")
}), "/")
## SMOKE_LEGACY=1 runs every case on the frozen generators; SMOKE_ONLY selects
## cases by regular expression; SMOKE_PROGRESS names a progress file.
Sys.setenv(CPPDE_PY_DIR = file.path(root, "inst", "python"))
legacy_dir <- file.path(root, "dev", "python", "legacy")

`%||%` <- function(a, b) if (is.null(a)) b else a

cases <- list()
pairs <- list()
## A case returns list(err = named errors, val = raw values); `vs_legacy` adds
## a twin run with the installed generators and compares the raw values.
case <- function(name, fn, limit = 1e-8, vs_legacy = FALSE, pair_limit = 1e-10) {
  cases[[name]] <<- list(fn = fn, legacy = FALSE, limit = limit)
  if (vs_legacy) {
    twin <- paste(name, "[legacy]")
    cases[[twin]] <<- list(fn = fn, legacy = TRUE, limit = limit)
    pairs[[name]] <<- list(a = name, b = twin, limit = pair_limit)
  }
}

## ---------------------------------------------------------------------------
## Models and helpers (shipped to the workers)
## ---------------------------------------------------------------------------

S1 <- list(
  rhs = c(A = "-k1*A*B + k2*C - piecewise(kd*A, time - ts > 0, 0.1*kd*A)",
          B = "-k1*A*B + k2*C + u*kin",
          C = "k1*A*B - k2*C - exp(-time/tau)*C"),
  pars = c(A = 1, B = 0.8, C = 0.1, k1 = 0.9, k2 = 0.3, kd = 0.5, ts = 1.5,
           kin = 0.4, tau = 2),
  forcings = list(u = data.frame(time = c(0, 2, 5), value = c(0.1, 0.3, 0.05))),
  times = seq(0, 5, length.out = 11),
  ## the switch time is not an event: its sensitivity misses the jump
  no_fd = "ts")

chain <- function(n) {
  rhs <- character(n); names(rhs) <- paste0("x", seq_len(n))
  rhs[1] <- "-k1*x1 + k0/(1 + x1)"
  for (i in 2:n)
    rhs[i] <- sprintf("k%d*x%d - k%d*x%d^2", i, i - 1, i + 1, i)
  rhs
}
S2 <- list(rhs = chain(40), times = seq(0, 4, length.out = 9))
S2$pars <- c(setNames(rep(0.1, 40), names(S2$rhs)),
             setNames(seq(0.5, 1.5, length.out = 42), paste0("k", 0:41)))

tol <- list(abstol = 1e-12, reltol = 1e-12, roottol = 1e-12)

contract <- function(tangent, W)
  vapply(seq_len(dim(W)[3]),
         function(k) apply(tangent * as.vector(W[, , k]), 3, sum),
         numeric(dim(tangent)[3]))

hess_forward <- function(res, W) {
  ns <- dim(res$hessian)[3]
  outer(seq_len(ns), seq_len(ns),
        Vectorize(function(a, b) sum(as.vector(W) * res$hessian[, , a, b])))
}

relerr <- function(a, b) max(abs(a - b)) / max(1, max(abs(b)))

cotangent_for <- function(n_t, n_x) { set.seed(1); array(rnorm(n_t * n_x), c(n_t, n_x, 1)) }

solve_with <- function(m, M, ...)
  do.call(cppDE::solveODE, c(list(m, M$times, M$pars), M$tol %||% tol,
                             if (!is.null(M$forcings)) list(forcings = M$forcings),
                             list(...)))

## reverse against the forward tangent (first order)
rev_vs_fwd <- function(M, method, outdir, tag, sparse = NULL) {
  mf <- cppDE::cppODE(M$rhs, events = M$events, forcings = names(M$forcings), method = method,
                      sparse = sparse, deriv = TRUE, fixed = M$fixed, outdir = outdir,
                      modelname = paste0(tag, "_f"))
  mr <- cppDE::cppODE(M$rhs, events = M$events, forcings = names(M$forcings), method = method,
                      sparse = sparse, derivMode = "reverse", fixed = M$fixed,
                      outdir = outdir, modelname = paste0(tag, "_r"))
  ff <- solve_with(mf, M)
  W <- cotangent_for(nrow(ff$variable), ncol(ff$variable))
  fr <- solve_with(mr, M, cotangent = W)
  ref <- contract(ff$tangent, W)[, 1]
  list(err = c(gradient = relerr(fr$cotangent[names(ref), 1], ref)),
       val = c(fr$cotangent[names(ref), 1]))
}

## forward-reverse against forward-forward (second order)
rev2_vs_fwd2 <- function(M, method, outdir, tag, sparse = NULL) {
  mf <- cppDE::cppODE(M$rhs, events = M$events, forcings = names(M$forcings), method = method,
                      sparse = sparse, derivMode = "forward-forward", fixed = M$fixed,
                      outdir = outdir, modelname = paste0(tag, "_ff"))
  mr <- cppDE::cppODE(M$rhs, events = M$events, forcings = names(M$forcings), method = method,
                      sparse = sparse, derivMode = "forward-reverse", fixed = M$fixed,
                      outdir = outdir, modelname = paste0(tag, "_fr"))
  ff <- solve_with(mf, M)
  W <- cotangent_for(nrow(ff$variable), ncol(ff$variable))
  W[!(ff$time %in% M$times), , ] <- 0
  fr <- solve_with(mr, M, cotangent = W)
  ref <- contract(ff$tangent, W)[, 1]
  hr <- unname(fr$curvature[names(ref), , 1])
  list(err = c(gradient = relerr(fr$cotangent[names(ref), 1], ref),
               hessian = relerr(hr, hess_forward(ff, W))),
       val = c(fr$cotangent[names(ref), 1], hr))
}


## Jacobian strategy of the generator for the rest of this worker.
set_jac <- function(strategy) {
  if (nzchar(Sys.getenv("CPPDE_PY_DIR")))
    reticulate::py_run_string(sprintf(
      "import os; os.environ['CPPDE_JAC'] = '%s'", strategy))
}

## Central differences of the trajectory (and of the tangent) in each parameter.
fd_sens <- function(m, M, pars, what = "variable", h = 1e-5) {
  base <- solve_with(m, M)[[what]]
  out <- array(0, c(dim(base), length(pars)))
  for (k in seq_along(pars)) {
    nm <- pars[k]
    hk <- h * max(1, abs(M$pars[[nm]]))
    Mp <- M; Mp$pars[[nm]] <- Mp$pars[[nm]] + hk
    Mm <- M; Mm$pars[[nm]] <- Mm$pars[[nm]] - hk
    d <- (solve_with(m, Mp)[[what]] - solve_with(m, Mm)[[what]]) / (2 * hk)
    if (length(dim(base)) == 2) out[, , k] <- d else out[, , , k] <- d
  }
  dimnames(out) <- c(dimnames(base), list(pars))
  out
}

## Forward sensitivities against finite differences; raw sens as values.
fwd_vs_fd <- function(M, method, outdir, tag, sparse = NULL, jac = "") {
  set_jac(jac)
  t0 <- proc.time()[["elapsed"]]
  m <- cppDE::cppODE(M$rhs, events = M$events, forcings = names(M$forcings), method = method,
                     sparse = sparse, deriv = TRUE, fixed = M$fixed, outdir = outdir,
                     modelname = tag, compile = FALSE)
  t1 <- proc.time()[["elapsed"]]
  cppDE::compile(m)
  t2 <- proc.time()[["elapsed"]]
  r <- solve_with(m, M)
  keep <- setdiff(dimnames(r$tangent)[[3]], M$no_fd)
  fd <- fd_sens(m, M, keep, h = M$fd_h %||% 1e-5)
  list(err = c(sens1 = relerr(r$tangent[, , keep], fd)), val = c(r$variable, r$tangent),
       note = sprintf("codegen %.1f s, compile %.1f s, %.0f kB source",
                      t1 - t0, t2 - t1, file.size(attr(m, "srcfile")) / 1e3))
}

## Second-order sensitivities against differences of the tangent.
fwd2_vs_fd <- function(M, method, outdir, tag, sparse = NULL, jac = "") {
  set_jac(jac)
  m <- cppDE::cppODE(M$rhs, events = M$events, forcings = names(M$forcings), method = method,
                     sparse = sparse, derivMode = "forward-forward",
                     outdir = outdir, modelname = tag)
  r <- solve_with(m, M)
  m1 <- cppDE::cppODE(M$rhs, events = M$events, forcings = names(M$forcings), method = method,
                      sparse = sparse, deriv = TRUE, outdir = outdir,
                      modelname = paste0(tag, "_1"))
  pars <- dimnames(r$hessian)[[4]]
  fd2 <- fd_sens(m1, M, pars, what = "tangent")
  fd2 <- fd2[, , dimnames(r$hessian)[[3]], , drop = FALSE]
  list(err = c(sens2 = relerr(r$hessian, fd2)), val = c(r$tangent, r$hessian))
}

## Small S2 for second order: 12 chain states.
S2s <- list(rhs = chain(12), times = seq(0, 3, length.out = 7))
S2s$pars <- c(setNames(rep(0.1, 12), names(S2s$rhs)),
              setNames(seq(0.5, 1.5, length.out = 14), paste0("k", 0:13)))
S1f <- S1
S1f$forcings <- NULL
S1f$rhs["B"] <- "-k1*A*B + k2*C + kin"


S3 <- list(
  ## The fixed event resets A, whose right-hand side reads the clock.
  rhs = c(A = "-k1*A + k2*B*time", B = "k1*A - k2*B^2"),
  events = data.frame(var = c("A", "B"), value = c("A + dose*B", "B*fac - 0.2"),
                      time = c("t_dose*2", NA),
                      root = c(NA, "A^2 - thr*B - 0.1*time"),
                      method = c("add", "replace"), stringsAsFactors = FALSE),
  pars = c(A = 1, B = 0.5, k1 = 0.8, k2 = 0.3, dose = 0.4, t_dose = 0.7,
           fac = 0.9, thr = 0.3),
  times = seq(0, 4, length.out = 9))

## Forward-reverse gradient of the seeded functional against differences.
rev_grad_fd <- function(M, method, outdir, tag) {
  mr <- cppDE::cppODE(M$rhs, events = M$events, method = method,
                      derivMode = "forward-reverse", outdir = outdir,
                      modelname = tag)
  m0 <- cppDE::cppODE(M$rhs, events = M$events, method = method, deriv = FALSE,
                      outdir = outdir, modelname = paste0(tag, "_v"))
  base <- solve_with(m0, M)
  grid <- base$time %in% M$times
  W <- cotangent_for(nrow(base$variable), ncol(base$variable))
  W[!grid, , ] <- 0
  fr <- solve_with(mr, M, cotangent = W)
  fun <- function(Mp) {
    r <- solve_with(m0, Mp)
    sum(r$variable[r$time %in% M$times, , drop = FALSE] * W[grid, , 1])
  }
  nm <- names(M$pars)
  fd <- vapply(nm, function(k) {
    h <- 1e-6 * max(1, abs(M$pars[[k]]))
    Mp <- M; Mp$pars[[k]] <- Mp$pars[[k]] + h
    Mm <- M; Mm$pars[[k]] <- Mm$pars[[k]] - h
    (fun(Mp) - fun(Mm)) / (2 * h)
  }, 0)
  list(err = c(gradient = relerr(fr$cotangent[nm, 1], fd)),
       val = c(fr$cotangent[nm, 1], fr$curvature[, , 1]))
}


## S5: LLG macrospins with dense dipolar coupling (dev/python/harness.py).
llg_dipole <- function(ns, drive = TRUE) {
  k <- 0:(ns - 1)
  pos <- cbind(cos(2 * pi * k / ns), sin(2 * pi * k / ns), 0.35 * k / ns)
  nm <- function(a, s) sprintf("m%s%d", c("x", "y", "z")[a], s - 1)
  tens <- function(s, t) {
    if (s == t) return(diag(-0.2, 3))
    r <- pos[t, ] - pos[s, ]
    d <- sqrt(sum(r^2))
    u <- r / d
    (3 * outer(u, u) - diag(3)) * 0.01 / d^3
  }
  ext <- c(if (drive) "Hx0 + h1*cos(phi)" else "Hx0", "Hy0", "Hz0")
  rhs <- character(0)
  for (s in seq_len(ns)) {
    H <- character(3)
    for (a in 1:3) {
      terms <- character(0)
      for (t in seq_len(ns)) {
        D <- tens(s, t)
        for (b in 1:3) if (D[a, b] != 0)
          terms <- c(terms, sprintf("%s*%s", format(D[a, b], digits = 17), nm(b, t)))
      }
      H[a] <- sprintf("(%s + Ms*(%s))", ext[a], paste(terms, collapse = " + "))
    }
    mx <- nm(1, s); my <- nm(2, s); mz <- nm(3, s)
    rhs[mx] <- sprintf("-gp*(%s*(-alpha*(%s^2 + %s^2)) + %s*(alpha*%s*%s - %s) + %s*(%s + alpha*%s*%s))",
                       H[1], my, mz, H[2], mx, my, mz, H[3], my, mx, mz)
    rhs[my] <- sprintf("-gp*(%s*(%s + alpha*%s*%s) + %s*(-alpha*(%s^2 + %s^2)) + %s*(alpha*%s*%s - %s))",
                       H[1], mz, mx, my, H[2], mx, mz, H[3], my, mz, mx)
    rhs[mz] <- sprintf("-gp*(%s*(alpha*%s*%s - %s) + %s*(%s + alpha*%s*%s) + %s*(-alpha*(%s^2 + %s^2)))",
                       H[1], mx, mz, my, H[2], mx, my, mz, H[3], mx, my)
  }
  if (drive) rhs["phi"] <- "omega"
  th <- 0.3 + 0.1 * k
  ph <- 0.2 * k
  init <- setNames(c(rbind(sin(th) * cos(ph), sin(th) * sin(ph), cos(th))),
                   names(rhs)[seq_len(3 * ns)])
  list(rhs = rhs, times = seq(0, 3, length.out = 7),
       pars = c(init, if (drive) c(phi = 0),
                gp = 1, alpha = 0.1, Ms = 1, Hx0 = 0.3, Hy0 = 0.1, Hz0 = 1,
                if (drive) c(h1 = 0.2, omega = 2)),
       no_fd = names(rhs))
}
S5 <- llg_dipole(20)
## Second order in the parameters only, over a short span.
S5s <- llg_dipole(6)
S5s$fixed <- names(S5s$rhs)
S5s$times <- seq(0, 1, length.out = 5)

## Row threshold of the linear map for the rest of this worker ("" = default).
set_lin <- function(value) {
  if (nzchar(Sys.getenv("CPPDE_PY_DIR")))
    reticulate::py_run_string(sprintf(
      "import os; os.environ['CPPDE_LINEAR'] = '%s'", value))
}

## S6: Brusselator on an nx x ny grid (dev/python/harness.py); the initial
## values are fixed, so the sensitivities are those of the four parameters.
brusselator2d <- function(nx, ny) {
  idx <- expand.grid(j = seq_len(ny) - 1, i = seq_len(nx) - 1)[, c("i", "j")]
  rhs <- character(0)
  for (s in c("u", "v")) for (r in seq_len(nrow(idx))) {
    i <- idx$i[r]; j <- idx$j[r]
    nb <- list(c(i - 1, j), c(i + 1, j), c(i, j - 1), c(i, j + 1))
    nb <- Filter(function(p) p[1] >= 0 && p[1] < nx && p[2] >= 0 && p[2] < ny, nb)
    cc <- sprintf("%s%d_%d", s, i, j)
    lap <- paste0(paste(vapply(nb, function(p) sprintf("%s%d_%d", s, p[1], p[2]), ""),
                        collapse = " + "), sprintf(" - %d*%s", length(nb), cc))
    u <- sprintf("u%d_%d", i, j); v <- sprintf("v%d_%d", i, j)
    rhs[cc] <- if (s == "u")
      sprintf("Du*(%s) + a - (b + 1)*%s + %s^2*%s", lap, u, u, v)
    else sprintf("Dv*(%s) + b*%s - %s^2*%s", lap, u, u, v)
  }
  x <- idx$i / max(1, nx - 1); y <- idx$j / max(1, ny - 1)
  init <- c(setNames(1 + 0.2 * sin(2 * pi * x) * cos(pi * y), sprintf("u%d_%d", idx$i, idx$j)),
            setNames(3 - 0.1 * cos(pi * x * y), sprintf("v%d_%d", idx$i, idx$j)))
  list(rhs = rhs, times = seq(0, 1, length.out = 5),
       pars = c(init[names(rhs)], Du = 0.002 * nx * nx, Dv = 0.001 * nx * nx, a = 1, b = 3),
       fixed = names(rhs), no_fd = character(0))
}
## Loop threshold of the generator for the rest of this worker ("" = default).
set_vec <- function(value) {
  if (nzchar(Sys.getenv("CPPDE_PY_DIR")))
    reticulate::py_run_string(sprintf(
      "import os; os.environ['CPPDE_VECTORISE'] = '%s'", value))
}

## CVODE against the native backend on the requested times: trajectory and
## sensitivities forward, the seeded gradient in reverse.
cv_vs_native <- function(M, outdir, tag, sparse = NULL, reverse = FALSE) {
  mn <- cppDE::cppODE(M$rhs, events = M$events, forcings = names(M$forcings),
                      fixed = M$fixed, deriv = TRUE, sparse = sparse, method = "bdf",
                      outdir = outdir, modelname = paste0(tag, "_n"))
  mc <- cppDE::cvode(M$rhs, events = M$events, forcings = names(M$forcings),
                     fixed = M$fixed, deriv = !reverse, sparse = sparse,
                     derivMode = if (reverse) "reverse" else "forward",
                     outdir = outdir, modelname = paste0(tag, "_c"))
  rn <- solve_with(mn, M)
  on <- match(M$times, rn$time)
  if (reverse) {
    W <- cotangent_for(nrow(rn$variable), ncol(rn$variable))
    W[-on, , ] <- 0
    rc <- solve_with(mc, M, cotangent = W[on, , , drop = FALSE])
    ref <- contract(rn$tangent, W)[, 1]
    return(list(err = c(gradient = relerr(rc$cotangent[names(ref), 1], ref)),
                val = rc$cotangent[names(ref), 1]))
  }
  rc <- solve_with(mc, M)
  oc <- match(M$times, rc$time)
  list(err = c(variable = relerr(rc$variable[oc, ], rn$variable[on, ]),
               sens1 = relerr(rc$tangent[oc, , ], rn$tangent[on, , ])),
       val = c(rc$variable[oc, ], rc$tangent[oc, , ]))
}

## cppFUN S4: jac/hess against differences, vjp against jac, its curvature
## against hess.
fun_s4 <- function(outdir, tag) {
  eq <- c(y1 = "piecewise(a*x^2, x > c, b*sqrt(x))",
          y2 = "abs(x - a)*pow(b, 2.5)",
          y3 = "exp(-k*x) + a*b*c")
  f <- cppDE::cppFUN(eq, variables = "x", parameters = c("a", "b", "c", "k"),
                     fixed = "c", deriv2 = TRUE, derivMode = c("forward", "reverse"),
                     compile = TRUE, outdir = outdir, modelname = tag,
                     convenient = FALSE)
  X <- matrix(c(0.3, 1.7), ncol = 1, dimnames = list(NULL, "x"))
  P <- c(a = 0.9, b = 1.3, c = 1.1, k = 0.4)
  J <- f$jac(X, P)
  H <- f$hess(X, P)
  syms <- dimnames(J)[[3]]
  val <- function(x, p) f$func(matrix(x, ncol = 1, dimnames = list(NULL, "x")), p)
  jfd <- array(0, dim(J))
  for (k in seq_along(syms)) {
    h <- 1e-6
    if (syms[k] == "x") {
      jfd[, , k] <- (val(X[, 1] + h, P) - val(X[, 1] - h, P)) / (2 * h)
    } else {
      Pp <- P; Pp[syms[k]] <- Pp[syms[k]] + h
      Pm <- P; Pm[syms[k]] <- Pm[syms[k]] - h
      jfd[, , k] <- (val(X[, 1], Pp) - val(X[, 1], Pm)) / (2 * h)
    }
  }
  jac_at <- function(x, p) f$jac(matrix(x, ncol = 1, dimnames = list(NULL, "x")), p)
  hfd <- array(0, dim(H))
  for (k in seq_along(syms)) {
    h <- 1e-5
    if (syms[k] == "x") {
      d <- (jac_at(X[, 1] + h, P) - jac_at(X[, 1] - h, P)) / (2 * h)
    } else {
      Pp <- P; Pp[syms[k]] <- Pp[syms[k]] + h
      Pm <- P; Pm[syms[k]] <- Pm[syms[k]] - h
      d <- (jac_at(X[, 1], Pp) - jac_at(X[, 1], Pm)) / (2 * h)
    }
    hfd[, , , k] <- d
  }
  set.seed(2)
  W <- matrix(rnorm(6), 2, 3)
  r <- f$vjp(X, P, W)
  jx <- vapply(1:2, function(o) sum(W[o, ] * J[o, , "x"]), 0)
  jp <- vapply(c("a", "b", "k"), function(nm) sum(W * J[, , nm]), 0)
  VX <- array(rnorm(2), c(2, 1, 1))
  VP <- matrix(rnorm(4), 4, 1, dimnames = list(names(P), NULL))
  VP["c", 1] <- 0
  r2 <- f$vjp(X, P, W, tangentX = VX, tangentP = VP)
  full <- array(0, c(2, 4, 4), list(NULL, c("x", "a", "b", "k"), c("x", "a", "b", "k")))
  for (o in 1:2) full[o, , ] <- W[o, 1] * H[o, 1, , ] + W[o, 2] * H[o, 2, , ] + W[o, 3] * H[o, 3, , ]
  ref_dx <- vapply(1:2, function(o) sum(full[o, "x", ] * c(VX[o, 1, 1], VP[c("a", "b", "k"), 1])), 0)
  ref_dp <- vapply(c("a", "b", "k"), function(nm)
    sum(vapply(1:2, function(o) sum(full[o, nm, ] * c(VX[o, 1, 1], VP[c("a", "b", "k"), 1])), 0)), 0)
  list(err = c(jac_fd = relerr(J, jfd), hess_fd = relerr(H, hfd),
               vjp_jac = relerr(c(r$cotangentX[, 1, 1], r$cotangentP[c("a", "b", "k"), 1]),
                                c(jx, jp)),
               curvature_hess = relerr(c(r2$curvatureX[, 1, 1, 1],
                                         r2$curvatureP[c(1, 2, 4), 1, 1]),
                                       c(ref_dx, ref_dp))),
       val = c(J, H, r$cotangentX, r$cotangentP, r2$curvatureX, r2$curvatureP))
}

## ---------------------------------------------------------------------------
## Cases per phase
## ---------------------------------------------------------------------------

if (phase == "2") {
  case("S1 reverse bdf", function(d) rev_vs_fwd(S1, "bdf", d, "s1bdf"),
       vs_legacy = TRUE)
  ## The time switch is not an event: ff and fr grids differ by O(tol) there.
  case("S1 forward-reverse rb4", function(d) rev2_vs_fwd2(S1, "rb4", d, "s1rb4"),
       limit = 1e-6, vs_legacy = TRUE)
  case("S2 reverse bdf sparse", function(d) rev_vs_fwd(S2, "bdf", d, "s2", sparse = TRUE),
       vs_legacy = TRUE)
}

if (phase == "3") {
  case("S1 forward bdf", function(d) fwd_vs_fd(S1, "bdf", d, "s1f"),
       limit = 1e-6, vs_legacy = TRUE, pair_limit = 1e-8)
  case("S1 tsit5", function(d) fwd_vs_fd(S1, "tsit5", d, "s1t"),
       limit = 1e-5, vs_legacy = TRUE, pair_limit = 1e-8)
  case("S2 ff rb4 sparse entries", function(d)
         fwd2_vs_fd(S2s, "rb4", d, "s2e", sparse = TRUE, jac = "entries"),
       limit = 1e-5, vs_legacy = TRUE, pair_limit = 1e-8)
  case("S2 ff rb4 sparse colour", function(d)
         fwd2_vs_fd(S2s, "rb4", d, "s2c", sparse = TRUE, jac = "colour"),
       limit = 1e-5)
  case("S2 forward bdf dense colour", function(d)
         fwd_vs_fd(S2, "bdf", d, "s2dc", sparse = FALSE, jac = "colour"),
       limit = 1e-6)
  case("S2 forward bdf dense entries", function(d)
         fwd_vs_fd(S2, "bdf", d, "s2de", sparse = FALSE, jac = "entries"),
       limit = 1e-6)
  pairs[["S2 colour vs entries (rb4 ff)"]] <- list(
    a = "S2 ff rb4 sparse colour", b = "S2 ff rb4 sparse entries", limit = 1e-10)
  pairs[["S2 colour vs entries (bdf)"]] <- list(
    a = "S2 forward bdf dense colour", b = "S2 forward bdf dense entries", limit = 1e-10)
}

if (phase == "4") {
  for (meth in c("bdf", "rb4", "tsit5")) local({
    mm <- meth
    case(paste("S3 ff vs fr", mm), function(d)
           rev2_vs_fwd2(S3, mm, d, paste0("s3", mm)),
         limit = 1e-5, vs_legacy = TRUE, pair_limit = 1e-8)
    case(paste("S3 fr gradient vs fd", mm), function(d)
           rev_grad_fd(S3, mm, d, paste0("s3g", mm)),
         limit = 1e-5, vs_legacy = TRUE, pair_limit = 1e-8)
  })
}

if (phase == "5") {
  case("S4 cppFUN forward+reverse deriv2", function(d) fun_s4(d, "s4"),
       limit = 1e-5, vs_legacy = TRUE, pair_limit = 1e-10)
}

if (phase == "6") {
  for (lin in c("E", "plain")) local({
    lv <- if (lin == "E") "" else "0"
    tg <- tolower(lin)
    case(paste("S5 forward rb4 dense", lin), function(d) {
      set_lin(lv)
      fwd_vs_fd(S5, "rb4", d, paste0("s5f", tg), sparse = FALSE)
    }, limit = 1e-5)
    case(paste("S5 reverse bdf sparse", lin), function(d) {
      set_lin(lv)
      rev_vs_fwd(S5, "bdf", d, paste0("s5r", tg), sparse = TRUE)
    }, limit = 1e-6)
    case(paste("S5s ff vs fr rb4", lin), function(d) {
      set_lin(lv)
      rev2_vs_fwd2(S5s, "rb4", d, paste0("s5s", tg))
    }, limit = 1e-6)
  })
  for (nm in c("S5 forward rb4 dense", "S5 reverse bdf sparse", "S5s ff vs fr rb4"))
    pairs[[paste(nm, "E vs plain")]] <- list(a = paste(nm, "E"), b = paste(nm, "plain"),
                                            limit = 1e-8)
}

if (phase == "7") {
  S6 <- brusselator2d(70, 70)
  ## Second order on a small grid.
  S6s <- brusselator2d(8, 8)
  S6s$tol <- list(abstol = 1e-9, reltol = 1e-9)
  ## 481 states with dense coupling over a short span.
  S6l <- llg_dipole(160)
  S6l$fixed <- names(S6l$rhs)
  S6l$times <- c(0, 0.25, 0.5)
  S6l$tol <- list(abstol = 1e-10, reltol = 1e-10)
  S6l$fd_h <- 1e-6
  case("S6 forward bdf sparse 2D", function(d)
         fwd_vs_fd(S6, "bdf", d, "s6f", sparse = TRUE), limit = 1e-5)
  case("S6 reverse rb4 sparse 2D", function(d)
         rev_vs_fwd(S6, "rb4", d, "s6r", sparse = TRUE), limit = 1e-6)
  case("S6 LLG 481 forward bdf dense", function(d)
         fwd_vs_fd(S6l, "bdf", d, "s6l", sparse = FALSE), limit = 1e-4)
  for (v in c("loops", "plain")) local({
    vv <- if (v == "loops") "" else "0"
    case(paste("S6s forward-reverse rb4", v), function(d) {
      set_vec(vv)
      rev2_vs_fwd2(S6s, "rb4", d, paste0("s6s", v), sparse = TRUE)
    }, limit = 1e-6)
  })
  pairs[["S6s loops vs plain"]] <- list(a = "S6s forward-reverse rb4 loops",
                                        b = "S6s forward-reverse rb4 plain",
                                        limit = 1e-8)
}

if (phase == "8") {
  S6c <- brusselator2d(12, 12)
  case("CV S3 events forward", function(d) cv_vs_native(S3, d, "cv3"),
       limit = 1e-6, vs_legacy = TRUE, pair_limit = 1e-8)
  case("CV S1 forcing forward", function(d) cv_vs_native(S1f, d, "cv1"),
       limit = 1e-6, vs_legacy = TRUE, pair_limit = 1e-8)
  case("CV S2 sparse forward", function(d) cv_vs_native(S2, d, "cv2", sparse = TRUE),
       limit = 1e-6, vs_legacy = TRUE, pair_limit = 1e-8)
  case("CV S2 reverse", function(d) cv_vs_native(S2, d, "cv2r", reverse = TRUE),
       limit = 1e-5, vs_legacy = TRUE, pair_limit = 1e-8)
  case("CV S5s map forward", function(d) cv_vs_native(S5s, d, "cv5"),
       limit = 1e-6, vs_legacy = TRUE, pair_limit = 1e-8)
  case("CV S6c loops reverse sparse", function(d)
         cv_vs_native(S6c, d, "cv6", sparse = TRUE, reverse = TRUE), limit = 1e-5)
}

## ---------------------------------------------------------------------------
## Run
## ---------------------------------------------------------------------------

only <- Sys.getenv("SMOKE_ONLY")
if (nzchar(only)) {
  cases <- cases[grepl(only, names(cases))]
  pairs <- pairs[vapply(pairs, function(p) all(c(p$a, p$b) %in% names(cases)), NA)]
}
if (nzchar(Sys.getenv("SMOKE_LEGACY")))
  for (nm in names(cases)) cases[[nm]]$legacy <- TRUE
progress <- Sys.getenv("SMOKE_PROGRESS")

outdir <- tempfile("cppde_smoke_")
dir.create(outdir)
cat(sprintf("phase %s: %d cases, %d workers, sources in %s\n",
            phase, length(cases), workers, outdir))
cl <- parallel::makePSOCKcluster(min(workers, length(cases)))
on.exit(parallel::stopCluster(cl), add = TRUE)
parallel::clusterExport(cl, setdiff(ls(), c("cl", "cases")))
t0 <- proc.time()[["elapsed"]]
res <- parallel::clusterApplyLB(cl, names(cases), function(nm, cases) {
  Sys.setenv(CPPDE_PY_DIR = if (cases[[nm]]$legacy) legacy_dir
                            else file.path(root, "inst", "python"))
  suppressPackageStartupMessages(library(cppDE))
  ## One worker runs many cases: every case imports its own generators.
  cache <- cppDE:::.cppde_py_cache
  rm(list = ls(cache), envir = cache)
  if (reticulate::py_available())
    reticulate::py_run_string(paste(
      "import os, sys",
      "for _v in ('CPPDE_JAC', 'CPPDE_LINEAR', 'CPPDE_VECTORISE'): os.environ.pop(_v, None)",
      "for _m in [m for m in sys.modules if m.startswith(('codegen_', 'cppde_', 'cppsympy', 'derivSymb'))]:",
      "    del sys.modules[_m]", sep = "\n"))
  d <- file.path(outdir, gsub("[^A-Za-z0-9]", "_", nm))
  dir.create(d, showWarnings = FALSE)
  t1 <- proc.time()[["elapsed"]]
  r <- tryCatch(cases[[nm]]$fn(d), error = function(e) conditionMessage(e))
  s <- proc.time()[["elapsed"]] - t1
  if (nzchar(progress))
    cat(sprintf("%-40s %s (%.0f s)\n", nm,
                if (is.character(r)) paste("ERROR", r)
                else paste(sprintf("%s %.1e", names(r$err), r$err), collapse = ", "), s),
        file = progress, append = TRUE)
  list(result = r, seconds = s)
}, cases = cases)
names(res) <- names(cases)

ok <- TRUE
for (nm in names(res)) {
  r <- res[[nm]]$result
  if (is.character(r)) {
    ok <- FALSE
    cat(sprintf("%-40s ERROR %s\n", nm, r))
    next
  }
  pass <- all(r$err <= cases[[nm]]$limit)
  ok <- ok && pass
  cat(sprintf("%-40s %s  %s  (%.0f s)%s\n", nm, if (pass) "ok  " else "FAIL",
              paste(sprintf("%s %.1e", names(r$err), r$err), collapse = ", "),
              res[[nm]]$seconds, if (is.null(r$note)) "" else paste0("  ", r$note)))
}
for (nm in names(pairs)) {
  p <- pairs[[nm]]
  a <- res[[p$a]]$result
  b <- res[[p$b]]$result
  if (is.character(a) || is.character(b)) next
  e <- relerr(a$val, b$val)
  pass <- e <= p$limit
  ok <- ok && pass
  cat(sprintf("%-40s %s  values differ by %.1e\n", nm, if (pass) "ok  " else "FAIL", e))
}
cat(sprintf("total %.0f s\n", proc.time()[["elapsed"]] - t0))
quit(status = if (ok) 0 else 1)
