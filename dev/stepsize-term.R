## =================================================================
## How much gradient does the step-size control carry?
##
## Forward sensitivities differentiate the solution at a frozen step
## sequence. The adaptive controller picks that sequence from theta, so
## the computed y_h(theta) additionally depends on theta through h, and
## that contribution is missing from every sensitivity cppDE returns
## today. This script measures whether it is worth recovering.
##
## The observable is the discretisation error e_h(theta) = y_h - y_ref.
## Its theta-derivative is exactly the missing term. A smooth e_h means
## the term is small and well defined; a sawtooth along the step-count
## switches means it interpolates between jumps and is an artefact.
## =================================================================
rm(list = ls(all.names = TRUE))

.workingDir <- file.path(tempdir(), "cppDE_stepsize_term")
dir.create(.workingDir, showWarnings = FALSE, recursive = TRUE)

library(cppDE)

RTOLS   <- c(1e-4, 1e-6, 1e-8, 1e-10)
RTOL_REF <- 1e-13
N_SWEEP <- 400

# --- Models ------------------------------------------------------------------
# Decay is non-stiff with a nearly constant step size, Robertson is stiff and
# spans many decades of h. The contrast is the point of running both.

eq_decay <- c(A = "-k1 * A", B = "k1 * A - k2 * B")
p_decay  <- c(A = 1, B = 0, k1 = 0.1, k2 = 0.2)
t_decay  <- c(0, 25, 50)

eq_rob <- c(y1 = "-k1*y1 + k2*y2*y3",
            y2 = "k1*y1 - k2*y2*y3 - k3*y2*y2",
            y3 = "k3*y2*y2")
p_rob  <- c(y1 = 1, y2 = 0, y3 = 0, k1 = 0.04, k2 = 1e4, k3 = 3e7)
t_rob  <- c(0, 10^seq(-5, 4, length.out = 40))

# --- Sweep -------------------------------------------------------------------
# One solve per theta value. `accepted` comes along because the step count is
# what makes the switching points visible.

sweep_one <- function(model, pars, times, pname, pvals, rtol, atol) {
  n <- length(pvals)
  y <- NULL
  steps <- integer(n)
  for (i in seq_len(n)) {
    p <- pars; p[[pname]] <- pvals[i]
    r <- solveODE(model, times, p, abstol = atol, reltol = rtol)
    v <- r$variable[nrow(r$variable), ]
    if (is.null(y)) y <- matrix(NA_real_, n, length(v), dimnames = list(NULL, names(v)))
    y[i, ] <- v
    steps[i] <- r$diagnostics$accepted
  }
  list(y = y, steps = steps)
}

# Forward sensitivity at the sweep centre, for the scale the error is judged on.
sens_at <- function(model, pars, times, pname, rtol, atol) {
  r <- solveODE(model, times, pars, abstol = atol, reltol = rtol)
  idx <- which(dimnames(r$sens1)$sens == pname)
  r$sens1[dim(r$sens1)[1], , idx]
}

# --- Analysis of one model ---------------------------------------------------
# Returns the sweep, the reference, and the relative size of the missing term:
# max |d e_h / d theta| measured across the sweep, against |dy/dtheta|.

analyse <- function(label, eqns, pars, times, pname, rel_span, atol, method = "bdf") {
  cat("\n=========================================================\n")
  cat(label, " sweeping ", pname, " over +/-", rel_span, " relative\n", sep = "")
  cat("=========================================================\n")

  mod_v <- cppODE(eqns, method = method, deriv = FALSE,
                  modelname = paste0("sst_v_", label), outdir = .workingDir)
  mod_s <- cppODE(eqns, method = method, deriv = TRUE,
                  modelname = paste0("sst_s_", label), outdir = .workingDir)

  p0    <- pars[[pname]]
  pvals <- p0 * (1 + seq(-rel_span, rel_span, length.out = N_SWEEP))
  dtheta <- diff(pvals)[1]

  ref <- sweep_one(mod_v, pars, times, pname, pvals, RTOL_REF, atol * 1e-3)

  res <- list()
  for (rt in RTOLS) {
    sw <- sweep_one(mod_v, pars, times, pname, pvals, rt, atol)
    e  <- sw$y - ref$y                       # discretisation error over theta
    de <- apply(e, 2, function(col) diff(col) / dtheta)   # its theta-derivative
    S  <- sens_at(mod_s, pars, times, pname, rt, atol)

    # Relative size of the missing term, per state, worst over the sweep.
    rel <- apply(abs(de), 2, max) / pmax(abs(S), .Machine$double.eps)

    res[[format(rt)]] <- list(rtol = rt, e = e, de = de, S = S, rel = rel,
                              steps = sw$steps, y = sw$y)
    cat(sprintf("  rtol %-8.0e  |e_h| max %10.3e   steps %4d..%4d   ",
                rt, max(abs(e)), min(sw$steps), max(sw$steps)))
    cat("rel. missing term ", paste(sprintf("%.2e", rel), collapse = " "), "\n", sep = "")
  }

  invisible(list(pvals = pvals, ref = ref, res = res, pname = pname))
}

# --- Run ---------------------------------------------------------------------

dec <- analyse("decay", eq_decay, p_decay, t_decay, "k1",
               rel_span = 1e-3, atol = 1e-10)

rob <- analyse("robertson", eq_rob, p_rob, t_rob, "k1",
               rel_span = 1e-3, atol = 1e-10)

# --- Plots -------------------------------------------------------------------
# Top: the discretisation error over theta. Smooth means the missing term is a
# real derivative; steps mean it is an artefact of where the controller switches.
# Bottom: the accepted step count, which locates those switches.

plot_case <- function(a, state, main) {
  op <- par(mfrow = c(2, 1), mar = c(4, 4.5, 2.5, 1))
  on.exit(par(op))
  cols <- seq_along(a$res)
  ymax <- max(sapply(a$res, function(r) max(abs(r$e[, state]))))
  plot(a$pvals, a$res[[1]]$e[, state], type = "n", ylim = c(-ymax, ymax),
       xlab = a$pname, ylab = paste0("e_h [", state, "]"), main = main)
  abline(h = 0, col = "grey70")
  for (i in cols) lines(a$pvals, a$res[[i]]$e[, state], col = i)
  legend("topright", legend = sprintf("rtol %.0e", sapply(a$res, `[[`, "rtol")),
         col = cols, lty = 1, bty = "n", cex = 0.8)
  plot(a$pvals, a$res[[1]]$steps, type = "s", xlab = a$pname,
       ylab = "accepted steps", col = 1,
       ylim = range(sapply(a$res, function(r) range(r$steps))))
  for (i in cols) lines(a$pvals, a$res[[i]]$steps, type = "s", col = i)
  invisible(NULL)
}

plot_case(dec, "A",  "decay: discretisation error over k1")
plot_case(rob, "y1", "robertson: discretisation error over k1")
