## =================================================================
## Does the adjoint's own per-step number predict where the gradient
## error sits?
##
## Stufe 9 of dev/adjoint-plan.md wants the step size chosen from lambda
## as well as from the state, under err = max(err_state, err_lambda). The
## weight it uses has to be a refinement indicator: a per-step number
## whose ranking says where refining buys the most gradient accuracy.
##
## Two candidates fall out of a sweep for free, and they are not the same
## quantity:
##
##   eta_k = lambda_{k+1}' e_k     e_k the step's embedded error estimate.
##                                 The dual-weighted residual. O(tol) per
##                                 step, and its sum estimates the error
##                                 in J.
##   wdt_k = dJ/dh_k               the transport derivative, lambda' f.
##                                 O(J), not O(error): it says how much J
##                                 moves when the step covers more time,
##                                 not how much accuracy the step costs.
##
## The plan named wdt as the indicator. This script is the check, and it
## is written so that a wrong guess shows up rather than hides.
##
## Measured, per model and over rtol:
##
##   1. sum |eta_k| against errJ, the error in the seeded functional
##      J = sum W_oi x_i(t_o) itself. This is what the dual-weighted
##      residual estimates, and the check is direct: the ratio should sit
##      near one and stay there over the decades.
##   2. sum |eta_k| against errG, the error in the gradient. That is the
##      quantity Stufe 9 actually wants to control, and eta does not
##      estimate it, which would take a second-order adjoint. What the
##      measurement can settle is whether the two fall together, so that
##      controlling the one controls the other.
##   3. the same for sum |wdt_k| h_k, the plan's original guess. If wdt is
##      the transport derivative, this ratio grows by a decade per decade
##      of tolerance, because the numerator does not fall and the
##      denominator does.
##   4. how concentrated each indicator is. One that is flat over the
##      trajectory has nothing to steer with, however well it sums.
##
## Reference values come from a solve at rtol_ref. Errors are measured per
## component with a floor: a gradient whose components span decades, as
## Robertson's does, is not summarised by a norm divided by its largest
## entry.
##
## The models are the ones the plan asks for: a stiff one and an event
## one, not only the decay toy.
## =================================================================
rm(list = ls(all.names = TRUE))

.workingDir <- file.path(tempdir(), "cppDE_lambda_indicator")
dir.create(.workingDir, showWarnings = FALSE, recursive = TRUE)

library(cppDE)

RTOLS    <- c(1e-4, 1e-6, 1e-8, 1e-10)
RTOL_REF <- 1e-13

# Robertson's y2 sits near 1e-5 for the whole run, so an absolute tolerance
# scaled to y1 and y3 is meaningless there and one scaled to y2 is unreachable
# for the others. It gets a fixed atol and a reference that stays inside what
# the method can deliver; the classic treatment.

# --- Models ------------------------------------------------------------------
# Decay is the control: smooth, non-stiff, nearly constant h. Robertson is
# stiff, so lambda spans decades. The root model is where the plan expects the
# first-order remainder to misbehave, a crossing making the map discontinuous
# in the state.

models <- list(
  decay = list(
    eq     = c(A = "-k1 * A", B = "k1 * A - k2 * B"),
    pars   = c(A = 1, B = 0, k1 = 0.1, k2 = 0.2),
    times  = seq(0, 50, by = 2.5),
    events = NULL, method = "tsit5"),

  robertson = list(
    eq     = c(y1 = "-k1*y1 + k2*y2*y3",
               y2 = "k1*y1 - k2*y2*y3 - k3*y2*y2",
               y3 = "k3*y2*y2"),
    pars   = c(y1 = 1, y2 = 0, y3 = 0, k1 = 0.04, k2 = 1e4, k3 = 3e7),
    times  = c(0, 10^seq(-5, 4, length.out = 30)),
    events = NULL, method = "bdf",
    atol = 1e-10, rtol_ref = 1e-12),

  rooted = list(
    eq     = c(A = "-k1 * A + k2 * B",
               B = "k1 * A - k2 * B - k3 * B * B"),
    pars   = c(A = 1.2, B = 0.4, k1 = 0.7, k2 = 0.35, k3 = 1.1),
    times  = seq(0, 12, by = 0.5),
    events = data.frame(var = "B", time = NA, value = 1.5, root = "A - 0.6",
                        method = "multiply", stringsAsFactors = FALSE),
    method = "tsit5")
)

# --- Build ------------------------------------------------------------------
# A value model beside the reverse one: its solve fixes the output row count,
# which a root event makes larger than length(times), and the seed's first
# dimension has to match it.

build <- function(name, m) {
  cat("compiling ", name, " ...\n", sep = "")
  m$value <- cppODE(m$eq, events = m$events, modelname = paste0("li_v_", name),
                    deriv = FALSE, method = m$method, outdir = .workingDir)
  m$rev   <- cppODE(m$eq, events = m$events, modelname = paste0("li_r_", name),
                    sweep = "reverse", method = m$method, outdir = .workingDir)

  v <- solveODE(m$value, m$times, m$pars,
                abstol = if (is.null(m$atol)) 1e-10 else m$atol, reltol = 1e-10)
  n_out <- length(v$time)
  n_x   <- ncol(v$variable)
  ## J = sum over observed times and states, weighted 1, 1/2, 1/3, ... over
  ## the states. A weight that is constant across states is annihilated by a
  ## linear conservation law. Robertson conserves y1+y2+y3, so an all-ones
  ## seed makes J the constant 1 per time point and its gradient exactly zero.
  ## Nothing about the study needs a residual-shaped objective; what matters is
  ## that lambda is not trivial.
  m$seed  <- array(rep(1 / seq_len(n_x), each = n_out),
                   dim = c(n_out, n_x, 1L))
  m$n_out <- n_out
  m
}

gradient <- function(m, rtol, grid = FALSE) {
  solveODE(m$rev, m$times, m$pars,
           abstol = if (is.null(m$atol)) rtol else m$atol, reltol = rtol,
           seed = m$seed, adjointGrid = grid)
}

# The functional the seed defines: J = sum over observed times and states of
# W * x. The reverse pass returns its gradient, so this is the value its
# adjoint belongs to.
functional <- function(r, m) sum(r$variable * m$seed[, , 1])

# Worst relative deviation over the components that carry something. A
# component three decades below the largest contributes nothing to any use of
# the gradient and would otherwise dominate a per-component ratio.
rel_worst <- function(g, g_ref) {
  keep <- abs(g_ref) > 1e-6 * max(abs(g_ref))
  if (!any(keep)) return(NA_real_)
  max(abs(g[keep] - g_ref[keep]) / abs(g_ref[keep]))
}

# Share of the total carried by the worst tenth of steps. 1 means one step
# carries everything, 0.1 means the indicator is flat.
concentration <- function(v) {
  s <- sum(v)
  if (!is.finite(s) || s <= 0) return(NA_real_)
  top <- sort(v, decreasing = TRUE)[seq_len(max(1L, length(v) %/% 10))]
  sum(top) / s
}

study <- function(name, m) {
  m <- build(name, m)

  ref   <- gradient(m, if (is.null(m$rtol_ref)) RTOL_REF else m$rtol_ref)
  g_ref <- as.numeric(ref$adjoint[, 1])
  J_ref <- functional(ref, m)

  out <- do.call(rbind, lapply(RTOLS, function(rt) {
    r <- gradient(m, rt, grid = TRUE)
    g <- as.numeric(r$adjoint[, 1])
    G <- r$adjointGrid
    eta <- abs(G$eta[, 1])
    wdt <- abs(G$wdt[, 1]) * G$h
    data.frame(rtol     = rt,
               steps    = length(G$h),
               errJ     = abs(functional(r, m) - J_ref) / max(abs(J_ref), 1),
               errG     = rel_worst(g, g_ref),
               sum_eta  = sum(eta),
               sum_wdth = sum(wdt),
               conc_eta = concentration(eta),
               conc_wdt = concentration(wdt))
  }))

  ## The ratios are the point: an indicator that tracks a quantity keeps a flat
  ## ratio to it over the decades, whatever its prefactor.
  out$eta_over_J <- out$sum_eta  / out$errJ
  out$eta_over_G <- out$sum_eta  / out$errG
  out$wdt_over_G <- out$sum_wdth / out$errG

  cat("\n=== ", name, " ===\n", sep = "")
  print(out, digits = 3, row.names = FALSE)
  spread <- function(v) {
    v <- v[is.finite(v) & v > 0]
    if (length(v) < 2) return(NA_real_)
    signif(max(v) / min(v), 3)
  }
  cat("\n  eta / errJ spreads by ", spread(out$eta_over_J),
      "   (this is what eta estimates)\n", sep = "")
  cat("  eta / errG spreads by ", spread(out$eta_over_G),
      "   (this is what Stufe 9 wants to control)\n", sep = "")
  cat("  (wdt*h) / errG spreads by ", spread(out$wdt_over_G), "\n", sep = "")
  cat("  worst tenth of steps carries ",
      signif(100 * mean(out$conc_eta), 3), "% of sum |eta|\n", sep = "")
  invisible(out)
}

for (nm in names(models))
  tryCatch(study(nm, models[[nm]]),
           error = function(e)
             cat("\n", nm, " failed: ", conditionMessage(e), "\n", sep = ""))
