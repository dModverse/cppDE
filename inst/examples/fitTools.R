## =================================================================
## Fitting helpers for the example_fit*.R scripts and the README
## =================================================================

## Each function works on `residuals(pars)`, which returns the scaled residuals
## `r` and their Jacobian `J` with columns named by parameter, or NULL where
## the model cannot be solved. A name starting with "logit10", "log10",
## "logit" or "log" marks a parameter on that scale; logit10 is
## log10(x / (1 - x)).

library(ggplot2)

toNatural <- function(name, value) {
  if (startsWith(name, "logit10")) return(1 / (1 + 10^(-value)))
  if (startsWith(name, "log10"))   return(10^value)
  if (startsWith(name, "logit"))   return(plogis(value))
  if (startsWith(name, "log"))     return(exp(value))
  value
}
natural <- function(p) mapply(toNatural, names(p), p)

## Gaussian prior as additional residuals (p - mean) / sd, flagged in `prior`.
withPrior <- function(residuals, mean, sd) {
  function(pars) {
    rj <- residuals(pars)
    if (is.null(rj)) return(NULL)
    nm <- names(mean)
    Jp <- matrix(0, length(nm), ncol(rj$J), dimnames = list(nm, colnames(rj$J)))
    Jp[cbind(nm, nm)] <- 1 / sd
    list(r = c(rj$r, (pars[nm] - mean) / sd), J = rbind(rj$J, Jp),
         prior = c(logical(length(rj$r)), rep(TRUE, length(nm))))
  }
}

## Wide model output to long form, one row per time and component.
asLong <- function(m, tt) {
  data.frame(time = rep(tt, ncol(m)),
             component = rep(colnames(m), each = nrow(m)),
             m = as.vector(m))
}

chi2Of <- function(residuals, pars) {
  rj <- residuals(pars)
  if (is.null(rj)) Inf else sum(rj$r^2)
}

## The part of the objective a prior contributes, zero without one.
priorOf <- function(residuals, pars) {
  rj <- residuals(pars)
  if (is.null(rj)) return(Inf)
  if (is.null(rj$prior)) 0 else sum(rj$r[rj$prior]^2)
}

## -----------------------------------------------------------------
## Gauss-Newton on a trust region
## -----------------------------------------------------------------

## The subproblem is solved exactly, More-Sorensen: the interior Newton step
## where it fits the radius, otherwise the boundary solution from Newton
## iteration on the secular equation. J'J is positive semidefinite, so the
## indefinite case cannot arise here.
subproblem <- function(B, g, delta) {
  chol2 <- function(lambda) tryCatch(chol(B + lambda * diag(nrow(B))),
                                     error = function(e) NULL)
  solveR <- function(R) drop(backsolve(R, forwardsolve(t(R), -g)))
  R <- chol2(0)
  if (!is.null(R)) {
    p <- solveR(R)
    if (sqrt(sum(p^2)) <= delta) return(p)
  }
  lambda <- sqrt(sum(g^2)) / delta
  for (k in 1:50) {
    R <- chol2(lambda)
    if (is.null(R)) { lambda <- 2 * lambda + 1e-10; next }
    p  <- solveR(R)
    np <- sqrt(sum(p^2))
    if (abs(np - delta) <= 0.01 * delta) break
    q      <- forwardsolve(t(R), p)
    lambda <- max(0, lambda + (np / sqrt(sum(q^2)))^2 * (np - delta) / delta)
  }
  p
}

## Hessian approximation J'J. Residuals and Jacobian of an accepted trial are
## reused in the next iteration. Parameters in `fixed` are held constant.
gaussNewton <- function(residuals, pars, fixed = NULL, maxit = 400,
                        tol = 1e-10, delta = 1, deltaMax = 10) {
  free <- setdiff(names(pars), fixed)
  rj   <- residuals(pars)
  if (is.null(rj)) stop("the solver failed at the start")
  value <- sum(rj$r^2)
  it    <- 0L
  repeat {
    it <- it + 1L
    if (it > maxit) break
    J <- rj$J[, free, drop = FALSE]
    g <- drop(2 * crossprod(J, rj$r))
    B <- 2 * crossprod(J)
    p <- subproblem(B, g, delta)

    trial       <- pars
    trial[free] <- pars[free] + p
    rjTrial     <- residuals(trial)
    vTrial      <- if (is.null(rjTrial)) Inf else sum(rjTrial$r^2)
    pred        <- -(sum(g * p) + 0.5 * sum(p * (B %*% p)))
    rho         <- if (pred > 0) (value - vTrial) / pred else -1

    if (rho < 0.25) delta <- delta / 4
    else if (rho > 0.75 && sqrt(sum(p^2)) > 0.99 * delta)
      delta <- min(2 * delta, deltaMax)

    if (rho > 0.1) {
      change <- value - vTrial
      pars   <- trial
      value  <- vTrial
      rj     <- rjTrial
      if (change < tol * (1 + abs(value))) break
    }
    if (delta < 1e-12) break
  }
  prior <- if (is.null(rj$prior)) 0 else sum(rj$r[rj$prior]^2)
  list(par = pars, value = value, prior = prior, iterations = it,
       converged = it <= maxit)
}

## -----------------------------------------------------------------
## Multistart
## -----------------------------------------------------------------

## One fit per row of `starts`, in parallel. Status "global": within `tol` of
## the best value; "local": converged elsewhere; "iteration limit": not
## converged.
multistart <- function(residuals, starts, cores = 1L, tol = 0.1) {
  fits <- parallel::mclapply(seq_len(nrow(starts)), function(i)
    tryCatch(gaussNewton(residuals, starts[i, ]), error = function(e) NULL),
    mc.cores = cores)
  value     <- vapply(fits, function(f) if (is.list(f)) f$value else NA_real_,
                      numeric(1))
  converged <- vapply(fits, function(f) is.list(f) && f$converged, logical(1))
  best      <- which.min(value)
  status    <- ifelse(value - value[best] < tol, "global",
                      ifelse(converged, "local", "iteration limit"))
  status[is.na(value)] <- "failed"
  list(fits = fits, value = value, status = status, best = fits[[best]])
}

## Sorted objective values minus the best, on a pseudo-log axis.
plotWaterfall <- function(ms) {
  o  <- order(ms$value, na.last = NA)
  df <- data.frame(rank = seq_along(o), delta = ms$value[o] - ms$value[o[1]],
                   status = factor(ms$status[o],
                                   c("global", "local", "iteration limit")))
  ## The axis spans at least two decades, so that fits within the noise of the
  ## optimum sit on one level.
  ggplot(df, aes(rank, delta)) +
    geom_step(colour = "grey70") +
    geom_point(aes(colour = status), size = 2) +
    expand_limits(y = c(0, 100)) +
    scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
                       breaks = c(0, 10^(0:8))) +
    scale_colour_manual(values = c(global = "#1b7837", local = "#b2182b",
                                   `iteration limit` = "grey55"),
                        drop = FALSE) +
    labs(x = "start, sorted by final objective", y = "objective - best",
         colour = NULL, title = "Waterfall",
         subtitle = sprintf("%d of %d starts reach the best optimum",
                            sum(ms$status == "global"), length(ms$value)))
}

## -----------------------------------------------------------------
## Noise model
## -----------------------------------------------------------------

## Moving-average noise of order q = length(w) on one series,
## u_t = eta_t + sum_j w_j eta_(t-j), with white innovations eta of scale
## exp(logSigma). Correlations beyond q samples are zero. The likelihood is the
## conditional sum of squares: innovations before the first sample are zero,
## which turns the whitening into a recursive filter. `u` are data-unit
## residuals, `Ju` their Jacobian; the log term enters as the pseudo-residual
## sqrt(n (C + 2 logSigma)), which needs C + 2 logSigma > 0.
## Box, Jenkins, Reinsel, Ljung, Time Series Analysis, 5th ed., Wiley (2015).
maResiduals <- function(u, Ju, logSigma, w, C = 100) {
  if (C + 2 * logSigma <= 0) return(NULL)
  n      <- length(u)
  sg     <- exp(logSigma)
  whiten <- function(x) unclass(stats::filter(x, -w, method = "recursive"))
  eta    <- as.vector(whiten(u))
  if (!all(is.finite(eta)) || max(abs(eta)) > 1e6 * max(abs(u))) return(NULL)
  dW <- vapply(seq_along(w), function(j)
    -as.vector(whiten(c(numeric(j), eta[seq_len(n - j)]))), numeric(n))
  a  <- sqrt(n * (C + 2 * logSigma))
  list(r      = c(eta / sg, a),
       J      = rbind(matrix(whiten(Ju), n, dimnames = dimnames(Ju)) / sg, 0),
       dSigma = c(-eta / sg, n / a),
       dW     = rbind(dW / sg, 0))
}

## The same objective, sum(eta^2) / sigma^2 + n (C + 2 logSigma), with its
## gradient in `u`, `logSigma` and `w` instead of a Jacobian. `du` is the seed a
## reverse sweep contracts with du/dtheta. `eta` are the scaled innovations.
maObjective <- function(u, logSigma, w, C = 100) {
  if (C + 2 * logSigma <= 0) return(NULL)
  n   <- length(u)
  sg2 <- exp(2 * logSigma)
  eta <- as.vector(stats::filter(u, -w, method = "recursive"))
  if (!all(is.finite(eta)) || max(abs(eta)) > 1e6 * max(abs(u))) return(NULL)
  v <- rev(as.vector(stats::filter(rev(2 * eta / sg2), -w,
                                   method = "recursive")))
  list(value  = sum(eta^2) / sg2 + n * (C + 2 * logSigma),
       du     = v,
       dSigma = 2 * n - 2 * sum(eta^2) / sg2,
       dW     = vapply(seq_along(w), function(j)
                  -sum(v[(j + 1):n] * eta[seq_len(n - j)]), numeric(1)),
       eta    = eta / sqrt(sg2))
}

## -----------------------------------------------------------------
## Profile likelihood
## -----------------------------------------------------------------

## One parameter is stepped away from the optimum and the others are
## re-optimised from the previous point. The step doubles, up to `stepCap`, while the objective rises
## by less than 0.1 per step; a step that rises by more than 0.5 is retried
## at half the size. Each side walks until the data part has risen past the
## threshold, `maxSteps` points are taken or the walk is `range` away from the
## optimum; a prior never ends a walk.
profile <- function(residuals, name, par, step, threshold, maxSteps, stepCap,
                    range) {
  ref  <- chi2Of(residuals, par)
  refP <- priorOf(residuals, par)
  walk <- function(direction) {
    p    <- par
    last <- ref
    h    <- step
    out  <- matrix(numeric(0), 0, 3,
                   dimnames = list(NULL, c("value", "chi2", "prior")))
    while (nrow(out) < maxSteps) {
      room <- range - abs(p[[name]] - par[[name]])
      if (room <= 1e-12 * range) break
      trial         <- p
      trial[[name]] <- p[[name]] + direction * min(h, room)
      fit <- tryCatch(gaussNewton(residuals, trial, fixed = name),
                      error = function(e) NULL)
      if (is.null(fit)) break
      rise <- fit$value - last
      if (rise > 0.5 && h > step / 64) {
        h <- h / 2
        next
      }
      p    <- fit$par
      last <- fit$value
      out  <- rbind(out, c(value = p[[name]], chi2 = fit$value,
                           prior = fit$prior))
      if ((fit$value - fit$prior) - (ref - refP) > threshold) break
      if (rise < 0.1) h <- min(2 * h, stepCap)
    }
    out
  }
  left <- walk(-1)
  out  <- rbind(left[rev(seq_len(nrow(left))), , drop = FALSE],
                c(value = par[[name]], chi2 = ref, prior = refP),
                walk(1))
  structure(out, centre = nrow(left) + 1L)
}

## The parameters in `which`, over forked workers. A walk starts with
## `stepSE` standard errors from the curvature, at most `stepCap` on the
## fitting scale, and reaches at most `range` from the optimum.
profiles <- function(residuals, par, which = names(par),
                     threshold = qchisq(0.95, df = 1), range = 5,
                     stepSE = 0.1, stepCap = range / 50, maxSteps = 100,
                     cores = 1L) {
  se   <- sqrt(diag(solve(crossprod(residuals(par)$J))))
  step <- pmin(stepSE * se, stepCap)
  out  <- parallel::mclapply(which, function(nm)
            profile(residuals, nm, par, step[[nm]], threshold, maxSteps,
                    stepCap, range),
          mc.cores = cores)
  names(out) <- which
  structure(out, reference = chi2Of(residuals, par),
            referencePrior = priorOf(residuals, par), threshold = threshold)
}

## Where the data part crosses the threshold, interpolated linearly. A side
## that does not reach the threshold gives -Inf or Inf.
profileIntervals <- function(prof) {
  ref <- attr(prof, "reference") - attr(prof, "referencePrior")
  thr <- attr(prof, "threshold")
  t(vapply(names(prof), function(nm) {
    pr   <- prof[[nm]]
    x    <- pr[, "value"]
    d    <- pr[, "chi2"] - pr[, "prior"] - ref
    edge <- function(idx, open) {
      above <- which(d[idx] > thr)
      if (!length(above)) return(open)
      j <- idx[above[1]]
      i <- idx[above[1] - 1]
      x[i] + (thr - d[i]) * (x[j] - x[i]) / (d[j] - d[i])
    }
    c0 <- attr(pr, "centre")
    c(lower = toNatural(nm, edge(c0:1, -Inf)),
      upper = toNatural(nm, edge(c0:nrow(pr), Inf)))
  }, numeric(2)))
}

## Data and prior contributions of each profile, plotted separately. Log and
## logit parameters are plotted as log10 of the natural value.
plotProfiles <- function(prof, truth = NULL,
                         labels = setNames(names(prof), names(prof)),
                         markLabel = "true value") {
  ref    <- attr(prof, "reference")
  refP   <- attr(prof, "referencePrior")
  thr    <- attr(prof, "threshold")
  axis   <- function(nm, v) {
    x <- toNatural(nm, v)
    if (startsWith(nm, "log")) log10(x) else x
  }
  title  <- function(nm) {
    if (startsWith(nm, "log")) paste("log10", labels[[nm]]) else labels[[nm]]
  }
  levels <- vapply(names(prof), title, character(1), USE.NAMES = FALSE)
  curves <- do.call(rbind, lapply(names(prof), function(nm) {
    pr    <- prof[[nm]]
    total <- pr[, "chi2"] - ref
    prior <- pr[, "prior"] - refP
    data.frame(parameter = title(nm),
               value = axis(nm, pr[, "value"]),
               delta = c(total - prior, prior),
               part  = rep(c("data", "prior"), each = nrow(pr)))
  }))
  if (all(curves$delta[curves$part == "prior"] == 0))
    curves <- curves[curves$part == "data", ]
  curves$parameter <- factor(curves$parameter, levels = levels)
  curves$part      <- factor(curves$part, c("data", "prior"))

  ## the optimum each walk started from
  best <- data.frame(parameter = factor(levels, levels = levels),
                     value = vapply(names(prof), function(nm) {
                       pr <- prof[[nm]]
                       axis(nm, pr[attr(pr, "centre"), "value"])
                     }, numeric(1)))

  subtitle <- "dashed grey: 95 percent threshold, dot: best fit"
  g <- ggplot(curves, aes(value, delta, colour = part, linetype = part)) +
    geom_hline(yintercept = thr, linetype = 2, colour = "grey40") +
    geom_line(linewidth = 0.7) +
    geom_point(data = best, aes(value, 0), inherit.aes = FALSE,
               colour = "#b2182b", size = 2) +
    scale_colour_manual(values = c(data = "black", prior = "#2166ac")) +
    scale_linetype_manual(values = c(data = "solid", prior = "longdash")) +
    facet_wrap(~parameter, scales = "free_x", ncol = 3) +
    coord_cartesian(ylim = c(max(min(0, curves$delta), -thr), 1.1 * thr)) +
    labs(x = NULL, y = expression(Delta * chi^2), colour = NULL,
         linetype = NULL, title = "Profile likelihood", subtitle = subtitle)
  if (is.null(truth)) return(g)

  marks <- data.frame(parameter = factor(vapply(names(truth), title,
                                                character(1)),
                                         levels = levels),
                      value = mapply(axis, names(truth), truth))
  g + geom_vline(data = marks, aes(xintercept = value), inherit.aes = FALSE,
                 colour = "#1b7837") +
    labs(subtitle = paste0(subtitle, ", green: ", markLabel))
}

## -----------------------------------------------------------------
## Residual diagnostics
## -----------------------------------------------------------------

## Per group of residuals: the lag-one autocorrelation and a Ljung-Box test
## over the first `lag` coefficients (Ljung, Box, Biometrika 65, 297, 1978).
residualCheck <- function(res, lag = 10) {
  t(vapply(res, function(r) {
    lb <- Box.test(r, lag = lag, type = "Ljung-Box")
    c(n = length(r), chi2_per_point = mean(r^2),
      lag1 = acf(r, lag.max = 1, plot = FALSE)$acf[2],
      LjungBox = unname(lb$statistic), p = lb$p.value)
  }, numeric(5)))
}

plotResidualACF <- function(res, lag.max = 40) {
  df <- do.call(rbind, lapply(names(res), function(nm) {
    a <- acf(res[[nm]], lag.max = lag.max, plot = FALSE)$acf[-1]
    data.frame(group = nm, lag = seq_along(a), acf = a,
               band = 1.96 / sqrt(length(res[[nm]])))
  }))
  df$group <- factor(df$group, names(res))
  ggplot(df, aes(lag, acf)) +
    geom_ribbon(aes(ymin = -band, ymax = band), fill = "grey85") +
    geom_hline(yintercept = 0, colour = "grey50") +
    geom_segment(aes(xend = lag, yend = 0), colour = "#2166ac") +
    facet_wrap(~group, ncol = 4) +
    labs(x = "lag [samples]", y = "autocorrelation",
         title = "Autocorrelation of the residuals",
         subtitle = "grey band: 95 percent range for white noise")
}
