## =================================================================
## Fitting measured magnetisation precession with cppDE
## =================================================================

## [Co/Pt]x5 multilayer with perpendicular anisotropy, field tilted by about
## 30 degrees from the film normal. An 80 fs pump pulse lowers the anisotropy,
## and the magnetisation precesses about the shifted equilibrium while the film
## cools. The probe measures the anomalous Hall conductivity, proportional to mz.
##
## Data: Zhang et al., Phys. Rev. Applied 24, 044079 (2025), Fig. 3;
## doi:10.5281/zenodo.17333041, CC BY 4.0, see inst/COPYRIGHTS.

rm(list = ls(all.names = TRUE))
setwd(tempdir())

library(cppDE)
library(ggplot2)

## Helpers: Gauss-Newton with More-Sorensen trust region, multistart, profile
## likelihood, moving-average noise model.
source(system.file("examples", "fitTools.R", package = "cppDE"))

theme_set(theme_bw(base_size = 11))
cores <- if (.Platform$OS.type == "windows") 1L else parallel::detectCores()

## -----------------------------------------------------------------
## 1. Data
## -----------------------------------------------------------------

raw    <- read.csv(system.file("extdata", "zhang2025_precession.csv",
                               package = "cppDE"))
fields <- sort(unique(raw$field_kOe))

## The first picosecond (ultrafast demagnetisation) is not described by the
## model and is left out of the fit.
obs      <- lapply(split(raw, raw$field_kOe), function(tr) tr[tr$time_ps >= 1, ])
measured <- do.call(rbind, obs)
kOe      <- as_labeller(function(x) paste(x, "kOe"))

## y axis linear near zero and logarithmic above.
signalAxis <- scale_y_continuous(
  trans  = scales::pseudo_log_trans(sigma = 0.005, base = 10),
  breaks = c(-0.01, 0, 0.01, 0.03, 0.1, 0.3, 1))

print(
  ggplot(raw, aes(time_ps, signal)) +
    geom_point(size = 0.5) +
    signalAxis +
    facet_wrap(~field_kOe, ncol = 4, labeller = kOe) +
    labs(x = "pump-probe delay [ps]", y = "normalised Hall signal",
         title = "Measured traces")
)

## -----------------------------------------------------------------
## 2. Model
## -----------------------------------------------------------------

muB <- 8.794e-5                      # mu_B / hbar in rad / (ps mT)

## LLG, m' = -g (m x h + alpha m x (m x h)), g = g-factor mu_B / hbar /
## (1 + alpha^2), field h = H (sin phi, 0, cos phi) + Ha (1 - K kappa) mz z.
## States xi = (m - m0) / K and kappa, with the equilibrium
## m0 = (sin th, 0, cos th). Since h(m0) = lam m0, the right-hand side is a
## polynomial in K without constant term. Positive parameters enter as log10.
a   <- "(10^log10Alpha)"
g   <- sprintf("(%g*10^log10G/(1 + %s^2))", muB, a)
xm  <- "(xi_x*sin(th) + xi_z*cos(th))"                  # xi . m0
xx  <- "(xi_x^2 + xi_y^2 + xi_z^2)"                     # xi . xi
eta <- "(10^log10Ha*(xi_z - kap*(cos(th) + K*xi_z)))"   # (h - lam m0) / K, z
a1  <- sprintf("(lam*%s + %s*cos(th))", xm, eta)
a2  <- sprintf("(%s*xi_z)", eta)
b1  <- sprintf("(2*%s)", xm)
b2  <- xx

## m x h / K
cx <- sprintf("(lam*xi_y*cos(th) + K*%s*xi_y)", eta)
cy <- sprintf("(lam*(xi_z*sin(th) - xi_x*cos(th)) - %s*sin(th) - K*%s*xi_x)",
              eta, eta)
cz <- "(-lam*xi_y*sin(th))"

## (m (m.h) - h (m.m)) / K
dx <- sprintf(paste0("(lam*xi_x + sin(th)*%s - lam*sin(th)*%s",
                     " + K*(xi_x*%s + sin(th)*%s - lam*sin(th)*%s)",
                     " + K^2*xi_x*%s)"),
              a1, b1, a1, a2, b2, a2)
dy <- sprintf("(lam*xi_y + K*xi_y*%s + K^2*xi_y*%s)", a1, a2)
dz <- sprintf(paste0("(lam*xi_z + cos(th)*%s - %s - lam*cos(th)*%s",
                     " + K*(xi_z*%s + cos(th)*%s - %s*%s - lam*cos(th)*%s)",
                     " + K^2*(xi_z*%s - %s*%s))"),
              a1, eta, b1, a1, a2, eta, b1, b2, a2, eta, b2)

llg <- c(
  xi_x = sprintf("-%s*(%s + %s*%s)", g, cx, a, dx),
  xi_y = sprintf("-%s*(%s + %s*%s)", g, cy, a, dy),
  xi_z = sprintf("-%s*(%s + %s*%s)", g, cz, a, dz),
  kap  = "-kap*10^(-log10TauK)"
)

## Initial states xi = 0 and kappa = 1 are not estimated. H and phi enter only
## through th and lam.
model <- cppODE(llg, fixed = c("xi_x", "xi_y", "xi_z", "kap"),
                modelname = "copt_precession", method = "bdf", deriv = TRUE)

## -----------------------------------------------------------------
## 3. Likelihood
## -----------------------------------------------------------------

## Equilibrium before the pump: H sin(th - phi) + Ha sin(th) cos(th) = 0 and
## lam = H cos(th - phi) + Ha cos(th)^2. Derivatives of th and lam with respect
## to Ha and phi from the implicit function theorem.
equilibrium <- function(Ha, H, phi) {
  torque <- function(th) H * sin(th - phi) + Ha * sin(th) * cos(th)
  th     <- uniroot(torque, c(0, phi), tol = 1e-12 * phi)$root
  sc     <- sin(th) * cos(th)
  den    <- H * cos(th - phi) + Ha * cos(2 * th)
  thHa   <- -sc / den
  thPhi  <- H * cos(th - phi) / den
  list(theta = th,
       lam   = H * cos(th - phi) + Ha * cos(th)^2,
       dHa   = c(th = thHa,  lam = cos(th)^2 - Ha * sc * thHa),
       dPhi  = c(th = thPhi, lam = H * sin(th - phi) - Ha * sc * thPhi))
}

## Signal: two-exponential background plus S (mz0 - mz) = -S K xi_z. Physics
## parameters and noise kernel are shared by all fields, background and noise
## level are per field. Positive parameters are fitted as log10, K and
## phi / 90 deg as log10(x / (1 - x)).
physics   <- c("log10Alpha", "log10G", "log10Ha", "log10TauK", "log10S",
               "logit10K", "logit10Phi")
kernel    <- paste0("w", 1:5)
bgNames   <- function(f) paste0(c("log10A1_", "log10Tau1_", "log10A2_",
                                  "log10Tau2_"), f)
sigmaName <- function(f) paste0("logSigma_", f)
ln10      <- log(10)

traceSignal <- function(p, k, tt) {
  H   <- 100 * fields[k]
  bn  <- bgNames(fields[k])
  K   <- toNatural("logit10K", p[["logit10K"]])
  S   <- 10^p[["log10S"]]
  Ha  <- 10^p[["log10Ha"]]
  u   <- toNatural("logit10Phi", p[["logit10Phi"]])
  phi <- pi / 2 * u                          # between 0 and 90 degrees
  eqm <- equilibrium(Ha, H, phi)
  pars <- c(xi_x = 0, xi_y = 0, xi_z = 0, kap = 1, K = K, th = eqm$theta,
            lam = eqm$lam, p[c("log10Alpha", "log10G", "log10Ha", "log10TauK")])
  sol <- tryCatch(solveODE(model, c(0, tt), pars, abstol = 1e-9, reltol = 1e-9),
                  error = function(e) NULL)
  if (is.null(sol)) return(NULL)

  xi      <- sol$variable[-1, "xi_z"]
  dxi     <- sol$tangent[-1, "xi_z", , drop = TRUE]
  viaRest <- function(d) dxi[, "th"] * d[["th"]] + dxi[, "lam"] * d[["lam"]]
  prec    <- -S * K * xi
  T1 <- 10^p[[bn[2]]]
  T2 <- 10^p[[bn[4]]]
  e1 <- 10^p[[bn[1]]] * exp(-tt / T1)
  e2 <- 10^p[[bn[3]]] * exp(-tt / T2)

  J <- matrix(0, length(tt), length(p), dimnames = list(NULL, names(p)))
  J[, "log10Alpha"] <- -S * K * dxi[, "log10Alpha"]
  J[, "log10G"]     <- -S * K * dxi[, "log10G"]
  J[, "log10Ha"]    <- -S * K * (dxi[, "log10Ha"] + ln10 * Ha * viaRest(eqm$dHa))
  J[, "log10TauK"]  <- -S * K * dxi[, "log10TauK"]
  J[, "log10S"]     <-  ln10 * prec
  J[, "logit10K"]   <- -S * (xi + K * dxi[, "K"]) * ln10 * K * (1 - K)
  J[, "logit10Phi"] <- -S * K * viaRest(eqm$dPhi) * pi / 2 * ln10 * u * (1 - u)
  J[, bn[1]] <- ln10 * e1
  J[, bn[2]] <- ln10 * e1 * tt / T1
  J[, bn[3]] <- ln10 * e2
  J[, bn[4]] <- ln10 * e2 * tt / T2
  list(y = e1 + e2 + prec, J = J)
}

## Noise: moving average of order 5, the lag at which the autocorrelation of
## the residuals decays. Kernel and noise levels are estimated. `innovations`
## are the whitened residuals. Moving-average models and the conditional sum
## of squares: Box, Jenkins, Reinsel, Ljung, Time Series Analysis, 5th ed.,
## Wiley (2015).
likelihood <- function(p) {
  r <- J <- innovations <- vector("list", length(fields))
  for (k in seq_along(fields)) {
    part <- traceSignal(p, k, obs[[k]]$time_ps)
    if (is.null(part)) return(NULL)
    sn <- sigmaName(fields[k])
    z  <- maResiduals(part$y - obs[[k]]$signal, part$J, p[[sn]], p[kernel])
    if (is.null(z)) return(NULL)
    z$J[, sn]     <- z$dSigma
    z$J[, kernel] <- z$dW
    r[[k]] <- z$r
    J[[k]] <- z$J
    innovations[[k]] <- z$r[seq_along(part$y)]
  }
  list(r = unlist(r), J = do.call(rbind, J),
       innovations = setNames(innovations, paste(fields, "kOe")))
}

## Gaussian prior on the physics parameters, standard deviation 3 decades,
## centred on order-of-magnitude values and the nominal field angle.
priorSd   <- 3
prior     <- c(log10Alpha = -1, log10G = log10(2), log10Ha = 3,
               log10TauK = log10(30), log10S = 0, logit10K = 0,
               logit10Phi = log10(0.5))
residuals <- withPrior(likelihood, prior, priorSd)

## -----------------------------------------------------------------
## 4. Multistart
## -----------------------------------------------------------------

set.seed(2)
nStart <- 100
perTrace <- do.call(cbind, lapply(fields, function(f) {
  tail <- obs[[as.character(f)]]$signal[obs[[as.character(f)]]$time_ps > 30]
  m <- cbind(runif(nStart, log10(0.3),  log10(3)),
             runif(nStart, log10(0.5),  log10(5)),
             runif(nStart, log10(0.05), log10(2)),
             runif(nStart, log10(2),    log10(30)),
             log(sd(diff(tail))) + runif(nStart, -0.5, 0.5))
  colnames(m) <- c(bgNames(f), sigmaName(f))
  m
}))

## Geometrically decaying kernels are invertible, so every start can be
## whitened.
decay <- runif(nStart, 0.2, 0.8)
kernelStarts <- outer(decay, seq_along(kernel), `^`)
colnames(kernelStarts) <- kernel
phiStart <- runif(nStart, 10, 60) / 90
starts <- cbind(log10Alpha = runif(nStart, log10(0.005), 0),
                log10G     = runif(nStart, log10(1.5),   log10(3)),
                log10Ha    = runif(nStart, 2,            log10(5000)),
                log10TauK  = runif(nStart, log10(2),     log10(500)),
                log10S     = runif(nStart, -1,           log10(50)),
                logit10K   = runif(nStart, -2, 1),
                logit10Phi = log10(phiStart / (1 - phiStart)),
                kernelStarts, perTrace)

ms <- multistart(residuals, starts, cores)

print(table(ms$status))
print(plotWaterfall(ms))

fit <- ms$best

## -----------------------------------------------------------------
## 5. Fit and residuals
## -----------------------------------------------------------------

fitted <- do.call(rbind, lapply(seq_along(fields), function(k) {
  tt <- seq(1, max(obs[[k]]$time_ps), length.out = 400)
  data.frame(field_kOe = fields[k], time_ps = tt,
             signal = traceSignal(fit$par, k, tt)$y)
}))

print(
  ggplot(raw, aes(time_ps, signal)) +
    geom_point(size = 0.5, alpha = 0.6) +
    geom_line(data = fitted, colour = "#b2182b", linewidth = 0.7) +
    signalAxis +
    facet_wrap(~field_kOe, ncol = 4, labeller = kOe) +
    labs(x = "pump-probe delay [ps]", y = "normalised Hall signal",
         title = "Fit", subtitle = "points measured, line fitted")
)

## Whitened residuals, Ljung-Box test (Ljung, Box, Biometrika 65, 297, 1978)
## and autocorrelation.
res <- likelihood(fit$par)$innovations

print(
  ggplot(data.frame(field_kOe = measured$field_kOe, time_ps = measured$time_ps,
                    r = unlist(res)),
         aes(time_ps, r)) +
    geom_hline(yintercept = c(-2, 2), colour = "grey70") +
    geom_hline(yintercept = 0, colour = "grey40") +
    geom_line(colour = "#2166ac", linewidth = 0.4) +
    facet_wrap(~field_kOe, ncol = 4, labeller = kOe) +
    labs(x = "pump-probe delay [ps]", y = "whitened residual",
         title = "Residuals")
)

print(signif(residualCheck(res), 3))
print(plotResidualACF(res))

labels <- c(log10Alpha = "alpha", log10G = "g", log10Ha = "Ha [mT]",
            log10TauK = "tauK [ps]", log10S = "S", logit10K = "anisotropy drop",
            logit10Phi = "phi / 90 deg")
estimates <- natural(fit$par[physics])
names(estimates) <- labels[physics]
print(signif(estimates, 4))

noiseLevel <- setNames(exp(fit$par[sigmaName(fields)]), paste(fields, "kOe"))
print(signif(noiseLevel, 3))
print(signif(fit$par[kernel], 3))

## -----------------------------------------------------------------
## 6. Profile likelihood
## -----------------------------------------------------------------

## For the profiles the prior is centred on the best fit with standard
## deviation 10 decades, and the fit is repeated.
priorSd   <- 10
prior[]   <- fit$par[names(prior)]
residuals <- withPrior(likelihood, prior, priorSd)
wide      <- gaussNewton(residuals, fit$par)

## At most 100 steps and 5 decades on each side of the best fit.
prof <- profiles(residuals, wide$par, which = physics, range = 4,
                 maxSteps = 100, cores = cores)
print(plotProfiles(prof, labels = labels))

## 95 percent intervals from the data part of the profiles.
intervals <- profileIntervals(prof)
rownames(intervals) <- labels[rownames(intervals)]
print(signif(intervals, 4))
cat(sprintf("field angle between %.0f and %.0f degrees\n",
            90 * intervals["phi / 90 deg", "lower"],
            90 * intervals["phi / 90 deg", "upper"]))

## -----------------------------------------------------------------
## 7. Precession frequencies compared with the paper
## -----------------------------------------------------------------

## Small-angle precession frequency (Smit-Beljers) from the fitted parameters,
## compared with the frequencies from damped-sine fits in the paper.
fmr <- function(H, Ha, alpha, g, phi) {
  th <- equilibrium(Ha, H, phi)$theta
  w  <- g * muB * sqrt((H * cos(th - phi) + Ha * cos(2 * th)) * H * sin(phi) / sin(th))
  1000 * w / (2 * pi * (1 + alpha^2))
}
paper <- read.csv(system.file("extdata", "zhang2025_frequencies.csv",
                              package = "cppDE"))
paper$model_GHz <- sapply(100 * paper$field_kOe, fmr,
                          Ha    = estimates[["Ha [mT]"]],
                          alpha = estimates[["alpha"]],
                          g     = estimates[["g"]],
                          phi   = pi / 2 * estimates[["phi / 90 deg"]])
print(paper)
