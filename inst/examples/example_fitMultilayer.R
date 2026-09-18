## =================================================================
## Fitting 150 parameters of a magnetic multilayer with reverse-mode gradients
## =================================================================

## A [Co/Pt] multilayer heated by a femtosecond laser pulse: electron and
## lattice temperature on 100 cells in depth, and the magnetisation of each Co
## layer by the microscopic three-temperature model (Koopmans et al., Nat.
## Mater. 9, 259, 2010). 248 states, 150 parameters, simulated data from three
## experiments. optim's L-BFGS-B fits them on reverse-mode gradients; section 4
## times the reverse against the forward gradient.

rm(list = ls(all.names = TRUE))
setwd(tempdir())

library(cppDE)
library(ggplot2)

## Helpers: waterfall plot.
source(system.file("examples", "fitTools.R", package = "cppDE"))

theme_set(theme_bw(base_size = 11))
cores <- if (.Platform$OS.type == "windows") 1L else min(10L, parallel::detectCores())
ln10  <- log(10)

## -----------------------------------------------------------------
## 1. Stack
## -----------------------------------------------------------------

## A 3 nm Pt cap on twelve bilayers of 1 nm Co and 2 nm Pt, four cells per
## layer. Interface L lies between layer L and layer L + 1.
material <- c("Pt", rep(c("Co", "Pt"), 12))
thick    <- ifelse(material == "Co", 1, 2)
thick[1] <- 3                                               # nm
nLayer   <- length(material)
coLayers <- which(material == "Co")

nSub    <- 4L
layerOf <- rep(seq_len(nLayer), each = nSub)
dz      <- rep(thick / nSub, each = nSub)                   # nm
zc      <- cumsum(dz) - dz / 2                              # cell centres, nm
nCell   <- length(dz)
coCells <- which(material[layerOf] == "Co")

## Held fixed: the electron conductivity at Te = Tp, in W/(m K), the lattice
## conductivity, and the absorbed fluences. The stack is pumped from the top at
## two fluences and once from the bottom, through the substrate.
kappaE <- c(Co = 100, Pt = 72)[material]
kappaP <- 5
pump   <- rbind(top_low  = c(top = 6, bottom = 0),
                top_high = c(top = 12, bottom = 0),
                bottom   = c(top = 0, bottom = 12))         # J/m2
experiments <- rownames(pump)
depth  <- sum(thick)

TeN <- sprintf("Te%03d", seq_len(nCell))
TpN <- sprintf("Tp%03d", seq_len(nCell))
lmN <- sprintf("lm%03d", coCells)
nm  <- function(name, i) sprintf("%s%02d", name, i)

## -----------------------------------------------------------------
## 2. Model
## -----------------------------------------------------------------

## Parameters, all as log10:
##   per layer      gamma  electron heat capacity Te*gamma, J/(m3 K2)
##                  Cp     lattice heat capacity, MJ/(m3 K)
##                  G      electron-lattice coupling, 1e17 W/(m3 K)
##   per Co layer   R      demagnetisation rate, 1/ps
##                  Tc     Curie temperature, K
##   per interface  He, Hp electron and lattice conductance, GW/(m2 K)
##   per material   Tau    thermalisation time, ps
##   global         Lambda absorption depth, nm
parNames <- c(nm("lgGamma", seq_len(nLayer)), nm("lgCp", seq_len(nLayer)),
              nm("lgG", seq_len(nLayer)), nm("lgR", coLayers),
              nm("lgTc", coLayers), nm("lgHe", seq_len(nLayer - 1)),
              nm("lgHp", seq_len(nLayer - 1)), "lgTauCo", "lgTauPt", "lgLambda")
cat(sprintf("%d parameters\n", length(parNames)))

## Conductance per area between neighbouring cells a and b, in W/(m2 K): two
## half cells in series, and the interface where they belong to different
## layers. The electron conductivity scales with Te/Tp.
conductance <- function(a, b, carrier) {
  La <- layerOf[a]
  Lb <- layerOf[b]
  if (carrier == "e") {
    ka <- sprintf("%g*%s/%s", kappaE[La], TeN[a], TpN[a])
    kb <- sprintf("%g*%s/%s", kappaE[Lb], TeN[b], TpN[b])
    h  <- sprintf("1e9*10^%s", nm("lgHe", La))
  } else {
    ka <- kb <- kappaP
    h  <- sprintf("1e9*10^%s", nm("lgHp", La))
  }
  cells <- sprintf("%g/(%s) + %g/(%s)", 0.5e-9 * dz[a], ka, 0.5e-9 * dz[b], kb)
  if (La == Lb) return(sprintf("1/(%s)", cells))
  sprintf("1/(%s + 1/(%s))", cells, h)
}

## Net heat flow into cell i from its neighbours, in W/m2; the top and the
## bottom of the stack are insulating.
inflow <- function(i, carrier, T) {
  below <- if (i < nCell)
    sprintf("%s*(%s - %s)", conductance(i, i + 1, carrier), T[i + 1], T[i]) else "0"
  above <- if (i > 1)
    sprintf("%s*(%s - %s)", conductance(i - 1, i, carrier), T[i], T[i - 1]) else "0"
  sprintf("(%s - %s)", below, above)
}

## The pulse is absorbed at time zero, decaying as exp(-z/Lambda) from the
## pumped side, and reaches Te within Tau. The source is largest at the start, so
## the solver cannot step over it. `top` and `bottom` are the fluences.
heating <- function(i) {
  tau <- paste0("10^lgTau", material[layerOf[i]])
  sprintf(paste0("(top*exp(-%g/10^lgLambda) + bottom*exp(-%g/10^lgLambda))",
                 "*1e21/10^lgLambda*exp(-time/%s)/%s"),
          zc[i], depth - zc[i], tau, tau)
}

## Electron and lattice temperature in K/ps: power densities in W/m3 over the
## heat capacities in J/(m3 K), times 1e-12 s/ps.
eqTe <- vapply(seq_len(nCell), function(i) {
  L <- layerOf[i]
  sprintf("1e-12*(%s/%g + %s - 1e17*10^%s*(%s - %s))/(10^%s*%s)",
          inflow(i, "e", TeN), 1e-9 * dz[i], heating(i), nm("lgG", L),
          TeN[i], TpN[i], nm("lgGamma", L), TeN[i])
}, "")
eqTp <- vapply(seq_len(nCell), function(i) {
  L <- layerOf[i]
  sprintf("1e-12*(%s/%g + 1e17*10^%s*(%s - %s))/(1e6*10^%s)",
          inflow(i, "p", TpN), 1e-9 * dz[i], nm("lgG", L), TeN[i], TpN[i],
          nm("lgCp", L))
}, "")

## Magnetisation of the Co cells, microscopic three-temperature model:
## dm/dt = R m (Tp/Tc) (1 - m coth(m Tc/Te)). The state is lm = log m, which
## keeps its precision and its sign when a layer demagnetises almost fully.
eqM <- vapply(seq_along(coCells), function(k) {
  i <- coCells[k]
  L <- layerOf[i]
  sprintf("10^%s*%s/10^%s*(1 - exp(%s)/tanh(exp(%s)*10^%s/%s))",
          nm("lgR", L), TpN[i], nm("lgTc", L), lmN[k], lmN[k],
          nm("lgTc", L), TeN[i])
}, "")

eqs <- c(setNames(eqTe, TeN), setNames(eqTp, TpN), setNames(eqM, lmN))

tBuild <- system.time(
  model <- cppODE(eqs, modelname = "multilayer", method = "bdf",
                  derivMode = "reverse")
)
cat(sprintf("%d states, code generation and compilation %.0f s\n",
            length(eqs), tBuild[["elapsed"]]))

## Before the pulse everything sits at 300 K, and each Co layer at its
## mean-field magnetisation m = tanh(m Tc/T). `dlm` is the derivative of its
## log in log10 Tc, by implicit differentiation.
T0 <- 300
magnetisation <- function(lgTc) {
  Tc <- 10^lgTc
  m  <- uniroot(function(m) m - tanh(m * Tc / T0), c(1e-3, 1), tol = 1e-14)$root
  s2 <- 1 - tanh(m * Tc / T0)^2
  list(lm = log(m), dlm = s2 / T0 / (1 - s2 * Tc / T0) * ln10 * Tc)
}

initialState <- function(p) {
  mag <- lapply(coLayers, function(L) magnetisation(p[[nm("lgTc", L)]]))
  lm0 <- vapply(mag, `[[`, 0, "lm")[match(layerOf[coCells], coLayers)]
  c(setNames(rep(T0, nCell), TeN), setNames(rep(T0, nCell), TpN),
    setNames(lm0, lmN))
}

## -----------------------------------------------------------------
## 3. Simulated data
## -----------------------------------------------------------------

## Literature values of Co and Pt; the true parameters scatter by 12 percent
## around them, per layer and per interface.
isCo    <- material == "Co"
isTc    <- startsWith(parNames, "lgTc")
nominal <- setNames(log10(c(
  ifelse(isCo, 700, 750), ifelse(isCo, 3.7, 2.85), ifelse(isCo, 4, 2.5),
  rep(25, length(coLayers)), rep(1000, length(coLayers)),
  rep(20, nLayer - 1), rep(0.5, nLayer - 1), 0.1, 0.2, 12)), parNames)
set.seed(1)
truth <- nominal + rnorm(length(nominal), 0, ifelse(isTc, 0.03, 0.05))

## Observed: the mean over the cells of a layer, for Te and Tp of every layer
## and m of every Co layer, on a logarithmic grid up to 100 ps.
times  <- c(0, 10^seq(log10(0.02), log10(100), length.out = 200))
traces <- c(paste0("Te", seq_len(nLayer)), paste0("Tp", seq_len(nLayer)),
            paste0("m", coLayers))
A <- matrix(0, length(traces), length(eqs), dimnames = list(traces, names(eqs)))
for (L in seq_len(nLayer)) {
  cells <- which(layerOf == L)
  A[paste0("Te", L), TeN[cells]] <- 1 / nSub
  A[paste0("Tp", L), TpN[cells]] <- 1 / nSub
  if (isCo[L]) A[paste0("m", L), lmN[match(cells, coCells)]] <- 1 / nSub
}
sigma <- ifelse(startsWith(traces, "m"), 0.002, 1)          # K; m has no unit

## The traces from a solution x, with m = exp(lm), and the derivative of
## sum(w * traces) in the states along x.
observe <- function(x) {
  x[, lmN] <- exp(x[, lmN])
  tcrossprod(x, A)
}
dObserve <- function(x, w) {
  d <- w %*% A
  d[, lmN] <- d[, lmN] * exp(x[, lmN])
  d
}

## A solve for values keeps the checkpoints a later reverse sweep can use. One
## that needs more than 1e4 steps, twenty times the usual, counts as failed.
parsOf <- function(p, e) c(initialState(p), p[parNames], pump[e, ])
solve1 <- function(p, e, tol)
  solveODE(model, times, parsOf(p, e), abstol = tol, reltol = tol,
           maxsteps = 1e4, keepStore = TRUE)

data <- sapply(experiments, function(e) {
  clean <- observe(solve1(truth, e, 1e-10)$variable[, names(eqs)])
  clean + matrix(rnorm(length(clean), 0, rep(sigma, each = nrow(clean))),
                 nrow(clean))
}, simplify = FALSE)

## -----------------------------------------------------------------
## 4. Objective and gradient
## -----------------------------------------------------------------

## Half the sum of squared residuals over the noise level. The gradient adds one
## reverse sweep per experiment to the forward solves, seeded with the residuals.
## L-BFGS-B stops on relative changes of about 2e-9, so the solves run at 1e-10.
tol <- 1e-10

evaluateAt <- function(p) {
  runs <- lapply(experiments, function(e) {
    sol <- tryCatch(solve1(p, e, tol), error = function(err) NULL)
    if (is.null(sol)) return(NULL)
    r <- (observe(sol$variable[, names(eqs)]) - data[[e]]) /
      rep(sigma, each = length(times))
    list(sol = sol, r = r)
  })
  if (any(vapply(runs, is.null, TRUE))) return(NULL)
  list(p = p, runs = setNames(runs, experiments),
       value = 0.5 * sum(vapply(runs, function(u) sum(u$r^2), 0)))
}

gradientAt <- function(o) {
  g <- setNames(numeric(length(parNames)), parNames)
  for (e in experiments) {
    u    <- o$runs[[e]]
    seed <- dObserve(u$sol$variable[, names(eqs)],
                     u$r / rep(sigma, each = length(times)))
    adj  <- solveODE(model, times, parsOf(o$p, e), abstol = tol, reltol = tol,
                     seed = seed, store = u$sol$store)$adjoint[, 1]
    g <- g + adj[parNames]
    for (L in coLayers) {
      cells <- lmN[layerOf[coCells] == L]
      g[[nm("lgTc", L)]] <- g[[nm("lgTc", L)]] +
        sum(adj[cells]) * magnetisation(o$p[[nm("lgTc", L)]])$dlm
    }
  }
  g
}

## optim asks for the value and the gradient at the same point; both share the
## forward solves. R does not see the size of the checkpoints, so those of the
## previous point are collected first.
cache <- new.env(parent = emptyenv())
evaluate <- function(p) {
  names(p) <- parNames
  if (!identical(cache$p, p)) {
    cache$out <- NULL
    invisible(gc())
    cache$p   <- p
    cache$out <- evaluateAt(p)
  }
  cache$out
}

## L-BFGS-B needs a finite value; a failed solve returns 1e20, above any start.
fn <- function(p) {
  o <- evaluate(p)
  if (is.null(o)) 1e20 else o$value
}
gr <- function(p) {
  o <- evaluate(p)
  if (is.null(o)) numeric(length(p)) else gradientAt(o)
}

## The same gradient from forward sensitivities, for comparison: one tangent
## per parameter, 150 in all, seeded so that the initial magnetisation moves
## with Tc.
modelF <- cppODE(eqs, modelname = "multilayerF", method = "bdf",
                 derivMode = "forward")

gradientForward <- function(p) {
  S <- matrix(0, length(parNames) + length(lmN), length(parNames),
              dimnames = list(c(parNames, lmN), parNames))
  S[cbind(parNames, parNames)] <- 1
  for (L in coLayers)
    S[lmN[layerOf[coCells] == L], nm("lgTc", L)] <-
      magnetisation(p[[nm("lgTc", L)]])$dlm
  g <- setNames(numeric(length(parNames)), parNames)
  for (e in experiments) {
    sol <- solveODE(modelF, times, parsOf(p, e), sens1ini = S,
                    abstol = tol, reltol = tol)
    x <- sol$variable[, names(eqs)]
    w <- dObserve(x, (observe(x) - data[[e]]) / rep(sigma^2, each = length(times)))
    g <- g + vapply(parNames, function(j) sum(w * sol$sens1[, names(eqs), j]), 0)
  }
  g
}

tFn  <- system.time(f0 <- fn(truth))[["elapsed"]]
tRev <- system.time(gRev <- gr(truth))[["elapsed"]]
tFwd <- system.time(gFwd <- gradientForward(truth))[["elapsed"]]
cat(sprintf("value %.2f s; gradient on top: reverse %.2f s, forward %.2f s\n",
            tFn, tRev, tFwd))
cat(sprintf("largest relative difference of the two gradients: %.2g\n",
            max(abs(gRev - gFwd)) / max(abs(gFwd))))

## -----------------------------------------------------------------
## 5. Multistart
## -----------------------------------------------------------------

## Each start moves every parameter by up to a factor of 2 from the truth. The
## Curie temperatures stay within 400 to 2500 K, above the 300 K of the sample;
## the other parameters are free.
nStart <- 20
starts <- t(replicate(nStart, truth + runif(length(truth), -0.3, 0.3)))
lower  <- ifelse(isTc, log10(400), -Inf)
upper  <- ifelse(isTc, log10(2500), Inf)

## L-BFGS-B keeps 20 steps for its curvature model instead of 5: the deep layers
## and the outer interfaces are determined weakly and converge slower with fewer.
control <- list(maxit = 10000, lmm = 20)

## The multistart result ships with the package and is read back unless
## `refit` is TRUE.
refit   <- FALSE
fitFile <- system.file("examples", "fitMultilayer.rds", package = "cppDE")
if (refit || !nzchar(fitFile)) {
  tFit <- system.time(
    fits <- parallel::mclapply(seq_len(nStart), function(s)
      tryCatch(optim(starts[s, ], fn, gr, method = "L-BFGS-B",
                     lower = lower, upper = upper,
                     control = control),
               error = function(e) NULL),
      mc.cores = cores, mc.preschedule = FALSE)
  )
  multistart <- list(starts = starts, fits = fits,
                     seconds = tFit[["elapsed"]], cores = cores)
} else {
  multistart <- readRDS(fitFile)
}
fits <- multistart$fits

values    <- vapply(fits, function(f) if (is.null(f)) NA_real_ else f$value,
                    numeric(1))
## L-BFGS-B reports 1 at the iteration limit. Codes 0 (the value stops dropping)
## and 52 (the line search finds no lower point) end the search by itself.
atLimit   <- vapply(fits, function(f) !is.null(f) && f$convergence == 1,
                    logical(1))
iBest     <- which.min(values)
status    <- ifelse(values - values[iBest] < 1, "global",
                    ifelse(atLimit, "iteration limit", "local"))
status[is.na(values)] <- "failed"

cat(sprintf("multistart: %d starts, %.0f s on %d cores\n", nStart,
            multistart$seconds, multistart$cores))
print(table(status))
print(plotWaterfall(list(value = values, status = status)))

fit <- fits[[iBest]]
names(fit$par) <- parNames

## -----------------------------------------------------------------
## 6. Recovered parameters and fit
## -----------------------------------------------------------------

## One panel per parameter group; the three global parameters share one.
group <- sub("[0-9]+$", "", parNames)
group[group %in% c("lgTauCo", "lgTauPt", "lgLambda")] <- "lgTau, lgLambda"
recovered <- data.frame(parameter = parNames, group = group,
                        truth = truth, start = starts[iBest, ],
                        estimate = fit$par)

print(
  ggplot(recovered, aes(truth, estimate)) +
    geom_abline(colour = "grey60") +
    geom_point(aes(y = start), colour = "grey75", size = 1) +
    geom_point(colour = "#b2182b", size = 1.2) +
    facet_wrap(~group, scales = "free") +
    labs(x = "true value", y = "estimate",
         title = "150 parameters recovered",
         subtitle = "grey: start, red: estimate, line: truth (log10 units)")
)
cat("largest |log10 error| per group:\n")
print(signif(tapply(abs(fit$par - truth), group, max), 2))

best <- evaluate(fit$par)
show <- c("Te1", "Te2", "Tp2", "m2", "Te24", "m24")
long <- do.call(rbind, lapply(experiments, function(e) {
  model_ <- data[[e]] + best$runs[[e]]$r * rep(sigma, each = length(times))
  do.call(rbind, lapply(show, function(tr)
    data.frame(experiment = e, trace = tr, time = times,
               measured = data[[e]][, tr], fitted = model_[, tr])))
}))

print(
  ggplot(long[long$time > 0, ], aes(time, measured, colour = experiment)) +
    geom_point(size = 0.5, alpha = 0.5) +
    geom_line(aes(y = fitted)) +
    scale_x_log10() +
    facet_wrap(~trace, scales = "free_y") +
    labs(x = "time after the pulse [ps]", y = NULL, colour = "pumped",
         title = "Simulated data and fit",
         subtitle = "layer 1: Pt cap, 2: first Co, 24: last Co")
)
