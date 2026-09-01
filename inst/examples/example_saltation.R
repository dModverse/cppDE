## =================================================================
## Saltation transport, checked against a symbolic solution
##
## Every model below is solved twice. cppODE() integrates it with first
## and second order sensitivities, and SymPy solves each smooth segment
## with dsolve(), glues the segments at the firing time and
## differentiates the composition twice. The two have to agree panel by
## panel, including the parameter dependence of the firing time itself.
##
## In every figure the symbolic solution is the solid line and what
## solveODE() returned is dashed on top of it.
## =================================================================
rm(list = ls(all.names = TRUE))
setwd(tempdir())

library(cppDE)
library(ggplot2)
sp <- reticulate::import("sympy")

## The base pdf() device has no glyph for the derivative sign in the panel
## titles. Interactively the plot pane draws it, a scripted run needs cairo.
if (!interactive()) cairo_pdf("saltation.pdf", onefile = TRUE, width = 10, height = 7)

## -----------------------------------------------------------------
## Symbolic layer
## -----------------------------------------------------------------
sym  <- function(name, ...) sp$Symbol(name, real = TRUE, ...)
tSym <- sym("t")

## Integer exponents only. A float power blocks the simplifications that
## turn the solution of an oscillator back into a sine.
pow <- function(base, exponent) sp$Pow(base, as.integer(exponent))

## One smooth segment of a first order equation, pinned by y(t0) = y0.
segment <- function(fun, rhs, t0, y0)
  sp$dsolve(sp$Eq(sp$Derivative(fun(tSym), tSym), rhs), fun(tSym),
            ics = reticulate::py_dict(list(fun(t0)), list(y0)))$rhs

## The same for a second order equation, pinned by y(t0) and y'(t0).
segment2 <- function(fun, rhs, t0, y0, dy0)
  sp$dsolve(sp$Eq(sp$Derivative(fun(tSym), tSym, 2L), rhs), fun(tSym),
            ics = reticulate::py_dict(
              list(fun(t0), sp$Derivative(fun(tSym), tSym)$subs(tSym, t0)),
              list(y0, dy0)))$rhs

## Segments glued at their firing times, one branch longer than fires.
composed <- function(branches, fires) {
  pieces <- Map(function(branch, fire) reticulate::tuple(branch, sp$Lt(tSym, fire)),
                branches[seq_along(fires)], fires)
  do.call(sp$Piecewise,
          c(pieces, list(reticulate::tuple(branches[[length(branches)]], sp$`true`))))
}

## A parameter point: the symbols, their values, and the grid to compare on.
## The names are the ones solveODE() labels its sensitivity slices with, where
## a state name stands for that state's initial value.
point <- function(syms, pars, times)
  list(syms = syms, names = names(pars), values = as.list(unname(pars)),
       times = times)

toNumeric <- function(expr)
  suppressWarnings(as.numeric(reticulate::py_str(sp$N(expr))))

## Every parameter replaced by its value, time left free.
atPoint <- function(expr, at) expr$subs(reticulate::py_dict(at$syms, at$values))

## SymPy's lambdify reads its caller's frame to build a docstring, and a call
## from R has no Python frame to read. Going through a Python function gives it
## one, and compiling once per expression turns one round trip per time point
## into one per panel.
sympyGrid <- reticulate::py_run_string(
"import sympy
def over_grid(sym, expr, grid):
    f = sympy.lambdify(sym, expr, 'math')
    return [f(v) for v in grid]
")

## The expression, or a derivative of it, over the whole time grid. A constant
## expression still yields one value per time point.
series <- function(expr, at, wrt = list()) {
  for (s in wrt) expr <- sp$diff(expr, s)
  as.numeric(sympyGrid$over_grid(tSym, atPoint(expr, at), at$times))
}

## Numeric value of an expression that no longer depends on time.
value <- function(expr, at) toNumeric(atPoint(expr, at))

## The first zero after `after`, over one condition or several at once, picked
## numerically among the symbolic roots. Complex branches drop out as NA.
crossing <- function(conditions, at, after = 0) {
  if (!is.list(conditions)) conditions <- list(conditions)
  roots <- unlist(lapply(conditions,
                         function(g) sp$solve(sp$Eq(g, sp$Integer(0L)), tSym)),
                  recursive = FALSE)
  when <- vapply(roots, value, 0, at = at)
  roots[[which.min(ifelse(!is.na(when) & when > after + 1e-9, when, Inf))]]
}

## -----------------------------------------------------------------
## Comparison and figures
## -----------------------------------------------------------------
## Panel labels carry the partial derivative sign, U+2202.
lab1 <- function(var, a) sprintf("\u2202%s / \u2202%s", var, a)
lab2 <- function(var, a, b) {
  if (a == b) return(sprintf("\u2202²%s / \u2202%s²", var, a))
  sprintf("\u2202²%s / \u2202%s\u2202%s", var, a, b)
}

## The last row at a requested time carries the post-event state, whether
## or not a localised root inserted rows of its own.
rowsAt <- function(res, times)
  vapply(times, function(s) max(which(abs(res$time - s) < 1e-9)), 1L)

## One panel per quantity: the state, its gradient, and each distinct
## second derivative, symbolic against what solveODE() returned.
comparison <- function(expr, at, res, var) {
  i  <- rowsAt(res, at$times)
  np <- length(at$syms)

  ## One entry per panel, each holding both curves over the whole grid.
  panel <- function(order, quantity, symbolic, numeric)
    list(order = order, quantity = quantity,
         symbolic = symbolic, numeric = as.numeric(numeric))

  panels <- list(panel(0L, var, series(expr, at), res$variable[i, var]))
  for (a in seq_len(np))
    panels[[length(panels) + 1L]] <-
      panel(1L, lab1(var, at$names[a]), series(expr, at, at$syms[a]),
            res$sens1[i, var, at$names[a]])
  for (a in seq_len(np)) for (b in a:np)
    panels[[length(panels) + 1L]] <-
      panel(2L, lab2(var, at$names[a], at$names[b]),
            series(expr, at, at$syms[c(a, b)]),
            res$sens2[i, var, at$names[a], at$names[b]])

  labels <- vapply(panels, `[[`, "", "quantity")
  data.frame(
    order    = rep(vapply(panels, `[[`, 0L, "order"), each = length(at$times)),
    quantity = factor(rep(labels, each = length(at$times)), levels = labels),
    time     = rep(at$times, length(panels)),
    symbolic = unlist(lapply(panels, `[[`, "symbolic"), use.names = FALSE),
    numeric  = unlist(lapply(panels, `[[`, "numeric"), use.names = FALSE))
}

## The numeric side of the figures, one line per derivative order.
report <- function(df, title) {
  cat("\n", title, "\n", sep = "")
  for (o in sort(unique(df$order))) {
    d <- df[df$order == o, ]
    cat(sprintf("  %-13s max |solveODE - symbolic| = %.2e\n",
                c("state", "first order", "second order")[o + 1L],
                max(abs(d$numeric - d$symbolic))))
  }
}

comparisonPlot <- function(df, orders, title)
  ggplot(df[df$order %in% orders, ], aes(x = time)) +
    ## A quantity that vanishes identically leaves the solver a cancellation
    ## residue of a few ulp, and a free y scale would blow that up into a
    ## picture of its own round-off. No panel is drawn tighter than this.
    expand_limits(y = c(-1e-10, 1e-10)) +
    geom_line(aes(y = symbolic, colour = "symbolic", linetype = "symbolic"),
              linewidth = 1) +
    geom_line(aes(y = numeric, colour = "solveODE", linetype = "solveODE"),
              linewidth = 1) +
    facet_wrap(~quantity, scales = "free_y") +
    scale_colour_manual(values = c(symbolic = "#2166ac", solveODE = "#d6604d")) +
    scale_linetype_manual(values = c(symbolic = "solid", solveODE = "dotdash")) +
    labs(title = title, x = "time", y = NULL, colour = NULL, linetype = NULL) +
    theme_bw(base_size = 9) +
    theme(legend.position = "top")

## =================================================================
## 1. Root event on a coupled pair
##
## S' = a and C' = -b C, with d added to C when S reaches c. The firing
## time t* = (c - S0)/a depends on three parameters, so every
## sensitivity of C after the event carries a saltation term. With a = 1
## the crossing lands exactly on a requested time.
## =================================================================
events_pair <- data.frame(var = "C", time = NA, value = "d", method = "add",
                          root = "S - c", stringsAsFactors = FALSE)
model_pair  <- cppODE(c(S = "a", C = "-b * C"), events = events_pair,
                      deriv = TRUE, deriv2 = TRUE, modelname = "saltation_pair")

pars_pair  <- c(S = 2, C = 4, a = 1, b = 0.1, c = 14, d = 3)
times_pair <- seq(0, 20, len = 1000)
## Tight tolerances throughout, so that what is left of the residual is the
## solver and not the saltation transport. roottol matters as much as the two
## others here: the firing time enters the state through the velocity at the
## event, so a loosely localised root shows up as a state error.
res_pair <- solveODE(model_pair, times_pair, pars_pair,
                     abstol = 1e-12, reltol = 1e-12, roottol = 1e-12)

S0 <- sym("S0"); C0 <- sym("C0"); a <- sym("a")
b  <- sym("b");  cc <- sym("c");  d <- sym("d")
Sfun <- sp$Function("S"); Cfun <- sp$Function("C")

point_pair <- point(list(S0, C0, a, b, cc, d),
                    pars_pair[c("S", "C", "a", "b", "c", "d")], times_pair)
S_of_t <- segment(Sfun, a, sp$Integer(0L), S0)
fire_pair <- crossing(S_of_t - cc, point_pair)
C_before  <- segment(Cfun, -b * Cfun(tSym), sp$Integer(0L), C0)
C_after   <- segment(Cfun, -b * Cfun(tSym), fire_pair,
                     C_before$subs(tSym, fire_pair) + d)
C_of_t    <- composed(list(C_before, C_after), list(fire_pair))

cmp_pair <- comparison(C_of_t, point_pair, res_pair, "C")
report(cmp_pair, "1. root event on a coupled pair")

plot_pair_1 <- comparisonPlot(cmp_pair, 0:1, "C and its gradient across a root event")
plot_pair_1
plot_pair_2 <- comparisonPlot(cmp_pair, 2L, "Second derivatives of C across a root event")
plot_pair_2

## =================================================================
## 2. Time-triggered event
##
## x' = -k x with v added at the parameter time te. Here dt*/dtheta is
## supplied directly by the user expression rather than reconstructed
## from a root condition, which is the other saltation path.
## =================================================================
events_dose <- data.frame(var = "x", time = "te", value = "v", method = "add",
                          root = NA, stringsAsFactors = FALSE)
model_dose  <- cppODE(c(x = "-k * x"), events = events_dose, deriv = TRUE,
                      deriv2 = TRUE, modelname = "saltation_dose")

pars_dose  <- c(x = 1, k = 0.3, v = 2, te = 4)
times_dose <- seq(0, 10, len = 1000)
res_dose   <- solveODE(model_dose, times_dose, pars_dose,
                       abstol = 1e-12, reltol = 1e-12)

x0 <- sym("x0"); k <- sym("k"); dose <- sym("v"); te <- sym("te")
xfun <- sp$Function("x")

point_dose <- point(list(x0, k, dose, te),
                    pars_dose[c("x", "k", "v", "te")], times_dose)
x_before <- segment(xfun, -k * xfun(tSym), sp$Integer(0L), x0)
x_after  <- segment(xfun, -k * xfun(tSym), te, x_before$subs(tSym, te) + dose)
x_dose   <- composed(list(x_before, x_after), list(te))

cmp_dose <- comparison(x_dose, point_dose, res_dose, "x")
report(cmp_dose, "2. time-triggered event")

plot_dose_1 <- comparisonPlot(cmp_dose, 0:1, "x and its gradient across a timed dose")
plot_dose_1
plot_dose_2 <- comparisonPlot(cmp_dose, 2L, "Second derivatives of x across a timed dose")
plot_dose_2

## =================================================================
## 3. A root event that fires three times
##
## x' = -k x with v added whenever x falls back to xc. Each firing time
## sits on top of the previous one, so the saltation corrections
## compose and the sensitivities pick up all three.
## =================================================================
events_pulse <- data.frame(var = "x", time = NA, value = "v", method = "add",
                           root = "x - xc", stringsAsFactors = FALSE)
model_pulse  <- cppODE(c(x = "-k * x"), events = events_pulse, deriv = TRUE,
                       deriv2 = TRUE, modelname = "saltation_pulse")

pars_pulse  <- c(x = 1, k = 0.3, v = 1, xc = 0.4)
times_pulse <- seq(0, 12, len = 1000)
res_pulse   <- solveODE(model_pulse, times_pulse, pars_pulse, maxroot = 3L,
                        abstol = 1e-12, reltol = 1e-12, roottol = 1e-12)

xc <- sym("xc")
point_pulse <- point(list(x0, k, dose, xc),
                     pars_pulse[c("x", "k", "v", "xc")], times_pulse)

branches_pulse <- list(segment(xfun, -k * xfun(tSym), sp$Integer(0L), x0))
fires_pulse    <- list()
for (j in 1:3) {
  previous <- if (j == 1) 0 else value(fires_pulse[[j - 1]], point_pulse)
  fires_pulse[[j]] <- crossing(branches_pulse[[j]] - xc, point_pulse, previous)
  branches_pulse[[j + 1]] <- segment(xfun, -k * xfun(tSym), fires_pulse[[j]],
                                     xc + dose)
}
x_pulse <- composed(branches_pulse, fires_pulse)

cmp_pulse <- comparison(x_pulse, point_pulse, res_pulse, "x")
report(cmp_pulse, "3. a root event that fires three times")

plot_pulse_1 <- comparisonPlot(cmp_pulse, 0:1, "x and its gradient over three firings")
plot_pulse_1
plot_pulse_2 <- comparisonPlot(cmp_pulse, 2L, "Second derivatives of x over three firings")
plot_pulse_2

## =================================================================
## 4. Explicit time dependence in the right-hand side
##
## x' = -k t x with v added when x falls to xc. The right-hand side
## depends on time itself, so the second order transport needs the
## time argument of the Heun legs and not only the shifted state.
## =================================================================
events_drift <- data.frame(var = "x", time = NA, value = "v", method = "add",
                           root = "x - xc", stringsAsFactors = FALSE)
model_drift  <- cppODE(c(x = "-k * time * x"), events = events_drift,
                       deriv = TRUE, deriv2 = TRUE, modelname = "saltation_drift")

pars_drift  <- c(x = 1, k = 0.2, v = 0.5, xc = 0.3)
times_drift <- seq(0, 6, len = 1000)
res_drift   <- solveODE(model_drift, times_drift, pars_drift,
                        abstol = 1e-12, reltol = 1e-12, roottol = 1e-12)

point_drift <- point(list(x0, k, dose, xc),
                     pars_drift[c("x", "k", "v", "xc")], times_drift)
drift_before <- segment(xfun, -k * tSym * xfun(tSym), sp$Integer(0L), x0)
fire_drift   <- crossing(drift_before - xc, point_drift)
drift_after  <- segment(xfun, -k * tSym * xfun(tSym), fire_drift, xc + dose)
x_drift      <- composed(list(drift_before, drift_after), list(fire_drift))

cmp_drift <- comparison(x_drift, point_drift, res_drift, "x")
report(cmp_drift, "4. explicit time dependence in the right-hand side")

plot_drift_1 <- comparisonPlot(cmp_drift, 0:1, "x and its gradient, time-dependent decay")
plot_drift_1
plot_drift_2 <- comparisonPlot(cmp_drift, 2L, "Second derivatives of x, time-dependent decay")
plot_drift_2

## =================================================================
## 5. Oscillator between two elastic walls
##
## x'' = -w^2 x written as the pair (x, v), with a wall at +L and one at
## -L. Each wall is a root event that turns the velocity around, and the
## two fire alternately. The position stays continuous at a bounce, so
## everything the event does to its sensitivities is saltation. The wall
## it just left also holds that root at exactly zero when the solver
## restarts, which must not count as the next crossing.
## =================================================================
events_wall <- data.frame(var = c("v", "v"), time = c(NA, NA),
                          value = c("-1", "-1"),
                          method = c("multiply", "multiply"),
                          root = c("x - L", "x + L"), stringsAsFactors = FALSE)
model_wall  <- cppODE(c(x = "v", v = "-w^2 * x"), events = events_wall,
                      deriv = TRUE, deriv2 = TRUE, modelname = "saltation_wall")

pars_wall  <- c(x = 0.2, v = 1.2, w = 1, L = 0.8)
times_wall <- seq(0, 11, len = 1000)
## Two firings per wall, so 8 bounces before the grid ends.
res_wall <- solveODE(model_wall, times_wall, pars_wall, maxroot = 4L,
                     abstol = 1e-12, reltol = 1e-12, roottol = 1e-12)

x_ini <- sym("x0"); v_ini <- sym("v0")
freq  <- sp$Symbol("w", real = TRUE, positive = TRUE)
wall  <- sp$Symbol("L", real = TRUE, positive = TRUE)

point_wall <- point(list(x_ini, v_ini, freq, wall),
                    pars_wall[c("x", "v", "w", "L")], times_wall)
oscillation <- function(t0, y0, dy0)
  segment2(xfun, -pow(freq, 2) * xfun(tSym), t0, y0, dy0)

## From the initial state to whichever wall comes first.
approach   <- oscillation(sp$Integer(0L), x_ini, v_ini)
first_hit  <- crossing(list(approach - wall, approach + wall), point_wall)
first_side <- sign(value(approach$subs(tSym, first_hit), point_wall))

## The energy survives the reflection, so the speed at a wall is the same at
## every bounce and every flight between the walls is one arc up to a sign and
## a time shift. Solved once it stays small, while composing it bounce by
## bounce through dsolve() nests one inverse tangent per wall.
energy <- sp$simplify(pow(sp$diff(approach, tSym), 2) + pow(freq, 2) * pow(approach, 2))
speed  <- sp$sqrt(energy - pow(freq, 2) * pow(wall, 2))
stopifnot(isTRUE(all.equal(
  value(speed, point_wall),
  first_side * value(sp$diff(approach, tSym)$subs(tSym, first_hit), point_wall))))

## One arc, from the wall it leaves to the opposite one, in the time since
## that bounce.
flight_arc  <- oscillation(sp$Integer(0L), wall, -speed)
flight_time <- crossing(flight_arc + wall, point_wall)

branches_wall <- list(approach)
fires_wall    <- list(first_hit)
side          <- first_side
for (j in 1:8) {
  branches_wall[[j + 1]] <- sp$Integer(as.integer(side)) *
    flight_arc$subs(tSym, tSym - fires_wall[[j]])
  if (j < 8) fires_wall[[j + 1]] <- fires_wall[[j]] + flight_time
  side <- -side
}
x_wall <- composed(branches_wall, fires_wall)
v_wall <- sp$diff(x_wall, tSym)

cmp_wall_x <- comparison(x_wall, point_wall, res_wall, "x")
cmp_wall_v <- comparison(v_wall, point_wall, res_wall, "v")
report(cmp_wall_x, "5. oscillator between elastic walls, position")
report(cmp_wall_v, "5. oscillator between elastic walls, velocity")

plot_wall_1 <- comparisonPlot(cmp_wall_x, 0:1, "Position and its gradient over 8 bounces")
plot_wall_1
plot_wall_2 <- comparisonPlot(cmp_wall_x, 2L, "Second derivatives of the position")
plot_wall_2
plot_wall_3 <- comparisonPlot(cmp_wall_v, 0:1, "Velocity and its gradient over 8 bounces")
plot_wall_3
plot_wall_4 <- comparisonPlot(cmp_wall_v, 2L, "Second derivatives of the velocity")
plot_wall_4

if (!interactive()) dev.off()
