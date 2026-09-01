## =====================================================================
##  plots.R -- ggplot2 figures for the benchmark results.
##
##  Five figures, each answering one question:
##    1  work-precision   how much time does a given accuracy cost?
##    2  speedup          how does cppDE compare to CVODE, per problem?
##    3  scaling          how does cost grow with system size?
##    4  sens overhead    what do first-order sensitivities cost?
##    5  summary          one number per mode, geometric mean over all cells
##
##  Every figure is also written as the CSV it was built from, so the
##  numbers can be read directly rather than measured off an axis.
## =====================================================================

for (p in c("ggplot2", "scales")) if (!requireNamespace(p, quietly = TRUE))
  stop("package '", p, "' is required for the benchmark plots")

library(ggplot2)

## Categorical palette: slots 1-4 of the validated default order
## (blue / orange / aqua / violet).  Verified for all-pairs colour-vision
## separation, so the same colours work in scatter and bar alike.  Shape
## is carried alongside colour throughout as a secondary encoding.
BENCH_COLS <- c(
  cppDE_ndf = "#2a78d6",
  CVODE_bdf  = "#eb6834",
  cppDE_bdf = "#1baf7a",
  cppDE_rb4 = "#4a3aa7")
BENCH_SHAPES <- c(cppDE_ndf = 16, CVODE_bdf = 17, cppDE_bdf = 15, cppDE_rb4 = 18)

## Keyed by backend, for the sparse-sweep figures.
BACKEND_COLS <- c(cppde = "#2a78d6", cvode = "#eb6834")
BACKEND_LABELS <- c(cppde = "cppDE", cvode = "CVODE")

MODE_LABELS <- c(nosens = "without sensitivities",
                 sens1  = "with 1st-order sensitivities",
                 sens2  = "with 2nd-order sensitivities")

theme_bench <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.3, colour = "grey88"),
      panel.border     = element_rect(fill = NA, colour = "grey85", linewidth = 0.4),
      axis.title       = element_text(colour = "grey20"),
      axis.text        = element_text(colour = "grey35"),
      strip.text       = element_text(face = "bold", size = base_size - 1,
                                      colour = "grey15"),
      strip.background = element_rect(fill = "grey96", colour = NA),
      plot.title       = element_text(face = "bold", size = base_size + 3),
      plot.subtitle    = element_text(colour = "grey35", size = base_size - 1),
      plot.caption     = element_text(colour = "grey45", size = base_size - 2,
                                      hjust = 0),
      legend.position  = "top",
      legend.title     = element_blank(),
      legend.key.height = unit(10, "pt"))
}

scale_solver <- function(...) list(
  scale_colour_manual(values = BENCH_COLS, ...),
  scale_shape_manual(values = BENCH_SHAPES, ...))

label_mode <- function(x) unname(MODE_LABELS[x])

## Modes as a factor ordered by MODE_LABELS; an unlabelled mode keeps its
## own name instead of becoming NA.
mode_factor <- function(x) {
  lab <- unname(MODE_LABELS[x])
  lab[is.na(lab)] <- x[is.na(lab)]
  factor(lab, levels = intersect(c(unname(MODE_LABELS), unique(lab)), unique(lab)))
}

## Log ticks that stay readable when the range spans many decades.
## Counts (states, parameters) get plain integer breaks -- a label such
## as 10^0.48 on an axis titled "number of states" is unreadable.
log_x <- scale_x_log10(labels = scales::label_log(digits = 2))
log_y <- scale_y_log10(labels = scales::label_number(drop0trailing = TRUE))

count_breaks <- function(x) {
  rng <- range(x[is.finite(x) & x > 0])
  cand <- c(1, 2, 3, 5, 8, 12, 20, 30, 50, 80, 125, 200, 300, 500, 1000, 2000)
  b <- cand[cand >= rng[1] * 0.95 & cand <= rng[2] * 1.05]
  if (length(b) < 2L) b <- unique(round(exp(seq(log(rng[1]), log(rng[2]),
                                               length.out = 4))))
  b
}
log_x_count <- scale_x_log10(breaks = count_breaks,
                             labels = scales::label_number(accuracy = 1))


## ---------------------------------------------------------------------
##  1  Work-precision
## ---------------------------------------------------------------------

## The standard ODE-solver diagram: achieved accuracy on x, cost on y,
## one point per tolerance.  Down-and-left is better; a curve lying below
## another is uniformly cheaper at equal accuracy.
plot_work_precision <- function(df, mode = "nosens", max_panels = 12L) {
  d <- df[df$ok & df$mode == mode & is.finite(df$err) & df$err > 0 &
          is.finite(df$time_ms), ]
  if (!nrow(d)) return(NULL)
  d$panel <- if (all(is.na(d$condition))) d$problem else
    ifelse(is.na(d$condition), d$problem, paste0(d$problem, "\n", d$condition))
  ## Keep the panels with the most complete solver coverage.
  keep <- names(sort(tapply(d$solver, d$panel, function(x) length(unique(x))),
                     decreasing = TRUE))[seq_len(min(max_panels,
                       length(unique(d$panel))))]
  d <- d[d$panel %in% keep, ]
  d$panel <- factor(d$panel, levels = keep)

  ggplot(d, aes(err, time_ms, colour = solver, shape = solver)) +
    geom_line(linewidth = 0.6, alpha = 0.85) +
    geom_point(size = 2.1) +
    facet_wrap(~ panel, scales = "free", ncol = 4) +
    scale_solver() + log_x + log_y +
    labs(title = "Work-precision diagram",
         subtitle = paste0(label_mode(mode),
                           ", lower and further left is better"),
         x = "achieved error (scaled max-norm vs reference solution)",
         y = "wall-clock time per solve [ms]",
         caption = paste("Each point is one (absolute, relative) tolerance pair.",
                         "Reference: CVODE at atol 1e-14 / rtol 1e-12.")) +
    theme_bench()
}


## ---------------------------------------------------------------------
##  2  Speedup against CVODE
## ---------------------------------------------------------------------

## One bar per problem: the CVODE time divided by the solver's time, at
## matched tolerance, aggregated over tolerances with a geometric mean.
plot_speedup <- function(df, baseline = "CVODE_bdf") {
  s <- speedup_table(df, baseline)
  if (is.null(s) || !nrow(s)) return(NULL)
  s <- s[s$solver != baseline, ]
  if (!nrow(s)) return(NULL)

  agg <- stats::aggregate(speedup ~ problem + mode + solver + nstates,
                          data = s, FUN = geo_mean)
  agg <- agg[is.finite(agg$speedup) & agg$speedup > 0, ]
  if (!nrow(agg)) return(NULL)
  ## coord_flip() puts the first factor level at the bottom, so reverse
  ## the size ordering to read smallest-at-top.
  ord <- unique(agg[order(agg$nstates, agg$problem), c("problem", "nstates")])
  agg$problem <- factor(agg$problem, levels = rev(ord$problem))
  agg$mode <- factor(label_mode(as.character(agg$mode)),
                     levels = unname(MODE_LABELS))
  ## Bars encode log2(ratio) so they grow from parity rather than from
  ## the bottom of a log axis, which would exaggerate every difference.
  agg$lr <- log2(agg$speedup)
  brk <- c(0.25, 0.5, 0.71, 1, 1.4, 2, 4, 8)
  brk <- brk[log2(brk) >= min(agg$lr) - 0.3 & log2(brk) <= max(agg$lr) + 0.3]

  ggplot(agg, aes(problem, lr, fill = solver)) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.5) +
    geom_col(position = position_dodge(width = 0.75), width = 0.68) +
    coord_flip() +
    facet_wrap(~ mode, nrow = 1) +
    scale_fill_manual(values = BENCH_COLS) +
    scale_y_continuous(breaks = log2(brk),
                       labels = sprintf("%.2g×", brk)) +
    labs(title = sprintf("Speed-up relative to %s", baseline),
         subtitle = "geometric mean over all tolerances; problems ordered by system size (smallest first)",
         x = NULL, y = sprintf("%s time / solver time  (right of the line = faster)", baseline),
         caption = "Identical right-hand side, analytic Jacobian, output grid and tolerances on both sides.") +
    theme_bench()
}


## ---------------------------------------------------------------------
##  3  Scaling with system size
## ---------------------------------------------------------------------

plot_scaling <- function(df) {
  d <- df[df$ok & is.finite(df$time_ms), ]
  if (!nrow(d)) return(NULL)
  agg <- stats::aggregate(time_ms ~ problem + nstates + solver + mode,
                          data = d, FUN = stats::median)
  agg$mode <- factor(label_mode(as.character(agg$mode)),
                     levels = unname(MODE_LABELS))

  ggplot(agg, aes(nstates, time_ms, colour = solver, shape = solver)) +
    geom_point(size = 2.2, alpha = 0.9) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                linewidth = 0.6, linetype = "22") +
    facet_wrap(~ mode, nrow = 1) +
    scale_solver() + log_x_count + log_y +
    labs(title = "Cost versus system size",
         subtitle = "median over tolerances; dashed lines are log-log fits",
         x = "number of states", y = "wall-clock time per solve [ms]",
         caption = "A steeper slope means worse asymptotic scaling in the number of states.") +
    theme_bench()
}


## ---------------------------------------------------------------------
##  4  Cost of first-order sensitivities
## ---------------------------------------------------------------------

## The quantity that matters for gradient-based fitting: how many plain
## solves does one solve-with-gradient cost?  The ideal for forward
## sensitivities is well below the finite-difference cost of M+1.
plot_sens_overhead <- function(df) {
  d <- df[df$ok & is.finite(df$time_ms), ]
  if (!nrow(d)) return(NULL)
  key <- c("problem", "condition", "solver", "atol", "rtol")
  a <- d[d$mode == "nosens", c(key, "time_ms")]
  b <- d[d$mode == "sens1",  c(key, "time_ms", "nsens")]
  names(a)[names(a) == "time_ms"] <- "t_plain"
  names(b)[names(b) == "time_ms"] <- "t_sens"
  m <- merge(a, b, by = key)
  m <- m[is.finite(m$t_plain) & m$t_plain > 0 & m$nsens > 0, ]
  if (!nrow(m)) return(NULL)
  m$overhead <- m$t_sens / m$t_plain
  agg <- stats::aggregate(overhead ~ problem + solver + nsens, data = m,
                          FUN = geo_mean)

  ggplot(agg, aes(nsens, overhead, colour = solver, shape = solver)) +
    geom_abline(slope = 1, intercept = 1, colour = "grey55",
                linetype = "22", linewidth = 0.5) +
    geom_point(size = 2.3, alpha = 0.9) +
    scale_solver() + log_x_count + log_y +
    labs(title = "Cost of first-order sensitivities",
         subtitle = "solve-with-gradient divided by plain solve, per problem",
         x = "number of sensitivity parameters M",
         y = "relative cost  (plain solve = 1)",
         caption = paste("Grey line: the M+1 cost of one-sided finite differences.",
                         "Points below it beat differencing.")) +
    theme_bench()
}


## ---------------------------------------------------------------------
##  5b  Second-order sensitivities
## ---------------------------------------------------------------------

## cppDE-only: CVODES provides no second-order sensitivities, so this is
## a cost curve, not a comparison.  What it shows is how the price of a
## full Hessian grows with the number of parameters -- forward-over-
## forward AD carries M(M+1)/2 second-order directions, so the
## expectation is quadratic.
plot_sens2_cost <- function(df) {
  d <- df[df$ok & is.finite(df$time_ms) & df$mode == "sens2", ]
  if (!nrow(d)) return(NULL)
  key <- c("problem", "condition", "solver", "atol", "rtol")
  base <- df[df$ok & df$mode == "nosens", c(key, "time_ms")]
  names(base)[names(base) == "time_ms"] <- "t_plain"
  m <- merge(d, base, by = key, all.x = TRUE)
  m$rel <- m$time_ms / m$t_plain
  agg <- stats::aggregate(cbind(time_ms, rel) ~ problem + nsens + nstates,
                          data = m[is.finite(m$rel), , drop = FALSE], FUN = geo_mean)
  if (!nrow(agg)) return(NULL)

  ggplot(agg, aes(nsens, rel)) +
    geom_point(colour = BENCH_COLS[["cppDE_ndf"]], size = 2.4, alpha = 0.9) +
    geom_text(aes(label = problem), hjust = -0.12, size = 2.7, colour = "grey35") +
    scale_x_log10(breaks = count_breaks,
                  labels = scales::label_number(accuracy = 1),
                  expand = expansion(mult = c(0.08, 0.35))) + log_y +
    labs(title = "Cost of second-order sensitivities",
         subtitle = "cppDE only -- CVODES has no second-order sensitivities",
         x = "number of sensitivity parameters M",
         y = "cost relative to a plain solve",
         caption = paste("Forward-over-forward AD carries M(M+1)/2 second-order",
                         "directions, so the expected growth is quadratic in M.")) +
    theme_bench()
}


## ---------------------------------------------------------------------
##  5  Headline summary
## ---------------------------------------------------------------------

plot_summary <- function(df, baseline = "CVODE_bdf") {
  s <- speedup_table(df, baseline)
  if (is.null(s) || !nrow(s)) return(NULL)
  s <- s[s$solver != baseline, ]
  if (!nrow(s)) return(NULL)
  agg <- stats::aggregate(cbind(speedup, fev_ratio) ~ solver + mode,
                          data = s, FUN = geo_mean)
  n <- stats::aggregate(speedup ~ solver + mode, data = s, FUN = length)
  names(n)[3] <- "n"
  agg <- merge(agg, n, by = c("solver", "mode"))
  agg$mode <- factor(label_mode(as.character(agg$mode)),
                     levels = unname(MODE_LABELS))

  ggplot(agg, aes(solver, speedup, fill = solver)) +
    geom_hline(yintercept = 1, colour = "grey40", linewidth = 0.5) +
    geom_col(width = 0.55) +
    geom_text(aes(label = sprintf("%.2f×  (n=%d)", speedup, n)),
              vjust = -0.6, size = 3.4, colour = "grey20") +
    facet_wrap(~ mode, nrow = 1) +
    scale_fill_manual(values = BENCH_COLS) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(title = sprintf("Overall speed-up versus %s", baseline),
         subtitle = "geometric mean over every matched problem / condition / tolerance cell",
         x = NULL, y = sprintf("%s time / solver time", baseline)) +
    theme_bench()
}


## ---------------------------------------------------------------------
##  8  Dense against sparse
## ---------------------------------------------------------------------

## One bar per problem and backend: dense time over sparse time at matched
## tolerance, aggregated over tolerances with a geometric mean.
plot_sparse_gain <- function(df) {
  s <- sparse_gain_table(df)
  if (is.null(s) || !nrow(s)) return(NULL)
  agg <- stats::aggregate(gain ~ problem + mode + backend + nstates,
                          data = s, FUN = geo_mean)
  agg <- agg[is.finite(agg$gain) & agg$gain > 0, ]
  if (!nrow(agg)) return(NULL)

  ## Reversed because coord_flip() puts the first level at the bottom.
  ord <- unique(agg[order(agg$nstates, agg$problem), c("problem", "nstates")])
  agg$problem <- factor(sprintf("%s (%d)", agg$problem, agg$nstates),
                        levels = rev(sprintf("%s (%d)", ord$problem, ord$nstates)))
  agg$mode <- mode_factor(agg$mode)
  ## Bars encode log2(ratio).
  agg$lr <- log2(agg$gain)
  ## Parity is always among the breaks.
  brk <- c(0.25, 0.5, 0.71, 1, 1.4, 2, 3, 4, 6, 8)
  brk <- sort(unique(c(1, brk[log2(brk) >= min(agg$lr) - 0.3 &
                             log2(brk) <= max(agg$lr) + 0.3])))

  ggplot(agg, aes(problem, lr, fill = backend)) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.5) +
    geom_col(position = position_dodge(width = 0.75), width = 0.68) +
    coord_flip() +
    facet_wrap(~ mode, nrow = 1) +
    scale_fill_manual(values = BACKEND_COLS, labels = BACKEND_LABELS) +
    scale_y_continuous(breaks = log2(brk), labels = sprintf("%.2g×", brk)) +
    labs(title = "Dense against sparse linear solver",
         subtitle = paste("Right of parity the sparse path was faster.",
                          "Problem size in states in brackets."),
         x = NULL, y = "dense time / sparse time", fill = NULL) +
    theme_bench()
}


## The same gains against system size rather than by name.
plot_sparse_crossover <- function(df) {
  s <- sparse_gain_table(df)
  if (is.null(s) || !nrow(s)) return(NULL)
  agg <- stats::aggregate(gain ~ problem + mode + backend + nstates,
                          data = s, FUN = geo_mean)
  agg <- agg[is.finite(agg$gain) & agg$gain > 0, ]
  if (!nrow(agg)) return(NULL)
  agg$mode <- mode_factor(agg$mode)
  agg$lr <- log2(agg$gain)
  ## Parity is always among the breaks.
  brk <- c(0.25, 0.5, 0.71, 1, 1.4, 2, 3, 4, 6, 8)
  brk <- sort(unique(c(1, brk[log2(brk) >= min(agg$lr) - 0.3 &
                             log2(brk) <= max(agg$lr) + 0.3])))

  ggplot(agg, aes(nstates, lr, colour = backend, shape = backend)) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.5) +
    geom_point(size = 2.4, alpha = 0.9) +
    facet_wrap(~ mode, nrow = 1) +
    scale_colour_manual(values = BACKEND_COLS, labels = BACKEND_LABELS) +
    scale_shape_manual(values = c(cppde = 16, cvode = 17),
                       labels = BACKEND_LABELS) +
    log_x_count +
    scale_y_continuous(breaks = log2(brk), labels = sprintf("%.2g×", brk)) +
    labs(title = "Where the sparse path starts to pay",
         subtitle = "one point per model; above parity the sparse path was faster",
         x = "number of states", y = "dense time / sparse time",
         colour = NULL, shape = NULL) +
    theme_bench()
}


## ---------------------------------------------------------------------
##  Write everything
## ---------------------------------------------------------------------

save_plots <- function(df, outdir, device = "png", width = 11, height = 7,
                       dpi = 160) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  n_prob <- length(unique(df$problem))
  ## The head-to-head figures read the auto-detected rows only, so they
  ## are drawn from that subset rather than from everything the sweep
  ## added alongside it.
  auto <- df[pinned_col(df) == "auto", , drop = FALSE]
  n_sweep <- length(unique(sparse_gain_table(df)$problem))
  ## A results.csv from a pinned-only run (the sweep as its own run, as it
  ## used to be) has no auto rows at all; those figures are then skipped
  ## rather than built from nothing.
  on_auto <- function(f, ...) if (nrow(auto)) f(auto, ...) else NULL
  figs <- list(
    `01-work-precision-nosens` = list(p = on_auto(plot_work_precision, "nosens"),
                                      h = height + 1),
    `02-work-precision-sens1`  = list(p = on_auto(plot_work_precision, "sens1"),
                                      h = height + 1),
    `03-speedup`               = list(p = on_auto(plot_speedup),
                                      h = max(4, 0.28 * n_prob + 2.5)),
    `04-scaling`               = list(p = on_auto(plot_scaling), h = 5),
    `05-sens-overhead`         = list(p = on_auto(plot_sens_overhead), h = 5),
    `06-summary`               = list(p = on_auto(plot_summary), h = 4.5),
    `07-sens2-cost`            = list(p = on_auto(plot_sens2_cost), h = 5),
    `08-sparse-gain`           = list(p = plot_sparse_gain(df),
                                      h = max(4, 0.30 * n_sweep + 2.5)),
    `09-sparse-crossover`      = list(p = plot_sparse_crossover(df), h = 5))
  written <- character(0)
  for (nm in names(figs)) {
    f <- figs[[nm]]
    if (is.null(f$p)) { message("  [skip] ", nm, " -- no data"); next }
    path <- file.path(outdir, paste0(nm, ".", device))
    ggsave(path, f$p, width = width, height = f$h, dpi = dpi,
           device = device, bg = "white", limitsize = FALSE)
    written <- c(written, path)
    message("  wrote ", path)
  }
  written
}
