#!/usr/bin/env Rscript
# Regenerate the figures of Section "Benchmarks" in Methods.Rmd:
#   fig-speedup.pdf   head-to-head ratio t_CVODE / t_cppDE against M
#   fig-gradient.pdf  cost of one gradient by forward sensitivities against M
#   fig-hessian.pdf   cost of one Hessian, forward-forward and forward-reverse
#   fig-sparse.pdf    gain of the sparse linear solver against n_x
#   fig-adjoint.pdf   CVODES adjoint against the discrete adjoint (adjoint.csv)
#
# Reads the results.csv of one benchmark run. Run it from the package root after
# a new run, then re-render the vignette with dev/render-methods.R.
#
#   Rscript dev/methods/make-bench-figure.R [<results-dir>]
#
# The default is the checked-in copy of the run the vignette text quotes;
# dev/methods/bench/README.md records what that run was.

suppressPackageStartupMessages(library(ggplot2))

RUN <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(RUN))
  RUN <- "dev/methods/bench"
DIR <- "dev/methods"

BLUE <- "#2a78d6"; ORANGE <- "#eb6834"; GREEN <- "#1baf7a"

gm <- function(x) exp(mean(log(x[is.finite(x) & x > 0])))

theme_fig <- function()
  theme_bw(base_size = 9, base_family = "serif") +
  theme(
    panel.grid       = element_blank(),
    panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.4),
    axis.ticks       = element_line(colour = "black", linewidth = 0.3),
    axis.ticks.length = unit(2.5, "pt"),
    axis.text        = element_text(colour = "black"),
    legend.position  = "inside",
    legend.position.inside = c(0.015, 0.985),
    legend.justification   = c(0, 1),
    legend.background = element_blank(),
    legend.key        = element_blank(),
    legend.key.spacing.y = unit(0, "pt"),
    legend.margin     = margin(0, 0, 0, 0),
    plot.margin       = margin(2, 4, 2, 2))

log2_axis <- function(breaks, suffix = "", ...)
  list(transform = "log2", breaks = breaks, labels = paste0(breaks, suffix), ...)

save <- function(p, name, width = 6.3, height = 3.0) {
  ggsave(file.path(DIR, name), p, width = width, height = height, device = cairo_pdf)
  message("Wrote ", file.path(DIR, name), " from ", RUN)
}

d <- read.csv(file.path(RUN, "results.csv"), stringsAsFactors = FALSE)
if (!"deriv" %in% names(d))
  d$deriv <- ifelse(d$mode == "nosens", "none", "forward")
ok   <- d[d$ok, ]
auto <- ok[ok$pinned == "auto", ]

## Plain solves per problem, backend and tolerance: the unit the derivative
## costs are expressed in.
plain <- auto[auto$mode == "nosens", c("problem", "backend", "rtol", "time_ms")]
names(plain)[4] <- "t_plain"


## ---------------------------------------------------------------------
##  Head-to-head, without and with forward sensitivities
## ---------------------------------------------------------------------

LABELS <- c(nosens = "without sensitivities", sens1 = "with sensitivities")
hh  <- auto[auto$mode %in% c("nosens", "sens1") & auto$deriv %in% c("none", "forward"), ]
key <- c("problem", "mode", "rtol", "nstates")
m   <- merge(hh[hh$backend == "cppde", c(key, "nsens", "time_ms")],
             hh[hh$backend == "cvode", c(key, "time_ms")],
             by = key, suffixes = c(".cp", ".cv"))
m$ratio <- m$time_ms.cv / m$time_ms.cp

## One point per problem and mode: geometric mean over the tolerances. The
## plain solve sits at the M of its sensitivity counterpart.
pp <- do.call(rbind, lapply(split(m, list(m$problem, m$mode), drop = TRUE),
  function(g) data.frame(problem = g$problem[1], mode = g$mode[1],
                         M = g$nsens[1], ratio = gm(g$ratio))))
Mof  <- setNames(pp$M[pp$mode == "sens1"], pp$problem[pp$mode == "sens1"])
pp$M <- unname(Mof[pp$problem])
pp   <- pp[is.finite(pp$M) & pp$M > 0, ]
pp$mode <- factor(pp$mode, levels = names(LABELS))

p <- ggplot(pp, aes(M, ratio, colour = mode, shape = mode)) +
  geom_hline(yintercept = 1, linewidth = 0.35, colour = "grey30") +
  geom_point(size = 1.7, stroke = 0.7, alpha = 0.9) +
  do.call(scale_x_continuous, log2_axis(c(1, 2, 4, 8, 16, 32),
                                        expand = expansion(mult = 0.04))) +
  do.call(scale_y_continuous, log2_axis(c(0.71, 1, 1.41, 2, 2.83), "×",
                                        expand = expansion(mult = 0.05))) +
  scale_colour_manual(values = c(nosens = BLUE, sens1 = ORANGE), labels = LABELS,
                      name = NULL) +
  scale_shape_manual(values = c(nosens = 1, sens1 = 2), labels = LABELS, name = NULL) +
  labs(x = expression("differentiated parameters" ~ italic(M)),
       y = expression(italic(t)["CVODE"] / italic(t)["cppDE"])) +
  theme_fig()
save(p, "fig-speedup.pdf")


## ---------------------------------------------------------------------
##  One gradient by forward sensitivities
## ---------------------------------------------------------------------

g <- merge(auto[auto$mode == "sens1" & auto$deriv == "forward",
                c("problem", "backend", "rtol", "nsens", "time_ms")],
           plain, by = c("problem", "backend", "rtol"))
g$cost <- g$time_ms / g$t_plain
g <- aggregate(cost ~ problem + backend + nsens, data = g, FUN = gm)
GLAB <- c(cppde = "cppDE", cvode = "CVODES")

## The forward cost is linear in M; one least-squares line per backend.
fit_line <- do.call(rbind, lapply(names(GLAB), function(b) {
  f <- lm(cost ~ nsens, data = g[g$backend == b, ])
  M <- exp(seq(log(1), log(max(g$nsens)), length.out = 200))
  data.frame(backend = b, M = M, cost = predict(f, data.frame(nsens = M)))
}))

p <- ggplot(g, aes(nsens, cost, colour = backend)) +
  geom_line(data = fit_line, aes(M, cost), linewidth = 0.45, show.legend = FALSE) +
  geom_point(aes(shape = backend), size = 1.7, stroke = 0.7, alpha = 0.9) +
  do.call(scale_x_continuous, log2_axis(c(1, 2, 4, 8, 16, 32),
                                        expand = expansion(mult = 0.04))) +
  do.call(scale_y_continuous, log2_axis(c(1, 2, 4, 8, 16, 32), "×",
                                        expand = expansion(mult = 0.05))) +
  scale_colour_manual(values = c(cppde = BLUE, cvode = ORANGE), labels = GLAB,
                      name = NULL) +
  scale_shape_manual(values = c(cppde = 1, cvode = 2), labels = GLAB, name = NULL) +
  labs(x = expression("differentiated parameters" ~ italic(M)),
       y = "gradient cost / plain solve") +
  theme_fig()
save(p, "fig-gradient.pdf")


## ---------------------------------------------------------------------
##  One Hessian: forward-forward against forward-reverse
## ---------------------------------------------------------------------

h <- merge(ok[ok$mode == "sens2", c("problem", "rtol", "deriv", "nsens", "time_ms")],
           plain[plain$backend == "cppde", c("problem", "rtol", "t_plain")],
           by = c("problem", "rtol"))
h$cost <- h$time_ms / h$t_plain
h <- aggregate(cost ~ problem + deriv + nsens, data = h, FUN = gm)
HLAB <- c(`forward-forward` = "forward-forward", `forward-reverse` = "forward-reverse")
h$deriv <- factor(h$deriv, levels = names(HLAB))

## Power laws fitted in log-log over M > 1, one per mode.
pl_line <- do.call(rbind, lapply(names(HLAB), function(dm) {
  s <- h[h$deriv == dm & h$nsens > 1, ]
  f <- lm(log(cost) ~ log(nsens), data = s)
  M <- exp(seq(log(1), log(max(s$nsens)), length.out = 100))
  data.frame(deriv = dm, M = M, cost = exp(predict(f, data.frame(nsens = M))),
             slope = coef(f)[2])
}))
pl_line$deriv <- factor(pl_line$deriv, levels = names(HLAB))
slopes <- tapply(pl_line$slope, pl_line$deriv, `[`, 1)
HLAB2 <- sprintf("%s, slope %.1f", HLAB, slopes[names(HLAB)])
names(HLAB2) <- names(HLAB)

p <- ggplot(h, aes(nsens, cost, colour = deriv)) +
  geom_line(data = pl_line, aes(M, cost), linewidth = 0.45, show.legend = FALSE) +
  geom_point(aes(shape = deriv), size = 1.7, stroke = 0.7, alpha = 0.9) +
  do.call(scale_x_continuous, log2_axis(c(1, 2, 4, 8, 16, 32),
                                        expand = expansion(mult = 0.04))) +
  do.call(scale_y_continuous, log2_axis(c(2, 4, 8, 16, 32, 64, 128, 256), "×",
                                        expand = expansion(mult = 0.05))) +
  scale_colour_manual(values = c(`forward-forward` = BLUE, `forward-reverse` = ORANGE),
                      labels = HLAB2, name = NULL) +
  scale_shape_manual(values = c(`forward-forward` = 1, `forward-reverse` = 2),
                     labels = HLAB2, name = NULL) +
  labs(x = expression("differentiated parameters" ~ italic(M)),
       y = "Hessian cost / plain solve") +
  theme_fig()
save(p, "fig-hessian.pdf")


## ---------------------------------------------------------------------
##  Dense against sparse LU, per backend
## ---------------------------------------------------------------------

sw <- ok[ok$pinned != "auto" & ok$mode %in% c("nosens", "sens1") &
         ok$deriv %in% c("none", "forward"), ]
sg <- merge(sw[sw$pinned == "dense",  c("problem", "backend", "mode", "rtol", "nstates", "time_ms")],
            sw[sw$pinned == "sparse", c("problem", "backend", "mode", "rtol", "time_ms")],
            by = c("problem", "backend", "mode", "rtol"), suffixes = c(".d", ".s"))
sg$gain <- sg$time_ms.d / sg$time_ms.s
sg <- aggregate(gain ~ problem + backend + nstates, data = sg, FUN = gm)
SLAB <- c(cppde = "cppDE", cvode = "CVODE")

p <- ggplot(sg, aes(nstates, gain, colour = backend, shape = backend)) +
  geom_hline(yintercept = 1, linewidth = 0.35, colour = "grey30") +
  geom_point(size = 1.7, stroke = 0.7, alpha = 0.9) +
  do.call(scale_x_continuous, log2_axis(c(8, 16, 32, 64, 128, 256, 512),
                                        expand = expansion(mult = 0.05))) +
  do.call(scale_y_continuous, log2_axis(c(1, 1.41, 2, 2.83, 4), "×",
                                        expand = expansion(mult = 0.05))) +
  scale_colour_manual(values = c(cppde = BLUE, cvode = ORANGE), labels = SLAB,
                      name = NULL) +
  scale_shape_manual(values = c(cppde = 1, cvode = 2), labels = SLAB, name = NULL) +
  labs(x = expression("states" ~ italic(n)[x]),
       y = expression(italic(t)["dense"] / italic(t)["sparse"])) +
  theme_fig()
save(p, "fig-sparse.pdf", width = 3.1, height = 2.3)


## ---------------------------------------------------------------------
##  One gradient in reverse: CVODES adjoint against the discrete adjoint
## ---------------------------------------------------------------------

## A separate run over the models with M >= 40, both backends in reverse.
adj_file <- file.path(RUN, "adjoint.csv")
if (file.exists(adj_file)) {
  a  <- read.csv(adj_file, stringsAsFactors = FALSE)
  a  <- a[a$ok, ]
  aw <- merge(a[a$backend == "cppde", c("problem", "rtol", "nsens", "nstates", "time_ms")],
              a[a$backend == "cvode", c("problem", "rtol", "time_ms")],
              by = c("problem", "rtol"), suffixes = c(".cp", ".cv"))
  aw$ratio <- aw$time_ms.cv / aw$time_ms.cp
  lab <- unique(a[a$backend == "cppde" & a$problem %in% aw$problem,
                  c("problem", "nsens", "nstates")])
  lab$label <- sprintf("%s  (M = %d, n = %d)", sub("_.*", "", lab$problem),
                       lab$nsens, lab$nstates)
  lab <- lab[order(lab$nsens), ]
  aw$label <- factor(lab$label[match(aw$problem, lab$problem)], levels = lab$label)
  aw$rtol  <- factor(format(aw$rtol, scientific = TRUE),
                     levels = format(sort(unique(aw$rtol), decreasing = TRUE),
                                     scientific = TRUE))
  RCOL <- c("#9ec5f3", "#2a78d6", "#0d3b73")
  names(RCOL) <- levels(aw$rtol)

  p <- ggplot(aw, aes(ratio, label, colour = rtol, shape = rtol)) +
    geom_vline(xintercept = 1, linewidth = 0.35, colour = "grey30") +
    geom_point(size = 2, stroke = 0.8) +
    do.call(scale_x_continuous, log2_axis(c(0.5, 1, 2, 4, 8, 16, 32), "×",
                                          expand = expansion(mult = 0.05))) +
    scale_colour_manual(values = RCOL, name = "rtol") +
    scale_shape_manual(values = c(1, 2, 0), name = "rtol") +
    labs(x = expression(italic(t)["CVODES adjoint"] / italic(t)["discrete adjoint"]),
         y = NULL) +
    theme_fig() +
    theme(legend.position.inside = c(0.3, 0.95), legend.justification = c(0, 1))
  save(p, "fig-adjoint.pdf", height = 2.4)
}
