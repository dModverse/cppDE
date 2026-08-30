#!/usr/bin/env Rscript
# Regenerate dev/methods/fig-speedup.pdf, the figure of Section "Benchmarks"
# in Methods.Rmd.
#
# Reads the results.csv of one benchmark run and plots the head-to-head
# ratio t_CVODE / t_cppDE of every problem against the number of
# differentiated parameters M, once without and once with first-order
# sensitivities.  Run it from the package root after a new benchmark run,
# then re-render the vignette with dev/render-methods.R.
#
#   Rscript dev/methods/make-bench-figure.R [<results-dir>]
#
# The default is the checked-in copy of the run the vignette text quotes;
# benchmarks/results/ is generated output and not in the repository.
# dev/methods/bench/README.md records what that run was.

suppressPackageStartupMessages(library(ggplot2))

RUN <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(RUN))
  RUN <- "dev/methods/bench"
OUT <- "dev/methods/fig-speedup.pdf"

SERIES <- c("#2a78d6", "#eb6834")
LABELS <- c("without sensitivities", "with sensitivities")

gm <- function(x) exp(mean(log(x[is.finite(x) & x > 0])))

d <- read.csv(file.path(RUN, "results.csv"), stringsAsFactors = FALSE)
d <- d[d$ok & d$pinned == "auto", ]

key <- c("problem", "mode", "rtol", "nstates")
cp  <- d[d$backend == "cppde", c(key, "nsens", "time_ms")]
cv  <- d[d$backend == "cvode",  c(key, "time_ms")]
m   <- merge(cp, cv, by = key, suffixes = c(".cp", ".cv"))
m$ratio <- m$time_ms.cv / m$time_ms.cp

## One point per problem and mode: geometric mean over the tolerances.
pp <- do.call(rbind, lapply(split(m, list(m$problem, m$mode), drop = TRUE),
  function(g) data.frame(problem = g$problem[1], mode = g$mode[1],
                         nstates = g$nstates[1], M = g$nsens[1],
                         ratio = gm(g$ratio))))

## The plain solve differentiates nothing; place it at the M its
## sensitivity counterpart carries, so both series share the abscissa.
Mof <- setNames(pp$M[pp$mode == "sens1"], pp$problem[pp$mode == "sens1"])
pp$M <- unname(Mof[pp$problem])
pp <- pp[is.finite(pp$M) & pp$M > 0, ]
pp$mode <- factor(pp$mode, levels = c("nosens", "sens1"))

brk_y <- c(0.71, 1, 2, 4, 8, 16)
brk_x <- c(1, 2, 4, 8, 16, 32, 64)

p <- ggplot(pp, aes(M, ratio, colour = mode, shape = mode)) +
  geom_hline(yintercept = 1, linewidth = 0.35, colour = "grey30") +
  geom_point(size = 1.7, stroke = 0.7, alpha = 0.9) +
  scale_x_continuous(transform = "log2", breaks = brk_x, labels = brk_x,
                     expand = expansion(mult = 0.04)) +
  scale_y_continuous(transform = "log2", breaks = brk_y,
                     labels = paste0(brk_y, "×"),
                     expand = expansion(mult = 0.05)) +
  scale_colour_manual(values = SERIES, labels = LABELS, name = NULL) +
  scale_shape_manual(values = c(1, 2), labels = LABELS, name = NULL) +
  labs(x = expression("differentiated parameters" ~ italic(M)),
       y = expression(italic(t)["CVODE"] / italic(t)["cppDE"])) +
  theme_bw(base_size = 9, base_family = "serif") +
  theme(
    panel.grid       = element_blank(),
    panel.border     = element_rect(colour = "black", fill = NA,
                                    linewidth = 0.4),
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

ggsave(OUT, p, width = 6.3, height = 3.0, device = cairo_pdf)
message("Wrote ", OUT, " from ", RUN)
