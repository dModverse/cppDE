## ---- longFormat
# states and sensitivities of a solveODE() result in long format, one panel per quantity
longFormat <- function(out, states = colnames(out$variable), sens = TRUE) {
  nt <- length(out$time)
  value <- data.frame(time = out$time, panel = rep(states, each = nt),
                      value = as.vector(out$variable[, states]))
  if (!sens) return(value)
  tangent <- out$tangent[, states, , drop = FALSE]
  grid <- expand.grid(time = out$time, state = states,
                      par = dimnames(tangent)[[3]], stringsAsFactors = FALSE)
  panels <- c(states, sprintf("∂%s/∂%s", grid$state, grid$par))
  long <- rbind(value, data.frame(time = grid$time,
                                  panel = sprintf("∂%s/∂%s", grid$state, grid$par),
                                  value = as.vector(tangent)))
  long$panel <- factor(long$panel, levels = unique(panels))
  long
}

## ---- style
# the plot style of the dMod2 site: framed panels without grid, ticks inside, white strips, and
# the Okabe-Ito colors for discrete scales
ggplot2::theme_set(
  ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      line = ggplot2::element_line(color = "black"),
      rect = ggplot2::element_rect(fill = "white", color = NA),
      text = ggplot2::element_text(color = "black"),
      axis.text = ggplot2::element_text(color = "black"),
      axis.text.x = ggplot2::element_text(margin = ggplot2::margin(4, 4, 0, 4, "mm")),
      axis.text.y = ggplot2::element_text(margin = ggplot2::margin(4, 4, 4, 0, "mm")),
      axis.ticks = ggplot2::element_line(color = "black"),
      axis.ticks.length = grid::unit(-2, "mm"),
      legend.key = ggplot2::element_rect(color = NA),
      panel.border = ggplot2::element_rect(color = "black"),
      panel.grid = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "white", color = NA),
      strip.text = ggplot2::element_text(size = ggplot2::rel(1))))
okabeIto <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00",
              "#CC79A7")
options(ggplot2.discrete.colour = okabeIto, ggplot2.discrete.fill = okabeIto)
