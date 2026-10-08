tF <- seq(0, 5, by = 0.5)
forcings <- list(F1 = data.frame(time = tF, value = sin(tF)),
                 F2 = data.frame(time = 0, value = 2))
times <- seq(-1, 7, by = 0.25)
v <- forcingValues(times, forcings)
matplot(times, v, type = "l", lty = 1, xlab = "time", ylab = "forcing")
points(tF, sin(tF))
