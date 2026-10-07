\donttest{
if (codegenAvailable()) {
  ## A fast forcing the state does not see, but its derivative in q does.
  rev <- cppODE(c(x = "-k*x + sin(100*time)*q"), derivMode = "reverse",
                modelname = "adjointControl_example")
  times <- seq(0, 2, 0.1)
  parms <- c(x = 1, k = 1, q = 0)
  W <- cbind(x = rep(1, length(times)))
  ref <- solveODE(rev, times, parms, cotangent = W, abstol = 1e-12,
                  reltol = 1e-12)$cotangent[, 1]

  ## The sweep on the grid of the value run, and checked against the
  ## tolerances, with the gradient held to gradtol as well.
  val <- solveODE(rev, times, parms, cotangent = W, reltol = 1e-6)
  chk <- solveODE(rev, times, parms, cotangent = W, reltol = 1e-6,
                  adjoint = adjointControl(refine = TRUE, gradtol = 1e-9,
                                           trace = TRUE))
  print(rbind(value_grid = val$cotangent[, 1] - ref, checked = chk$cotangent[, 1] - ref))
  sum(chk$adjoint$substeps)
}
}
