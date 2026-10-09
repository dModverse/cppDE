\donttest{
if (codegenAvailable()) {
  ## Comparisons and logic give 1 or 0: as a factor or as a condition.
  ##   u: a pulse between t1 and t2, `&&` joins two comparisons on time
  ##   x: decays fast below xc and slowly above, unless mode is 0
  ##   y: produced once x exceeds xc, with powers as `^` and `**`
  eqns <- c(
    u = "a*(time >= t1 && time < t2) - u/tau",
    x = "k*u - piecewise(kfast, mode == 0 || !(x > xc), kslow)*x",
    y = "kp*Heaviside(x - xc)*(1 - y^2) - y**2/(K + y)"
  )
  parms <- c(u = 0, x = 0, y = 0, a = 2, t1 = 1, t2 = 3, tau = 1, k = 1,
             kslow = 0.1, kfast = 1, xc = 0.5, mode = 1, kp = 1, K = 0.5)
  times <- 0:8

  ## The same equations on both backends, compiled together
  models <- list(cppODE = cppODE(eqns, deriv = FALSE, compile = FALSE))
  if (codegenAvailable("cvode"))
    models$cvode <- cvode(eqns, compile = FALSE)
  do.call(compile, c(unname(models), list(output = "expressions_example")))

  out <- lapply(models, function(m) solveODE(m, times, parms)$variable)
  cbind(time = times, out$cppODE)
  if (!is.null(out$cvode)) max(abs(out$cppODE - out$cvode))
}
}
