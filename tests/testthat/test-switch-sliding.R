# Sliding along switching surfaces (Filippov): where the branches of switches
# push the solution onto their surfaces g = 0, it stays there with f the mix
# that keeps every g' at 0, and leaves where one branch turns away. Values and
# sensitivities against closed forms. The blocks of one model slide at
# overlapping times, so the mix of several slides is exercised throughout.

skip_on_cran()

library(cppDE)

sl_rhs <- c(
  x = "kin - v*(x > 0)", y = "v*(x > 0)",
  e = "c*time - ve*(e > 0)",
  A = "k*(A*A < Acrit) - d*A",
  cx = "u - vc*(cx > 0) - w*(cz > 0)", cz = "s - q*(cz > 0)")
sl_p <- c(x = 1, y = 0, kin = 0.2, v = 0.8,
          e = 0.5, c = 0.5, ve = 1,
          A = 0, k = 1, d = 0.5, Acrit = 1,
          cx = 1, cz = 0.5, u = 1, vc = 1.5, w = 0.5, s = 0.2, q = 1)
sl_t <- c(0, 0.3, 0.9, 1.5, 2.5, 3.2)

sl_closed <- function(p, t = sl_t) {
  ts <- p[["x"]] / (p[["v"]] - p[["kin"]])
  te <- (p[["ve"]] - sqrt(p[["ve"]]^2 - 2 * p[["c"]] * p[["e"]])) / p[["c"]]
  t1 <- p[["ve"]] / p[["c"]]
  Ainf <- p[["k"]] / p[["d"]]
  ta <- -log((Ainf - sqrt(p[["Acrit"]])) / (Ainf - p[["A"]])) / p[["d"]]
  tz <- p[["cz"]] / (p[["q"]] - p[["s"]])
  r1 <- p[["u"]] - p[["vc"]] - p[["w"]]
  r2 <- p[["u"]] - p[["vc"]] - p[["w"]] * p[["s"]] / p[["q"]]
  xz <- p[["cx"]] + r1 * tz
  tx <- tz - xz / r2
  cbind(x = ifelse(t < ts, p[["x"]] - (p[["v"]] - p[["kin"]]) * t, 0),
        y = p[["y"]] + ifelse(t < ts, p[["v"]] * t, p[["x"]] + p[["kin"]] * t),
        e = ifelse(t < te, p[["e"]] + p[["c"]] * t^2 / 2 - p[["ve"]] * t,
                   ifelse(t < t1, 0, p[["c"]] * (t^2 - t1^2) / 2 - p[["ve"]] * (t - t1))),
        A = ifelse(t < ta, Ainf + (p[["A"]] - Ainf) * exp(-p[["d"]] * t), sqrt(p[["Acrit"]])),
        cx = ifelse(t < tz, p[["cx"]] + r1 * t, ifelse(t < tx, xz + r2 * (t - tz), 0)),
        cz = ifelse(t < tz, p[["cz"]] - (p[["q"]] - p[["s"]]) * t, 0))
}
sl_d1 <- function(p, n, h = 1e-6) {
  up <- dn <- p
  up[n] <- up[n] + h
  dn[n] <- dn[n] - h
  (sl_closed(up) - sl_closed(dn)) / (2 * h)
}
sl_d2 <- function(n, m, h = 1e-4) {
  up <- dn <- sl_p
  up[m] <- up[m] + h
  dn[m] <- dn[m] - h
  (sl_d1(up, n) - sl_d1(dn, n)) / (2 * h)
}
sl_pn <- names(sl_p)

sl <- function(nm, ...) cppODE(sl_rhs, modelname = nm, compile = FALSE, ...)
sl_mods <- list(
  d2   = sl("sl_d2", deriv2 = TRUE),
  rev  = sl("sl_rev", derivMode = "reverse"),
  rrb4 = sl("sl_rrb4", method = "rb4", derivMode = "reverse"),
  fr   = sl("sl_fr", derivMode = "forward-reverse"),
  mm   = cppODE(c(x = "kin - v*x/(Km + x)"), modelname = "sl_mm", deriv = FALSE,
                compile = FALSE)
)
do.call(compile, c(unname(sl_mods), list(output = "test_switch_sliding",
                                         cores = test_cores())))
has_cvode <- isTRUE(cvodeConfig$available)
if (has_cvode) {
  sl_cv <- list(
    fwd = cvode(sl_rhs, modelname = "sl_cv", deriv = TRUE, compile = FALSE),
    rev = cvode(sl_rhs, modelname = "sl_cv_rev", derivMode = "reverse", compile = FALSE)
  )
  do.call(compile, c(unname(sl_cv), list(output = "test_switch_sliding_cv",
                                         cores = test_cores())))
}

tight <- function(m, ...)
  solveODE(m, sl_t, sl_p, abstol = 1e-10, reltol = 1e-10, roottol = 1e-10, ...)
set.seed(7)
sl_W <- array(rnorm(length(sl_closed(sl_p))), c(dim(sl_closed(sl_p)), 1))
sl_grad <- sapply(sl_pn, function(n) sum(sl_d1(sl_p, n) * sl_W[, , 1]))

test_that("the slides hold their surfaces, first and second order", {
  out <- tight(sl_mods$d2)
  expect_equal(unname(out$variable), unname(sl_closed(sl_p)), tolerance = 1e-8)
  for (n in sl_pn) {
    expect_equal(unname(out$tangent[, , n]), unname(sl_d1(sl_p, n)), tolerance = 1e-6,
                 info = n)
    for (m in sl_pn)
      expect_equal(unname(out$hessian[, , n, m]), unname(sl_d2(n, m)), tolerance = 1e-4,
                   info = paste(n, m))
  }
  expect_gt(max(abs(out$hessian[sl_t > 2, "e", "c", "c"])), 1)
  expect_gt(max(abs(out$hessian[sl_t > 1.5, "A", "Acrit", "Acrit"])), 0.1)
})

test_that("the reverse gradient and the curvature are those of the closed form", {
  for (nm in c("rev", "rrb4")) {
    out <- tight(sl_mods[[nm]], cotangent = sl_W)
    expect_equal(unname(out$cotangent[sl_pn, 1]), unname(sl_grad), tolerance = 1e-7, info = nm)
  }
  out <- tight(sl_mods$rev, cotangent = sl_W, adjoint = adjointControl(refine = TRUE))
  expect_equal(unname(out$cotangent[sl_pn, 1]), unname(sl_grad), tolerance = 1e-7)
  H <- tight(sl_mods$d2)$hessian
  curv <- apply(H * as.vector(sl_W[, , 1]), c(3, 4), sum)[sl_pn, sl_pn]
  fr <- tight(sl_mods$fr, cotangent = sl_W)
  expect_equal(unname(fr$curvature[sl_pn, sl_pn, 1]), unname(curv), tolerance = 1e-6)
})

test_that("Michaelis-Menten with a small Km runs into the slide", {
  out <- solveODE(sl_mods$mm, sl_t, c(x = 1, kin = 0.2, v = 0.8, Km = 1e-7),
                  abstol = 1e-10, reltol = 1e-10)
  expect_equal(unname(out$variable[, "x"]), unname(sl_closed(sl_p)[, "x"]), tolerance = 1e-5)
})

test_that("cvode() slides, enters and leaves, forward and reverse", {
  skip_if_not(has_cvode, "CVODE not available")
  out <- tight(sl_cv$fwd)
  expect_equal(unname(out$variable), unname(sl_closed(sl_p)), tolerance = 1e-8)
  for (n in sl_pn)
    expect_equal(unname(out$tangent[, , n]), unname(sl_d1(sl_p, n)), tolerance = 1e-6,
                 info = n)
  # The quadrature of the reverse sweep is error-controlled where f reads the
  # time explicitly.
  rev <- tight(sl_cv$rev, cotangent = sl_W, adjoint = adjointControl(gradtol = 1e-10))
  expect_equal(unname(rev$cotangent[sl_pn, 1]), unname(sl_grad), tolerance = 1e-6)
})
