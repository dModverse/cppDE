## =====================================================================
##  problems-classic.R -- the standard stiff IVP test problems.
##
##  The PEtab collection covers systems-biology models: moderately
##  stiff, many parameters, short horizons.  These add the failure modes
##  it does not exercise -- extreme stiffness, scale separation over
##  eleven decades, limit cycles, and a sparse Jacobian that grows with
##  the discretisation.
##
##  Sources: Hairer & Wanner, *Solving Ordinary Differential Equations
##  II*, 2nd ed. (1996), and the accompanying Bari/CWI IVP test set.
## =====================================================================

## `atol` pins the absolute tolerance for problems where the swept value
## would be meaningless -- see E5, whose solution components decay to
## ~1e-30, so any atol from the sweep sits far above the solution itself
## and the step size underflows.
bench_problem <- function(id, name, rhs, parms, times, sens,
                          source = "classic", notes = character(0),
                          traits = character(0), atol = NULL) {
  stopifnot(all(sens %in% names(parms)))
  p <- list(
    id = id, name = name, source = source,
    rhs = rhs, parms = parms, times = times, atol = atol,
    sens = sens, sens_available = sens,
    fixed = setdiff(names(parms), sens),
    condition = NA_character_,
    nstates = length(rhs), npars = length(parms) - length(rhs),
    nsens = length(sens), nevents = 0L, usable = TRUE, notes = notes
  )
  p$traits <- unique(c(traits, infer_traits(p)))
  structure(p, class = "bench_problem")
}


classic_problems <- function(which = NULL) {
  p <- list()

  ## -- Robertson (1966): three-species autocatalysis, the canonical
  ##    stiff test.  Rate constants span nine decades.
  p$robertson <- bench_problem(
    "rob", "Robertson",
    c(y1 = "-k1*y1 + k2*y2*y3",
      y2 = " k1*y1 - k2*y2*y3 - k3*y2^2",
      y3 = " k3*y2^2"),
    c(y1 = 1, y2 = 0, y3 = 0, k1 = 0.04, k2 = 1e4, k3 = 3e7),
    c(0, 10^seq(-5, 7, length.out = 200)),
    sens = c("k1", "k2", "k3"),
    traits = c("stiff-extreme", "mass-action"))

  ## -- HIRES: High Irradiance RESponse of plant photomorphogenesis
  ##    (Schaefer 1975); 8 states, moderately stiff.
  p$hires <- bench_problem(
    "hires", "HIRES",
    c(y1 = "-1.71*y1 + 0.43*y2 + 8.32*y3 + c1",
      y2 = " 1.71*y1 - 8.75*y2",
      y3 = "-10.03*y3 + 0.43*y4 + 0.035*y5",
      y4 = " 8.32*y2 + 1.71*y3 - 1.12*y4",
      y5 = "-1.745*y5 + 0.43*y6 + 0.43*y7",
      y6 = "-c2*y6*y8 + 0.69*y4 + 1.71*y5 - 0.43*y6 + 0.69*y7",
      y7 = " c2*y6*y8 - 1.81*y7",
      y8 = "-c2*y6*y8 + 1.81*y7"),
    c(y1 = 1, y2 = 0, y3 = 0, y4 = 0, y5 = 0, y6 = 0, y7 = 0, y8 = 0.0057,
      c1 = 0.0007, c2 = 280),
    seq(0, 321.8122, length.out = 200),
    sens = c("c1", "c2"),
    traits = "stiff-moderate")

  ## -- Oregonator (Field & Noyes 1974): Belousov-Zhabotinsky
  ##    oscillator, stiff limit cycle.
  p$orego <- bench_problem(
    "orego", "OREGO",
    c(y1 = "s * (y2 + y1*(1 - q*y1 - y2))",
      y2 = "(y3 - (1 + y1)*y2) / s",
      y3 = "w * (y1 - y3)"),
    c(y1 = 1, y2 = 2, y3 = 3, s = 77.27, q = 8.375e-6, w = 0.161),
    seq(0, 360, length.out = 400),
    sens = c("s", "q", "w"),
    traits = c("oscillatory", "stiff-moderate"))

  ## -- E5: chemical pyrolysis; solution components differ by twelve
  ##    orders of magnitude and the horizon spans 1e-5 .. 1e11.
  p$e5 <- bench_problem(
    "e5", "E5",
    c(y1 = "-A*y1 - B*y1*y3",
      y2 = " A*y1 - M*C*y2*y3",
      y3 = " A*y1 - B*y1*y3 - M*C*y2*y3 + C*y4",
      y4 = " B*y1*y3 - C*y4"),
    c(y1 = 1.76e-3, y2 = 0, y3 = 0, y4 = 0,
      A = 7.89e-10, B = 1.1e7, C = 1.13e3, M = 1e6),
    c(0, 10^seq(-5, 11, length.out = 200)),
    sens = c("A", "B", "C"),
    traits = c("stiff-extreme", "mass-action"),
    ## The IVP test set runs E5 at atol = 1.7e-24: y2..y4 decay to about
    ## 1e-30, so a swept atol of 1e-6..1e-12 exceeds the solution by
    ## eighteen orders of magnitude and both solvers stall on step-size
    ## underflow long before t = 1e11.
    atol = 1.7e-24,
    notes = "absolute tolerance pinned to 1.7e-24 (IVP test set convention)")

  ## -- Van der Pol at mu = 1000: relaxation oscillator whose stiffness
  ##    switches on and off along the trajectory.
  p$vdp <- bench_problem(
    "vdp", "VanDerPol_mu1000",
    c(x = "y", y = "mu * (1 - x^2) * y - x"),
    c(x = 2, y = 0, mu = 1000),
    seq(0, 3000, length.out = 500),
    sens = "mu",
    traits = c("relaxation", "oscillatory", "stiff-moderate"))

  ## -- Pollution: 20 species, 25 reactions of atmospheric chemistry
  ##    (Verwer 1994).  Sparse but small.
  p$pollution <- local({
    k <- c(.35, .266e2, .123e5, .86e-3, .82e-3, .15e5, .13e-3, .24e5, .165e5, .9e4,
           .22e-1, .12e5, .188e1, .163e5, .48e7, .35e-3, .175e-1, .1e9, .444e12, .124e4,
           .21e1, .578e1, .474e-1, .178e4, .312e1)
    y0 <- c(y1=0, y2=.2, y3=0, y4=.04, y5=0, y6=0, y7=.1, y8=.3, y9=.01, y10=0,
            y11=0, y12=0, y13=0, y14=0, y15=0, y16=0, y17=.007, y18=0, y19=0, y20=0)
    rhs <- c(
      y1  = "-k1*y1 - k10*y11*y1 - k14*y1*y6 - k23*y1*y4 - k24*y19*y1 + k2*y2*y4 + k3*y5*y2 + k9*y11*y2 + k11*y13 + k12*y10*y2 + k22*y19 + k25*y20",
      y2  = "-k2*y2*y4 - k3*y5*y2 - k9*y11*y2 - k12*y10*y2 + k1*y1 + k21*y19",
      y3  = "-k15*y3 + k1*y1 + k17*y4 + k19*y16 + k22*y19",
      y4  = "-k2*y2*y4 - k16*y4 - k17*y4 - k23*y1*y4 + k15*y3",
      y5  = "-k3*y5*y2 + 2*k4*y7 + k6*y7*y6 + k7*y9 + k13*y14 + k20*y17*y6",
      y6  = "-k6*y7*y6 - k8*y9*y6 - k14*y1*y6 - k20*y17*y6 + k3*y5*y2 + 2*k18*y16",
      y7  = "-k4*y7 - k5*y7 - k6*y7*y6 + k13*y14",
      y8  = "k4*y7 + k5*y7 + k6*y7*y6 + k7*y9",
      y9  = "-k7*y9 - k8*y9*y6",
      y10 = "-k12*y10*y2 + k7*y9 + k9*y11*y2",
      y11 = "-k9*y11*y2 - k10*y11*y1 + k8*y9*y6 + k11*y13",
      y12 = "k9*y11*y2",
      y13 = "-k11*y13 + k10*y11*y1",
      y14 = "-k13*y14 + k12*y10*y2",
      y15 = "k14*y1*y6",
      y16 = "-k18*y16 - k19*y16 + k16*y4",
      y17 = "-k20*y17*y6",
      y18 = "k20*y17*y6",
      y19 = "-k21*y19 - k22*y19 - k24*y19*y1 + k23*y1*y4 + k25*y20",
      y20 = "-k25*y20 + k24*y19*y1")
    bench_problem("poll", "Pollution", rhs,
                  c(y0, stats::setNames(k, paste0("k", seq_along(k)))),
                  seq(0, 60, length.out = 200),
                  sens = paste0("k", 1:25),
                  traits = c("stiff-moderate", "mass-action"))
  })

  ## -- Brusselator, 1D method of lines.  The one problem here whose
  ##    size is a free knob: N grid points give 2N states with a banded
  ##    Jacobian, which is what the sparse path is for.
  ##    The small variants exist so that the `tiny` tier still covers the
  ##    sparse path without paying for 128 states.
  p$brusselator_small <- build_brusselator(24L)
  p$brusselator       <- build_brusselator(64L)

  ## -- FitzHugh-Nagumo chain: N diffusively coupled neurons, stiffness
  ##    set by eps.  Complements Brusselator with a cubic nonlinearity.
  p$fhn_small <- build_fhn_chain(24L)
  p$fhn       <- build_fhn_chain(64L)

  if (!is.null(which)) p <- p[intersect(which, names(p))]
  p
}


## N cells -> 2N states.  Periodic boundary, so the Jacobian is banded
## with two corner entries.
build_brusselator <- function(N = 64L) {
  N <- as.integer(N)
  u <- sprintf("u%03d", seq_len(N)); v <- sprintf("v%03d", seq_len(N))
  wrap <- function(i) ((i - 1L) %% N) + 1L
  rhs <- stats::setNames(character(2L * N), c(u, v))
  for (i in seq_len(N)) {
    lapu <- sprintf("(%s + %s - 2*%s)", u[wrap(i - 1L)], u[wrap(i + 1L)], u[i])
    lapv <- sprintf("(%s + %s - 2*%s)", v[wrap(i - 1L)], v[wrap(i + 1L)], v[i])
    rhs[u[i]] <- sprintf("A + %s^2*%s - (B + 1)*%s + alpha*%s", u[i], v[i], u[i], lapu)
    rhs[v[i]] <- sprintf("B*%s - %s^2*%s + alpha*%s", u[i], u[i], v[i], lapv)
  }
  x <- (seq_len(N) - 0.5) / N
  parms <- c(stats::setNames(1 + sin(2 * pi * x), u),
             stats::setNames(rep(3, N), v),
             A = 1, B = 3, alpha = 0.02 * N^2)
  bench_problem(sprintf("bruss%d", N), sprintf("Brusselator1D_N%d", N),
                rhs, parms, seq(0, 10, length.out = 100),
                sens = c("A", "B", "alpha"),
                traits = c("sparse", "stiff-moderate"),
                notes = sprintf("method of lines, %d cells", N))
}

build_fhn_chain <- function(N = 64L, eps = 0.02) {
  N <- as.integer(N)
  v <- sprintf("v%03d", seq_len(N)); w <- sprintf("w%03d", seq_len(N))
  rhs <- stats::setNames(character(2L * N), c(v, w))
  for (i in seq_len(N)) {
    nb <- c(if (i > 1L) i - 1L, if (i < N) i + 1L)
    coup <- if (length(nb))
      sprintf(" + D*(%s)", paste(sprintf("(%s - %s)", v[nb], v[i]), collapse = " + ")) else ""
    rhs[v[i]] <- sprintf("(%s - %s^3/3 - %s + Iext)/eps%s", v[i], v[i], w[i], coup)
    rhs[w[i]] <- sprintf("%s + a - b*%s", v[i], w[i])
  }
  parms <- c(stats::setNames(rep(-1, N), v), stats::setNames(rep(1, N), w),
             a = 0.7, b = 0.8, D = 1, Iext = 0.5, eps = eps)
  bench_problem(sprintf("fhn%d", N), sprintf("FitzHughNagumo_N%d", N),
                rhs, parms, seq(0, 100, length.out = 100),
                sens = c("a", "b", "D", "Iext", "eps"),
                traits = c("sparse", "relaxation", "stiff-moderate"),
                notes = sprintf("%d coupled neurons, eps = %g", N, eps))
}
