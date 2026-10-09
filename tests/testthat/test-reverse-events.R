# The reverse mode on events, forcings, a run stopped at equilibrium, a wide
# model and sparse Jacobians, each checked against one forward oracle on bdf.
#
# The oracle is not sharp: the forward solve adapts under sensitivities and the
# reverse one in plain double, so the two discretisations differ by O(tol).

skip_on_cran()

eqns <- c(A = "-k1 * A + k2 * B",
          B = "k1 * A - k2 * B - k3 * B * B")
pars  <- c(A = 1.2, B = 0.4, k1 = 0.7, k2 = 0.35, k3 = 1.1)
times <- c(0, 0.3, 1, 2.5, 5)
tol   <- list(abstol = 1e-10, reltol = 1e-10)
methods <- c("bdf", "adams", "rb4", "tsit5")

# A dose and a root-triggered jump on a forced variant of eqns.
ev_eqns <- c(A = "-k1 * A + k2 * B + u",
             B = "k1 * A - k2 * B - k3 * B * B")
ev_def  <- data.frame(var    = c("A", "B"),
                      time   = c("t_dose", NA),
                      value  = c("d_amt", 1.5),
                      root   = c(NA, "A - 0.6"),
                      method = c("add", "multiply"),
                      stringsAsFactors = FALSE)
# A forcing that multiplies the states it acts on.
fc_eqns <- c(A = "u * A - k1 * A", B = "k1 * A - k2 * B * u")
# A cascade that settles, for rootfunc = "equilibrate".
eq_eqns <- c(R  = "k_act - k_deact * R",
             A  = "-k1 * A * R + k2 * pA",
             pA = "k1 * A * R - k2 * pA")

# A model wider than the Nordsieck history is deep. Every other model in this
# file has fewer states than a step has slots, so a buffer sized by one and used
# for the other fits, and a contraction writing per state stays inside it.
wide  <- local({
  n <- 10L
  v <- paste0("A", seq_len(n))
  k <- paste0("k", seq_len(n))
  eq <- setNames(character(n), v)
  eq[1] <- paste0("-", k[1], "*", v[1])
  for (i in 2:n)
    eq[i] <- paste0(k[i - 1], "*", v[i - 1], " - ", k[i], "*", v[i], "*", v[i])
  eq
})
wpars <- c(setNames(seq(1.5, 0.6, length.out = 10), paste0("A", 1:10)),
           setNames(seq(0.9, 0.2, length.out = 10), paste0("k", 1:10)))

# Models nest in lists; compile() takes them flat.
flat <- function(x) if (is.list(x)) do.call(c, lapply(unname(x), flat)) else list(x)

# One uncompiled cppODE model per method, named by method.
per_method <- function(rhs, prefix, ms = methods, ...)
  lapply(setNames(nm = ms), function(m)
    cppODE(rhs, method = m, modelname = paste0(prefix, m), compile = FALSE, ...))

# A forward oracle on bdf per configuration and the reverse models it checks.
# wide still turns sparse by auto-detection when KLU is present.
oracle <- function(rhs, nm, ...)
  cppODE(rhs, modelname = nm, deriv = TRUE, compile = FALSE, ...)
models <- list(
  fwd      = oracle(eqns, "reve_f"),
  ev_fwd   = oracle(ev_eqns, "reve_ev_f", events = ev_def, forcings = "u"),
  ev_rev   = per_method(ev_eqns, "reve_ev_r_", events = ev_def,
                        forcings = "u", derivMode = "reverse"),
  fc_fwd   = oracle(fc_eqns, "reve_fc_f", forcings = "u"),
  fc_rev   = per_method(fc_eqns, "reve_fc_r_", c("rb4", "bdf"),
                        forcings = "u", derivMode = "reverse"),
  eq_val   = per_method(eq_eqns, "reve_eq_v_", c("bdf", "rb4"),
                        rootfunc = "equilibrate", deriv = FALSE),
  eq_rev   = per_method(eq_eqns, "reve_eq_r_", c("bdf", "rb4"),
                        rootfunc = "equilibrate", derivMode = "reverse"),
  wide_fwd = oracle(wide, "reve_wide_f"),
  wide_rev = per_method(wide, "reve_wide_r_", c("bdf", "adams"),
                        derivMode = "reverse"))
do.call(compile, c(flat(models),
                   list(output = "test_reverse_events", cores = test_cores())))

# The sparse models need KLU. They all have the Jacobian of eqns, so the KLU
# defines each one records coincide and one shared object serves them all.
if (isTRUE(cvodeConfig$klu_available)) {
  sp_ev <- data.frame(var = "A", time = "t_dose", value = "d_amt",
                      method = "add", stringsAsFactors = FALSE)
  sparse_models <- list(
    ev_fwd = oracle(eqns, "reve_sp_ev_f", events = sp_ev, sparse = TRUE),
    ev_rev = per_method(eqns, "reve_sp_ev_r_", c("bdf", "rb4"), events = sp_ev,
                        sparse = TRUE, derivMode = "reverse"),
    rev    = cppODE(eqns, method = "rb4", modelname = "reve_sp_r", sparse = TRUE,
                    derivMode = "reverse", compile = FALSE))
  do.call(compile, c(flat(sparse_models),
                     list(output = "test_reverse_events_klu", cores = test_cores())))
}

# w' S contracted over times and states, the quantity the sweep returns.
contract <- function(tangent, W) {
  vapply(seq_len(dim(W)[3]),
         function(k) apply(tangent * as.vector(W[, , k]), 3, sum),
         numeric(dim(tangent)[3]))
}

seed_for <- function(res, n_seed = 1L, seed = 1L) {
  set.seed(seed)
  n <- nrow(res$variable); p <- ncol(res$variable)
  array(rnorm(n * p * n_seed), c(n, p, n_seed))
}

test_that("the reverse mode handles events, roots and forcings", {
  p  <- c(A = 1.2, B = 0.4, k1 = 0.7, k2 = 0.35, k3 = 1.1,
          d_amt = 0.4, t_dose = 1.0)
  fc <- list(u = data.frame(time = c(0, 2, 5), value = c(0.1, 0.25, 0.05)))

  fwd <- do.call(solveODE, c(list(models$ev_fwd, times, p, forcings = fc), tol))
  # Both jumps have to be in the run, or the test proves nothing.
  expect_gt(nrow(fwd$variable), length(times))
  W   <- seed_for(fwd)
  ref <- contract(fwd$tangent, W)[, 1]

  for (m in methods) {
    rev <- do.call(solveODE, c(list(models$ev_rev[[m]], times, p, forcings = fc,
                                    cotangent = W), tol))
    expect_equal(rev$variable, fwd$variable, tolerance = 1e-6, info = m)
    expect_equal(unname(rev$cotangent[names(ref), 1]), unname(ref),
                 tolerance = 1e-5, info = m)
  }
})

test_that("a sparse model jumps backwards too", {
  skip_if_not(isTRUE(cvodeConfig$klu_available), "KLU not available")
  pe <- c(pars, t_dose = 1.0, d_amt = 0.4)
  fwd <- do.call(solveODE, c(list(sparse_models$ev_fwd, times, pe), tol))
  W   <- seed_for(fwd)
  ref <- contract(fwd$tangent, W)[, 1]
  for (m in c("bdf", "rb4")) {
    rv  <- do.call(solveODE, c(list(sparse_models$ev_rev[[m]], times, pe,
                                    cotangent = W), tol))
    expect_equal(unname(rv$cotangent[names(ref), 1]), unname(ref),
                 tolerance = 1e-5, info = m)
  }
})

test_that("a written Rosenbrock adjoint handles a multiplicative forcing", {
  # Four of the six stages add a multiple of df/dt, so its derivative in the
  # state and in the parameters is part of the adjoint. A forcing reaches df/dt
  # through a chain term the Jacobian emitter appends, and multiplicatively is
  # the way that term keeps a state in it.
  pf <- c(A = 1.1, B = 0.3, k1 = 0.8, k2 = 0.45)
  fc <- list(u = data.frame(time = c(0, 1, 3, 5), value = c(0.2, 0.5, 0.1, 0.3)))
  fwd <- do.call(solveODE, c(list(models$fc_fwd, times, pf, forcings = fc), tol))
  W   <- seed_for(fwd)
  ref <- contract(fwd$tangent, W)[, 1]
  for (m in c("rb4", "bdf")) {
    rv  <- do.call(solveODE, c(list(models$fc_rev[[m]], times, pf, forcings = fc,
                                    cotangent = W), tol))
    expect_equal(unname(rv$cotangent[names(ref), 1]), unname(ref),
                 tolerance = 1e-5, info = m)
  }
})

test_that("an equilibrated run goes backwards", {
  # rootfunc stops the run when the state stops moving, so the output grid is
  # the run's own. A sensitivity run does not produce the same one: it adapts
  # on the tangents as well and reaches the root elsewhere, so there is no
  # forward gradient on this grid to compare against. What is asserted is that
  # the backward pass answers on the grid its own forward pass produced, and
  # that this pass is the value run to the last bit.
  pe <- c(R = 1, A = 1, pA = 0, k_act = 0.1, k_deact = 0.7, k1 = 0.1, k2 = 0.05)
  tt <- seq(0, 1e3, length.out = 50)
  for (m in c("bdf", "rb4")) {
    mv <- models$eq_val[[m]]
    mr <- models$eq_rev[[m]]
    val <- do.call(solveODE, c(list(mv, tt, pe), tol))
    expect_lt(nrow(val$variable), length(tt))   # it really did stop early
    W   <- seed_for(val)
    rv  <- do.call(solveODE, c(list(mr, tt, pe, cotangent = W), tol))
    expect_identical(rv$variable, val$variable)

    # With it, the rows come from the reverse run's own values.
    fp  <- do.call(solveODE, c(list(mr, tt, pe, keepStore = TRUE), tol))
    expect_lt(nrow(fp$variable), length(tt))
    rv  <- do.call(solveODE, c(list(mr, tt, pe, cotangent = seed_for(fp),
                                    store = fp$store), tol))
    expect_equal(dim(rv$cotangent), c(length(pe), 1L), info = m)
    expect_true(all(is.finite(rv$cotangent)), info = m)
    expect_gt(max(abs(rv$cotangent)), 1e-6)
  }
})

test_that("rb4 goes backwards on a sparse Jacobian", {
  skip_if_not(isTRUE(cvodeConfig$klu_available), "KLU not available")
  fwd <- do.call(solveODE, c(list(models$fwd, times, pars), tol))
  W   <- seed_for(fwd)
  rv  <- do.call(solveODE, c(list(sparse_models$rev, times, pars, cotangent = W), tol))
  ref <- contract(fwd$tangent, W)[, 1]
  expect_equal(unname(rv$cotangent[names(ref), 1]), unname(ref), tolerance = 1e-5)
})

test_that("the reverse mode handles a model wider than its history is deep", {
  fwd <- do.call(solveODE, c(list(models$wide_fwd, times, wpars), tol))
  W   <- seed_for(fwd)
  ref <- contract(fwd$tangent, W)[, 1]
  for (m in c("bdf", "adams")) {
    rv  <- do.call(solveODE, c(list(models$wide_rev[[m]], times, wpars,
                                    cotangent = W), tol))
    expect_equal(unname(rv$cotangent[names(ref), 1]), unname(ref),
                 tolerance = 1e-5, info = m)
  }
})
