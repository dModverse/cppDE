# solveODEBatch(): results must match the serial path exactly, not merely
# closely -- each condition runs the same steps on thread-local state.
# cores is pinned to 2 throughout: testthat runs files in parallel processes,
# so letting each grab detectCores() would oversubscribe badly.

skip_on_cran()

decay <- c(A = "-k * A", B = "k * A")
tt    <- seq(0, 3, by = 0.25)
ks    <- c(0.3, 0.7, 1.4, 2.2)
conds <- lapply(ks, function(k) list(parms = c(A = 1, B = 0, k = k)))
names(conds) <- paste0("k", seq_along(ks))

serial_ref <- function(model, times = tt) {
  lapply(ks, function(k) solveODE(model, times = times, parms = c(A = 1, B = 0, k = k)))
}

expect_batch_identical <- function(bat, ser, what = c("time", "variable")) {
  expect_equal(length(bat), length(ser))
  for (i in seq_along(ser))
    for (el in what)
      expect_identical(bat[[i]][[el]], ser[[i]][[el]],
                       info = paste0("condition ", i, ", element ", el))
}

# -- Bit-identical against the serial path ------------------------------------

test_that("solveODEBatch matches solveODE exactly, with sensitivities", {
  m <- cppODE(decay, modelname = "batch_sens", deriv = TRUE)
  ser <- serial_ref(m)
  bat <- solveODEBatch(m, conds, times = tt, cores = 2)

  expect_named(bat, names(conds))
  expect_batch_identical(bat, ser, c("time", "variable", "sens1"))
})

test_that("solveODEBatch matches solveODE exactly without sensitivities", {
  m <- cppODE(decay, modelname = "batch_nosens", deriv = FALSE)
  ser <- serial_ref(m)
  bat <- solveODEBatch(m, conds, times = tt, cores = 2)
  expect_batch_identical(bat, ser)
})

test_that("second-order sensitivities survive the batch path", {
  m <- cppODE(decay, modelname = "batch_d2", deriv = TRUE, deriv2 = TRUE, nStack = 3L)
  ser <- serial_ref(m)
  bat <- solveODEBatch(m, conds, times = tt, cores = 2)
  expect_batch_identical(bat, ser, c("time", "variable", "sens1", "sens2"))
})

# The arena is thread-local and pops when solve_impl returns, so heap AD is
# the case where a result that was not flattened in time would show up.
test_that("heap AD (nStack = Inf) batches correctly", {
  m <- cppODE(decay, modelname = "batch_heap", deriv = TRUE, nStack = Inf)
  ser <- serial_ref(m)
  bat <- solveODEBatch(m, conds, times = tt, cores = 2)
  expect_batch_identical(bat, ser, c("time", "variable", "sens1"))
})

# -- Thread count must not change the answer ----------------------------------

test_that("results are invariant in the number of threads", {
  m <- cppODE(decay, modelname = "batch_threads_a", deriv = TRUE)
  one  <- solveODEBatch(m, conds, times = tt, cores = 1)
  many <- solveODEBatch(m, conds, times = tt, cores = 4)
  for (i in seq_along(conds)) {
    expect_identical(one[[i]]$variable, many[[i]]$variable)
    expect_identical(one[[i]]$sens1,    many[[i]]$sens1)
  }
})

# -- One failing condition must not take the others down ----------------------

test_that("a failing condition is isolated and reported", {
  m <- cppODE(decay, modelname = "batch_fail", deriv = FALSE)
  cs <- list(ok1 = list(parms = c(A = 1, B = 0, k = 0.5)),
             bad = list(parms = c(A = 1, B = 0, k = 0.5)),
             ok2 = list(parms = c(A = 1, B = 0, k = 1.5)))

  # maxsteps applies batch-wide, so starve every condition and check that the
  # per-condition return codes come back rather than an error being thrown.
  expect_warning(solveODEBatch(m, cs, times = tt, maxsteps = 2L, cores = 2),
                 "did not complete")
  starved <- suppressWarnings(
    solveODEBatch(m, cs, times = tt, maxsteps = 2L, cores = 2))
  expect_length(starved, 3L)
  expect_true(all(vapply(starved, function(r) r$diagnostics$return_code, integer(1)) != 0L))

  # silent suppresses the warning but keeps the codes
  quiet <- solveODEBatch(m, cs, times = tt, maxsteps = 2L, cores = 2,
                         onFailure = "silent")
  expect_length(quiet, 3L)

  # and a healthy batch warns about nothing
  fine <- solveODEBatch(m, cs, times = tt, cores = 2)
  expect_true(all(vapply(fine, function(r) r$diagnostics$return_code, integer(1)) == 0L))
})

# -- Per-condition overrides and argument handling ----------------------------

test_that("per-condition arguments override the batch-wide ones", {
  m <- cppODE(decay, modelname = "batch_over", deriv = FALSE)
  cs <- list(short = list(parms = c(A = 1, B = 0, k = 0.5), times = seq(0, 1, 0.5)),
             long  = list(parms = c(A = 1, B = 0, k = 0.5)))
  bat <- solveODEBatch(m, cs, times = tt, cores = 2)

  expect_equal(bat$short$time, seq(0, 1, 0.5), tolerance = 1e-12)
  expect_equal(bat$long$time,  tt,             tolerance = 1e-12)
})

test_that("solveODEBatch rejects malformed input", {
  m <- cppODE(decay, modelname = "batch_input", deriv = FALSE)
  expect_error(solveODEBatch(m, list(), times = tt), "non-empty")
  expect_error(solveODEBatch(m, list(list(nope = 1)), times = tt), "unknown per-condition")
  expect_error(solveODEBatch(m, list(list(parms = c(A = 1, B = 0, k = 1)))),
               "no 'times'")
  expect_error(
    solveODEBatch(m, conds, times = tt, cores = 0),
    "positive integer")
})

# -- Optional backends --------------------------------------------------------

test_that("the CVODE backend batches like the native one", {
  skip_if_not(isTRUE(cvodeConfig$available), "CVODE backend not available")
  m <- cvode(decay, modelname = "batch_cv", deriv = TRUE)
  ser <- serial_ref(m)
  bat <- solveODEBatch(m, conds, times = tt, cores = 2)
  expect_batch_identical(bat, ser, c("time", "variable", "sens1"))
})

test_that("the sparse KLU path batches correctly", {
  skip_if_not(isTRUE(cvodeConfig$klu_available), "KLU not available")
  n <- 10L
  nm <- paste0("S", seq_len(n))
  eqs <- setNames(c("-r * S1", paste0("r * S", seq_len(n - 1L), " - r * S", seq_len(n - 1L) + 1L)), nm)
  eqs[n] <- paste0("r * S", n - 1L)
  m <- cppODE(eqs, modelname = "batch_sparse", deriv = TRUE, sparse = TRUE)
  expect_true(isTRUE(attr(m, "sparse")))

  p0 <- setNames(c(1, rep(0, n - 1L)), nm)
  cs <- lapply(c(0.4, 0.9), function(r) list(parms = c(p0, r = r)))
  ser <- lapply(c(0.4, 0.9), function(r)
    solveODE(m, times = tt, parms = c(p0, r = r)))
  bat <- solveODEBatch(m, cs, times = tt, cores = 2)
  for (i in seq_along(ser)) {
    expect_identical(bat[[i]]$variable, ser[[i]]$variable)
    expect_identical(bat[[i]]$sens1,    ser[[i]]$sens1)
  }
})


## ---- prepared batch handles ----------------------------------------------

test_that("solveBatch on a prepared handle equals solveODEBatch", {
  m <- cppODE(decay, modelname = "batch_prep_a", deriv = TRUE)
  h <- prepareBatch(m, conds, times = tt)
  a <- solveBatch(h, cores = 2L)
  b <- solveODEBatch(m, conds, times = tt, cores = 2L)
  expect_named(a, names(conds))
  for (i in seq_along(conds)) {
    expect_identical(a[[i]]$variable, b[[i]]$variable)
    expect_identical(a[[i]]$sens1,    b[[i]]$sens1)
  }
})


test_that("solveBatch re-solves with new parameters", {
  m  <- cppODE(decay, modelname = "batch_prep2", deriv = TRUE)
  h  <- prepareBatch(m, conds, times = tt)
  np <- lapply(conds, function(c) { p <- c$parms; p["k"] <- p["k"] * 1.5; p })
  a  <- solveBatch(h, parms = np, cores = 2L)
  b  <- solveODEBatch(m, lapply(np, function(p) list(parms = p)),
                      times = tt, cores = 2L)
  for (i in seq_along(conds))
    expect_identical(a[[i]]$variable, b[[i]]$variable)

  # the handle is not consumed: solving again with the originals still works
  expect_identical(solveBatch(h, cores = 1L)[[1]]$variable,
                   solveODEBatch(m, conds, times = tt, cores = 1L)[[1]]$variable)
})


test_that("solveBatch rejects input it cannot reuse", {
  m <- cppODE(decay, modelname = "batch_prep", deriv = TRUE)
  h <- prepareBatch(m, conds, times = tt)
  expect_error(solveBatch(h, parms = conds[[1]]$parms), "list with one element")
  expect_error(solveBatch(h, parms = lapply(conds, function(c) unname(c$parms))),
               "must be named")
  expect_error(solveBatch("nope"), "prepareBatch")
})


test_that("batchAvailable reports why a batch would be serial", {
  m <- cppODE(decay, modelname = "batch_avail", deriv = TRUE)
  a <- batchAvailable(m)
  expect_type(a, "list")
  expect_named(a, c("symbol", "openmp", "modelOpenmp", "parallel"))
  expect_true(a$symbol)
  expect_identical(a$parallel, a$symbol && a$openmp && a$modelOpenmp)
})


test_that("solveODEBatch reports the thread count it used", {
  m <- cppODE(decay, modelname = "batch_threads", deriv = TRUE)
  skip_if_not(isTRUE(batchAvailable(m)$parallel),
              "batch falls back to a serial loop")
  out <- solveODEBatch(m, conds, times = tt, cores = 2L)
  expect_identical(attr(out, "threads"), 2L)
})


## ---- per-condition inputs -------------------------------------------------

test_that("conditions may carry their own sens1ini labels", {
  # dMod's normal case: each condition depends on a different outer parameter
  # set, so the batch cannot hand one shared dimnames pair to the C++ side.
  m <- cppODE(decay, modelname = "batch_hetsens", deriv = TRUE)
  mk <- function(k, lab) {
    s <- matrix(0, 3, 2, dimnames = list(c("A", "B", "k"), c("shared", lab)))
    s["A", 1] <- 1; s["k", 2] <- 1
    list(parms = c(A = 1, B = 0, k = k), sens1ini = s)
  }
  cs <- list(a = mk(0.3, "pa"), b = mk(0.7, "pb"), c = mk(1.4, "pc"))
  bat <- solveODEBatch(m, cs, times = tt, cores = 2L)

  for (i in seq_along(cs)) {
    ser <- solveODE(m, times = tt, parms = cs[[i]]$parms,
                    sens1ini = cs[[i]]$sens1ini)
    expect_identical(bat[[i]]$variable, ser$variable)
    expect_identical(bat[[i]]$sens1, ser$sens1)
    expect_identical(dimnames(bat[[i]]$sens1)[[3]], colnames(cs[[i]]$sens1ini))
  }
})


test_that("conditions may fix different parameters", {
  m <- cppODE(decay, modelname = "batch_hetfixed", deriv = TRUE)
  cs <- list(free  = list(parms = c(A = 1, B = 0, k = 0.5)),
             fixk  = list(parms = c(A = 1, B = 0, k = 0.5), fixed = "k"),
             fixA  = list(parms = c(A = 1, B = 0, k = 0.5), fixed = "A"))
  bat <- solveODEBatch(m, cs, times = tt, cores = 2L)

  widths <- vapply(bat, function(b) dim(b$sens1)[3], 0L)
  expect_true(widths[["free"]] > widths[["fixk"]])
  for (i in seq_along(cs)) {
    ser <- solveODE(m, times = tt, parms = cs[[i]]$parms, fixed = cs[[i]]$fixed)
    expect_identical(bat[[i]]$sens1, ser$sens1)
  }
})


test_that("conditions may carry their own forcings", {
  m <- cppODE(c(A = "u - k * A"), forcings = "u",
              modelname = "batch_forcings", deriv = TRUE)
  mkf <- function(a) list(u = data.frame(time = c(0, 1, 2, 3),
                                         value = a * c(0, 1, 1, 0)))
  cs <- list(lo = list(parms = c(A = 0, k = 0.5), forcings = mkf(0.5)),
             hi = list(parms = c(A = 0, k = 0.5), forcings = mkf(2.0)))
  bat <- solveODEBatch(m, cs, times = tt, cores = 2L)

  for (i in seq_along(cs)) {
    ser <- solveODE(m, times = tt, parms = cs[[i]]$parms,
                    forcings = cs[[i]]$forcings)
    expect_identical(bat[[i]]$variable, ser$variable)
  }
  expect_gt(max(bat$hi$variable[, "A"]), max(bat$lo$variable[, "A"]))
})


test_that("conditions may carry their own solver options", {
  m <- cppODE(decay, modelname = "batch_hetopts", deriv = FALSE)
  cs <- list(loose = list(parms = c(A = 1, B = 0, k = 0.5), abstol = 1e-3,
                          reltol = 1e-3),
             tight = list(parms = c(A = 1, B = 0, k = 0.5), abstol = 1e-10,
                          reltol = 1e-10))
  bat <- solveODEBatch(m, cs, times = tt, cores = 2L)
  for (i in seq_along(cs)) {
    ser <- do.call(solveODE, c(list(m, times = tt), cs[[i]]))
    expect_identical(bat[[i]]$variable, ser$variable)
  }
  expect_gt(bat$tight$diagnostics$accepted, bat$loose$diagnostics$accepted)
})


test_that("one failing condition does not take the others down", {
  m <- cppODE(decay, modelname = "batch_partialfail", deriv = FALSE)
  cs <- list(ok1  = list(parms = c(A = 1, B = 0, k = 0.5)),
             bad  = list(parms = c(A = 1, B = 0, k = 0.5), maxsteps = 2L),
             ok2  = list(parms = c(A = 1, B = 0, k = 1.2)))
  expect_warning(solveODEBatch(m, cs, times = tt, cores = 2L), "did not complete")
  bat <- solveODEBatch(m, cs, times = tt, cores = 2L, onFailure = "silent")
  rc <- vapply(bat, function(b) b$diagnostics$return_code, 0L)
  expect_identical(unname(rc[c("ok1", "ok2")]), c(0L, 0L))
  expect_true(rc[["bad"]] != 0L)
  expect_identical(bat$ok1$variable,
                   solveODE(m, times = tt, parms = cs$ok1$parms)$variable)
})


test_that("a batch with events matches the serial path", {
  eqns <- c(A = "-k1 * A")
  evt  <- data.frame(var = "A", time = "t_e", value = "dose", method = "add",
                     root = NA, stringsAsFactors = FALSE)
  m <- cppODE(eqns, events = evt, modelname = "batch_events", deriv = TRUE)
  te <- seq(0, 50, length.out = 60)
  cs <- lapply(c(0.2, 0.6), function(d)
    list(parms = c(A = 1, k1 = 0.1, t_e = 25, dose = d)))
  bat <- solveODEBatch(m, cs, times = te, cores = 2L)
  for (i in seq_along(cs)) {
    ser <- solveODE(m, times = te, parms = cs[[i]]$parms)
    expect_identical(bat[[i]]$variable, ser$variable)
    expect_identical(bat[[i]]$sens1, ser$sens1)
  }
})


test_that("a time event takes the preallocated path when it is on the grid", {
  # A time event does not change n_out as long as its time is a requested
  # output time, so the batch can size the results up front. When it is not,
  # the extra row makes the prediction wrong and the sink must decline rather
  # than write out of bounds -- both paths have to match the serial solve.
  eqns <- c(A = "-k1 * A")
  evt  <- data.frame(var = "A", time = "t_e", value = "dose", method = "add",
                     root = NA, stringsAsFactors = FALSE)
  m <- cppODE(eqns, events = evt, modelname = "batch_evt_grid", deriv = TRUE)
  off <- seq(0, 50, length.out = 60)          # 25 is not a grid point
  on  <- sort(unique(c(off, 25)))
  cs  <- lapply(c(0.2, 5), function(d)
    list(parms = c(A = 1, k1 = 0.1, t_e = 25, dose = d)))

  expect_length(solveODE(m, times = on,  parms = cs[[1]]$parms)$time, length(on))
  expect_length(solveODE(m, times = off, parms = cs[[1]]$parms)$time, length(off) + 1L)

  for (tt in list(on, off)) {
    bat <- solveODEBatch(m, cs, times = tt, cores = 2L)
    for (i in seq_along(cs)) {
      ser <- solveODE(m, times = tt, parms = cs[[i]]$parms)
      expect_identical(bat[[i]]$time,     ser$time)
      expect_identical(bat[[i]]$variable, ser$variable)
      expect_identical(bat[[i]]$sens1,    ser$sens1)
    }
  }

  # The event time is a parameter, so it differs per condition; the batch has to
  # evaluate it per condition to predict each grid.
  cs2 <- lapply(c(10.3, 25.7, 40.1), function(te)
    list(parms = c(A = 1, k1 = 0.1, t_e = te, dose = 3)))
  b2 <- solveODEBatch(m, cs2, times = off, cores = 3L)
  for (i in seq_along(cs2)) {
    ser <- solveODE(m, times = off, parms = cs2[[i]]$parms)
    expect_identical(b2[[i]]$time,     ser$time)
    expect_identical(b2[[i]]$variable, ser$variable)
    expect_identical(b2[[i]]$sens1,    ser$sens1)
    expect_length(b2[[i]]$time, length(off) + 1L)
  }
})


test_that("a root event keeps the dynamic path and still matches the serial solve", {
  # A root event's firing time is not known before the solve, so the batch must
  # decline to size the output up front rather than guess.
  evt <- data.frame(var = "A", time = NA, root = "A - 0.5", value = "0.9",
                    method = "multiply", stringsAsFactors = FALSE)
  m <- cppODE(c(A = "-k1 * A"), events = evt, modelname = "batch_evt_root",
              deriv = TRUE)
  tt <- seq(0, 50, length.out = 60)
  cs <- lapply(c(0.08, 0.12), function(k) list(parms = c(A = 1, k1 = k)))
  bat <- solveODEBatch(m, cs, times = tt, cores = 2L)
  for (i in seq_along(cs)) {
    ser <- solveODE(m, times = tt, parms = cs[[i]]$parms)
    expect_identical(bat[[i]]$time,     ser$time)
    expect_identical(bat[[i]]$variable, ser$variable)
    expect_identical(bat[[i]]$sens1,    ser$sens1)
  }
  expect_gt(length(bat[[1]]$time), length(tt))   # the root inserted rows
})


test_that("the reported thread count is the one actually used", {
  m <- cppODE(decay, modelname = "batch_nt", deriv = FALSE)
  expect_identical(attr(solveODEBatch(m, conds, times = tt, cores = 1L),
                        "threads"), 1L)
  # capped by the number of conditions, never above it
  nt <- attr(solveODEBatch(m, conds, times = tt, cores = 99L), "threads")
  expect_lte(nt, length(conds))
})


test_that("both backends put an event time into the output", {
  skip_if_not(isTRUE(cvodeConfig$available), "CVODE backend not available")
  # An event fires whether or not its time was requested, and the time it fires
  # at becomes an output row carrying the post-event state. The native and the
  # CVODE backend have to agree on that grid, otherwise the same model returns
  # different rows depending on `backend`.
  eqns <- c(A = "-k1 * A")
  evt  <- data.frame(var = "A", time = "t_e", value = "dose", method = "add",
                     root = NA, stringsAsFactors = FALSE)
  mn <- cppODE(eqns, events = evt, modelname = "grid_native", deriv = FALSE)
  mc <- cvode(eqns, events = evt, modelname = "grid_cvode",  deriv = FALSE)
  p  <- c(A = 1, k1 = 0.1, t_e = 25, dose = 5)
  off <- seq(0, 50, length.out = 60)          # 25 is not a grid point
  on  <- sort(unique(c(off, 25)))

  rn <- solveODE(mn, times = off, parms = p)
  rc <- solveODE(mc, times = off, parms = p)
  expect_length(rn$time, length(off) + 1L)
  expect_length(rc$time, length(off) + 1L)
  expect_true(25 %in% rn$time)
  expect_true(25 %in% rc$time)
  expect_equal(rn$time, rc$time, tolerance = 1e-12)
  # Two different integrators at their own default tolerances: the grid has to
  # match exactly, the trajectory only to solver accuracy.
  expect_equal(rn$variable, rc$variable, tolerance = 1e-4)

  # already a requested time: one row, not two
  r2 <- solveODE(mc, times = on, parms = p)
  expect_length(r2$time, length(on))
  expect_false(as.logical(anyDuplicated(r2$time)))

  # the post-event state is what lands in that row
  expect_equal(unname(rc$variable[rc$time == 25, 1]), exp(-0.1 * 25) + 5,
               tolerance = 1e-4)
})


test_that("the CVODE batch preallocates when the grid is fixed", {
  skip_if_not(isTRUE(cvodeConfig$available), "CVODE backend not available")
  eqns <- c(A = "-k1 * A", B = "k1 * A - k2 * B")
  m  <- cvode(eqns, modelname = "cv_prealloc", deriv = TRUE)
  tt <- seq(0.5, 50, length.out = 60)
  si <- diag(4); dimnames(si) <- list(NULL, c("A", "B", "k1", "k2"))
  cs <- lapply(c(0.08, 0.12, 0.2), function(k)
    list(parms = c(A = 1, B = 0, k1 = k, k2 = 0.05), sens1ini = si))
  bat <- solveODEBatch(m, cs, times = tt, cores = 3L)
  for (i in seq_along(cs)) {
    ser <- solveODE(m, times = tt, parms = cs[[i]]$parms, sens1ini = si)
    expect_identical(bat[[i]]$time,     ser$time)
    expect_identical(bat[[i]]$variable, ser$variable)
    expect_identical(bat[[i]]$sens1,    ser$sens1)
  }
  expect_identical(dimnames(bat[[1]]$variable)[[2]], c("A", "B"))
})


test_that("a root event puts its before/after pair into the output", {
  skip_if_not(isTRUE(cvodeConfig$available), "CVODE backend not available")
  # The crossing is not a requested time, so both the state just before the
  # event and the state just after it become rows. Both backends have to do it.
  evt <- data.frame(var = "A", time = NA, root = "A - 0.5", value = "0.9",
                    method = "multiply", stringsAsFactors = FALSE)
  mn <- cppODE(c(A = "-k1 * A"), events = evt, modelname = "rootpair_n", deriv = FALSE)
  mc <- cvode(c(A = "-k1 * A"),  events = evt, modelname = "rootpair_c", deriv = FALSE)
  tt <- seq(0, 50, length.out = 60)
  p  <- c(A = 1, k1 = 0.1)
  rn <- solveODE(mn, times = tt, parms = p)
  rc <- solveODE(mc, times = tt, parms = p)

  expect_length(rn$time, length(tt) + 2L)
  expect_length(rc$time, length(tt) + 2L)
  extra <- function(r) which(!round(r$time, 8) %in% round(tt, 8))
  en <- extra(rn); ec <- extra(rc)
  expect_length(en, 2L)
  expect_length(ec, 2L)
  # pre-event 0.5 (the root), post-event 0.45 (times 0.9)
  expect_equal(unname(rn$variable[en, 1]), c(0.5, 0.45), tolerance = 1e-6)
  expect_equal(unname(rc$variable[ec, 1]), c(0.5, 0.45), tolerance = 1e-4)
  expect_equal(rn$time[en], rc$time[ec], tolerance = 1e-4)
})


test_that("cvode(includeTimeZero) matches the native grid", {
  skip_if_not(isTRUE(cvodeConfig$available), "CVODE backend not available")
  mn <- cppODE(c(A = "-k1 * A"), modelname = "tz_native", deriv = FALSE)
  mc <- cvode(c(A = "-k1 * A"),  modelname = "tz_cv",     deriv = FALSE)
  m0 <- cvode(c(A = "-k1 * A"),  modelname = "tz_cv0",    deriv = FALSE,
              includeTimeZero = FALSE)
  tt <- c(1, 2, 5, 10)
  p  <- c(A = 1, k1 = 0.1)
  expect_equal(solveODE(mn, times = tt, parms = p)$time, c(0, tt))
  expect_equal(solveODE(mc, times = tt, parms = p)$time, c(0, tt))
  expect_equal(solveODE(m0, times = tt, parms = p)$time, tt)
  # the batch has to predict the injected grid, not the requested one
  b <- solveODEBatch(mc, list(list(parms = p), list(parms = c(A = 2, k1 = 0.2))),
                     times = tt, cores = 2L)
  for (i in 1:2) {
    ser <- solveODE(mc, times = tt, parms = if (i == 1) p else c(A = 2, k1 = 0.2))
    expect_identical(b[[i]]$time,     ser$time)
    expect_identical(b[[i]]$variable, ser$variable)
  }
})


test_that("the batch writes one trace file per condition", {
  m <- cppODE(c(A = "-k1 * A", B = "k1 * A - k2 * B"), modelname = "batch_trace",
              deriv = FALSE, stepTrace = TRUE)
  cs <- list(c1 = list(parms = c(A = 1, B = 0, k1 = 0.1, k2 = 0.05)),
             c2 = list(parms = c(A = 2, B = 0, k1 = 0.3, k2 = 0.05)))
  tt <- seq(0, 20, length.out = 30)
  d  <- file.path(tempdir(), "batch_trace_out")
  dir.create(d, showWarnings = FALSE)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)

  # one path used as a template: the condition name goes before the extension
  solveODEBatch(m, cs, times = tt, cores = 2L,
                traceFile = file.path(d, "trace.csv"))
  expect_setequal(basename(list.files(d)), c("trace_c1.csv", "trace_c2.csv"))
  expect_true(all(vapply(list.files(d, full.names = TRUE),
                         function(f) nrow(utils::read.csv(f)) > 0L, TRUE)))

  # or one path per condition
  unlink(list.files(d, full.names = TRUE))
  solveODEBatch(m, cs, times = tt, cores = 2L,
                traceFile = file.path(d, c("a.csv", "b.csv")))
  expect_setequal(basename(list.files(d)), c("a.csv", "b.csv"))

  # results are the same with and without a trace
  a <- solveODEBatch(m, cs, times = tt, cores = 2L)
  b <- solveODEBatch(m, cs, times = tt, cores = 2L,
                     traceFile = file.path(d, "t.csv"))
  expect_identical(lapply(a, `[[`, "variable"), lapply(b, `[[`, "variable"))
  expect_error(solveODEBatch(m, cs, times = tt, traceFile = c("a", "b", "c")),
               "one per condition")
})
