## =====================================================================
##  harness.R -- compile, solve, time and score benchmark problems.
##
##  The comparison is deliberately apples-to-apples: every solver gets
##  the same C++ right-hand side, the same analytic Jacobian, the same
##  output grid and the same tolerances.  Nothing is evaluated through
##  an R callback on either side.
## =====================================================================

bench_source <- function(dir = "benchmarks/R") {
  for (f in c("sbml.R", "tiers.R", "petab.R", "problems-classic.R",
              "shards.R", "plots.R", "report.R")) {
    p <- file.path(dir, f)
    if (file.exists(p)) source(p)
  }
  invisible(TRUE)
}


## ---------------------------------------------------------------------
##  Timing
## ---------------------------------------------------------------------

## Repeats the call until it has run for at least `min_time` seconds so
## that sub-millisecond solves are not measured as clock noise, then
## reports the median over `nrep` such batches.
##
## `min_time` has to grow when the run is parallel.  Measured on the tiny
## tier, going from 1 to 8 workers left the aggregate sens1 ratio within
## 2 % but moved the nosens ratio by 12 %, with individual cells off by
## up to 96 % -- those solves take well under a millisecond, so at 50 ms
## batches the scheduler jitter under load swamps the signal.  Averaging
## over a longer window is what buys the ratio back.
bench_time <- function(fn, nrep = 5L, min_time = 0.05) {
  fn()                                              # warm-up / page-in
  inner <- 1L
  repeat {
    t0 <- proc.time()[["elapsed"]]
    for (i in seq_len(inner)) fn()
    el <- proc.time()[["elapsed"]] - t0
    if (el >= min_time || inner >= 4096L) break
    inner <- inner * 2L
  }
  reps <- vapply(seq_len(nrep), function(k) {
    t0 <- proc.time()[["elapsed"]]
    for (i in seq_len(inner)) fn()
    (proc.time()[["elapsed"]] - t0) / inner
  }, 0)
  list(median = stats::median(reps), min = min(reps), max = max(reps),
       inner = inner, nrep = nrep)
}


## ---------------------------------------------------------------------
##  Memory guards
## ---------------------------------------------------------------------

## None of this makes a run faster.  It exists so that a run asking for
## more memory than the machine has dies in one worker instead of taking
## the machine with it.
##
## The peak is not the integration, it is the build: a sweep over a
## 124-state model with sensitivities hands g++ a very large translation
## unit, and with one worker per problem the compilers run concurrently.
## Under a scheduler the allocation bounds their sum; on a plain ssh host
## nothing does, and enough of them at once put the kernel into reclaim
## until it stops answering -- which is how a 12-core box was lost to a
## two-job, eight-worker sweep on 2026-08-01.
##
## Three guards, because none of them covers what the others do:
##
##   compile_sem_*()         how many compilers run at once, counted
##                           across workers *and* across jobs on one host
##   bench_limit_compiler()  how large one compiler may grow, so a single
##                           pathological model fails to build instead of
##                           pushing the host into swap
##   bench_limit_process()   how large one worker may grow, covering the
##                           R side of the peak (needs the `unix` package)
##
## Together they bound a host's benchmark footprint at roughly
## `compile_slots * max_compile_gb` while building, which is the number
## that has to fit next to whatever else the machine is doing.

BENCH_SEM <- new.env(parent = emptyenv())
BENCH_SEM$dir   <- NULL
BENCH_SEM$slots <- 0L

## A slot is a directory.  dir.create() is the one portable filesystem
## call that both creates and reports losing the race, so a successful
## create *is* the acquisition and nothing has to be read back.
##
## The semaphore lives under TMPDIR because that is host-local, and that
## is exactly the scope wanted: two jobs on one machine have to share the
## count, jobs on different machines must not.  Putting it in the home
## directory would do the opposite of the right thing on a cluster whose
## /home is a shared mount -- one lock there would serialise every host
## in the pool against every other.
compile_sem_init <- function(slots, id = "cppde-bench-compile") {
  slots <- suppressWarnings(as.integer(slots))
  if (is.na(slots) || slots < 1L) {
    BENCH_SEM$dir <- NULL; BENCH_SEM$slots <- 0L
    return(invisible(NULL))
  }
  user <- unname(Sys.info()[["user"]])
  if (is.na(user) || !nzchar(user)) user <- "user"
  d <- file.path(Sys.getenv("TMPDIR", unset = "/tmp"), sprintf("%s-%s", id, user))
  dir.create(d, showWarnings = FALSE, recursive = TRUE, mode = "0700")
  BENCH_SEM$dir   <- d
  BENCH_SEM$slots <- slots
  invisible(d)
}

## Blocks until a slot is free.  Returns NULL when no semaphore is
## configured, which is what makes every call site unconditional.
compile_sem_acquire <- function(poll = 0.5, stale = 3600) {
  if (is.null(BENCH_SEM$dir)) return(NULL)
  slots <- file.path(BENCH_SEM$dir,
                     sprintf("slot-%02d", seq_len(BENCH_SEM$slots)))
  repeat {
    for (s in slots) {
      if (isTRUE(suppressWarnings(dir.create(s, showWarnings = FALSE)))) {
        writeLines(c(as.character(Sys.getpid()),
                     unname(Sys.info()[["nodename"]])),
                   file.path(s, "owner"))
        return(s)
      }
    }
    for (s in slots) compile_sem_reap(s, stale)
    Sys.sleep(poll)
  }
}

compile_sem_release <- function(slot) {
  if (!is.null(slot)) unlink(slot, recursive = TRUE, force = TRUE)
  invisible(NULL)
}

## A worker killed while holding a slot -- by the OOM killer, or by the
## very limit these guards impose -- would leak it for good, and a leaked
## slot is a permanently smaller semaphore.  The owner's pid is only
## trusted when the recorded nodename matches, because the semaphore path
## is the same string on every host and a pid from another machine says
## nothing here.  The age fallback covers a host without /proc.
##
## Two reapers can in principle free the same slot and both take it; that
## needs a crash first, and the cost is one slot of overshoot, so it is
## left as a race rather than paid for with a second lock.
compile_sem_reap <- function(slot, stale = 3600) {
  if (!dir.exists(slot)) return(invisible(NULL))
  own <- tryCatch(readLines(file.path(slot, "owner"), warn = FALSE),
                  error = function(e) character(0), warning = function(w) character(0))
  dead <- length(own) >= 2L &&
    identical(own[[2L]], unname(Sys.info()[["nodename"]])) &&
    !dir.exists(file.path("/proc", own[[1L]]))
  age <- tryCatch(as.numeric(difftime(Sys.time(), file.info(slot)$mtime,
                                      units = "secs")),
                  error = function(e) 0)
  if (dead || (isTRUE(is.finite(age)) && age > stale))
    unlink(slot, recursive = TRUE, force = TRUE)
  invisible(NULL)
}

## Caps one compiler process.  cppDE builds a model through R CMD SHLIB
## (R/tools.R), which takes its compiler from R's Makeconf with
## R_MAKEVARS_USER layered on top -- so pointing the CXX* variables at a
## wrapper that calls `ulimit -v` first bounds every build the run
## starts, with no package, no privilege and no change to cppDE.
##
## The wrapper execs whatever `R CMD config` reports rather than a bare
## `g++`, because that string is not the same everywhere: this R answers
## CXX20 with the compiler alone and keeps the standard in CXX20STD,
## while other builds carry it inline as "g++ -std=gnu++17".  Passing the
## configured command through covers both; hard-coding a compiler name
## would silently drop the standard on the second kind.
##
## `ulimit -v` is RLIMIT_AS, address space rather than resident set, so
## it sits well above the true footprint; the point is a ceiling that
## stops runaway growth, not an accurate accounting -- set it generously.
## g++ reports hitting it as "virtual memory exhausted: Cannot allocate
## memory" and exits nonzero, which cppDE turns into a failed compile
## and run_all_problems() into one skipped model.
bench_limit_compiler <- function(gb, dir = tempdir()) {
  gb <- suppressWarnings(as.numeric(gb))
  if (!isTRUE(is.finite(gb)) || gb <= 0 || .Platform$OS.type != "unix")
    return(invisible(FALSE))
  ## A limit too small to start a compiler is not a smaller limit, it is
  ## a build that cannot succeed; measured here, a trivial translation
  ## unit still needs somewhere between 20 and 128 MB before g++ even
  ## reports an error.
  if (gb < 0.5)
    stop(sprintf("a %g GB compile cap cannot build anything -- no compiler starts in that", gb))
  if (gb < 2)
    warning(sprintf("a %g GB compile cap is below what a mid-size model needs", gb),
            call. = FALSE, immediate. = TRUE)
  kb   <- format(round(gb * 1024 * 1024), scientific = FALSE)
  rbin <- file.path(R.home("bin"), "R")
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)

  ## The toolchain has to be probed *unwrapped*.  `R CMD config` reads
  ## R_MAKEVARS_USER as well, so probing while an earlier call's file is
  ## still in place answers with the wrapper -- and a wrapper whose exec
  ## target is itself is an infinite exec loop that hangs the build
  ## instead of bounding it.  Unsetting first also makes repeat calls
  ## idempotent, which is what a driver that limits twice relies on.
  ## R CMD config still merges ~/.R/Makevars on its own, so a compiler the
  ## user chose there is what gets wrapped.
  old <- Sys.getenv("R_MAKEVARS_USER", unset = NA)
  Sys.unsetenv("R_MAKEVARS_USER")
  probe <- function(v) suppressWarnings(tryCatch(
    trimws(paste(system(paste(shQuote(rbin), "CMD config", v), intern = TRUE,
                        ignore.stderr = TRUE), collapse = " ")),
    error = function(e) ""))
  vars <- c("CC", "CXX", "CXX11", "CXX14", "CXX17", "CXX20", "CXX23")
  cmds <- vapply(vars, probe, "")
  if (!is.na(old)) Sys.setenv(R_MAKEVARS_USER = old)

  lines <- character(0)
  for (v in vars) {
    cmd <- cmds[[v]]
    if (!nzchar(cmd)) next          # a standard this R does not define
    if (grepl("cppde-limit-", cmd, fixed = TRUE)) {
      warning("refusing to wrap ", v, ": it already points at a limit wrapper",
              call. = FALSE, immediate. = TRUE)
      next
    }
    w <- file.path(dir, paste0("cppde-limit-", tolower(v)))
    writeLines(c("#!/bin/sh",
                 sprintf("ulimit -v %s 2>/dev/null || true", kb),
                 sprintf("exec %s \"$@\"", cmd)), w)
    Sys.chmod(w, "0755")
    lines <- c(lines, sprintf("%s = %s", v, w))
  }
  if (!length(lines)) return(invisible(FALSE))
  ## R_MAKEVARS_USER *replaces* ~/.R/Makevars rather than adding to it,
  ## so a user who has one would silently lose their flags.
  home_mk <- path.expand("~/.R/Makevars")
  if (file.exists(home_mk)) lines <- c(sprintf("include %s", home_mk), lines)
  mk <- file.path(dir, "Makevars-cppde-limit")
  writeLines(lines, mk)
  Sys.setenv(R_MAKEVARS_USER = mk)
  invisible(TRUE)
}

## Caps the R worker itself.  Base R has no setrlimit, so this one is
## best-effort: it returns NA when the `unix` package is missing, and the
## callers report that rather than failing, because the compiler cap
## above is the guard that always applies.
##
## rlimits are inherited across fork(), so setting this before mclapply()
## covers every worker and everything a worker starts.  The limit is
## per process, not per tree: eight workers under an 8 GB limit can still
## reach 64 GB between them, which is what the semaphore is for.
bench_limit_process <- function(gb) {
  gb <- suppressWarnings(as.numeric(gb))
  if (!isTRUE(is.finite(gb)) || gb <= 0) return(invisible(FALSE))
  if (!requireNamespace("unix", quietly = TRUE)) return(invisible(NA))
  invisible(tryCatch({ unix::rlimit_as(gb * 1024^3); TRUE },
                     error = function(e) FALSE))
}


## ---------------------------------------------------------------------
##  Model compilation, with caching
## ---------------------------------------------------------------------

## Compiling dominates the wall time of a benchmark run, and the same
## right-hand side is reused across tolerances and often across
## conditions, so every distinct (rhs, backend, deriv, method, fixed)
## tuple is built exactly once per session.
## `tag` keeps generated file and symbol names distinct across parallel
## workers, which otherwise race on the same .cpp in a shared build dir.
new_model_cache <- function(outdir, tag = "") {
  e <- new.env(parent = emptyenv())
  e$store <- list()
  e$outdir <- outdir
  e$counter <- 0L
  e$tag <- tag
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  e
}

cache_key <- function(...) paste(vapply(list(...), function(x)
  paste(as.character(x), collapse = "\x1f"), ""), collapse = "\x1e")

## `sparse`: NULL auto-detects, TRUE/FALSE pins the linear solver.
get_model <- function(cache, prob, backend, deriv, deriv2 = FALSE, method = "bdf",
                      useNDF = TRUE, sparse = NULL, verbose = FALSE) {
  ev <- prob$events
  key <- cache_key(paste(prob$rhs, collapse = ";"), paste(names(prob$rhs), collapse = ";"),
                   backend, deriv, deriv2, method, useNDF,
                   if (is.null(sparse)) "auto" else sparse,
                   paste(sort(prob$fixed), collapse = ";"),
                   if (is.null(ev)) "" else paste(unlist(ev), collapse = ";"))
  if (!is.null(cache$store[[key]])) return(cache$store[[key]])

  cache$counter <- cache$counter + 1L
  nm <- sprintf("bm%s%04d_%s", cache$tag %||% "", cache$counter,
                substr(prob$id, 1, 24))
  fixed <- if (deriv) prob$fixed else NULL

  ## Only the build is metered, and only the build waits: a cached model
  ## has already returned above, and the timing starts after the slot is
  ## in hand so that queueing never lands in compile_seconds.
  slot <- compile_sem_acquire()
  on.exit(compile_sem_release(slot), add = TRUE)

  t0 <- proc.time()[["elapsed"]]
  m <- if (backend == "cvode") {
    if (deriv2) stop("CVODES has no second-order sensitivities")
    cvode(prob$rhs, events = ev, deriv = deriv, fixed = fixed, outdir = cache$outdir,
          modelname = nm, method = method, sparse = sparse,
          compile = TRUE, verbose = verbose)
  } else {
    ## deriv2 needs a compile-time AD width; NULL selects exactly the
    ## number of active sensitivity parameters.
    cppODE(prob$rhs, events = ev, deriv = deriv, deriv2 = deriv2, fixed = fixed,
           nStack = if (deriv2) NULL else Inf,
           outdir = cache$outdir, sparse = sparse,
           modelname = nm, method = method, useNDF = useNDF,
           compile = TRUE, verbose = verbose)
  }
  attr(m, "compile_seconds") <- proc.time()[["elapsed"]] - t0
  cache$store[[key]] <- m
  m
}


## R caps the number of simultaneously loaded DLLs (100 by default), and
## a full sweep compiles several per problem.  Releasing a problem's
## models before moving on keeps the run within that budget.  Names are
## never reused, so nothing is ever re-loaded under an old name.
release_cache <- function(cache) {
  for (m in cache$store) {
    so <- sub("\\.cpp$", .Platform$dynlib.ext, attr(m, "srcfile"))
    if (!is.null(so) && file.exists(so)) try(dyn.unload(so), silent = TRUE)
  }
  cache$store <- list()
  invisible(cache)
}


## ---------------------------------------------------------------------
##  Accuracy
## ---------------------------------------------------------------------

## cppDE reports an extra output row at each event time while CVODE
## reports only the requested grid, so results are matched on the time
## column before being compared rather than by row position.
align_rows <- function(res, times) {
  if (is.null(res) || is.null(res$time)) return(NULL)
  idx <- match(times, res$time)
  keep <- !is.na(idx)
  if (!any(keep)) return(NULL)
  idx <- idx[keep]
  list(time = res$time[idx],
       variable = res$variable[idx, , drop = FALSE],
       sens1 = if (!is.null(res$sens1)) res$sens1[idx, , , drop = FALSE] else NULL)
}

## Each state is normalised by its own trajectory magnitude, so a
## species living at 1e-9 counts as much as one at 1e5 -- without that,
## the error of a stiff biological model is decided by one large state.
traj_error <- function(x, ref) {
  if (is.null(x) || is.null(ref)) return(NA_real_)
  n <- min(nrow(x), nrow(ref))
  if (n < 2L || ncol(x) != ncol(ref)) return(NA_real_)
  x <- x[seq_len(n), , drop = FALSE]; ref <- ref[seq_len(n), , drop = FALSE]
  scale <- apply(abs(ref), 2L, max, na.rm = TRUE)
  scale[!is.finite(scale) | scale <= 0] <- 1
  err <- abs(x - ref) / rep(scale, each = n)
  if (!any(is.finite(err))) return(NA_real_)
  max(err[is.finite(err)])
}

## Same normalisation, applied over the [time, state, parameter] cube.
sens_error <- function(s, ref) {
  if (is.null(s) || is.null(ref)) return(NA_real_)
  if (!identical(dim(s), dim(ref))) return(NA_real_)
  scale <- apply(abs(ref), 2L, max, na.rm = TRUE)
  scale[!is.finite(scale) | scale <= 0] <- 1
  err <- abs(s - ref) / array(rep(scale, each = dim(s)[1L]), dim(s))
  if (!any(is.finite(err))) return(NA_real_)
  max(err[is.finite(err)])
}


## ---------------------------------------------------------------------
##  Solver configurations under test
## ---------------------------------------------------------------------

## The head-to-head the suite is built around: cppDE's own NDF/BDF
## multistep solver against SUNDIALS CVODE(S), which is what the R
## ecosystem otherwise reaches for.  `extra = TRUE` adds cppDE's plain
## BDF (CVODE's own formula, so differences isolate the NDF change) and
## the Rosenbrock4 one-step solver.
solver_configs <- function(extra = FALSE) {
  cfg <- list(
    list(label = "cppDE_ndf", backend = "cppde", method = "bdf", useNDF = TRUE),
    list(label = "CVODE_bdf",  backend = "cvode",  method = "bdf", useNDF = NA)
  )
  if (extra) cfg <- c(cfg, list(
    list(label = "cppDE_bdf", backend = "cppde", method = "bdf", useNDF = FALSE),
    list(label = "cppDE_rb4", backend = "cppde", method = "rb4", useNDF = NA)
  ))
  cfg
}


## Every backend twice, with the linear solver pinned dense and sparse
## instead of auto-detected.  These run *alongside* solver_configs() on
## the models marked by mark_sparse_sweep(), so one run answers both "is
## cppDE faster than CVODE" and "was auto-detection right".
sparse_sweep_configs <- function() {
  list(
    list(label = "cppDE_dense",  backend = "cppde", method = "bdf",
         useNDF = TRUE, sparse = FALSE),
    list(label = "cppDE_sparse", backend = "cppde", method = "bdf",
         useNDF = TRUE, sparse = TRUE),
    list(label = "CVODE_dense",   backend = "cvode",  method = "bdf",
         useNDF = NA,   sparse = FALSE),
    list(label = "CVODE_sparse",  backend = "cvode",  method = "bdf",
         useNDF = NA,   sparse = TRUE)
  )
}


## What a config asked of the linear solver, which is not the same as the
## `lu` column: `lu` is what the model compiled to, `pinned` is what was
## requested.  Separating them is what lets the dense/sparse comparison
## key on the pinned pairs and leave the auto-detected rows -- the
## head-to-head -- out of it.
pin_label <- function(sparse)
  if (is.null(sparse)) "auto" else if (isTRUE(sparse)) "sparse" else "dense"


## Marks the cases whose linear solver is worth pinning both ways, and
## reports the ones it leaves out.  Two gates, for two different reasons:
##
##   nstates >= min_states  below this the codegen never auto-selects
##                          sparse (decide_sparse() in codegen_cppODE.py),
##                          so both pinnings would measure a decision
##                          that is not taken
##   density  <= max_density  a structurally dense Jacobian has no sparse
##                          path to measure; the proxy overestimates, so
##                          a model kept here may still turn out dense
##
## Raising --max-density past the codegen's own 0.4 cutoff is what turns
## the sweep from "does sparse pay where we use it" into "is 0.4 the
## right place to switch".
mark_sparse_sweep <- function(problems, max_density = 0.25, min_states = 8L,
                              verbose = TRUE) {
  kept <- character(0)
  for (pname in names(problems)) {
    p0 <- problems[[pname]][[1L]]
    d  <- jac_density_proxy(p0)
    why <- if (p0$nstates < min_states)
             sprintf("%d states < %d", p0$nstates, min_states)
           else if (is.na(d)) "Jacobian density unknown"
           else if (d > max_density)
             sprintf("Jacobian density %.2f > %.2f", d, max_density)
           else NA_character_
    ok <- is.na(why)
    if (ok) kept <- c(kept, pname)
    else if (verbose)
      cat(sprintf("  [no sweep] %-32s %s\n", pname, why))
    for (i in seq_along(problems[[pname]])) problems[[pname]][[i]]$sweep <- ok
  }
  if (verbose)
    cat(sprintf("  sparse sweep: %d of %d model(s)\n", length(kept), length(problems)))
  attr(problems, "sweep_models") <- kept
  problems
}


## ---------------------------------------------------------------------
##  One problem, all solvers x modes x tolerances
## ---------------------------------------------------------------------

## `ref_tol` must be strictly tighter than every benchmarked tolerance.
## It is produced by CVODE, so that cppDE is never scored against
## itself -- but that cuts both ways: if the reference used the same
## setting as the tightest swept cell, CVODE would be compared against
## itself there, score an error of exactly zero, and drop out of the
## work-precision plot as if it could not reach that accuracy.  Hence
## two decades of headroom below the tightest sweep (atol 1e-13).
## `max_sens2` is separate from the first-order cap because the cost of
## second-order forward AD grows with M^2: a width that is routine for
## sens1 makes sens2 the whole run.
##
## `sweep_configs` are appended for a case that mark_sparse_sweep() kept,
## so the dense/sparse comparison rides along in the same run and on the
## same worker as the head-to-head it has to be read next to.
run_problem <- function(prob, cache, tolerances, modes = c("nosens", "sens1"),
                        configs = solver_configs(), nrep = 5L,
                        ref_tol = c(atol = 1e-14, rtol = 1e-12),
                        verbose = TRUE, ref_backend = "cvode",
                        max_sens2 = 10L, min_time = 0.05,
                        sweep_configs = NULL) {
  if (isTRUE(prob$sweep) && length(sweep_configs))
    configs <- c(configs, sweep_configs)
  rows <- list()
  say <- function(...) if (verbose) cat(...)
  emit <- function(...) rows[[length(rows) + 1L]] <<- data.frame(..., stringsAsFactors = FALSE)

  say(sprintf("\n=== %s%s | %d states, %d params, %d sens, %d out ===\n",
              prob$name,
              if (is.null(prob$condition) || is.na(prob$condition)) ""
              else paste0(" [", prob$condition, "]"),
              prob$nstates, prob$npars, prob$nsens, length(prob$times)))

  ## Second order is a cppDE-only capability -- CVODES does not provide
  ## it -- so there is nothing to compare against and no reference is
  ## computed: the sens2 rows report cost only.
  narrow_for_sens2 <- function(p) {
    if (length(p$sens) <= max_sens2) return(p)
    p$sens  <- p$sens[seq_len(max_sens2)]
    p$fixed <- setdiff(names(p$parms), p$sens)
    p$nsens <- length(p$sens)
    p
  }

  ## -- reference trajectories -------------------------------------------
  ref <- list()
  for (mode in setdiff(modes, "sens2")) {
    deriv <- identical(mode, "sens1")
    r <- tryCatch({
      m <- get_model(cache, prob, ref_backend, deriv = deriv, method = "bdf")
      align_rows(solveODE(m, prob$times, prob$parms,
                          abstol = min(ref_tol[["atol"]], prob$atol %||% Inf),
                          reltol = ref_tol[["rtol"]], onFailure = "stop"),
                 prob$times)
    }, error = function(e) { say("  reference (", mode, ") failed: ",
                                 sub("\n.*", "", conditionMessage(e)), "\n"); NULL })
    ref[[mode]] <- r
  }

  for (mode in modes) {
    deriv2 <- identical(mode, "sens2")
    deriv  <- deriv2 || identical(mode, "sens1")
    prob_m <- if (deriv2) narrow_for_sens2(prob) else prob
    if (deriv && prob_m$nsens == 0L) {
      say(sprintf("  [%s skipped: no free parameters]\n", mode)); next
    }

    for (cfg in configs) {
      if (deriv2 && cfg$backend == "cvode") next   # CVODES cannot do this
      m <- tryCatch(get_model(cache, prob_m, cfg$backend, deriv = deriv,
                              deriv2 = deriv2, method = cfg$method,
                              useNDF = if (is.na(cfg$useNDF)) TRUE else cfg$useNDF,
                              sparse = cfg$sparse),
                    error = function(e) { say(sprintf("  %-11s %-6s COMPILE FAILED: %s\n",
                                              cfg$label, mode,
                                              sub("\n.*", "", conditionMessage(e)))); NULL })
      if (is.null(m)) next
      ## What the model compiled to, not what was requested.
      lu <- if (isTRUE(attr(m, "sparse"))) "sparse" else "dense"
      pinned <- pin_label(cfg$sparse)

      for (ti in seq_len(nrow(tolerances))) {
        rtol <- tolerances$rtol[ti]
        ## A problem may pin its absolute tolerance; the relative one
        ## is still swept, which is how the IVP test set treats E5.
        atol <- prob$atol %||% tolerances$atol[ti]
        run <- function() solveODE(m, prob$times, prob$parms, abstol = atol,
                                   reltol = rtol, onFailure = "silent")
        res <- tryCatch(run(), error = function(e) NULL)
        ok <- !is.null(res) && !is.null(res$diagnostics) &&
              res$diagnostics$return_code == 0L
        if (!ok) {
          say(sprintf("  %-11s %-6s rtol=%-7.0e FAILED\n", cfg$label, mode, rtol))
          emit(problem = prob$name, condition = prob$condition %||% NA_character_,
               source = prob$source, nstates = prob$nstates, npars = prob$npars,
               nsens = if (deriv) prob_m$nsens else 0L, nout = length(prob$times),
               solver = cfg$label, backend = cfg$backend, lu = lu,
               pinned = pinned, mode = mode,
               atol = atol, rtol = rtol, ok = FALSE, time_ms = NA_real_,
               time_min_ms = NA_real_, time_max_ms = NA_real_,
               solves_per_batch = NA_integer_, batches = NA_integer_,
               accepted = NA_integer_, rejected = NA_integer_, fevals = NA_integer_,
               jevals = NA_integer_, setups = NA_integer_,
               err = NA_real_, err_sens = NA_real_,
               compile_s = unname(attr(m, "compile_seconds")))
          next
        }
        tm <- bench_time(run, nrep = nrep, min_time = min_time)
        ## sens2 has no reference (nothing else computes Hessians here),
        ## so its rows carry timing and step counts only.
        al  <- if (deriv2) NULL else align_rows(res, prob$times)
        e_x <- if (deriv2) NA_real_ else traj_error(al$variable, ref[[mode]]$variable)
        e_s <- if (deriv2 || !deriv) NA_real_
               else sens_error(al$sens1, ref[[mode]]$sens1)
        d <- res$diagnostics
        say(sprintf("  %-11s %-6s rtol=%-7.0e %8.2f ms  steps=%-6s fev=%-7s%s\n",
                    cfg$label, mode, rtol, tm$median * 1000,
                    d$accepted %||% NA, d$fevals %||% NA,
                    if (deriv2) sprintf("  M=%d", prob_m$nsens)
                    else sprintf(" err=%.1e%s", e_x,
                                 if (deriv) sprintf(" errS=%.1e", e_s) else "")))
        emit(problem = prob$name, condition = prob$condition %||% NA_character_,
             source = prob$source, nstates = prob$nstates, npars = prob$npars,
             nsens = if (deriv) prob_m$nsens else 0L, nout = length(prob$times),
             solver = cfg$label, backend = cfg$backend, lu = lu,
             pinned = pinned, mode = mode,
             atol = atol, rtol = rtol, ok = TRUE, time_ms = tm$median * 1000,
             time_min_ms = tm$min * 1000, time_max_ms = tm$max * 1000,
             solves_per_batch = tm$inner, batches = tm$nrep,
             accepted = d$accepted %||% NA_integer_, rejected = d$rejected %||% NA_integer_,
             fevals = d$fevals %||% NA_integer_, jevals = d$jevals %||% NA_integer_,
             setups = d$setups %||% NA_integer_,
             err = e_x, err_sens = e_s,
             compile_s = unname(attr(m, "compile_seconds")))
      }
    }
  }
  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}

`%||%` <- function(a, b) if (is.null(a)) b else a


## ---------------------------------------------------------------------
##  Running many problems, optionally in parallel
## ---------------------------------------------------------------------

## Parallelism is applied at the granularity of a whole *problem*: one
## worker runs every solver, mode and tolerance of that problem, in
## sequence.  That is what keeps the results meaningful.
##
## Wall-clock measurements taken under load are inflated -- shared L3,
## memory bandwidth, and above all turbo clocking, which drops the core
## frequency once many cores are busy.  Those effects are largely
## *uniform*, so the suite's headline number, a ratio between two solvers
## on the same cell, survives them: both solvers of a matched cell ran in
## the same worker under the same neighbourhood.  Splitting a problem's
## cells across workers would break exactly that property.
##
## Absolute milliseconds from a run with cores > 1 are therefore not
## comparable with those from a serial run.  The `cores` column in
## results.csv records which it was, so the two are never mixed silently.
run_all_problems <- function(problems, builddir, tolerances, modes, configs,
                             nrep, cores = 1L, max_sens2 = 10L,
                             max_sens = Inf, max_states = Inf, min_time = NULL,
                             compile_slots = 0L, sweep_configs = NULL,
                             on_skip = function(name, why) invisible()) {
  ## Longer batches under load -- see the note on bench_time().
  if (is.null(min_time)) min_time <- if (cores > 1L) 0.25 else 0.05

  ## Set up before the fork so the workers inherit the configuration, and
  ## torn down after so a serial caller in the same session is not left
  ## with a semaphore it never asked for.  The slots themselves outlive
  ## the process by design -- that is what lets a second job on the same
  ## host count against the same ceiling.
  compile_sem_init(compile_slots)
  on.exit(compile_sem_init(0L), add = TRUE)

  keep <- vapply(problems, function(cases) {
    p0 <- cases[[1L]]
    if (p0$nstates > max_states) {
      on_skip(p0$name, sprintf("%d states > max-states %d", p0$nstates, max_states))
      return(FALSE)
    }
    TRUE
  }, NA)
  problems <- problems[keep]
  if (!length(problems)) return(NULL)

  one <- function(i) {
    pname <- names(problems)[i]
    cache <- new_model_cache(file.path(builddir, sprintf("w%02d", i)),
                             tag = sprintf("%02d", i %% 100L))
    on.exit(release_cache(cache), add = TRUE)
    out <- list()
    for (case in problems[[i]]) {
      if (length(case$sens) > max_sens) {
        case$sens  <- case$sens[seq_len(max_sens)]
        case$fixed <- setdiff(names(case$parms), case$sens)
        case$nsens <- length(case$sens)
      }
      r <- tryCatch(
        run_problem(case, cache, tolerances = tolerances, modes = modes,
                    configs = configs, nrep = nrep, max_sens2 = max_sens2,
                    verbose = cores == 1L, min_time = min_time,
                    sweep_configs = sweep_configs),
        error = function(e) { message("  [error] ", case$name, ": ",
                                      sub("\n.*", "", conditionMessage(e))); NULL })
      if (!is.null(r)) out[[length(out) + 1L]] <- r
    }
    if (!length(out)) return(NULL)
    do.call(rbind, out)
  }

  res <- if (cores > 1L && .Platform$OS.type == "unix") {
    ## Nested threading would silently oversubscribe every worker.
    old <- Sys.getenv(c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                        "MKL_NUM_THREADS"), unset = NA)
    Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
               MKL_NUM_THREADS = "1")
    ## Sys.setenv() takes named *arguments*, not a named vector, so the
    ## restore has to go through do.call.  The bug hides on a machine
    ## where these variables were unset to begin with -- only the
    ## unsetenv branch runs there -- and surfaces on a cluster, where
    ## the module system sets OMP_NUM_THREADS before R starts.
    on.exit({
      for (n in names(old)) {
        if (is.na(old[[n]])) Sys.unsetenv(n)
        else do.call(Sys.setenv, stats::setNames(list(old[[n]]), n))
      }
    }, add = TRUE)
    parallel::mclapply(seq_along(problems), function(i)
      tryCatch(one(i), error = function(e) NULL),
      mc.cores = cores, mc.preschedule = FALSE)
  } else {
    lapply(seq_along(problems), one)
  }

  res <- Filter(function(x) is.data.frame(x) && nrow(x), res)
  if (!length(res)) return(NULL)
  df <- do.call(rbind, res)
  df$cores <- cores
  df
}


## ---------------------------------------------------------------------
##  Summaries
## ---------------------------------------------------------------------

## The `pinned` column, or what it would have been for a results.csv
## written before it existed: there, the sweep was a separate run whose
## solver labels carry the pinning, and everything else was auto.
pinned_col <- function(df) {
  if ("pinned" %in% names(df)) return(df$pinned)
  ifelse(grepl("_dense$", df$solver), "dense",
         ifelse(grepl("_sparse$", df$solver), "sparse", "auto"))
}


## Speed-up of every solver against CVODE on the identical
## (problem, condition, mode, tolerance) cell.  Ratios are combined with
## a geometric mean, which is the correct average for relative measures.
##
## Restricted to the auto-detected rows.  The sweep's pinned rows are the
## same two backends a second and third time, and letting them in would
## count those models repeatedly in a mean that is supposed to weight
## every model once.
speedup_table <- function(df, baseline = "CVODE_bdf") {
  df <- df[pinned_col(df) == "auto", , drop = FALSE]
  if (!nrow(df)) return(NULL)
  ok <- df[df$ok & is.finite(df$time_ms), ]
  if (!nrow(ok)) return(NULL)
  base <- ok[ok$solver == baseline,
             c("problem", "condition", "mode", "atol", "rtol", "time_ms", "fevals")]
  if (!nrow(base)) return(NULL)
  names(base)[names(base) == "time_ms"] <- "base_ms"
  names(base)[names(base) == "fevals"]  <- "base_fev"
  m <- merge(ok, base, by = c("problem", "condition", "mode", "atol", "rtol"))
  m$speedup <- m$base_ms / m$time_ms
  m$fev_ratio <- m$fevals / m$base_fev
  m
}

## Dense/sparse time ratio per backend and cell.
##
## Keyed on `pinned`, not on `lu`: the auto-detected head-to-head rows
## also carry an `lu`, and pairing one of those with a pinned row would
## compare a model against itself.
sparse_gain_table <- function(df) {
  if (!nrow(df) || !"lu" %in% names(df)) return(NULL)
  df <- df[pinned_col(df) %in% c("dense", "sparse"), , drop = FALSE]
  if (!nrow(df)) return(NULL)
  ok <- df[df$ok & is.finite(df$time_ms), ]
  ok$pinned <- pinned_col(ok)
  key <- c("problem", "condition", "mode", "atol", "rtol", "backend")
  d <- ok[ok$pinned == "dense",  c(key, "nstates", "nsens", "time_ms")]
  s <- ok[ok$pinned == "sparse", c(key, "time_ms")]
  if (!nrow(d) || !nrow(s)) return(NULL)
  names(d)[names(d) == "time_ms"] <- "dense_ms"
  names(s)[names(s) == "time_ms"] <- "sparse_ms"
  m <- merge(d, s, by = key)
  if (!nrow(m)) return(NULL)
  m$gain <- m$dense_ms / m$sparse_ms
  m[order(m$backend, m$nstates, m$problem), ]
}

## TRUE when some cell exists in both a dense and a sparse variant.
is_sparse_sweep <- function(df) {
  g <- sparse_gain_table(df)
  !is.null(g) && nrow(g) > 0L
}


geo_mean <- function(x) {
  x <- x[is.finite(x) & x > 0]
  if (!length(x)) return(NA_real_)
  exp(mean(log(x)))
}
