## =====================================================================
##  harness.R: compile, solve, time and score benchmark problems.
## =====================================================================

##  Every solver gets the same C++ right-hand side, the same analytic Jacobian,
##  the same output grid and the same tolerances. Nothing runs through an R
##  callback on either side.

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

## Repeats the call until it has run for at least `min_time` seconds, so a
## sub-millisecond solve is not measured as clock noise, then reports the
## median over `nrep` such batches.

## `min_time` has to grow with the worker count: at 50 ms batches the scheduler
## jitter under load swamps a sub-millisecond solve, and only a longer window
## brings the ratio back.
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

## None of this makes a run faster. It bounds the peak, which is the build and
## not the integration: one worker per problem runs the compilers concurrently,
## and enough at once put the kernel into reclaim until the host stops answering.

## Three guards, none covering what the others do: compile_sem_*() how many
## compilers run at once across workers and jobs, bench_limit_compiler() how
## large one compiler may grow, bench_limit_process() how large one worker may.

## Together they bound a host's footprint at roughly
## `compile_slots * max_compile_gb` while building.

BENCH_SEM <- new.env(parent = emptyenv())
BENCH_SEM$dir   <- NULL
BENCH_SEM$slots <- 0L

## A slot is a directory: dir.create() is the one portable call that both
## creates and reports losing the race, so a successful create is the
## acquisition and nothing has to be read back.

## The semaphore lives under TMPDIR because that is host-local, which is the
## scope wanted: two jobs on one machine share the count, jobs on different
## machines must not. A shared /home would serialise the whole pool.
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

## A worker killed while holding a slot would leak it for good, and a leaked
## slot is a permanently smaller semaphore. The owner's pid is trusted only
## when the nodename matches; the age fallback covers a host without /proc.

## Two reapers can free the same slot and both take it. That needs a crash
## first and costs one slot of overshoot, so it is left as a race.
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

## Caps one compiler process. R CMD SHLIB takes its compiler from Makeconf with
## R_MAKEVARS_USER on top, so pointing CXX* at a wrapper that calls `ulimit -v`
## bounds every build the run starts, with no privilege and no change to cppDE.

## The wrapper execs what `R CMD config` reports, not a bare `g++`: some builds
## answer CXX20 with the compiler alone and keep the standard in CXX20STD,
## others carry it inline, and hard-coding would drop it on the latter.

## `ulimit -v` is RLIMIT_AS, address space rather than resident set, so set it
## generously: it is a ceiling against runaway growth, not an accounting. g++
## exits nonzero on it, which becomes a failed compile and one skipped model.
bench_limit_compiler <- function(gb, dir = tempdir()) {
  gb <- suppressWarnings(as.numeric(gb))
  if (!isTRUE(is.finite(gb)) || gb <= 0 || .Platform$OS.type != "unix")
    return(invisible(FALSE))
  ## A limit too small to start a compiler is not a smaller limit, it is a build
  ## that cannot succeed: a trivial translation unit needs tens of MB before g++
  ## even reports an error.
  if (gb < 0.5)
    stop(sprintf("a %g GB compile cap cannot build anything, no compiler starts in that", gb))
  if (gb < 2)
    warning(sprintf("a %g GB compile cap is below what a mid-size model needs", gb),
            call. = FALSE, immediate. = TRUE)
  kb   <- format(round(gb * 1024 * 1024), scientific = FALSE)
  rbin <- file.path(R.home("bin"), "R")
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)

  ## The toolchain has to be probed unwrapped: `R CMD config` reads
  ## R_MAKEVARS_USER too, and a wrapper whose exec target is itself is an infinite
  ## exec loop. Unsetting first also makes a repeated limit call idempotent.

  ## R CMD config still merges ~/.R/Makevars, so a compiler chosen there is the
  ## one that gets wrapped.
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

## Caps the R worker itself, best-effort: base R has no setrlimit, so this
## returns NA without the `unix` package and callers report that rather than
## fail, the compiler cap above being the guard that always applies.

## rlimits are inherited across fork(), so setting this before mclapply() covers
## every worker. It is per process, not per tree: eight workers under 8 GB can
## still reach 64 GB between them, which is what the semaphore is for.
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

## Compiling dominates a run's wall time and the same right-hand side is reused
## across tolerances and conditions, so every distinct (rhs, backend, deriv,
## method, fixed) tuple is built once per session.

## `tag` keeps generated file and symbol names distinct across parallel workers,
## which otherwise race on the same .cpp in a shared build dir.
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


## R caps the number of simultaneously loaded DLLs, and a full sweep compiles
## several per problem, so a problem's models are released before the next.
## Names are never reused, so nothing is re-loaded under an old name.
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
## species living at 1e-9 counts as much as one at 1e5, without that,
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

## The head-to-head the suite is built around: cppDE's NDF/BDF multistep solver
## against SUNDIALS CVODE(S). `extra = TRUE` adds cppDE's plain BDF, CVODE's own
## formula, so a difference isolates the NDF change, and Rosenbrock4.
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


## Every backend twice, with the linear solver pinned dense and sparse instead
## of auto-detected. These run alongside solver_configs() on the models marked
## by mark_sparse_sweep(), so one run answers both questions.
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


## What a config asked of the linear solver: `lu` is what the model compiled to,
## `pinned` what was requested. Separating them lets the dense/sparse comparison
## key on the pinned pairs and leave the auto-detected head-to-head rows out.
pin_label <- function(sparse)
  if (is.null(sparse)) "auto" else if (isTRUE(sparse)) "sparse" else "dense"


## Marks the cases whose linear solver is worth pinning both ways, and reports
## the ones left out: below `min_states` the codegen never auto-selects sparse,
## and above `max_density` there is no sparse path to measure.

## Raising --max-density past the codegen's own 0.4 cutoff turns the sweep from
## "does sparse pay where we use it" into "is 0.4 the right place to switch".
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

## `ref_tol` must be strictly tighter than every benchmarked tolerance, two
## decades below the tightest sweep: at equal settings CVODE would score exactly
## zero error against itself and drop out of the work-precision plot.

## `max_sens2` is separate from the first-order cap because second-order forward
## AD grows with M^2, so a width that is routine for sens1 makes sens2 the run.

## `sweep_configs` are appended for a case mark_sparse_sweep() kept, so the
## dense/sparse comparison rides in the same run and worker as the head-to-head.
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

  ## Second order is a cppDE-only capability, CVODES does not provide
  ## it, so there is nothing to compare against and no reference is
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

## Parallelism is applied per problem: one worker runs every solver, mode and
## tolerance of that problem in sequence, so both solvers of a matched cell run
## in the same worker under the same neighbourhood.

## Wall-clock under load is inflated by shared L3, memory bandwidth and turbo
## clocking, but largely uniformly, so a ratio survives it and an absolute
## millisecond does not. The `cores` column records which run it came from.
run_all_problems <- function(problems, builddir, tolerances, modes, configs,
                             nrep, cores = 1L, max_sens2 = 10L,
                             max_sens = Inf, max_states = Inf, min_time = NULL,
                             compile_slots = 0L, sweep_configs = NULL,
                             on_skip = function(name, why) invisible()) {
  ## Longer batches under load, see the note on bench_time().
  if (is.null(min_time)) min_time <- if (cores > 1L) 0.25 else 0.05

  ## Set up before the fork so the workers inherit it, torn down after so a serial
  ## caller in the same session keeps no semaphore it never asked for. The slots
  ## outlive the process, which is what lets a second job share the ceiling.
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
    ## Sys.setenv() takes named arguments, not a named vector, so the restore goes
    ## through do.call. Only the unsetenv branch runs where these variables were
    ## unset to begin with, so a cluster module setting OMP_NUM_THREADS exposes it.
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


## Speed-up of every solver against CVODE on the identical (problem, condition,
## mode, tolerance) cell, combined with a geometric mean, the correct average
## for relative measures.

## Restricted to the auto-detected rows: the sweep's pinned rows are the same
## two backends again, and would count those models repeatedly in a mean meant
## to weight every model once.
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

## Dense/sparse time ratio per backend and cell. Keyed on `pinned`, not on
## `lu`: the auto-detected head-to-head rows also carry an `lu`, and pairing
## one of those with a pinned row would compare a model against itself.
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
