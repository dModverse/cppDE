#!/usr/bin/env Rscript

## =====================================================================
##  run-benchmarks.R: cppDE versus SUNDIALS CVODE(S).
## =====================================================================

##      Rscript benchmarks/run-benchmarks.R --quick
##      Rscript benchmarks/run-benchmarks.R --suite all --conditions all
##      Rscript benchmarks/run-benchmarks.R --models Boehm,Fujita --tol wp

##  Run `--help` for the full option list. Results land in
##  benchmarks/results/<timestamp>/ as results.csv plus six figures.

suppressPackageStartupMessages({
  library(cppDE)
})

ROOT <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) normalizePath(file.path(dirname(f[1]), "..")) else getwd()
})

source(file.path(ROOT, "benchmarks", "R", "harness.R"))
bench_source(file.path(ROOT, "benchmarks", "R"))


## ---------------------------------------------------------------------
##  Options
## ---------------------------------------------------------------------

OPTS <- list(
  tier         = "medium",    # tiny | medium | full
  suite        = "all",       # all | petab | classic
  models       = "",          # comma-separated substrings; "" = every model
  skip         = "",          # comma-separated substrings to leave out
  conditions   = "1",         # integer, or "all"
  modes        = "nosens,sens1",
  `max-sens2`  = "10",
  tol          = "default",   # default | wp | loose | tight
  nrep         = "5",
  `min-time`   = "",           # seconds per timing batch; "" = auto
  cores        = "1",
  ## Memory guards (R/harness.R).  Serial runs never trip them; they
  ## start to matter at --cores > 1, where one compiler per worker is
  ## what a large model actually costs.
  `compile-slots` = "",        # "" = min(cores, 4); 0 disables
  `max-compile-gb` = "8",      # ulimit -v per compiler process
  `max-worker-gb` = "",        # RLIMIT_AS per worker; needs `unix`
  shard        = "",           # "i/n": run only every n-th problem
  `max-states` = "400",
  `max-sens`   = "32",
  `min-points` = "25",
  outdir       = file.path(ROOT, "benchmarks", "results"),
  `petab-root` = file.path(ROOT, "benchmarks", "cache", "petab", "Benchmark-Models"),
  builddir     = "",          # default: a fresh temp directory
  `extra-solvers` = "FALSE",
  `sparse-sweep` = "TRUE",
  `max-density` = "0.25",
  `min-sweep-states` = "8",
  `include-excluded` = "FALSE",
  plots        = "TRUE",
  quick        = "FALSE",
  help         = "FALSE")

## Returns the option list with a `.set` attribute naming the options the
## caller passed explicitly, so a tier's defaults never override them.
parse_args <- function(args, opts) {
  seen <- character(0)
  i <- 1L
  while (i <= length(args)) {
    a <- args[i]
    if (!grepl("^--", a)) { warning("ignoring argument: ", a); i <- i + 1L; next }
    key <- sub("^--", "", a)
    if (grepl("=", key)) {
      kv <- strsplit(key, "=", fixed = TRUE)[[1L]]
      key <- kv[1L]; val <- paste(kv[-1L], collapse = "=")
    } else if (key %in% c("quick", "help", "extra-solvers", "sparse-sweep",
                          "include-excluded", "no-plots", "no-sparse-sweep")) {
      val <- "TRUE"
    } else {
      i <- i + 1L; val <- if (i <= length(args)) args[i] else ""
    }
    if (key == "no-plots") { opts$plots <- "FALSE"; seen <- c(seen, "plots") }
    else if (key == "no-sparse-sweep") {
      opts$`sparse-sweep` <- "FALSE"; seen <- c(seen, "sparse-sweep")
    }
    else if (!key %in% names(opts)) stop("unknown option --", key)
    else { opts[[key]] <- val; seen <- c(seen, key) }
    i <- i + 1L
  }
  attr(opts, ".set") <- unique(seen)
  opts
}

OPTS <- parse_args(commandArgs(trailingOnly = TRUE), OPTS)
user_opts <- stats::setNames(vector("list", length(attr(OPTS, ".set"))),
                             attr(OPTS, ".set"))

if (isTRUE(as.logical(OPTS$help))) {
  cat("
cppDE benchmark suite, cppDE vs SUNDIALS CVODE(S)

  --tier <tiny|medium|full>     benchmark depth                  [medium]
        tiny    ~5 min,  one representative of every problem trait,
                the basis for judging a small change
        medium  ~30 min, at least two per trait, wider size range
        full    everything the size limit allows
  --suite <all|petab|classic>   which problem sources to run     [all]
  --models <a,b,c>              only models whose name contains one of these
  --skip <a,b,c>                leave out models whose name contains one of
                                these; applied after --models
  --conditions <n|all>          PEtab conditions per model       [1]
  --modes <nosens,sens1,sens2>  integration modes                [nosens,sens1]
        sens2 measures second-order (Hessian) runtimes.  CVODES has no
        second-order sensitivities, so those rows are cppDE-only and
        carry timings without a cross-implementation comparison.
  --max-sens2 <n>               cap on sensitivity parameters for sens2,
                                whose cost grows with M^2            [10]
  --tol <default|wp|loose|tight>  tolerance set; 'wp' is a 5-point sweep
                                for work-precision diagrams      [default]
  --nrep <n>                    timing batches per cell          [5]
  --min-time <s>                seconds per batch; the solve is repeated
                                until a batch lasts this long, so short
                                solves are averaged instead of measured
                                once.  Default 0.05 serial, 0.25 with
                                --cores > 1.  The reported time is the
                                median over the batches; results.csv
                                records solves_per_batch and batches.
  --shard <i/n>                 split the problem list into n shards and
                                run only shard i (1-based).  Problems are
                                assigned round-robin over a deterministic
                                ordering, so the shards are reproducible
                                and cover the list exactly once.  Intended
                                for cluster array jobs; merge the shards'
                                results.csv afterwards.  Absolute times are
                                only comparable within a shard unless every
                                shard ran on identical hardware.
  --cores <n>                   run n problems concurrently      [1]
        Parallelism is per problem: one worker runs all of a problem's
        solvers, modes and tolerances, so both sides of every matched
        cell see the same machine load and the speed-up ratios stay
        valid.  Absolute times are inflated by shared cache, memory
        bandwidth and turbo clocking, use --cores 1 for those, and
        never compare absolute times across different --cores values.
        The value is recorded in results.csv and run-info.txt.
  --compile-slots <n>           concurrent compilers            [min(cores,4)]
        The peak of a run is the build, not the solve: a large model with
        sensitivities is a very large translation unit and g++ answers it
        in gigabytes.  This bounds how many of those run at once, counted
        per *machine*, a second benchmark process on the same box shares
        the same slots.  0 disables the semaphore.
  --max-compile-gb <g>          address-space cap per compiler   [8]
        A model needing more fails to build and is skipped, instead of
        pushing the machine into swap.
  --max-worker-gb <g>           address-space cap per worker     [off]
        Needs the `unix` package; without it the run warns and continues.
  --max-states <n>              skip models with more states     [400]
  --max-sens <n>                cap on sensitivity parameters    [32]
  --min-points <n>              minimum output grid size         [25]
  --extra-solvers               also run cppDE BDF and Rosenbrock4
  --no-sparse-sweep             leave out the dense/sparse sweep
        By default every model where the linear solver is actually a
        choice is additionally run with it pinned dense and pinned
        sparse, on both backends, four extra cells per tolerance, next
        to the two auto-detected ones.  That is what validates the
        auto-detection instead of only using it, and it lands in the same
        results.csv (column `pinned`) and the same figures as the
        head-to-head.  Turn it off for a cheap run.
  --max-density <x>             sweep density cutoff             [0.25]
        Only models with a structurally sparser Jacobian are swept.  The
        codegen switches to KLU at 0.4, so raising this past 0.4 is what
        probes where the switch belongs rather than only what it costs.
  --min-sweep-states <n>        smallest system worth sweeping   [8]
        Below this the codegen never picks sparse, so both pinnings would
        measure a decision that is never taken.
  --include-excluded            also run models flagged unusable
  --outdir <dir>                results root      [benchmarks/results]
  --petab-root <dir>            PEtab collection  [benchmarks/cache/petab/...]
  --no-plots                    write the CSV only
  --quick                       tiny smoke run (a few small models)

Fetch the PEtab collection first:  Rscript benchmarks/fetch-models.R
")
  quit(save = "no")
}

as_bool <- function(x) isTRUE(as.logical(x))
QUICK <- as_bool(OPTS$quick)

TIER <- BENCH_TIERS[[OPTS$tier]]
if (is.null(TIER)) stop("unknown --tier '", OPTS$tier,
                        "'; expected one of ", paste(names(BENCH_TIERS), collapse = ", "))

## Shared with the two remote submitters (benchmarks/R/shards.R), so a
## local and a shipped run sweep the same tolerances.
TOLSETS <- BENCH_TOLSETS
## An explicit --tol always wins; otherwise the tier picks the set.
tol_name <- if ("tol" %in% names(user_opts)) OPTS$tol else TIER$tol
tolerances <- TOLSETS[[tol_name]]
if (is.null(tolerances)) stop("unknown --tol '", tol_name, "'")
if (QUICK) tolerances <- TOLSETS$loose

modes <- trimws(strsplit(OPTS$modes, ",")[[1L]])
nrep  <- if (QUICK) 3L else as.integer(OPTS$nrep)
max_states <- if (QUICK) 40L else
  if ("max-states" %in% names(user_opts)) as.integer(OPTS$`max-states`) else TIER$max_states
max_sens <- if ("max-sens" %in% names(user_opts)) as.integer(OPTS$`max-sens`) else
  TIER$max_sens
cores <- max(1L, as.integer(OPTS$cores))
if (cores > 1L && .Platform$OS.type != "unix") {
  warning("--cores > 1 needs a unix fork; falling back to serial",
          call. = FALSE, immediate. = TRUE)
  cores <- 1L
}
## --shard i/n: round-robin over a deterministic ordering, so the union
## of all shards is exactly the unsharded list, with no overlap.
shard_i <- NA_integer_; shard_n <- NA_integer_
if (nzchar(OPTS$shard)) {
  pp <- as.integer(strsplit(OPTS$shard, "/", fixed = TRUE)[[1L]])
  if (length(pp) != 2L || anyNA(pp) || pp[1L] < 1L || pp[2L] < 1L || pp[1L] > pp[2L])
    stop("--shard must look like i/n with 1 <= i <= n")
  shard_i <- pp[1L]; shard_n <- pp[2L]
}
shard_keep <- function(nms) {
  if (is.na(shard_n)) return(rep(TRUE, length(nms)))
  (rank(nms, ties.method = "first") - 1L) %% shard_n == (shard_i - 1L)
}
SPARSE_SWEEP <- as_bool(OPTS$`sparse-sweep`)
configs <- solver_configs(extra = as_bool(OPTS$`extra-solvers`))
sweep_configs <- if (SPARSE_SWEEP) sparse_sweep_configs() else NULL
model_filter <- trimws(strsplit(OPTS$models, ",")[[1L]])
model_filter <- model_filter[nzchar(model_filter)]
model_skip <- trimws(strsplit(OPTS$skip, ",")[[1L]])
model_skip <- model_skip[nzchar(model_skip)]

n_conditions <- if (identical(OPTS$conditions, "all")) "all" else
                  as.integer(OPTS$conditions)

builddir <- if (nzchar(OPTS$builddir)) OPTS$builddir else
  file.path(tempdir(), "cppDE_bench_build")
dir.create(builddir, showWarnings = FALSE, recursive = TRUE)

## Memory guards.  Applied here, before any model is built and before the
## workers are forked, so the caps are inherited rather than re-derived.
if (!nzchar(OPTS$`compile-slots`))
  OPTS$`compile-slots` <- as.character(min(cores, 4L))
compile_slots  <- as.integer(OPTS$`compile-slots`)
max_compile_gb <- suppressWarnings(as.numeric(OPTS$`max-compile-gb`))
max_worker_gb  <- if (nzchar(OPTS$`max-worker-gb`))
  suppressWarnings(as.numeric(OPTS$`max-worker-gb`)) else 0
if (is.na(compile_slots) || compile_slots < 0L)
  stop("--compile-slots must be a non-negative integer")
if (is.na(max_compile_gb) || is.na(max_worker_gb))
  stop("--max-compile-gb and --max-worker-gb must be numbers")
bench_limit_compiler(max_compile_gb, dir = builddir)
if (max_worker_gb > 0 && !isTRUE(bench_limit_process(max_worker_gb)))
  warning("no per-worker memory limit set: install the `unix` package",
          call. = FALSE, immediate. = TRUE)

## The directory name carries what the run was, not only when: a results tree is
## browsed far more often than opened, and "tiny_c8_nosens-sens1" answers at a
## glance. The timestamp stays first so the tree still sorts by time.
stamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
tag <- run_tag(stamp, OPTS$tier, cores, OPTS$modes, suite = OPTS$suite,
               models = model_filter,
               all_conditions = identical(OPTS$conditions, "all"))
if (!is.na(shard_n)) tag <- paste0(tag, sprintf("_shard%02dof%02d", shard_i, shard_n))
outdir <- file.path(OPTS$outdir, tag)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)


## ---------------------------------------------------------------------
##  Environment check
## ---------------------------------------------------------------------

cvode_ok <- isTRUE(cppDE:::cvodeConfig$available)
if (!cvode_ok)
  stop("The CVODE backend is unavailable, so there is nothing to compare against.\n",
       "  Install SUNDIALS and re-install cppDE, or run:\n",
       "      cppDE::install_libs(\"sundials\")")

cat(sprintf("cppDE benchmark suite\n  results  : %s\n  build    : %s\n", outdir, builddir))
cat(sprintf("  R %s | cppDE %s | KLU %s\n",
            getRversion(), utils::packageVersion("cppDE"),
            if (isTRUE(cppDE:::cvodeConfig$klu_available)) "yes" else "no"))


## ---------------------------------------------------------------------
##  Collect problems
## ---------------------------------------------------------------------

problems <- list(); skipped <- list()

if (OPTS$suite %in% c("all", "classic")) {
  cl <- classic_problems()
  if (!is.null(TIER$classic)) cl <- cl[intersect(TIER$classic, names(cl))]
  cl <- cl[shard_keep(names(cl))]
  if (QUICK) cl <- cl[intersect(c("robertson", "hires", "orego"), names(cl))]
  problems <- c(problems, lapply(cl, list))   # one "case" per problem
}

if (OPTS$suite %in% c("all", "petab")) {
  if (!dir.exists(OPTS$`petab-root`)) {
    warning("PEtab collection not found at ", OPTS$`petab-root`,
            "\n  run: Rscript benchmarks/fetch-models.R")
  } else {
    idx <- petab_list(OPTS$`petab-root`)
    if (!is.null(TIER$petab)) {
      missing <- setdiff(TIER$petab, idx$name)
      if (length(missing))
        warning("tier '", OPTS$tier, "' lists model(s) not in the collection: ",
                paste(missing, collapse = ", "), call. = FALSE, immediate. = TRUE)
      idx <- idx[idx$name %in% TIER$petab, , drop = FALSE]
    }
    if (length(model_filter))
      idx <- idx[Reduce(`|`, lapply(model_filter, function(f)
        grepl(f, idx$name, fixed = TRUE))), , drop = FALSE]
    if (length(model_skip))
      idx <- idx[!Reduce(`|`, lapply(model_skip, function(f)
        grepl(f, idx$name, fixed = TRUE))), , drop = FALSE]
    if (QUICK) idx <- idx[idx$name %in% c("Boehm_JProteomeRes2014",
                                          "Crauste_CellSystems2017",
                                          "Elowitz_Nature2000"), , drop = FALSE]
    ## Apply the size filter before parsing: the largest models in the
    ## collection take minutes to translate and would then be skipped.
    nsp <- vapply(idx$yaml, petab_species_count, 0L)
    too_big <- !is.na(nsp) & nsp > max_states
    for (i in which(too_big)) {
      skipped[[idx$name[i]]] <- sprintf("~%d species > --max-states %d",
                                        nsp[i], max_states)
      cat(sprintf("  [skip] %-34s %s\n", idx$name[i], skipped[[idx$name[i]]]))
    }
    idx <- idx[!too_big, , drop = FALSE]
    idx <- idx[shard_keep(idx$name), , drop = FALSE]
    cat(sprintf("  PEtab    : %d model(s), parsing with %d core(s)\n",
                nrow(idx), cores))

    ## Translating SBML happens before anything is timed, so it can be
    ## parallelised without affecting a single measurement.
    build_case <- function(i) tryCatch(
      petab_problems(idx$yaml[i],
                     conditions = if (identical(n_conditions, "all")) "all" else NULL,
                     max_conditions = if (identical(n_conditions, "all")) Inf
                                      else n_conditions,
                     max_sens = max_sens, min_points = as.integer(OPTS$`min-points`)),
      error = function(e) { message("  [skip] ", idx$name[i], ": ",
                                    sub("\n.*", "", conditionMessage(e))); list() })
    built <- if (cores > 1L && .Platform$OS.type == "unix" && nrow(idx) > 1L)
      parallel::mclapply(seq_len(nrow(idx)), build_case, mc.cores = cores)
    else lapply(seq_len(nrow(idx)), build_case)
    for (i in seq_along(built))
      if (length(built[[i]])) problems[[idx$name[i]]] <- built[[i]]
  }
}

## Match either the registry key ("e5") or the problem's display name
## ("E5"), case-insensitively, the two differ for the classic problems.
if (length(model_filter)) {
  disp <- vapply(problems, function(cases) cases[[1L]]$name, "")
  hit <- Reduce(`|`, lapply(model_filter, function(f)
    grepl(f, names(problems), ignore.case = TRUE) |
    grepl(f, disp, ignore.case = TRUE)))
  problems <- problems[hit]
}
if (length(model_skip)) {
  disp <- vapply(problems, function(cases) cases[[1L]]$name, "")
  out <- Reduce(`|`, lapply(model_skip, function(f)
    grepl(f, names(problems), ignore.case = TRUE) |
    grepl(f, disp, ignore.case = TRUE)))
  problems <- problems[!out]
}

if (!length(problems)) stop("no problems selected")

## The sweep narrows which models get the extra pinned cells; it never
## narrows the run itself, so trait coverage below stays the tier's.
if (SPARSE_SWEEP)
  problems <- mark_sparse_sweep(problems,
                                max_density = as.numeric(OPTS$`max-density`),
                                min_states  = as.integer(OPTS$`min-sweep-states`))

## A tier is only a usable regression basis if it still touches every
## trait; this warns loudly when a gap opens up.
coverage <- tier_coverage(problems, OPTS$tier)
print_tier_coverage(coverage, OPTS$tier)
cat(sprintf("\n%d model(s), %d case(s) selected\n",
            length(problems), sum(lengths(problems))))


## ---------------------------------------------------------------------
##  Run
## ---------------------------------------------------------------------

t_start <- proc.time()[["elapsed"]]

## Models flagged unusable are dropped before the run so that the reason
## reaches run-info.txt rather than a worker's stderr.
for (pname in names(problems)) {
  p0 <- problems[[pname]][[1L]]
  if (!isTRUE(p0$usable) && !as_bool(OPTS$`include-excluded`)) {
    skipped[[pname]] <- paste(p0$notes[grepl("EXCLUDED", p0$notes)], collapse = "; ")
    cat(sprintf("[skip] %-34s %s\n", pname, skipped[[pname]]))
  }
}
runnable <- problems[!names(problems) %in% names(skipped)]

df <- run_all_problems(
  runnable, builddir = builddir, tolerances = tolerances, modes = modes,
  configs = configs, sweep_configs = sweep_configs, nrep = nrep, cores = cores,
  max_sens2 = as.integer(OPTS$`max-sens2`), max_sens = max_sens,
  max_states = max_states, compile_slots = compile_slots,
  min_time = if (nzchar(OPTS$`min-time`)) as.numeric(OPTS$`min-time`) else NULL,
  on_skip = function(name, why) {
    skipped[[name]] <<- why
    cat(sprintf("[skip] %-34s %s\n", name, why))
  })

if (is.null(df) || !nrow(df)) stop("no successful runs")
elapsed <- proc.time()[["elapsed"]] - t_start


## ---------------------------------------------------------------------
##  Report
## ---------------------------------------------------------------------

csv <- file.path(outdir, "results.csv")
utils::write.csv(df, csv, row.names = FALSE)

write_run_readme(df, outdir, list(
  tag = tag, tier = OPTS$tier, modes = OPTS$modes, cores = cores,
  coverage = coverage, skipped = skipped,
  nrep = nrep, elapsed = sprintf("%.1f s", elapsed),
  date = format(Sys.time()),
  tolerances = paste(sprintf("atol %.0e / rtol %.0e",
                             tolerances$atol, tolerances$rtol), collapse = "; "),
  versions = sprintf("%s / %s", getRversion(), utils::packageVersion("cppDE")),
  cpu = tryCatch(sub(".*: ", "", grep("model name",
        readLines("/proc/cpuinfo", warn = FALSE), value = TRUE)[1L]),
        error = function(e) NULL),
  klu = isTRUE(cppDE:::cvodeConfig$klu_available),
  options = paste(sprintf("%s=%s", names(OPTS), unlist(OPTS)), collapse = " ")))

cat("\n\n==============================================================\n")
cat("SUMMARY, geometric mean speed-up versus CVODE_bdf\n")
cat("==============================================================\n")
s <- speedup_table(df)
if (!is.null(s)) {
  for (md in unique(s$mode)) {
    sub <- s[s$mode == md & s$solver != "CVODE_bdf", ]
    if (!nrow(sub)) next
    cat(sprintf("\n-- %s --\n", MODE_LABELS[md]))
    for (sv in unique(sub$solver)) {
      x <- sub[sub$solver == sv, ]
      cat(sprintf("  %-12s  time x%.2f   rhs-evals x%.2f   (%d cells, %d problems)\n",
                  sv, geo_mean(x$speedup), geo_mean(x$fev_ratio),
                  nrow(x), length(unique(x$problem))))
    }
  }
}
sg <- sparse_gain_table(df)
if (!is.null(sg) && nrow(sg)) {
  cat("\n==============================================================\n")
  cat("DENSE vs SPARSE LU, geometric mean, > 1 means sparse was faster\n")
  cat("==============================================================\n")
  for (bk in unique(sg$backend)) {
    x <- sg[sg$backend == bk, ]
    cat(sprintf("  %-8s x%.2f   (%d cells, %d problems)\n",
                bk, geo_mean(x$gain), nrow(x), length(unique(x$problem))))
  }
}

nfail <- sum(!df$ok)
cat(sprintf("\n%d rows, %d failed cells, %.1f s total\n", nrow(df), nfail, elapsed))
cat(sprintf("results: %s\n", csv))

if (as_bool(OPTS$plots)) {
  cat("\nwriting figures...\n")
  invisible(save_plots(df, outdir))
}
cat("done.\n")
