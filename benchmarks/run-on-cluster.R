#!/usr/bin/env Rscript
## =====================================================================
##  run-on-cluster.R -- submit the benchmark to a SLURM cluster (helix)
##  through dMod2's distributedComputing().
##
##      Rscript benchmarks/run-on-cluster.R --dry-run          # check locally
##      Rscript benchmarks/run-on-cluster.R --machine you@helix --submit
##      Rscript benchmarks/run-on-cluster.R --collect --jobname cppde_full
##
##  How the work is split
##  ---------------------
##  The problems are parsed *locally* and travel to the cluster inside
##  the transferred workspace, so the compute nodes need neither the
##  PEtab collection nor this repository -- only cppDE (with SUNDIALS
##  and KLU) installed for the R that `module load math/R` provides, and
##  a C++ compiler.  Each SLURM array task gets one shard of the problem
##  list and returns its rows; the shards are merged here.
##
##  What this does NOT fix
##  ----------------------
##  A benchmark measures wall-clock time, and array tasks land on
##  whatever nodes SLURM has free.  Consequently:
##
##    * ratios (cppDE vs CVODE) stay valid -- both solvers of a cell run
##      in the same task on the same node;
##    * absolute milliseconds are comparable *within* a shard only,
##      unless the partition is hardware-homogeneous and tasks got whole
##      nodes.  Every returned row therefore carries the node name and
##      CPU model so this can be checked afterwards rather than assumed.
##
##  Request whole nodes (`--cores` = the node's core count) if the
##  absolute numbers matter; otherwise a co-tenant job will show up in
##  them.
## =====================================================================

suppressPackageStartupMessages({
  library(cppDE)
  if (!requireNamespace("dMod2", quietly = TRUE))
    stop("dMod2 is required for cluster submission")
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

OPT <- list(
  jobname      = "cppde_bench",
  machine      = "helix",
  partition    = "cpu-single",
  cores        = "64",            # cores requested from SLURM (a whole node)
  `bench-cores` = "12",           # benchmark workers actually used
  `mem-per-core` = "2",           # GB per requested core
  walltime     = "08:00:00",
  shards       = "8",
  tier         = "full",
  modes        = "nosens,sens1,sens2",
  conditions   = "1",
  `max-states` = "400",
  `max-sens`   = "32",
  `max-sens2`  = "10",
  tol          = "default",
  nrep         = "5",
  `min-time`   = "0.25",
  `sparse-sweep` = "TRUE",
  `max-density` = "0.25",
  `min-sweep-states` = "8",
  ## Memory guards (R/harness.R).  SLURM already reserves the node's
  ## memory, so these are a second line rather than the only one: what
  ## they buy here is that a model too large for the allocation fails to
  ## build and is skipped, instead of the whole task being cancelled for
  ## exceeding it.
  `compile-slots` = "",           # "" = min(bench-cores, 4); 0 disables
  `max-compile-gb` = "8",
  `max-worker-gb` = "",           # needs the `unix` package on the nodes
  outdir       = file.path(ROOT, "benchmarks", "results"),
  `petab-root` = file.path(ROOT, "benchmarks", "cache", "petab", "Benchmark-Models"),
  `ssh-passwd` = "",
  `dry-run`    = "FALSE",
  submit       = "FALSE",
  collect      = "FALSE",
  help         = "FALSE")

args <- commandArgs(trailingOnly = TRUE)
SET <- character(0)
i <- 1L
while (i <= length(args)) {
  k <- sub("^--", "", args[i])
  if (identical(k, "no-sparse-sweep")) {
    OPT$`sparse-sweep` <- "FALSE"; SET <- c(SET, "sparse-sweep"); i <- i + 1L; next }
  if (k %in% c("dry-run", "submit", "collect", "sparse-sweep", "help")) {
    OPT[[k]] <- "TRUE"; SET <- c(SET, k); i <- i + 1L; next }
  if (grepl("=", k)) { kv <- strsplit(k, "=", fixed = TRUE)[[1L]]
    OPT[[kv[1L]]] <- paste(kv[-1L], collapse = "="); SET <- c(SET, kv[1L])
    i <- i + 1L; next }
  if (!k %in% names(OPT)) stop("unknown option --", k)
  OPT[[k]] <- args[i + 1L]; SET <- c(SET, k); i <- i + 2L
}
tf <- function(x) isTRUE(as.logical(x))

if (tf(OPT$help)) {
  cat("\nSubmit the cppDE benchmark to a SLURM cluster via dMod2.\n\n")
  cat(sprintf("  --%-14s [%s]\n", names(OPT), unlist(OPT)), sep = "")
  cat("\n  --dry-run   build the shards locally and report, submit nothing",
      "\n  --submit    build, transfer and submit",
      "\n  --collect   fetch results of an earlier --submit and merge them\n\n")
  quit(save = "no")
}

n_shards <- as.integer(OPT$shards)
cores    <- as.integer(OPT$cores)          # SLURM allocation
bench_cores <- as.integer(OPT$`bench-cores`)

## Requesting the whole node but running fewer workers is deliberate.
## Booking all 64 cores keeps co-tenants off the node, which is what
## makes a wall-clock measurement mean anything; running only a dozen
## workers keeps the job inside the node's memory, because generating
## and compiling a large model peaks in the gigabytes -- a 124-state
## model drove the local R process to 7.4 GB.  A full node at 2 GB per
## core holds ~128 GB, so a dozen workers is the right order and 64
## would not be.
if (bench_cores > cores) stop("--bench-cores cannot exceed --cores")
mem_per_worker <- cores * as.numeric(OPT$`mem-per-core`) / bench_cores
if (mem_per_worker < 6)
  warning(sprintf("only %.1f GB per benchmark worker; large models peak near 8 GB",
                  mem_per_worker), call. = FALSE, immediate. = TRUE)

if (!nzchar(OPT$`compile-slots`))
  OPT$`compile-slots` <- as.character(min(bench_cores, 4L))
compile_slots  <- as.integer(OPT$`compile-slots`)
max_compile_gb <- suppressWarnings(as.numeric(OPT$`max-compile-gb`))
max_worker_gb  <- if (nzchar(OPT$`max-worker-gb`))
  suppressWarnings(as.numeric(OPT$`max-worker-gb`)) else 0
if (is.na(compile_slots) || compile_slots < 0L)
  stop("--compile-slots must be a non-negative integer")
if (is.na(max_compile_gb) || is.na(max_worker_gb))
  stop("--max-compile-gb and --max-worker-gb must be numbers")

## One task owns the node here, so the ceiling the guards impose has to
## fit the allocation rather than the machine.
if (compile_slots > 0 && max_compile_gb > 0 &&
    compile_slots * max_compile_gb > cores * as.numeric(OPT$`mem-per-core`))
  warning(sprintf(paste("%d compile slot(s) at %g GB exceed the %g GB requested",
                        "from SLURM; the allocation, not the guard, will decide"),
                  compile_slots, max_compile_gb,
                  cores * as.numeric(OPT$`mem-per-core`)),
          call. = FALSE, immediate. = TRUE)


## ---------------------------------------------------------------------
##  Collect mode: fetch an earlier submission and merge
## ---------------------------------------------------------------------

JOBDIR <- file.path(ROOT, "benchmarks", "results", ".jobs")

if (tf(OPT$collect)) {
  ## The submission's own options are replayed from disk, so --collect
  ## needs nothing but the job name.  Anything given on the command line
  ## still wins, for the case where the record is missing or stale.
  rec <- file.path(JOBDIR, paste0(OPT$jobname, ".rds"))
  if (file.exists(rec)) {
    saved <- readRDS(rec)
    for (k in names(saved)) if (!k %in% SET) OPT[[k]] <- saved[[k]]
    n_shards <- as.integer(OPT$shards)
    bench_cores <- as.integer(OPT$`bench-cores`)
    cat(sprintf("replaying submission record: tier %s, modes %s, %d shard(s)\n",
                OPT$tier, OPT$modes, n_shards))
  } else {
    cat("no submission record found; using the options given on the command line.\n",
        "If --shards differs from the submission, the fetch will not line up.\n",
        sep = "")
  }
  ## dMod2 validates var_values/no_rep before it looks at `recover`, so
  ## the array shape has to be restated even when only fetching.
  job <- dMod2::distributedComputing({ NULL }, jobname = OPT$jobname,
                                    machine = OPT$machine, recover = TRUE,
                                    var_values = list(seq_len(n_shards)))
  if (!job$check()) message("note: the cluster reports the job as incomplete; ",
                            "merging whatever has arrived")
  job$get()
  res <- get("cluster_result", envir = .GlobalEnv)
  parts <- Filter(is.data.frame, if (is.data.frame(res)) list(res) else res)
  if (!length(parts)) stop("no data frames came back from the cluster")
  df <- do.call(rbind, parts)
  tag <- paste0(format(Sys.time(), "%Y%m%d-%H%M%S"), "_", OPT$tier,
                "_cluster_c", bench_cores, "_", gsub(",", "-", OPT$modes))
  out <- file.path(OPT$outdir, tag)
  dir.create(out, showWarnings = FALSE, recursive = TRUE)
  utils::write.csv(df, file.path(out, "results.csv"), row.names = FALSE)
  write_run_readme(df, out, list(
    tag = tag, tier = OPT$tier, modes = OPT$modes, cores = bench_cores,
    nrep = OPT$nrep, date = format(Sys.time()),
    options = paste(sprintf("%s=%s", names(OPT), unlist(OPT)), collapse = " ")))
  suppressWarnings(save_plots(df, out))
  cat(sprintf("\nmerged %d rows from %d shard(s) into %s\n",
              nrow(df), length(parts), out))
  if ("node" %in% names(df)) {
    cat("\nnodes used (absolute times are only comparable within one):\n")
    print(table(df$node))
  }
  quit(save = "no")
}


## ---------------------------------------------------------------------
##  Build the problems locally, then shard them
## ---------------------------------------------------------------------

TIER <- BENCH_TIERS[[OPT$tier]]
if (is.null(TIER)) stop("unknown --tier '", OPT$tier, "'")

## Same rule as run-benchmarks.R: an explicit --tol wins, otherwise the
## tier decides.  Without this the cluster run would sweep a different
## tolerance set than the identical local command, and the two would not
## be comparable.
tol_name <- if ("tol" %in% SET) OPT$tol else TIER$tol
tolerances <- BENCH_TOLSETS[[tol_name]]
if (is.null(tolerances)) stop("unknown --tol '", tol_name, "'")
max_states <- as.integer(OPT$`max-states`)
max_sens   <- as.integer(OPT$`max-sens`)

cat(sprintf("building problems for tier '%s' ...\n", OPT$tier))
problems <- bench_problems_for_tier(OPT$tier, OPT$`petab-root`,
                                    conditions = OPT$conditions,
                                    max_states = max_states, max_sens = max_sens)
if (!length(problems)) stop("no problems selected")

## Marks which models additionally get the pinned dense/sparse cells; the
## problem list itself is unchanged, so the shard balance below already
## sees the extra weight.
if (isTRUE(as.logical(OPT$`sparse-sweep`)))
  problems <- mark_sparse_sweep(problems,
                                max_density = as.numeric(OPT$`max-density`),
                                min_states  = as.integer(OPT$`min-sweep-states`))

sh <- balance_shards(problems, n_shards)
shards <- sh$shards

cov <- tier_coverage(problems, OPT$tier)
print_tier_coverage(cov, OPT$tier)
cat(sprintf("\n%d problems, %d cases -> %d shard(s)\n",
            length(problems), sum(lengths(problems)), n_shards))
print_shard_plan(shards, sh$load)

sz <- format(utils::object.size(shards), units = "MB")
cat(sprintf("\nworkspace payload: %s\n", sz))
cat(sprintf("tolerances (%s): %s\n",
            if ("tol" %in% SET) "explicit --tol" else paste("from tier", OPT$tier),
            paste(sprintf("atol %.0e / rtol %.0e", tolerances$atol, tolerances$rtol),
                  collapse = "; ")))
cat(sprintf("equivalent local run:\n  Rscript benchmarks/run-benchmarks.R --tier %s --modes %s --cores %d\n",
            OPT$tier, OPT$modes, bench_cores))

if (!tf(OPT$submit)) {
  cat("\nDry run -- nothing submitted.\n")
  cat("Before submitting, verify on the cluster:\n",
      " 1. cppDE installs and loads for `module load math/R`\n",
      " 2. cppDE:::cvodeConfig$available is TRUE there (CVODE backend)\n",
      " 3. a C++ compiler is on PATH on the compute nodes -- the benchmark\n",
      "    generates and compiles C++ at run time\n",
      " 4. $TMPDIR is node-local and writable; that is where models are built\n\n",
      "Then re-run with --machine user@host --submit\n", sep = "")
  quit(save = "no")
}


## ---------------------------------------------------------------------
##  Submit
## ---------------------------------------------------------------------

if (!nzchar(OPT$machine)) stop("--submit needs --machine user@host")

modes    <- trimws(strsplit(OPT$modes, ",")[[1L]])
cfgs     <- solver_configs()
sweep_cfgs <- if (isTRUE(as.logical(OPT$`sparse-sweep`))) sparse_sweep_configs() else NULL
nrep     <- as.integer(OPT$nrep)
min_time <- as.numeric(OPT$`min-time`)
maxs2    <- as.integer(OPT$`max-sens2`)

cat(sprintf("\nsubmitting '%s' to %s\n  partition %s, %d cores allocated, %d workers, walltime %s\n",
            OPT$jobname, OPT$machine, OPT$partition, cores, bench_cores, OPT$walltime))

job <- dMod2::distributedComputing(
  {
    ## Runs on a compute node.  `var_1` is this task's shard index; the
    ## benchmark functions and `shards` arrive in the workspace.
    ##
    ## The workspace carries function *objects*, not attached packages,
    ## so cppDE has to be loaded here -- without this every model fails
    ## with "could not find function cppDE" and the shard returns
    ## nothing.
    library(cppDE)
    stopifnot(isTRUE(cppDE:::cvodeConfig$available))
    Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
    builddir <- file.path(tempdir(), paste0("cppde_bench_", var_1))
    dir.create(builddir, showWarnings = FALSE, recursive = TRUE)

    ## Memory guards, before the fork so the workers inherit them.  The
    ## Makevars goes under builddir, which is $TMPDIR and node-local --
    ## on a cluster with a shared home that is the difference between one
    ## file per task and every task fighting over the same one.
    bench_limit_compiler(max_compile_gb, dir = builddir)
    if (max_worker_gb > 0 && !isTRUE(bench_limit_process(max_worker_gb)))
      message("note: could not set a per-worker memory limit on ",
              Sys.info()[["nodename"]])

    out <- run_all_problems(shards[[var_1]], builddir = builddir,
                            tolerances = tolerances, modes = modes,
                            configs = cfgs, sweep_configs = sweep_cfgs,
                            nrep = nrep, cores = bench_cores,
                            max_sens2 = maxs2, min_time = min_time,
                            compile_slots = compile_slots)
    ## Stamp the hardware: absolute times are only comparable across rows
    ## that ran on the same node.
    if (!is.null(out) && nrow(out)) {
      out$shard <- var_1
      out$node  <- unname(Sys.info()[["nodename"]])
      out$cpu   <- tryCatch(sub(".*: ", "", grep("model name",
                     readLines("/proc/cpuinfo", warn = FALSE), value = TRUE)[1L]),
                     error = function(e) NA_character_)
    }
    out
  },
  jobname      = OPT$jobname,
  partition    = OPT$partition,
  cores        = cores,
  nodes        = 1,
  mem_per_core = as.numeric(OPT$`mem-per-core`),
  walltime     = OPT$walltime,
  machine      = OPT$machine,
  ssh_passwd   = if (nzchar(OPT$`ssh-passwd`)) OPT$`ssh-passwd` else NULL,
  var_values   = list(seq_len(n_shards)),
  recover      = FALSE)

dir.create(JOBDIR, showWarnings = FALSE, recursive = TRUE)
saveRDS(OPT[setdiff(names(OPT), c("dry-run", "submit", "collect", "help"))],
        file.path(JOBDIR, paste0(OPT$jobname, ".rds")))

cat("\nsubmitted. Poll and fetch with:\n",
    sprintf("  Rscript benchmarks/run-on-cluster.R --collect --jobname %s\n",
            OPT$jobname), sep = "")
