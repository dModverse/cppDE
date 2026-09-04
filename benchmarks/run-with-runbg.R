#!/usr/bin/env Rscript

## =====================================================================
##  run-with-runbg.R: run the benchmark on plain ssh hosts via runbg().
## =====================================================================

##      Rscript benchmarks/run-with-runbg.R --dry-run
##      Rscript benchmarks/run-with-runbg.R --machines you@box --submit
##      Rscript benchmarks/run-with-runbg.R --collect --jobname cppde_bg

##  Same tiers and cost-balanced shards as run-on-cluster.R, so results from
##  the two routes are comparable. The problems are parsed here and travel in
##  the workspace; a host needs only cppDE and a C++ compiler.

##  runbg() has no scheduler: `--machines` is the entire placement policy, and
##  naming a host n times starts n concurrent R processes on it. The preflight
##  refuses a placement that does not fit and the memory guards bound the rest.

##  Nothing reserves the machine, so ratios stay valid while absolute
##  milliseconds do not; every row carries node and cpu. Authentication has to
##  be non-interactive, including for `--machines localhost`.

##  runbg's own `walltime` is not exposed: in dMod2 it splices a second
##  statement into the try() call and the remote script no longer parses.
##  Bound the run with the tier and --max-states instead.

suppressPackageStartupMessages({
  library(cppDE)
  if (!requireNamespace("dMod2", quietly = TRUE))
    stop("dMod2 is required for runbg submission")
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
  jobname      = "cppde_bg",
  machines     = "localhost",     # comma-separated host pool for ssh
  shards       = "",              # jobs to start; default: one per host
  `bench-cores` = "4",            # benchmark workers inside one job
  tier         = "tiny",
  modes        = "nosens,sens1",
  conditions   = "1",
  `max-states` = "400",
  `max-sens`   = "32",
  `max-sens2`  = "10",
  tol          = "default",
  nrep         = "5",
  `min-time`   = "0.25",
  `log-lines`  = "40",
  `sparse-sweep` = "TRUE",
  `max-density` = "0.25",
  `min-sweep-states` = "8",
  ## Memory guards; see the "Memory guards" section of R/harness.R for
  ## what each one bounds and why one alone is not enough.
  `compile-slots` = "",           # "" = min(bench-cores, 4); 0 disables
  `max-compile-gb` = "8",         # ulimit -v per compiler process
  `max-worker-gb` = "",           # RLIMIT_AS per worker; "" = off
  `mem-per-worker` = "6",         # preflight floor, GB available per worker
  `mem-budget`  = "0.5",          # fraction of a host's free RAM we may plan for
  outdir       = file.path(ROOT, "benchmarks", "results"),
  `petab-root` = file.path(ROOT, "benchmarks", "cache", "petab", "Benchmark-Models"),
  `dry-run`    = "FALSE",
  submit       = "FALSE",
  wait         = "FALSE",
  check        = "FALSE",
  collect      = "FALSE",
  log          = "FALSE",
  purge        = "FALSE",
  terminate    = "FALSE",
  force        = "FALSE",
  help         = "FALSE")

FLAGS <- c("dry-run", "submit", "wait", "check", "collect", "log",
           "purge", "terminate", "sparse-sweep", "force", "help")

args <- commandArgs(trailingOnly = TRUE)
SET <- character(0)
i <- 1L
while (i <= length(args)) {
  k <- sub("^--", "", args[i])
  if (identical(k, "no-sparse-sweep")) {
    OPT$`sparse-sweep` <- "FALSE"; SET <- c(SET, "sparse-sweep"); i <- i + 1L; next }
  if (k %in% FLAGS) { OPT[[k]] <- "TRUE"; SET <- c(SET, k); i <- i + 1L; next }
  if (grepl("=", k)) { kv <- strsplit(k, "=", fixed = TRUE)[[1L]]
    OPT[[kv[1L]]] <- paste(kv[-1L], collapse = "="); SET <- c(SET, kv[1L])
    i <- i + 1L; next }
  if (!k %in% names(OPT)) stop("unknown option --", k)
  OPT[[k]] <- args[i + 1L]; SET <- c(SET, k); i <- i + 2L
}
tf <- function(x) isTRUE(as.logical(x))

if (tf(OPT$help)) {
  cat("\nRun the cppDE benchmark on ssh-reachable hosts via dMod2::runbg().\n\n")
  cat(sprintf("  --%-14s [%s]\n", names(OPT), unlist(OPT)), sep = "")
  cat("\n  --dry-run    build and balance the shards locally, send nothing",
      "\n  --submit     probe the hosts, transfer and start the jobs",
      "\n  --wait       with --submit: block until done, then merge",
      "\n               (runbg runs the hosts in sequence when it waits,",
      "\n                so this is meant for a single shard)",
      "\n  --check      ask whether the results are on the hosts yet",
      "\n  --log        tail each host's R CMD BATCH output",
      "\n  --collect    fetch an earlier --submit and merge it",
      "\n  --purge      delete the remote job folders and local scratch",
      "\n  --terminate  kill the remote R processes of this job",
      "\n  --no-sparse-sweep  leave out the dense/sparse sweep.  On by",
      "\n               default: every model where the linear solver is a",
      "\n               real choice (>= --min-sweep-states states, denser",
      "\n               than --max-density) is additionally run with it",
      "\n               pinned both ways on both backends, in the same",
      "\n               shard and the same results.csv as the head-to-head\n")
  cat("\n  Memory guards, a host without a scheduler has nothing else",
      "\n  standing between a sweep and the kernel's OOM path:",
      "\n  --compile-slots <n>   concurrent compilers per *host*, shared by",
      "\n                        every job on it; the build is the peak, so",
      "\n                        this is what bounds a host's footprint.",
      "\n                        Default min(--bench-cores, 4); 0 disables.",
      "\n  --max-compile-gb <g>  ulimit -v per compiler; a model that needs",
      "\n                        more fails to build instead of swapping",
      "\n  --max-worker-gb <g>   RLIMIT_AS per worker (needs the `unix`",
      "\n                        package on the hosts); off by default",
      "\n  --mem-per-worker <g>  preflight refuses to submit below this much",
      "\n                        available RAM per worker",
      "\n  --mem-budget <f>      fraction of a host's free RAM the build",
      "\n                        ceiling may claim; above it, refuse  [0.5]",
      "\n  --force               downgrade the preflight's refusals to",
      "\n                        warnings and submit anyway\n\n")
  quit(save = "no")
}

## Written by --submit, replayed by every later mode, so --collect needs
## nothing but the job name.  Anything given on the command line still
## wins, for a record that is missing or stale.
JOBDIR <- file.path(ROOT, "benchmarks", "results", ".jobs")
REC    <- file.path(JOBDIR, paste0(OPT$jobname, ".runbg.rds"))

if (any(c("check", "collect", "log", "purge", "terminate") %in% SET) &&
    file.exists(REC)) {
  saved <- readRDS(REC)
  for (k in names(saved)) if (!k %in% SET) OPT[[k]] <- saved[[k]]
  cat(sprintf("replaying submission record: tier %s, modes %s, %s shard(s) on %s\n",
              OPT$tier, OPT$modes, OPT$shards, OPT$machines))
}

## Paths are resolved before anything changes directory: runbg writes its
## scratch into, and fetches results into, the working directory.
OPT$outdir       <- normalizePath(OPT$outdir, mustWork = FALSE)
OPT$`petab-root` <- normalizePath(OPT$`petab-root`, mustWork = FALSE)

machines <- trimws(strsplit(OPT$machines, ",")[[1L]])
machines <- machines[nzchar(machines)]
if (!length(machines)) stop("--machines needs at least one host")
n_shards <- if (nzchar(OPT$shards)) as.integer(OPT$shards) else length(machines)
if (is.na(n_shards) || n_shards < 1L) stop("--shards must be a positive integer")
## One runbg job per shard; the host pool is dealt round-robin, so a
## single host named once and --shards 4 gives four jobs on that host.
placement <- machines[((seq_len(n_shards) - 1L) %% length(machines)) + 1L]
bench_cores <- as.integer(OPT$`bench-cores`)
OPT$shards <- as.character(n_shards)

## The semaphore counts per host, not per job, so the default is a property of
## the host. Four concurrent builds at the 8 GB ceiling is ~32 GB of compiler
## on a machine, which is the budget being defended.
if (!nzchar(OPT$`compile-slots`))
  OPT$`compile-slots` <- as.character(min(bench_cores, 4L))
compile_slots  <- as.integer(OPT$`compile-slots`)
max_compile_gb <- suppressWarnings(as.numeric(OPT$`max-compile-gb`))
max_worker_gb  <- if (nzchar(OPT$`max-worker-gb`))
  suppressWarnings(as.numeric(OPT$`max-worker-gb`)) else 0
if (is.na(compile_slots) || compile_slots < 0L)
  stop("--compile-slots must be a non-negative integer")
if (is.na(max_compile_gb)) stop("--max-compile-gb must be a number")
if (is.na(max_worker_gb)) stop("--max-worker-gb must be a number")

## runbg() names its files after `filename` and works in getwd(), scp'ing every
## *.c/*.cpp/*.o/*.so next to them and fetching results back there. A scratch
## directory keeps it from shipping unrelated sources into the repo.
STAGE <- file.path(ROOT, "benchmarks", "results", ".runbg")
dir.create(STAGE, showWarnings = FALSE, recursive = TRUE)


## ---------------------------------------------------------------------
##  Talking to the hosts
## ---------------------------------------------------------------------

## BatchMode=yes turns a missing key into an error instead of a password prompt
## that would hang an Rscript forever. A host that answered but could not report
## a field prints "?", so a gap reads as a gap and not as a number.
fmt_num <- function(x, digits = 0) {
  if (!isTRUE(is.finite(x))) "?" else sprintf(paste0("%.", digits, "f"), x)
}

ssh_run <- function(host, cmd) {
  out <- suppressWarnings(system2(
    "ssh", c("-o", "BatchMode=yes", "-o", "ConnectTimeout=10",
             shQuote(host), shQuote(cmd)),
    stdout = TRUE, stderr = TRUE))
  list(status = attr(out, "status") %||% 0L, out = out)
}

## Everything the job needs from a host, in one round trip and with the same R
## runbg will start (`R`, not `Rscript`, and --vanilla, so ~/.Renviron is out of
## the picture on both sides alike).

## The last field generates, compiles and solves a one-state model, because
## `cvodeConfig$available` is an install-time flag: a SUNDIALS displaced on the
## loader path links and loads, then aborts the process at the first solve.

## The static fields are flushed before the solve, so a reply that stops short
## is distinguishable from a host that never answered. The reply is tagged,
## because compiling prints the compiler command line into the same stream.

## The capacity fields are read on the host rather than guessed: the pool is
## heterogeneous and `MemAvailable` moves with whoever else is logged in.
host_report <- function(host) {
  rcode <- paste0(
    'mi <- tryCatch(readLines("/proc/meminfo", warn = FALSE),',
    ' error = function(e) character(0));',
    'gb <- function(k) { l <- grep(paste0("^", k, ":"), mi, value = TRUE);',
    ' if (length(l)) as.numeric(gsub("[^0-9]", "", l[[1]])) / 1048576 else NA };',
    'la <- tryCatch(as.numeric(strsplit(readLines("/proc/loadavg",',
    ' warn = FALSE)[[1]], " ")[[1]][[1]]), error = function(e) NA);',
    'v <- c(as.character(getRversion()),',
    ' isTRUE(try(requireNamespace("cppDE", quietly = TRUE), silent = TRUE)),',
    ' isTRUE(try(cppDE:::cvodeConfig$available, silent = TRUE)),',
    ' isTRUE(try(cppDE:::cvodeConfig$klu_available, silent = TRUE)),',
    ' nzchar(Sys.which("g++")) || nzchar(Sys.which("clang++")),',
    ' parallel::detectCores(), sprintf("%.1f", gb("MemTotal")),',
    ' sprintf("%.1f", gb("MemAvailable")), sprintf("%.2f", la),',
    ' requireNamespace("unix", quietly = TRUE));',
    'cat("CPPDEPROBEa,", paste(v, collapse = ","), "\\n", sep = ""); flush(stdout());',
    'ok <- tryCatch({',
    ' m <- cppDE::cvode(c(A = "-k*A"), modelname = "cppde_hostprobe", compile = TRUE);',
    ' r <- cppDE::solveODE(m, 0:1, c(A = 1, k = 1));',
    ' isTRUE(all(is.finite(r$variable)))',
    '}, error = function(e) FALSE);',
    'cat("CPPDEPROBEb,", ok, "\\n", sep = "")')
  r <- ssh_run(host, paste("R --vanilla --no-echo -e", shQuote(rcode)))
  tail3 <- function() paste(utils::tail(r$out, 3L), collapse = " | ")
  pick <- function(tag) {
    ln <- grep(tag, r$out, fixed = TRUE, value = TRUE)
    if (!length(ln)) return(NULL)
    strsplit(sub(paste0(".*", tag, ","), "", utils::tail(ln, 1L)), ",", fixed = TRUE)[[1L]]
  }
  f <- pick("CPPDEPROBEa")
  if (is.null(f) || length(f) < 10L) return(list(ok = FALSE, why = tail3()))
  s <- pick("CPPDEPROBEb")
  num <- function(x) suppressWarnings(as.numeric(x))
  list(ok = TRUE, rversion = f[1L], cppde = f[2L] == "TRUE",
       cvode = f[3L] == "TRUE", klu = f[4L] == "TRUE", cxx = f[5L] == "TRUE",
       cores = num(f[6L]), mem_total = num(f[7L]), mem_avail = num(f[8L]),
       load1 = num(f[9L]), unixpkg = f[10L] == "TRUE",
       ## The static fields are flushed before the solve, so a reply that stops after
       ## them means the solve took the R process down rather than raising a
       ## condition, which is what a displaced SUNDIALS does and no flag can report.
       solve = if (is.null(s)) NA else s[1L] == "TRUE",
       why = if (is.null(s)) paste("the CVODE solve aborted the R process:",
                                   tail3()) else "")
}

## Probing every distinct host beats discovering on collect that the jobs died
## an hour ago. A missing CVODE is a fact and always refuses; a placement too
## big for the box is an estimate and can be overridden with --force.
preflight <- function(hosts) {
  bad <- character(0)
  over <- character(0)
  floor_gb <- suppressWarnings(as.numeric(OPT$`mem-per-worker`))
  budget_frac <- suppressWarnings(as.numeric(OPT$`mem-budget`))
  if (!isTRUE(is.finite(budget_frac)) || budget_frac <= 0 || budget_frac > 1)
    stop("--mem-budget must be a fraction in (0, 1]")
  for (h in unique(hosts)) {
    rep <- host_report(h)
    if (!rep$ok) {
      cat(sprintf("  %-24s UNREACHABLE / no R: %s\n", h, rep$why))
      bad <- c(bad, h); next
    }
    cat(sprintf("  %-24s R %s | cppDE %s | CVODE %s | KLU %s | C++ %s | solve %s\n",
                h, rep$rversion,
                if (rep$cppde) "yes" else "NO",
                if (rep$cvode) "yes" else "NO",
                if (rep$klu) "yes" else "no",
                if (rep$cxx) "yes" else "NO",
                if (isTRUE(rep$solve)) "yes" else if (is.na(rep$solve)) "CRASH" else "NO"))
    if (nzchar(rep$why %||% "")) cat(sprintf("  %-24s %s\n", "", rep$why))
    if (!rep$cppde || !rep$cvode || !rep$cxx || !isTRUE(rep$solve)) bad <- c(bad, h)

    ## What this submission is about to ask of this particular box.
    jobs    <- sum(hosts == h)
    workers <- jobs * bench_cores
    cat(sprintf("  %-24s %s cores, %s GB RAM (%s GB available), load %s -> %d job(s) x %d = %d worker(s)\n",
                "", fmt_num(rep$cores, 0), fmt_num(rep$mem_total, 0),
                fmt_num(rep$mem_avail, 0), fmt_num(rep$load1, 2),
                jobs, bench_cores, workers))
    if (isTRUE(rep$cores > 0) && workers > rep$cores) {
      cat(sprintf("  %-24s TOO MANY WORKERS: %d workers on %s cores\n",
                  "", workers, fmt_num(rep$cores, 0)))
      over <- c(over, h)
    }
    if (isTRUE(is.finite(rep$mem_avail)) && isTRUE(is.finite(floor_gb)) &&
        rep$mem_avail / workers < floor_gb) {
      cat(sprintf("  %-24s TOO LITTLE MEMORY: %.1f GB per worker, floor is %.1f\n",
                  "", rep$mem_avail / workers, floor_gb))
      over <- c(over, h)
    }

    ## The guards' own ceiling, checked against the machine rather than assumed to
    ## fit it. The semaphore counts per host, so it bounds the concurrent compilers,
    ## but it can never exceed the workers either.

    ## Planning for only part of the free memory is deliberate: MemAvailable is a
    ## snapshot of a box nothing reserves, and the next co-tenant has to fit too.
    eff_slots <- if (compile_slots > 0L) min(compile_slots, workers) else workers
    ceiling_gb <- if (max_compile_gb > 0) eff_slots * max_compile_gb else Inf
    cat(sprintf("  %-24s build ceiling %s vs %s GB budget (%.0f%% of free)\n", "",
                if (is.finite(ceiling_gb))
                  sprintf("%s GB (%d compiler(s) x %s GB)", format(ceiling_gb),
                          eff_slots, format(max_compile_gb))
                else "UNBOUNDED",
                fmt_num(budget_frac * rep$mem_avail, 0), 100 * budget_frac))
    if (!is.finite(ceiling_gb)) {
      cat(sprintf("  %-24s NO COMPILE CEILING: --max-compile-gb is off, so nothing\n", ""))
      cat(sprintf("  %-24s   bounds what the builds may allocate on this host\n", ""))
      over <- c(over, h)
    } else if (isTRUE(is.finite(rep$mem_avail)) &&
               ceiling_gb > budget_frac * rep$mem_avail) {
      cat(sprintf("  %-24s CEILING EXCEEDS BUDGET: %s GB planned, %.0f%% of %s GB free is %s GB\n",
                  "", format(ceiling_gb), 100 * budget_frac,
                  fmt_num(rep$mem_avail, 0), fmt_num(budget_frac * rep$mem_avail, 0)))
      over <- c(over, h)
    }
    if (max_worker_gb > 0 && !isTRUE(rep$unixpkg)) {
      cat(sprintf("  %-24s --max-worker-gb asked for, but the `unix` package is missing here;\n",
                  ""))
      cat(sprintf("  %-24s   the per-worker cap will NOT be enforced (install.packages(\"unix\"))\n",
                  ""))
    }
  }

  ## A broken host is reported before an overloaded one: it is the harder
  ## fact, and --force must never be able to talk past it.
  if (length(bad))
    stop("host(s) not ready: ", paste(unique(bad), collapse = ", "),
         "\n  runbg starts `R CMD BATCH --vanilla` over a non-login ssh shell,",
         "\n  so anything a module system or ~/.profile sets up must be reachable",
         "\n  from ~/.bashrc. The benchmark also compiles C++ at run time and",
         "\n  needs cppDE's CVODE backend (cppDE::install_libs(\"sundials\")).",
         "\n  A failing or crashing solve with everything else green usually",
         "\n  means another SUNDIALS on the loader's search path is displacing",
         "\n  the one cppDE was installed against; run",
         "\n  `ssh <host> 'R --vanilla --no-echo' < dev/diagnose-cvode.R` to see",
         "\n  which libraries the model actually loads.")

  if (length(over)) {
    msg <- paste0(
      "placement exceeds the host(s): ", paste(unique(over), collapse = ", "),
      "\n  Lower --shards or --bench-cores, spread over more --machines,",
      "\n  or wait for the box to empty out. A sweep that oversubscribes a",
      "\n  machine without a scheduler does not merely run slowly: enough",
      "\n  concurrent compilers drive the kernel into reclaim until it stops",
      "\n  answering, and the whole job is lost along with the host.",
      "\n  Override with --force once you are sure the numbers are wrong.")
    if (tf(OPT$force)) warning(msg, call. = FALSE, immediate. = TRUE)
    else stop(msg)
  }
}

## The control handles of an earlier submission, without re-sending it.
runbg_handle <- function() dMod2::runbg({ NULL }, machine = placement,
                                       filename = OPT$jobname, recover = TRUE)


## ---------------------------------------------------------------------
##  Merging what came back
## ---------------------------------------------------------------------

merge_results <- function(res) {
  res <- if (is.data.frame(res)) list(res) else res
  for (k in seq_along(res))
    if (inherits(res[[k]], "try-error"))
      cat(sprintf("  shard %d failed on %s: %s\n", k, placement[k],
                  sub("\n.*", "", conditionMessage(attr(res[[k]], "condition")))))
  parts <- Filter(function(x) is.data.frame(x) && nrow(x), res)
  if (!length(parts))
    stop("no data frames came back; try --log to see what the hosts did")
  df <- do.call(rbind, parts)
  tag <- paste0(format(Sys.time(), "%Y%m%d-%H%M%S"), "_", OPT$tier,
                "_runbg_c", bench_cores, "_", gsub(",", "-", OPT$modes))
  out <- file.path(OPT$outdir, tag)
  dir.create(out, showWarnings = FALSE, recursive = TRUE)
  utils::write.csv(df, file.path(out, "results.csv"), row.names = FALSE)
  write_run_readme(df, out, list(
    tag = tag, tier = OPT$tier, modes = OPT$modes, cores = bench_cores,
    nrep = OPT$nrep, date = format(Sys.time()),
    options = paste(sprintf("%s=%s", names(OPT), unlist(OPT)), collapse = " ")))
  suppressWarnings(save_plots(df, out))
  cat(sprintf("\nmerged %d rows from %d of %d shard(s) into %s\n",
              nrow(df), length(parts), length(res), out))
  if ("node" %in% names(df)) {
    cat("\nnodes used (absolute times are only comparable within one):\n")
    print(table(df$node))
  }
  invisible(df)
}


## ---------------------------------------------------------------------
##  Modes that only need the job name
## ---------------------------------------------------------------------

## R CMD BATCH names its transcript after the basename of the script and ssh
## drops us in the home directory, so the log sits beside the job folder rather
## than inside it. Both are tried: an unfindable log is an undebuggable job.
if (tf(OPT$log)) {
  n <- as.integer(OPT$`log-lines`)
  for (k in seq_along(placement)) {
    f <- sprintf("%s_%d", OPT$jobname, k)
    cat(sprintf("\n===== shard %d on %s =====\n", k, placement[k]))
    r <- ssh_run(placement[k], sprintf(
      "tail -n %d %s.Rout 2>/dev/null || tail -n %d %s_folder/%s.Rout", n, f, n, f, f))
    cat(paste(r$out, collapse = "\n"), "\n")
  }
  quit(save = "no")
}

if (tf(OPT$terminate)) { runbg_handle()$terminate(); quit(save = "no") }

if (tf(OPT$purge)) {
  setwd(STAGE)
  runbg_handle()$purge()
  unlink(REC)
  cat("purged remote folders and local scratch for '", OPT$jobname, "'\n", sep = "")
  quit(save = "no")
}

if (tf(OPT$check)) {
  setwd(STAGE)
  invisible(runbg_handle()$check())
  quit(save = "no")
}

if (tf(OPT$collect)) {
  setwd(STAGE)
  job <- runbg_handle()
  if (!job$check())
    message("note: not every shard has reported; merging whatever has arrived")
  job$get()
  merge_results(get(".runbgOutput", envir = .GlobalEnv))
  quit(save = "no")
}


## ---------------------------------------------------------------------
##  Build the problems locally, then shard them
## ---------------------------------------------------------------------

TIER <- BENCH_TIERS[[OPT$tier]]
if (is.null(TIER)) stop("unknown --tier '", OPT$tier, "'")

## An explicit --tol wins, otherwise the tier decides, the same rule as
## run-benchmarks.R and run-on-cluster.R, so the identical command run
## the other way sweeps the identical tolerances.
tol_name <- if ("tol" %in% SET) OPT$tol else TIER$tol
tolerances <- BENCH_TOLSETS[[tol_name]]
if (is.null(tolerances)) stop("unknown --tol '", tol_name, "'")

cat(sprintf("building problems for tier '%s' ...\n", OPT$tier))
problems <- bench_problems_for_tier(OPT$tier, OPT$`petab-root`,
                                    conditions = OPT$conditions,
                                    max_states = as.integer(OPT$`max-states`),
                                    max_sens = as.integer(OPT$`max-sens`))
if (!length(problems)) stop("no problems selected")

## The sweep only adds cells to the models it marks; the problem list,
## and with it the shard balance and the trait coverage, is the tier's
## either way.
if (tf(OPT$`sparse-sweep`))
  problems <- mark_sparse_sweep(problems,
                                max_density = as.numeric(OPT$`max-density`),
                                min_states  = as.integer(OPT$`min-sweep-states`))

sh <- balance_shards(problems, n_shards)
shards <- sh$shards

cov <- tier_coverage(problems, OPT$tier)
print_tier_coverage(cov, OPT$tier)
cat(sprintf("\n%d problems, %d cases -> %d shard(s) over %d host(s)\n",
            length(problems), sum(lengths(problems)), n_shards,
            length(unique(placement))))
print_shard_plan(shards, sh$load)
for (h in unique(placement))
  cat(sprintf("  %-24s %d job(s) x %d worker(s) = %d concurrent solves\n",
              h, sum(placement == h), bench_cores,
              sum(placement == h) * bench_cores))

## The compile ceiling is per host and not per job, so it is reported
## that way, with two jobs on a box the number below is still the whole
## budget, which is the property that makes it safe.
cat(sprintf("memory guards: %s concurrent compiler(s) per host x %s GB = %s GB peak build%s\n",
            if (compile_slots > 0L) as.character(compile_slots) else "unlimited",
            if (max_compile_gb > 0) format(max_compile_gb) else "unlimited",
            if (compile_slots > 0L && max_compile_gb > 0)
              format(compile_slots * max_compile_gb) else "?",
            if (max_worker_gb > 0)
              sprintf("; %s GB per worker", format(max_worker_gb)) else ""))

modes    <- trimws(strsplit(OPT$modes, ",")[[1L]])
cfgs     <- solver_configs()
sweep_cfgs <- if (tf(OPT$`sparse-sweep`)) sparse_sweep_configs() else NULL
nrep     <- as.integer(OPT$nrep)
min_time <- as.numeric(OPT$`min-time`)
maxs2    <- as.integer(OPT$`max-sens2`)

## runbg() saves one workspace and copies it to every job, so the payload is
## paid per shard. Naming what travels keeps `problems` out of it: R's
## serialiser writes the shards' cases twice if the undivided list is present.

## BENCH_SEM is an environment, not a function, and the guards resolve it by
## name in .GlobalEnv at call time, so it has to travel explicitly.
payload <- c(Filter(function(n) is.function(get(n, envir = .GlobalEnv)),
                    ls(.GlobalEnv)),
             "shards", "tolerances", "modes", "cfgs", "sweep_cfgs",
             "nrep", "bench_cores", "maxs2", "min_time", "BENCH_SEM",
             "compile_slots", "max_compile_gb", "max_worker_gb")

cat(sprintf("\nworkspace payload: %s per shard\n",
            format(utils::object.size(mget(payload, envir = .GlobalEnv)),
                   units = "MB")))
cat(sprintf("tolerances (%s): %s\n",
            if ("tol" %in% SET) "explicit --tol" else paste("from tier", OPT$tier),
            paste(sprintf("atol %.0e / rtol %.0e", tolerances$atol, tolerances$rtol),
                  collapse = "; ")))
cat(sprintf("equivalent local run:\n  Rscript benchmarks/run-benchmarks.R --tier %s --modes %s --cores %d%s\n",
            OPT$tier, OPT$modes, bench_cores,
            if (tf(OPT$`sparse-sweep`))
              sprintf(" --max-density %s --min-sweep-states %s",
                      OPT$`max-density`, OPT$`min-sweep-states`)
            else " --no-sparse-sweep"))

if (!tf(OPT$submit)) {
  cat("\nDry run, nothing sent.\n")
  if ("machines" %in% SET) { cat("host probe:\n"); try(preflight(placement)) }
  else cat("Add --machines user@host to probe the hosts from here.\n")
  cat("\nThen re-run with --machines user@host --submit\n")
  quit(save = "no")
}


## ---------------------------------------------------------------------
##  Submit
## ---------------------------------------------------------------------

cat("\nhost probe:\n")
preflight(placement)

WAIT <- tf(OPT$wait)
if (WAIT && n_shards > 1L)
  warning("runbg() starts the jobs sequentially when it waits, so --wait ",
          "with ", n_shards, " shards runs them one after another",
          call. = FALSE, immediate. = TRUE)

cat(sprintf("\nstarting '%s': %d shard(s), %d worker(s) each\n",
            OPT$jobname, n_shards, bench_cores))

dir.create(JOBDIR, showWarnings = FALSE, recursive = TRUE)
saveRDS(OPT[setdiff(names(OPT), FLAGS)], REC)

setwd(STAGE)
invisible(dMod2::runbg(
  {
    ## Runs on a host. `.node` is this job's shard index, set by runbg after the
    ## workspace is loaded; `shards` and the benchmark functions arrive in it.

    ## runbg only re-attaches the packages attached here, so cppDE is loaded
    ## explicitly or every model fails and the shard returns nothing.
    library(cppDE)
    stopifnot(isTRUE(cppDE:::cvodeConfig$available))
    Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
               MKL_NUM_THREADS = "1")
    builddir <- file.path(tempdir(), paste0("cppde_bench_", .node))
    dir.create(builddir, showWarnings = FALSE, recursive = TRUE)

    ## The memory guards, set up before the workers are forked so they inherit
    ## both. The compiler cap is written into a Makevars under builddir, node-local
    ## ($TMPDIR): on a shared home it must not be a file other jobs can see.

    ## The per-worker cap is best-effort and reports back rather than failing, so a
    ## host without the `unix` package still runs; the preflight has said so.
    bench_limit_compiler(max_compile_gb, dir = builddir)
    if (max_worker_gb > 0 && !isTRUE(bench_limit_process(max_worker_gb)))
      message("note: could not set a per-worker memory limit on ",
              Sys.info()[["nodename"]])

    out <- run_all_problems(shards[[.node]], builddir = builddir,
                            tolerances = tolerances, modes = modes,
                            configs = cfgs, sweep_configs = sweep_cfgs,
                            nrep = nrep, cores = bench_cores,
                            max_sens2 = maxs2, min_time = min_time,
                            compile_slots = compile_slots)
    ## Stamp the hardware: absolute times are only comparable across rows
    ## that ran on the same host, and here nothing guaranteed they did.
    if (!is.null(out) && nrow(out)) {
      out$shard <- .node
      out$node  <- unname(Sys.info()[["nodename"]])
      out$cpu   <- tryCatch(sub(".*: ", "", grep("model name",
                     readLines("/proc/cpuinfo", warn = FALSE), value = TRUE)[1L]),
                     error = function(e) NA_character_)
    }
    out
  },
  machine  = placement,
  filename = OPT$jobname,
  input    = payload,
  wait     = WAIT))

if (WAIT) {
  merge_results(get(".runbgOutput", envir = .GlobalEnv))
} else {
  cat("\nstarted. Poll, inspect and fetch with:\n",
      sprintf("  Rscript benchmarks/run-with-runbg.R --check   --jobname %s\n", OPT$jobname),
      sprintf("  Rscript benchmarks/run-with-runbg.R --log     --jobname %s\n", OPT$jobname),
      sprintf("  Rscript benchmarks/run-with-runbg.R --collect --jobname %s\n", OPT$jobname),
      sep = "")
}
