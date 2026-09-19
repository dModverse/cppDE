## Logs the generator calls of both test suites and of the benchmark models
## (up to max_states states) to dev/python/tap_log/.
##
##   Rscript dev/python/collect_corpus.R [max_states]

args <- commandArgs(trailingOnly = TRUE)
max_states <- if (length(args)) as.integer(args[1]) else 150L

root <- normalizePath(local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  file.path(dirname(f[1]), "..", "..")
}), "/")
logdir <- file.path(root, "dev", "python", "tap_log")
unlink(logdir, recursive = TRUE)
dir.create(logdir, recursive = TRUE)

Sys.setenv(CPPDE_PY_DIR = file.path(root, "dev", "python", "tap"),
           CPPDE_TAP_REAL = system.file("python", package = "cppDE"),
           CPPDE_TAP_LOG = logdir,
           NOT_CRAN = "true")

run_suite <- function(dir, pkg) {
  cat("== tests of", pkg, "\n")
  old <- setwd(dir); on.exit(setwd(old))
  t0 <- proc.time()[["elapsed"]]
  try(testthat::test_dir("tests/testthat", package = pkg,
                         load_package = "installed", reporter = "summary",
                         stop_on_failure = FALSE))
  cat(sprintf("   %.0f s\n", proc.time()[["elapsed"]] - t0))
}

run_suite(root, "cppDE")
run_suite(file.path(dirname(root), "dMod2"), "dMod2")

## Benchmark models, through codegen only.
cat("== benchmark models up to", max_states, "states\n")
suppressPackageStartupMessages(library(cppDE))
ROOT <- root
source(file.path(root, "benchmarks", "R", "harness.R"))
bench_source(file.path(root, "benchmarks", "R"))
petab_root <- file.path(root, "benchmarks", "cache", "petab", "Benchmark-Models")
out <- file.path(tempdir(), "corpus_models")
dir.create(out, showWarnings = FALSE)
if (dir.exists(petab_root)) {
  lst <- petab_list(petab_root)
  for (i in seq_len(nrow(lst))) {
    n <- petab_species_count(lst$yaml[i])
    if (is.na(n) || n > max_states) next
    prob <- tryCatch(petab_problem(lst$yaml[i]), error = function(e) NULL)
    if (is.null(prob)) next
    cat(sprintf("   %-36s %4d states\n", lst$name[i], length(prob$rhs)))
    try(cppODE(prob$rhs, events = prob$events, fixed = prob$fixed,
               deriv = TRUE, compile = FALSE, outdir = out,
               modelname = paste0("corpus_", i)))
  }
}
for (prob in tryCatch(classic_problems(), error = function(e) list())) {
  if (length(prob$rhs) > max_states) next
  try(cppODE(prob$rhs, events = prob$events, deriv = TRUE, compile = FALSE,
             outdir = out, modelname = paste0("corpus_", prob$id)))
}
cat("log in", logdir, "\n")
