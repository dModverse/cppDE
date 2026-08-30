#!/usr/bin/env Rscript
## =====================================================================
##  fetch-models.R -- download (or update) the PEtab benchmark collection.
##
##  Hass et al. (2019), "Benchmark problems for dynamic modeling of
##  intracellular processes", Bioinformatics 35(17):3073-3082.
##  https://doi.org/10.1093/bioinformatics/btz020
##
##  The collection is fetched into benchmarks/cache/petab, which is not
##  tracked by git -- the models stay under their own licence and are
##  never vendored into this repository.
## =====================================================================

ROOT <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) normalizePath(file.path(dirname(f[1]), "..")) else getwd()
})

REPO  <- "https://github.com/Benchmarking-Initiative/Benchmark-Models-PEtab.git"
DEST  <- file.path(ROOT, "benchmarks", "cache", "petab")

if (Sys.which("git") == "")
  stop("git is required to fetch the benchmark collection")

dir.create(dirname(DEST), showWarnings = FALSE, recursive = TRUE)

if (dir.exists(file.path(DEST, ".git"))) {
  cat("updating existing clone in", DEST, "\n")
  st <- system2("git", c("-C", shQuote(DEST), "pull", "--ff-only", "--quiet"))
} else {
  cat("cloning", REPO, "\n  into", DEST, "\n")
  st <- system2("git", c("clone", "--depth", "1", "--quiet", REPO, shQuote(DEST)))
}
if (st != 0L) stop("git failed with status ", st)

models <- list.dirs(file.path(DEST, "Benchmark-Models"), recursive = FALSE)
cat(sprintf("\n%d model directories available in %s\n",
            length(models), file.path(DEST, "Benchmark-Models")))
cat("next:  Rscript benchmarks/run-benchmarks.R --quick\n")
