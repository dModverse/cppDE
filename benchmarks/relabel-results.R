#!/usr/bin/env Rscript
## =====================================================================
##  relabel-results.R -- retro-fit older result folders.
##
##  Folders written before the naming change carry a bare timestamp and
##  a plain-text run-info.txt.  This reads what each run actually was
##  out of that file, renames the folder to the descriptive scheme, and
##  writes the same README.md a fresh run produces.
##
##      Rscript benchmarks/relabel-results.R            # show the plan
##      Rscript benchmarks/relabel-results.R --apply    # do it
##
##  Nothing is deleted: only directory names change, and README.md is
##  added.  run-info.txt is left in place.
## =====================================================================

suppressPackageStartupMessages(library(cppDE))

ROOT <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) normalizePath(file.path(dirname(f[1]), "..")) else getwd()
})
source(file.path(ROOT, "benchmarks", "R", "harness.R"))
bench_source(file.path(ROOT, "benchmarks", "R"))

args  <- commandArgs(trailingOnly = TRUE)
APPLY <- "--apply" %in% args
RES   <- file.path(ROOT, "benchmarks", "results")

## run-info.txt keeps the option string on one line; pull the fields the
## folder name needs back out of it.
parse_info <- function(path) {
  if (!file.exists(path)) return(list())
  txt <- readLines(path, warn = FALSE)
  field <- function(key) {
    i <- grep(paste0("^", key, "\\s*:"), txt)
    if (!length(i)) return(NULL)
    trimws(sub("^[^:]*:", "", txt[i[1L]]))
  }
  opt <- field("options")
  o <- list()
  if (!is.null(opt)) {
    kv <- strsplit(strsplit(opt, " +")[[1L]], "=", fixed = TRUE)
    for (p in kv) if (length(p) >= 1L)
      o[[p[1L]]] <- if (length(p) >= 2L) paste(p[-1L], collapse = "=") else ""
  }
  cov <- local({
    i <- grep("^trait coverage:", txt)
    if (!length(i)) return(NULL)
    rows <- txt[(i + 1L):length(txt)]
    rows <- rows[grepl("^\\s+\\S+\\s+\\d+\\s*$", rows)]
    if (!length(rows)) return(NULL)
    m <- do.call(rbind, strsplit(trimws(rows), "\\s+"))
    data.frame(trait = m[, 1L], n = as.integer(m[, 2L]),
               stringsAsFactors = FALSE)
  })
  skip <- local({
    i <- grep("^skipped:", txt)
    if (!length(i)) return(list())
    rows <- trimws(txt[(i + 1L):length(txt)])
    rows <- rows[nzchar(rows) & rows != "(none)"]
    if (!length(rows)) return(list())
    nm <- sub("\\s.*", "", rows)
    stats::setNames(as.list(trimws(sub("^\\S+\\s+", "", rows))), nm)
  })
  list(opts = o, coverage = cov, skipped = skip,
       date = field("date"), elapsed = field("elapsed"),
       cpu = field("cpu"), klu = field("KLU"),
       versions = if (!is.null(field("R")))
         sprintf("%s / %s", field("R"), field("cppDE")) else NULL,
       options = opt)
}

dirs <- list.dirs(RES, recursive = FALSE)
dirs <- Filter(function(d) file.exists(file.path(d, "results.csv")), dirs)
if (!length(dirs)) stop("no result folders with a results.csv under ", RES)

cat(sprintf("%d result folder(s) in %s\n\n", length(dirs), RES))
planned <- character(0)

for (d in dirs) {
  base <- basename(d)
  df <- utils::read.csv(file.path(d, "results.csv"), stringsAsFactors = FALSE)
  info <- parse_info(file.path(d, "run-info.txt"))
  o <- info$opts

  ## Fall back to what the data itself reveals when run-info.txt is old
  ## or absent: the modes are in the CSV, and so is the core count.
  stamp <- sub("_.*", "", base)
  cores <- if (!is.null(o$cores)) as.integer(o$cores) else
           if (!is.null(df$cores)) as.integer(df$cores[1L]) else 1L
  modes <- if (!is.null(o$modes)) o$modes else
           paste(unique(df$mode), collapse = ",")
  ## Runs made before --tier existed were driven by explicit limits, so
  ## "custom" describes them; "unknown" would suggest information lost.
  tier  <- if (!is.null(o$tier)) o$tier else "custom"
  models <- if (!is.null(o$models) && nzchar(o$models))
    trimws(strsplit(o$models, ",")[[1L]]) else character(0)

  newname <- run_tag(stamp, tier, cores, modes,
                     suite = if (is.null(o$suite)) "all" else o$suite,
                     models = models,
                     all_conditions = identical(o$conditions, "all"))
  target <- file.path(RES, newname)

  cat(sprintf("%-24s -> %s\n", base, newname))
  if (identical(base, newname)) { cat("    (name already current)\n") }
  else if (dir.exists(target) && !identical(normalizePath(target), normalizePath(d))) {
    cat("    !! target exists, leaving this folder alone\n"); next
  }
  planned <- c(planned, newname)

  if (APPLY) {
    dest <- d
    if (!identical(base, newname)) {
      if (!file.rename(d, target)) { cat("    !! rename failed\n"); next }
      dest <- target
    }
    write_run_readme(df, dest, c(list(tag = newname, tier = tier, modes = modes,
                                      cores = cores, nrep = o$nrep),
                                 info[c("coverage", "skipped", "date", "elapsed",
                                        "cpu", "klu", "versions", "options")]))
    cat("    wrote README.md\n")
  }
}

if (!APPLY)
  cat("\nDry run. Re-run with --apply to rename the folders and write the READMEs.\n")
