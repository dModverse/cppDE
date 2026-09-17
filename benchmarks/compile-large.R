#!/usr/bin/env Rscript

## =====================================================================
##  compile-large.R: code generation and compile time of large models.
## =====================================================================

##      Rscript benchmarks/compile-large.R
##      Rscript benchmarks/compile-large.R --min-states 100 --models Chen,LLG

##  For each model and derivative mode (forward, reverse):
##    codegen   cppODE(compile = FALSE)
##    cppDE     cppDE::compile() of that model
##    dMod2     dMod2::odemodel(compile = FALSE) + dMod2::compile(cores = 1),
##              which builds the value and the derivative model
##  then one solve over the first five output times, checked against central
##  differences (forward) or against the forward sensitivities (reverse).
##  Everything runs serially.
##  Results land in benchmarks/results/compile-large-<stamp>.csv.

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

OPTS <- list(`min-states` = "100", models = "", modes = "forward,reverse",
             dmod2 = "TRUE", sizes = "500",
             `petab-root` = file.path(ROOT, "benchmarks", "cache", "petab",
                                      "Benchmark-Models"),
             outdir = file.path(ROOT, "benchmarks", "results"))
args <- commandArgs(trailingOnly = TRUE)
for (i in seq_len(length(args) %/% 2L)) {
  key <- sub("^--", "", args[2L * i - 1L])
  if (!key %in% names(OPTS)) stop("unknown option --", key)
  OPTS[[key]] <- args[2L * i]
}
min_states <- as.integer(OPTS$`min-states`)
modes <- strsplit(OPTS$modes, ",")[[1L]]
use_dmod2 <- as.logical(OPTS$dmod2) && requireNamespace("dMod2", quietly = TRUE)
pick <- if (nzchar(OPTS$models)) strsplit(OPTS$models, ",")[[1L]] else character(0)
tol <- list(abstol = 1e-10, reltol = 1e-10)
chk_times <- function(prob) head(prob$times, 5L)

## ---------------------------------------------------------------------
##  Problems
## ---------------------------------------------------------------------

problems <- list()
sizes <- as.integer(strsplit(OPTS$sizes, ",")[[1L]])
for (n in sizes) local({
  nn <- n
  problems[[length(problems) + 1L]] <<- function() build_brusselator(nn)
  problems[[length(problems) + 1L]] <<- function() build_fhn_chain(nn)
})
problems[[length(problems) + 1L]] <- function() build_llg_dipole(160L)
if (dir.exists(OPTS$`petab-root`)) {
  pl <- petab_list(OPTS$`petab-root`)
  nsp <- vapply(pl$yaml, petab_species_count, 0L)
  for (y in pl$yaml[!is.na(nsp) & nsp >= min_states])
    local({ yy <- y; problems[[length(problems) + 1L]] <<- function() petab_problem(yy) })
} else {
  message("no PEtab collection at ", OPTS$`petab-root`, "; run fetch-models.R")
}

## ---------------------------------------------------------------------
##  One model
## ---------------------------------------------------------------------

stamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
build <- file.path(tempdir(), paste0("compile-large-", stamp))
dir.create(build, recursive = TRUE, showWarnings = FALSE)
dir.create(OPTS$outdir, recursive = TRUE, showWarnings = FALSE)
out_csv <- file.path(OPTS$outdir, paste0("compile-large-", stamp, ".csv"))

elapsed <- function(expr) {
  t0 <- proc.time()[["elapsed"]]
  force(expr)
  proc.time()[["elapsed"]] - t0
}

relerr <- function(a, b) max(abs(a - b)) / max(1, max(abs(b)))

jac_note <- function(src) {
  l <- grep("^// Jacobian \\(", readLines(src, warn = FALSE), value = TRUE)[1L]
  if (is.na(l)) return(c(strategy = NA, nnz = NA))
  c(strategy = sub("^// Jacobian \\([a-z]+, ([a-z]+),.*$", "\\1", l),
    nnz = sub("^.* ([0-9]+) nonzeros.*$", "\\1", l))
}

check_forward <- function(m, prob) {
  tt <- chk_times(prob)
  r <- do.call(solveODE, c(list(m, tt, prob$parms), tol))
  ps <- head(intersect(dimnames(r$sens1)[[3]], prob$sens), 2L)
  err <- 0
  for (p in ps) {
    h <- 1e-6 * max(1, abs(prob$parms[[p]]))
    pp <- prob$parms; pp[[p]] <- pp[[p]] + h
    pm <- prob$parms; pm[[p]] <- pm[[p]] - h
    fd <- (do.call(solveODE, c(list(m, tt, pp), tol))$variable -
           do.call(solveODE, c(list(m, tt, pm), tol))$variable) / (2 * h)
    err <- max(err, relerr(r$sens1[, , p], fd))
  }
  c(finite = all(is.finite(r$variable)), check = err)
}

check_reverse <- function(m, mf, prob) {
  tt <- chk_times(prob)
  rf <- do.call(solveODE, c(list(mf, tt, prob$parms), tol))
  set.seed(1)
  W <- array(rnorm(prod(dim(rf$variable))), c(dim(rf$variable), 1L))
  rr <- do.call(solveODE, c(list(m, tt, prob$parms, seed = W), tol))
  ref <- apply(rf$sens1 * as.vector(W[, , 1L]), 3, sum)
  c(finite = all(is.finite(rr$adjoint)),
    check = relerr(rr$adjoint[names(ref), 1L], ref))
}

run_one <- function(prob) {
  rows <- list()
  fwd <- NULL
  for (mode in modes) {
    tag <- sprintf("%s_%s", substr(gsub("[^A-Za-z0-9]", "", prob$id), 1, 20),
                   substr(mode, 1, 3))
    row <- list(model = prob$name, states = prob$nstates, sens = prob$nsens,
                mode = mode)
    m <- NULL
    row$codegen_s <- elapsed(m <- cppODE(prob$rhs, events = prob$events,
                                         derivMode = mode, fixed = prob$fixed,
                                         compile = FALSE, outdir = build,
                                         modelname = tag))
    src <- attr(m, "srcfile")
    row$source_kB <- round(file.size(src) / 1e3)
    jn <- jac_note(src)
    row$strategy <- jn[["strategy"]]
    row$nnz <- jn[["nnz"]]
    row$cppDE_compile_s <- elapsed(cppDE::compile(m))
    chk <- tryCatch(
      if (mode == "forward") check_forward(m, prob)
      else if (!is.null(fwd)) check_reverse(m, fwd, prob)
      else c(finite = NA, check = NA),
      error = function(e) { message("  check failed: ", conditionMessage(e));
                            c(finite = FALSE, check = NA) })
    row$finite <- as.logical(chk[["finite"]])
    row$check <- signif(chk[["check"]], 3)
    if (mode == "forward") fwd <- m
    if (use_dmod2) {
      dtag <- paste0("d", tag)
      row$dMod2_compile_s <- tryCatch({
        om <- dMod2::odemodel(prob$rhs, deriv = TRUE, derivMode = mode,
                              fixed = prob$fixed, events = prob$events,
                              modelname = dtag, outdir = build, compile = FALSE)
        elapsed(dMod2::compile(om, cores = 1))
      }, error = function(e) { message("  dMod2 failed: ", conditionMessage(e)); NA })
    }
    rows[[length(rows) + 1L]] <- as.data.frame(row, stringsAsFactors = FALSE)
    cat(sprintf("  %-8s codegen %7.1f s  cppDE %7.1f s  dMod2 %7.1f s  %6d kB  %s %s  check %.1e\n",
                mode, row$codegen_s, row$cppDE_compile_s,
                if (is.null(row$dMod2_compile_s)) NA_real_ else row$dMod2_compile_s,
                as.integer(row$source_kB), row$strategy, row$nnz, row$check))
  }
  do.call(rbind, rows)
}

## ---------------------------------------------------------------------
##  Run
## ---------------------------------------------------------------------

cat(sprintf("compile-large: %d candidate models, modes %s, dMod2 %s\n",
            length(problems), paste(modes, collapse = "/"), use_dmod2))
results <- NULL
for (mk in problems) {
  prob <- tryCatch(mk(), error = function(e) { message("load failed: ",
                                                       conditionMessage(e)); NULL })
  if (is.null(prob)) next
  if (prob$nstates < min_states) next
  if (length(pick) && !any(vapply(pick, grepl, NA, x = prob$name, fixed = TRUE))) next
  if (!isTRUE(prob$usable)) { cat(prob$name, ": excluded\n"); next }
  cat(sprintf("%s: %d states, %d sensitivities\n", prob$name, prob$nstates, prob$nsens))
  res <- tryCatch(run_one(prob), error = function(e) {
    message("  failed: ", conditionMessage(e)); NULL })
  if (!is.null(res)) {
    results <- rbind(results, res)
    utils::write.csv(results, out_csv, row.names = FALSE)
  }
}
cat("results in", out_csv, "\n")
