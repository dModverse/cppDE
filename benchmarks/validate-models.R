#!/usr/bin/env Rscript
## =====================================================================
##  validate-models.R -- check the SBML translation against ground truth.
##
##  The PEtab collection ships a `simulatedData_*.tsv` for every model:
##  the observables as produced by an independent, established
##  implementation at the nominal parameter values.  Reproducing those
##  numbers exercises the whole chain at once -- MathML parsing,
##  stoichiometry, compartment scaling, assignment rules, initial
##  assignments, condition overrides, code generation and the
##  integrator.  Agreement between cppDE and CVODE cannot show this,
##  because both are fed the same translated right-hand side.
##
##  Usage:
##      Rscript benchmarks/validate-models.R
##      Rscript benchmarks/validate-models.R --models Boehm,Fujita
##
##  Models whose PEtab problem specifies a steady-state
##  pre-equilibration are reported as NOT COMPARABLE rather than as
##  failures: the benchmark starts from the SBML initials by design, so
##  a mismatch there says nothing about the translation.
## =====================================================================

suppressPackageStartupMessages(library(cppDE))

ROOT <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) normalizePath(file.path(dirname(f[1]), "..")) else getwd()
})
source(file.path(ROOT, "benchmarks", "R", "harness.R"))
bench_source(file.path(ROOT, "benchmarks", "R"))

args <- commandArgs(trailingOnly = TRUE)
get_opt <- function(name, default = "") {
  i <- which(args == paste0("--", name))
  if (length(i) && i[1] < length(args)) return(args[i[1] + 1L])
  kv <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(kv)) return(sub(paste0("^--", name, "="), "", kv[1]))
  default
}
PETAB <- get_opt("petab-root",
                 file.path(ROOT, "benchmarks", "cache", "petab", "Benchmark-Models"))
filt <- trimws(strsplit(get_opt("models"), ",")[[1L]]); filt <- filt[nzchar(filt)]
TOL  <- as.numeric(get_opt("tol", "1e-3"))

builddir <- file.path(tempdir(), "cppDE_validate")
dir.create(builddir, showWarnings = FALSE, recursive = TRUE)


## ---------------------------------------------------------------------
##  Observable formulas
## ---------------------------------------------------------------------

## PEtab writes observable formulas in infix notation over the original
## SBML ids, which may not be valid R names, so the same renaming the
## model went through is applied here.  `observableParameterN_<id>`
## placeholders are filled from the measurement table.
prepare_formula <- function(formula, ids, obs_params) {
  f <- gsub("**", "^", formula, fixed = TRUE)
  ren <- ids[san_id(ids) != ids]
  for (old in ren[order(-nchar(ren))])
    f <- gsub(paste0("(?<![A-Za-z0-9._])", old, "(?![A-Za-z0-9._])"),
              san_id(old), f, perl = TRUE)
  for (i in seq_along(obs_params))
    f <- gsub(paste0("(?<![A-Za-z0-9._])observableParameter", i, "_[A-Za-z0-9_]+"),
              sprintf("(%s)", obs_params[i]), f, perl = TRUE)
  f
}

validate_one <- function(yaml_path) {
  loaded <- petab_load(yaml_path)
  pt <- loaded$pt
  sim <- pt$simulated
  if (is.null(sim)) return(list(status = "no simulatedData"))
  val_col <- intersect(c("simulation", "simulatedData", "measurement"), names(sim))
  if (!length(val_col)) return(list(status = "no simulation column"))
  val_col <- val_col[1L]

  ## Pre-equilibration is out of scope for the benchmark (see README).
  preeq <- "preequilibrationConditionId" %in% names(sim) &&
           any(!is.na(sim$preequilibrationConditionId) &
               nzchar(sim$preequilibrationConditionId))
  cond <- petab_conditions(loaded)$conditionId[1L]
  prob <- petab_case(loaded, condition = if (is.na(cond)) NULL else cond)
  if (!isTRUE(prob$usable)) return(list(status = "excluded (state-dependent switch)"))

  rows <- sim
  if (!is.na(cond) && "simulationConditionId" %in% names(sim))
    rows <- sim[sim$simulationConditionId == cond, , drop = FALSE]
  rows$.t <- suppressWarnings(as.numeric(rows$time))
  rows <- rows[is.finite(rows$.t), , drop = FALSE]
  if (!nrow(rows)) return(list(status = "steady-state data only"))

  ## Solve on the union of the model grid and the data times.
  times <- sort(unique(c(0, prob$times, rows$.t)))
  cache <- new_model_cache(builddir)
  on.exit(release_cache(cache), add = TRUE)
  m <- get_model(cache, prob, "cppde", deriv = FALSE)
  res <- solveODE(m, times, prob$parms, abstol = 1e-11, reltol = 1e-9,
                  onFailure = "silent")
  if (is.null(res$diagnostics) || res$diagnostics$return_code != 0L)
    return(list(status = "solve failed"))

  obs_tab <- pt$observables
  ids <- c(loaded$sbml$species$id, loaded$sbml$parameters$id,
           loaded$sbml$compartments$id)
  ## san_id() was already applied on read, so recover the originals from
  ## the file to build the rename map.
  raw_ids <- unlist(lapply(c("species", "parameter", "compartment"), function(k)
    xml2::xml_attr(xml2::xml_find_all(
      { d <- xml2::read_xml(pt$sbml_file); xml2::xml_ns_strip(d); d },
      paste0(".//", k)), "id")))
  raw_ids <- raw_ids[!is.na(raw_ids)]

  ## Observable formulas may use parameters that never appear in the
  ## ODE (scaling and offset factors), so they are not in prob$parms.
  ## Layer the PEtab nominal values and the SBML defaults underneath.
  sbml_defaults <- stats::setNames(loaded$sbml$parameters$value,
                                   loaded$sbml$parameters$id)
  env0 <- c(as.list(prob$parms), as.list(loaded$nominal),
            as.list(sbml_defaults))
  env0 <- env0[!duplicated(names(env0))]
  worst <- 0; n_cmp <- 0L; failures <- character(0)

  for (oid in unique(rows$observableId)) {
    orow <- obs_tab[obs_tab$observableId == oid, , drop = FALSE]
    if (!nrow(orow)) next
    sub <- rows[rows$observableId == oid, , drop = FALSE]
    op <- if ("observableParameters" %in% names(sub)) sub$observableParameters[1L] else NA
    op <- if (is.na(op) || !nzchar(op)) character(0) else
            trimws(strsplit(op, ";", fixed = TRUE)[[1L]])
    op <- vapply(op, function(x) {
      v <- suppressWarnings(as.numeric(x))
      if (!is.na(v)) num_str(v) else san_id(x)
    }, "")
    f <- prepare_formula(orow$observableFormula[1L], raw_ids, op)
    e <- tryCatch(str2lang(f), error = function(err) NULL)
    if (is.null(e)) { failures <- c(failures, paste0(oid, ": unparseable")); next }

    idx <- match(sub$.t, res$time)
    keep <- !is.na(idx)
    if (!any(keep)) next
    st <- as.data.frame(res$variable[idx[keep], , drop = FALSE])
    pred <- tryCatch(eval(e, c(as.list(st), env0)), error = function(err) NULL)
    if (is.null(pred) || length(pred) != sum(keep)) {
      failures <- c(failures, paste0(oid, ": not evaluable")); next
    }
    truth <- suppressWarnings(as.numeric(sub[[val_col]][keep]))
    ok <- is.finite(pred) & is.finite(truth)
    if (!any(ok)) next
    scale <- max(abs(truth[ok]), 1e-30)
    worst <- max(worst, max(abs(pred[ok] - truth[ok])) / scale)
    n_cmp <- n_cmp + sum(ok)
  }

  list(status = if (n_cmp == 0L) "nothing comparable"
               else if (preeq) "pre-equilibration (not comparable)"
               else if (worst <= TOL) "PASS" else "MISMATCH",
       worst = worst, n = n_cmp, preeq = preeq,
       failures = failures)
}


## ---------------------------------------------------------------------
##  Run
## ---------------------------------------------------------------------

if (!dir.exists(PETAB))
  stop("PEtab collection not found at ", PETAB,
       "\n  run: Rscript benchmarks/fetch-models.R")

idx <- petab_list(PETAB)
if (length(filt))
  idx <- idx[Reduce(`|`, lapply(filt, function(f)
    grepl(f, idx$name, fixed = TRUE))), , drop = FALSE]

cat(sprintf("Validating %d model(s) against simulatedData (tol %.0e)\n\n",
            nrow(idx), TOL))

tally <- character(0)
for (i in seq_len(nrow(idx))) {
  r <- tryCatch(validate_one(idx$yaml[i]),
                error = function(e) list(status = paste("ERROR:",
                                          sub("\n.*", "", conditionMessage(e)))))
  tally <- c(tally, r$status)
  detail <- if (!is.null(r$worst) && r$n > 0)
    sprintf("  max rel. dev %.2e over %d points", r$worst, r$n) else ""
  cat(sprintf("%-38s %-32s%s\n", idx$name[i], r$status, detail))
  for (f in r$failures) cat("      note: ", f, "\n", sep = "")
  flush(stdout())
}

cat("\n--- summary ---\n")
for (s in names(sort(table(tally), decreasing = TRUE)))
  cat(sprintf("  %-34s %d\n", s, sum(tally == s)))
