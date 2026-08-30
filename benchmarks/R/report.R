## =====================================================================
##  report.R -- naming and per-run documentation of a results folder.
##
##  A results tree is browsed far more often than it is opened, so the
##  directory name carries what the run *was* -- depth, parallelism,
##  modes -- not only when it happened.  The README inside carries what
##  it found.  Both are generated from the same code here so that a
##  fresh run and a retro-fitted old one describe themselves alike.
## =====================================================================

## <timestamp>_<tier>[_<suite>]_c<cores>_<modes>[_<models>][_allcond]
run_tag <- function(stamp, tier, cores, modes, suite = "all",
                    models = character(0), all_conditions = FALSE) {
  paste0(
    stamp, "_", tier,
    if (!identical(suite, "all")) paste0("_", suite) else "",
    "_c", cores,
    "_", gsub(",", "-", modes),
    if (length(models))
      paste0("_", substr(gsub("[^A-Za-z0-9]", "", paste(models, collapse = "")),
                         1L, 20L)) else "",
    if (isTRUE(all_conditions)) "_allcond" else "")
}

FIGURE_LEGEND <- c(
  "`01`/`02-work-precision` | achieved error against cost, one panel per problem; down-and-left is better",
  "`03-speedup` | per-problem ratio against CVODE, bars growing from parity",
  "`04-scaling` | cost against number of states, log-log",
  "`05-sens-overhead` | gradient cost in units of plain solves; grey line is the finite-difference cost",
  "`06-summary` | one geometric-mean number per mode",
  "`07-sens2-cost` | Hessian cost against M (cppDE only)")

## Appended when the run also pinned the linear solver both ways.
SPARSE_FIGURE_LEGEND <- c(
  "`08-sparse-gain` | dense/sparse ratio per model and backend, bars growing from parity",
  "`09-sparse-crossover` | the same ratios against system size")


## `info` carries whatever is known about the run; every field is
## optional so that a folder retro-fitted from an old run still gets a
## README rather than an error.
write_run_readme <- function(df, outdir, info = list()) {
  g <- function(k, default = NULL) if (is.null(info[[k]])) default else info[[k]]
  line <- function(label, value) if (is.null(value)) NULL else
    sprintf("| %s | %s |", label, value)

  ## -- headline ---------------------------------------------------------
  head_tbl <- local({
    s <- tryCatch(speedup_table(df), error = function(e) NULL)
    if (is.null(s) || !nrow(s)) return("_no matched cells_")
    out <- c("| mode | solver | time vs CVODE | rhs-evals vs CVODE | cells | problems |",
             "|---|---|---:|---:|---:|---:|")
    for (md in unique(s$mode)) for (sv in setdiff(unique(s$solver), "CVODE_bdf")) {
      x <- s[s$mode == md & s$solver == sv, , drop = FALSE]
      if (!nrow(x)) next
      out <- c(out, sprintf("| %s | %s | **%.2f×** | %.2f× | %d | %d |",
                            MODE_LABELS[md] %||% md, sv, geo_mean(x$speedup),
                            geo_mean(x$fev_ratio), nrow(x),
                            length(unique(x$problem))))
    }
    out
  })

  ## -- sparse sweep: dense against sparse, per backend ------------------
  ## Empty unless the run pinned the linear solver both ways.
  ## Aggregated over tolerances -- per-cell would be several hundred rows
  ## on the full tier, and the per-cell numbers are in results.csv.
  ## `chose` is what auto-detection actually did on that model, so the
  ## table can be read as a verdict on the choice and not only as a cost.
  sparse_tbl <- local({
    sg <- tryCatch(sparse_gain_table(df), error = function(e) NULL)
    if (is.null(sg) || !nrow(sg)) return(NULL)
    agg <- stats::aggregate(gain ~ problem + nstates + backend + mode,
                            data = sg, FUN = geo_mean)
    agg <- agg[order(agg$nstates, agg$problem, agg$backend, agg$mode), ]
    auto <- df[pinned_col(df) == "auto" & df$ok, , drop = FALSE]
    chose_of <- local({
      tab <- if (nrow(auto))
        tapply(auto$lu, paste(auto$problem, auto$backend), function(x) x[1L])
        else character(0)
      function(p, b) { v <- tab[paste(p, b)]; if (is.na(v)) "?" else unname(v) }
    })
    out <- c("", "## Dense vs sparse LU", "",
             sprintf(paste("%d model(s) run a second and third time with the linear",
                           "solver pinned, next to the auto-detected head-to-head."),
                     length(unique(agg$problem))),
             "Gain > 1 means the sparse path was faster; `chose` is what",
             "auto-detection picked, so a gain > 1 next to `dense` is a",
             "mis-detection.", "",
             "| problem | states | backend | mode | chose | gain |",
             "|---|---:|---|---|---|---:|")
    for (i in seq_len(nrow(agg)))
      out <- c(out, sprintf("| %s | %d | %s | %s | %s | **%.2f×** |",
                            agg$problem[i], agg$nstates[i], agg$backend[i],
                            MODE_LABELS[agg$mode[i]] %||% agg$mode[i],
                            chose_of(agg$problem[i], agg$backend[i]),
                            agg$gain[i]))
    out
  })

  ## -- second order has no baseline: report cost, not a ratio -----------
  sens2_tbl <- local({
    s2 <- df[df$mode == "sens2" & df$ok, , drop = FALSE]
    if (!nrow(s2)) return(NULL)
    pl <- stats::aggregate(time_ms ~ problem,
                           data = df[df$mode == "nosens" & df$ok, , drop = FALSE],
                           FUN = stats::median)
    names(pl)[2] <- "plain"
    m <- merge(stats::aggregate(cbind(time_ms, nsens) ~ problem, data = s2,
                                FUN = stats::median), pl, by = "problem")
    if (!nrow(m)) return(NULL)
    m <- m[order(-m$time_ms / m$plain), ]
    c("", "**Second order** — cppDE only, CVODES has no second-order",
      "sensitivities, so this is a cost, not a comparison:", "",
      "| problem | M | plain [ms] | Hessian [ms] | factor |",
      "|---|---:|---:|---:|---:|",
      sprintf("| %s | %d | %.2f | %.1f | %.0f× |", m$problem, as.integer(m$nsens),
              m$plain, m$time_ms, m$time_ms / m$plain))
  })

  ## -- what ran ---------------------------------------------------------
  prob_tbl <- local({
    d <- unique(df[, c("problem", "source", "nstates", "npars", "nout")])
    cnd <- tapply(df$condition, df$problem, function(x) length(unique(x)))
    sns <- tapply(df$nsens, df$problem, max)
    d$conditions <- as.integer(cnd[d$problem])
    d$sens <- as.integer(sns[d$problem])
    d <- d[order(d$nstates, d$problem), ]
    c("| problem | source | states | params | sens | conditions | out |",
      "|---|---|---:|---:|---:|---:|---:|",
      sprintf("| %s | %s | %d | %d | %d | %d | %d |", d$problem, d$source,
              d$nstates, d$npars, d$sens, d$conditions, d$nout))
  })

  cov <- g("coverage")
  cov_tbl <- if (is.null(cov)) NULL else c(
    "", "## Trait coverage", "",
    "A tier is only a usable regression basis while every trait still has",
    "a representative.", "",
    "| trait | problems |", "|---|---:|",
    sprintf("| %s | %d%s |", cov$trait, cov$n,
            ifelse(cov$n == 0L, " **(gap)**", "")))

  skipped <- g("skipped", list())
  cores <- g("cores", unique(df$cores)[1L])

  md <- c(
    sprintf("# Benchmark run — `%s`", g("tag", basename(outdir))), "",
    sprintf("cppDE against SUNDIALS CVODE(S). %d rows, %d failed cells.",
            nrow(df), sum(!df$ok)),
    if (!is.null(cores) && !is.na(cores) && cores > 1L)
      paste("\n> Run with", cores, "workers: **absolute times are inflated**",
            "by shared cache and turbo clocking. The ratios below are the",
            "quantity to read; treat differences under ~10 % as noise.") else NULL,
    "", "## Result", "", head_tbl, sparse_tbl, sens2_tbl, "",
    "## Configuration", "", "| | |", "|---|---|",
    line("tier", g("tier")), line("modes", g("modes")),
    line("cores", if (is.null(cores)) NULL else
      sprintf("%d%s", cores, if (cores > 1L) " (parallel)" else " (serial)")),
    line("tolerances", g("tolerances")),
    line("repetitions", g("nrep")),
    line("elapsed", g("elapsed")),
    line("date", g("date")),
    line("R / cppDE", g("versions")),
    line("CPU", g("cpu")),
    line("KLU", g("klu")),
    if (!is.null(g("options"))) c("", sprintf("Full options: `%s`", g("options"))) else NULL,
    "", "## Problems", "", prob_tbl,
    cov_tbl,
    "", "## Skipped", "",
    if (length(skipped)) sprintf("- `%s` — %s", names(skipped), unlist(skipped))
    else "_none_",
    "", "## Figures", "", "| file | shows |", "|---|---|",
    sprintf("| %s |", c(if (any(pinned_col(df) == "auto")) FIGURE_LEGEND,
                        if (is_sparse_sweep(df)) SPARSE_FIGURE_LEGEND)), "",
    "---", "",
    "Raw numbers are in `results.csv`: one row per problem / condition /",
    "solver / mode / tolerance.")

  writeLines(md[!vapply(md, is.null, NA)], file.path(outdir, "README.md"))
  invisible(file.path(outdir, "README.md"))
}
