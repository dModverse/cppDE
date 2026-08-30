## =====================================================================
##  tiers.R -- three benchmark depths, and the coverage rule that makes
##  the smallest one trustworthy.
##
##  The point of `tiny` is to answer "did this change help or hurt?" in a
##  couple of minutes.  That only works if every *kind* of problem is
##  represented: a tier that happens to contain no sparse system cannot
##  detect a regression in the sparse linear solver, and one with no
##  event cannot detect a regression in event handling.
##
##  Tier membership is therefore an explicit, stable list -- so numbers
##  stay comparable across runs -- and `tier_coverage()` checks that the
##  list still touches every trait, failing loudly when it does not.
## =====================================================================

## The traits are the axes along which solver performance actually
## differs.  Each one exercises a distinct code path or numerical regime.
BENCH_TRAITS <- c(
  "stiff-extreme"  = "scale separation over many decades (Robertson, E5)",
  "stiff-moderate" = "stiff, but without extreme scale separation",
  "oscillatory"    = "sustained limit cycle -- error accumulates as phase drift",
  "relaxation"     = "stiffness switches on and off along the trajectory",
  "sparse"         = "sparse or banded Jacobian -- exercises the KLU path",
  "large"          = ">= 100 states",
  "events"         = "timed discontinuities -- stop, jump, restart",
  "many-sens"      = ">= 20 sensitivity parameters -- wide AD",
  "few-sens"       = "<= 5 sensitivity parameters -- narrow AD",
  "rational"       = "Michaelis-Menten / Hill kinetics (division by states)",
  "transcendental" = "exp / log / trigonometric terms in the RHS",
  "log-horizon"    = "output grid spanning >= 6 decades in time",
  "systems-biology"= "published model with realistic parameter values")


## ---------------------------------------------------------------------
##  Trait inference
## ---------------------------------------------------------------------

## Derived from the assembled problem, so PEtab models are classified
## without a hand-maintained table.  The regime traits (stiffness,
## oscillation) cannot be read off the equations and are set explicitly
## where they apply.
infer_traits <- function(p) {
  tr <- character(0)
  rhs <- p$rhs
  states <- names(rhs)

  if (p$nstates >= 100L) tr <- c(tr, "large")
  if (isTRUE(p$nevents > 0L)) tr <- c(tr, "events")
  if (p$nsens >= 20L) tr <- c(tr, "many-sens")
  if (p$nsens > 0L && p$nsens <= 5L) tr <- c(tr, "few-sens")
  if (identical(p$source, "petab")) tr <- c(tr, "systems-biology")

  tt <- p$times[p$times > 0]
  if (length(tt) > 1L && max(tt) / min(tt) >= 1e6) tr <- c(tr, "log-horizon")

  ## Structural traits, read off the expression trees once.
  langs <- lapply(unname(rhs), function(s) tryCatch(str2lang(s), error = function(e) NULL))
  langs <- Filter(Negate(is.null), langs)
  if (length(langs)) {
    calls <- unique(unlist(lapply(langs, lang_call_heads)))
    if (any(c("exp", "log", "log10", "sin", "cos", "tan", "tanh") %in% calls))
      tr <- c(tr, "transcendental")
    ## A division whose denominator mentions a state is a saturating
    ## term; a division by parameters alone is just a rescaling.
    if (any(vapply(langs, function(e) lang_divides_by(e, states), NA)))
      tr <- c(tr, "rational")

    ## The codegen switches to a sparse matrix well below 10 %.
    d <- jac_density_proxy(p)
    if (p$nstates >= 20L && !is.na(d) && d < 0.1) tr <- c(tr, "sparse")
  }
  unique(tr)
}


## Structural density of the Jacobian, counted as "equation i mentions
## state j".  An upper bound on the true density.  NA if the equations do
## not parse.
jac_density_proxy <- function(p) {
  states <- names(p$rhs)
  langs <- lapply(unname(p$rhs), function(s) tryCatch(str2lang(s), error = function(e) NULL))
  langs <- Filter(Negate(is.null), langs)
  if (!length(langs) || !p$nstates) return(NA_real_)
  refs <- vapply(langs, function(e) sum(lang_symbols(e) %in% states), 0L)
  sum(refs) / (p$nstates^2)
}

lang_call_heads <- function(e) {
  out <- character(0)
  walk <- function(x) {
    if (!is.call(x)) return(invisible())
    if (is.symbol(x[[1L]])) out[[length(out) + 1L]] <<- as.character(x[[1L]])
    for (i in seq_along(x)) if (i > 1L && !is.null(x[[i]])) walk(x[[i]])
    invisible()
  }
  walk(e); unique(out)
}

lang_divides_by <- function(e, states) {
  found <- FALSE
  walk <- function(x) {
    if (found || !is.call(x)) return(invisible())
    if (is.symbol(x[[1L]]) && identical(as.character(x[[1L]]), "/") && length(x) == 3L)
      if (any(lang_symbols(x[[3L]]) %in% states)) { found <<- TRUE; return(invisible()) }
    for (i in seq_along(x)) if (i > 1L && !is.null(x[[i]])) walk(x[[i]])
    invisible()
  }
  walk(e); found
}


## ---------------------------------------------------------------------
##  The tiers
## ---------------------------------------------------------------------

## Names are keys of classic_problems() and PEtab directory names.
## Entries are chosen as the *cheapest* representative of the traits they
## contribute, so `tiny` stays inside a few minutes.
BENCH_TIERS <- list(

  ## ~5 min.  One representative per trait -- the regression basis.
  ## The members were chosen as the cheapest carrier of the traits they
  ## contribute; the per-solve costs in the comments come from a measured
  ## run, so a future edit can keep the tier inside its budget.
  tiny = list(
    classic = c("robertson",        # stiff-extreme, log-horizon, few-sens
                "orego",            # oscillatory                      2.3 ms
                "vdp",              # relaxation oscillator            2.0 ms
                "brusselator_small",# sparse, banded Jacobian          2.0 ms
                "brusselator"),     # + large (128 states)             2.2 ms
    petab   = c("Boehm_JProteomeRes2014",    # transcendental          0.5 ms
                "Crauste_CellSystems2017",   # numerically hard canary 0.6 ms
                "Sneyd_PNAS2002",            # 6-state signalling      0.4 ms
                "Borghans_BiophysChem1997",  # rational + many-sens    1.7 ms
                "Beer_MolBioSystems2014"),   # timed event             0.4 ms
    tol = "loose_tight",
    max_states = 200L, max_sens = 24L),

  ## ~30 min.  At least two representatives per trait, plus mid-size
  ## models and a wider sensitivity range.
  medium = list(
    classic = c("robertson", "hires", "orego", "e5", "vdp", "pollution",
                "brusselator_small", "brusselator", "fhn_small"),
    petab   = c("Boehm_JProteomeRes2014", "Crauste_CellSystems2017",
                "Sneyd_PNAS2002", "Beer_MolBioSystems2014",
                "Borghans_BiophysChem1997",
                "Fujita_SciSignal2010", "Elowitz_Nature2000",
                "Bachmann_MSB2011", "Isensee_JCB2018",
                "Weber_BMC2015", "Raia_CancerResearch2011",
                "Lucarelli_CellSystems2018"),
    tol = "default",
    max_states = 200L, max_sens = 32L),

  ## Everything that fits the size limit, every trait many times over.
  full = list(
    classic = NULL,     # all
    petab   = NULL,     # all
    tol = "default",
    max_states = 400L, max_sens = 64L)
)


## ---------------------------------------------------------------------
##  Coverage
## ---------------------------------------------------------------------

## Called after the tier's problems are assembled.  A trait with zero
## representatives means the tier cannot detect a regression in that
## area, which is a defect in the tier definition, not in the run.
tier_coverage <- function(problems, tier_name = "") {
  flat <- unlist(problems, recursive = FALSE)
  tr <- lapply(flat, function(p) p$traits %||% character(0))
  n <- vapply(names(BENCH_TRAITS), function(t) sum(vapply(tr, function(x) t %in% x, NA)), 0L)
  data.frame(trait = names(BENCH_TRAITS), n = unname(n),
             description = unname(BENCH_TRAITS), stringsAsFactors = FALSE)
}

print_tier_coverage <- function(cov, tier_name = "") {
  cat(sprintf("\nTrait coverage%s:\n", if (nzchar(tier_name)) paste0(" [", tier_name, "]") else ""))
  for (i in seq_len(nrow(cov)))
    cat(sprintf("  %-16s %s %-3d %s\n", cov$trait[i],
                if (cov$n[i] == 0L) "!!" else "  ", cov$n[i], cov$description[i]))
  gaps <- cov$trait[cov$n == 0L]
  if (length(gaps))
    warning("tier covers no problem with trait(s): ", paste(gaps, collapse = ", "),
            "\n  a regression in that area would go unnoticed -- extend BENCH_TIERS",
            call. = FALSE, immediate. = TRUE)
  invisible(cov)
}
