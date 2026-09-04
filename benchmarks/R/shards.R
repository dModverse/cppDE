## =====================================================================
##  shards.R: the problem set a *remote* run works on.
## =====================================================================

##  Both submitters parse the problems locally and ship them in the transferred
##  workspace. They have to select and split identically, or their results are
##  comparable neither with each other nor with a local run-benchmarks.R.

## The tolerance sets every driver sweeps.  A tier names one of these;
## an explicit --tol overrides it.
BENCH_TOLSETS <- list(
  loose       = data.frame(atol = 1e-8,  rtol = 1e-6),
  tight       = data.frame(atol = 1e-10, rtol = 1e-8),
  loose_tight = data.frame(atol = c(1e-8, 1e-11), rtol = c(1e-6, 1e-9)),
  default     = data.frame(atol = c(1e-6, 1e-9, 1e-12),
                           rtol = c(1e-4, 1e-7, 1e-10)),
  wp          = data.frame(atol = 10^-c(5, 7, 9, 11, 13),
                           rtol = 10^-c(3, 5, 7, 9, 11)))


## ---------------------------------------------------------------------
##  Building the problem list
## ---------------------------------------------------------------------

## Parses the tier's classic and PEtab problems, then drops what the runner
## would drop anyway, flagged unusable or over the state cap, so that a shard
## never carries a case that dies on arrival.

## `conditions` is the string the drivers take on the command line: "all", or a
## count per PEtab model.
bench_problems_for_tier <- function(tier, petab_root, conditions = "1",
                                    max_states = 400L, max_sens = 32L,
                                    cores = max(1L, parallel::detectCores() - 2L),
                                    verbose = TRUE) {
  TIER <- BENCH_TIERS[[tier]]
  if (is.null(TIER)) stop("unknown tier '", tier, "'")
  all_cond <- identical(conditions, "all")

  problems <- list()
  cl <- classic_problems()
  if (!is.null(TIER$classic)) cl <- cl[intersect(TIER$classic, names(cl))]
  problems <- c(problems, lapply(cl, list))

  if (dir.exists(petab_root)) {
    idx <- petab_list(petab_root)
    if (!is.null(TIER$petab)) idx <- idx[idx$name %in% TIER$petab, , drop = FALSE]
    ## Filter on size before parsing: the largest models in the
    ## collection take minutes to translate and would then be skipped.
    nsp <- vapply(idx$yaml, petab_species_count, 0L)
    idx <- idx[is.na(nsp) | nsp <= max_states, , drop = FALSE]
    built <- parallel::mclapply(seq_len(nrow(idx)), function(k) tryCatch(
      petab_problems(idx$yaml[k],
                     conditions = if (all_cond) "all" else NULL,
                     max_conditions = if (all_cond) Inf else as.integer(conditions),
                     max_sens = max_sens),
      error = function(e) list()),
      mc.cores = cores)
    for (k in seq_along(built))
      if (length(built[[k]])) problems[[idx$name[k]]] <- built[[k]]
  } else warning("PEtab collection not found at ", petab_root, call. = FALSE)

  usable <- vapply(problems, function(cs) isTRUE(cs[[1L]]$usable) &&
                     cs[[1L]]$nstates <= max_states, NA)
  if (verbose && any(!usable))
    cat(sprintf("  excluded: %s\n", paste(names(problems)[!usable], collapse = ", ")))
  problems[usable]
}


## ---------------------------------------------------------------------
##  Splitting it
## ---------------------------------------------------------------------

## Balance the shards by expected cost rather than by count: the in-process
## --shard option has to decide before the models are parsed and can only
## round-robin over names, whereas here their size is known.

## Longest-processing-time first, and deterministic because the sort is, so the
## same problem set always splits the same way and two runs of different width
## stay comparable shard by shard.

## The cost proxy is states x conditions x sensitivities. A model marked for the
## sparse sweep carries three times the solver configs and is weighted for it,
## or the shard holding the sparse models finishes long after the others.
balance_shards <- function(problems, n) {
  cost <- vapply(problems, function(cs)
    cs[[1L]]$nstates * length(cs) * max(1L, cs[[1L]]$nsens) *
      (if (isTRUE(cs[[1L]]$sweep)) 3 else 1), 0)
  ord <- order(-cost, names(problems))
  load <- numeric(n)
  shard_of <- integer(length(problems))
  for (k in ord) {
    pick <- which.min(load)
    shard_of[k] <- pick
    load[pick] <- load[pick] + cost[k]
  }
  list(shards = lapply(seq_len(n), function(k) problems[shard_of == k]),
       load = load)
}

## The per-shard line both submitters print before they send anything.
print_shard_plan <- function(shards, load) {
  for (k in seq_along(shards))
    cat(sprintf("  shard %2d: %2d problems, %6.0f cost, %5d states  (%s)\n", k,
                length(shards[[k]]), load[k],
                sum(vapply(shards[[k]], function(cs) cs[[1L]]$nstates, 0L)),
                paste(utils::head(names(shards[[k]]), 3L), collapse = ", ")))
  cat(sprintf("  balance: heaviest / lightest shard = %.2f\n",
              max(load) / max(min(load), 1)))
}
