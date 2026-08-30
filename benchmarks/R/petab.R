## =====================================================================
##  petab.R -- turn a PEtab problem into a cppDE benchmark problem.
##
##  A "benchmark problem" is the minimum a solver comparison needs:
##
##      rhs    named character   dx/dt as R-syntax strings
##      parms  named numeric     state initials followed by parameters
##      times  numeric           output grid
##      sens   character         parameters to differentiate w.r.t.
##
##  PEtab carries more than that (observables, noise models, multiple
##  conditions, pre-equilibration).  What is used and what is dropped is
##  recorded in the problem's `notes`, and printed by the runner, so no
##  simplification is silent.
## =====================================================================

## Requires sbml.R to have been sourced first (see bench_source() in
## harness.R, which loads the files in dependency order).
if (!exists("sbml_read")) stop("source benchmarks/R/sbml.R before petab.R")
if (!requireNamespace("yaml", quietly = TRUE))
  stop("package 'yaml' is required to read PEtab problems")


## ---------------------------------------------------------------------
##  PEtab file discovery
## ---------------------------------------------------------------------

## The collection follows `<Dir>/<Dir>.yaml`, with a handful of extra
## yaml files for model variants sitting in sub-directories.
## `variants = TRUE` also returns the alternative formulations that some
## models keep in sub-directories; the default is the one canonical
## problem per publication.
petab_list <- function(root, variants = FALSE) {
  dirs <- list.dirs(root, recursive = FALSE)
  yamls <- unlist(lapply(dirs, function(d) {
    y <- list.files(d, pattern = "\\.yaml$", recursive = variants,
                    full.names = TRUE)
    if (!length(y)) return(NULL)
    ## <Dir>/<Dir>.yaml is the canonical entry point.
    canon <- file.path(d, paste0(basename(d), ".yaml"))
    if (!variants) return(if (file.exists(canon)) canon else y[1L])
    unique(c(if (file.exists(canon)) canon, y))
  }))
  if (!length(yamls)) stop("no PEtab problems found under ", root)
  keep <- vapply(yamls, function(y) {
    txt <- readLines(y, warn = FALSE)
    any(grepl("sbml_files", txt)) && any(grepl("parameter_file", txt))
  }, NA)
  yamls <- yamls[keep]
  ## Name problems after their directory -- that is the publication key
  ## ("Bertozzi_PNAS2020"), whereas the yaml is sometimes "problem.yaml".
  nm <- basename(dirname(yamls))
  sub <- nm[duplicated(nm)]
  if (length(sub))
    nm[nm %in% sub] <- paste0(nm[nm %in% sub], "/",
                              basename(tools::file_path_sans_ext(yamls[nm %in% sub])))
  data.frame(name = nm, yaml = unname(yamls), stringsAsFactors = FALSE)
}

## Counting `<species` in the file costs milliseconds, whereas fully
## parsing Froehlich_CellSystems2018 (1228 states, 2686 reactions) costs
## minutes.  The runner uses this to apply --max-states before paying
## for a model it is going to skip.
petab_species_count <- function(yaml_path) {
  d <- dirname(yaml_path)
  xml <- list.files(d, pattern = "^model_.*\\.xml$", full.names = TRUE)
  if (!length(xml)) return(NA_integer_)
  sum(vapply(readLines(xml[1L], warn = FALSE),
             function(l) lengths(regmatches(l, gregexpr("<species ", l, fixed = TRUE))),
             0L))
}

read_tsv <- function(path) {
  if (!file.exists(path)) return(NULL)
  utils::read.delim(path, sep = "\t", header = TRUE, quote = "",
                    check.names = FALSE, stringsAsFactors = FALSE,
                    colClasses = "character", na.strings = c("", "NA"),
                    comment.char = "")
}

petab_read <- function(yaml_path) {
  y <- yaml::read_yaml(yaml_path)
  dir <- dirname(yaml_path)
  p1 <- y$problems[[1L]]
  pick <- function(x) if (is.null(x)) NULL else file.path(dir, x[[1L]])
  list(
    name        = basename(tools::file_path_sans_ext(yaml_path)),
    dir         = dir,
    sbml_file   = pick(p1$sbml_files),
    parameters  = read_tsv(file.path(dir, y$parameter_file)),
    conditions  = read_tsv(pick(p1$condition_files)),
    measurements= read_tsv(pick(p1$measurement_files)),
    observables = read_tsv(pick(p1$observable_files)),
    simulated   = read_tsv(file.path(dir, paste0("simulatedData_", basename(dir), ".tsv"))))
}


## ---------------------------------------------------------------------
##  Initial-value resolution
## ---------------------------------------------------------------------

## initialAssignment expressions may reference parameters *and* other
## species' initial values.  Substitute until everything is numeric.
resolve_inits <- function(init_exprs, env, max_iter = 25L) {
  states <- names(init_exprs)
  cur <- init_exprs
  for (iter in seq_len(max_iter)) {
    vals <- suppressWarnings(lapply(cur, function(s)
      tryCatch(eval(str2lang(s), env), error = function(e) NULL)))
    done <- !vapply(vals, is.null, NA)
    if (all(done)) return(stats::setNames(unlist(vals), states))
    ## Feed the already-numeric initials back in and retry.
    known <- stats::setNames(
      lapply(num_str_vec(unlist(vals[done])), identity), states[done])
    if (!length(known)) break
    cur[!done] <- subst_str(cur[!done], known)
  }
  bad <- states[vapply(cur, function(s)
    is.null(tryCatch(eval(str2lang(s), env), error = function(e) NULL)), NA)]
  stop("could not evaluate initial value(s) for: ",
       paste(utils::head(bad, 5L), collapse = ", "))
}


## ---------------------------------------------------------------------
##  PEtab problem -> benchmark problem
## ---------------------------------------------------------------------

## Reading and sympifying the SBML dominates the cost for the large
## models, so it is done once per problem and shared by all conditions.
petab_load <- function(yaml_path) {
  pt   <- petab_read(yaml_path)
  sbml <- sbml_read(pt$sbml_file)
  ode  <- sbml_to_ode(sbml)

  par_tab <- pt$parameters
  par_tab$parameterId <- san_id(par_tab$parameterId)
  nominal <- stats::setNames(as.numeric(par_tab$nominalValue), par_tab$parameterId)

  list(pt = pt, sbml = sbml, ode = ode, nominal = nominal,
       estimated = par_tab$parameterId[par_tab$estimate == "1"],
       yaml = yaml_path)
}

## Conditions ranked by how much trajectory they actually carry.  PEtab
## marks steady-state measurements with time = inf; a condition made up
## entirely of those has nothing to integrate.
petab_conditions <- function(loaded) {
  meas <- loaded$pt$measurements
  if (is.null(meas) || !"simulationConditionId" %in% names(meas))
    return(data.frame(conditionId = NA_character_, n_times = 0L,
                      stringsAsFactors = FALSE))
  tt  <- suppressWarnings(as.numeric(meas$time))
  fin <- is.finite(tt)
  ids <- unique(meas$simulationConditionId)
  n_times <- vapply(ids, function(id)
    length(unique(tt[fin & meas$simulationConditionId == id])), 0L)
  df <- data.frame(conditionId = ids, n_times = n_times,
                   stringsAsFactors = FALSE)
  df[order(-df$n_times, df$conditionId), ]
}

petab_problem <- function(yaml_path, condition = NULL, ...)
  petab_case(petab_load(yaml_path), condition = condition, ...)

## All conditions of a problem, most informative first.  `max_conditions`
## caps the fan-out for collections such as Froehlich (9570 conditions).
petab_problems <- function(yaml_path, conditions = NULL,
                           max_conditions = Inf, ...) {
  loaded <- petab_load(yaml_path)
  cond <- petab_conditions(loaded)
  ids <- if (is.null(conditions)) cond$conditionId[1L]
         else if (identical(conditions, "all")) cond$conditionId
         else intersect(conditions, cond$conditionId)
  if (!length(ids)) ids <- NA_character_
  if (length(ids) > max_conditions) ids <- ids[seq_len(max_conditions)]
  out <- lapply(ids, function(id)
    tryCatch(petab_case(loaded, condition = if (is.na(id)) NULL else id, ...),
             error = function(e) NULL))
  Filter(Negate(is.null), out)
}

petab_case <- function(loaded, condition = NULL, max_sens = 64L,
                       min_points = 25L) {
  pt <- loaded$pt; sbml <- loaded$sbml; ode <- loaded$ode
  nominal <- loaded$nominal; estimated <- loaded$estimated
  notes <- character(0)

  vals <- ode$par_values
  hit <- intersect(names(vals), names(nominal))
  vals[hit] <- nominal[hit]

  meas <- pt$measurements
  cond_tab <- pt$conditions
  cond_id <- NULL
  if (!is.null(meas) && "simulationConditionId" %in% names(meas)) {
    cond_id <- if (!is.null(condition)) condition
               else petab_conditions(loaded)$conditionId[1L]
    if (!is.na(cond_id) && !cond_id %in% meas$simulationConditionId)
      stop("condition '", cond_id, "' has no measurements")
  }

  ## Condition columns override model parameters / species initials.
  cond_map <- list()
  if (!is.null(cond_tab) && !is.null(cond_id) && nrow(cond_tab)) {
    row <- cond_tab[cond_tab$conditionId == cond_id, , drop = FALSE]
    if (nrow(row) == 1L) {
      cols <- setdiff(names(cond_tab), c("conditionId", "conditionName"))
      for (cl in cols) {
        v <- row[[cl]][1L]
        if (is.na(v)) next
        target <- san_id(cl)
        num <- suppressWarnings(as.numeric(v))
        if (!is.na(num)) {
          cond_map[[target]] <- num
        } else {
          ## the cell names another parameter -- use its nominal value
          ref <- san_id(v)
          if (!is.na(nominal[ref])) {
            cond_map[[target]] <- unname(nominal[ref])
            if (ref %in% estimated) estimated <- c(estimated, target)
          }
        }
      }
    }
  }
  for (nm in names(cond_map)) vals[nm] <- cond_map[[nm]]

  ## Parameters with no value anywhere default to 0 and are reported.
  if (any(is.na(vals))) {
    notes <- c(notes, sprintf("%d parameter(s) had no value in SBML or PEtab; set to 0",
                              sum(is.na(vals))))
    vals[is.na(vals)] <- 0
  }

  ## -- output time grid --------------------------------------------------
  times <- NULL
  if (!is.null(meas) && "time" %in% names(meas)) {
    tt <- suppressWarnings(as.numeric(
      meas$time[is.null(cond_id) | meas$simulationConditionId == cond_id]))
    tt <- sort(unique(tt[is.finite(tt)]))
    if (length(tt)) times <- tt
  }
  if (is.null(times) || length(times) < 2L) {
    times <- seq(0, 100, length.out = min_points)
    notes <- c(notes, "no usable measurement times; using seq(0, 100)")
  }
  if (times[1L] > 0) times <- c(0, times)
  ## A handful of output points would make the timing mostly call
  ## overhead, so refine sparse grids while keeping the original span.
  if (length(times) < min_points) {
    times <- sort(unique(c(times,
      seq(min(times), max(times), length.out = min_points))))
  }
  tspan <- range(times)

  ## -- resolve piecewise switches over this window -----------------------
  ## Switches that fall inside the window become timed events on
  ## auxiliary states; switches outside it fold to a constant.
  env <- as.list(vals)
  pw <- resolve_piecewise(ode$rhs, env, tspan)
  ode$rhs <- pw$exprs
  events <- pw$events
  if (length(pw$switch_times))
    notes <- c(notes, sprintf(
      "%d discontinuous time switch(es) at t = %s, integrated as timed events",
      length(pw$switch_times),
      paste(format(pw$switch_times, digits = 6), collapse = ", ")))
  if (isTRUE(pw$state_switch))
    notes <- c(notes, "piecewise switches on a state variable -- EXCLUDED")
  ode$rhs <- fold_constants(ode$rhs)

  ## Auxiliary switch states are constant between events.
  aux <- pw$aux_init
  if (length(aux)) {
    ode$rhs  <- c(ode$rhs, stats::setNames(rep("0", length(aux)), names(aux)))
    ode$init <- c(ode$init, unlist(aux))
  }
  ## Initial-value expressions are evaluated at t = 0 only.
  init_pw <- resolve_piecewise(ode$init, env, c(tspan[1L], tspan[1L]), ngrid = 1L)
  ode$init <- init_pw$exprs

  ## -- initial values -----------------------------------------------------
  y0 <- resolve_inits(ode$init, env)
  ## Condition overrides may target a species initial directly.
  for (nm in intersect(names(cond_map), names(y0))) y0[nm] <- cond_map[[nm]]

  ## Species initials are model inputs, not free parameters.
  vals <- vals[setdiff(names(vals), names(y0))]
  parms <- c(y0, vals)

  ## -- sensitivity parameter set --------------------------------------------
  sens_all <- intersect(names(vals), unique(estimated))
  ## Event times are emitted as numeric constants, so the saltation term
  ## that a parameter controlling a switch time would contribute is not
  ## in the model.  Rather than report a silently wrong gradient, drop
  ## those parameters from the sensitivity set.
  switch_pars <- intersect(pw$cond_symbols, sens_all)
  if (length(pw$switch_times) && length(switch_pars)) {
    sens_all <- setdiff(sens_all, switch_pars)
    notes <- c(notes, sprintf(
      "excluded from sensitivities (control an event time): %s",
      paste(switch_pars, collapse = ", ")))
  }
  sens <- sens_all
  if (length(sens) > max_sens) {
    sens <- sens[seq_len(max_sens)]
    notes <- c(notes, sprintf("%d of %d estimated parameters used for sensitivities",
                              max_sens, length(sens_all)))
  }

  model_name <- basename(dirname(loaded$yaml))
  structure(list(
    ## `id` becomes part of a DLL name, so it must be unique per case
    ## and contain nothing a C identifier would reject.
    id      = gsub("[^A-Za-z0-9]", "", paste0(model_name, "_",
                     if (is.null(cond_id) || is.na(cond_id)) "" else cond_id)),
    name    = model_name,
    source  = "petab",
    rhs     = ode$rhs,
    parms   = parms,
    times   = times,
    events  = events,
    sens    = sens,
    sens_available = sens_all,
    fixed   = setdiff(names(parms), sens),
    condition = cond_id,
    nstates = length(ode$rhs),
    npars   = length(vals),
    nsens   = length(sens),
    nevents = if (is.null(events)) 0L else nrow(events),
    usable  = !isTRUE(pw$state_switch),
    notes   = notes,
    petab   = pt,
    sbml    = sbml,
    ode     = ode
  ), class = "bench_problem") -> out
  out$traits <- infer_traits(out)
  out
}

print.bench_problem <- function(x, ...) {
  cat(sprintf("<bench_problem> %s [%s]\n", x$name, x$source))
  cat(sprintf("  states %d | parameters %d | sensitivities %d | out %d | t in [%g, %g]\n",
              x$nstates, x$npars, x$nsens, length(x$times),
              min(x$times), max(x$times)))
  if (!is.null(x$condition)) cat(sprintf("  condition: %s\n", x$condition))
  for (n in x$notes) cat("  note: ", n, "\n", sep = "")
  invisible(x)
}
