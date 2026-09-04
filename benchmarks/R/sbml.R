## =====================================================================
##  sbml.R: an SBML reader for the character-vector ODE form cppODE() consumes.
## =====================================================================

##  Scope is limited to what the Benchmark-Models-PEtab collection uses:
##  compartments, species, parameters, initialAssignment, assignmentRule,
##  functionDefinition, reaction/kineticLaw, and MathML piecewise.

##  algebraicRule, rateRule, stoichiometryMath, delay and events do not occur
##  in the collection and raise an explicit error. Everything is produced as
##  R-syntax strings, the interchange format cppDE uses end to end.

if (!requireNamespace("xml2", quietly = TRUE))
  stop("package 'xml2' is required to read SBML models")


## ---------------------------------------------------------------------
##  Small language-object utilities
## ---------------------------------------------------------------------

## deparse() defaults to 15 significant digits, which silently truncates
## rate constants.  Every deparse in this file goes through here.
deparse_expr <- function(e)
  paste(deparse(e, width.cutoff = 500L,
                control = c("digits17", "keepInteger")), collapse = "")

num_str <- function(x) {
  if (length(x) != 1L || !is.finite(x)) stop("non-finite numeric constant")
  sprintf("%.17g", x)
}

## SBML ids are far more permissive than R names: the collection contains ids
## such as `_Nuc__FOXM1_DBD_Thr600p_`. Every id and every <ci> goes through
## here so that the same SBML symbol always maps to the same R symbol.
.py_keywords <- c(
  "False","None","True","and","as","assert","async","await","break","class",
  "continue","def","del","elif","else","except","finally","for","from","global",
  "if","import","in","is","lambda","nonlocal","not","or","pass","raise","return",
  "try","while","with","yield")

san_id <- function(x) {
  if (!length(x)) return(x)
  y <- make.names(x)
  ## cppDE rejects a Python keyword as a symbol name, SymPy parsing through
  ## Python's own parser. C++ keywords need no folding: the generator
  ## substitutes a slot, so `default` never reaches the source as a name.
  bad <- y %in% .py_keywords
  y[bad] <- paste0(y[bad], "X")
  y[x == "time"] <- "timeX"          # `time` is reserved by cppDE
  y
}

## Recursively replace symbols by language objects.  The function
## position of a call is left alone so that a parameter named e.g. `exp`
## cannot rewrite a call to exp().
subst_lang <- function(e, map) {
  if (is.symbol(e)) {
    nm <- as.character(e)
    if (!is.null(map[[nm]])) return(map[[nm]])
    return(e)
  }
  if (is.call(e)) {
    start <- if (is.symbol(e[[1L]])) 2L else 1L
    if (length(e) >= start)
      for (i in seq.int(start, length(e)))
        if (!is.null(e[[i]])) e[[i]] <- subst_lang(e[[i]], map)
    return(e)
  }
  e
}

## Symbol names occurring in a language object, minus call heads. The large
## models carry expressions with tens of thousands of nodes, so these walk the
## tree directly: deparsing and re-parsing per node is quadratic.
lang_symbols <- function(e) {
  out <- character(0)
  walk <- function(x) {
    if (is.symbol(x)) { out[[length(out) + 1L]] <<- as.character(x); return(invisible()) }
    if (is.call(x)) {
      start <- if (is.symbol(x[[1L]])) 2L else 1L
      if (length(x) >= start)
        for (i in seq.int(start, length(x))) if (!is.null(x[[i]])) walk(x[[i]])
    }
    invisible()
  }
  walk(e)
  unique(out)
}

lang_has_symbol <- function(e) {
  if (is.symbol(e)) return(TRUE)
  if (!is.call(e)) return(FALSE)
  start <- if (is.symbol(e[[1L]])) 2L else 1L
  if (length(e) >= start)
    for (i in seq.int(start, length(e)))
      if (!is.null(e[[i]]) && lang_has_symbol(e[[i]])) return(TRUE)
  FALSE
}

## All symbol names occurring in expression strings.
expr_symbols <- function(txt) {
  txt <- txt[nzchar(txt)]
  if (!length(txt)) return(character(0))
  unique(unlist(lapply(txt, function(s) lang_symbols(str2lang(s)))))
}

## Evaluate every symbol-free subexpression up front. Besides shrinking the
## generated C++ this avoids a codegen corner: CSE types a hoisted constant as
## the AD scalar, which has no cppde::log(double) to build it from.
fold_constants <- function(txt) {
  fold <- function(e) {
    if (!is.call(e)) return(e)
    fname <- if (is.symbol(e[[1L]])) as.character(e[[1L]]) else ""
    if (length(e) >= 2L)
      for (i in seq.int(2L, length(e))) if (!is.null(e[[i]])) e[[i]] <- fold(e[[i]])
    if (fname == ".pw") return(e)          # conditions must stay symbolic
    if (!lang_has_symbol(e)) {
      v <- tryCatch(eval(e, baseenv()), error = function(x) NULL)
      if (is.numeric(v) && length(v) == 1L && is.finite(v))
        return(str2lang(num_str(v)))
    }
    e
  }
  out <- vapply(txt, function(s) deparse_expr(fold(str2lang(s))), "")
  names(out) <- names(txt)
  out
}

## Substitute a named map of expression *strings* into an expression
## string, repeating until a fixed point is reached.  SBML assignment
## rules may reference one another and are not required to be ordered.
subst_str <- function(txt, map, max_iter = 50L) {
  if (!length(map) || !length(txt)) return(txt)
  lang_map <- lapply(map, str2lang)
  keys <- names(map)
  ## Parse once and iterate on the trees: re-parsing per round is what
  ## makes this dominate the run time on the 1228-state models.
  langs <- lapply(txt, str2lang)
  converged <- FALSE
  for (iter in seq_len(max_iter)) {
    hit <- vapply(langs, function(e) any(keys %in% lang_symbols(e)), NA)
    if (!any(hit)) { converged <- TRUE; break }
    langs[hit] <- lapply(langs[hit], subst_lang, map = lang_map)
  }
  if (!converged) stop("assignment rules did not converge, circular reference?")
  out <- vapply(langs, deparse_expr, "")
  names(out) <- names(txt)
  out
}


## ---------------------------------------------------------------------
##  MathML -> R
## ---------------------------------------------------------------------

## A piecewise is emitted as a call to the placeholder `.pw`, taking the value
## and condition pairs and an otherwise branch. It stays in the tree until
## resolve_piecewise() folds it against a concrete interval (see petab.R).
mathml_to_r <- function(node, funcs = list()) {
  nm <- xml2::xml_name(node)

  if (nm == "math" || nm == "lambda") {
    kids <- mathml_children(node)
    kids <- kids[!xml2::xml_name(kids) %in% c("annotation", "bvar", "semantics")]
    if (length(kids) != 1L)
      stop("<", nm, "> must wrap exactly one expression")
    return(mathml_to_r(kids[[1L]], funcs))
  }

  switch(nm,
    apply     = mathml_apply(node, funcs),
    ci        = san_id(trimws(xml2::xml_text(node))),
    cn        = mathml_cn(node),
    piecewise = mathml_piecewise(node, funcs),
    csymbol   = {
      url <- xml2::xml_attr(node, "definitionURL")
      if (!is.na(url) && grepl("symbols/time$", url)) "time"
      else if (!is.na(url) && grepl("symbols/avogadro$", url)) num_str(6.02214179e23)
      else stop("unsupported <csymbol> definitionURL: ", url)
    },
    `true`        = "TRUE",
    `false`       = "FALSE",
    pi            = num_str(pi),
    exponentiale  = num_str(exp(1)),
    infinity      = "Inf",
    notanumber    = "NaN",
    stop("unsupported MathML element <", nm, ">")
  )
}

mathml_children <- function(node) {
  kids <- xml2::xml_children(node)
  if (!length(kids)) return(kids)
  kids[!xml2::xml_name(kids) %in% c("annotation", "annotation-xml")]
}

mathml_cn <- function(node) {
  type <- xml2::xml_attr(node, "type")
  if (is.na(type)) type <- "real"
  parts <- strsplit(trimws(xml2::xml_text(node)), "[[:space:]]+")[[1L]]
  parts <- parts[nzchar(parts)]
  switch(type,
    "e-notation" = {
      if (length(parts) != 2L) stop("malformed e-notation <cn>")
      num_str(as.numeric(parts[1L]) * 10^as.numeric(parts[2L]))
    },
    "rational" = {
      if (length(parts) != 2L) stop("malformed rational <cn>")
      num_str(as.numeric(parts[1L]) / as.numeric(parts[2L]))
    },
    {
      if (length(parts) != 1L) stop("malformed <cn>: ", xml2::xml_text(node))
      num_str(as.numeric(parts[1L]))
    })
}

mathml_piecewise <- function(node, funcs) {
  pieces <- mathml_children(node)
  vals <- character(0); conds <- character(0); otherwise <- NULL
  for (i in seq_along(pieces)) {
    p  <- pieces[[i]]
    pn <- xml2::xml_name(p)
    kid <- mathml_children(p)
    if (pn == "piece") {
      if (length(kid) != 2L) stop("<piece> needs value and condition")
      vals  <- c(vals,  mathml_to_r(kid[[1L]], funcs))
      conds <- c(conds, mathml_to_r(kid[[2L]], funcs))
    } else if (pn == "otherwise") {
      if (length(kid) != 1L) stop("<otherwise> needs exactly one value")
      otherwise <- mathml_to_r(kid[[1L]], funcs)
    } else stop("unexpected <", pn, "> inside <piecewise>")
  }
  if (is.null(otherwise)) otherwise <- "NaN"   # SBML: undefined outside pieces
  args <- character(0)
  for (i in seq_along(vals)) args <- c(args, vals[i], conds[i])
  sprintf(".pw(%s)", paste(c(args, otherwise), collapse = ", "))
}

mathml_apply <- function(node, funcs) {
  kids <- mathml_children(node)
  if (!length(kids)) stop("empty <apply>")
  ## Some models wrap a sub-expression in a bare <apply> purely as
  ## parentheses (Beer_MolBioSystems2014).  Unwrap it.
  if (length(kids) == 1L && xml2::xml_name(kids[[1L]]) == "apply")
    return(mathml_to_r(kids[[1L]], funcs))
  op     <- kids[[1L]]
  opname <- xml2::xml_name(op)
  rest   <- if (length(kids) > 1L) kids[-1L] else kids[0]

  ## <degree> / <logbase> are qualifiers, not operands.
  qual_names <- vapply(seq_along(rest), function(i) xml2::xml_name(rest[[i]]), "")
  degree  <- if ("degree"  %in% qual_names)
    mathml_to_r(mathml_children(rest[[which(qual_names == "degree")[1L]]])[[1L]], funcs)
  logbase <- if ("logbase" %in% qual_names)
    mathml_to_r(mathml_children(rest[[which(qual_names == "logbase")[1L]]])[[1L]], funcs)
  operands <- rest[!qual_names %in% c("degree", "logbase")]

  a <- vapply(seq_along(operands),
              function(i) mathml_to_r(operands[[i]], funcs), "")
  p <- function(x) paste0("(", x, ")")
  infix  <- function(sep) p(paste(p(a), collapse = sep))
  unary1 <- function(fn) { stopifnot(length(a) == 1L); paste0(fn, p(a)) }

  ## A <ci> in operator position calls a functionDefinition.
  if (opname == "ci") {
    fname <- san_id(trimws(xml2::xml_text(op)))
    fd <- funcs[[fname]]
    if (is.null(fd)) stop("call to undefined function '", fname, "'")
    if (length(fd$args) != length(a))
      stop("function '", fname, "' called with ", length(a),
           " argument(s), expected ", length(fd$args))
    map <- stats::setNames(lapply(a, function(s) str2lang(p(s))), fd$args)
    return(p(deparse_expr(subst_lang(fd$body, map))))
  }

  switch(opname,
    plus   = if (!length(a)) "0" else infix(" + "),
    times  = if (!length(a)) "1" else infix(" * "),
    minus  = if (length(a) == 1L) p(paste0("-", p(a))) else infix(" - "),
    divide = { stopifnot(length(a) == 2L); p(paste0(p(a[1L]), " / ", p(a[2L]))) },
    power  = { stopifnot(length(a) == 2L); p(paste0(p(a[1L]), "^", p(a[2L]))) },
    root   = {
      d <- if (is.null(degree)) "2" else degree
      p(paste0(p(a[1L]), "^(1 / ", p(d), ")"))
    },
    log = {
      if (is.null(logbase)) paste0("log10", p(a[1L]))
      else p(paste0("log", p(a[1L]), " / log", p(logbase)))
    },
    abs = unary1("abs"), exp = unary1("exp"), ln = unary1("log"),
    floor = unary1("floor"), ceiling = unary1("ceiling"),
    sin = unary1("sin"), cos = unary1("cos"), tan = unary1("tan"),
    arcsin = unary1("asin"), arccos = unary1("acos"), arctan = unary1("atan"),
    sinh = unary1("sinh"), cosh = unary1("cosh"), tanh = unary1("tanh"),
    eq  = infix(" == "), neq = infix(" != "),
    gt  = infix(" > "),  lt  = infix(" < "),
    geq = infix(" >= "), leq = infix(" <= "),
    and = infix(" & "),  or  = infix(" | "),
    not = unary1("!"),
    stop("unsupported MathML operator <", opname, ">")
  )
}


## ---------------------------------------------------------------------
##  SBML document -> structured list
## ---------------------------------------------------------------------

sbml_read <- function(path) {
  doc <- xml2::read_xml(path)
  xml2::xml_ns_strip(doc)
  model <- xml2::xml_find_first(doc, ".//model")
  if (inherits(model, "xml_missing")) stop("no <model> element in ", path)

  find <- function(xpath) xml2::xml_find_all(model, xpath)
  attr_chr <- function(nodes, a) xml2::xml_attr(nodes, a)
  attr_num <- function(nodes, a) suppressWarnings(as.numeric(xml2::xml_attr(nodes, a)))
  attr_lgl <- function(nodes, a, default = FALSE) {
    v <- xml2::xml_attr(nodes, a)
    ifelse(is.na(v), default, v == "true")
  }

  ## -- unsupported constructs: fail loudly ---------------------------
  for (bad in c("algebraicRule", "stoichiometryMath")) {
    if (length(find(paste0(".//", bad))))
      stop(basename(path), ": <", bad, "> is not supported")
  }

  ## -- function definitions ------------------------------------------
  funcs <- list()
  for (fd in find(".//functionDefinition")) {
    id  <- san_id(xml2::xml_attr(fd, "id"))
    lam <- xml2::xml_find_first(fd, ".//lambda")
    bvars <- san_id(vapply(xml2::xml_find_all(lam, "./bvar/ci"),
                           function(n) trimws(xml2::xml_text(n)), ""))
    ## Bodies may call previously declared functions.
    funcs[[id]] <- list(args = bvars,
                        body = str2lang(mathml_to_r(lam, funcs)))
  }

  ## -- compartments ---------------------------------------------------
  cmp <- find(".//compartment")
  compartments <- data.frame(
    id       = san_id(attr_chr(cmp, "id")),
    size     = ifelse(is.na(attr_num(cmp, "size")), 1, attr_num(cmp, "size")),
    constant = attr_lgl(cmp, "constant", TRUE),
    stringsAsFactors = FALSE)

  ## -- species --------------------------------------------------------
  sp <- find(".//species")
  conc <- attr_num(sp, "initialConcentration")
  amt  <- attr_num(sp, "initialAmount")
  species <- data.frame(
    id          = san_id(attr_chr(sp, "id")),
    compartment = san_id(attr_chr(sp, "compartment")),
    init        = ifelse(is.na(conc), ifelse(is.na(amt), 0, amt), conc),
    is_amount   = attr_lgl(sp, "hasOnlySubstanceUnits", FALSE),
    boundary    = attr_lgl(sp, "boundaryCondition", FALSE),
    constant    = attr_lgl(sp, "constant", FALSE),
    stringsAsFactors = FALSE)
  ## An initialAmount on a concentration species must be divided by the
  ## compartment size to become the value the symbol denotes.
  needs_scale <- !is.na(amt) & is.na(conc) & !species$is_amount
  if (any(needs_scale)) {
    vol <- compartments$size[match(species$compartment[needs_scale], compartments$id)]
    species$init[needs_scale] <- species$init[needs_scale] / vol
  }

  ## -- global parameters ----------------------------------------------
  pr <- find(".//listOfParameters/parameter")
  parameters <- data.frame(
    id    = san_id(attr_chr(pr, "id")),
    value = attr_num(pr, "value"),
    stringsAsFactors = FALSE)

  ## -- initial assignments ---------------------------------------------
  init_assign <- list()
  for (ia in find(".//initialAssignment")) {
    sym <- san_id(xml2::xml_attr(ia, "symbol"))
    init_assign[[sym]] <- mathml_to_r(xml2::xml_find_first(ia, ".//math"), funcs)
  }

  ## -- assignment rules -------------------------------------------------
  assign_rules <- list()
  for (ar in find(".//assignmentRule")) {
    v <- san_id(xml2::xml_attr(ar, "variable"))
    assign_rules[[v]] <- mathml_to_r(xml2::xml_find_first(ar, ".//math"), funcs)
  }

  ## -- rate rules, dx/dt given directly rather than via reactions -----
  rate_rules <- list()
  for (rr in find(".//rateRule")) {
    v <- san_id(xml2::xml_attr(rr, "variable"))
    rate_rules[[v]] <- mathml_to_r(xml2::xml_find_first(rr, ".//math"), funcs)
  }

  ## -- reactions ---------------------------------------------------------
  reactions <- list()
  for (rx in find(".//reaction")) {
    refs <- function(which) {
      nodes <- xml2::xml_find_all(rx, paste0("./listOf", which, "/speciesReference"))
      if (!length(nodes)) return(stats::setNames(numeric(0), character(0)))
      st <- suppressWarnings(as.numeric(xml2::xml_attr(nodes, "stoichiometry")))
      st[is.na(st)] <- 1
      ids <- san_id(xml2::xml_attr(nodes, "species"))
      ## the same species may appear twice in one list
      tapply(st, ids, sum)
    }
    kl <- xml2::xml_find_first(rx, "./kineticLaw/math")
    if (inherits(kl, "xml_missing"))
      stop("reaction ", xml2::xml_attr(rx, "id"), " has no kinetic law")
    ## Local parameters inside a kineticLaw are inlined by value.
    local_pars <- xml2::xml_find_all(rx, "./kineticLaw//listOfParameters/parameter")
    law <- mathml_to_r(kl, funcs)
    if (length(local_pars)) {
      lp <- stats::setNames(
        as.list(num_str_vec(suppressWarnings(as.numeric(xml2::xml_attr(local_pars, "value"))))),
        san_id(xml2::xml_attr(local_pars, "id")))
      law <- subst_str(law, lp)
    }
    reactions[[length(reactions) + 1L]] <- list(
      id        = xml2::xml_attr(rx, "id"),
      reactants = refs("Reactants"),
      products  = refs("Products"),
      law       = law)
  }

  list(id = xml2::xml_attr(model, "id"), file = path,
       compartments = compartments, species = species, parameters = parameters,
       init_assign = init_assign, assign_rules = assign_rules,
       rate_rules = rate_rules, reactions = reactions, functions = funcs)
}

num_str_vec <- function(x) vapply(x, num_str, "")


## ---------------------------------------------------------------------
##  SBML -> ODE system
## ---------------------------------------------------------------------

## Returns rhs (dx/dt, named by state symbol), init (initial-value expressions
## per state), par_values (default value of every free parameter) and assign
## (resolved assignment rules, for observables).
sbml_to_ode <- function(sbml) {
  species      <- sbml$species
  compartments <- sbml$compartments

  ## Species fixed by an assignment rule are not states.  Constant and
  ## boundary species are not driven by reactions, but a rateRule makes
  ## one a state again.
  ruled  <- names(sbml$assign_rules)
  rated  <- names(sbml$rate_rules)
  is_state <- (!species$constant & !species$boundary & !(species$id %in% ruled)) |
              (species$id %in% rated)
  states <- species$id[is_state]
  ## A rateRule may also target a plain parameter, which then becomes a state.
  extra_rated <- setdiff(rated, species$id)
  states <- c(states, extra_rated)
  if (!length(states)) stop("model has no dynamic species")

  ## -- stoichiometric assembly ----------------------------------------
  terms <- stats::setNames(vector("list", length(states)), states)
  for (rx in sbml$reactions) {
    involved <- union(names(rx$reactants), names(rx$products))
    for (s in intersect(involved, states)) {
      nu <- (if (is.na(rx$products[s]))  0 else rx$products[s]) -
            (if (is.na(rx$reactants[s])) 0 else rx$reactants[s])
      if (nu == 0) next
      terms[[s]] <- c(terms[[s]],
                      if (nu == 1) paste0("(", rx$law, ")")
                      else sprintf("%s * (%s)", num_str(nu), rx$law))
    }
  }

  rhs <- vapply(states, function(s) {
    if (!is.null(sbml$rate_rules[[s]])) return(sbml$rate_rules[[s]])
    tm <- terms[[s]]
    if (!length(tm)) return("0")
    body <- paste(tm, collapse = " + ")
    row <- species[match(s, species$id), ]
    ## Kinetic laws are extents (amount/time).  A species in
    ## concentration units therefore needs division by its volume.
    if (row$is_amount) body else sprintf("(%s) / (%s)", body, row$compartment)
  }, "")

  ## -- fold assignment rules into the RHS ------------------------------
  rules <- unlist(sbml$assign_rules)
  if (length(rules)) {
    rules <- subst_str(rules, as.list(rules))   # rules may nest
    rhs   <- subst_str(rhs, as.list(rules))
  }

  ## -- free parameters --------------------------------------------------
  used <- expr_symbols(c(rhs, unlist(sbml$init_assign), rules))
  free <- setdiff(used, c(states, "time", "TRUE", "FALSE", "Inf", "NaN", "pi"))

  par_values <- stats::setNames(rep(NA_real_, length(free)), free)
  m <- match(free, sbml$parameters$id)
  par_values[!is.na(m)] <- sbml$parameters$value[m[!is.na(m)]]
  m <- match(free, compartments$id)
  par_values[!is.na(m)] <- compartments$size[m[!is.na(m)]]
  m <- match(free, species$id)                    # constant/boundary species
  par_values[!is.na(m)] <- species$init[m[!is.na(m)]]

  ## -- initial values as expressions -------------------------------------
  init <- stats::setNames(
    vapply(states, function(s) {
      ia <- sbml$init_assign[[s]]
      if (!is.null(ia)) return(ia)
      i <- match(s, species$id)
      if (!is.na(i)) return(num_str(species$init[i]))
      j <- match(s, sbml$parameters$id)        # rateRule on a parameter
      num_str(if (!is.na(j) && !is.na(sbml$parameters$value[j]))
                sbml$parameters$value[j] else 0)
    }, ""), states)

  list(rhs = rhs, init = init, par_values = par_values,
       assign = rules,
       init_assign = sbml$init_assign,
       states = states, id = sbml$id)
}


## ---------------------------------------------------------------------
##  Piecewise resolution
## ---------------------------------------------------------------------

## Every `piecewise` in the collection switches on `time` alone, so over a
## concrete window it folds to a constant, becomes a timed event at the jump,
## or is state-dependent and the caller excludes the problem.

## A timed event is the numerically correct treatment: the solver stops at t*,
## applies the jump and restarts instead of stepping across a discontinuity.

## Events act on state variables, so each switch gets an auxiliary state with
## dx/dt = 0 carrying the switched value, or indicator states where a branch is
## itself a function of time. Identical switches share one auxiliary state.
resolve_piecewise <- function(exprs, env, tspan, ngrid = 4001L,
                              make_events = TRUE, prefix = "pw") {
  switches <- numeric(0); state_switch <- FALSE
  cond_symbols <- character(0)
  aux_init <- list()                     # aux state -> initial value expression
  ev <- list()                           # rows of the event table
  registry <- list()                     # .pw text -> replacement expression
  counter <- 0L
  tgrid <- if (ngrid > 1L) seq(tspan[1L], tspan[2L], length.out = ngrid) else tspan[1L]
  ngrid <- length(tgrid)

  eval_at <- function(e, t) {
    v <- try(eval(e, c(as.list(env), list(time = t))), silent = TRUE)
    if (inherits(v, "try-error")) return(NULL)
    if (length(v) == 1L && length(t) > 1L) v <- rep(v, length(t))
    if (length(v) != length(t)) return(NULL)
    v
  }

  branch_grid <- function(conds) {
    b <- rep(length(conds) + 1L, ngrid)  # the <otherwise> branch
    for (k in rev(seq_along(conds))) {
      cv <- eval_at(conds[[k]], tgrid)
      if (is.null(cv)) { state_switch <<- TRUE; cv <- rep(FALSE, ngrid) }
      b[which(as.logical(cv))] <- k
    }
    b
  }
  branch_at <- function(conds, t) {
    b <- length(conds) + 1L
    for (k in rev(seq_along(conds))) {
      cv <- eval_at(conds[[k]], t)
      if (!is.null(cv) && isTRUE(as.logical(cv))) b <- k
    }
    b
  }
  ## Bisect down to the exact instant the active branch changes.
  refine <- function(conds, lo, hi) {
    blo <- branch_at(conds, lo)
    for (i in seq_len(60L)) {
      mid <- (lo + hi) / 2
      if (!(mid > lo && mid < hi)) break
      if (branch_at(conds, mid) == blo) lo <- mid else hi <- mid
    }
    hi
  }

  walk <- function(e) {
    if (!is.call(e)) return(e)
    if (is.symbol(e[[1L]]) && identical(as.character(e[[1L]]), ".pw")) {
      args <- lapply(as.list(e)[-1L], walk)
      n <- length(args)
      if (n %% 2L != 1L) stop("malformed .pw() call")
      vals  <- args[c(seq(1L, n - 1L, by = 2L), n)]      # branch values, last = otherwise
      conds <- args[seq(2L, n - 1L, by = 2L)]
      key <- deparse_expr(as.call(c(list(as.symbol(".pw")), args)))
      if (!is.null(registry[[key]])) return(registry[[key]])

      for (cd in conds)
        cond_symbols <<- c(cond_symbols, setdiff(expr_symbols(deparse_expr(cd)), "time"))

      b <- branch_grid(conds)
      chg <- which(diff(b) != 0L)
      if (!length(chg)) return(vals[[b[1L]]])           # constant: fold

      ## -- segment boundaries and the branch active on each segment,
      tsw <- vapply(chg, function(i) refine(conds, tgrid[i], tgrid[i + 1L]), 0)
      seg_branch <- c(b[1L], b[chg + 1L])
      switches <<- c(switches, tsw)
      if (!make_events) return(vals[[b[1L]]])

      counter <<- counter + 1L
      time_dep <- vapply(vals[unique(seg_branch)],
                         function(v) "time" %in% expr_symbols(deparse_expr(v)), NA)

      repl <- if (!any(time_dep)) {
        ## One auxiliary state carrying the switched value.
        nm <- sprintf("%ssw%d", prefix, counter)
        aux_init[[nm]] <<- deparse_expr(vals[[seg_branch[1L]]])
        for (j in seq_along(tsw))
          ev[[length(ev) + 1L]] <<- data.frame(
            var = nm, time = num_str(tsw[j]),
            value = deparse_expr(vals[[seg_branch[j + 1L]]]),
            method = "replace", root = NA_character_, stringsAsFactors = FALSE)
        str2lang(nm)
      } else {
        ## Indicator states selecting between time-dependent branches.
        used <- unique(seg_branch)
        nms <- sprintf("%sind%d_%d", prefix, counter, seq_along(used))
        for (u in seq_along(used)) {
          aux_init[[nms[u]]] <<- if (used[u] == seg_branch[1L]) "1" else "0"
          for (j in seq_along(tsw))
            ev[[length(ev) + 1L]] <<- data.frame(
              var = nms[u], time = num_str(tsw[j]),
              value = if (used[u] == seg_branch[j + 1L]) "1" else "0",
              method = "replace", root = NA_character_, stringsAsFactors = FALSE)
        }
        str2lang(paste(sprintf("%s * (%s)", nms,
                               vapply(vals[used], deparse_expr, "")),
                       collapse = " + "))
      }
      registry[[key]] <<- repl
      return(repl)
    }
    start <- if (is.symbol(e[[1L]])) 2L else 1L
    if (length(e) >= start)
      for (i in seq.int(start, length(e)))
        if (!is.null(e[[i]])) e[[i]] <- walk(e[[i]])
    e
  }

  out <- vapply(exprs, function(s) deparse_expr(walk(str2lang(s))), "")
  names(out) <- names(exprs)
  list(exprs        = out,
       switch_times = sort(unique(switches)),
       state_switch = state_switch,
       aux_init     = aux_init,
       events       = if (length(ev)) do.call(rbind, ev) else NULL,
       cond_symbols = unique(cond_symbols))
}
