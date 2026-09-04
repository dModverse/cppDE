#' Generate and Compile a C++ ODE Solver
#'
#' Generates C++ source code for an ODE system \eqn{\dot{x} = f(x, p)},
#' compiles it via `R CMD SHLIB`, and returns a handle for use with
#' [solveODE()]. The four available methods (`"bdf"`, `"adams"`,
#' `"rb4"`, `"tsit5"`) and their selection criteria, the dual-number
#' sensitivity framework, and the event handling are described in
#' `vignette("Methods", package = "cppDE")`.
#'
#' Events are specified as a `data.frame` with columns `var` (affected
#' variable), `value` (expression to apply), `method` (`"replace"`,
#' `"add"`, or `"multiply"`), and exactly one of `time` (event time) or
#' `root` (root expression). `rootfunc` terminates integration:
#' `"equilibrate"` stops at steady state; a character vector of
#' expressions stops at the first zero crossing.
#'
#' @param rhs Named character vector of ODE right-hand sides. Names are
#'   the state variables.
#' @param events Optional event `data.frame`. See Details.
#' @param rootfunc Optional integration-termination root: `"equilibrate"`
#'   or a character vector of expressions.
#' @param fixed Optional character vector of state or parameter names
#'   excluded from sensitivities at compile time.
#' @param forcings Optional character vector of forcing-function names
#'   referenced by `rhs`.
#' @param compile Logical. Compile and load the generated C++ code.
#' @param modelname Optional base name for the generated source file
#'   and exported C symbols. A random identifier is used when `NULL`.
#' @param outdir Directory for the generated C++ source. Default
#'   `tempdir()`.
#' @param deriv Logical. Compute first-order parameter sensitivities.
#' @param deriv2 Logical. Compute second-order parameter sensitivities;
#'   implies `deriv = TRUE`.
#' @param nStack Compile-time AD slab width. `Inf` (default) selects
#'   heap-allocated AD; the per-call sensitivity dimension is taken at
#'   run time from `ncol(sens1ini)` in [solveODE()]. A positive integer
#'   `K` fixes the width to `K` at compile time. `NULL` selects
#'   `length(c(variables, parameters)) - length(fixed)`. `deriv2 = TRUE`
#'   requires a finite width and silently demotes `Inf` to `NULL`.
#' @param includeTimeZero Logical. Ensure that `0` is part of the
#'   integration times.
#' @param useDenseOutput Logical. Use Hermite dense output for
#'   user-requested time points. Applies to `method = "rb4"` and
#'   `"tsit5"`, whose interpolant is of lower order than the method;
#'   `FALSE` makes the solver land on each requested time exactly.
#'   Ignored with a warning for `"bdf"` and `"adams"`.
#' @param sparse Logical or `NULL`. `NULL` auto-selects sparse vs.
#'   dense LU from the Jacobian sparsity; `TRUE` forces sparse, `FALSE`
#'   forces dense. Sparse LU requires KLU at install time.
#' @param method Integration method: `"bdf"` (default), `"adams"`,
#'   `"rb4"`, or `"tsit5"`.
#' @param useNDF Logical. Use Klopfenstein-Shampine NDF coefficients in
#'   the BDF corrector. Applies to `method = "bdf"`; ignored otherwise.
#' @param profile Logical. Compile with profiling counters.
#' @param stepTrace Logical. Compile to record per-step diagnostics,
#'   returned as `$trace` from [solveODE()].
#' @param verbose Logical. Print progress messages.
#'
#' @return The compiled model name (character) carrying the attributes
#'   required by [solveODE()]: `equations`, `srcfile`, `variables`,
#'   `parameters`, `forcings`, `events`, `rootfunc`, `fixed`, `jacobian`
#'   (with components `f.x` and `f.time`), `deriv`, `deriv2`, `nStack`,
#'   `sparse`, `method`, `useNDF`, `dimNames`, `compileArgs`, `backend`.
#'
#' @example inst/examples/example_ODE.R
#' @importFrom stats setNames
#' @seealso [solveODE()] for integration; [cvode()] for the
#'   SUNDIALS-backed alternative; [funCpp()] for algebraic functions;
#'   `vignette("Methods", package = "cppDE")` for behaviour.
#' @export
cppODE <- function(rhs, events = NULL, rootfunc = NULL, fixed = NULL, forcings = NULL,
                   compile = TRUE, modelname = NULL, outdir = tempdir(),
                   deriv = TRUE, deriv2 = FALSE,
                   nStack = Inf,
                   includeTimeZero = TRUE, useDenseOutput = TRUE,
                   sparse = NULL,
                   method = c("bdf", "adams", "rb4", "tsit5"),
                   useNDF = TRUE,
                   profile = FALSE, stepTrace = FALSE, verbose = FALSE) {

  # --- Validate arguments ---
  if (deriv2 && !deriv) {
    warning("deriv2 = TRUE requires deriv = TRUE. Setting deriv = TRUE automatically.")
    deriv <- TRUE
  }
  method <- match.arg(method)
  if (method == "rb4") method <- "rosenbrock4"
  # Both multistep methods ("bdf", "adams") instantiate cppde::multistepper,
  # selected at compile time by the multistep_method enum; the single-step ones
  # use onestep_controller / onestep_dense_output.
  is_multistep <- function(m) m %in% c("bdf", "adams")
  is_explicit  <- function(m) m %in% c("tsit5")

  if (!useDenseOutput && is_multistep(method)) {
    warning("'useDenseOutput = FALSE' is ignored for method = \"", method,
            "\"; multistep methods always use Nordsieck dense output",
            call. = FALSE)
    useDenseOutput <- TRUE
  }

  # --- Clean up ODE definitions ---
  rhs <- unclass(rhs)
  rhs <- gsub("\n", "", rhs)

  # Every expression the model is built from. Numeric event columns hold no
  # identifiers, and "equilibrate" is a keyword of rootfunc, not an expression.
  checkSymbolNames(rhs, forcings, fixed)
  if (!is.null(events)) {
    for (col in c("var", "value", "time", "root"))
      if (col %in% names(events) && is.character(events[[col]]))
        checkSymbolNames(events[[col]])
  }
  if (!is.null(rootfunc) && !identical(tolower(rootfunc), "equilibrate"))
    checkSymbolNames(rootfunc)

  # --- Extract variable and parameter names ---
  variables <- names(rhs)

  # Collect all expressions for symbol extraction
  all_expressions <- rhs
  if (!is.null(events)) {
    bad <- which(!xor(!is.na(events$time), !is.na(events$root)))
    if (length(bad) > 0) {
      stop(
        sprintf(
          "Each event must define exactly one of 'time' or 'root'. Invalid event(s): %s",
          paste(bad, collapse = ", ")
        )
      )
    }
    all_expressions <- c(all_expressions,
                         events$value,
                         if ("time" %in% names(events)) events$time,
                         if ("root" %in% names(events)) events$root)
  }
  # Include rootfunc expressions (but not "equilibrate" which is a keyword)
  if (!is.null(rootfunc) && !identical(tolower(rootfunc), "equilibrate")) {
    all_expressions <- c(all_expressions, rootfunc)
  }

  symbols <- getSymbols(all_expressions)

  # --- Validate forcings ---
  if (is.null(forcings)) forcings <- character(0)

  # Forcings must be symbols in rhs but NOT state names
  if (length(forcings) > 0) {
    unknown_forcings <- setdiff(forcings, symbols)
    if (length(unknown_forcings) > 0) {
      stop("Unknown forcing symbols: ", paste(unknown_forcings, collapse = ", "))
    }
    forcing_states <- intersect(forcings, variables)
    if (length(forcing_states) > 0) {
      stop("Forcing names cannot be state variables: ", paste(forcing_states, collapse = ", "))
    }
  }
  n_forcings <- length(forcings)

  # Parameters are symbols that are not variables, forcings, or time

  params <- setdiff(symbols, c(variables, forcings, "time"))

  # --- Handle fixed initial conditions and parameters ---
  if (is.null(fixed)) fixed <- character(0)
  fixed_initials <- if (deriv) intersect(fixed, variables) else character(0)
  fixed_params <- if (deriv) intersect(fixed, params) else character(0)
  sens_initials  <- if (deriv) setdiff(variables, fixed_initials) else character(0)
  sens_params  <- if (deriv) setdiff(params, fixed_params) else character(0)

  # Index maps
  variable_idx0 <- setNames(seq_along(variables) - 1L, variables)
  param_idx0 <- setNames(seq_along(params) - 1L, params)
  fixed_initial_idx  <- variable_idx0[fixed_initials]
  fixed_param_idx  <- param_idx0[fixed_params]

  # --- Calculate dimensions ---
  n_variables <- length(variables)
  n_params <- length(params)
  n_sens_initials <- length(sens_initials)
  n_sens_params <- length(sens_params)
  n_total_sens <- n_sens_initials + n_sens_params

  # --- Resolve nStack (compile-time AD slab width) ---
  # Inf heap-allocates and takes the width from ncol(sens1ini) at run time,
  # NULL stacks it at n_total_sens, K fixes it and every call has to fit K.
  if (!deriv && is.numeric(nStack) && length(nStack) == 1L && is.infinite(nStack)) {
    nStack <- NULL  # heap AD is meaningful only with deriv = TRUE
  }
  if (is.null(nStack)) {
    nStack_width <- as.integer(n_total_sens)
    is_heap <- FALSE
  } else if (is.numeric(nStack) && length(nStack) == 1L && is.infinite(nStack) && nStack > 0) {
    if (!deriv) stop("'nStack = Inf' requires deriv = TRUE")
    nStack_width <- 0L  # routes codegen to <double, 0> (heap spec)
    is_heap <- TRUE
  } else if (is.numeric(nStack) && length(nStack) == 1L && is.finite(nStack) && nStack >= 0 &&
             nStack == as.integer(nStack)) {
    if (!deriv) stop("'nStack' is only meaningful when deriv = TRUE")
    nStack_width <- as.integer(nStack)
    is_heap <- FALSE
  } else {
    stop("'nStack' must be NULL, a non-negative integer, or Inf")
  }
  # Codegen helper: under heap AD, every diff() seeding call must pass the
  # runtime size as a second arg so the tangent slab is allocated. Stack AD
  # uses the static-N spec where diff(idx) takes no size arg.
  dyn_arg <- if (is_heap) ", n_sens" else ""

  # --- Generate unique model name ---
  if (is.null(modelname)) {
    modelname <- paste(c("x", sample(c(letters, 0:9), 8, TRUE)), collapse = "")
  }
  modelname <- unique_modelname(modelname)

  # Lazy import
  codegen <- get_codegen_cppODE_py()

  if (verbose) message("Generating ODE and Jacobian code...")

  # Single Python call generates everything
  if (deriv2) {
    numType <- "AD2"  # F<F<double>>
  } else if (deriv) {
    numType <- "AD"   # F<double>
  } else {
    numType <- "double"
  }

  # Auto-selected sparse without KLU -> dense. An explicit sparse = TRUE
  # is left alone and reaches compile().
  sparse_for_codegen <- sparse
  if (is.null(sparse) && !isTRUE(cvodeConfig$klu_available)) {
    sparse_for_codegen <- FALSE
  }

  codegen_result <- codegen$generate_ode_cpp(
    rhs_dict = as.list(setNames(rhs, variables)),
    params_list = params,
    num_type = numType,
    fixed_states = fixed_initials,
    fixed_params = fixed_params,
    forcings_list = forcings,
    sparse = sparse_for_codegen,
    skip_jacobian = is_explicit(method)
  )

  ode_code <- codegen_result$ode_code
  jac_code <- codegen_result$jac_code
  time_derivs_str <- codegen_result$time_derivs

  if (verbose) message("  \u2713 ODE and Jacobian generated")

  # --- Sparse LU decision ---
  # use_sparse must match what the codegen generated: the Jacobian functor
  # signature is tied to the stepper's matrix type.
  use_sparse <- isTRUE(codegen_result$use_sparse)

  if (use_sparse) {
    stats <- codegen_result$sparsity_stats
    message(sprintf("Sparse Jacobian detected (%dx%d, %d nnz, %.3f%% sparse)",
                    stats$n, stats$n, stats$jac_nnz, stats$jac_zeros_pct))
  }

  # --- Generate event code if needed ---
  event_code <- ""
  if (!is.null(events)) {
    if (verbose) message("Generating event code...")

    event_lines <- codegen$generate_event_code(
      events_df = events,
      states_list = variables,
      params_list = params,
      n_states = n_variables,
      num_type = numType,
      forcings_list = forcings,
      rhs_dict = as.list(setNames(rhs, variables))
    )

    event_code <- paste(event_lines, collapse = "\n")

    ## Fixed-event times as plain-double expressions over the flat [states, params]
    ## vector, so the batch entry can size its output exactly. NULL means it cannot:
    ## a root event, or a time built from a forcing.
    event_time_exprs <- codegen$fixed_event_time_exprs(
      events_df = events, states_list = variables, params_list = params,
      n_states = n_variables, forcings_list = forcings)
  }

  # --- Generate rootfunc code if needed ---
  rootfunc_code <- ""
  if (!is.null(rootfunc)) {
    if (verbose) message("Generating rootfunc code...")

    rootfunc_lines <- codegen$generate_rootfunc_code(
      rootfunc = rootfunc,
      states_list = variables,
      params_list = params,
      n_states = n_variables,
      num_type = numType,
      forcings_list = forcings
    )

    rootfunc_code <- paste(rootfunc_lines, collapse = "\n")
    if (verbose) message("  \u2713 rootfunc generated")
  }

  # --- Generate forcing initialization code ---
  forcing_init_code <- paste(codegen$generate_forcing_init_code(n_forcings, numType), collapse = "\n")

  # --- C++ includes ---
  includings <- c(
    "#define R_NO_REMAP",
    "#include <R.h>",
    "#include <Rinternals.h>",
    "#include <algorithm>",
    "#include <vector>",
    "#include <cmath>",
    "#include <climits>",
    "#include <cstdio>",
    "#include <string>",
    "#include <cppde/cppde.hpp>",
    "#include <cppde/cppde_r_batch.hpp>"
  )

  # --- Using declarations ---
  # N is nStack_width, so the AD type is cppde::dual<double, N>, arena-allocated
  # at N = 0; deriv2 uses cppde::dual2nd<double, N> with the same layout.
  if (deriv2) {
    usings <- c(
      "using namespace cppde;",
      sprintf("using AD = cppde::dual<double, %d>;", nStack_width),
      sprintf("using AD2 = cppde::dual2nd<double, %d>;", nStack_width)
    )
  } else if (deriv) {
    usings <- c(
      "using namespace cppde;",
      sprintf("using AD = cppde::dual<double, %d>;", nStack_width)
    )
  } else {
    usings <- c(
      "using namespace cppde;"
    )
  }

  # --- Observer ---
  observer_lines <- c(
    "// Observer: stores trajectory values in vectors",
    "struct observer {",
    sprintf("  std::vector<%s>& times;", numType),
    sprintf("  std::vector<%s>& y;", numType),
    "",
    sprintf("  explicit observer(std::vector<%s>& t, std::vector<%s>& y_)", numType, numType),
    "    : times(t), y(y_) {}",
    "",
    sprintf("  void operator()(const cppde::vector_t<%s>& x, const %s& t) {", numType, numType),
    "    times.push_back(t);",
    "    for (size_t i = 0; i < x.size(); ++i) y.push_back(x[i]);",
    "  }",
    "};"
  )
  observer_code <- paste(observer_lines, collapse = "\n")

  # --- Solver function (externC) ---
  externC <- c(
    "// Pure C++ solve: no R API, no escaping exceptions.  Results are read out",
    "// of the AD types into plain double here, while the arena scope is still",
    "// open, tangents live in the arena and die when it pops.",
    "static int solve_impl(const cppde::rbatch::solve_args& args,",
    "                      cppde::rbatch::solve_result& res) noexcept {",
    "try {",
    "",
    "  cppde::dual_arena::scope _cppde_arena_scope;",
    "  cppde::ndf_detail::trace_scope _cppde_trace_scope(res.trace);",
    "",
    "  StepChecker checker(args.maxprogress, args.maxsteps);",
    "",
    sprintf("  cppde::vector_t<%s> x(%d);", numType, n_variables),
    sprintf("  cppde::vector_t<%s> full_params(%d);", numType, n_variables + n_params),
    ""
  )

  # --- Custom sensitivity initial values ---
  if (deriv) {
    externC <- c(
      externC,
      "  // Custom sensitivity initial values (length validated after n_sens is computed)",
      "  bool has_sens1ini = (args.sens1ini != nullptr);",
      "  const double* sens1ini = has_sens1ini ? args.sens1ini : nullptr;"
    )

    if (deriv2) {
      externC <- c(
        externC,
        "  bool has_sens2ini = (args.sens2ini != nullptr);",
        "  const double* sens2ini = has_sens2ini ? args.sens2ini : nullptr;"
      )
    } else {
      externC <- c(
        externC,
        "  if (args.sens2ini != nullptr)",
        "    return res.fail(cppde::RC_ILL_INPUT, \"sens2ini supplied but deriv2 = FALSE\");"
      )
    }

    externC <- c(externC, "")
  }

  # --- Runtime fixed parameters (must be defined before n_sens and active_idx) ---
  if (deriv) {
    externC <- c(
      externC,
      "  // Runtime fixed parameters - O(1) lookup via boolean vector",
      sprintf("  std::vector<bool> is_runtime_fixed(%d, false);  // size = n_sens_total (compile-time)", n_total_sens),
      "  if ((args.n_fixed > 0)) {",
      "    const int* fixed_ptr = args.fixed;",
      "    int n_fixed = args.n_fixed;",
      "    for (int i = 0; i < n_fixed; ++i) {",
      sprintf("      if (fixed_ptr[i] >= 0 && fixed_ptr[i] < %d) {", n_total_sens),
      "        is_runtime_fixed[fixed_ptr[i]] = true;",
      "      }",
      "    }",
      "  }",
      "  int n_runtime_fixed = std::count(is_runtime_fixed.begin(), is_runtime_fixed.end(), true);",
      ""
    )
  }

  # --- Sensitivity dimensions and index helpers (depends on is_runtime_fixed) ---
  if (deriv) {
    externC <- c(
      externC,
      sprintf("  const int n_states     = %d;", n_variables),
      sprintf("  const int n_params_all = %d;", n_params),
      sprintf("  const int n_phi_rows   = %d;  // n_states + n_params (full Phi' row count)", n_variables + n_params),
      sprintf("  const int n_sens_total = %d;  // compile-time total (excl. compile-time fixed)", n_total_sens),
      sprintf("  const int n_stack_max  = %s;  // compile-time AD slab width (INT_MAX under heap AD)",
              if (is_heap) "INT_MAX" else as.character(nStack_width)),
      "  // Per-call active sens dimension: from sens1ini's column count when supplied",
      "  // (which is the auto-extended full Phi' shape on the R side, legacy",
      "  // [n_states, n_active] input is padded with an identity block on the param",
      "  // rows before reaching here), or n_sens_total - n_runtime_fixed when not",
      "  // supplied (identity-fallback seeding).",
      "  const int n_sens = has_sens1ini",
      "      ? args.n_sens1_cols",
      "      : (n_sens_total - n_runtime_fixed);",
      "  if (n_sens > n_stack_max) {",
      "    char _m[256];",
      "    snprintf(_m, sizeof(_m), \"sens1ini has %d columns but the model's compile-time nStack is %d\", n_sens, n_stack_max);",
      "    return res.fail(cppde::RC_ILL_INPUT, _m);",
      "  }",
      "",
      "  // active_idx[i]: compile-time sens index i -> active index in [0, n_sens), or -1 if runtime-fixed",
      "  std::vector<int> active_idx(n_sens_total, -1);",
      "  {",
      "    int k = 0;",
      "    for (int i = 0; i < n_sens_total; ++i) {",
      "      if (!is_runtime_fixed[i]) active_idx[i] = k++;",
      "    }",
      "  }",
      "",
      "  // IDX1 / IDX2 index into sens1ini / sens2ini, which are passed as the full",
      "  // Phi'(theta) / Phi''(theta) with first dim = n_phi_rows (state rows + param rows).",
      "  // R-side coerce auto-extends legacy [n_states, n_active] input with an identity",
      "  // block on param rows, so C++ always sees the full shape here.",
      "  auto IDX1 = [n_phi_rows](int g, int v) {",
      "    return g + n_phi_rows * v;",
      "  };"
    )

    if (deriv2) {
      externC <- c(
        externC,
        "  auto IDX2 = [n_phi_rows, n_sens](int g, int v1, int v2) {",
        "    return g + n_phi_rows * (v1 + n_sens * v2);",
        "  };"
      )
    }

    # Build global_to_sens compile-time literal:
    # global index v -> sens index in [0, n_sens_total), or -1 if compile-time-fixed
    n_global <- n_variables + n_params
    fixed_global_idx <- c(fixed_initial_idx, n_variables + fixed_param_idx)
    sens_counter <- 0L
    global_to_sens_vals <- integer(n_global)
    for (v in seq_len(n_global) - 1L) {
      if (v %in% fixed_global_idx) {
        global_to_sens_vals[v + 1L] <- -1L
      } else {
        global_to_sens_vals[v + 1L] <- sens_counter
        sens_counter <- sens_counter + 1L
      }
    }
    global_to_sens_literal <- paste(global_to_sens_vals, collapse = ", ")

    externC <- c(
      externC,
      "",
      "  // global_to_sens[v]: global index v -> compile-time sens index, or -1 if compile-time-fixed",
      sprintf("  static const int global_to_sens_arr[%d] = {%s};", n_global, global_to_sens_literal),
      sprintf("  auto global_to_sens = [](int v) -> int { return global_to_sens_arr[v]; };"),
      "",
      "  // Validate sens1ini / sens2ini length against full Phi'/Phi'' shape",
      "  // [n_phi_rows, n_sens] and [n_phi_rows, n_sens, n_sens] respectively.",
      "  if (has_sens1ini && args.n_sens1 != n_phi_rows * n_sens) {",
      "    char _m[256];",
      "    snprintf(_m, sizeof(_m), \"sens1ini has wrong length: expected n_phi_rows * n_sens = %d * %d = %d, got %d\", n_phi_rows, n_sens, n_phi_rows * n_sens, args.n_sens1);",
      "    return res.fail(cppde::RC_ILL_INPUT, _m);",
      "  }"
    )

    if (deriv2) {
      externC <- c(
        externC,
        "  if (has_sens2ini && args.n_sens2 != n_phi_rows * n_sens * n_sens) {",
        "    char _m[256];",
        "    snprintf(_m, sizeof(_m), \"sens2ini has wrong length: expected n_phi_rows * n_sens^2 = %d * %d^2 = %d, got %d\", n_phi_rows, n_sens, n_phi_rows * n_sens * n_sens, args.n_sens2);",
        "    return res.fail(cppde::RC_ILL_INPUT, _m);",
        "  }",
        "  if (has_sens2ini && !has_sens1ini)",
        "    return res.fail(cppde::RC_ILL_INPUT, \"sens2ini requires sens1ini (Phi'' without Phi' is inconsistent)\");"
      )
    }

    externC <- c(externC, "")
  }

  # --- Slab-bind state and parameter tangents -------------------------------
  # On the heap-AD path every state and parameter tangent points into one
  # contiguous slab block per vector; for every other value type it is a stub.
  if (deriv) {
    externC <- c(
      externC,
      sprintf("  cppde::detail::tangent_slab<%s> _x_slab;", numType),
      sprintf("  cppde::detail::tangent_slab<%s> _full_params_slab;", numType),
      "  if (n_sens > 0) {",
      sprintf("    _x_slab.prime(x, static_cast<unsigned>(%d), static_cast<unsigned>(n_sens));",
              n_variables),
      sprintf("    _full_params_slab.prime(full_params, static_cast<unsigned>(%d), static_cast<unsigned>(n_sens));",
              n_variables + n_params),
      "  }",
      ""
    )
  }

  # --- initialize states ---
  externC <- c(
    externC,
    "  // initialize variables",
    sprintf("  for (int i = 0; i < %d; ++i) {", n_variables),
    "    bool is_fixed = false;",
    "    (void)is_fixed;  // suppress unused warning"
  )

  if (deriv && length(fixed_initial_idx) > 0) {
    externC <- c(
      externC,
      sprintf(
        "    is_fixed = (%s);",
        paste(sprintf("i == %d", fixed_initial_idx), collapse = " || ")
      )
    )
  }

  if (deriv2) {
    externC <- c(
      externC,
      "    x[i].x().x() = args.params[i];",
      "    if (!is_fixed) {",
      "      int si = global_to_sens(i);  // compile-time sens index (-1 if compile-time-fixed)",
      "      int ai = (si >= 0) ? active_idx[si] : -1;  // active index (-1 if any-fixed)",
      "      if (ai >= 0) {",
      "        // First-order sensitivities (inner layer): seed from Phi' or identity",
      "        if (has_sens1ini) {",
      "          if (n_sens > 0) {  // M=0 under reparam: leave F default (no propagation)",
      sprintf("            x[i].x().diff(0%s);  // allocate n_sens components", dyn_arg),
      "            for (int av = 0; av < n_sens; ++av)",
      "              x[i].x().d(av) = sens1ini[IDX1(i, av)];",
      "          }",
      "        } else {",
      sprintf("          x[i].x().diff(ai%s);  // identity: d(ai) = 1", dyn_arg),
      "        }",
      "        // Second-order sensitivities (outer layer): seed from Phi'' or zeros.",
      "        // Note: .diff(idx) is DESTRUCTIVE (zeros tangents, sets [idx]=1,",
      "        // depend=true). We therefore call diff(0) ONCE per layer to arm",
      "        // depend, then use .d() (non-destructive accessor) to write values.",
      "        if (has_sens2ini) {",
      "          if (n_sens > 0) {",
      sprintf("            x[i].diff(0%s);  // arm outer m_depend", dyn_arg),
      "            for (int av1 = 0; av1 < n_sens; ++av1) {",
      sprintf("              x[i].d(av1).diff(0%s);  // arm inner m_depend of m_diff[av1]", dyn_arg),
      "              x[i].d(av1).x() = sens1ini[IDX1(i, av1)];  // first-order value",
      "              for (int av2 = 0; av2 < n_sens; ++av2)",
      "                x[i].d(av1).d(av2) = sens2ini[IDX2(i, av1, av2)];",
      "            }",
      "          }",
      "        } else {",
      sprintf("          x[i].diff(ai%s);  // identity/allocate outer (inner m_depend stays false)", dyn_arg),
      "        }",
      "      }",
      "    }"
    )
  } else if (deriv) {
    externC <- c(
      externC,
      "    x[i] = args.params[i];",
      "    if (!is_fixed) {",
      "      int si = global_to_sens(i);  // compile-time sens index (-1 if compile-time-fixed)",
      "      int ai = (si >= 0) ? active_idx[si] : -1;  // active index (-1 if any-fixed)",
      "      if (ai >= 0) {",
      "        if (has_sens1ini) {",
      "          if (n_sens > 0) {  // M=0 under reparam: leave F default",
      "            // Seed from Phi'(theta): row i of sens1ini",
      sprintf("            x[i].diff(0%s);  // allocate n_sens components", dyn_arg),
      "            for (int av = 0; av < n_sens; ++av)",
      "              x[i].d(av) = sens1ini[IDX1(i, av)];",
      "          }",
      "        } else {",
      sprintf("          x[i].diff(ai%s);  // identity: d(ai) = 1", dyn_arg),
      "        }",
      "      }",
      "    }"
    )
  } else {
    externC <- c(externC, "    x[i] = args.params[i];")
  }

  externC <- c(
    externC,
    "    full_params[i] = x[i];",
    "  }",
    "",
    "  // initialize parameters",
    sprintf("  for (int i = 0; i < %d; ++i) {", n_params),
    sprintf("    int param_index = %d + i;", n_variables),
    "    bool is_fixed = false;",
    "    (void)is_fixed;  // suppress unused warning"
  )

  if (deriv && length(fixed_param_idx) > 0) {
    externC <- c(
      externC,
      sprintf(
        "    is_fixed = (%s);",
        paste(sprintf("i == %d", fixed_param_idx), collapse = " || ")
      )
    )
  }

  if (deriv2) {
    externC <- c(
      externC,
      "    int global_idx = n_states + i;",
      "    full_params[param_index].x().x() = args.params[param_index];",
      "    if (!is_fixed) {",
      "      int si = global_to_sens(global_idx);  // compile-time sens index",
      "      int ai = (si >= 0) ? active_idx[si] : -1;  // active index",
      "      if (ai >= 0) {",
      "        // First-order (inner layer): seed from Phi' or identity fallback",
      "        if (has_sens1ini) {",
      "          if (n_sens > 0) {  // M=0 under reparam: leave F default",
      sprintf("            full_params[param_index].x().diff(0%s);  // allocate n_sens components", dyn_arg),
      "            for (int av = 0; av < n_sens; ++av)",
      "              full_params[param_index].x().d(av) = sens1ini[IDX1(global_idx, av)];",
      "          }",
      "        } else {",
      sprintf("          full_params[param_index].x().diff(ai%s);  // identity: dp_i/dp_j = delta_ij", dyn_arg),
      "        }",
      "        // Second-order (outer layer): see note in state-seed block.",
      "        if (has_sens2ini) {",
      "          if (n_sens > 0) {",
      sprintf("            full_params[param_index].diff(0%s);  // arm outer m_depend", dyn_arg),
      "            for (int av1 = 0; av1 < n_sens; ++av1) {",
      sprintf("              full_params[param_index].d(av1).diff(0%s);  // arm inner m_depend", dyn_arg),
      "              full_params[param_index].d(av1).x() = sens1ini[IDX1(global_idx, av1)];",
      "              for (int av2 = 0; av2 < n_sens; ++av2)",
      "                full_params[param_index].d(av1).d(av2) = sens2ini[IDX2(global_idx, av1, av2)];",
      "            }",
      "          }",
      "        } else {",
      sprintf("          full_params[param_index].diff(ai%s);  // identity/allocate outer", dyn_arg),
      "        }",
      "      }",
      "    }"
    )
  } else if (deriv) {
    externC <- c(
      externC,
      "    int global_idx = n_states + i;",
      "    full_params[param_index] = args.params[param_index];",
      "    if (!is_fixed) {",
      "      int si = global_to_sens(global_idx);  // compile-time sens index",
      "      int ai = (si >= 0) ? active_idx[si] : -1;  // active index",
      "      if (ai >= 0) {",
      "        if (has_sens1ini) {",
      "          if (n_sens > 0) {  // M=0 under reparam: leave F default",
      "            // Seed from Phi'(theta): row n_states + i of sens1ini",
      sprintf("            full_params[param_index].diff(0%s);  // allocate n_sens components", dyn_arg),
      "            for (int av = 0; av < n_sens; ++av)",
      "              full_params[param_index].d(av) = sens1ini[IDX1(global_idx, av)];",
      "          }",
      "        } else {",
      "          // Identity fallback: dp_i/dp_j = delta_{ij}",
      sprintf("          full_params[param_index].diff(ai%s);", dyn_arg),
      "        }",
      "      }",
      "    }"
    )
  } else {
    externC <- c(externC, "    full_params[param_index] = args.params[param_index];")
  }

  externC <- c(externC, "  }", "")
  # --- Forcing Initialization ---
  externC <- c(externC, forcing_init_code, "")

  externC <- c(externC,
               "  // --- Copy integration times ---",
               "  std::vector<double> times_dbl(args.times, args.times + args.n_times);",
               ""
  )

  if (includeTimeZero) {
    externC <- c(externC,
                 "  // ensure time zero is included",
                 "  if (std::find(times_dbl.begin(), times_dbl.end(), 0.0) == times_dbl.end()) {",
                 "    times_dbl.push_back(0.0);",
                 "  }",
                 ""
    )
  }

  externC <- c(externC,
               "  // sort times ascending and remove duplicates",
               "  std::sort(times_dbl.begin(), times_dbl.end());",
               "  times_dbl.erase(std::unique(times_dbl.begin(), times_dbl.end()), times_dbl.end());",
               "",
               "  // convert to AD vector",
               sprintf("  std::vector<%s> times;", numType),
               "  times.reserve(times_dbl.size());",
               "  for (double tval : times_dbl) {",
               "    times.emplace_back(tval);",
               "  }",
               "",
               "  // storage for results",
               sprintf("  std::vector<%s> result_times;", numType),
               sprintf("  std::vector<%s> y;", numType),
               "  result_times.reserve(times_dbl.size());",
               sprintf("  y.reserve(times_dbl.size() * %d);", n_variables),
               "",
               "  // --- Event containers ---",
               sprintf("  std::vector<FixedEvent<cppde::vector_t<%s>, %s>> fixed_events;", numType, numType),
               sprintf("  std::vector<RootEvent<cppde::vector_t<%s>, %s>> root_events;", numType, numType)
  )

  # Insert event code from Python
  if (event_code != "") {
    externC <- c(externC, event_code)
  }

  # Note: rootfunc_code is inserted later, after sys is defined

  # --- Integration setup ---
  # Dense or sparse LU, Rosenbrock4 or one of the two multistep methods, which
  # instantiate the same cppde::multistepper. NDF against BDF is a runtime flag.
  resizer_tag   <- "cppde::initially_resizer"

  # Generate the C++ stepper type for value type V and LU pattern J.
  # Only meaningful for is_multistep(method), callers on the rb4 path
  # never invoke this helper.
  make_stepper_type <- function(V, J) {
    method_enum <- switch(method,
      "bdf"   = "cppde::multistep_method::bdf",
      "adams" = "cppde::multistep_method::adams",
      stop(sprintf("internal: unhandled multistep method '%s'", method))
    )
    sprintf("cppde::multistepper<%s, %s, %s, %s>",
            method_enum, V, J, resizer_tag)
  }

  # make_stepper_type() errors on a method name that is not multistep, so the
  # stepper type strings are built lazily, only for is_multistep(method).
  ms_double <- ms_AD <- ms_AD2 <- NULL
  if (use_sparse) {
    rb4_double <- "rosenbrock4<double, sparse_lu_tag>"
    rb4_AD     <- "rosenbrock4<AD, sparse_lu_tag>"
    rb4_AD2    <- "rosenbrock4<AD2, sparse_lu_tag>"
    if (is_multistep(method)) {
      ms_double <- make_stepper_type("double", "sparse_lu_tag")
      ms_AD     <- make_stepper_type("AD",     "sparse_lu_tag")
      ms_AD2    <- make_stepper_type("AD2",    "sparse_lu_tag")
    }
  } else {
    rb4_double <- "rosenbrock4<double>"
    rb4_AD     <- "rosenbrock4<AD>"
    rb4_AD2    <- "rosenbrock4<AD2>"
    if (is_multistep(method)) {
      ms_double <- make_stepper_type("double", "cppde::dense_lu_tag")
      ms_AD     <- make_stepper_type("AD",     "cppde::dense_lu_tag")
      ms_AD2    <- make_stepper_type("AD2",    "cppde::dense_lu_tag")
    }
  }
  # Tsit5 stepper types (no Jacobian pattern, explicit method)
  tsit5_double <- "cppde::tsit5<double>"
  tsit5_AD     <- "cppde::tsit5<AD>"
  tsit5_AD2    <- "cppde::tsit5<AD2>"

  if (is_multistep(method)) {
    # ---- Multistep stepper (bdf / adams) ----
    # cppde::multistepper_controller is templated on the stepper type and works
    # with any multistepper instantiation.
    ms_type <- if (deriv2) ms_AD2 else if (deriv) ms_AD else ms_double
    stepper_setup_lines <- c(
      sprintf("  controlledStepper.stepper().set_use_ndf_kappa(%s);",
              if (useNDF) "true" else "false"),
      # Slab priming for any AD path (heap dual<T,0> or static-N dual<T,N>),
      # a no-op for non-AD and nested-AD types. Emitted before the
      # std::move into denseStepper.
      if (deriv) {
        "  controlledStepper.prepare_sensitivities(static_cast<unsigned>(n_sens));"
      } else {
        character()
      }
    )
    # Termination argument for equilibrate (empty string when not used)
    is_equilibrate <- identical(tolower(rootfunc), "equilibrate")
    termination_arg <- if (is_equilibrate) ", ss_termination" else ""

    # Multistep methods always use the Nordsieck dense output.
    stepper_line <- paste(
      c(sprintf("  auto controlledStepper = cppde::multistepper_controller<%s>(abstol, reltol);",
                ms_type),
        stepper_setup_lines,
        "  auto denseStepper = cppde::multistepper_dense_output<decltype(controlledStepper)>(std::move(controlledStepper));"),
      collapse = "\n"
    )
    integrate_line <- sprintf("  integrate_times_dense(denseStepper, std::make_pair(sys, jac), x, times.begin(), times.end(), dt, obs, fixed_events, root_events, checker, root_tol, maxroot, dt_est%s);", termination_arg)
  } else {
    # ---- Single-step stepper (rb4 or tsit5) ----
    # Both use the generic onestep_controller with Gustafsson PI control.

    # Select the stepper type string based on method and AD level.
    if (method == "tsit5") {
      os_type <- if (deriv2) tsit5_AD2 else if (deriv) tsit5_AD else tsit5_double
    } else {
      # rosenbrock4
      os_type <- if (deriv2) rb4_AD2 else if (deriv) rb4_AD else rb4_double
    }

    # Termination argument for equilibrate (empty string when not used)
    is_equilibrate <- identical(tolower(rootfunc), "equilibrate")
    termination_arg <- if (is_equilibrate) ", cppde::detail::no_dt_estimator{}, ss_termination" else ""

    # Slab priming for the single-step path, as in the multistep branch: the call
    # lands on the controller before the move into the dense wrapper, so the slabs
    # reachable through it are primed for the whole solve. Gated on the value type.
    onestep_prep_line <- if (deriv) {
      "  controlledStepper.prepare_sensitivities(static_cast<unsigned>(n_sens));"
    } else {
      character()
    }

    if (useDenseOutput) {
      stepper_line <- paste(
        c(sprintf("  auto controlledStepper = cppde::onestep_controller<%s>(abstol, reltol);", os_type),
          onestep_prep_line,
          "  auto denseStepper = cppde::onestep_dense_output<decltype(controlledStepper)>(std::move(controlledStepper));"),
        collapse = "\n"
      )
      integrate_line <- sprintf("  integrate_times_dense(denseStepper, std::make_pair(sys, jac), x, times.begin(), times.end(), dt, obs, fixed_events, root_events, checker, root_tol, maxroot%s);", termination_arg)
    } else {
      stepper_line <- paste(
        c(sprintf("  auto controlledStepper = cppde::onestep_controller<%s>(abstol, reltol);", os_type),
          onestep_prep_line),
        collapse = "\n"
      )
      integrate_line <- sprintf("  integrate_times(controlledStepper, std::make_pair(sys, jac), x, times.begin(), times.end(), dt, obs, fixed_events, root_events, checker, root_tol, maxroot%s);", termination_arg)
    }
  }

  # Both dense and sparse paths use the same Jacobian functor name.
  # Codegen produces exactly one 'struct jacobian' matching the path.
  externC <- c(externC, "",
               sprintf("  // --- Solver setup (%s LU) ---",
                       if (use_sparse) "sparse" else "dense"),
               "  double abstol = args.abstol;",
               "  double reltol = args.reltol;",
               "  double root_tol = args.root_tol;",
               "  double hini = args.hini;",
               "  int maxroot = args.maxroot;",
               "  ode_system sys(full_params, F);",
               "  jacobian jac(full_params, F);",
               "  observer obs(result_times, y);")

  # Insert rootfunc code from Python (after sys is defined, needed for equilibrate)
  if (rootfunc_code != "") {
    externC <- c(externC, rootfunc_code)
  }

  # --- Initial step size estimation ---
  # The multistep methods use the cvHin port, the single-step ones
  # estimate_initial_dt. order = 1 suits all of them, see cppde_utils.hpp.
  method_order <- 1L

  if (is_multistep(method)) {
    # cppde_hin, needs only `sys`; no compute_ydd closure, no order param.
    estimate_dt_block <- c(
      sprintf("  // --- Determine initial dt (%s, CVODES cvHin port) ---", method),
      sprintf("  %s dt;", numType),
      "  if (hini == 0.0) {",
      sprintf("    dt = odeint_utils::cppde_hin<%s>(", numType),
      "      sys, x, times.front(),",
      "      odeint_utils::scalar_value(times.back()),",
      "      abstol, reltol);",
      "  } else {",
      "    dt = hini;",
      "  }"
    )
  } else {
    if (is_explicit(method)) {
      compute_ydd_lines <- c(
        sprintf("  // --- compute_ydd (FD fallback for explicit method) ---"),
        "  auto compute_ydd = odeint_utils::make_fd_ydd(sys);"
      )
    } else if (use_sparse) {
      compute_ydd_lines <- c(
        "  // --- compute_ydd (analytic: dfdt + sparse J*f) ---",
        "  auto compute_ydd = [&](const auto& x_, auto t_, const auto& f_,",
        "                          double /*h_trial*/, auto& ydd_) {",
        sprintf("    cppde::csc_matrix<%s> J_init;", numType),
        sprintf("    std::vector<%s> dfdt_(x_.size());", numType),
        "    jac(x_, J_init, t_, dfdt_);",
        "    ydd_ = dfdt_;",
        "    cppde::csc_matvec_add(J_init, f_, ydd_);",
        "  };"
      )
    } else {
      compute_ydd_lines <- c(
        "  // --- compute_ydd (analytic: dfdt + dense J*f) ---",
        "  auto compute_ydd = [&](const auto& x_, auto t_, const auto& f_,",
        "                          double /*h_trial*/, auto& ydd_) {",
        sprintf("    cppde::dense_matrix<%s> J_init(x_.size(), x_.size());", numType),
        sprintf("    std::vector<%s> dfdt_(x_.size());", numType),
        "    jac(x_, J_init, t_, dfdt_);",
        "    ydd_ = dfdt_;",
        "    for (std::size_t i = 0; i < x_.size(); ++i)",
        "      for (std::size_t j = 0; j < x_.size(); ++j)",
        "        ydd_[i] += J_init(i,j) * f_[j];",
        "  };"
      )
    }

    estimate_dt_block <- c(
      compute_ydd_lines,
      sprintf("  // --- Determine initial dt (%s, order=%d) ---", method, method_order),
      sprintf("  %s dt;", numType),
      "  if (hini == 0.0) {",
      sprintf("    dt = odeint_utils::estimate_initial_dt<%s>(", numType),
      "      sys, compute_ydd, x, times.front(),",
      "      odeint_utils::scalar_value(times.back()),",
      sprintf("      abstol, reltol, /*order=*/%d);", method_order),
      "  } else {",
      "    dt = hini;",
      "  }"
    )
  }

  # --- Multistep methods: dt re-estimator lambda for event restarts ---
  # Handed to integrate_times as the DtEstimator, called on every restart.
  # times.back() is the upper-bound hint that keeps the remaining window in view.
  dt_est_block <- character(0)
  if (is_multistep(method)) {
    dt_est_block <- c(
      "  // --- Multistep event restart: re-estimate dt from post-event state ---",
      "  auto dt_est = [&](auto& x_ev, auto t_ev) {",
      sprintf("    return odeint_utils::cppde_hin<%s>(", numType),
      "      sys, x_ev, t_ev,",
      "      odeint_utils::scalar_value(times.back()),",
      "      abstol, reltol);",
      "  };"
    )
  }

  externC <- c(externC,
               stepper_line, "",
               estimate_dt_block,
               dt_est_block,
               "",
               "  // --- Integration (catch recoverable errors for partial results) ---",
               "  std::string solver_message;",
               "  try {",
               paste0("    ", integrate_line),
               "  } catch (const cppde::no_progress_error& e) {",
               "    // StepChecker sets RC_TOO_MUCH_WORK / RC_CONV_FAILURE before",
               "    // throwing, but the dense-output wrappers' failed_step_checker",
               "    // throws the same exception type without touching StepChecker.",
               "    // Guarantee a non-zero return code in that case.",
               "    if (checker.return_code() == cppde::RC_SUCCESS)",
               "      checker.set_return_code(cppde::RC_CONV_FAILURE);",
               "    solver_message = e.what();",
               "  } catch (const std::runtime_error& e) {",
               "    // KLU factor/analyze failures and similar linear-solver issues.",
               "    checker.set_return_code(cppde::RC_LSETUP_FAIL);",
               "    solver_message = e.what();",
               "  } catch (const std::exception& e) {",
               "    // Anything else that derives from std::exception (bad_alloc,",
               "    // overflow_error from a wild RHS, ...), classify as unclassified",
               "    // but still return partial results rather than aborting R.",
               "    checker.set_return_code(cppde::RC_UNRECOGNIZED_ERR);",
               "    solver_message = e.what();",
               "  }",
               "",
               "  // --- Populate diagnostics from available state ---",
               "  if (!result_times.empty()) {",
               sprintf("    checker.set_t_reached(static_cast<double>(%s));",
                       if (deriv2) "result_times.back().x().x()"
                       else if (deriv) "result_times.back().x()"
                       else "result_times.back()"),
               "  }",
               "  // last_dt is set by the integration loop (process_dense/process_controlled)",
               "  // and reflects the true controller step size, not the output grid spacing.",
               "",
               "  const int n_out = static_cast<int>(result_times.size());",
               "  res.n_out       = n_out;",
               "  res.return_code = checker.return_code();",
               "  res.message     = solver_message;",
               "  res.accepted    = checker.n_accepted();",
               "  res.rejected    = checker.n_rejected();",
               "  res.fevals      = checker.n_fevals();",
               "  res.jevals      = checker.n_jevals();",
               "  res.setups      = checker.n_setups();",
               "  res.last_dt     = checker.last_dt();",
               "  res.last_order  = checker.last_order();",
               "  res.t_reached   = checker.t_reached();",
               "",
               "  if (n_out <= 0) {",
               "    if (res.return_code == cppde::RC_SUCCESS)",
               "      res.fail(cppde::RC_UNRECOGNIZED_ERR, \"Integration produced no output\");",
               "    return res.return_code;",
               "  }",
               ""
  )

  # --- Flatten AD results to plain double.  Must happen here: the arena
  # scope is still open, and tangents point into it.
  if (!deriv) {
    externC <- c(externC,
                 sprintf("  res.prepare_out(n_out, 0, %d, false, false);", n_variables),
                 "  for (int i = 0; i < n_out; ++i) {",
                 "    res.t_out[i] = result_times[i];",
                 sprintf("    for (int s = 0; s < %d; ++s)", n_variables),
                 sprintf("      res.v_out[i + (size_t)n_out * s] = y[i * %d + s];", n_variables),
                 "  }")
  } else if (!deriv2) {
    externC <- c(externC,
                 "  const int n_sens_out = n_sens;",
                 sprintf("  res.prepare_out(n_out, n_sens_out, %d, true, false);", n_variables),
                 "  for (int i = 0; i < n_out; ++i) {",
                 "    res.t_out[i] = result_times[i].x();",
                 sprintf("    for (int s = 0; s < %d; ++s) {", n_variables),
                 sprintf("      %s& xi = y[i * %d + s];", numType, n_variables),
                 "      res.v_out[i + (size_t)n_out * s] = xi.x();",
                 "      for (int av = 0; av < n_sens_out; ++av)",
                 sprintf("        res.s1_out[i + (size_t)n_out * (s + %d * av)] = xi.d(av);", n_variables),
                 "    }",
                 "  }")
  } else {
    externC <- c(externC,
                 "  const int n_sens_out = n_sens;",
                 sprintf("  res.prepare_out(n_out, n_sens_out, %d, true, true);", n_variables),
                 "  for (int i = 0; i < n_out; ++i) {",
                 "    res.t_out[i] = result_times[i].x().x();",
                 sprintf("    for (int s = 0; s < %d; ++s) {", n_variables),
                 sprintf("      %s& xi = y[i * %d + s];", numType, n_variables),
                 "      res.v_out[i + (size_t)n_out * s] = xi.x().x();",
                 "      for (int av1 = 0; av1 < n_sens_out; ++av1) {",
                 sprintf("        res.s1_out[i + (size_t)n_out * (s + %d * av1)] = xi.d(av1).x();", n_variables),
                 "        for (int av2 = 0; av2 < n_sens_out; ++av2)",
                 sprintf("          res.s2_out[i + (size_t)n_out * (s + %d * (av1 + n_sens_out * av2))] = xi.d(av1).d(av2);", n_variables),
                 "      }",
                 "    }",
                 "  }")
  }

  externC <- c(externC, "  return res.return_code;")

  # --- End try/catch ---
  externC <- c(externC,
               "  } catch (const std::exception& e) {",
               "    return res.fail(cppde::RC_UNRECOGNIZED_ERR,",
               "                    std::string(\"ODE solver failed: \") + e.what());",
               "  } catch (...) {",
               "    return res.fail(cppde::RC_UNRECOGNIZED_ERR,",
               "                    \"ODE solver failed: unknown C++ exception\");",
               "  }",
               "}",
               "")

  # --- extern "C" entry points: single condition, and the OpenMP batch ---
  dflag  <- if (deriv)  "true" else "false"
  d2flag <- if (deriv2) "true" else "false"
  ## Whether `times` alone fixes the output grid, so the batch can size its
  ## results up front. An event time off the grid adds a row, which pre_acquire
  ## declines; a root event or a rootfunc is dynamic and never qualifies.
  has_root_events <- !is.null(events) && "root" %in% names(events) &&
    any(!is.na(events$root))
  fixed_grid <- !has_root_events && rootfunc_code == ""
  if (!exists("event_time_exprs", inherits = FALSE)) event_time_exprs <- NULL
  event_time_exprs <- if (is.null(event_time_exprs)) NULL else
    as.character(unlist(event_time_exprs))
  ## A model with fixed events needs its event times before the solve; without
  ## them the predicted row count is wrong for every off-grid event and
  ## pre_acquire declines, which costs the staging copy.
  if (fixed_grid && !is.null(events) && is.null(event_time_exprs))
    fixed_grid <- FALSE
  n_ev_times <- length(event_time_exprs)
  if (fixed_grid && n_ev_times > 0L) {
    externC <- c(
      externC,
      sprintf("static void %s_fixed_event_times(const double* params, double* out) {",
              modelname),
      sprintf("  out[%d] = %s;", seq_len(n_ev_times) - 1L, event_time_exprs),
      "}",
      "")
  }
  externC <- c(
    externC,
    sprintf('extern "C" SEXP solve_%s(SEXP timesSEXP, SEXP paramsSEXP, SEXP sens1iniSEXP, SEXP sens2iniSEXP, SEXP fixedSEXP, SEXP abstolSEXP, SEXP reltolSEXP, SEXP maxprogressSEXP, SEXP maxstepsSEXP, SEXP hiniSEXP, SEXP root_tolSEXP, SEXP maxrootSEXP, SEXP forcingTimesSEXP, SEXP forcingValuesSEXP, SEXP dimnamesSEXP) {',
            modelname),
    "  cppde::rbatch::solve_args a = cppde::rbatch::read_solve_args(",
    "      timesSEXP, paramsSEXP, sens1iniSEXP, sens2iniSEXP, fixedSEXP,",
    "      abstolSEXP, reltolSEXP, maxprogressSEXP, maxstepsSEXP, hiniSEXP,",
    "      root_tolSEXP, maxrootSEXP, forcingTimesSEXP, forcingValuesSEXP);",
    sprintf("  return cppde::rbatch::solve_one(a, %d, %s, %s, &solve_impl, dimnamesSEXP);",
            n_variables, dflag, d2flag),
    "}",
    "",
    sprintf('extern "C" SEXP solve_%s_batch(SEXP condsSEXP, SEXP nthreadsSEXP, SEXP dimnamesSEXP) {', modelname),
    "  const int K  = Rf_length(condsSEXP);",
    "  const int nt = Rf_isNull(nthreadsSEXP) ? 0 : INTEGER(nthreadsSEXP)[0];",
    if (fixed_grid) "" else "  (void)dimnamesSEXP;",
    "",
    "  // Phase A: every R read and every R allocation, serially.",
    "  std::vector<cppde::rbatch::solve_args> A;",
    "  A.reserve(K);",
    "  for (int k = 0; k < K; ++k)",
    "    A.push_back(cppde::rbatch::read_cond_args(VECTOR_ELT(condsSEXP, k)));",
    "  std::vector<cppde::rbatch::solve_result> R_(K);",
    "  SEXP out = PROTECT(Rf_allocVector(VECSXP, K));",
    if (fixed_grid) paste0(
      "\n  // No events and no root finding, so n_out follows from `times`:\n",
      "  // allocate the results now and let the workers gather straight into\n",
      "  // them.  Saves the staging copy and puts the first-touch page faults\n",
      "  // inside the parallel region.\n",
      if (n_ev_times > 0L) sprintf(
        paste0("  std::vector<std::vector<double> > EV(K, std::vector<double>(%d));\n",
               "  for (int k = 0; k < K; ++k) %s_fixed_event_times(A[k].params, EV[k].data());\n"),
        n_ev_times, modelname) else "",
      sprintf("  std::vector<cppde::rbatch::pre_ctx> P = cppde::rbatch::prealloc_batch(\n      out, A, %d, %d, %s, %s, %s, dimnamesSEXP%s);\n",
              n_variables, n_total_sens, dflag, d2flag,
              if (includeTimeZero) "true" else "false",
              if (n_ev_times > 0L) ", &EV" else ""),
      "  for (int k = 0; k < K; ++k) {\n",
      "    R_[k].acquire = &cppde::rbatch::pre_acquire;\n",
      "    R_[k].acquire_ctx = &P[k];\n",
      "  }") else "",
    "",
    "  // Phase B: no R API, no R allocation.",
    "  cppde::rbatch::run_batch(K, nt, [&](int k) noexcept { solve_impl(A[k], R_[k]); });",
    "",
    "  // Phase C: back on the R side.",
    "  for (int k = 0; k < K; ++k) {",
    if (fixed_grid) sprintf(
      "    if (R_[k].used_sink) { cppde::rbatch::finish_prealloc(out, k, R_[k], %s, %s); continue; }",
      dflag, d2flag) else "",
    sprintf("    SEXP e = PROTECT(cppde::rbatch::build_result_sexp(R_[k], %d, %s, %s));",
            n_variables, dflag, d2flag),
    "    SET_VECTOR_ELT(out, k, e);",
    "    UNPROTECT(1);",
    "  }",
    "  cppde::rbatch::set_threads_attr(out, K, nt);",
    "  UNPROTECT(1);",
    "  return out;",
    "}")

  externC <- paste(externC, collapse = "\n")

  # --- Write C++ file ---
  if (!dir.exists(outdir)) {
    stop("outdir does not exist: ", outdir)
  }

  filename <- file.path(outdir, paste0(modelname, ".cpp"))

  # Warn if file already exists
  if (file.exists(filename)) {
    message("Overwriting existing file: ", normalizePath(filename, winslash = "/", mustWork = FALSE))
  }

  cpp_text <- c(
    paste0("/** Code auto-generated by cppDE ", as.character(utils::packageVersion("cppDE")), " **/"),
    "", includings, "", usings, "", "namespace {",
    ode_code, "", jac_code,
    "", observer_code,
    "", "}", "", externC
  )

  writeLines(cpp_text, filename, useBytes = TRUE)

  if (verbose) message("Wrote: ", normalizePath(filename, winslash = "/", mustWork = FALSE))

  # --- Attach attributes ---
  # Reconstruct character Jacobian from sparse triplet format
  jac_matrix_R <- matrix("0", nrow = n_variables, ncol = n_variables,
                         dimnames = list(variables, variables))
  jac_rows <- codegen_result$jac_nnz_rows
  jac_cols <- codegen_result$jac_nnz_cols
  jac_exprs <- codegen_result$jac_nnz_exprs
  if (length(jac_rows) > 0L) {
    jac_matrix_R[cbind(as.integer(jac_rows) + 1L,
                       as.integer(jac_cols) + 1L)] <- as.character(jac_exprs)
  }

  attr(modelname, "equations")     <- rhs
  attr(modelname, "srcfile")       <- normalizePath(filename, winslash = "/", mustWork = FALSE)
  attr(modelname, "variables")     <- variables
  attr(modelname, "parameters")    <- params
  attr(modelname, "forcings")      <- forcings
  attr(modelname, "events")        <- events
  attr(modelname, "rootfunc")      <- rootfunc
  attr(modelname, "fixed")         <- c(fixed_initials, fixed_params)
  attr(modelname, "jacobian")      <- list(f.x = jac_matrix_R, f.time = time_derivs_str)
  attr(modelname, "deriv")         <- deriv
  attr(modelname, "deriv2")        <- deriv2
  attr(modelname, "nStack")        <- if (is_heap) Inf else as.numeric(nStack_width)
  attr(modelname, "sparse")        <- use_sparse
  attr(modelname, "method")        <- method
  attr(modelname, "useNDF")        <- useNDF
  # Dimension names: under reparametrization the sens columns are theta slots.
  # The sens dim defaults to model-parameter names; solveODE() overrides it per
  # call when sens1ini carries a full Phi' shape.
  if (deriv) {
    attr(modelname, "dimNames") <- list(
      time = "time",
      variable = variables,
      sens = c(sens_initials, sens_params)
    )
  } else {
    attr(modelname, "dimNames") <- list(
      time = "time",
      variable = variables
    )
  }

  # --- Build compile arguments ---
  compile_args <- character(0)
  if (isTRUE(profile))   compile_args <- c(compile_args, "-DCPPDE_PROFILE")
  if (isTRUE(stepTrace)) compile_args <- c(compile_args, "-DCPPDE_STEP_TRACE")

  # KLU auto-tuning: pass codegen-determined settings as compile-time defines
  klu_settings <- codegen_result$klu_settings
  if (!is.null(klu_settings)) {
    klu_btf <- if (isTRUE(klu_settings$use_btf)) 1L else 0L
    klu_ord <- as.integer(klu_settings$ordering)
    compile_args <- c(compile_args,
                      sprintf("-DKLUBTF=%d", klu_btf),
                      sprintf("-DKLUAMD=%d", klu_ord))
    if (verbose) {
      message(sprintf("  KLU settings (codegen): BTF=%s, ordering=%s (nblocks=%d, cv=%.2f)",
                      if (klu_btf) "on" else "off",
                      klu_settings$ordering_name,
                      as.integer(klu_settings$nblocks),
                      klu_settings$cv_row_degree))
    }
  }

  # Record OpenMP on the model, not just inside compile(): a consumer that
  # builds the shared object itself (dMod) drives its flags from these
  # attributes, and without them the batch entry compiles to a serial loop.
  if (isTRUE(cvodeConfig$openmp_available)) {
    compile_args <- c(compile_args, cvodeConfig$openmp_cxxflags)
    attr(modelname, "linkArgs") <- cvodeConfig$openmp_libs
  }

  attr(modelname, "compileArgs") <- paste(compile_args, collapse = " ")

  if (compile) {
    compile(modelname, verbose = verbose)
  }
  return(modelname)
}
