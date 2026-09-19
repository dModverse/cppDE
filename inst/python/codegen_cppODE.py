"""C++ generators of the native cppDE solvers.

Public API: generate_ode_cpp, generate_event_code, generate_rootfunc_code,
generate_forcing_init_code, fixed_event_time_exprs, decide_sparse and
analyze_klu_settings. Expressions and derivatives come from cppde_graph via
cppde_model.
"""

import math
import time

import cppde_model

# =====================================================================
# Main ODE generator
# =====================================================================

def generate_ode_cpp(
    rhs_dict,
    params_list,
    num_type="double",
    fixed_states=None,
    fixed_params=None,
    forcings_list=None,
    sparse=None,
    skip_jacobian=False,
    ad_level=0,
    arena=False,
    emit_contractions=False,
    emit_jvp=False):
    """C++ of `ode_system`, `jacobian` and (reverse) `adjoint_terms`.

    Args:
        rhs_dict: state name -> right-hand side string.
        params_list, forcings_list: symbol names.
        num_type, ad_level, arena: the C++ scalar (see cppde_model.Scalar).
        sparse: None (decide), True or False.
        skip_jacobian: emit a Jacobian stub (explicit methods).
        emit_contractions: emit adjoint_terms; emit_jvp adds its jvp_* and
            dfdt_* pairs (reverse mode with Rosenbrock4).
        fixed_states, fixed_params: unused.

    Returns:
        dict with ode_code, jac_code, adj_code, data_code, jac_nnz_rows,
        jac_nnz_cols, states, params, forcings, use_sparse, sparsity_stats,
        klu_settings and codegen_stats.
    """
    t0 = time.perf_counter()
    scalar = cppde_model.Scalar(num_type, ad_level, arena)
    model = cppde_model.model_for(rhs_dict, params_list, forcings_list)
    n = model.n
    stats = {"graph_nodes": 0, "parse_seconds": time.perf_counter() - t0}

    ode_lines = cppde_model.cpp_ode_system(model, scalar)

    pairs = [] if skip_jacobian else [
        (i, j) for i, row in enumerate(model.pattern()) for j in row]
    rows = [p[0] for p in pairs]
    cols = [p[1] for p in pairs]
    nnz = len(pairs)
    use_sparse = decide_sparse(sparse, n, nnz, has_jacobian=not skip_jacobian)
    n2 = n * n
    sparsity_stats = {
        "n": n,
        "jac_nnz": nnz,
        "jac_zeros_pct": 100.0 * (1.0 - nnz / n2) if n2 else 0,
        "jac_pattern": sorted(pairs),
    }

    if skip_jacobian:
        jac_lines = cppde_model.cpp_noop_jacobian(scalar)
    else:
        strategy = model.jacobian_strategy()
        stats["strategy"] = strategy
        stats["entry_work"], stats["rhs_nodes"] = model.entry_work()
        if strategy == "colour":
            stats["colours"] = model.colouring()[1]
        jac_lines = cppde_model.cpp_jacobian(model, scalar, use_sparse, strategy)

    adj_lines = []
    if emit_contractions:
        adj_lines = cppde_model.cpp_adjoint_terms(model, scalar, emit_jvp)

    klu_settings = analyze_klu_settings(n, rows, cols) if use_sparse else None
    data_lines = model.linmap.cpp() if model.linmap is not None else []
    stats["lin_rows"] = model.nlin
    stats["lin_nnz"] = model.linmap.nnz() if model.linmap is not None else 0
    stats["graph_nodes"] = len(model.g)
    stats["seconds"] = time.perf_counter() - t0

    return {
        "ode_code": "\n".join(ode_lines),
        "jac_code": "\n".join(jac_lines),
        "adj_code": "\n".join(adj_lines),
        "data_code": "\n".join(data_lines),
        "jac_nnz_rows": rows,
        "jac_nnz_cols": cols,
        "states": model.states,
        "params": model.params,
        "forcings": model.forcings,
        "use_sparse": use_sparse,
        "sparsity_stats": sparsity_stats,
        "klu_settings": klu_settings,
        "codegen_stats": stats,
    }


def generate_forcing_init_code(n_forcings, num_type=None):
    """Generate C++ code to initialize PchipForcing objects from R raw data."""
    return [
        "",
        "  // --- Initialize forcings (PCHIP interpolation) ---",
        f"  const int n_forcings = static_cast<int>(args.flen.size());",
        f"  std::vector<cppde::PchipForcing<{num_type}>> forcing_storage(n_forcings);",
        f"  std::vector<const cppde::PchipForcing<{num_type}>*> F(n_forcings);",
        "",
        "  for (int fi = 0; fi < n_forcings; ++fi) {",
        "    const int n_points = args.flen[fi];",
        "",
        "    std::vector<double> ftimes(args.ftimes[fi], args.ftimes[fi] + n_points);",
        "    std::vector<double> fvalues(args.fvalues[fi], args.fvalues[fi] + n_points);",
        "",
        "    forcing_storage[fi].initialize(ftimes, fvalues);",
        "    F[fi] = &forcing_storage[fi];",
        "  }",
        "",
    ]


# =====================================================================
# Events and root functions
# =====================================================================

def _event_model(states_list, params_list, forcings_list, rhs_dict=None):
    states = cppde_model.as_list(states_list)
    rhs = ({s: rhs_dict[s] for s in states} if rhs_dict is not None
           else {s: "0" for s in states})
    return cppde_model.model_for(rhs, params_list, forcings_list)


def fixed_event_time_exprs(events_df, states_list, params_list, n_states,
                           forcings_list=None):
    """Double expressions of the fixed-event times over the flat vector
    `params`, one per fixed event.

    Returns:
        list of str, or None when the model has no fixed event, has a root
        event, or an event time reads a state, the clock or a forcing.
    """
    rows = cppde_model.event_rows(events_df)
    if not rows:
        return None
    model = _event_model(states_list, params_list, forcings_list)
    return cppde_model.fixed_event_times(model, rows)


def generate_event_code(events_df, states_list, params_list, n_states,
                        num_type="double", forcings_list=None, rhs_dict=None,
                        ad_level=0, arena=False, emit_adjoint=False):
    """C++ lines filling `fixed_events` and `root_events`, or with
    `emit_adjoint` the struct `event_adjoint_terms`.

    Args:
        events_df: table with var, value, time or root, method, terminal,
            direction.
        rhs_dict: right-hand side; without it root events get no G_tt.
        num_type, ad_level, arena: the C++ scalar (see cppde_model.Scalar).

    Returns:
        list of str.
    """
    scalar = cppde_model.Scalar(num_type, ad_level, arena)
    rows = cppde_model.event_rows(events_df)
    if not rows:
        if emit_adjoint:
            return cppde_model.empty_event_adjoint_terms(
                n_states, len(cppde_model.as_list(params_list)), scalar.name)
        return []
    model = _event_model(states_list, params_list, forcings_list, rhs_dict)
    code = cppde_model.EventCode(model, scalar, has_rhs=rhs_dict is not None)
    return code.adjoint(rows) if emit_adjoint else code.forward(rows)


def generate_rootfunc_code(rootfunc, states_list, params_list, n_states,
                           num_type="double", forcings_list=None,
                           ad_level=0, arena=False):
    """C++ lines of terminal root functions ('equilibrate' or expressions)."""
    scalar = cppde_model.Scalar(num_type, ad_level, arena)
    if rootfunc is None:
        return []
    model = _event_model(states_list, params_list, forcings_list)
    return cppde_model.rootfunc_code(model, rootfunc, scalar)


# =====================================================================
# Sparse LU pattern code generation
# =====================================================================

def decide_sparse(sparse, n_states, jac_nnz, has_jacobian=True):
    """Sparse (True) or dense (False) linear solver, for both backends.

    `sparse` is None (auto) or pinned by R, which also handles a missing KLU.
    Auto picks sparse from 8 states at a Jacobian density <= 0.4 (bounds from
    benchmarks/run-benchmarks.R --sparse-sweep); without a Jacobian, dense.
    """
    if sparse is not None:
        return bool(sparse)
    if not has_jacobian or n_states <= 0:
        return False
    density = jac_nnz / float(n_states * n_states)
    return (n_states >= 8) and (density <= 0.4)


def analyze_klu_settings(n, jac_nnz_rows, jac_nnz_cols):
    """KLU settings from the Jacobian pattern: BTF when its digraph has more
    than one strongly connected component; COLAMD when the coefficient of
    variation of the row degrees exceeds 0.5, else AMD.

    Returns:
        dict with use_btf, ordering (0 = AMD, 1 = COLAMD), nblocks,
        ordering_name, cv_row_degree and mean_row_degree; only the first two
        for an empty pattern.
    """
    rows = list(jac_nnz_rows)
    cols = list(jac_nnz_cols)

    if n == 0 or len(rows) == 0:
        return {"use_btf": False, "ordering": 0}

    # --- BTF: strongly connected components via Tarjan's algorithm ---
    # Build adjacency list from (row, col) pairs
    adj = [[] for _ in range(n)]
    for r, c in zip(rows, cols):
        if r != c:  # skip self-loops for SCC analysis
            adj[r].append(c)

    # Iterative Tarjan's SCC
    index_counter = [0]
    stack = []
    on_stack = [False] * n
    index = [-1] * n
    lowlink = [0] * n
    n_scc = [0]

    def _strongconnect_iter(v):
        """Iterative Tarjan's SCC."""
        work_stack = [(v, 0)]  # (node, neighbor_index)
        index[v] = lowlink[v] = index_counter[0]
        index_counter[0] += 1
        stack.append(v)
        on_stack[v] = True

        while work_stack:
            v, ni = work_stack[-1]
            if ni < len(adj[v]):
                work_stack[-1] = (v, ni + 1)
                w = adj[v][ni]
                if index[w] == -1:
                    index[w] = lowlink[w] = index_counter[0]
                    index_counter[0] += 1
                    stack.append(w)
                    on_stack[w] = True
                    work_stack.append((w, 0))
                elif on_stack[w]:
                    lowlink[v] = min(lowlink[v], index[w])
            else:
                # Done with v's neighbors
                if lowlink[v] == index[v]:
                    # Pop SCC
                    while True:
                        w = stack.pop()
                        on_stack[w] = False
                        if w == v:
                            break
                    n_scc[0] += 1
                work_stack.pop()
                if work_stack:
                    w = v
                    v = work_stack[-1][0]
                    lowlink[v] = min(lowlink[v], lowlink[w])

    for v in range(n):
        if index[v] == -1:
            _strongconnect_iter(v)

    nblocks = n_scc[0]
    use_btf = nblocks > 1

    # --- Ordering: row degree variance heuristic ---
    row_degrees = [0] * n
    for r in rows:
        row_degrees[r] += 1

    mean_deg = sum(row_degrees) / n
    var_deg = sum((d - mean_deg) ** 2 for d in row_degrees) / n
    cv = math.sqrt(var_deg) / mean_deg if mean_deg > 0 else 0.0

    # High CV means an irregular pattern (pathway models with hub nodes): COLAMD
    # Low CV means a uniform pattern (PDE stencils): AMD
    ordering = 1 if cv > 0.5 else 0
    ordering_name = "COLAMD" if ordering == 1 else "AMD"

    return {
        "use_btf": use_btf,
        "ordering": ordering,
        "nblocks": nblocks,
        "ordering_name": ordering_name,
        "cv_row_degree": float(cv),
        "mean_row_degree": float(mean_deg),
    }
