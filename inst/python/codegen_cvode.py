"""CVODE / CVODES backend of cppDE.

`generate_cvode_cpp` writes a C++ source that compiles against SUNDIALS and
exposes `solve_<name>` with the SEXP signature of the native models, so that
`solveODE()` calls it unchanged. Supported: BDF and Adams, dense or KLU
Jacobian, forward sensitivities (CVODES) and their adjoint (ASA), PCHIP
forcings, a terminating rootfunc and time- and root-triggered events.
Second-order sensitivities are not supported.

All derivatives come from the expression graph (cppde_graph); a sensitivity
right-hand side is one Jacobian-vector product.
"""

import os
import time

import cppde_model
from codegen_cppODE import analyze_klu_settings, decide_sparse


def _as_list(x):
    if x is None:
        return []
    if isinstance(x, str):
        return [x]
    return list(x)


def _valid(val):
    """False for missing event entries: None, NaN, R's NA (a bool through
    reticulate) and the placeholder strings."""
    if val is None or isinstance(val, bool):
        return False
    if isinstance(val, float) and val != val:
        return False
    if isinstance(val, str) and val.strip().lower() in {"", "na", "nan", "none"}:
        return False
    return True


def generate_cvode_cpp(
    rhs_dict,
    params_list,
    modelname,
    outdir,
    deriv=False,
    reverse=False,
    asa_checkpoints=200,
    fixed_states=None,
    fixed_params=None,
    sparse=None,
    lapack=None,
    method="bdf",
    forcings_list=None,
    events=None,
    rootfunc=None,
    include_time_zero=True,
    version="unknown",
):
    """Write `<outdir>/<modelname>.cpp` for the CVODE backend.

    Returns:
        dict with srcfile, variables, parameters, forcings, sens_names,
        jac_nnz_rows, jac_nnz_cols, use_sparse, use_lapack, compile_defs and
        codegen_stats.
    """
    t_start = time.perf_counter()
    fixed_states = _as_list(fixed_states)
    fixed_params = _as_list(fixed_params)
    params_list = _as_list(params_list)
    forcings_list = _as_list(forcings_list)

    model = cppde_model.model_for(rhs_dict, params_list, forcings_list)
    states_list = model.states
    n_states = model.n
    n_params = len(params_list)
    n_global = n_states + n_params
    n_forcings = len(forcings_list)
    g = model.g

    if method not in ("bdf", "adams"):
        raise ValueError(f"method must be 'bdf' or 'adams', got {method!r}")
    cv_method = "CV_BDF" if method == "bdf" else "CV_ADAMS"

    pattern = model.pattern()
    jac_rows = [i for i, row in enumerate(pattern) for _ in row]
    jac_cols = [j for row in pattern for j in row]
    use_sparse = decide_sparse(sparse, n_states, len(jac_rows))
    klu_settings = (analyze_klu_settings(n_states, jac_rows, jac_cols)
                    if use_sparse else None)
    use_lapack = (not use_sparse) and bool(lapack)
    strategy = model.jacobian_strategy()

    sens_ic_names = [s for s in states_list if s not in fixed_states]
    sens_pr_names = [p for p in params_list if p not in fixed_params]
    sens_names = sens_ic_names + sens_pr_names
    sens_params = {params_list.index(p) for p in sens_pr_names}

    def lines(stmts, indent="    "):
        return "\n".join(cppde_model.cvode_lines(model, stmts, indent))

    ode_body = lines(cppde_model.cvode_rhs_statements(model))
    jac_stmts = cppde_model.jacobian_statements(model, use_sparse, strategy, cvode=True)
    jac_body = lines(jac_stmts)
    sens_body = lines(cppde_model.cvode_sens_statements(
        model, sens_params)) if deriv else ""
    nnz_total = len(jac_rows)

    # --- rootfunc ---
    rootfunc_mode = "none"
    rootfunc_nodes = []
    if rootfunc is not None:
        if isinstance(rootfunc, str):
            if rootfunc.strip().lower() == "equilibrate":
                rootfunc_mode = "equilibrate"
                rootfunc_list = []
            else:
                rootfunc_list = [rootfunc]
                rootfunc_mode = "user"
        elif isinstance(rootfunc, (list, tuple)):
            rootfunc_list = list(rootfunc)
            rootfunc_mode = "user"
        else:
            raise ValueError(
                f"rootfunc must be 'equilibrate' or a character vector, got {type(rootfunc)}")
        for expr_str in rootfunc_list:
            s = str(expr_str).strip()
            if s:
                rootfunc_nodes.append(model.parse(s))

    # --- events ---
    # Each event has the post-event value g of its state, the method
    # folded in, with its partials in x, p and t; a time event adds the time
    # and its dt/dp, a root event the condition r and its partials.
    ev = cppde_model.CvodeEvent(model)
    time_events = []
    root_events = []
    if events is not None:
        ev_dict = events.to_dict("list") if hasattr(events, "to_dict") else events
        list_lens = [len(v) for v in ev_dict.values() if isinstance(v, (list, tuple))]
        n_ev = max(list_lens) if list_lens else 1

        def _get(key, i):
            v = ev_dict.get(key)
            if isinstance(v, (list, tuple)):
                return v[i] if i < len(v) else None
            return v

        for i in range(n_ev):
            var_raw = _get("var", i)
            if var_raw is None:
                continue
            var_name = str(var_raw)
            if var_name not in states_list:
                raise ValueError(f"Event {i}: unknown state variable '{var_name}'")
            var_idx = states_list.index(var_name)
            time_raw = _get("time", i)
            root_raw = _get("root", i)
            if not _valid(time_raw) and not _valid(root_raw):
                raise ValueError(f"Event {i}: either 'time' or 'root' is required")
            if _valid(time_raw) and _valid(root_raw):
                raise ValueError(f"Event {i}: specify exactly one of 'time' or 'root'")
            value_raw = _get("value", i)
            if not _valid(value_raw):
                raise ValueError(f"Event {i}: 'value' is required")
            method_raw = _get("method", i)
            ev_method = str(method_raw).lower() if _valid(method_raw) else "replace"
            if ev_method not in ("replace", "add", "multiply"):
                raise ValueError(
                    f"Event {i}: method must be replace/add/multiply, got {ev_method!r}")

            h = model.parse(str(value_raw))
            x_var = g.state(var_idx)
            gn = {"replace": h, "add": g.add(x_var, h),
                  "multiply": g.mul(x_var, h)}[ev_method]
            item = {
                "var_idx": var_idx,
                "g": lines(ev.value(gn), " " * 6),
                "dg_dx": [(j, lines(b, " " * 10)) for j, b in ev.cases_x(gn)],
                "dg_dp": [(k, lines(b, " " * 10)) for k, b in ev.cases_p(gn)],
                "dg_dt": lines(ev.partial_t(gn), " " * 6),
            }
            if _valid(time_raw):
                tn = model.parse(str(time_raw))
                item["t"] = lines(ev.value(tn), " " * 6)
                item["dt_dp"] = [(k, lines(b, " " * 10)) for k, b in ev.cases_p(tn)]
                time_events.append(item)
                continue
            rn = model.parse(str(root_raw))
            direction_raw = _get("direction", i)
            try:
                direction = int(direction_raw) if _valid(direction_raw) else 0
            except (ValueError, TypeError):
                direction = 0
            if direction not in (-1, 0, 1):
                raise ValueError(
                    f"Event {i}: direction must be -1, 0, or 1, got {direction}")
            terminal = str(_get("terminal", i)).lower() == "true"
            item.update({
                "r": rn,
                "dr_dx": [(j, lines(b, " " * 10)) for j, b in ev.cases_x(rn)],
                "dr_dp": [(k, lines(b, " " * 10)) for k, b in ev.cases_p(rn)],
                "dr_dt": lines(ev.partial_t(rn), " " * 6),
                "direction": direction,
                "terminal": terminal,
            })
            root_events.append(item)

    # The state switches follow the root events of the table: no reset, no
    # limit on how often they fire, and the mode they set (see
    # cppde_model.hold_switches).
    n_table_roots = len(root_events)
    for k, (rn, closed) in enumerate(model.switches):
        root_events.append({
            "var_idx": -1, "mode": k, "closed": closed, "r": rn,
            "r_val": lines(ev.value(rn), " " * 6),
            "dr_dx": [(j, lines(b, " " * 10)) for j, b in ev.cases_x(rn)],
            "dr_dp": [(q, lines(b, " " * 10)) for q, b in ev.cases_p(rn)],
            "dr_dt": lines(ev.partial_t(rn), " " * 6),
            "direction": 0, "terminal": False,
        })

    # root_fn: user roots, then event roots; an exhausted event root reads +1.
    root_stores = [(("vec", "gout", k), n, "=") for k, n in enumerate(rootfunc_nodes)]
    n_user = len(rootfunc_nodes)
    for j, e in enumerate(root_events):
        if "mode" in e:
            root_stores.append((("vec", "gout", n_user + j), e["r"], "="))
            continue
        root_stores.append((("call", {
            "cpp": "gout[{0}] = (ud->root_fired[{1}] >= ud->maxroot) ? 1.0 : (%s);",
            "py": "gout[{0}] = 1.0 if ud.root_fired[{1}] >= ud.maxroot else (%s)"},
            (n_user + j, j)), e["r"], "="))
    root_body = ""
    if root_stores:
        root_body = lines(cppde_model.stores_prelude(model, root_stores, "double")
                          + [cppde_model.block(g, root_stores, "_r")], "  ")

    # --- adjoint ---
    adj_rhs_body = adj_quad_body = ""
    if reverse:
        xb, qb = cppde_model.cvode_adjoint_statements(model)
        adj_rhs_body, adj_quad_body = lines(xb), lines(qb)

    compile_defs = ["-DCVODE_KLU"] if use_sparse else []
    data_code = "\n".join(model.linmap.cpp()) if model.linmap is not None else ""

    src = _render_source(
        modelname=modelname, version=version,
        n_states=n_states, n_global=n_global,
        ode_body=ode_body,
        jac_body=jac_body,
        nnz_total=nnz_total,
        sens_body=sens_body,
        reverse=reverse,
        adj_rhs_body=adj_rhs_body,
        adj_quad_body=adj_quad_body,
        asa_checkpoints=int(asa_checkpoints),
        cv_method=cv_method,
        deriv=deriv,
        use_sparse=use_sparse,
        use_lapack=use_lapack,
        klu_settings=klu_settings,
        n_forcings=n_forcings,
        rootfunc_mode=rootfunc_mode,
        n_user_rootfunc=n_user,
        root_body=root_body,
        time_events=time_events,
        root_events=root_events,
        n_table_roots=n_table_roots,
        n_modes=len(model.switches),
        n_params=n_params,
        states_list=states_list,
        params_list=params_list,
        forcings_list=forcings_list,
        include_time_zero=include_time_zero,
        data_code=data_code,
    )

    srcfile = os.path.join(outdir, f"{modelname}.cpp")
    if os.path.exists(srcfile):
        print(f"Overwriting existing file: {srcfile}")
    with open(srcfile, "w", encoding="utf-8") as f:
        f.write(src)

    return {
        "srcfile": srcfile,
        "variables": states_list,
        "parameters": params_list,
        "forcings": forcings_list,
        "sens_names": sens_names,
        "jac_nnz_rows": jac_rows,
        "jac_nnz_cols": jac_cols,
        "use_sparse": use_sparse,
        "compile_defs": compile_defs,
        "use_lapack": use_lapack,
        "codegen_stats": {"strategy": strategy, "lin_rows": model.nlin,
                          "graph_nodes": len(g),
                          "seconds": time.perf_counter() - t_start},
    }


# =====================================================================
# C++ source template
# =====================================================================

def _fail_main(code, msg):
    """Failure statement of the entry point's setup."""
    return '{ cleanup(); return res.fail(cppde::%s, "%s"); }' % (code, msg)


def _fail_seg(code, msg):
    """Failure statement of the ASA restart, which returns -1 to its caller."""
    return '{ asa_msg = "%s"; return -1; }' % msg


def _switch_lambda(cases, var):
    """Lambda (x, t, var) -> double over index cases [(index, body)]."""
    out = ["[params, &F](const double* x, double t, int %s) -> double {" % var,
           "      (void)x; (void)t;"]
    if cases:
        out.append("      switch (%s) {" % var)
        for idx, body in cases:
            out += ["        case %d: {" % idx, body, "        }"]
        out += ["        default: return 0.0;", "      }"]
    else:
        out.append("      (void)%s; return 0.0;" % var)
    out.append("    }")
    return "\n".join(out)


def _render_source(
    modelname, version,
    n_states, n_global,
    ode_body, jac_body, nnz_total,
    sens_body,
    cv_method, deriv, use_sparse,
    reverse=False, adj_rhs_body="", adj_quad_body="", asa_checkpoints=200,
    use_lapack=False,
    klu_settings=None,
    n_forcings=0,
    rootfunc_mode="none",
    n_user_rootfunc=0,
    root_body="",
    time_events=None,
    root_events=None,
    n_table_roots=None,
    n_modes=0,
    n_params=0,
    states_list=None,
    params_list=None,
    forcings_list=None,
    include_time_zero=True,
    data_code="",
):
    deriv_flag = "true" if deriv else "false"
    # The modes of the state switches live as long as the solve and its sweep.
    mode_scope = ((
        "  signed char _switch_mode[%d] = {0};\n"
        "  cppde::switch_mode_scope _cppde_mode_scope(_switch_mode, %d);\n")
        % (n_modes, n_modes)) if n_modes else ""
    linmap_include = "#include <cppde/cppde_linmap.hpp>\n" if data_code else ""
    data_block = "\n" + data_code + "\n" if data_code else ""
    has_forcings = n_forcings > 0
    if time_events is None:
        time_events = []
    if root_events is None:
        root_events = []
    if n_table_roots is None:
        n_table_roots = len(root_events)
    has_time_events = len(time_events) > 0
    has_root_events = len(root_events) > 0
    has_events      = has_time_events or has_root_events
    # ASA over events: one CVODES memory per stretch between jumps.
    asa_seg = reverse and has_events

    def reinit(t):
        if asa_seg:
            return f"asa_restart({t}) < 0"
        return f"CVodeReInit(cvode_mem, {t}, y) < 0"

    # UserData holds the PchipForcing storage and F whenever forcings or event
    # lambdas read it; the event lambdas capture F even without forcings.
    need_pchip = has_forcings
    if has_events and not has_forcings:
        need_pchip = True

    if need_pchip:
        forcing_include = "#include <cppde/cppde_pchip_forcing.hpp>\n"
        forcing_ud_members = (
            "  std::vector<cppde::PchipForcing<double>> forcing_storage;\n"
            "  std::vector<const cppde::PchipForcing<double>*> F;\n"
        )
        forcing_local = "  const auto& F = ud->F;\n  (void)F;"
    else:
        forcing_include = ""
        forcing_ud_members = ""
        forcing_local = ""

    # Root-event fire counters live in UserData, where root_fn reads them and
    # the main loop increments them.
    events_ud_members = (
        "  std::vector<int> root_fired;           // per event-root fire count\n"
        "  int              maxroot = 1;          // cap from solveODE(maxroot=)\n"
    )

    # UserData holds Phi_prime (the auto-extended Phi'(theta) flat,
    # column-major) so the sens rhs and event saltation lambdas can apply
    # the chain rule. R always supplies a full-shape sens1ini at solve time.
    if deriv:
        reparam_ud_members = (
            "  std::vector<double> Phi_prime;         // (NEQ + n_params) * Ns_active, column-major\n"
        )
    else:
        reparam_ud_members = ""

    if has_forcings:
        forcing_init_block = f"""  // --- Initialize forcings (PCHIP interpolation) ---
  {{
    const int n_f_in = static_cast<int>(args.flen.size());
    if (n_f_in != {n_forcings}) {{
      char _m[160];
      snprintf(_m, sizeof(_m), "Forcing count mismatch: model expects %d, got %d", {n_forcings}, n_f_in);
      return res.fail(cppde::RC_ILL_INPUT, _m);
    }}
    ud.forcing_storage.resize({n_forcings});
    ud.F.resize({n_forcings});
    for (int fi = 0; fi < {n_forcings}; ++fi) {{
      const int np = args.flen[fi];
      std::vector<double> ft(args.ftimes[fi],  args.ftimes[fi]  + np);
      std::vector<double> fv(args.fvalues[fi], args.fvalues[fi] + np);
      ud.forcing_storage[fi].initialize(ft, fv);
      ud.F[fi] = &ud.forcing_storage[fi];
    }}
  }}
"""
    else:
        forcing_init_block = ""

    # --- Rootfunc / event-root registration ---

    # root_fn writes the user rootfunc into gout[0..n_user-1] and the event roots
    # after it. On CV_ROOT_RETURN a user root terminates, an event root applies its
    # state change and saltation, then reinits and continues unless terminal.
    n_event_roots   = len(root_events)
    n_total_roots   = n_user_rootfunc + n_event_roots
    has_user_root   = rootfunc_mode == "user"
    need_cvode_root = (n_total_roots > 0)

    if need_cvode_root:
        root_gout_body = root_body or "  (void)gout;"
        rootfunc_decl = ("static int root_fn(sunrealtype t, N_Vector y, "
                         "sunrealtype* gout, void* ud_vp);")
        rootfunc_impl = f"""
static int root_fn(sunrealtype t, N_Vector y, sunrealtype* gout, void* ud_vp) {{
  (void)t;
  // An exception must not unwind through SUNDIALS' C frames (UB); these
  // bodies allocate, so bad_alloc is reachable.
  try {{
  UserData* ud = static_cast<UserData*>(ud_vp);
  const double* params = ud->params.data();
  const double* x = N_VGetArrayPointer(y);
  (void)x; (void)params;
{forcing_local}
{root_gout_body}
  return 0;
  }} catch (...) {{
    return -1;   // unrecoverable
  }}
}}
"""
        # Build direction array for CVodeSetRootDirection.  User roots
        # get direction 0 (any crossing); event roots use their spec.
        direction_entries = ["0"] * n_user_rootfunc + [
            str(e["direction"]) for e in root_events
        ]
        direction_list = ", ".join(direction_entries)

        # A state switch whose new branch holds the state on its surface leaves
        # its root at zero, which is no reason for a warning.
        quiet = "    CVodeSetNoInactiveRootWarn(cvode_mem);\n" if n_modes else ""

        def rootfunc_init(fail):
            return (
                f"  {{ int rd[{n_total_roots}] = {{ {direction_list} }};\n"
                f"    if (CVodeRootInit(cvode_mem, {n_total_roots}, root_fn) < 0) "
                + fail("RC_LINIT_FAIL", "CVodeRootInit failed") + "\n"
                f"    if (CVodeSetRootDirection(cvode_mem, rd) < 0) "
                + fail("RC_LINIT_FAIL", "CVodeSetRootDirection failed") + "\n"
                + quiet + "  }\n"
            )
    else:
        rootfunc_decl = ""
        rootfunc_impl = ""

        def rootfunc_init(fail):
            return ""
    rootfunc_init_block = rootfunc_init(_fail_main)

    # --- Events (time and root) ---

    # Time events land in a sorted std::vector<TimeEvent> traversed in the main
    # loop, root events in std::vector<RootEvent> dispatched from CV_ROOT_RETURN.
    event_block_includes = "#include <functional>\n" if has_events else ""
    if n_modes:
        event_block_includes += "#include <cppde/cppde_switch_modes.hpp>\n"
    if has_time_events:
        event_block_includes += "#include <cppde/cppde_event_window.hpp>\n"
    event_struct = ""
    if has_time_events:
        event_struct += f"""
struct TimeEvent {{
  double time;
  int    var_idx;
  // All lambdas capture `params` (pointer) and `F` (ref) from the builder.
  std::function<double(const double* x, double t)> g_fn;
  // dg/dx_i as a single dispatch lambda: returns 0 outside range.
  std::function<double(const double* x, double t, int i)> dg_dx_fn;
  // dg/dp_k (explicit partial, at fixed t_e).
  std::function<double(const double* x, double t, int k)> dg_dp_fn;
  // dg/dt (explicit time partial of g, at fixed x/p).
  std::function<double(const double* x, double t)> dg_dt_fn;
  // dt_e/dp_k (saltation correction for parameterized event times).
  std::function<double(const double* x, double t, int k)> dt_dp_fn;
}};

static std::vector<TimeEvent> build_time_events(const double* params,
                                                const std::vector<const cppde::PchipForcing<double>*>& F) {{
  (void)params; (void)F;
  std::vector<TimeEvent> ev;
"""
        ev_body = []
        for i, e in enumerate(time_events):
            ev_body.append(f"""  {{
    TimeEvent e;
    e.time    = [&]() -> double {{
{e['t']}
    }}();
    e.var_idx = {e['var_idx']};
    e.g_fn    = [params, &F](const double* x, double t) -> double {{
      (void)x; (void)t;
{e['g']}
    }};
    e.dg_dx_fn = {_switch_lambda(e['dg_dx'], 'i')};
    e.dg_dp_fn = {_switch_lambda(e['dg_dp'], 'k')};
    e.dg_dt_fn = [params, &F](const double* x, double t) -> double {{
      (void)x; (void)t;
{e['dg_dt']}
    }};
    e.dt_dp_fn = {_switch_lambda(e['dt_dp'], 'k')};
    ev.push_back(std::move(e));
  }}""")
        event_struct += "\n".join(ev_body) + """
  std::sort(ev.begin(), ev.end(),
            [](const TimeEvent& a, const TimeEvent& b) { return a.time < b.time; });
  return ev;
}
"""

    if has_root_events:
        event_struct += f"""
struct RootEvent {{
  int  var_idx;
  int  direction;
  bool terminal;
  // Event map g(x, t, p) and its partials (for saltation)
  std::function<double(const double* x, double t)> g_fn;
  std::function<double(const double* x, double t, int i)> dg_dx_fn;
  std::function<double(const double* x, double t, int k)> dg_dp_fn;
  std::function<double(const double* x, double t)> dg_dt_fn;
  // Root condition r(x, t, p) and its partials (for the IFT of dt_root/dp)
  std::function<double(const double* x, double t, int i)> dr_dx_fn;
  std::function<double(const double* x, double t, int k)> dr_dp_fn;
  std::function<double(const double* x, double t)> dr_dt_fn;
  // A state switch: the mode its root sets, which holds where r > 0 (r >= 0
  // when closed); var_idx is -1 and the event map unset. -1 for an event.
  int  mode = -1;
  bool closed = false;
  std::function<double(const double* x, double t)> r_fn;
}};

static std::vector<RootEvent> build_root_events(const double* params,
                                                const std::vector<const cppde::PchipForcing<double>*>& F) {{
  (void)params; (void)F;
  std::vector<RootEvent> ev;
"""
        re_body = []
        for i, e in enumerate(root_events):
            terminal = "true" if e["terminal"] else "false"
            if "mode" in e:
                re_body.append(f"""  {{
    RootEvent e;
    e.var_idx   = -1;
    e.direction = 0;
    e.terminal  = false;
    e.mode      = {e['mode']};
    e.closed    = {'true' if e['closed'] else 'false'};
    e.r_fn    = [params, &F](const double* x, double t) -> double {{
      (void)x; (void)t;
{e['r_val']}
    }};
    e.dr_dx_fn = {_switch_lambda(e['dr_dx'], 'i')};
    e.dr_dp_fn = {_switch_lambda(e['dr_dp'], 'k')};
    e.dr_dt_fn = [params, &F](const double* x, double t) -> double {{
      (void)x; (void)t;
{e['dr_dt']}
    }};
    ev.push_back(std::move(e));
  }}""")
                continue
            re_body.append(f"""  {{
    RootEvent e;
    e.var_idx   = {e['var_idx']};
    e.direction = {e['direction']};
    e.terminal  = {terminal};
    e.g_fn    = [params, &F](const double* x, double t) -> double {{
      (void)x; (void)t;
{e['g']}
    }};
    e.dg_dx_fn = {_switch_lambda(e['dg_dx'], 'i')};
    e.dg_dp_fn = {_switch_lambda(e['dg_dp'], 'k')};
    e.dg_dt_fn = [params, &F](const double* x, double t) -> double {{
      (void)x; (void)t;
{e['dg_dt']}
    }};
    e.dr_dx_fn = {_switch_lambda(e['dr_dx'], 'i')};
    e.dr_dp_fn = {_switch_lambda(e['dr_dp'], 'k')};
    e.dr_dt_fn = [params, &F](const double* x, double t) -> double {{
      (void)x; (void)t;
{e['dr_dt']}
    }};
    ev.push_back(std::move(e));
  }}""")
        event_struct += "\n".join(re_body) + """
  return ev;
}
"""

    # ---- Conditional sens code snippets ----
    if deriv:
        phi_rows = n_states + n_params
        sens_init_block = f"""  if (Ns_active > 0) {{
    yS = N_VCloneVectorArray(Ns_active, y);
    Ns_alloc = Ns_active;
    if (!yS) {{ cleanup(); return res.fail(cppde::RC_NO_MALLOC, "N_VCloneVectorArray failed"); }}
    // sens1ini is the auto-extended Phi'(theta) of shape
    // [NEQ + n_params, Ns_active] (column-major); copy to UserData for the
    // sens RHS callback and event-saltation chain rule.
    const int phi_rows_rt = {phi_rows};
    const int expected_len = phi_rows_rt * Ns_active;
    if (args.n_sens1 != expected_len) {{
      char _m[160];
      snprintf(_m, sizeof(_m), "tangent length %d != expected %d (phi_rows * Ns_active)",
               args.n_sens1, expected_len);
      cleanup();
      return res.fail(cppde::RC_ILL_INPUT, _m);
    }}
    ud.Phi_prime.assign(args.sens1ini, args.sens1ini + expected_len);
    // yS0[iS][i] = Phi'[i, iS] for state-row i
    for (int iS = 0; iS < Ns_active; ++iS) {{
      double* v = N_VGetArrayPointer(yS[iS]);
      for (int i = 0; i < NEQ; ++i) v[i] = ud.Phi_prime[i + phi_rows_rt * iS];
    }}
    if (CVodeSensInit1(cvode_mem, Ns_active, CV_STAGGERED, sens_rhs1_fn, yS) < 0) {{
      cleanup(); return res.fail(cppde::RC_LINIT_FAIL, "CVodeSensInit1 failed");
    }}
    if (CVodeSensEEtolerances(cvode_mem) < 0) {{
      cleanup(); return res.fail(cppde::RC_LINIT_FAIL, "CVodeSensEEtolerances failed");
    }}
    if (CVodeSetSensErrCon(cvode_mem, SUNTRUE) < 0) {{
      cleanup(); return res.fail(cppde::RC_LINIT_FAIL, "CVodeSetSensErrCon failed");
    }}
  }}
"""
    else:
        sens_init_block = ""

    if deriv:
        sens_t0_block = """    for (int iS = 0; iS < Ns_active; ++iS) {
      double* v = N_VGetArrayPointer(yS[iS]);
      for (int i = 0; i < NEQ; ++i) out_s.push_back(v[i]);
    }"""
        sens_get_block = """    if (Ns_active > 0) {
      int flag_gs = CVodeGetSens(cvode_mem, &tret, yS);
      if (flag_gs < 0) {
        return_code = flag_gs;
        solver_msg = "CVodeGetSens failed";
        break;
      }
    }"""
        sens_store_block = """    for (int iS = 0; iS < Ns_active; ++iS) {
      double* v = N_VGetArrayPointer(yS[iS]);
      for (int i = 0; i < NEQ; ++i) out_s.push_back(v[i]);
    }"""
        sens_fevals_block = """  n_fe += cvc.nfeS;"""
    else:
        sens_t0_block = ""
        sens_get_block = ""
        sens_store_block = ""
        sens_fevals_block = ""

    # Event helpers emitted into the entry point. `apply_event` updates y[var_idx]
    # and, with deriv, yS[*][var_idx] from the precomputed partials dg/dx, dg/dp.

    # The saltation body is templated on the event handle: TimeEvent uses
    # ev.dt_dp_fn, RootEvent the IFT-computed dt_root_dp passed in.
    event_sens_reinit = ""
    time_event_apply_lambda = ""
    root_event_apply_lambda = ""

    if has_events:
        event_builder_block = (
            "  f_buf = N_VClone(y);\n"
            "  if (!f_buf) { cleanup(); return res.fail(cppde::RC_NO_MALLOC, \"N_VClone failed for event scratch\"); }\n"
        )
        if has_time_events:
            event_builder_block += (
                "  auto time_events = build_time_events(ud.params.data(), ud.F);\n"
            )
        if has_root_events:
            event_builder_block += (
                "  auto event_roots = build_root_events(ud.params.data(), ud.F);\n"
            )
        if n_modes:
            event_builder_block += (
                "  // The modes of the state switches: read off the state, or set by\n"
                "  // the crossing of their root, after which a reset among the events\n"
                "  // reads every other mode afresh.\n"
                "  auto cv_read_modes = [&](const double* x, double t) {\n"
                "    for (const auto& e : event_roots)\n"
                "      if (e.mode >= 0)\n"
                "        cppde::detail::set_switch_mode(\n"
                "            e.mode, cppde::detail::mode_holds(e.r_fn(x, t), e.closed));\n"
                "  };\n"
                "  auto cv_switch_after = [&](const std::vector<int>& trig,\n"
                "                             const std::vector<int>& dir,\n"
                "                             const double* x, double t) {\n"
                "    bool jumped = false;\n"
                "    for (size_t q = 0; q < trig.size(); ++q) {\n"
                "      const auto& e = event_roots[trig[q]];\n"
                "      if (e.mode >= 0) cppde::detail::set_switch_mode(e.mode, dir[q] > 0);\n"
                "      else if (!e.terminal) jumped = true;\n"
                "    }\n"
                "    if (!jumped) return;\n"
                "    for (size_t j = 0; j < event_roots.size(); ++j) {\n"
                "      const auto& e = event_roots[j];\n"
                "      if (e.mode < 0 ||\n"
                "          std::find(trig.begin(), trig.end(), (int)j) != trig.end()) continue;\n"
                "      cppde::detail::set_switch_mode(\n"
                "          e.mode, cppde::detail::mode_holds(e.r_fn(x, t), e.closed));\n"
                "    }\n"
                "  };\n"
                "  cv_read_modes(N_VGetArrayPointer(y), t0);\n"
                + ("  asa_segs[0].modes = cppde::detail::switch_mode_snapshot();\n"
                   if asa_seg else "")
            )
        if deriv:
            event_sens_reinit = ("      if (Ns_active > 0 && CVodeSensReInit(cvode_mem, CV_STAGGERED, yS) < 0) "
                                 "{ hard_fail = true; return_code = cppde::RC_LSETUP_FAIL; "
                                 "solver_msg = \"CVodeSensReInit failed\"; }\n")
    else:
        event_builder_block = ""

    if has_time_events:
        # First-order saltation for a parameterised time event (Barton & Lee), with
        # Gdot = sum_i dg/dx_i f_old[i] + dg/dt - f_new[var]:
        #   S_new[var] = sum_i dg/dx_i S_old[i] + dg/dp_k + Gdot dt_e/dp_k

        # Every other state shifts by (f_old[j] - f_new[j]) dt_e/dp_k. iS indexes theta
        # slots, and the chain rule dt_e/dtheta = sum_k dt_e/dp_k M[NEQ+k, iS] is
        # applied at the call site, with M from ud.Phi_prime.
        if deriv:
            phi_rows_e = n_states + n_params
            time_event_apply_lambda = f"""  auto apply_time_event = [&](const TimeEvent& ev) {{
    double* y_arr = N_VGetArrayPointer(y);
    double t_e = ev.time;
    std::vector<double> x_old(NEQ);
    for (int i = 0; i < NEQ; ++i) x_old[i] = y_arr[i];

    rhs_fn(t_e, y, f_buf, &ud);
    std::vector<double> f_old(NEQ);
    {{ const double* fb = N_VGetArrayPointer(f_buf);
      for (int i = 0; i < NEQ; ++i) f_old[i] = fb[i]; }}

    double new_v = ev.g_fn(x_old.data(), t_e);
    y_arr[ev.var_idx] = new_v;
@READ_MODES@
    rhs_fn(t_e, y, f_buf, &ud);
    std::vector<double> f_new(NEQ);
    {{ const double* fb = N_VGetArrayPointer(f_buf);
      for (int i = 0; i < NEQ; ++i) f_new[i] = fb[i]; }}

    const double dg_dt = ev.dg_dt_fn(x_old.data(), t_e);
    const int phi_rows_rt = {phi_rows_e};

    for (int iS = 0; iS < Ns_active; ++iS) {{
      double* yS_k = N_VGetArrayPointer(yS[iS]);
      std::vector<double> S_old(NEQ);
      for (int i = 0; i < NEQ; ++i) S_old[i] = yS_k[i];

      // Chain-rule over model-param slots via Phi'(theta) param block.
      double dt_dp = 0.0;
      double dg_dp_k = 0.0;
      for (int pk = 0; pk < {n_params}; ++pk) {{
        double M = ud.Phi_prime[(NEQ + pk) + phi_rows_rt * iS];
        if (M == 0.0) continue;
        dt_dp   += ev.dt_dp_fn(x_old.data(), t_e, pk) * M;
        dg_dp_k += ev.dg_dp_fn(x_old.data(), t_e, pk) * M;
      }}

      double sum_gx_S = 0.0, sum_gx_f = 0.0;
      for (int i = 0; i < NEQ; ++i) {{
        double gx_i = ev.dg_dx_fn(x_old.data(), t_e, i);
        sum_gx_S += gx_i * S_old[i];
        sum_gx_f += gx_i * f_old[i];
      }}

      yS_k[ev.var_idx] = sum_gx_S + dg_dp_k
                       + (sum_gx_f + dg_dt - f_new[ev.var_idx]) * dt_dp;

      if (dt_dp != 0.0) {{
        for (int j = 0; j < NEQ; ++j) {{
          if (j == ev.var_idx) continue;
          yS_k[j] = S_old[j] + (f_old[j] - f_new[j]) * dt_dp;
        }}
      }}
    }}
  }};
"""
        elif reverse:
            # Under ASA the jump is applied as without derivatives and recorded
            # for the backward sweep, with f on both sides of it.
            time_event_apply_lambda = """  auto apply_time_event = [&](const TimeEvent& ev) {
    double* y_arr = N_VGetArrayPointer(y);
    AsaJump J;
    J.kind = 0;
    J.ev = (int)(&ev - time_events.data());
    J.t = ev.time;
    J.x_old.assign(y_arr, y_arr + NEQ);
    rhs_fn(J.t, y, f_buf, &ud);
    { const double* fb = N_VGetArrayPointer(f_buf); J.f_old.assign(fb, fb + NEQ); }
    y_arr[ev.var_idx] = ev.g_fn(J.x_old.data(), J.t);
@READ_MODES_J@    rhs_fn(J.t, y, f_buf, &ud);
    { const double* fb = N_VGetArrayPointer(f_buf); J.f_new.assign(fb, fb + NEQ); }
    asa_jumps.push_back(std::move(J));
  };
"""
        else:
            time_event_apply_lambda = """  auto apply_time_event = [&](const TimeEvent& ev) {
    double* y_arr = N_VGetArrayPointer(y);
    double t_e = ev.time;
    std::vector<double> x_old(NEQ);
    for (int i = 0; i < NEQ; ++i) x_old[i] = y_arr[i];
    y_arr[ev.var_idx] = ev.g_fn(x_old.data(), t_e);
@READ_MODES@
  };
"""
        time_event_apply_lambda = (
            time_event_apply_lambda
            .replace("@READ_MODES@\n", "    cv_read_modes(y_arr, t_e);\n" if n_modes else "")
            .replace("@READ_MODES_J@", "    cv_read_modes(y_arr, J.t);\n" if n_modes else ""))
        # A time event fires inside the grid's window only, the rule of
        # cppde_event_window.hpp that the native backend applies: one before t0
        # is skipped, one at t0 is applied before the solve starts, and the t0
        # row is written again with the state after it, as every row at an
        # event time has. The main loop stops short of the last time.
        # Under ASA the t0 row written below belongs to the stretch after them.
        asa_t0_rows = ("      for (auto& s_ : asa_segs) s_.row0 = s_.row1 = 0;\n"
                       if asa_seg else "")
        event_pre_t0_block = f"""  const double t_last = times[n_times - 1];
  size_t ev_idx = 0;
  while (ev_idx < time_events.size() && time_events[ev_idx].time < t0 &&
         !cppde::detail::same_event_time(time_events[ev_idx].time, t0)) ev_idx++;
  {{
    size_t applied_pre = 0;
    while (ev_idx < time_events.size() &&
           cppde::detail::same_event_time(time_events[ev_idx].time, t0)) {{
      apply_time_event(time_events[ev_idx]);
      ev_idx++;
      applied_pre++;
    }}
    if (applied_pre > 0) {{
      cv_harvest();
      if ({reinit("t0")})
        {{ cleanup(); return res.fail(cppde::RC_LINIT_FAIL, "CVodeReInit failed (events at t0)"); }}
{asa_t0_rows}{event_sens_reinit}      cv_rebase();
      out_y.clear();
      out_s.clear();
      {{
        const double* y0 = N_VGetArrayPointer(y);
        for (int i = 0; i < NEQ; ++i) out_y.push_back(y0[i]);
{sens_t0_block}
      }}
    }}
  }}
"""
    else:
        event_pre_t0_block = ""

    if has_root_events:
        # Batched root-event application: every event at the same t_e resolves against
        # one pre-event snapshot and shares one dt/dp per sens slot, so the flow
        # correction (f_old - f_new) dt_dp lands at most once per state.

        # dt_dp comes from the first non-terminal triggered event: simultaneous roots
        # cross at the same t by construction, so their dt agrees to first order.
        if deriv:
            phi_rows_r = n_states + n_params
            root_event_apply_lambda = f"""  auto apply_root_events_batch = [&](const std::vector<int>& triggered_idx,
                                     const std::vector<int>& trig_dir,
                                     double t_e) -> bool {{
    (void)trig_dir;
    double* y_arr = N_VGetArrayPointer(y);
    std::vector<double> x_old(NEQ);
    for (int i = 0; i < NEQ; ++i) x_old[i] = y_arr[i];

    rhs_fn(t_e, y, f_buf, &ud);
    std::vector<double> f_old(NEQ);
    {{ const double* fb = N_VGetArrayPointer(f_buf);
      for (int i = 0; i < NEQ; ++i) f_old[i] = fb[i]; }}

    bool any_terminal = false;
    std::vector<char> is_modified(NEQ, 0);
    for (int j : triggered_idx) {{
      const auto& ev = event_roots[j];
      if (ev.terminal) {{ any_terminal = true; continue; }}
      if (ev.var_idx < 0) continue;
      y_arr[ev.var_idx] = ev.g_fn(x_old.data(), t_e);
      is_modified[ev.var_idx] = 1;
    }}
@SWITCH_AFTER@
    rhs_fn(t_e, y, f_buf, &ud);
    std::vector<double> f_new(NEQ);
    {{ const double* fb = N_VGetArrayPointer(f_buf);
      for (int i = 0; i < NEQ; ++i) f_new[i] = fb[i]; }}

    int ref_j = -1;
    for (int j : triggered_idx) {{
      if (!event_roots[j].terminal) {{ ref_j = j; break; }}
    }}
    if (ref_j < 0) return any_terminal;
    const auto& ref = event_roots[ref_j];

    double g_dot = ref.dr_dt_fn(x_old.data(), t_e);
    std::vector<double> rx(NEQ);
    for (int i = 0; i < NEQ; ++i) {{
      rx[i] = ref.dr_dx_fn(x_old.data(), t_e, i);
      g_dot += rx[i] * f_old[i];
    }}
    if (std::fabs(g_dot) < 1e-14) {{
      Rf_warning("Root event: tangential crossing (g_dot ~ 0): sensitivities may be unreliable at t=%.6e", t_e);
    }}

    const int phi_rows_rt = {phi_rows_r};

    for (int iS = 0; iS < Ns_active; ++iS) {{
      double* yS_k = N_VGetArrayPointer(yS[iS]);
      std::vector<double> S_old(NEQ);
      for (int i = 0; i < NEQ; ++i) S_old[i] = yS_k[i];

      // dr_num = Σ rx_i · S_old_i + Σ_pk dr/dp_pk · M[NEQ+pk, iS]
      double dr_num = 0.0;
      for (int i = 0; i < NEQ; ++i) dr_num += rx[i] * S_old[i];
      for (int pk = 0; pk < {n_params}; ++pk) {{
        double M = ud.Phi_prime[(NEQ + pk) + phi_rows_rt * iS];
        if (M == 0.0) continue;
        dr_num += ref.dr_dp_fn(x_old.data(), t_e, pk) * M;
      }}
      double dt_dp = (g_dot != 0.0) ? -dr_num / g_dot : 0.0;

      if (dt_dp != 0.0) {{
        for (int i = 0; i < NEQ; ++i) {{
          if (is_modified[i]) continue;
          yS_k[i] = S_old[i] + (f_old[i] - f_new[i]) * dt_dp;
        }}
      }}

      for (int j : triggered_idx) {{
        const auto& ev = event_roots[j];
        if (ev.terminal || ev.var_idx < 0) continue;
        int v = ev.var_idx;

        double sum_gx_S = 0.0, sum_gx_f = 0.0;
        for (int i = 0; i < NEQ; ++i) {{
          double gx_i = ev.dg_dx_fn(x_old.data(), t_e, i);
          sum_gx_S += gx_i * S_old[i];
          sum_gx_f += gx_i * f_old[i];
        }}
        double dg_dt = ev.dg_dt_fn(x_old.data(), t_e);
        double dg_dp_k = 0.0;
        for (int pk = 0; pk < {n_params}; ++pk) {{
          double M = ud.Phi_prime[(NEQ + pk) + phi_rows_rt * iS];
          if (M == 0.0) continue;
          dg_dp_k += ev.dg_dp_fn(x_old.data(), t_e, pk) * M;
        }}

        yS_k[v] = sum_gx_S + dg_dp_k
                + (sum_gx_f + dg_dt - f_new[v]) * dt_dp;
      }}
    }}

    return any_terminal;
  }};
"""
        elif reverse:
            root_event_apply_lambda = """  auto apply_root_events_batch = [&](const std::vector<int>& triggered_idx,
                                     const std::vector<int>& trig_dir,
                                     double t_e) -> bool {
    (void)trig_dir;
    double* y_arr = N_VGetArrayPointer(y);
    AsaJump J;
    J.kind = 1;
    J.trig = triggered_idx;
    J.t = t_e;
    J.x_old.assign(y_arr, y_arr + NEQ);
    rhs_fn(t_e, y, f_buf, &ud);
    { const double* fb = N_VGetArrayPointer(f_buf); J.f_old.assign(fb, fb + NEQ); }
    bool any_terminal = false;
    for (int j : triggered_idx) {
      const auto& ev = event_roots[j];
      if (ev.terminal) { any_terminal = true; continue; }
      if (ev.var_idx < 0) continue;
      y_arr[ev.var_idx] = ev.g_fn(J.x_old.data(), t_e);
    }
@SWITCH_AFTER@
    rhs_fn(t_e, y, f_buf, &ud);
    { const double* fb = N_VGetArrayPointer(f_buf); J.f_new.assign(fb, fb + NEQ); }
    asa_jumps.push_back(std::move(J));
    return any_terminal;
  };
"""
        else:
            root_event_apply_lambda = """  auto apply_root_events_batch = [&](const std::vector<int>& triggered_idx,
                                     const std::vector<int>& trig_dir,
                                     double t_e) -> bool {
    (void)trig_dir;
    double* y_arr = N_VGetArrayPointer(y);
    std::vector<double> x_old(NEQ);
    for (int i = 0; i < NEQ; ++i) x_old[i] = y_arr[i];
    bool any_terminal = false;
    for (int j : triggered_idx) {
      const auto& ev = event_roots[j];
      if (ev.terminal) { any_terminal = true; continue; }
      if (ev.var_idx < 0) continue;
      y_arr[ev.var_idx] = ev.g_fn(x_old.data(), t_e);
    }
@SWITCH_AFTER@
    return any_terminal;
  };
"""
        root_event_apply_lambda = root_event_apply_lambda.replace(
            "@SWITCH_AFTER@\n",
            "    cv_switch_after(triggered_idx, trig_dir, y_arr, t_e);\n" if n_modes else "")

    event_apply_lambda = time_event_apply_lambda + root_event_apply_lambda

    if rootfunc_mode == "equilibrate":
        # Post-step check: every |ydot| <= root_tol |y| + abstol.
        equilibrate_check_block = """    {
      N_Vector ydot_tmp = N_VClone(y);
      rhs_fn(t_reached, y, ydot_tmp, &ud);
      const double* yd = N_VGetArrayPointer(ydot_tmp);
      const double* yv = N_VGetArrayPointer(y);
      bool steady = true;
      for (int i = 0; i < NEQ && steady; ++i)
        steady = std::fabs(yd[i]) <= root_tol * std::fabs(yv[i]) + abstol;
      N_VDestroy(ydot_tmp);
      if (steady) {
        solver_msg = "Terminated: steady state reached (equilibrate)";
        root_terminated = true;
        break;
      }
    }
"""
    else:
        equilibrate_check_block = ""

    # --- do_cvode_step: single step with CV_ROOT_RETURN dispatch ---

    # Returns 0 for the target reached, 1 for a user rootfunc stop, 2 for a terminal
    # root event, -1 for an error, with return_code holding the raw CVODE flag.
    # On 1 and 2 the caller pushes the final output row.
    if has_events:
        if deriv:
            sens_get_lambda = """      {
        int flag_gs = CVodeGetSens(cvode_mem, &tret, yS);
        if (flag_gs < 0) {
          return_code = flag_gs;
          solver_msg = "CVodeGetSens failed";
          return -1;
        }
      }
"""
        else:
            sens_get_lambda = ""

        # Root-dispatch body: only meaningful when CVodeRootInit was called.
        if need_cvode_root:
            user_check = ""
            if n_user_rootfunc > 0:
                user_check = f"""      for (int i = 0; i < {n_user_rootfunc}; ++i) {{
        if (rinfo[i] != 0) {{
          solver_msg = "Terminated: rootfunc crossed zero";
          root_terminated = true;
          return 1;
        }}
      }}
"""
            event_apply_sec = ""
            if has_root_events:
                event_apply_sec = f"""      std::vector<int> triggered_idx, trig_dir;
      bool rows = false;   // a state switch alone writes no row
      for (int j = 0; j < {n_event_roots}; ++j) {{
        if (rinfo[{n_user_rootfunc} + j] != 0 &&
            (ud.root_fired[j] < ud.maxroot || event_roots[j].mode >= 0)) {{
          triggered_idx.push_back(j);
          trig_dir.push_back(rinfo[{n_user_rootfunc} + j]);
          rows = rows || event_roots[j].mode < 0;
        }}
      }}
      bool any_terminal = false;
      if (!triggered_idx.empty()) {{
        // The native backend emits the state just before the event and the
        // state just after it, both at the root time. Skip a row that the
        // requested-times loop is about to write anyway.
        if (rows && std::abs(target - ((double)tret - 1e-15)) >= 1e-14) {{
          out_t.push_back((double)tret - 1e-15);
          {{ const double* y_arr = N_VGetArrayPointer(y);
            for (int i = 0; i < NEQ; ++i) out_y.push_back(y_arr[i]); }}
{sens_store_block}
        }}
        any_terminal = apply_root_events_batch(triggered_idx, trig_dir, (double)tret);
        for (int j : triggered_idx) ud.root_fired[j]++;
#ifdef CVODE_STEP_TRACE
        {{
          N_Vector _ele_r = N_VClone(y);
          N_Vector _ewt_r = N_VClone(y);
          cvode_emit_trace_row(cvode_mem, _ele_r, _ewt_r, nullptr, nullptr,
                               Ns_active, (double)tret, "CVODE_event",
                               cvc.nst - cvc.b_nst, cvc.nfe - cvc.b_nfe,
                               cvc.nje - cvc.b_nje, cvc.nsetups - cvc.b_nsetups);
          N_VDestroy(_ele_r); N_VDestroy(_ewt_r);
        }}
#endif
        cv_harvest();
        if ({reinit("tret")}) {{
          return_code = cppde::RC_ILL_INPUT;
          solver_msg = "CVodeReInit after root event failed";
          return -1;   // no cleanup: the caller still reads y
        }}
{event_sens_reinit}        if (hard_fail) return -1;
        cv_rebase();
        // A terminal event stops here and the caller writes the closing row.
        if (rows && !any_terminal && std::abs(target - (double)tret) >= 1e-14) {{
          out_t.push_back((double)tret);
          {{ const double* y_arr = N_VGetArrayPointer(y);
            for (int i = 0; i < NEQ; ++i) out_y.push_back(y_arr[i]); }}
{sens_store_block}
        }}
      }}
      if (any_terminal) {{
        solver_msg = "Terminated: terminal root event fired";
        root_terminated = true;
        return 2;
      }}
"""
            root_dispatch_block = f"""      std::vector<int> rinfo({n_total_roots});
      CVodeGetRootInfo(cvode_mem, rinfo.data());
{user_check}{event_apply_sec}"""
        else:
            root_dispatch_block = ("      solver_msg = \"unexpected CV_ROOT_RETURN\";\n"
                                   "      return_code = CV_UNRECOGNIZED_ERR; return -1;\n")

        if asa_seg:
            step_call = ("      int ncheck_ = 0;\n"
                         "      int flag = CVodeF(cvode_mem, target, y, &tret, CV_NORMAL, &ncheck_);\n"
                         "      asa_segs.back().used = true;\n")
        else:
            step_call = "      int flag = CVode(cvode_mem, target, y, &tret, CV_NORMAL);\n"
        do_cvode_step_lambda = f"""  auto do_cvode_step = [&](double target, double& tret_out) -> int {{
    while (true) {{
      sunrealtype tret;
{step_call}      if (flag < 0) {{
        return_code = flag;
        char buf[160];
        std::snprintf(buf, sizeof(buf),
                      "CVode failed at t=%.6e with flag %d", target, flag);
        solver_msg = buf;
        return -1;
      }}
{sens_get_lambda}      tret_out = (double)tret;
      if (flag != CV_ROOT_RETURN) return 0;
{root_dispatch_block}      if (tret_out >= target) return 0;
    }}
  }};
"""
    else:
        do_cvode_step_lambda = ""

    # --- Zero-copy sink: the batch entry sizes the results before the solve when
    # the grid is fixed (no root event, no rootfunc, no ASA adjoint). Time-event
    # rows are counted per condition through <model>_fixed_event_times.
    cv_fixed_grid = (n_table_roots == 0 and rootfunc_mode == "none"
                     and not reverse)
    n_cv_ev = len(time_events) if cv_fixed_grid else 0
    if cv_fixed_grid:
        _lines = "".join(
            f"  out[{i}] = [&]() -> double {{\n{e['t']}\n  }}();\n"
            for i, e in enumerate(time_events))
        cv_event_times_fn = (
            f"static void {modelname}_fixed_event_times(const double* params,\n"
            f"                                          double* out) {{\n"
            f"  (void)params; (void)out;\n{_lines}}}\n\n")
    else:
        cv_event_times_fn = ""

    cv_zero_block = (
        "  // ensure time zero is included, as the native backend does\n"
        "  if (std::find(times.begin(), times.end(), 0.0) == times.end())\n"
        "    times.push_back(0.0);\n") if include_time_zero else ""

    if cv_fixed_grid:
        _ev = ""
        if n_cv_ev:
            _ev = (f"  std::vector<std::vector<double> > EV(K, std::vector<double>({n_cv_ev}));\n"
                   f"  for (int k = 0; k < K; ++k)\n"
                   f"    {modelname}_fixed_event_times(A[k].params, EV[k].data());\n")
        cv_prealloc_block = (
            _ev +
            "  // CVODE needs sens1ini when deriv is on, so sens_width() reads the\n"
            "  // requested column count straight off the condition.\n"
            "  std::vector<cppde::rbatch::pre_ctx> P = cppde::rbatch::prealloc_batch(\n"
            f"      out, A, NEQ, 0, {'true' if deriv else 'false'}, false, "
            f"{'true' if include_time_zero else 'false'},\n"
            f"      dimnamesSEXP{', &EV' if n_cv_ev else ''});\n"
            "  for (int k = 0; k < K; ++k) {\n"
            "    R_[k].acquire = &cppde::rbatch::pre_acquire;\n"
            "    R_[k].acquire_ctx = &P[k];\n"
            "  }\n")
        cv_finish_block = (
            f"    if (R_[k].used_sink) {{ cppde::rbatch::finish_prealloc(out, k, R_[k], "
            f"{'true' if deriv else 'false'}, false); continue; }}\n")
    else:
        cv_prealloc_block = "  (void)dimnamesSEXP;   // root handling: the grid is dynamic\n"
        cv_finish_block = ""

    # --- Main integration loop body ---
    # Under ASA the forward pass runs CVodeF, which stores the checkpoints;
    # ncheck_ is not read.
    if reverse:
        cv_fwd_call = ("    int ncheck_ = 0;\n"
                       "    int flag = CVodeF(cvode_mem, times[k], y, &tret, "
                       "CV_NORMAL, &ncheck_);")
    else:
        cv_fwd_call = ("    int flag = CVode(cvode_mem, times[k], y, &tret, "
                       "CV_NORMAL);")

    if has_events:
        if has_time_events:
            time_interleave_block = f"""    // Time-event interleave: integrate to each event time, apply, reinit.
    // Events are sorted, so the first one at or past the last time ends it.
    while (ev_idx < time_events.size() && time_events[ev_idx].time <= times[k] &&
           cppde::detail::fixed_event_in_window(time_events[ev_idx].time, t0, t_last)) {{
      double t_e = time_events[ev_idx].time;
      if (t_e > t_reached) {{
        double tret_loc;
        int rc = do_cvode_step(t_e, tret_loc);
        if (rc < 0) {{ stop = true; break; }}
        t_reached = tret_loc;
        if (rc >= 1) {{
          out_t.push_back(t_reached);
          {{ const double* y_arr = N_VGetArrayPointer(y);
            for (int i = 0; i < NEQ; ++i) out_y.push_back(y_arr[i]); }}
{sens_store_block}
          stop = true; break;
        }}
      }}
      apply_time_event(time_events[ev_idx]);
#ifdef CVODE_STEP_TRACE
      cvode_emit_trace_row(cvode_mem, _ele_buf, _ewt_buf, nullptr, nullptr,
                           Ns_active, t_e, "CVODE_event",
                           cvc.nst - cvc.b_nst, cvc.nfe - cvc.b_nfe,
                           cvc.nje - cvc.b_nje, cvc.nsetups - cvc.b_nsetups);
#endif
      cv_harvest();
      if ({reinit("t_e")}) {{
        char _m[128];
        snprintf(_m, sizeof(_m), "CVodeReInit failed at t=%.6e", t_e);
        return_code = cppde::RC_ILL_INPUT; solver_msg = _m;
        stop = true; break;
      }}
{event_sens_reinit}      if (hard_fail) {{ stop = true; break; }}
      cv_rebase();
      // An event time that is not itself a requested output time still gets a
      // row, holding the post-event state, same contract as the native
      // backend. When it coincides with times[k] the branch below emits it.
      if (t_e < times[k]) {{
        out_t.push_back(t_e);
        {{ const double* y_arr = N_VGetArrayPointer(y);
          for (int i = 0; i < NEQ; ++i) out_y.push_back(y_arr[i]); }}
{sens_store_block}
      }}
      ev_idx++;
    }}
    if (stop) break;
"""
        else:
            time_interleave_block = ""

        main_loop_body = f"""  bool stop = false;
#ifdef CVODE_STEP_TRACE
  // Trace with events: CV_NORMAL, one row per output point and per applied
  // event; `nst` is cumulative, `h` and `q` are the last accepted step and order.
  N_Vector _ele_buf = N_VClone(y);
  N_Vector _ewt_buf = N_VClone(y);
#endif
  for (int k = 1; k < n_times && !stop; ++k) {{
{time_interleave_block}
    // Integrate to the next output time.
    if (times[k] > t_reached) {{
      double tret_loc;
      int rc = do_cvode_step(times[k], tret_loc);
      if (rc < 0) break;
      t_reached = tret_loc;
      out_t.push_back(t_reached);
      {{ const double* y_arr = N_VGetArrayPointer(y);
        for (int i = 0; i < NEQ; ++i) out_y.push_back(y_arr[i]); }}
{sens_store_block}
#ifdef CVODE_STEP_TRACE
      cvode_emit_trace_row(cvode_mem, _ele_buf, _ewt_buf, nullptr, nullptr,
                           Ns_active, t_reached, "CVODE",
                           cvc.nst - cvc.b_nst, cvc.nfe - cvc.b_nfe,
                           cvc.nje - cvc.b_nje, cvc.nsetups - cvc.b_nsetups);
#endif
      if (rc >= 1) {{ stop = true; break; }}
{equilibrate_check_block}
    }} else {{
      // Reached times[k] already via an event: record snapshot.
      out_t.push_back(t_reached);
      {{ const double* y_arr = N_VGetArrayPointer(y);
        for (int i = 0; i < NEQ; ++i) out_y.push_back(y_arr[i]); }}
{sens_store_block}
    }}
  }}
#ifdef CVODE_STEP_TRACE
  N_VDestroy(_ele_buf);
  N_VDestroy(_ewt_buf);
#endif
"""
    else:
        main_loop_body = f"""
#ifdef CVODE_STEP_TRACE
  {{
    // Trace without events: CV_ONE_STEP, one row per accepted step, outputs
    // interpolated by CVodeGetDky / CVodeGetSensDky; no rootfunc check.
    N_Vector _ele_buf    = N_VClone(y);
    N_Vector _ewt_buf    = N_VClone(y);
    N_Vector _y_interp   = N_VClone(y);
    N_Vector* _eleS_buf  = (Ns_active > 0) ? N_VCloneVectorArray(Ns_active, y) : nullptr;
    N_Vector* _ewtS_buf  = (Ns_active > 0) ? N_VCloneVectorArray(Ns_active, y) : nullptr;
    N_Vector* _yS_interp = (Ns_active > 0) ? N_VCloneVectorArray(Ns_active, y) : nullptr;

    const double t_final = (double)times[n_times - 1];
    CVodeSetStopTime(cvode_mem, t_final);

    int out_idx = 1;  // times[0] = t0 already written
    sunrealtype tret = t0;
    while (out_idx < n_times) {{
      int flag = CVode(cvode_mem, t_final, y, &tret, CV_ONE_STEP);
      if (flag < 0) {{
        return_code = flag;
        char buf[160];
        std::snprintf(buf, sizeof(buf),
                      "CVode failed at t=%.6e with flag %d", (double)tret, flag);
        solver_msg = buf;
        break;
      }}

      cvode_emit_trace_row(cvode_mem, _ele_buf, _ewt_buf,
                           _eleS_buf, _ewtS_buf, Ns_active, (double)tret);

      while (out_idx < n_times && (double)times[out_idx] <= (double)tret) {{
        {{
          int flag_dky = CVodeGetDky(cvode_mem, times[out_idx], 0, _y_interp);
          if (flag_dky != 0) {{
            return_code = flag_dky;
            solver_msg = "CVodeGetDky failed";
            out_idx = n_times;
            break;
          }}
        }}
        out_t.push_back((double)times[out_idx]);
        {{
          const double* y_arr = N_VGetArrayPointer(_y_interp);
          for (int i = 0; i < NEQ; ++i) out_y.push_back(y_arr[i]);
        }}
        if (Ns_active > 0 && _yS_interp != nullptr) {{
          int flag_sdky = CVodeGetSensDky(cvode_mem, times[out_idx], 0, _yS_interp);
          if (flag_sdky != 0) {{
            return_code = flag_sdky;
            solver_msg = "CVodeGetSensDky failed";
            out_idx = n_times;
            break;
          }}
          for (int j = 0; j < Ns_active; ++j) {{
            const double* yS_arr = N_VGetArrayPointer(_yS_interp[j]);
            for (int i = 0; i < NEQ; ++i) out_s.push_back(yS_arr[i]);
          }}
        }}
        ++out_idx;
      }}
      t_reached = (double)tret;
      if (flag == CV_TSTOP_RETURN) break;
    }}

    N_VDestroy(_ele_buf);
    N_VDestroy(_ewt_buf);
    N_VDestroy(_y_interp);
    if (_eleS_buf)  N_VDestroyVectorArray(_eleS_buf,  Ns_active);
    if (_ewtS_buf)  N_VDestroyVectorArray(_ewtS_buf,  Ns_active);
    if (_yS_interp) N_VDestroyVectorArray(_yS_interp, Ns_active);
  }}
#else
  for (int k = 1; k < n_times; ++k) {{
    sunrealtype tret;
{cv_fwd_call}
    if (flag < 0) {{
      return_code = flag;
      char buf[160];
      std::snprintf(buf, sizeof(buf),
                    "CVode failed at t=%.6e with flag %d", (double)times[k], flag);
      solver_msg = buf;
      break;
    }}
{sens_get_block}
    t_reached = (double)tret;
    out_t.push_back(t_reached);
    {{
      const double* y_arr = N_VGetArrayPointer(y);
      for (int i = 0; i < NEQ; ++i) out_y.push_back(y_arr[i]);
    }}
{sens_store_block}
    if (flag == CV_ROOT_RETURN) {{
      solver_msg = "Terminated: rootfunc crossed zero";
      root_terminated = true;
      break;
    }}
{equilibrate_check_block}
  }}
#endif
"""

    # --- ASA: checkpoint allocation before the forward pass, sweep after it ---
    # Per seed column, lambda' = -J'lambda and q' = -(df/dp)'lambda run from T
    # to t0, and lambda jumps by the seed row W_o at each output time.
    # A run cut into stretches interpolates the forward state by Hermite: the
    # polynomial one reads a wrong state at the very start of a stretch that a
    # root opened, where the backward problem ends.
    asa_interp = "CV_HERMITE" if asa_seg else "CV_POLYNOMIAL"
    if reverse:
        asa_init_block = """
  // --- adjoint sensitivity analysis: checkpoint allocation ---
  // @INTERP@ interpolation of the forward state between checkpoints.
  if (args.seed == nullptr) {
    cleanup();
    return res.fail(cppde::RC_ILL_INPUT,
                    "a model compiled with derivMode = reverse needs a cotangent");
  }
  if (CVodeAdjInit(cvode_mem, ASA_CHECKPOINTS, @INTERP@) < 0) {
    cleanup(); return res.fail(cppde::RC_LINIT_FAIL, "CVodeAdjInit failed");
  }
""".replace("@INTERP@", asa_interp)
        asa_sweep_block = """
  // --- the backward sweep ---
  // The result is indexed like sens1ini: state rows from lambda(t0), then
  // parameter rows from the quadrature.
  if (return_code == 0) {
    const int n_seed = args.n_seed_cols;
    const int n_out_b = (int)out_t.size();
    if (args.n_seed_rows != n_out_b) {
      char m[192];
      std::snprintf(m, sizeof(m),
                    "cotangent has %d rows but the run produced %d output rows",
                    args.n_seed_rows, n_out_b);
      cleanup(); return res.fail(cppde::RC_ILL_INPUT, m);
    }
    if (args.n_seed_states != NEQ) {
      cleanup(); return res.fail(cppde::RC_ILL_INPUT,
                                 "cotangent has the wrong state count");
    }

    const int n_phi_rows = NEQ + NPAR_ADJ;
    res.n_adj_rows = n_phi_rows;
    res.n_adj_cols = n_seed;
    res.adjoint.assign((size_t)n_phi_rows * n_seed, 0.0);

@ROOT_OFF@    int indexB = -1;
    N_Vector yB = N_VNew_Serial(NEQ, ctx);
    N_Vector qB = N_VNew_Serial(NPAR_ADJ > 0 ? NPAR_ADJ : 1, ctx);
    SUNMatrix       AB  = nullptr;
    SUNLinearSolver LSB = nullptr;
    auto cleanupB = [&]() {
      if (LSB) { SUNLinSolFree(LSB); LSB = nullptr; }
      if (AB)  { SUNMatDestroy(AB);  AB = nullptr; }
      if (ud.Jscratch) { SUNMatDestroy(ud.Jscratch); ud.Jscratch = nullptr; }
      if (qB)  { N_VDestroy(qB); qB = nullptr; }
      if (yB)  { N_VDestroy(yB); yB = nullptr; }
    };
    if (!yB || !qB) {
      cleanupB(); cleanup();
      return res.fail(cppde::RC_NO_MALLOC, "N_VNew_Serial (adjoint) failed");
    }

    for (int c = 0; c < n_seed && return_code == 0; ++c) {
      // lambda at T is the seed row for the final observation; every earlier
      // row is added when the sweep reaches its time.
      {
        double* lam = N_VGetArrayPointer(yB);
        for (int i = 0; i < NEQ; ++i)
          lam[i] = args.seed[(n_out_b - 1) + (size_t)n_out_b * i +
                             (size_t)n_out_b * NEQ * c];
      }
      { double* q = N_VGetArrayPointer(qB);
        for (int k = 0; k < NPAR_ADJ; ++k) q[k] = 0.0; }

      if (indexB < 0) {
        if (CVodeCreateB(cvode_mem, CV_BDF, &indexB) < 0 ||
            CVodeInitB(cvode_mem, indexB, adj_rhs_fn, out_t.back(), yB) < 0 ||
            CVodeSStolerancesB(cvode_mem, indexB, reltol, abstol) < 0 ||
            CVodeSetUserDataB(cvode_mem, indexB, &ud) < 0) {
          return_code = cppde::RC_LINIT_FAIL;
          solver_msg = "adjoint initialisation failed";
          break;
        }
        // The backward problem gets the caller's maxsteps, not the CVODES
        // default of 500.
        CVodeSetMaxNumStepsB(cvode_mem, indexB, maxsteps);
@LSB_SETUP@        ud.Jscratch = SUNMatClone(A);
        if (!AB || !LSB || !ud.Jscratch ||
            CVodeSetLinearSolverB(cvode_mem, indexB, LSB, AB) < 0 ||
            CVodeSetJacFnB(cvode_mem, indexB, jacB_fn) < 0) {
          return_code = cppde::RC_LINIT_FAIL;
          solver_msg = "adjoint linear solver failed";
          break;
        }
        if (NPAR_ADJ > 0) {
          const double qatol = args.gradtol > 0.0 ? args.gradtol : abstol;
          if (CVodeQuadInitB(cvode_mem, indexB, adj_quad_fn, qB) < 0 ||
              CVodeQuadSStolerancesB(cvode_mem, indexB, reltol, qatol) < 0) {
            return_code = cppde::RC_LINIT_FAIL;
            solver_msg = "adjoint quadrature failed";
            break;
          }
          // The quadrature, the gradient, enters the error test only under a
          // gradient tolerance: entries of very different scale would stall
          // the step against reltol alone.
          CVodeSetQuadErrConB(cvode_mem, indexB,
                              args.gradtol > 0.0 ? SUNTRUE : SUNFALSE);
        }
      } else {
        // Later columns re-initialise the one backward problem: CVodeB advances
        // every backward problem that exists.
        if (CVodeReInitB(cvode_mem, indexB, out_t.back(), yB) < 0) {
          return_code = cppde::RC_UNRECOGNIZED_ERR;
          solver_msg = "CVodeReInitB failed between seed columns"; break;
        }
        if (NPAR_ADJ > 0 && CVodeQuadReInitB(cvode_mem, indexB, qB) < 0) {
          return_code = cppde::RC_UNRECOGNIZED_ERR;
          solver_msg = "CVodeQuadReInitB failed between seed columns"; break;
        }
      }

      // Backwards over the output grid, adding each seed row on arrival. The
      // last row is already in lambda; the first is added after the sweep, at
      // t0, where there is nothing left to integrate.
      for (int k = n_out_b - 2; k >= 0 && return_code == 0; --k) {
        int fb = CVodeB(cvode_mem, out_t[k], CV_NORMAL);
        if (fb < 0) {
          return_code = fb; solver_msg = "CVodeB failed"; break;
        }
        sunrealtype tB;
        if (CVodeGetB(cvode_mem, indexB, &tB, yB) < 0) {
          return_code = cppde::RC_UNRECOGNIZED_ERR;
          solver_msg = "CVodeGetB failed"; break;
        }
        double* lam = N_VGetArrayPointer(yB);
        bool jumped = false;
        for (int i = 0; i < NEQ; ++i) {
          const double w = args.seed[k + (size_t)n_out_b * i +
                                     (size_t)n_out_b * NEQ * c];
          if (w != 0.0) { lam[i] += w; jumped = true; }
        }
        // A seeded jump is a new initial condition, so CVODES is re-initialised
        // there, and only there: a re-init restarts the method at order one.
        if (jumped && k > 0) {
          // Reading the quadrature out and re-initialising it with that value
          // restarts its history and keeps what it has integrated.
          if (NPAR_ADJ > 0) {
            sunrealtype tq_;
            if (CVodeGetQuadB(cvode_mem, indexB, &tq_, qB) < 0) {
              return_code = cppde::RC_UNRECOGNIZED_ERR;
              solver_msg = "CVodeGetQuadB failed at a seeded jump"; break;
            }
          }
          if (CVodeReInitB(cvode_mem, indexB, out_t[k], yB) < 0) {
            return_code = cppde::RC_UNRECOGNIZED_ERR;
            solver_msg = "CVodeReInitB failed"; break;
          }
          if (NPAR_ADJ > 0 && CVodeQuadReInitB(cvode_mem, indexB, qB) < 0) {
            return_code = cppde::RC_UNRECOGNIZED_ERR;
            solver_msg = "CVodeQuadReInitB failed"; break;
          }
        }
      }
      if (return_code != 0) break;

      if (NPAR_ADJ > 0) {
        sunrealtype tq;
        if (CVodeGetQuadB(cvode_mem, indexB, &tq, qB) < 0) {
          return_code = cppde::RC_UNRECOGNIZED_ERR;
          solver_msg = "CVodeGetQuadB failed"; break;
        }
      }

      {
        const double* lam = N_VGetArrayPointer(yB);
        const double* q   = N_VGetArrayPointer(qB);
        for (int i = 0; i < NEQ; ++i)
          res.adjoint[i + (size_t)n_phi_rows * c] = lam[i];
        // CVODES integrates the quadrature from T down to t0 with xi(T) = 0, so
        // xi(t0) = -int_{t0}^{T} fQB dt; fQB includes the adjoint's minus sign,
        // and the two cancel to int lambda' (df/dp) dt.
        for (int k = 0; k < NPAR_ADJ; ++k)
          res.adjoint[(NEQ + k) + (size_t)n_phi_rows * c] = q[k];
      }

    }
    cleanupB();
  }
"""
    else:
        asa_init_block = ""
        asa_sweep_block = ""

    # --- ASA over events ---
    # The forward run is cut at every jump into stretches, each on its own
    # CVODES memory and checkpoints. The backward sweep runs them last to first
    # and passes lambda and the quadrature through each jump's adjoint.
    asa_structs = asa_decl = asa_cleanup = ""
    if asa_seg:
        asa_structs = """
// One stretch of the forward run between jumps: its CVODES memory (null when
// nothing was integrated on it), its backward problem, and what it spans.
struct AsaSeg {
  void* mem = nullptr;
  SUNMatrix A = nullptr;
  SUNLinearSolver LS = nullptr;
  SUNMatrix AB = nullptr;
  SUNLinearSolver LSB = nullptr;
  int indexB = -1;
  bool used = false;
  double t0 = 0.0, t1 = 0.0;
  int row0 = 0, row1 = 0;     // output rows [row0, row1)
  int jump0 = 0, jump1 = 0;   // jumps applied at its start [jump0, jump1)
  std::vector<signed char> modes;   // of the state switches, none without
};

// A jump of the forward run: one time event, or the root events of one crossing.
struct AsaJump {
  int kind = 0;               // 0 time event, 1 root events
  int ev = -1;                // index into the sorted time events
  std::vector<int> trig;      // triggered root events
  double t = 0.0;
  std::vector<double> x_old, f_old, f_new;
};
"""
        asa_decl = """  std::vector<AsaSeg> asa_segs(1);
  asa_segs[0].t0 = t0;
  std::vector<AsaJump> asa_jumps;
  int asa_jump_mark = 0;
  std::string asa_msg;
"""
        asa_cleanup = """    for (auto& s_ : asa_segs) {
      if (s_.mem) CVodeFree(&s_.mem);
      if (s_.LSB) SUNLinSolFree(s_.LSB);
      if (s_.AB)  SUNMatDestroy(s_.AB);
      if (s_.LS)  SUNLinSolFree(s_.LS);
      if (s_.A)   SUNMatDestroy(s_.A);
      s_ = AsaSeg();
    }
"""

        jump_parts = []
        if has_time_events:
            jump_parts.append("""      if (J.kind == 0) {
        const TimeEvent& ev = time_events[J.ev];
        const int v = ev.var_idx;
        const double wv = w[v];
        // W is the cotangent of the event time.
        double W = wv * (ev.dg_dt_fn(xo, J.t) - fn[v]);
        for (int i = 0; i < NEQ; ++i) {
          const double gx = ev.dg_dx_fn(xo, J.t, i);
          W += wv * gx * fo[i];
          if (i != v) W += w[i] * (fo[i] - fn[i]);
          lam[i] = (i == v ? 0.0 : w[i]) + wv * gx;
        }
        for (int k = 0; k < NPAR_ADJ; ++k)
          q[k] += wv * ev.dg_dp_fn(xo, J.t, k) + W * ev.dt_dp_fn(xo, J.t, k);
        return;
      }
""")
        if has_root_events:
            jump_parts.append("""      // The first non-terminal event defines the crossing time, as forward.
      std::vector<int> writer(NEQ, -1);
      int ref = -1;
      for (int j : J.trig) {
        const auto& e = event_roots[j];
        if (e.terminal) continue;
        if (e.var_idx >= 0) writer[e.var_idx] = j;
        if (ref < 0) ref = j;
      }
      if (ref < 0) return;
      const auto& R = event_roots[ref];
      std::vector<double> rx(NEQ);
      double g_dot = R.dr_dt_fn(xo, J.t);
      for (int i = 0; i < NEQ; ++i) {
        rx[i] = R.dr_dx_fn(xo, J.t, i);
        g_dot += rx[i] * fo[i];
      }
      double W = 0.0;
      for (int i = 0; i < NEQ; ++i) lam[i] = (writer[i] < 0) ? w[i] : 0.0;
      for (int i = 0; i < NEQ; ++i) {
        if (writer[i] < 0) { W += w[i] * (fo[i] - fn[i]); continue; }
        const auto& e = event_roots[writer[i]];
        const double wi = w[i];
        if (wi == 0.0) continue;
        W += wi * (e.dg_dt_fn(xo, J.t) - fn[i]);
        for (int m = 0; m < NEQ; ++m) {
          const double gx = e.dg_dx_fn(xo, J.t, m);
          W += wi * gx * fo[m];
          lam[m] += wi * gx;
        }
        for (int k = 0; k < NPAR_ADJ; ++k) q[k] += wi * e.dg_dp_fn(xo, J.t, k);
      }
      // dt_e = -(r_x dx + r_p dp) / g_dot, zero on a tangential crossing.
      if (g_dot != 0.0) {
        const double c = -W / g_dot;
        for (int m = 0; m < NEQ; ++m) lam[m] += c * rx[m];
        for (int k = 0; k < NPAR_ADJ; ++k) q[k] += c * R.dr_dp_fn(xo, J.t, k);
      }
""")
        jump_adjoint = ("""    // The adjoint of one jump x+ = h(x-, p, t_e): lam- = h_x' lam+ and
    // q += h_p' lam+, plus the event time's cotangent W through dt_e/dp or,
    // for a root, through dt_e/dx- and dt_e/dp.
    auto asa_jump_adjoint = [&](const AsaJump& J, double* lam, double* q) {
      const double* xo = J.x_old.data();
      const double* fo = J.f_old.data();
      const double* fn = J.f_new.data();
      const std::vector<double> w(lam, lam + NEQ);
""" + "".join(jump_parts) + "    };\n")

        asa_sweep_block = """
  // --- the backward sweep, stretch by stretch from the last ---
  // The result is indexed like sens1ini: state rows from lambda(t0), then
  // parameter rows from the quadrature.
  if (return_code == 0) {
    const int n_seed = args.n_seed_cols;
    const int n_out_b = (int)out_t.size();
    if (args.n_seed_rows != n_out_b) {
      char m[192];
      std::snprintf(m, sizeof(m),
                    "cotangent has %d rows but the run produced %d output rows",
                    args.n_seed_rows, n_out_b);
      cleanup(); return res.fail(cppde::RC_ILL_INPUT, m);
    }
    if (args.n_seed_states != NEQ) {
      cleanup(); return res.fail(cppde::RC_ILL_INPUT,
                                 "cotangent has the wrong state count");
    }

    const int n_phi_rows = NEQ + NPAR_ADJ;
    res.n_adj_rows = n_phi_rows;
    res.n_adj_cols = n_seed;
    res.adjoint.assign((size_t)n_phi_rows * n_seed, 0.0);

    {
      AsaSeg& L = asa_segs.back();
      L.t1 = std::max(L.t0, out_t.back());
      L.row1 = n_out_b;
    }
    // The replay from a checkpoint must not stop at the roots again.
    if (cvode_mem) CVodeRootInit(cvode_mem, 0, nullptr);

    const int n_q = NPAR_ADJ > 0 ? NPAR_ADJ : 1;
    std::vector<double> lam_all((size_t)NEQ * n_seed, 0.0);
    std::vector<double> q_all((size_t)n_q * n_seed, 0.0);
    N_Vector yB = N_VNew_Serial(NEQ, ctx);
    N_Vector qB = N_VNew_Serial(n_q, ctx);
    ud.Jscratch = SUNMatClone(A);
    auto cleanupB = [&]() {
      if (ud.Jscratch) { SUNMatDestroy(ud.Jscratch); ud.Jscratch = nullptr; }
      if (qB)  { N_VDestroy(qB); qB = nullptr; }
      if (yB)  { N_VDestroy(yB); yB = nullptr; }
    };
    if (!yB || !qB || !ud.Jscratch) {
      cleanupB(); cleanup();
      return res.fail(cppde::RC_NO_MALLOC, "N_VNew_Serial (adjoint) failed");
    }

    // Rows at one time within rounding, such as the row before a root and the
    // root time itself, are added together.
    auto near = [](double a, double b) {
      return std::abs(a - b) <= 1e-13 * std::max(1.0, std::abs(b));
    };
    auto add_seed = [&](int k, int c, double* lam) -> bool {
      bool any = false;
      for (int i = 0; i < NEQ; ++i) {
        const double w = args.seed[k + (size_t)n_out_b * i +
                                   (size_t)n_out_b * NEQ * c];
        if (w != 0.0) { lam[i] += w; any = true; }
      }
      return any;
    };
@JUMP_ADJOINT@
    const int n_seg = (int)asa_segs.size();
    for (int s = n_seg - 1; s >= 0 && return_code == 0; --s) {
      AsaSeg& S = asa_segs[s];
@USE_MODES@
      void* mem = (s == n_seg - 1) ? cvode_mem : S.mem;
      const bool flows = mem != nullptr && S.used && S.t1 > S.t0 && !near(S.t0, S.t1);
      SUNMatrix&       AB  = S.AB;
      SUNLinearSolver& LSB = S.LSB;
      for (int c = 0; c < n_seed && return_code == 0; ++c) {
        double* lam_c = lam_all.data() + (size_t)NEQ * c;
        double* q_c   = q_all.data() + (size_t)n_q * c;
        int k = S.row1 - 1;
        double tB = S.t1;
        while (k >= S.row0 && near(out_t[k], tB)) { add_seed(k, c, lam_c); --k; }

        if (!flows) {
          // Nothing to integrate: every row here is at the stretch's one time.
          for (; k >= S.row0; --k) add_seed(k, c, lam_c);
        } else {
          double* lam = N_VGetArrayPointer(yB);
          double* q   = N_VGetArrayPointer(qB);
          for (int i = 0; i < NEQ; ++i) lam[i] = lam_c[i];
          for (int i = 0; i < n_q; ++i) q[i] = q_c[i];
          if (S.indexB < 0) {
            if (CVodeCreateB(mem, CV_BDF, &S.indexB) < 0 ||
                CVodeInitB(mem, S.indexB, adj_rhs_fn, tB, yB) < 0 ||
                CVodeSStolerancesB(mem, S.indexB, reltol, abstol) < 0 ||
                CVodeSetUserDataB(mem, S.indexB, &ud) < 0) {
              return_code = cppde::RC_LINIT_FAIL;
              solver_msg = "adjoint initialisation failed";
              break;
            }
            CVodeSetMaxNumStepsB(mem, S.indexB, maxsteps);
@LSB_SETUP@            if (!AB || !LSB ||
                CVodeSetLinearSolverB(mem, S.indexB, LSB, AB) < 0 ||
                CVodeSetJacFnB(mem, S.indexB, jacB_fn) < 0) {
              return_code = cppde::RC_LINIT_FAIL;
              solver_msg = "adjoint linear solver failed";
              break;
            }
            if (NPAR_ADJ > 0) {
              const double qatol = args.gradtol > 0.0 ? args.gradtol : abstol;
              if (CVodeQuadInitB(mem, S.indexB, adj_quad_fn, qB) < 0 ||
                  CVodeQuadSStolerancesB(mem, S.indexB, reltol, qatol) < 0) {
                return_code = cppde::RC_LINIT_FAIL;
                solver_msg = "adjoint quadrature failed";
                break;
              }
              CVodeSetQuadErrConB(mem, S.indexB,
                                  args.gradtol > 0.0 ? SUNTRUE : SUNFALSE);
            }
          } else {
            if (CVodeReInitB(mem, S.indexB, tB, yB) < 0 ||
                (NPAR_ADJ > 0 && CVodeQuadReInitB(mem, S.indexB, qB) < 0)) {
              return_code = cppde::RC_UNRECOGNIZED_ERR;
              solver_msg = "CVodeReInitB failed between seed columns"; break;
            }
          }

          // Backwards to each earlier row and on to the stretch's start; a
          // seeded row re-initialises the backward problem, as without events.
          while (return_code == 0) {
            const double tk = (k >= S.row0) ? std::max(out_t[k], S.t0) : S.t0;
            if (tk < tB && !near(tk, tB)) {
              int fb = CVodeB(mem, tk, CV_NORMAL);
              if (fb < 0) { return_code = fb; solver_msg = "CVodeB failed"; break; }
              sunrealtype tret_b;
              if (CVodeGetB(mem, S.indexB, &tret_b, yB) < 0) {
                return_code = cppde::RC_UNRECOGNIZED_ERR;
                solver_msg = "CVodeGetB failed"; break;
              }
              tB = tk;
            }
            if (k < S.row0) break;
            bool jumped = false;
            while (k >= S.row0 && near(std::max(out_t[k], S.t0), tB)) {
              jumped = add_seed(k, c, lam) || jumped;
              --k;
            }
            if (jumped && !near(tB, S.t0)) {
              sunrealtype tq_;
              if ((NPAR_ADJ > 0 && CVodeGetQuadB(mem, S.indexB, &tq_, qB) < 0) ||
                  CVodeReInitB(mem, S.indexB, tB, yB) < 0 ||
                  (NPAR_ADJ > 0 && CVodeQuadReInitB(mem, S.indexB, qB) < 0)) {
                return_code = cppde::RC_UNRECOGNIZED_ERR;
                solver_msg = "CVodeReInitB failed"; break;
              }
            }
          }
          if (return_code != 0) break;
          if (NPAR_ADJ > 0) {
            sunrealtype tq;
            if (CVodeGetQuadB(mem, S.indexB, &tq, qB) < 0) {
              return_code = cppde::RC_UNRECOGNIZED_ERR;
              solver_msg = "CVodeGetQuadB failed"; break;
            }
          }
          for (int i = 0; i < NEQ; ++i) lam_c[i] = lam[i];
          for (int i = 0; i < n_q; ++i) q_c[i] = q[i];
        }

        for (int j = S.jump1 - 1; j >= S.jump0; --j)
          asa_jump_adjoint(asa_jumps[j], lam_c, q_c);
      }
    }

    if (return_code == 0) {
      for (int c = 0; c < n_seed; ++c) {
        for (int i = 0; i < NEQ; ++i)
          res.adjoint[i + (size_t)n_phi_rows * c] = lam_all[i + (size_t)NEQ * c];
        for (int k = 0; k < NPAR_ADJ; ++k)
          res.adjoint[(NEQ + k) + (size_t)n_phi_rows * c] = q_all[k + (size_t)n_q * c];
      }
    }
    cleanupB();
  }
""".replace("@JUMP_ADJOINT@", jump_adjoint).replace(
            "@USE_MODES@\n",
            "      // The backward problem reads f as the stretch ran it.\n"
            "      cppde::detail::use_switch_modes(S.modes);\n" if n_modes else "")
    else:
        # A user rootfunc stops the run; the replay from a checkpoint must not
        # stop at it again.
        root_off = ("    CVodeRootInit(cvode_mem, 0, nullptr);\n"
                    if reverse and need_cvode_root else "")
        asa_sweep_block = asa_sweep_block.replace("@ROOT_OFF@", root_off)

    # A seed reaches the adjoint under reverse and means nothing otherwise, so
    # each direction refuses the other rather than ignoring the argument.
    if reverse:
        seed_guard = (
            '  if (a.seed == nullptr)\n'
            '    Rf_error("a model compiled with derivMode = reverse needs a cotangent");')
    else:
        seed_guard = (
            '  if (a.seed != nullptr)\n'
            '    Rf_error("this cvode model was compiled with derivMode = forward; '
            'recompile with derivMode = reverse");')

    # Linear-solver setup
    if use_sparse:
        ls_includes = (
            "#include <sunmatrix/sunmatrix_sparse.h>\n"
            "#include <sunlinsol/sunlinsol_klu.h>\n"
        )
        jac_decl = ("static int jac_fn(sunrealtype t, N_Vector y, N_Vector fy, "
                    "SUNMatrix J, void* ud_vp, N_Vector, N_Vector, N_Vector);")
        jac_impl = f"""
static int jac_fn(sunrealtype t, N_Vector y, N_Vector fy,
                  SUNMatrix J, void* ud_vp,
                  N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {{
  (void)t; (void)fy; (void)tmp1; (void)tmp2; (void)tmp3;
  // An exception must not unwind through SUNDIALS' C frames (UB); these
  // bodies allocate, so bad_alloc is reachable.
  try {{
  UserData* ud = static_cast<UserData*>(ud_vp);
  const double* params = ud->params.data();
  const double* x = N_VGetArrayPointer(y);
  (void)x; (void)params;
{forcing_local}
{jac_body}
  return 0;
  }} catch (...) {{
    return -1;   // unrecoverable
  }}
}}
"""
        # CVODES assumes the replay from a checkpoint takes the steps the run
        # took. KLU's refactorisation depends on earlier pivots, a full
        # factorisation at every setup on the matrix alone.
        if reverse:
            jac_impl += """
static int klu_setup_fresh(SUNLinearSolver S, SUNMatrix A) {
  int flag = SUNLinSol_KLUReInit(S, A, SUNSparseMatrix_NNZ(A), SUNKLU_REINIT_PARTIAL);
  if (flag != 0) return flag;
  return SUNLinSolSetup_KLU(S, A);
}
"""
        klu_fresh = "  LS->ops->setup = klu_setup_fresh;\n" if reverse else ""
        klu_btf = 1 if klu_settings and klu_settings["use_btf"] else 0
        klu_ord = int(klu_settings["ordering"]) if klu_settings else 0
        def ls_setup_for(fail):
            return f"""  A = SUNSparseMatrix(NEQ, NEQ, {nnz_total}, CSC_MAT, ctx);
  if (!A) {fail("RC_NO_MALLOC", "SUNSparseMatrix failed")}
  LS = SUNLinSol_KLU(y, A, ctx);
  if (!LS) {fail("RC_LINIT_FAIL", "SUNLinSol_KLU failed")}
  {{
    sun_klu_common* klu_c = SUNLinSol_KLUGetCommon(LS);
    if (klu_c) {{ klu_c->btf = {klu_btf}; klu_c->ordering = {klu_ord}; }}
  }}
{klu_fresh}  if (CVodeSetLinearSolver(cvode_mem, LS, A) < 0) {fail("RC_LINIT_FAIL", "CVodeSetLinearSolver failed")}
  if (CVodeSetJacFn(cvode_mem, jac_fn) < 0) {fail("RC_LINIT_FAIL", "CVodeSetJacFn failed")}
"""
    else:
        ls_includes = ("#include <sunlinsol/sunlinsol_lapackdense.h>\n"
                       if use_lapack else "")
        jac_decl = ("static int jac_fn(sunrealtype t, N_Vector y, N_Vector fy, "
                    "SUNMatrix J, void* ud_vp, N_Vector, N_Vector, N_Vector);")
        jac_impl = f"""
static int jac_fn(sunrealtype t, N_Vector y, N_Vector fy,
                  SUNMatrix J, void* ud_vp,
                  N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {{
  (void)t; (void)fy; (void)tmp1; (void)tmp2; (void)tmp3;
  // An exception must not unwind through SUNDIALS' C frames (UB); these
  // bodies allocate, so bad_alloc is reachable.
  try {{
  UserData* ud = static_cast<UserData*>(ud_vp);
  const double* params = ud->params.data();
  const double* x = N_VGetArrayPointer(y);
  (void)x; (void)params;
{forcing_local}
  SUNMatZero(J);
{jac_body}
  return 0;
  }} catch (...) {{
    return -1;   // unrecoverable
  }}
}}
"""
        ls_ctor = "SUNLinSol_LapackDense" if use_lapack else "SUNLinSol_Dense"
        def ls_setup_for(fail):
            return f"""  A = SUNDenseMatrix(NEQ, NEQ, ctx);
  if (!A) {fail("RC_NO_MALLOC", "SUNDenseMatrix failed")}
  LS = {ls_ctor}(y, A, ctx);
  if (!LS) {fail("RC_LINIT_FAIL", ls_ctor + " failed")}
  if (CVodeSetLinearSolver(cvode_mem, LS, A) < 0) {fail("RC_LINIT_FAIL", "CVodeSetLinearSolver failed")}
  if (CVodeSetJacFn(cvode_mem, jac_fn) < 0) {fail("RC_LINIT_FAIL", "CVodeSetJacFn failed")}
"""

    ls_setup = ls_setup_for(_fail_main)

    # A jump closes the stretch on the current memory; the next one starts on a
    # fresh memory, or on the same one when nothing was integrated on it yet.
    asa_restart_lambda = ""
    if asa_seg:
        asa_restart_lambda = f"""  auto asa_restart = [&](double t_r) -> int {{
    const size_t cur = asa_segs.size() - 1;
    asa_segs[cur].t1 = t_r;
    asa_segs[cur].row1 = (int)out_t.size();
    AsaSeg nxt;
    nxt.t0 = t_r;
    nxt.row0 = (int)out_t.size();
    nxt.jump0 = asa_jump_mark;
    nxt.jump1 = (int)asa_jumps.size();
    asa_jump_mark = nxt.jump1;
@SNAP_MODES@    if (!asa_segs[cur].used) {{
      asa_segs.push_back(nxt);
      return CVodeReInit(cvode_mem, t_r, y);
    }}
    // The replay from a checkpoint must not stop at the roots again.
    CVodeRootInit(cvode_mem, 0, nullptr);
    asa_segs[cur].mem = cvode_mem;
    asa_segs[cur].A = A;
    asa_segs[cur].LS = LS;
    cvode_mem = nullptr; A = nullptr; LS = nullptr;
    asa_segs.push_back(nxt);
    cvode_mem = CVodeCreate({cv_method}, ctx);
    if (!cvode_mem) {_fail_seg("", "CVodeCreate failed")}
    if (CVodeInit(cvode_mem, rhs_fn, t_r, y) < 0) {_fail_seg("", "CVodeInit failed")}
    if (CVodeSStolerances(cvode_mem, reltol, abstol) < 0) {_fail_seg("", "CVodeSStolerances failed")}
    if (CVodeSetUserData(cvode_mem, &ud) < 0) {_fail_seg("", "CVodeSetUserData failed")}
    CVodeSetMaxNumSteps(cvode_mem, maxsteps);
    if (hini > 0.0) CVodeSetInitStep(cvode_mem, hini);
{ls_setup_for(_fail_seg)}{rootfunc_init(_fail_seg)}    if (CVodeAdjInit(cvode_mem, ASA_CHECKPOINTS, {asa_interp}) < 0) {_fail_seg("", "CVodeAdjInit failed")}
    return 0;
  }};
""".replace("@SNAP_MODES@", "    nxt.modes = cppde::detail::switch_mode_snapshot();\n"
                    if n_modes else "")

    # The backward problem lambda' = -J' lambda gets the forward solver kind and
    # its Jacobian -J', formed from jac_fn into a scratch matrix of the forward
    # pattern.
    asa_ud_members = ""
    jacB_decl = jacB_impl = lsB_setup = ""
    if reverse:
        asa_ud_members = ("  SUNMatrix Jscratch = nullptr;              // J for jacB_fn\n"
                          "  std::vector<sunindextype> jacB_next;       // CSC transpose cursor\n")
        jacB_decl = ("static int jacB_fn(sunrealtype t, N_Vector y, N_Vector yB, "
                     "N_Vector fyB, SUNMatrix JB, void* ud_vp, N_Vector, N_Vector, "
                     "N_Vector);")
        if use_sparse:
            jacB_fill = """  // -J' in CSC is J in CSR: count per row, then scatter.
  const sunindextype* Jp = SUNSparseMatrix_IndexPointers(J);
  const sunindextype* Ji = SUNSparseMatrix_IndexValues(J);
  const sunrealtype*  Jx = SUNSparseMatrix_Data(J);
  sunindextype* Bp = SUNSparseMatrix_IndexPointers(JB);
  sunindextype* Bi = SUNSparseMatrix_IndexValues(JB);
  sunrealtype*  Bx = SUNSparseMatrix_Data(JB);
  for (sunindextype c = 0; c <= NEQ; ++c) Bp[c] = 0;
  for (sunindextype k = 0; k < Jp[NEQ]; ++k) ++Bp[Ji[k] + 1];
  for (sunindextype c = 0; c < NEQ; ++c) Bp[c + 1] += Bp[c];
  ud->jacB_next.assign(Bp, Bp + NEQ);
  for (sunindextype j = 0; j < NEQ; ++j)
    for (sunindextype k = Jp[j]; k < Jp[j + 1]; ++k) {
      const sunindextype pos = ud->jacB_next[Ji[k]]++;
      Bi[pos] = j;
      Bx[pos] = -Jx[k];
    }"""
            lsB_setup = f"""        AB  = SUNSparseMatrix(NEQ, NEQ, {nnz_total}, CSC_MAT, ctx);
        LSB = AB ? SUNLinSol_KLU(yB, AB, ctx) : nullptr;
        if (LSB) {{
          sun_klu_common* klu_cB = SUNLinSol_KLUGetCommon(LSB);
          if (klu_cB) {{ klu_cB->btf = {klu_btf}; klu_cB->ordering = {klu_ord}; }}
        }}
"""
        else:
            jacB_fill = """  for (sunindextype j = 0; j < NEQ; ++j)
    for (sunindextype i = 0; i < NEQ; ++i)
      SM_ELEMENT_D(JB, i, j) = -SM_ELEMENT_D(J, j, i);"""
            lsB_setup = f"""        AB  = SUNDenseMatrix(NEQ, NEQ, ctx);
        LSB = AB ? {ls_ctor}(yB, AB, ctx) : nullptr;
"""
        jacB_impl = f"""
// ---- Backward Jacobian: d(lambda')/d(lambda) = -J(x,p)' ----
static int jacB_fn(sunrealtype t, N_Vector y, N_Vector yB, N_Vector fyB,
                   SUNMatrix JB, void* ud_vp,
                   N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B) {{
  (void)yB; (void)fyB; (void)tmp1B; (void)tmp2B; (void)tmp3B;
  try {{
  UserData* ud = static_cast<UserData*>(ud_vp);
  SUNMatrix J = ud->Jscratch;
  if (jac_fn(t, y, nullptr, J, ud_vp, nullptr, nullptr, nullptr) != 0) return -1;
{jacB_fill}
  return 0;
  }} catch (...) {{
    return -1;
  }}
}}
"""
    asa_sweep_block = asa_sweep_block.replace("@LSB_SETUP@", lsB_setup)

    sens_decl = ""
    sens_impl = ""
    if deriv:
        sens_decl = ("static int sens_rhs1_fn(int Ns, sunrealtype t, N_Vector y, "
                     "N_Vector ydot, int iS, N_Vector yS, N_Vector ySdot, "
                     "void* ud_vp, N_Vector tmp1, N_Vector tmp2);")
        phi_rows = n_states + n_params
        sens_impl = f"""
static int sens_rhs1_fn(int Ns, sunrealtype t,
                        N_Vector y, N_Vector ydot,
                        int iS, N_Vector yS, N_Vector ySdot,
                        void* ud_vp,
                        N_Vector tmp1, N_Vector tmp2) {{
  (void)Ns; (void)t; (void)ydot; (void)tmp1; (void)tmp2;
  // An exception must not unwind through SUNDIALS' C frames (UB); these
  // bodies allocate, so bad_alloc is reachable.
  try {{
  UserData* ud = static_cast<UserData*>(ud_vp);
  const double* params = ud->params.data();
  const double* x = N_VGetArrayPointer(y);
  const double* yS_arr  = N_VGetArrayPointer(yS);
  double*       ySdot_arr = N_VGetArrayPointer(ySdot);
  // parameter rows of column iS of Phi'(theta)
  const double* _Mp = ud->Phi_prime.data() + NEQ + {phi_rows} * iS;
  (void)x; (void)params; (void)_Mp;
{forcing_local}
  // ySdot = J(x,p) yS + (df/dp) _Mp
{sens_body}
  return 0;
  }} catch (...) {{
    return -1;   // unrecoverable
  }}
}}
"""

    # --- Adjoint entry points, emitted only under reverse ---
    if reverse:
        adj_decl = (
            "static int adj_rhs_fn(sunrealtype t, N_Vector y, N_Vector yB,\n"
            "                      N_Vector yBdot, void* ud_vp);\n"
            "static int adj_quad_fn(sunrealtype t, N_Vector y, N_Vector yB,\n"
            "                       N_Vector qBdot, void* ud_vp);")
        adj_impl = f"""
// ---- Adjoint right-hand side: lambda' = -J(x,p)' lambda ----
// y is the forward state as CVODES interpolates it from the checkpoints.
static int adj_rhs_fn(sunrealtype t, N_Vector y, N_Vector yB,
                      N_Vector yBdot, void* ud_vp) {{
  (void)t;
  try {{
  UserData* ud = static_cast<UserData*>(ud_vp);
  const double* params = ud->params.data();
  const double* x      = N_VGetArrayPointer(y);
  const double* lam    = N_VGetArrayPointer(yB);
  double*       lamdot = N_VGetArrayPointer(yBdot);
  (void)x; (void)params;
{forcing_local}
{adj_rhs_body}
  return 0;
  }} catch (...) {{
    return -1;
  }}
}}

// ---- Adjoint quadrature: q' = -(df/dp)' lambda ----
// Integrated backwards from T to t0, so what it accumulates is the parameter
// half of the gradient. The state half needs no quadrature: it is lambda(t0).
static int adj_quad_fn(sunrealtype t, N_Vector y, N_Vector yB,
                       N_Vector qBdot, void* ud_vp) {{
  (void)t;
  try {{
  UserData* ud = static_cast<UserData*>(ud_vp);
  const double* params = ud->params.data();
  const double* x    = N_VGetArrayPointer(y);
  const double* lam  = N_VGetArrayPointer(yB);
  double*       qdot = N_VGetArrayPointer(qBdot);
  (void)x; (void)params;
{forcing_local}
{adj_quad_body}
  return 0;
  }} catch (...) {{
    return -1;
  }}
}}
"""
    else:
        adj_decl = ""
        adj_impl = ""

    return f"""/** Code auto-generated by cppDE {version} (CVODE backend) **/

#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include <cvodes/cvodes.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
{ls_includes}#include <sundials/sundials_types.h>
#include <sundials/sundials_context.h>
#include <cppde/cppde_scalar_ops.hpp>
{linmap_include}#include <cppde/cppde_step_trace.hpp>
#include <cppde/cppde_r_batch.hpp>
{forcing_include}{event_block_includes}

namespace {{
{data_block}
constexpr int NEQ    = {n_states};
constexpr int NPARMS = {n_global};           // n_states + n_params (flat layout)
// Rows the adjoint quadrature covers: the dynamic parameters. The state half
// of the answer is lambda(t0) and needs no quadrature.
constexpr int NPAR_ADJ = NPARMS - NEQ;
// Accepted forward steps between checkpoints (Nd of CVodeAdjInit).
constexpr int ASA_CHECKPOINTS = {asa_checkpoints};

struct UserData {{
  std::vector<double> params;               // length NPARMS
{reparam_ud_members}{forcing_ud_members}{events_ud_members}{asa_ud_members}}};

static int rhs_fn(sunrealtype t, N_Vector y, N_Vector ydot, void* ud_vp);
{jac_decl}
{sens_decl}
{adj_decl}
{jacB_decl}
{rootfunc_decl}

// ---- RHS ----
static int rhs_fn(sunrealtype t, N_Vector y, N_Vector ydot, void* ud_vp) {{
  (void)t;
  // An exception must not unwind through SUNDIALS' C frames (UB); these
  // bodies allocate, so bad_alloc is reachable.
  try {{
  UserData* ud = static_cast<UserData*>(ud_vp);
  const double* params = ud->params.data();
  const double* x = N_VGetArrayPointer(y);
  double* ydot_arr = N_VGetArrayPointer(ydot);
  (void)x; (void)params;
{forcing_local}
{ode_body}
  return 0;
  }} catch (...) {{
    return -1;   // unrecoverable
  }}
}}
{jac_impl}
{sens_impl}
{adj_impl}
{jacB_impl}
{rootfunc_impl}
{event_struct}{asa_structs}

// ---- Step trace (CVODE_STEP_TRACE) ----
// Appends trace rows to the buffer of cppde_step_trace.hpp, which reaches R as
// the `trace` element; fields the CVODES API does not expose are NaN.
#ifdef CVODE_STEP_TRACE
static void cvode_emit_trace_row(void* cvode_mem,
                                 N_Vector ele_buf, N_Vector ewt_buf,
                                 N_Vector* eleS_buf, N_Vector* ewtS_buf,
                                 int ns_active, double t_reached,
                                 const char* row_mode = "CVODE",
                                 long nst_carry = 0, long nfe_carry = 0,
                                 long nje_carry = 0, long nsetups_carry = 0) {{
  long n_steps = 0, n_fe = 0, n_je = 0, n_setups = 0;
  int  last_order = 0;
  double last_h = 0.0;
  CVodeGetNumSteps(cvode_mem, &n_steps);
  CVodeGetNumRhsEvals(cvode_mem, &n_fe);
  CVodeGetNumJacEvals(cvode_mem, &n_je);
  CVodeGetNumLinSolvSetups(cvode_mem, &n_setups);
  CVodeGetLastOrder(cvode_mem, &last_order);
  CVodeGetLastStep(cvode_mem, &last_h);

  // dsm = WRMS(ele * ewt) over the states: CVodeGetEstLocalErrors returns the
  // scaled estimate acor * tq[2]. The public API exposes neither acor, tq[2]
  // nor the sensitivity part of the error test, so those fields are NaN.
  (void)eleS_buf; (void)ewtS_buf; (void)ns_active;
  double sumsq = 0.0;
  long   N_eff = 0;
  if (CVodeGetEstLocalErrors(cvode_mem, ele_buf) == 0 &&
      CVodeGetErrWeights(cvode_mem, ewt_buf)    == 0) {{
    const double* e = N_VGetArrayPointer(ele_buf);
    const double* w = N_VGetArrayPointer(ewt_buf);
    for (int i = 0; i < NEQ; ++i) {{
      double r = e[i] * w[i];
      sumsq += r * r;
      ++N_eff;
    }}
  }}
  double dsm_reconstructed =
    (N_eff > 0) ? std::sqrt(sumsq / static_cast<double>(N_eff)) : 0.0;

  double nan_val = std::nan("");
  auto* _tbp = cppde::ndf_detail::trace_sink();
  if (!_tbp) return;
  auto& tb = *_tbp;
  n_steps += nst_carry; n_fe += nfe_carry;
  n_je += nje_carry;    n_setups += nsetups_carry;
  tb.nst.push_back(static_cast<int>(n_steps));
  tb.t.push_back(t_reached);
  tb.h.push_back(last_h);
  tb.q.push_back(last_order);
  tb.dsm.push_back(dsm_reconstructed);
  tb.acnrm.push_back(nan_val);        // raw acor not exposed by CVODES API
  tb.acnrm_state.push_back(nan_val);  // see comment above `dsm_reconstructed`
  tb.tq2.push_back(nan_val);
  tb.gamma.push_back(nan_val);
  tb.gamrat.push_back(nan_val);
  tb.newton_conv.push_back(1);
  tb.mode.emplace_back(row_mode);
  tb.nfe.push_back(static_cast<int>(n_fe));
  tb.njev.push_back(static_cast<int>(n_je));
  tb.nsetups.push_back(static_cast<int>(n_setups));
  tb.setup_reason.emplace_back("");
}}
#endif  // CVODE_STEP_TRACE

}} // anonymous namespace

// =====================================================================
// R entry points, with the signature of the native solve_<model>
// =====================================================================

// Pure C++ solve: no R API, no escaping exceptions.
static int solve_impl(const cppde::rbatch::solve_args& args,
                      cppde::rbatch::solve_result& res) noexcept
{{
try {{
  cppde::ndf_detail::trace_scope _cppde_trace_scope(res.trace);

  if (args.sens2ini != nullptr)
    return res.fail(cppde::RC_ILL_INPUT, "hessian not supported by CVODES backend");

  constexpr bool deriv = {deriv_flag};

  const int n_times_in = args.n_times;
  if (n_times_in < 1) return res.fail(cppde::RC_ILL_INPUT, "times must be non-empty");
  const double abstol = args.abstol;
  const double reltol = args.reltol;
  const double hini   = args.hini;
  const double root_tol = args.root_tol; (void)root_tol;
  const long   maxsteps = static_cast<long>(args.maxsteps);

  UserData ud;
  ud.params.assign(args.params, args.params + NPARMS);
{mode_scope}
  {{
    const int mr = args.maxroot;
    ud.maxroot = (mr > 0) ? mr : 1;
    ud.root_fired.assign({n_event_roots}, 0);
  }}

{forcing_init_block}
  // With deriv, R supplies sens1ini as [NEQ + n_params, M], zero rows for fixed
  // entries; Ns_active = M, and M = 0 skips the sensitivity integration.
  int Ns_active_tmp = 0;
  if (deriv) {{
    if (args.sens1ini == nullptr)
      return res.fail(cppde::RC_ILL_INPUT, "tangent is required when deriv = TRUE");
    Ns_active_tmp = args.n_sens1_cols;
  }}
  const int Ns_active = Ns_active_tmp;

  // --- times ---
  std::vector<double> times(args.times, args.times + n_times_in);
{cv_zero_block}  std::sort(times.begin(), times.end());
  times.erase(std::unique(times.begin(), times.end()), times.end());
  const int n_times = static_cast<int>(times.size());
  const double t0 = times[0];

  // --- output buffers ---
  std::vector<double> out_t;  out_t.reserve(n_times);
  std::vector<double> out_y;  out_y.reserve(static_cast<size_t>(n_times) * NEQ);
  std::vector<double> out_s;
  if (deriv) out_s.reserve(static_cast<size_t>(n_times) * NEQ * Ns_active);

  // --- SUNDIALS resources ---
  SUNContext ctx = nullptr;
  void* cvode_mem = nullptr;
  N_Vector y  = nullptr;
  N_Vector* yS = nullptr;
  N_Vector f_buf = nullptr;  // scratch for rhs evaluation at events
  SUNMatrix A = nullptr;
  SUNLinearSolver LS = nullptr;
  int Ns_alloc = 0;

{asa_decl}  auto cleanup = [&]() {{
{asa_cleanup}    if (yS) {{ N_VDestroyVectorArray(yS, Ns_alloc); yS = nullptr; }}
    if (f_buf) {{ N_VDestroy(f_buf); f_buf = nullptr; }}
    if (cvode_mem) {{ CVodeFree(&cvode_mem); cvode_mem = nullptr; }}
    if (LS) {{ SUNLinSolFree(LS); LS = nullptr; }}
    if (A)  {{ SUNMatDestroy(A); A = nullptr; }}
    if (y)  {{ N_VDestroy(y); y = nullptr; }}
    if (ctx) {{ SUNContext_Free(&ctx); ctx = nullptr; }}
  }};

  // SUNContext_Create signature changed in SUNDIALS 7: void* -> SUNComm.
  // SUN_COMM_NULL is defined in sundials_types.h on v7 only.
#if SUNDIALS_VERSION_MAJOR >= 7
  if (SUNContext_Create(SUN_COMM_NULL, &ctx) < 0)
#else
  if (SUNContext_Create(nullptr, &ctx) < 0)
#endif
    return res.fail(cppde::RC_NO_MALLOC, "SUNContext_Create failed");

  y = N_VNew_Serial(NEQ, ctx);
  if (!y) {{ cleanup(); return res.fail(cppde::RC_NO_MALLOC, "N_VNew_Serial failed"); }}
  {{
    double* y0 = N_VGetArrayPointer(y);
    for (int i = 0; i < NEQ; ++i) y0[i] = ud.params[i];
  }}

  cvode_mem = CVodeCreate({cv_method}, ctx);
  if (!cvode_mem) {{ cleanup(); return res.fail(cppde::RC_NO_MALLOC, "CVodeCreate failed"); }}
  if (CVodeInit(cvode_mem, rhs_fn, t0, y) < 0) {{ cleanup(); return res.fail(cppde::RC_LINIT_FAIL, "CVodeInit failed"); }}
  if (CVodeSStolerances(cvode_mem, reltol, abstol) < 0) {{ cleanup(); return res.fail(cppde::RC_LINIT_FAIL, "CVodeSStolerances failed"); }}
  if (CVodeSetUserData(cvode_mem, &ud) < 0) {{ cleanup(); return res.fail(cppde::RC_LINIT_FAIL, "CVodeSetUserData failed"); }}
  CVodeSetMaxNumSteps(cvode_mem, maxsteps);
  if (hini > 0.0) CVodeSetInitStep(cvode_mem, hini);
{ls_setup}
{rootfunc_init_block}
  // --- sensitivities ---
{sens_init_block}
{asa_init_block}

  // --- write t0 row ---
  out_t.push_back(t0);
  {{
    const double* y0 = N_VGetArrayPointer(y);
    for (int i = 0; i < NEQ; ++i) out_y.push_back(y0[i]);
{sens_t0_block}
  }}

  // --- time-step loop ---
  std::string solver_msg;
  int return_code = 0;
  bool hard_fail = false;   // set by event_sens_reinit, checked at both splice sites
  double t_reached = t0;
  bool root_terminated = false; (void)root_terminated;

  // --- cumulative step counters across CVodeReInit, which zeroes them ---
  // `cv_counters` holds a per-counter accumulator plus the raw value the
  // counter had right after the last re-init; `harvest` banks the delta
  // since that baseline, `rebase` records a new one.
  struct cv_counters {{
    long nst = 0, nfe = 0, nje = 0, netf = 0, nsetups = 0, nfeS = 0;
    long b_nst = 0, b_nfe = 0, b_nje = 0, b_netf = 0, b_nsetups = 0, b_nfeS = 0;
  }} cvc;
  auto cv_raw = [&](long& nst, long& nfe, long& nje, long& netf,
                    long& nsetups, long& nfeS) {{
    nst = nfe = nje = netf = nsetups = nfeS = 0;
    CVodeGetNumSteps(cvode_mem, &nst);
    CVodeGetNumRhsEvals(cvode_mem, &nfe);
    CVodeGetNumJacEvals(cvode_mem, &nje);
    CVodeGetNumErrTestFails(cvode_mem, &netf);
    CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
    if (Ns_active > 0) CVodeGetSensNumRhsEvals(cvode_mem, &nfeS);
  }};
  auto cv_harvest = [&]() {{
    long a, b, c, d, e, f;
    cv_raw(a, b, c, d, e, f);
    cvc.nst     += a - cvc.b_nst;      cvc.nfe  += b - cvc.b_nfe;
    cvc.nje     += c - cvc.b_nje;      cvc.netf += d - cvc.b_netf;
    cvc.nsetups += e - cvc.b_nsetups;  cvc.nfeS += f - cvc.b_nfeS;
  }};
  auto cv_rebase = [&]() {{
    cv_raw(cvc.b_nst, cvc.b_nfe, cvc.b_nje, cvc.b_netf,
           cvc.b_nsetups, cvc.b_nfeS);
  }};
  (void)cv_harvest; (void)cv_rebase;

{asa_restart_lambda}{event_builder_block}{event_apply_lambda}{do_cvode_step_lambda}{event_pre_t0_block}{main_loop_body}

{asa_sweep_block}
  // --- diagnostics ---
  long n_steps = 0, n_fe = 0, n_je = 0, n_etf = 0, n_setups = 0;
  int  last_order = 0;
  double last_h = 0.0;
  cv_harvest();
  n_steps = cvc.nst; n_fe = cvc.nfe; n_je = cvc.nje;
  n_etf = cvc.netf;  n_setups = cvc.nsetups;
  CVodeGetLastOrder(cvode_mem, &last_order);
  CVodeGetLastStep(cvode_mem, &last_h);
{sens_fevals_block}

  const int n_out = static_cast<int>(out_t.size());

  // --- fill solve_result (plain double, no R API) ---
  // Internal buffers are state-first (out_y: [NEQ, n_out],
  // out_s: [NEQ, Ns_active, n_out]); the R-facing layout is time-first.
  res.prepare_out(n_out, deriv ? Ns_active : 0, NEQ, deriv, false);
  std::memcpy(res.t_out, out_t.data(), sizeof(double) * (size_t)n_out);
  for (int i = 0; i < n_out; ++i)
    for (int j = 0; j < NEQ; ++j)
      res.v_out[i + (size_t)n_out * j] = out_y[(size_t)i * NEQ + j];

  if (deriv && res.s1_out) {{
    if (out_s.empty()) {{
      std::fill(res.s1_out, res.s1_out + (size_t)n_out * NEQ * Ns_active, 0.0);
    }} else {{
      for (int i = 0; i < n_out; ++i)
        for (int iS = 0; iS < Ns_active; ++iS)
          for (int j = 0; j < NEQ; ++j)
            res.s1_out[i + (size_t)n_out * (j + (size_t)NEQ * iS)] =
              out_s[((size_t)i * Ns_active + iS) * NEQ + j];
    }}
  }}

  res.return_code = return_code;
  res.message     = solver_msg;
  res.accepted    = (int)(n_steps - n_etf);
  res.rejected    = (int)n_etf;
  res.fevals      = (int)n_fe;
  res.jevals      = (int)n_je;
  res.setups      = (int)n_setups;
  res.last_dt     = last_h;
  res.last_order  = last_order;
  res.t_reached   = t_reached;

  cleanup();
  return res.return_code;

}} catch (const std::exception& e) {{
  return res.fail(cppde::RC_UNRECOGNIZED_ERR,
                  std::string("CVODE solver failed: ") + e.what());
}} catch (...) {{
  return res.fail(cppde::RC_UNRECOGNIZED_ERR,
                  "CVODE solver failed: unknown C++ exception");
}}
}}

extern "C" SEXP solve_{modelname}(
    SEXP timesSEXP, SEXP paramsSEXP, SEXP sens1iniSEXP, SEXP sens2iniSEXP,
    SEXP fixedSEXP, SEXP abstolSEXP, SEXP reltolSEXP, SEXP maxprogressSEXP,
    SEXP maxstepsSEXP, SEXP hiniSEXP, SEXP root_tolSEXP, SEXP maxrootSEXP,
    SEXP forcingTimesSEXP, SEXP forcingValuesSEXP, SEXP seedSEXP,
    SEXP dimnamesSEXP)
{{
  cppde::rbatch::solve_args a = cppde::rbatch::read_solve_args(
      timesSEXP, paramsSEXP, sens1iniSEXP, sens2iniSEXP, fixedSEXP,
      abstolSEXP, reltolSEXP, maxprogressSEXP, maxstepsSEXP, hiniSEXP,
      root_tolSEXP, maxrootSEXP, forcingTimesSEXP, forcingValuesSEXP, seedSEXP);
  // Nothing is protected yet, so the longjmp is safe here.
{seed_guard}
  return cppde::rbatch::solve_one(a, NEQ, {deriv_flag}, false, &solve_impl, dimnamesSEXP);
}}

{cv_event_times_fn}extern "C" SEXP solve_{modelname}_batch(SEXP condsSEXP, SEXP nthreadsSEXP, SEXP dimnamesSEXP)
{{
  const int K  = Rf_length(condsSEXP);
  const int nt = Rf_isNull(nthreadsSEXP) ? 0 : INTEGER(nthreadsSEXP)[0];

  // Phase A: all R reads, serially.
  std::vector<cppde::rbatch::solve_args> A;
  A.reserve(K);
  for (int k = 0; k < K; ++k)
    A.push_back(cppde::rbatch::read_cond_args(VECTOR_ELT(condsSEXP, k)));
  std::vector<cppde::rbatch::solve_result> R_(K);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, K));
{cv_prealloc_block}
  // Phase B: no R API.  Each condition builds its own SUNContext.
  cppde::rbatch::run_batch(K, nt, [&](int k) noexcept {{ solve_impl(A[k], R_[k]); }});

  // Phase C: back on the R side.
  for (int k = 0; k < K; ++k) {{
{cv_finish_block}    SEXP e = PROTECT(cppde::rbatch::build_result_sexp(R_[k], NEQ, {deriv_flag}, false));
    SET_VECTOR_ELT(out, k, e);
    UNPROTECT(1);
  }}
  UNPROTECT(1);
  return out;
}}
"""
