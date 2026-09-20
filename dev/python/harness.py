"""Harness of the graph-AD code generators: Python emission, sampling,
comparison, the legacy SymPy reference (legacy/, e55c599) and the model
collection."""

import importlib.util
import math
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PYDIR = os.path.join(os.path.dirname(os.path.dirname(HERE)), "inst", "python")
for _p in (PYDIR, HERE):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import numpy as np  # noqa: E402
import sympy as sp  # noqa: E402

import cppde_emit as em  # noqa: E402
import cppde_graph as cg  # noqa: E402
import cppde_model  # noqa: E402
import pyrt  # noqa: E402

LEGACY = os.path.join(HERE, "legacy")


# ===========================================================================
# Emission, sampling, comparison
# ===========================================================================


TOL = 1e-13


def as_list(x):
    if x is None:
        return []
    if isinstance(x, (list, tuple)):
        return list(x)
    return [x]


# ---------------------------------------------------------------------------
# Models
# ---------------------------------------------------------------------------

OdeModel = cppde_model.OdeModel


def py_function(g, stores, slot, args, prefix="_t", linmap=None):
    """Compile stores into a Python function def f(*args)."""
    stmts, names = em.schedule(g, stores, prefix=prefix)
    pr = em.Printer(g, slot, style="py", names=names)
    body = em.render_py(stmts, pr)
    src = "def f(%s):\n%s\n    return None\n" % (", ".join(args), "\n".join(body))
    scope = {"R": pyrt, "linmap_": linmap}
    exec(compile(src, "<generated>", "exec"), scope)
    f = scope["f"]
    f.source = src
    f.stmts = stmts
    return f


# ---------------------------------------------------------------------------
# Sampling
# ---------------------------------------------------------------------------

def draw(rng, kind):
    if kind == 0:
        return rng.uniform(0.2, 3.0)
    if kind == 1:
        return rng.uniform(-3.0, 3.0)
    if kind == 2:
        return math.exp(rng.uniform(math.log(1e-3), math.log(1e3)))
    return rng.choice((0.0, 1.0, -1.0, 0.5, 2.0))


class Point:
    """One evaluation point of an ODE model."""

    def __init__(self, model, rng, kind):
        n = model.n
        self.x = [draw(rng, kind) for _ in range(n)]
        self.p = [draw(rng, kind) for _ in range(n + len(model.params))]
        self.t = draw(rng, kind) * (10 if kind != 3 else 1)
        self.fv = [draw(rng, kind) for _ in model.forcings]
        self.fr = [draw(rng, 1) for _ in model.forcings]
        self.F = [pyrt.Forcing((lambda v: (lambda t: v))(v),
                               (lambda r: (lambda t: r))(r))
                  for v, r in zip(self.fv, self.fr)]
        self.lam = [draw(rng, 1) for _ in range(n)]
        self.v = [draw(rng, 1) for _ in range(n)]
        self.sc = draw(rng, 0)

    def env(self, model):
        """name -> float, as the legacy printer reads the names."""
        e = {}
        n = model.n
        for j, f in enumerate(model.forcings):
            e[f] = self.fv[j]
        for k, p in enumerate(model.params):
            e[p] = self.p[n + k]
        for i, s in enumerate(model.states):
            e[s] = self.x[i]
        for name, i in model.init_names.items():
            e[name] = self.p[i]
        e["time"] = self.t
        for i in range(n):
            e["_lam%d" % i] = self.lam[i]
            e["_v%d" % i] = self.v[i]
        for j in range(len(model.forcings)):
            e["_fd%d" % j] = self.fr[j]
        return e

    def leaf_env(self, model):
        """leaf attr -> float, for cppde_graph.Graph.evaluate."""
        e = {}
        for i in range(model.n):
            e[(cg.STATE, i)] = self.x[i]
            e[(cg.INIT, i)] = self.p[i]
        for k in range(len(model.params)):
            e[(cg.PARAM, k)] = self.p[model.n + k]
        for j in range(len(model.forcings)):
            e[(cg.FORCING, j)] = self.fv[j]
            e[(cg.FRATE, j)] = self.fr[j]
        e[(cg.TIME,)] = self.t
        for r in range(model.nlin):
            e[(cg.LINROW, r)] = pyrt.row_dot(model.linmap, r, self.x)
        return e


# Intermediate magnitude beyond which the last digits are cancellation noise.
ILL = 1e12


def ill_conditioned(model, pt):
    """True if the right-hand side is not finite at pt or a finite
    intermediate value exceeds ILL in magnitude."""
    g = model.g
    env = pt.leaf_env(model)
    val = {}
    for n in g.topo(model.rhs_plain):
        try:
            v = g._eval_node(n, val, env)
        except (ArithmeticError, ValueError):
            return True
        val[n] = v
        if not isinstance(v, bool) and math.isfinite(v) and abs(v) > ILL:
            return True
    return not all(math.isfinite(val[r]) for r in model.rhs_plain)


def conditions(g, roots):
    return [n for n in g.topo(roots) if g.op[n] == cg.CMP]


def sample_points(model, roots, rng, base=8, extra=48):
    """Points covering both outcomes of every reachable condition."""
    conds = conditions(model.g, roots)
    seen = {c: set() for c in conds}
    pts = []
    for k in range(base + extra):
        pt = Point(model, rng, k % 4)
        if k >= base:
            vals = model.g.evaluate(conds, pt.leaf_env(model)) if conds else []
            new = any(bool(v) not in seen[c] for c, v in zip(conds, vals))
            if not new:
                continue
        else:
            vals = model.g.evaluate(conds, pt.leaf_env(model)) if conds else []
        for c, v in zip(conds, vals):
            seen[c].add(bool(v))
        pts.append(pt)
        if k >= base and all(len(s) == 2 for s in seen.values()):
            break
    covered = sum(len(s) == 2 for s in seen.values())
    return pts, (covered, len(conds))


# ---------------------------------------------------------------------------
# Comparison
# ---------------------------------------------------------------------------

def close(a, b, tol=TOL):
    a = float(a)
    b = float(b)
    if math.isnan(a) or math.isnan(b):
        return math.isnan(a) and math.isnan(b)
    if math.isinf(a) or math.isinf(b):
        return a == b
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


class Tally:
    """Counts of comparisons and the first failures."""

    def __init__(self, name):
        self.name = name
        self.n = 0
        self.bad = 0
        self.examples = []
        self.skipped = []
        self.singular = 0

    def check(self, ok, what):
        self.n += 1
        if not ok:
            self.bad += 1
            if len(self.examples) < 12:
                self.examples.append(what)

    def skip(self, what):
        self.skipped.append(what)

    def compare(self, a, b, tol, what, singular=False):
        """Check a against reference b; a non-finite reference or a
        `singular` point is counted, not compared."""
        if singular or not math.isfinite(float(b)):
            self.singular += 1
            return
        self.check(close(a, b, tol), what)

    def report(self):
        status = "ok" if self.bad == 0 else "FAIL"
        print("%-44s %6d checks  %4d failed  %s" % (self.name, self.n, self.bad, status))
        if self.singular:
            print("    %d comparisons at a singular or ill-conditioned point"
                  % self.singular)
        for e in self.examples:
            print("    " + e)
        if self.skipped:
            print("    skipped %d: %s" % (len(self.skipped), "; ".join(self.skipped[:5])))
        return self.bad == 0


# ===========================================================================
# Legacy SymPy reference
# ===========================================================================

LEGACY = os.path.join(HERE, "legacy")


def _load(name, alias):
    spec = importlib.util.spec_from_file_location(
        alias, os.path.join(LEGACY, name + ".py"))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[alias] = mod
    spec.loader.exec_module(mod)
    return mod


def _load_all():
    keep = {k: sys.modules.get(k) for k in ("cppsympy", "codegen_cppODE")}
    try:
        sys.modules["cppsympy"] = _load("cppsympy", "_legacy_cppsympy")
        ode = _load("codegen_cppODE", "_legacy_codegen_cppODE")
        sys.modules["codegen_cppODE"] = ode
        cvode = _load("codegen_cvode", "_legacy_codegen_cvode")
        fun = _load("codegen_cppFUN", "_legacy_codegen_cppFUN")
    finally:
        for k, v in keep.items():
            if v is None:
                sys.modules.pop(k, None)
            else:
                sys.modules[k] = v
    return ode, cvode, fun



class _Legacy:
    """Legacy modules, loaded on first use."""

    def __getattr__(self, name):
        mods = dict(zip(("ODE", "CVODE", "FUN"), _load_all()))
        self.__dict__.update(mods)
        return mods[name]


legacy = _Legacy()



def replace_dirac(e):
    """DiracDelta as the discrete indicator; d floor and d ceiling as 0."""
    e = legacy.ODE._replace_dirac_delta(e)
    if isinstance(e, sp.Basic) and e.has(sp.Derivative):
        e = e.replace(lambda d: isinstance(d, sp.Derivative)
                      and isinstance(d.expr, (sp.floor, sp.ceiling)),
                      lambda d: sp.Integer(0))
    return e


def sym_table(names):
    return {n: sp.Symbol(n, real=True) for n in names}


def legacy_parse(text, local):
    return legacy.ODE._safe_sympify(str(text), local)


def ref_jacobian(exprs, syms):
    J = legacy.ODE._compute_ode_jacobian_serial(exprs, syms, set(syms))
    return [[replace_dirac(e) for e in row] for row in J]


def ref_dfdp(exprs, syms):
    """[(row, k, expr)] of the nonzero d exprs[row] / d syms[k]."""
    return [(i, k, replace_dirac(e))
            for i, k, e in legacy.ODE._compute_ode_dfdp(exprs, syms, set(syms))]


def ref_diff(expr, sym):
    return replace_dirac(sp.diff(expr, sym))


def ref_dfdt(exprs, t, forcings, rates):
    """df/dt with the forcing chain terms; `rates[j]` stands for F_j'(t)."""
    out = []
    for e in exprs:
        d = replace_dirac(sp.diff(e, t))
        for fsym, rsym in zip(forcings, rates):
            dF = sp.diff(e, fsym)
            if dF != 0:
                d = d + dF * rsym
        out.append(d)
    return out


# ---------------------------------------------------------------------------
# Numeric evaluation
# ---------------------------------------------------------------------------

def _heaviside(x, h0=0.5):
    return 0.0 if x < 0 else (h0 if x == 0 else 1.0)


def _sign(x):
    return 1.0 if x > 0 else (-1.0 if x < 0 else 0.0)


def _dirac(x, *k):
    return 1.0 if x == 0 else 0.0


NUMERIC = {
    "Heaviside": _heaviside, "sign": _sign, "DiracDelta": _dirac,
    "Min": min, "Max": max, "Abs": abs, "floor": math.floor,
    "ceiling": math.ceil, "erf": math.erf, "erfc": math.erfc,
    "gamma": math.gamma, "loggamma": math.lgamma,
}


class Evaluator:
    """IEEE float evaluation of SymPy expressions over named arguments."""

    def __init__(self, exprs, names):
        self.names = list(names)
        syms = [sp.Symbol(n, real=True) for n in self.names]
        self.fns = [sp.lambdify(syms, e, modules=[NUMERIC, "numpy"])
                    for e in exprs]

    def __call__(self, env):
        args = [np.float64(env[n]) for n in self.names]
        out = []
        with np.errstate(all="ignore"):
            for f in self.fns:
                try:
                    v = float(f(*args))
                except (ValueError, ZeroDivisionError, TypeError):
                    v = math.nan
                except OverflowError:
                    v = math.inf
                out.append(v)
        return out


# ===========================================================================
# Model collection
# ===========================================================================

def chain(n, k="k", sparse=True):
    rhs = {}
    for i in range(n):
        prev = "x%d" % (i - 1)
        cur = "x%d" % i
        if i == 0:
            rhs[cur] = "-%s1*%s + u*%s0/(1 + %s)" % (k, cur, k, cur)
        else:
            rhs[cur] = "%s%d*%s - %s%d*%s^2" % (k, i, prev, k, i + 1, cur)
    return rhs


def params_of(rhs, forcings=(), states=None):
    import re
    states = set(rhs) if states is None else set(states)
    names = []
    for e in rhs.values():
        for m in re.finditer(r"[A-Za-z_][A-Za-z0-9_]*", str(e)):
            s = m.group(0)
            if s in states or s in forcings or s == "time" or s in names:
                continue
            if s in ("exp", "log", "sqrt", "abs", "sign", "min", "max",
                     "piecewise", "Heaviside", "sin", "cos", "tanh",
                     "floor", "ceiling", "pow", "exp10", "log10"):
                continue
            names.append(s)
    return names


def ode_spec(name, rhs, forcings=(), events=None):
    return {"name": name, "rhs_dict": rhs, "params_list": params_of(rhs, forcings),
            "forcings_list": list(forcings), "events": events}


MODELS = [
    # S1: three states, forcing, time-dependent piecewise
    ode_spec("S1", {
        "A": "-k1*A*B + k2*C - piecewise(kd*A, time - ts > 0, 0.1*kd*A)",
        "B": "-k1*A*B + k2*C + u*kin",
        "C": "k1*A*B - k2*C - exp(-time/tau)*C",
    }, forcings=["u"]),
    # S2: 40-state chain, sparse
    ode_spec("S2", chain(40)),
    ode_spec("switches", {
        "x": "-abs(x - a)*k + sign(y)*b - min(x, y)*c + max(x, 2*y)*d",
        "y": "Heaviside(x - a)*k - y*floor(c) + ceiling(d)*x",
    }),
    ode_spec("init_values", {
        "A": "-k*A + k*A_0*exp(-time)",
        "B": "k*A - B/(1 + B_0^2)",
    }),
    ode_spec("no_params", {"x": "-x^3 + y", "y": "-y + x*sin(time)"}),
    ode_spec("explicit_time", {
        "x": "-a*x + b*time^2 - c*cos(w*time)*x",
        "y": "x*exp(-time*d) - y*log(1 + time)",
    }),
    ode_spec("forcings", {
        "x": "-k*x*F1 + F2^2",
        "y": "k*x*F1 - log(F2)*y",
    }, forcings=["F1", "F2"]),
    ode_spec("logic", {
        "x": "-k*x*piecewise(1, time > t1 && time <= t2, 0)",
        "y": "piecewise(k*x, x > y || !(time > t1), -k*y, x < 0, 0)",
    }),
    ode_spec("rational", {
        "S": "-vmax*S/(Km + S) + kin",
        "P": "vmax*S/(Km + S) - kout*P^n/(Kd^n + P^n)",
    }),
    ode_spec("powers", {
        "x": "sqrt(x)*a - x^(1/3)*b + x^(-2)*c + exp10(y)*d",
        "y": "log10(x)*a - y^p + 2^y",
    }),
]


def event_spec(name, rhs, events, forcings=()):
    texts = dict(rhs)
    for col in ("value", "time", "root"):
        for k, v in enumerate(events.get(col, [])):
            if v is not None:
                texts["_%s%d" % (col, k)] = str(v)
    params = params_of(texts, forcings, states=list(rhs))
    return {"name": name, "rhs_dict": rhs, "states_list": list(rhs),
            "params_list": params, "forcings_list": list(forcings),
            "events_df": events}


EVENT_MODELS = [
    # S3: non-terminal, time-dependent nonlinear root, fixed event at a
    # parameter-valued time
    event_spec("S3", {"A": "-k1*A + k2*B*time", "B": "k1*A - k2*B^2"}, {
        "var": ["A", "B"],
        "value": ["A + dose*B", "B*fac - A_0"],
        "time": ["t_dose*2", None],
        "root": [None, "A^2 - thr*B - 0.1*time"],
        "method": ["add", "replace"],
    }),
    event_spec("forced_root", {"x": "-k*x + u", "y": "x - y"}, {
        "var": ["y"], "value": ["y + u*x"], "time": [None],
        "root": ["x*u - thr"], "method": ["replace"],
    }, forcings=["u"]),
    event_spec("heaviside_root", {"x": "-k*x", "y": "x*Heaviside(x - a)"}, {
        "var": ["x"], "value": ["x + 1"], "time": [None],
        "root": ["Heaviside(y - a) - 0.5 + x*b"], "method": ["add"],
    }),
]


FUN_MODELS = [
    ("S4", {"y1": "piecewise(a*x^2, x > c, b*sqrt(x))", "y2": "abs(x - a)*pow(b, 2.5)",
            "y3": "exp(-k*x) + max(x, a)"}, ["x"], ["a", "b", "c", "k"]),
    ("no_variables", {"A": "exp(la)", "B": "la^2 + lb", "C": "1.5"}, [], ["la", "lb"]),
    ("no_parameters", {"y": "x^2 + 3*x", "z": "log(1 + x*w)"}, ["x", "w"], []),
    ("constants", {"A": "1.0", "B": "2.5"}, [], ["dummy"]),
]


# ---------------------------------------------------------------------------
# Models with long linear sums (E)
# ---------------------------------------------------------------------------

def _fmt(c):
    return repr(float(c))


def llg_dipole(n_spins, drive=True):
    """Landau-Lifshitz-Gilbert macrospins on a helix with dense dipolar
    coupling: 3 * n_spins states, plus a drive phase `phi` if `drive`.

    The field of spin s is H_a + Ms * sum_t D_st m_t with numeric D; each
    field component appears once per equation.
    """
    pos = [(math.cos(2 * math.pi * s / n_spins), math.sin(2 * math.pi * s / n_spins),
            0.35 * s / n_spins) for s in range(n_spins)]
    comp = "xyz"
    names = [["m%s%d" % (a, s) for a in comp] for s in range(n_spins)]

    def tensor(s, t):
        if s == t:
            return [[-0.2 if a == b else 0.0 for b in range(3)] for a in range(3)]
        r = [pos[t][k] - pos[s][k] for k in range(3)]
        d = math.sqrt(sum(v * v for v in r))
        u = [v / d for v in r]
        return [[(3 * u[a] * u[b] - (1.0 if a == b else 0.0)) * 0.01 / d ** 3
                 for b in range(3)] for a in range(3)]

    D = [[tensor(s, t) for t in range(n_spins)] for s in range(n_spins)]
    ext = {"x": "Hx0 + h1*cos(phi)" if drive else "Hx0", "y": "Hy0", "z": "Hz0"}
    rhs = {}
    for s in range(n_spins):
        H = []
        for a in range(3):
            terms = []
            for t in range(n_spins):
                for b in range(3):
                    c = D[s][t][a][b]
                    if c != 0.0:
                        terms.append("%s*%s" % (_fmt(c), names[t][b]))
            H.append("(%s + Ms*(%s))" % (ext[comp[a]], " + ".join(terms)))
        mx, my, mz = names[s]
        Hx, Hy, Hz = H
        rhs[mx] = ("-gp*(%s*(-alpha*(%s^2 + %s^2)) + %s*(alpha*%s*%s - %s) + %s*(%s + alpha*%s*%s))"
                   % (Hx, my, mz, Hy, mx, my, mz, Hz, my, mx, mz))
        rhs[my] = ("-gp*(%s*(%s + alpha*%s*%s) + %s*(-alpha*(%s^2 + %s^2)) + %s*(alpha*%s*%s - %s))"
                   % (Hx, mz, mx, my, Hy, mx, mz, Hz, my, mz, mx))
        rhs[mz] = ("-gp*(%s*(alpha*%s*%s - %s) + %s*(%s + alpha*%s*%s) + %s*(-alpha*(%s^2 + %s^2)))"
                   % (Hx, mx, mz, my, Hy, mx, my, mz, Hz, mx, my))
    params = ["gp", "alpha", "Ms", "Hx0", "Hy0", "Hz0"]
    if drive:
        rhs["phi"] = "omega"
        params += ["h1", "omega"]
    return rhs, params


def coupled_chain(n):
    """Dense linear coupling with parameter, time and forcing factors, a
    repeated row, a sum inside a condition and a nonlinear remainder."""
    rhs = {}
    for i in range(n):
        row = " + ".join("%s*x%d" % (_fmt(0.1 * ((i * 7 + j * 3) % 5 - 2) + 0.05), j)
                         for j in range(n))
        alt = " - ".join("%s*x%d" % (_fmt(0.5 + 0.01 * j), j) for j in range(n))
        rep = " + ".join("%s*x%d" % (_fmt(1.0 / (j + 1)), j) for j in range(n))
        rhs["x%d" % i] = (
            "-k*x%d + g*(%s) + exp(-time/tau)*(%s) + u*(%s) - x%d^2*(%s)"
            " + piecewise(a, %s > c, b) + (%s + x%d*x%d)"
            % (i, row, alt, rep, i, rep, "(" + rep + ")", row, i, (i + 1) % n))
    return rhs, ["k", "g", "tau", "a", "b", "c"], ["u"]


def _linear_models():
    out = [(label, rhs, params, []) for label, (rhs, params) in (
        ("llg4", llg_dipole(4)), ("llg3-static", llg_dipole(3, drive=False)))]
    rhs, params, forcings = coupled_chain(5)
    return out + [("coupled5", rhs, params, forcings)]


LINEAR_MODELS = _linear_models()


# ---------------------------------------------------------------------------
# Models of regular structure (D)
# ---------------------------------------------------------------------------

def rd1d(n, periodic=False):
    """Fisher-KPP on a line; Neumann ends unless periodic."""
    rhs = {}
    for i in range(n):
        left, right = i - 1, i + 1
        if periodic:
            left, right = left % n, right % n
        if left < 0:
            lap = "u%d - u%d" % (right, i)
        elif right >= n:
            lap = "u%d - u%d" % (left, i)
        else:
            lap = "u%d - 2*u%d + u%d" % (left, i, right)
        rhs["u%d" % i] = "D*(%s) + k*u%d*(1 - u%d/K)" % (lap, i, i)
    return rhs, ["D", "k", "K"]


def brusselator2d(nx, ny):
    """Brusselator on a grid, 5-point Laplacian with Neumann boundaries."""
    rhs = {}
    for s in ("u", "v"):
        for i in range(nx):
            for j in range(ny):
                nb = [(i + di, j + dj) for di, dj in ((-1, 0), (1, 0), (0, -1), (0, 1))
                      if 0 <= i + di < nx and 0 <= j + dj < ny]
                c = "%s%d_%d" % (s, i, j)
                lap = " + ".join("%s%d_%d" % (s, a, b) for a, b in nb) + " - %d*%s" % (len(nb), c)
                u, v = "u%d_%d" % (i, j), "v%d_%d" % (i, j)
                if s == "u":
                    rhs[c] = "Du*(%s) + a - (b + 1)*%s + %s^2*%s" % (lap, u, u, v)
                else:
                    rhs[c] = "Dv*(%s) + b*%s - %s^2*%s" % (lap, u, u, v)
    return rhs, ["Du", "Dv", "a", "b"]


def kuramoto(n):
    """Kuramoto oscillators with individual frequencies and a forcing."""
    rhs = {}
    for i in range(n):
        s = " + ".join("sin(th%d - th%d)" % (j, i) for j in range(n) if j != i)
        rhs["th%d" % i] = "w%d + K/%d*(%s) + drive*sin(time - th%d)" % (i, n, s, i)
    return rhs, ["w%d" % i for i in range(n)] + ["K"], ["drive"]


def _vector_models():
    out = []
    for label, (rhs, params) in (("rd1d", rd1d(40)), ("rd1d-periodic", rd1d(40, True)),
                                 ("brusselator2d", brusselator2d(8, 8)),
                                 ("llg20", llg_dipole(20))):
        out.append((label, rhs, params, []))
    rhs, params, forcings = kuramoto(20)
    out.append(("kuramoto", rhs, params, forcings))
    return out


VECTOR_MODELS = _vector_models()
