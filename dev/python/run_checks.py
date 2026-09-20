"""Checks of the graph-AD code generators against the legacy SymPy ones.

Usage:
    uv run --no-project --python 3.12 --with sympy==1.14 \
        python dev/python/run_checks.py [check ...]

Without arguments every check runs. Exit status 1 if any check fails.
"""

import json
import os
import random
import subprocess
import sys
import math
import time

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import harness as H  # noqa: E402

import cppde_graph as cg  # noqa: E402

CORPUS = os.path.join(HERE, "corpus.json")


def load_corpus():
    """corpus.json plus the models of models.py."""
    corpus = {}
    if os.path.exists(CORPUS):
        with open(CORPUS, encoding="utf-8") as f:
            corpus = json.load(f)
    corpus["ode"] = [dict(m) for m in H.MODELS] + corpus.get("ode", [])
    return corpus


# ---------------------------------------------------------------------------
# Expression sets of the corpus
# ---------------------------------------------------------------------------

def expression_sets(corpus):
    """(label, OdeModel-like spec, [expr strings], legacy parser) per entry."""
    out = []
    for i, e in enumerate(corpus.get("ode", [])):
        out.append(("ode[%d]" % i, e["rhs_dict"], e.get("params_list"),
                    e.get("forcings_list"), list(e["rhs_dict"].values()), "ode"))
    for i, e in enumerate(corpus.get("cvode", [])):
        out.append(("cvode[%d]" % i, e["rhs_dict"], e.get("params_list"),
                    e.get("forcings_list"), list(e["rhs_dict"].values()), "ode"))
    for i, e in enumerate(corpus.get("event", [])):
        ev = e.get("events_df") or {}
        exprs = []
        for col in ("value", "time", "root"):
            for v in H.as_list(ev.get(col)):
                if v is None or isinstance(v, bool):
                    continue
                s = str(v).strip()
                if s.lower() in ("", "na", "nan", "none", "true", "false"):
                    continue
                exprs.append(s)
        states = H.as_list(e.get("states_list"))
        rhs = e.get("rhs_dict") or {s: "0" for s in states}
        out.append(("event[%d]" % i, rhs, e.get("params_list"),
                    e.get("forcings_list"), exprs, "ode"))
    for i, e in enumerate(corpus.get("rootfunc", [])):
        rf = H.as_list(e.get("rootfunc"))
        rf = [r for r in rf if str(r).strip().lower() != "equilibrate"]
        states = H.as_list(e.get("states_list"))
        out.append(("rootfunc[%d]" % i, {s: "0" for s in states},
                    e.get("params_list"), e.get("forcings_list"), rf, "ode"))
    for i, e in enumerate(corpus.get("fun", [])):
        exprs = e.get("exprs")
        exprs = list(exprs.values()) if isinstance(exprs, dict) else H.as_list(exprs)
        names = H.as_list(e.get("variables")) + H.as_list(e.get("parameters"))
        out.append(("fun[%d]" % i, {}, names, [], exprs, "fun"))
    return out


def legacy_parse(kind, text, local, names):
    if kind == "fun":
        ctx = H.legacy.FUN.CodeGenContext([], names)
        return ctx._parse_expr(str(text))
    return H.legacy_parse(text, local)


# ---------------------------------------------------------------------------
# Checks
# ---------------------------------------------------------------------------

# Term order of long sums with cancellation moves the last digits.
PARSE_TOL = 1e-12


def check_parse(corpus):
    """Parser and Python backend against lambdify on the corpus."""
    t = H.Tally("parse: corpus values")
    rng = random.Random(1)
    cov = [0, 0]
    for label, rhs, params, forcings, exprs, kind in expression_sets(corpus):
        if not exprs:
            continue
        try:
            m = H.OdeModel(rhs, params, forcings)
            roots = [m.parser.parse(x) for x in exprs]
        except Exception as e:
            try:
                local = H.sym_table(H.OdeModel({}, params, forcings).names()
                                    + list(rhs))
                [legacy_parse(kind, x, local, H.as_list(params)) for x in exprs]
            except Exception:
                t.skip("%s: rejected by both parsers" % label)
                continue
            t.check(False, "%s: new parser: %s" % (label, e))
            continue
        local = H.sym_table(m.names())
        try:
            ref = [legacy_parse(kind, x, local, H.as_list(params)) for x in exprs]
        except Exception as e:
            t.skip("%s: legacy parser: %s" % (label, str(e)[:60]))
            continue
        ev = H.Evaluator(ref, m.names())
        stores = [(("vec", "out", i), r, "=") for i, r in enumerate(roots)]
        f = H.py_function(m.g, stores, m.slot("py"), ["x", "params", "t", "F", "out"])
        pts, c = H.sample_points(m, roots, rng)
        cov[0] += c[0]
        cov[1] += c[1]
        for pt in pts:
            out = [None] * len(roots)
            f(pt.x, pt.p, pt.t, pt.F, out)
            want = ev(pt.env(m))
            for k, (a, b) in enumerate(zip(out, want)):
                t.compare(a, b, PARSE_TOL, "%s expr %d %r: new %r legacy %r"
                          % (label, k, exprs[k][:60], a, b))
    print("  branch coverage: %d of %d conditions seen both ways" % tuple(cov))
    return t.report()


def check_fallback(corpus):
    """Expressions only SymPy reads."""
    t = H.Tally("parse: SymPy fallback")
    m = H.OdeModel({"x": "0"}, ["a", "b"], [])
    cases = [("x! + a", lambda x, a, b: __import__("math").gamma(x + 1) + a),
             ("re(a) * x", lambda x, a, b: a * x),
             ("beta(a, b)", None),
             ("a * piecewise(1, x > 2, 0, x <= 2)", lambda x, a, b: a * (x > 2)),
             ("Piecewise((b, x < 2), (a, x >= 2))", lambda x, a, b: b if x < 2 else a)]
    for text, fn in cases:
        try:
            n = m.parser.parse(text)
        except Exception as e:
            t.check(False, "%s: %s" % (text, e))
            continue
        if fn is None:
            ok = m.g.op[n] == cg.CALL and m.g.attr[n] == "beta"
            t.check(ok, "%s: got %s" % (text, m.g.describe(n)))
            continue
        for xv in (2.5, 1.5):
            env = {(cg.STATE, 0): xv, (cg.PARAM, 0): 0.7, (cg.PARAM, 1): 1.3}
            v = m.g.evaluate([n], env)[0]
            t.check(H.close(v, fn(xv, 0.7, 1.3)), "%s: %r" % (text, v))
    for text in ("piecewise(1, x > 2, 0, x < 2)", "Piecewise((a, x > 2))"):
        try:
            m.parser.parse(text)
            t.check(False, "%s: accepted without a default branch" % text)
        except ValueError as e:
            t.check("default branch" in str(e), "%s: %s" % (text, e))
    return t.report()


def emit_sample_text(corpus):
    """C++ text of the corpus right-hand sides, for the determinism check."""
    import cppde_emit as em
    parts = []
    for label, rhs, params, forcings, exprs, kind in expression_sets(corpus)[:200]:
        if not exprs:
            continue
        try:
            m = H.OdeModel(rhs, params, forcings)
            roots = [m.parser.parse(x) for x in exprs]
        except Exception:
            continue
        stores = [(("vec", "dxdt", i), r, "=") for i, r in enumerate(roots)]
        stmts, names = em.schedule(m.g, stores)
        pr = em.Printer(m.g, m.slot("cpp"), style="ad", names=names)
        try:
            parts += em.render_cpp(stmts, pr, "cppde::dual<double, 0>")
        except em.EmitError:
            continue
    return "\n".join(parts)


def check_determinism(corpus):
    """Emitted text does not depend on PYTHONHASHSEED."""
    t = H.Tally("emit: PYTHONHASHSEED independence")
    code = ("import sys, hashlib; sys.path.insert(0, %r); import run_checks as R; "
            "print(hashlib.sha256(R.emit_sample_text(R.load_corpus()).encode()).hexdigest())"
            % HERE)
    digests = []
    for seed in ("1", "12345"):
        env = dict(os.environ, PYTHONHASHSEED=seed)
        r = subprocess.run([sys.executable, "-c", code], env=env,
                           capture_output=True, text=True)
        digests.append(r.stdout.strip() or r.stderr[-200:])
    t.check(digests[0] == digests[1] and len(digests[0]) == 64,
            "digests differ: %s" % digests)
    return t.report()


SCALARS = (("double", "double", 0),
           ("ad", "cppde::dual<double, 0>", 1),
           ("ad", "cppde::dual2nd<double, 0>", 2))

CXX_HEAD = """#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <vector>
#include <cmath>
#include <limits>
#include <cppde/cppde.hpp>
using namespace cppde;
"""


def include_dir():
    """cppDE headers: CPPDE_INCLUDE, else inst/include next to this harness."""
    return os.environ.get("CPPDE_INCLUDE") or os.path.join(
        os.path.dirname(os.path.dirname(HERE)), "inst", "include")


def tool_paths():
    """(g++, R include dir), from CXX / R_INCLUDE or the Rtools defaults."""
    cxx = os.environ.get("CXX") or r"C:\rtools45\x86_64-w64-mingw32.static.posix\bin\g++.exe"
    rinc = os.environ.get("R_INCLUDE") or r"C:\Program Files\R\R-4.6.1\include"
    return cxx, rinc


def syntax_units(corpus, limit=80):
    """C++ functions of the corpus right-hand sides in all three scalars."""
    import cppde_emit as em
    units = []
    for label, rhs, params, forcings, exprs, kind in expression_sets(corpus):
        if kind != "ode" or not exprs or len(units) >= limit:
            continue
        try:
            m = H.OdeModel(rhs, params, forcings)
            roots = [m.parser.parse(x) for x in exprs]
        except Exception:
            continue
        stores = [(("vec", "out", i), r, "=") for i, r in enumerate(roots)]
        stmts, names = em.schedule(m.g, stores)
        body = []
        for style, scalar, lvl in SCALARS:
            pr = em.Printer(m.g, m.slot("cpp"), style=style, names=names,
                            ad_level=lvl)
            try:
                lines = em.render_cpp(stmts, pr, scalar)
            except em.EmitError:
                continue
            fname = "f%d" % lvl
            body.append(
                "void %s(const std::vector<%s>& x, std::vector<%s>& out,\n"
                "        const %s& t, const std::vector<%s>& params,\n"
                "        const std::vector<const cppde::PchipForcing<%s>*>& F) {\n"
                "  (void)x; (void)t; (void)params; (void)F;\n%s\n}\n"
                % (fname, scalar, scalar, scalar, scalar, scalar, "\n".join(lines)))
        units.append("namespace u%d { // %s\n%s}\n" % (len(units), label, "".join(body)))
    return units


def check_syntax(corpus):
    """g++ -fsyntax-only over the emitted corpus in double, dual, dual2nd."""
    import tempfile
    t = H.Tally("emit: C++ syntax (double, dual, dual2nd)")
    cxx, rinc = tool_paths()
    if not os.path.exists(cxx):
        t.skip("no compiler at " + cxx)
        return t.report()
    inc = include_dir()
    units = syntax_units(corpus)
    d = tempfile.mkdtemp(prefix="cppde_syntax_")
    src = os.path.join(d, "corpus.cpp")
    with open(src, "w", encoding="utf-8") as f:
        f.write(CXX_HEAD + "\n".join(units))
    r = subprocess.run([cxx, "-std=gnu++17", "-fsyntax-only", "-I", inc,
                        "-I", rinc, src], capture_output=True, text=True)
    errs = [l for l in r.stderr.splitlines() if "error" in l]
    t.check(r.returncode == 0, "g++ exit %d, %d error lines: %s"
            % (r.returncode, len(errs), " | ".join(errs[:3])))
    print("  %d units in %s" % (len(units), src))
    return t.report()


def check_model_syntax(corpus):
    """g++ -fsyntax-only over whole generated models with maps and loops, in
    double, dual and dual2nd, forward (dense, sparse) and reverse (rb4)."""
    import tempfile
    import codegen_cppODE as CO
    t = H.Tally("emit: C++ syntax of models with maps and loops")
    cxx, rinc = tool_paths()
    if not os.path.exists(cxx):
        t.skip("no compiler at " + cxx)
        return t.report()
    inc = include_dir()
    rhs_f, params_f, forcings_f = H.kuramoto(20)
    models = [("rd1d", H.rd1d(40), []), ("bruss", H.brusselator2d(6, 6), []),
              ("llg", H.llg_dipole(8), []), ("kuramoto", (rhs_f, params_f), forcings_f)]
    units = []
    for label, (rhs, params), forcings in models:
        for style, scalar, lvl in SCALARS:
            for kw in (dict(sparse=False), dict(sparse=True),
                       dict(sparse=False, emit_contractions=True, emit_jvp=True)):
                for strat in ("entries", "colour"):
                    os.environ["CPPDE_JAC"] = strat
                    try:
                        r = CO.generate_ode_cpp(rhs, params, num_type=scalar,
                                                ad_level=lvl, arena=lvl > 0,
                                                forcings_list=forcings, **kw)
                    finally:
                        os.environ.pop("CPPDE_JAC", None)
                    units.append("namespace u%d { // %s %s %s %s\n%s\n%s\n%s\n%s\n}\n"
                                 % (len(units), label, scalar, sorted(kw), strat,
                                    r["data_code"], r["ode_code"], r["jac_code"],
                                    r["adj_code"]))
    d = tempfile.mkdtemp(prefix="cppde_models_")
    src = os.path.join(d, "models.cpp")
    with open(src, "w", encoding="utf-8") as f:
        f.write(CXX_HEAD + "\n".join(units))
    r = subprocess.run([cxx, "-std=gnu++17", "-fsyntax-only", "-I", inc,
                        "-I", rinc, src], capture_output=True, text=True)
    errs = [l for l in r.stderr.splitlines() if "error" in l]
    t.check(r.returncode == 0, "g++ exit %d, %d error lines: %s"
            % (r.returncode, len(errs), " | ".join(errs[:3])))
    print("  %d units in %s" % (len(units), src))
    return t.report()


def check_cvode_syntax(corpus):
    """g++ -fsyntax-only over generated CVODE sources: dense and KLU, with
    sensitivities, adjoint, events, a linear map and loops."""
    import tempfile
    import codegen_cvode as CV
    t = H.Tally("cvode: C++ syntax of generated sources")
    cxx, rinc = tool_paths()
    sund = os.environ.get("SUNDIALS_INCLUDE") or r"C:\rtools45\ucrt64\include"
    if not os.path.exists(cxx) or not os.path.exists(os.path.join(sund, "cvodes")):
        t.skip("no compiler or SUNDIALS headers")
        return t.report()
    d = tempfile.mkdtemp(prefix="cppde_cvode_")
    ev = {"var": ["A", "B"], "value": ["A + dose*B", "B*fac"], "time": ["t_dose*2", None],
          "root": [None, "A^2 - thr*B - 0.1*time"], "method": ["add", "replace"],
          "terminal": [None, False], "direction": [None, 0]}
    rhs_e = {"A": "-k1*A + k2*B*u", "B": "k1*A - k2*B^2"}
    llg, llg_p = H.llg_dipole(8)
    rd, rd_p = H.rd1d(40)
    runs = [
        ("events", rhs_e, ["k1", "k2", "dose", "t_dose", "fac", "thr"], ["u"],
         dict(deriv=True, events=ev, sparse=False)),
        ("rootfunc", rhs_e, ["k1", "k2"], ["u"],
         dict(deriv=True, rootfunc=["A - 0.5"], sparse=True)),
        ("adjoint", rd, rd_p, [], dict(deriv=True, reverse=True, sparse=True)),
        ("llg_dense", llg, llg_p, [], dict(deriv=True, reverse=True, sparse=False)),
        ("llg_klu", llg, llg_p, [], dict(deriv=True, sparse=True)),
    ]
    for name, rhs, params, forcings, kw in runs:
        for strat in ("entries", "colour"):
            os.environ["CPPDE_JAC"] = strat
            try:
                r = CV.generate_cvode_cpp(rhs, params, "cv_%s_%s" % (name, strat), d,
                                          forcings_list=forcings, **kw)
            finally:
                os.environ.pop("CPPDE_JAC", None)
            args = [cxx, "-std=gnu++17", "-fsyntax-only", "-I", include_dir(),
                    "-I", rinc, "-I", sund, "-I", os.path.join(sund, "suitesparse")]
            args += list(r["compile_defs"]) + [r["srcfile"]]
            out = subprocess.run(args, capture_output=True, text=True)
            errs = [l for l in out.stderr.splitlines() if "error" in l]
            t.check(out.returncode == 0, "%s %s: g++ exit %d: %s"
                    % (name, strat, out.returncode,
                       " | ".join(errs[:3]) or out.stderr.strip()[:300]))
    print("  sources in %s" % d)
    return t.report()


# ---------------------------------------------------------------------------
# Phase 2: adjoint_terms
# ---------------------------------------------------------------------------

ARGS = ["x", "v", "lam", "t", "params", "F", "sc", "out"]
FD_TOL = 1e-5


def ode_models(corpus, limit=None):
    out = []
    for i, e in enumerate(corpus.get("ode", [])):
        label = e.get("name") or "ode[%d]" % i
        out.append((label, e["rhs_dict"], e.get("params_list"), e.get("forcings_list")))
    return out[:limit] if limit else out


def compile_functions(m, emit_jvp=True):
    """name -> (kind, python function) of adjoint_terms."""
    import cppde_model as CM
    fns = {}
    for name, kind, stores in m.contraction_stores(emit_jvp):
        stmts = CM.stores_prelude(m, stores, "double") + [CM.block(m.g, stores, "_t")]
        fns[name] = (kind, py_body(m, stmts, ARGS))
    return fns


def call(fn, kind, m, pt, x=None, lam=None, v=None, p=None):
    x = pt.x if x is None else x
    lam = pt.lam if lam is None else lam
    v = pt.v if v is None else v
    p = pt.p if p is None else p
    if kind == "dot":
        return [fn(x, v, lam, pt.t, p, pt.F, pt.sc, None)]
    if kind == "x":
        out = [0.0] * m.n
    else:
        out = [0.0] * (m.n + len(m.params))
    fn(x, v, lam, pt.t, p, pt.F, pt.sc, out)
    return out


def near(a, b, floor):
    return math.isfinite(a) and math.isfinite(b) and abs(a - b) <= floor


def cancel_floor(values):
    """Rounding left by cancelling terms: 1e-12 of the largest finite output
    of the same function."""
    return 1e-12 * max([abs(v) for v in values if math.isfinite(v)] + [0.0])


def reference(m):
    """name -> (kind, [sympy expr]) from the legacy derivatives."""
    local = H.sym_table(m.names())
    exprs = [H.legacy_parse(e, local) for e in m_exprs(m)]
    xs = [local[s] for s in m.states]
    t = local["time"]
    lam = [sp.Symbol("_lam%d" % i, real=True) for i in range(m.n)]
    v = [sp.Symbol("_v%d" % i, real=True) for i in range(m.n)]
    fs = [local[f] for f in m.forcings]
    rates = [sp.Symbol("_fd%d" % j, real=True) for j in range(len(m.forcings))]
    ps = [(k, local[p]) for k, p in enumerate(m.params) if p not in m.init_names]

    def pull(fs_):
        jx = [sum((sp.diff(f, xj) * lam[i] for i, f in enumerate(fs_)), sp.Integer(0))
              for xj in xs]
        jp = {k: sum((sp.diff(f, pk) * lam[i] for i, f in enumerate(fs_)),
                     sp.Integer(0)) for k, pk in ps}
        return [H.replace_dirac(e) for e in jx], {k: H.replace_dirac(e) for k, e in jp.items()}

    J = H.ref_jacobian(exprs, xs)
    ref = {}
    ref["jac_t_vec"] = ("x", [sum((J[i][j] * lam[i] for i in range(m.n)), sp.Integer(0))
                              for j in range(m.n)])
    dfdp = H.ref_dfdp(exprs, [pk for _, pk in ps])
    kmap = [k for k, _ in ps]
    acc = {}
    for i, kk, e in dfdp:
        acc[kmap[kk]] = acc.get(kmap[kk], 0) + e * lam[i]
    ref["dfdp_t_vec_axpy"] = ("p", acc)
    jv = [sum((J[i][k] * v[k] for k in range(m.n)), sp.Integer(0)) for i in range(m.n)]
    jx, jp = pull(jv)
    ref["jvp_x_t_vec"] = ("x", jx)
    ref["jvp_p_t_vec_axpy"] = ("p", jp)
    dd = H.ref_dfdt(exprs, t, fs, rates)
    px, pp = pull(dd)
    ref["dfdt_x_t_vec"] = ("x", px)
    ref["dfdt_p_t_vec_axpy"] = ("p", pp)
    ref["dfdt_dot"] = ("dot", [sum((d * lam[i] for i, d in enumerate(dd)), sp.Integer(0))])
    names = m.names() + [str(s) for s in lam + v + rates]
    return ref, names


def m_exprs(m):
    return m.exprs


def check_contractions(corpus, models_limit=None):
    """The seven contractions against the legacy SymPy ones."""
    t = H.Tally("adjoint: contractions vs SymPy")
    rng = random.Random(7)
    for label, rhs, params, forcings in ode_models(corpus, models_limit):
        try:
            m = H.OdeModel(rhs, params, forcings)
            m.exprs = [str(rhs[s]) for s in m.states]
            fns = compile_functions(m)
            ref, names = reference(m)
            evs = {}
            for name, (kind, exprs) in ref.items():
                if kind == "p":
                    keys = sorted(exprs)
                    evs[name] = (kind, keys, H.Evaluator([exprs[k] for k in keys], names))
                else:
                    evs[name] = (kind, None, H.Evaluator(exprs, names))
        except Exception as e:
            t.skip("%s: %s" % (label, str(e)[:80]))
            continue
        pts, _ = H.sample_points(m, m.rhs, rng, base=6, extra=24)
        for pt in pts:
            env = pt.env(m)
            ill = H.ill_conditioned(m, pt)
            for name, (kind, fn) in fns.items():
                got = call(fn, kind, m, pt)
                rk, keys, ev = evs[name]
                want = ev(env)
                if kind == "p":
                    want_full = [0.0] * len(got)
                    for k, w in zip(keys, want):
                        want_full[m.n + k] = pt.sc * w
                    want = want_full
                floor = cancel_floor(want)
                for i, (a, b) in enumerate(zip(got, want)):
                    t.compare(a, a if near(a, b, floor) else b, 1e-12,
                              "%s %s[%d]: new %r legacy %r" % (label, name, i, a, b),
                              singular=ill)
    return t.report()


def check_fd(corpus, models_limit=None):
    """jac_t_vec and dfdp against finite differences; the dual instantiation
    of jac_t_vec against finite differences of the double one."""
    t = H.Tally("adjoint: finite differences, dual")
    rng = random.Random(11)
    for label, rhs, params, forcings in ode_models(corpus, models_limit):
        try:
            m = H.OdeModel(rhs, params, forcings)
            fns = compile_functions(m)
        except Exception as e:
            t.skip("%s: %s" % (label, str(e)[:80]))
            continue
        f = rhs_function(m)
        for _ in range(3):
            pt = H.Point(m, rng, 0)
            base = call(f, "x", m, pt)
            if not all(math.isfinite(v) for v in base):
                continue
            if H.ill_conditioned(m, pt):
                continue
            lam_f = lambda **kw: sum(a * l for a, l in zip(call(f, "x", m, pt, **kw), pt.lam))
            # J' lam by differences of lam . f
            kind, jt = fns["jac_t_vec"]
            got = call(jt, kind, m, pt)
            for j in range(m.n):
                fd = fd5(lambda v: lam_f(x=_with(pt.x, j, v)), pt.x[j])
                if fd is None:
                    continue
                t.check(abs(fd - got[j]) <= FD_TOL * max(1.0, abs(fd)),
                        "%s jac_t_vec[%d]: %r vs fd %r" % (label, j, got[j], fd))
            # (df/dp)' lam
            kind, dp = fns["dfdp_t_vec_axpy"]
            gotp = call(dp, kind, m, pt)
            for k in range(len(m.params)):
                if m.params[k] in m.init_names:
                    continue
                s = m.n + k
                fd = fd5(lambda v: lam_f(p=_with(pt.p, s, v)), pt.p[s])
                if fd is None:
                    continue
                fd *= pt.sc
                t.check(abs(fd - gotp[s]) <= FD_TOL * max(1.0, abs(fd)),
                        "%s dfdp[%d]: %r vs fd %r" % (label, k, gotp[s], fd))
            # jac_t_vec over a dual in x: d/dx_j (J' lam) against differences
            xd = [H.pyrt.Dual(xv, [1.0 if k == j else 0.0 for k in range(m.n)])
                  for j, xv in enumerate(pt.x)]
            outd = call(jt, "x", m, pt, x=xd)
            for j in range(min(m.n, 6)):
                for i in range(m.n):
                    fd = fd5(lambda v: call(jt, "x", m, pt, x=_with(pt.x, j, v))[i], pt.x[j])
                    if fd is None:
                        continue
                    ad = H.pyrt.tangent(outd[i], j)
                    t.check(abs(fd - ad) <= FD_TOL * max(1.0, abs(fd)),
                            "%s dual jac_t_vec[%d]/dx%d: %r vs fd %r" % (label, i, j, ad, fd))
    return t.report()


def _with(v, i, x):
    out = list(v)
    out[i] = x
    return out


def fd5(fn, x0, rel=1e-4):
    """Five-point central difference of a scalar function at x0, or None
    where two step sizes disagree (a kink or jump nearby)."""
    def est(h):
        return (fn(x0 - 2 * h) - 8 * fn(x0 - h) + 8 * fn(x0 + h) - fn(x0 + 2 * h)) / (12 * h)
    h = rel * max(1.0, abs(x0))
    a, b = est(h), est(h / 2)
    if not (math.isfinite(a) and math.isfinite(b)) or abs(a - b) > 1e-6 * max(1.0, abs(b)):
        return None
    return b


# ---------------------------------------------------------------------------
# Phase 3: ode_system and jacobian
# ---------------------------------------------------------------------------

JAC_ARGS = ["x", "J", "W", "t", "dfdt", "params", "F"]


def rhs_function(m):
    """f(x, v, lam, t, params, F, sc, out) writing the right-hand side."""
    import cppde_model as CM
    stores = [(("vec", "out", i), r, "=") for i, r in enumerate(m.rhs)]
    return py_body(m, CM.stores_prelude(m, stores, "double")
                   + [CM.block(m.g, stores, "_t")], ARGS)


def py_body(m, stmts, args):
    """Python function def f(*args) from statements."""
    body = H.cppde_model.py_statements_raw(m, stmts)
    src = "def f(%s):\n%s\n    return None\n" % (", ".join(args), "\n".join(body))
    scope = {"R": H.pyrt, "linmap_": m.linmap}
    exec(compile(src, "<generated>", "exec"), scope)
    f = scope["f"]
    f.source = src
    return f


def jac_reference(m):
    """(f, J, dfdt) of the legacy generator, as SymPy."""
    local = H.sym_table(m.names())
    exprs = [H.legacy_parse(str(e), local) for e in m.exprs]
    xs = [local[s] for s in m.states]
    fs = [local[f] for f in m.forcings]
    rates = [sp.Symbol("_fd%d" % j, real=True) for j in range(len(m.forcings))]
    J = H.ref_jacobian(exprs, xs)
    dd = H.ref_dfdt(exprs, local["time"], fs, rates)
    names = m.names() + [str(r) for r in rates]
    return exprs, J, dd, names


def run_jac(fn, m, pt, sparse, n, W=None, x=None, p=None, t=None):
    """-J as dense rows and dfdt of one functor call."""
    x = pt.x if x is None else x
    p = pt.p if p is None else p
    t = pt.t if t is None else t
    dfdt = [0.0] * n
    if sparse:
        W = W if W is not None else H.pyrt.CscMatrix()
        fn(x, None, W, t, dfdt, p, pt.F)
        return W.dense(), dfdt, W
    J = H.pyrt.DenseMatrix(n)
    fn(x, J, None, t, dfdt, p, pt.F)
    return J.data, dfdt, None


def check_jacobian(corpus, models_limit=None):
    """ode_system and jacobian (entries, colour; dense, sparse) against the
    legacy SymPy derivatives; colour against entries."""
    import cppde_model as CM
    t = H.Tally("jacobian: functors vs SymPy")
    tp = H.Tally("jacobian: pattern, decide_sparse, KLU")
    rng = random.Random(21)
    for label, rhs, params, forcings in ode_models(corpus, models_limit):
        try:
            m = H.OdeModel(rhs, params, forcings)
            m.exprs = [str(rhs[s]) for s in m.states]
            exprs, J, dd, names = jac_reference(m)
            ev_f = H.Evaluator(exprs, names)
            flatJ = [(i, j) for i in range(m.n) for j in range(m.n) if J[i][j] != 0]
            ev_J = H.Evaluator([J[i][j] for i, j in flatJ], names)
            ev_d = H.Evaluator(dd, names)
            f = py_body(m, CM.ode_statements(m), ["x", "dxdt", "t", "params", "F"])
            fns = {(s, sp_): py_body(m, CM.jacobian_statements(m, sp_, s), JAC_ARGS)
                   for s in ("entries", "colour") for sp_ in (False, True)}
        except Exception as e:
            t.skip("%s: %s" % (label, str(e)[:100]))
            continue
        # pattern and the decisions taken from it
        old = sorted(flatJ)
        new = sorted((i, j) for i, row in enumerate(m.pattern()) for j in row)
        tp.check(old == new, "%s pattern differs: %d legacy, %d new"
                 % (label, len(old), len(new)))
        if old == new:
            r = [p_[0] for p_ in old]
            c = [p_[1] for p_ in old]
            tp.check(H.legacy.ODE.analyze_klu_settings(m.n, r, c)
                     == cppde_model_klu(m.n, r, c), "%s klu settings" % label)
        pts, _ = H.sample_points(m, m.rhs, rng, base=5, extra=20)
        for pt in pts:
            env = pt.env(m)
            ill = H.ill_conditioned(m, pt)
            want_f = ev_f(env)
            got_f = [0.0] * m.n
            f(pt.x, got_f, pt.t, pt.p, pt.F)
            for i, (a, b) in enumerate(zip(got_f, want_f)):
                t.compare(a, b, H.TOL, "%s f[%d]: %r vs %r" % (label, i, a, b))
            wantJ = dict(zip(flatJ, ev_J(env)))
            want_d = ev_d(env)
            # a coloured row is only as finite as its worst entry (0 * inf)
            bad_rows = {i for (i, j), b in wantJ.items() if not math.isfinite(b)}
            for (strat, sparse), fn in fns.items():
                negJ, dfdt, _ = run_jac(fn, m, pt, sparse, m.n)
                for i in range(m.n):
                    for j in range(m.n):
                        b = wantJ.get((i, j), 0.0)
                        if strat == "colour" and i in bad_rows:
                            b = math.nan
                        t.compare(-negJ[i][j], b, 1e-12, "%s %s/%s J[%d,%d]: %r vs %r"
                                  % (label, strat, "csc" if sparse else "dense",
                                     i, j, -negJ[i][j], b), singular=ill)
                for i, (a, b) in enumerate(zip(dfdt, want_d)):
                    t.compare(a, b, 1e-12, "%s %s dfdt[%d]: %r vs %r"
                              % (label, strat, i, a, b), singular=ill)
    ok = t.report()
    return tp.report() and ok


def cppde_model_klu(n, rows, cols):
    import codegen_cppODE
    return codegen_cppODE.analyze_klu_settings(n, rows, cols)


def check_jacobian_reuse(corpus, models_limit=None):
    """A reused CSC matrix: clock-dependent entries follow t, dense entries
    are reset between calls."""
    t = H.Tally("jacobian: reuse across calls")
    rng = random.Random(5)
    specs = [("time-only entry", {"x": "-k*exp(-time)*x + y", "y": "x - y*time"},
              ["k"], [])]
    for label, rhs, params, forcings in specs + ode_models(corpus, models_limit)[:20]:
        try:
            m = H.OdeModel(rhs, params, forcings)
        except Exception as e:
            t.skip("%s: %s" % (label, e))
            continue
        import cppde_model as CM
        for strat in ("entries", "colour"):
            for sparse in (True, False):
                fn = py_body(m, CM.jacobian_statements(m, sparse, strat), JAC_ARGS)
                pt = H.Point(m, rng, 0)
                J1, _, W = run_jac(fn, m, pt, sparse, m.n)
                J2, _, _ = run_jac(fn, m, pt, sparse, m.n, W=W, t=pt.t + 1.7,
                                   x=[v * 1.3 for v in pt.x])
                J3, _, _ = run_jac(fn, m, pt, sparse, m.n, t=pt.t + 1.7,
                                   x=[v * 1.3 for v in pt.x])
                for i in range(m.n):
                    for j in range(m.n):
                        t.check(H.close(J2[i][j], J3[i][j], 1e-14),
                                "%s %s %s J[%d,%d] reused %r fresh %r"
                                % (label, strat, sparse, i, j, J2[i][j], J3[i][j]))
    return t.report()


def check_jacobian_dual(corpus, models_limit=None):
    """Tangents of the jacobian functor under Dual and nested Dual parameters
    against SymPy d/dp of the legacy Jacobian."""
    import cppde_model as CM
    t = H.Tally("jacobian: dual and nested dual tangents")
    rng = random.Random(9)
    for label, rhs, params, forcings in ode_models(corpus, models_limit)[:25]:
        try:
            m = H.OdeModel(rhs, params, forcings)
            m.exprs = [str(rhs[s]) for s in m.states]
            exprs, J, dd, names = jac_reference(m)
            local = H.sym_table(m.names())
            pidx = [k for k, p_ in enumerate(m.params) if p_ not in m.init_names][:3]
            psyms = [local[m.params[k]] for k in pidx]
            flatJ = [(i, j) for i in range(m.n) for j in range(m.n) if J[i][j] != 0]
            dJ = [[H.replace_dirac(sp.diff(J[i][j], ps)) for ps in psyms] for i, j in flatJ]
            ev = H.Evaluator([e for row in dJ for e in row], names)
            d2 = [H.replace_dirac(sp.diff(J[i][j], psyms[0], psyms[-1])) for i, j in flatJ] if psyms else []
            ev2 = H.Evaluator(d2, names) if psyms else None
        except Exception as e:
            t.skip("%s: %s" % (label, str(e)[:100]))
            continue
        if not pidx:
            continue
        for strat in ("entries", "colour"):
            fn = py_body(m, CM.jacobian_statements(m, False, strat), JAC_ARGS)
            for _ in range(2):
                pt = H.Point(m, rng, 0)
                np_ = len(pidx)
                pd = list(pt.p)
                for q, k in enumerate(pidx):
                    pd[m.n + k] = H.pyrt.Dual(pt.p[m.n + k], [1.0 if r == q else 0.0 for r in range(np_)])
                pd = [v if isinstance(v, H.pyrt.Dual) else H.pyrt.Dual(v, [0.0] * np_) for v in pd]
                negJ, _, _ = run_jac(fn, m, pt, False, m.n, p=pd)
                want = ev(pt.env(m))
                for e_i, (i, j) in enumerate(flatJ):
                    for q in range(np_):
                        a = -H.pyrt.tangent(negJ[i][j], q)
                        b = want[e_i * np_ + q]
                        t.compare(a, b, 1e-11, "%s %s dJ[%d,%d]/dp%d: %r vs %r"
                                  % (label, strat, i, j, q, a, b))
                # nested: d2/dp0 dp_last
                k0, k1 = pidx[0], pidx[-1]
                pn = []
                for idx, v in enumerate(pt.p):
                    inner = H.pyrt.Dual(v, [1.0 if idx == m.n + k1 else 0.0])
                    outer_t = [H.pyrt.Dual(1.0 if idx == m.n + k0 else 0.0, [0.0])]
                    pn.append(H.pyrt.Dual(inner, outer_t))
                negJ, _, _ = run_jac(fn, m, pt, False, m.n, p=pn)
                want2 = ev2(pt.env(m))
                for e_i, (i, j) in enumerate(flatJ):
                    a = -H.pyrt.hessian(negJ[i][j], 0, 0)
                    t.compare(a, want2[e_i], 1e-10, "%s %s d2J[%d,%d]: %r vs %r"
                              % (label, strat, i, j, a, want2[e_i]))
    return t.report()


# ---------------------------------------------------------------------------
# Phase 4: events and roots
# ---------------------------------------------------------------------------

EV_ARGS = ["x", "t", "params", "full_params", "F", "sc", "out"]


def event_specs(corpus):
    """(label, rhs, params, forcings, rows, has_rhs) per event set."""
    import cppde_model as CM
    out = []
    for i, e in enumerate(H.EVENT_MODELS + corpus.get("event", [])):
        if e.get("emit_adjoint"):
            continue
        rows = CM.event_rows(e.get("events_df"))
        if not rows:
            continue
        states = H.as_list(e.get("states_list"))
        rhs = e.get("rhs_dict")
        has_rhs = rhs is not None
        rhs = {s: (rhs[s] if has_rhs else "0") for s in states}
        out.append((e.get("name") or "event[%d]" % i, rhs, e.get("params_list"),
                    e.get("forcings_list"), rows, has_rhs))
    return out


def ev_fn(m, stores):
    return H.py_function(m.g, stores, m.slot("py"), EV_ARGS, prefix="_e")


def ev_call(fn, m, pt, n_out, x=None):
    out = [0.0] * n_out
    r = fn(pt.x if x is None else x, pt.t, pt.p, pt.p, pt.F, pt.sc, out)
    return out if r is None else [r]


def check_events(corpus):
    """Event lambdas and event_adjoint_terms cases against the legacy
    event math, dg/dx against finite differences."""
    import cppde_model as CM
    t = H.Tally("events: lambdas and adjoint cases vs SymPy")
    tf = H.Tally("events: dg/dx vs finite differences")
    rng = random.Random(31)
    for label, rhs, params, forcings, rows, has_rhs in event_specs(corpus):
        try:
            m = H.OdeModel(rhs, params, forcings)
            code = CM.EventCode(m, CM.Scalar("double"), has_rhs=has_rhs)
            local = H.sym_table(m.names())
            f_ref = [H.legacy_parse(str(rhs[s]), local) for s in m.states]
            xs = [local[s] for s in m.states]
            tt = local["time"]
            fs = [local[f] for f in m.forcings]
            rates = [sp.Symbol("_fd%d" % j, real=True) for j in range(len(m.forcings))]
            names = m.names() + [str(r) for r in rates] + ["_sc"]
            sc = sp.Symbol("_sc", real=True)
            pk = [(k, local[p]) for k, p in enumerate(m.params) if p not in m.init_names]
            inits = [(i, local[s + "_0"]) for i, s in enumerate(m.states)
                     if s + "_0" in local]
        except Exception as e:
            t.skip("%s: %s" % (label, str(e)[:100]))
            continue

        def dt_forcing(e):
            d = H.replace_dirac(sp.diff(e, tt))
            for fsym, r in zip(fs, rates):
                d = d + H.replace_dirac(sp.diff(e, fsym)) * r
            return d

        def p_slots(e):
            out = {}
            for k, s_ in pk:
                d = H.replace_dirac(sp.diff(e, s_))
                if d != 0:
                    out[m.n + k] = sc * d
            for i, s_ in inits:
                d = H.replace_dirac(sp.diff(e, s_))
                if d != 0:
                    out[i] = out.get(i, 0) + sc * d
            return out

        checks = []   # (what, stores or node, ref exprs by slot, n_out)
        for r in rows:
            i = r["index"]
            if not CM.valid_value(r["value"]):
                continue
            h = m.parse(str(r["value"]))
            h_ref = H.legacy_parse(str(r["value"]), local)
            checks.append(("h", [(("ret",), h, "=")], {0: h_ref}, 0))
            checks.append(("dh_dx", code.case_x_stores(h),
                           {j: H.replace_dirac(sp.diff(h_ref, x)) for j, x in enumerate(xs)}, m.n))
            checks.append(("dh_dp", code.case_p_stores(h), p_slots(h_ref),
                           m.n + len(m.params)))
            checks.append(("dh_dt", code.case_t_stores(h), {0: sc * dt_forcing(h_ref)}, 1))
            if CM.valid_value(r["time"]):
                tn = m.parse(str(r["time"]))
                t_ref = H.legacy_parse(str(r["time"]), local)
                checks.append(("time", [(("ret",), tn, "=")], {0: t_ref}, 0))
                checks.append(("dtime_dp", code.case_p_stores(tn), p_slots(t_ref),
                               m.n + len(m.params)))
            elif CM.valid_value(r["root"]):
                gn = m.parse(str(r["root"]))
                g_ref = H.legacy_parse(str(r["root"]), local)
                checks.append(("g", [(("ret",), gn, "=")], {0: g_ref}, 0))
                checks.append(("dg_dx", code.dg_dx_stores(gn),
                               {j: H.replace_dirac(sp.diff(g_ref, x)) for j, x in enumerate(xs)}, m.n))
                checks.append(("dg_dt", [(("ret",), code.dg_dt_node(gn), "=")],
                               {0: dt_forcing(g_ref)}, 0))
                checks.append(("dg_dp", code.case_p_stores(gn), p_slots(g_ref),
                               m.n + len(m.params)))
                gdot = dt_forcing(g_ref) + sum(
                    (sp.diff(g_ref, x) * fi for x, fi in zip(xs, f_ref)), sp.Integer(0))
                gd = code.gdot_node(gn)
                checks.append(("gdot_dx", code.case_x_stores(gd),
                               {j: H.replace_dirac(sp.diff(gdot, x)) for j, x in enumerate(xs)}, m.n))
                checks.append(("gdot_dp", code.case_p_stores(gd), p_slots(gdot),
                               m.n + len(m.params)))
                if has_rhs and not m.forcings:
                    g1 = sp.diff(g_ref, tt) + sum(
                        (sp.diff(g_ref, x) * fi for x, fi in zip(xs, f_ref)), sp.Integer(0))
                    gtt = sp.diff(g1, tt) + sum(
                        (sp.diff(g1, x) * fi for x, fi in zip(xs, f_ref)), sp.Integer(0))
                    checks.append(("G_tt", [(("ret",), code.gtt_node(gn), "=")],
                                   {0: H.replace_dirac(gtt)}, 0))
                    # finite differences of g in x
                    tf_item = (gn, code.dg_dx_stores(gn))
                else:
                    tf_item = (gn, code.dg_dx_stores(gn))
                fg = ev_fn(m, [(("ret",), tf_item[0], "=")])
                fdx = ev_fn(m, tf_item[1])
                for _ in range(2):
                    pt = H.Point(m, rng, 0)
                    got = ev_call(fdx, m, pt, m.n)
                    for j in range(m.n):
                        hj = 1e-6 * max(1.0, abs(pt.x[j]))
                        xp = list(pt.x); xm = list(pt.x)
                        xp[j] += hj; xm[j] -= hj
                        fd = (ev_call(fg, m, pt, 0, xp)[0] - ev_call(fg, m, pt, 0, xm)[0]) / (2 * hj)
                        if math.isfinite(fd):
                            tf.check(abs(fd - got[j]) <= 1e-6 * max(1.0, abs(fd)),
                                     "%s event %d dg/dx[%d]: %r vs fd %r"
                                     % (label, i, j, got[j], fd))
        try:
            compiled = []
            for what, stores, ref, n_out in checks:
                keys = sorted(ref)
                ev = H.Evaluator([ref[k] for k in keys], names)
                compiled.append((what, ev_fn(m, stores), keys, ev, n_out))
        except Exception as e:
            t.skip("%s: %s" % (label, str(e)[:100]))
            continue
        pts, _ = H.sample_points(m, m.rhs, rng, base=4, extra=12)
        for pt in pts:
            env = pt.env(m)
            env["_sc"] = pt.sc
            for what, fn, keys, ev, n_out in compiled:
                got = ev_call(fn, m, pt, n_out)
                want = dict(zip(keys, ev(env)))
                for k in range(max(len(got), 1)):
                    a = got[k]
                    b = want.get(k, 0.0)
                    t.compare(a, b, 1e-12, "%s %s[%d]: %r vs %r" % (label, what, k, a, b))
    ok = t.report()
    return tf.report() and ok


# ---------------------------------------------------------------------------
# Phase 5: cppFUN
# ---------------------------------------------------------------------------

def fun_specs(corpus):
    out = list(H.FUN_MODELS)
    for i, e in enumerate(corpus.get("fun", [])):
        exprs = e.get("exprs")
        if isinstance(exprs, list):
            exprs = {"f%d" % (k + 1): v for k, v in enumerate(exprs)}
        out.append(("fun[%d]" % i, exprs, H.as_list(e.get("variables")),
                    H.as_list(e.get("parameters"))))
    return out


def fun_py(model, stmts, args, pre=""):
    import cppde_emit as em
    pr = em.Printer(model.g, model.slot, style="py")
    body = em.render_py(stmts, pr)
    src = "def f(%s):\n%s%s\n    return None\n" % (", ".join(args), pre, "\n".join(body))
    scope = {"R": H.pyrt}
    exec(compile(src, "<generated>", "exec"), scope)
    return scope["f"]


def check_fun(corpus):
    """cppFUN eval and vjp against SymPy; the dual vjp against finite
    differences of the double one."""
    import codegen_cppFUN as CF
    t = H.Tally("cppFUN: eval and vjp vs SymPy")
    tf = H.Tally("cppFUN: dual vjp vs finite differences")
    rng = random.Random(41)
    for label, exprs, variables, parameters in fun_specs(corpus):
        try:
            m = CF.FunModel(exprs, variables, parameters)
            names = variables + parameters
            ctx = H.legacy.FUN.CodeGenContext(variables, parameters)
            ref = [ctx._parse_expr(str(exprs[n])) for n in m.out_names]
            syms = [ctx.all_symbols[s] for s in names]
            ws = [sp.Symbol("_w%d" % i, real=True) for i in range(len(ref))]
            adj = [H.replace_dirac(sum((sp.diff(e, s) * w for e, w in zip(ref, ws)),
                                       sp.Integer(0))) for s in syms]
            ev_y = H.Evaluator(ref, names)
            ev_a = H.Evaluator(adj, names + [str(w) for w in ws])
            f_eval = fun_py(m, CF.eval_statements(m), ["x_obs", "p", "y_local"])
            outer, inner = CF.vjp_statements(m, cast="")
            nv, npar = len(variables), len(parameters)
            pre = "    obs = 0\n    n_obs = 1\n    s = 0\n    n_vars = %d\n    n_params = %d\n" % (nv, npar)
            f_vjp = fun_py(m, outer + inner, ["x_obs", "p", "_w", "wx", "wp"], pre)
        except Exception as e:
            t.skip("%s: %s" % (label, str(e)[:100]))
            continue
        for k in range(12):
            kind = k % 4
            vals = [H.draw(rng, kind) for _ in names]
            wv = [H.draw(rng, 1) for _ in ref]
            env = dict(zip(names, vals))
            x, p = vals[:nv], vals[nv:]
            y = [0.0] * len(ref)
            f_eval(x, p, y)
            for i, (a, b) in enumerate(zip(y, ev_y(env))):
                t.compare(a, b, 1e-12, "%s y[%d]: %r vs %r" % (label, i, a, b))
            wx = [0.0] * nv
            wp = [0.0] * npar
            f_vjp(x, p, wv, wx, wp)
            env.update({"_w%d" % i: w for i, w in enumerate(wv)})
            for i, (a, b) in enumerate(zip(wx + wp, ev_a(env))):
                t.compare(a, b, 1e-12, "%s adj[%d]: %r vs %r" % (label, i, a, b))
            if kind != 0 or not names:
                continue
            # vjp2: tangents of the dual vjp against differences of the vjp
            nd = len(names)
            xd = [H.pyrt.Dual(v, [1.0 if q == j else 0.0 for q in range(nd)])
                  for j, v in enumerate(vals)]
            wxd = [0.0] * nv
            wpd = [0.0] * npar
            f_vjp(xd[:nv], xd[nv:], wv, wxd, wpd)
            out_d = wxd + wpd
            for j in range(nd):
                h = 1e-6 * max(1.0, abs(vals[j]))
                vp = list(vals); vm = list(vals)
                vp[j] += h; vm[j] -= h
                ap_x, ap_p = [0.0] * nv, [0.0] * npar
                am_x, am_p = [0.0] * nv, [0.0] * npar
                f_vjp(vp[:nv], vp[nv:], wv, ap_x, ap_p)
                f_vjp(vm[:nv], vm[nv:], wv, am_x, am_p)
                for i, (a1, a2) in enumerate(zip(ap_x + ap_p, am_x + am_p)):
                    fd = (a1 - a2) / (2 * h)
                    if not math.isfinite(fd):
                        continue
                    ad = H.pyrt.tangent(out_d[i], j)
                    tf.check(abs(fd - ad) <= 1e-5 * max(1.0, abs(fd)),
                             "%s d adj[%d]/d%s: %r vs fd %r" % (label, i, names[j], ad, fd))
    ok = t.report()
    return tf.report() and ok



# ---------------------------------------------------------------------------
# Phase 6: linear maps
# ---------------------------------------------------------------------------

LIN_TOL = 1e-13


def _lin_functions(m):
    import cppde_model as CM
    fs = {"f": ("x", rhs_function(m))}
    fs.update(compile_functions(m))
    for strat in ("entries", "colour"):
        for sparse in (False, True):
            fs[(strat, sparse)] = ("jac", py_body(
                m, CM.jacobian_statements(m, sparse, strat), JAC_ARGS))
    return fs


def check_linear(corpus):
    """Every ODE function with the linear map against the same function
    without it; the Jacobian also across reused CSC calls and the adjoint
    under a Dual state."""
    t = H.Tally("linear map: functions with and without E")
    rng = random.Random(31)
    rows = 0
    for label, rhs, params, forcings in list(H.LINEAR_MODELS) + ode_models(corpus):
        try:
            m0 = H.OdeModel(rhs, params, forcings, linear=0)
            m1 = H.OdeModel(rhs, params, forcings, linear=3)
            if not m1.nlin:
                continue
            f0 = _lin_functions(m0)
            f1 = _lin_functions(m1)
        except Exception as e:
            t.skip("%s: %s" % (label, str(e)[:100]))
            continue
        rows += 1
        n = m0.n
        t.check(m0.pattern() == m1.pattern(), "%s: pattern differs" % label)
        pts, _ = H.sample_points(m0, m0.rhs, rng, base=4, extra=12)
        for pt in pts:
            ill = H.ill_conditioned(m0, pt)
            for key, (kind, fn0) in f0.items():
                fn1 = f1[key][1]
                if kind == "jac":
                    strat, sparse = key
                    a, da, W = run_jac(fn1, m1, pt, sparse, n)
                    b, db, _ = run_jac(fn0, m0, pt, sparse, n)
                    got = [v for row in a for v in row] + da
                    want = [v for row in b for v in row] + db
                    if sparse:
                        xs = [v * 1.3 for v in pt.x]
                        a2, _, _ = run_jac(fn1, m1, pt, sparse, n, W=W, x=xs, t=pt.t + 0.7)
                        b2, _, _ = run_jac(fn0, m0, pt, sparse, n, x=xs, t=pt.t + 0.7)
                        got += [v for row in a2 for v in row]
                        want += [v for row in b2 for v in row]
                else:
                    got = call(fn1, kind, m1, pt)
                    want = call(fn0, kind, m0, pt)
                floor = cancel_floor(want)
                for i, (x, y) in enumerate(zip(got, want)):
                    t.compare(x, x if near(x, y, floor) else y, LIN_TOL,
                              "%s %s[%d]: E %r plain %r" % (label, key, i, x, y),
                              singular=ill)
            # jac_t_vec under a Dual state: tangents
            xd = [H.pyrt.Dual(v, [1.0 if k == j else 0.0 for k in range(n)])
                  for j, v in enumerate(pt.x)]
            got = call(f1["jac_t_vec"][1], "x", m1, pt, x=xd)
            want = call(f0["jac_t_vec"][1], "x", m0, pt, x=xd)
            for i, (x, y) in enumerate(zip(got, want)):
                for k in range(n):
                    t.compare(H.pyrt.tangent(x, k), H.pyrt.tangent(y, k), 1e-12,
                              "%s dual jac_t_vec[%d]/dx%d" % (label, i, k), singular=ill)
    print("  %d models with map rows" % rows)
    return t.report()


def check_linear_codegen(corpus):
    """Code generation of the 481-state LLG dipole model with the map."""
    import codegen_cppODE as CO
    import cppde_model as CM
    t = H.Tally("linear map: codegen, 481-state LLG")
    rhs, params = H.llg_dipole(160)
    t0 = time.time()
    m = CM.model_for(rhs, params, [])
    print("  parse and map: %.2f s, %d rows, %d coefficients, W_ent %d"
          % (time.time() - t0, m.nlin, m.linmap.nnz(), m.entry_work()[0]))
    runs = [("forward dual, dense", dict(num_type="cppde::dual<double, 0>",
                                         ad_level=1, arena=True, sparse=False)),
            ("reverse rb4, double", dict(num_type="double", sparse=False,
                                         emit_contractions=True, emit_jvp=True)),
            ("forward double, sparse", dict(num_type="double", sparse=True))]
    for label, kw in runs:
        t0 = time.time()
        res = CO.generate_ode_cpp(rhs, params, **kw)
        dt = time.time() - t0
        size = sum(len(res[k]) for k in ("ode_code", "jac_code", "adj_code", "data_code"))
        st = res["codegen_stats"]
        print("  %s: %.2f s, %.0f kB (%.0f kB tables), strategy %s"
              % (label, dt, size / 1e3, len(res["data_code"]) / 1e3, st.get("strategy")))
        t.check(dt < 30, "%s took %.1f s" % (label, dt))
    return t.report()



# ---------------------------------------------------------------------------
# Phase 7: loops over statements of equal structure
# ---------------------------------------------------------------------------

def _set_vectorise(value):
    os.environ["CPPDE_VECTORISE"] = str(value)


def _model_statements(m, emit_jvp=True):
    """name -> (kind, statements, args) of every model function."""
    import cppde_model as CM
    out = {"f": ("x", CM.ode_statements(m), ["x", "dxdt", "t", "params", "F"])}
    for name, kind, stores in m.contraction_stores(emit_jvp):
        out[name] = (kind, CM.stores_prelude(m, stores, "double")
                     + [CM.block(m.g, stores, "_t")], ARGS)
    for strat in ("entries", "colour"):
        for sparse in (False, True):
            out[(strat, sparse)] = ("jac", CM.jacobian_statements(m, sparse, strat),
                                    JAC_ARGS)
    return out


def _count(stmts):
    n = 0
    for st in stmts:
        if st[0] == "table":
            continue
        n += 1
        if st[0] in ("block", "if"):
            n += _count(st[2])
        elif st[0] == "for":
            n += _count(st[3])
        elif st[0] == "loop":
            n += _count(st[4])
    return n


def _loops(stmts):
    n = 0
    for st in stmts:
        if st[0] == "for" and st[1] == "_q":
            n += 1
        if st[0] in ("block", "if"):
            n += _loops(st[2])
        elif st[0] in ("for", "loop"):
            n += _loops(st[-1])
    return n


def _run(m, key, kind, fn, pt):
    if kind == "jac":
        a, d, W = run_jac(fn, m, pt, key[1], m.n)
        vals = [v for row in a for v in row] + d
        if key[1]:
            xs = [v * 1.3 for v in pt.x]
            a2, _, _ = run_jac(fn, m, pt, True, m.n, W=W, x=xs, t=pt.t + 0.7)
            vals += [v for row in a2 for v in row]
        return vals
    if kind == "x" and key == "f":
        out = [0.0] * m.n
        fn(pt.x, out, pt.t, pt.p, pt.F)
        return out
    return call(fn, kind, m, pt)


def check_vectorise(corpus):
    """Every model function with and without loops, bit for bit."""
    t = H.Tally("loops: functions with and without, identical")
    ta = H.Tally("loops: with sums as loops, 1e-13")
    rng = random.Random(51)
    looped = 0
    try:
        for label, rhs, params, forcings in list(H.VECTOR_MODELS) + ode_models(corpus):
            try:
                m = H.OdeModel(rhs, params, forcings)
                _set_vectorise(0)
                plain = _model_statements(m)
                _set_vectorise("")
                vec = _model_statements(m)
            except Exception as e:
                t.skip("%s: %s" % (label, str(e)[:100]))
                continue
            nl = sum(_loops(v[1]) for v in vec.values())
            if not nl:
                continue
            looped += 1
            fns = {}
            for key in plain:
                kind, st0, args = plain[key]
                fns[key] = (kind, py_body(m, st0, args), py_body(m, vec[key][1], args))
            pts, _ = H.sample_points(m, m.rhs, rng, base=3, extra=6)
            for pt in pts:
                ill = H.ill_conditioned(m, pt)
                for key, (kind, f0, f1) in fns.items():
                    a = _run(m, key, kind, f0, pt)
                    b = _run(m, key, kind, f1, pt)
                    # a sum run as a loop, or reordered '+=' rows, round
                    # differently
                    exact = not _reordered(vec[key][1])
                    floor = cancel_floor(a)
                    for i, (x, y) in enumerate(zip(b, a)):
                        if exact:
                            same = (x == y) or (x != x and y != y)
                            t.check(same, "%s %s[%d]: loops %r plain %r"
                                    % (label, key, i, x, y))
                        else:
                            ta.compare(x, x if near(x, y, floor) else y, 1e-13,
                                       "%s %s[%d]: loops %r plain %r"
                                       % (label, key, i, x, y), singular=ill)
            print("  %-16s %4d states, %3d loops" % (label, m.n, nl))
    finally:
        _set_vectorise("")
    print("  %d models with loops" % looped)
    ok = t.report()
    return ta.report() and ok


def _reordered(stmts):
    """True if statements hold a sum loop or '+=' row calls."""
    for st in stmts:
        if st[0] == "loop" and st[1] == "_p":
            return True
        if st[0] == "store" and st[1][0] == "call":
            return True
        if st[0] in ("block", "if") and _reordered(st[2]):
            return True
        if st[0] in ("for", "loop") and _reordered(st[-1]):
            return True
    return False


def check_vector_size(corpus):
    """Statements outside tables at N and 2N."""
    t = H.Tally("loops: statement count independent of N")
    pairs = [("rd1d", H.rd1d(64), H.rd1d(128)),
             ("rd1d-periodic", H.rd1d(64, True), H.rd1d(128, True)),
             ("brusselator2d", H.brusselator2d(24, 24), H.brusselator2d(48, 48)),
             ("llg", H.llg_dipole(24), H.llg_dipole(48))]
    for label, (r1, p1), (r2, p2) in pairs:
        counts = []
        for rhs, params in ((r1, p1), (r2, p2)):
            t0 = time.time()
            m = H.OdeModel(rhs, params, [])
            st = _model_statements(m)
            counts.append({k: _count(v[1]) for k, v in st.items()})
            print("  %-14s %5d states: %.1f s" % (label, m.n, time.time() - t0))
        for k in counts[0]:
            t.check(counts[0][k] == counts[1][k], "%s %s: %d vs %d statements"
                    % (label, k, counts[0][k], counts[1][k]))
    return t.report()


# ---------------------------------------------------------------------------
# Phase 8: CVODE
# ---------------------------------------------------------------------------

CV_ARGS = ["x", "t", "params", "F", "J", "ydot_arr", "yS_arr", "_Mp", "ySdot_arr",
           "lam", "lamdot", "qdot", "gout", "ud"]


def _cv_call(fn, pt, x=None, **kw):
    a = {k: None for k in CV_ARGS}
    a.update(x=pt.x if x is None else x, t=pt.t, params=pt.p, F=pt.F)
    a.update(kw)
    r = fn(*[a[k] for k in CV_ARGS])
    return r


def check_cvode(corpus):
    """CVODE callbacks against SymPy: f, +J, the sensitivity right-hand side,
    the adjoint and its quadrature, event and root partials."""
    import cppde_model as CM
    t = H.Tally("cvode: callbacks vs SymPy")
    rng = random.Random(61)
    specs = [(l, r, p, f) for l, r, p, f in ode_models(corpus)]
    specs += [(l, r, p, f) for l, r, p, f in H.LINEAR_MODELS]
    for label, rhs, params, forcings in specs[:80] + specs[-3:]:
        try:
            m = H.OdeModel(rhs, params, forcings)
            m.exprs = [str(rhs[s]) for s in m.states]
            exprs, J, dd, names = jac_reference(m)
            local = H.sym_table(m.names())
            np_ = len(m.params)
            psyms = [local[p_] for p_ in m.params]
            ev_f = H.Evaluator(exprs, names)
            flatJ = [(i, j) for i in range(m.n) for j in range(m.n) if J[i][j] != 0]
            ev_J = H.Evaluator([J[i][j] for i, j in flatJ], names)
            dfdp = [[H.replace_dirac(sp.diff(e, ps)) for ps in psyms] for e in exprs]
            ev_P = H.Evaluator([d for row in dfdp for d in row], names)
            f = py_body(m, CM.cvode_rhs_statements(m), CV_ARGS)
            jac = {(s_, sp_): py_body(m, CM.jacobian_statements(m, sp_, s_, cvode=True),
                                      CV_ARGS)
                   for s_ in ("entries", "colour") for sp_ in (False, True)}
            sens_p = set(range(np_))
            sens = py_body(m, CM.cvode_sens_statements(m, sens_p), CV_ARGS)
            xb, qb = CM.cvode_adjoint_statements(m)
            adj = py_body(m, xb, CV_ARGS)
            quad = py_body(m, qb, CV_ARGS)
            nnz = sum(len(r) for r in m.pattern())
        except Exception as e:
            t.skip("%s: %s" % (label, str(e)[:100]))
            continue
        pts, _ = H.sample_points(m, m.rhs, rng, base=3, extra=9)
        for pt in pts:
            env = pt.env(m)
            ill = H.ill_conditioned(m, pt)
            fv = ev_f(env)
            Jv = dict(zip(flatJ, ev_J(env)))
            Pv = ev_P(env)
            got = [0.0] * m.n
            _cv_call(f, pt, ydot_arr=got)
            for i in range(m.n):
                t.compare(got[i], fv[i], 1e-12, "%s f[%d]" % (label, i), singular=ill)
            for (s_, sp_), fn in jac.items():
                if sp_:
                    Jm = H.pyrt.SunSparse(m.n, nnz)
                    _cv_call(fn, pt, J=Jm)
                    dense = Jm.dense()
                else:
                    Jm = H.pyrt.DenseMatrix(m.n)
                    _cv_call(fn, pt, J=Jm)
                    dense = Jm.data
                for i in range(m.n):
                    for j in range(m.n):
                        b = Jv.get((i, j), 0.0)
                        if s_ == "colour" and any(not math.isfinite(Jv.get((i, c), 0.0))
                                                  for c in range(m.n)):
                            b = math.nan
                        t.compare(dense[i][j], b, 1e-12, "%s %s/%s J[%d,%d]: %r vs %r"
                                  % (label, s_, "csc" if sp_ else "dense", i, j,
                                     dense[i][j], b), singular=ill)
            yS = [H.draw(rng, 1) for _ in range(m.n)]
            Mp = [H.draw(rng, 1) for _ in range(np_)]
            got = [0.0] * m.n
            _cv_call(sens, pt, yS_arr=yS, _Mp=Mp, ySdot_arr=got)
            for i in range(m.n):
                want = sum(Jv.get((i, j), 0.0) * yS[j] for j in range(m.n))
                want += sum(Pv[i * np_ + k] * Mp[k] for k in range(np_))
                t.compare(got[i], want, 1e-11, "%s sens[%d]: %r vs %r"
                          % (label, i, got[i], want), singular=ill)
            lamdot = [0.0] * m.n
            qdot = [0.0] * np_
            _cv_call(adj, pt, lam=pt.lam, lamdot=lamdot)
            _cv_call(quad, pt, lam=pt.lam, qdot=qdot)
            for j in range(m.n):
                want = -sum(Jv.get((i, j), 0.0) * pt.lam[i] for i in range(m.n))
                t.compare(lamdot[j], want, 1e-11, "%s lamdot[%d]: %r vs %r"
                          % (label, j, lamdot[j], want), singular=ill)
            for k in range(np_):
                want = -sum(Pv[i * np_ + k] * pt.lam[i] for i in range(m.n))
                t.compare(qdot[k], want, 1e-11, "%s qdot[%d]: %r vs %r"
                          % (label, k, qdot[k], want), singular=ill)
    return t.report() and _check_cvode_events(corpus)


def _check_cvode_events(corpus):
    import cppde_model as CM
    t = H.Tally("cvode: event partials vs SymPy")
    rng = random.Random(62)
    for label, rhs, params, forcings, rows, has_rhs in event_specs(corpus):
        try:
            m = H.OdeModel(rhs, params, forcings)
            local = H.sym_table(m.names())
            ev = CM.CvodeEvent(m)
            names = m.names()
            xs = [local[s_] for s_ in m.states]
            psyms = [local[p_] for p_ in m.params]
            tsym = local["time"]
            items = []
            for r in rows:
                if not CM.valid_value(r["value"]):
                    continue
                var = str(r["var"])
                vi = m.states.index(var)
                h = m.parse(str(r["value"]))
                meth = str(r["method"]).lower() if CM.valid_value(r["method"]) else "replace"
                gn = {"replace": h, "add": m.g.add(m.g.state(vi), h),
                      "multiply": m.g.mul(m.g.state(vi), h)}.get(meth, h)
                hs = H.legacy_parse(str(r["value"]), local)
                gs = {"replace": hs, "add": xs[vi] + hs, "multiply": xs[vi] * hs}.get(meth, hs)
                items.append((gn, gs))
                for col in ("time", "root"):
                    if CM.valid_value(r[col]):
                        items.append((m.parse(str(r[col])),
                                      H.legacy_parse(str(r[col]), local)))
            fns = []
            for node, sym in items:
                refs = [sym] + [H.replace_dirac(sp.diff(sym, v)) for v in xs + psyms]
                refs.append(H.replace_dirac(sp.diff(sym, tsym)))
                blocks = [ev.value(node)]
                cx = dict(ev.cases_x(node))
                cp = dict(ev.cases_p(node))
                blocks += [cx.get(j) for j in range(m.n)]
                blocks += [cp.get(k) for k in range(len(m.params))]
                blocks.append(ev.partial_t(node))
                fns.append(([py_body(m, b, CV_ARGS) if b else None for b in blocks],
                            H.Evaluator(refs, names)))
        except Exception as e:
            t.skip("%s: %s" % (label, str(e)[:100]))
            continue
        for _ in range(4):
            pt = H.Point(m, rng, 0)
            env = pt.env(m)
            for bodies, evr in fns:
                want = evr(env)
                for q, (fn, w) in enumerate(zip(bodies, want)):
                    got = 0.0 if fn is None else _cv_call(fn, pt)
                    t.compare(got, w, 1e-11, "%s event part %d: %r vs %r"
                              % (label, q, got, w))
    return t.report()


CHECKS = {
    "parse": check_parse,
    "fallback": check_fallback,
    "determinism": check_determinism,
    "syntax": check_syntax,
    "adjoint": lambda c: check_contractions(c),
    "adjoint-fd": lambda c: check_fd(c),
    "jacobian": lambda c: check_jacobian(c),
    "jacobian-reuse": lambda c: check_jacobian_reuse(c),
    "jacobian-dual": lambda c: check_jacobian_dual(c),
    "events": check_events,
    "fun": check_fun,
    "linear": check_linear,
    "linear-codegen": check_linear_codegen,
    "vectorise": check_vectorise,
    "vector-size": check_vector_size,
    "model-syntax": check_model_syntax,
    "cvode": check_cvode,
    "cvode-syntax": check_cvode_syntax,
}


def main(argv):
    names = argv or list(CHECKS)
    corpus = load_corpus()
    if not corpus:
        print("no corpus at", CORPUS)
    ok = True
    for name in names:
        t0 = time.time()
        ok = CHECKS[name](corpus) and ok
        print("  (%s: %.1f s)" % (name, time.time() - t0))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
