"""ODE model graphs and the stores of the generated model functions.

Flat parameter layout: params[i] is the initial value of state i,
params[n_states + k] parameter k. With a constant linear map (see
cppde_struct), row r of C x is `_lin[r]` and of C v is `_linv[r]`.
"""

import os

import cppde_emit as em
import cppde_graph as cg
import cppde_struct as cs


def _bit_list(b):
    out = []
    while b:
        low = b & -b
        out.append(low.bit_length() - 1)
        b ^= low
    return out


def as_list(x):
    if x is None:
        return []
    if isinstance(x, str):
        return [x]
    return list(x)


class OdeModel:
    """Parsed right-hand side of an ODE model.

    Args:
        rhs: dict state name -> expression string (order = state order).
        params: parameter names; a name `<state>_0` reads the initial value.
        forcings: forcing names.
        linear: row threshold of the linear map (cppde_struct.linear_min_terms).

    Attributes:
        g: cppde_graph.Graph; ad: cppde_graph.AD.
        rhs_plain: root node per state as parsed.
        rhs: the same with long linear sums as map rows.
        linmap: cppde_struct.LinMap or None; nlin its number of rows.
    """

    def __init__(self, rhs, params, forcings=None, linear=None):
        self.states = list(rhs.keys())
        self.params = as_list(params)
        self.forcings = as_list(forcings)
        self.n = len(self.states)
        g = self.g = cg.Graph()
        sym = {}
        for j, f in enumerate(self.forcings):
            sym[f] = g.forcing(j)
        for k, p in enumerate(self.params):
            sym[p] = g.param(k)
        for i, s in enumerate(self.states):
            sym[s] = g.state(i)
        self.init_names = {}
        for i, s in enumerate(self.states):
            sym[s + "_0"] = g.init(i)
            self.init_names[s + "_0"] = i
        sym["time"] = g.time()
        self.sym = sym
        self.parser = cg.Parser(g, sym)
        self.ad = cg.AD(g)
        self.rhs_plain = [self.parser.parse(str(rhs[s]), label=s) for s in self.states]
        self.rhs, self.linmap, self.lin_expand = cs.linear_rows(
            g, self.rhs_plain, self.n, cs.linear_min_terms(linear))
        self.nlin = len(self.linmap) if self.linmap is not None else 0
        self._cache = {}

    def parse(self, text, label=None):
        return self.parser.parse(text, label=label)

    def names(self):
        """Names the expressions may read, without duplicates."""
        out = []
        for s in (self.states + self.params + self.forcings
                  + list(self.init_names) + ["time"]):
            if s not in out:
                out.append(s)
        return out

    def slot(self, style="cpp", **kw):
        return em.model_slots(self.n, style=style, **kw)

    def lam(self):
        return [self.g.vec("lam", i) for i in range(self.n)]

    def vvec(self):
        return [self.g.vec("v", k) for k in range(self.n)]

    def sc(self):
        return self.g.vec("sc", None)

    def linrows(self):
        return [self.g.leaf(cg.LINROW, r) for r in range(self.nlin)]

    def linv(self):
        return [self.g.vec("_linv", r) for r in range(self.nlin)]

    # -- derived expressions -------------------------------------------------

    def cached(self, key, fn):
        hit = self._cache.get(key)
        if hit is None:
            hit = self._cache[key] = fn()
        return hit

    def dfdt(self):
        """df_i/dt including forcing rates."""
        return self.cached("dfdt", lambda: self.ad.time_derivative(self.rhs))

    def jv(self):
        """(J v)_i with v = VEC('v')."""
        v = self.vvec()
        seeds = {self.g.state(k): v[k] for k in range(self.n)}
        seeds.update(zip(self.linrows(), self.linv()))
        return self.cached("jv", lambda: self.ad.jvp(self.rhs, seeds))

    def pullback(self, roots):
        """(state, parameter, map row) adjoints of roots against lam. The
        state adjoints lack C^T times the row adjoints."""
        g = self.g
        adj = self.ad.vjp(roots, self.lam(), cg.F_STATE | cg.F_PARAM | cg.F_LINROW)
        out = {cg.STATE: {}, cg.PARAM: {}, cg.LINROW: {}}
        for leaf, a in adj.items():
            k = g.attr[leaf]
            out[k[0]][k[1]] = a
        return out[cg.STATE], out[cg.PARAM], out[cg.LINROW]

    # -- Jacobian -------------------------------------------------------------

    def pattern(self):
        """Per state, the sorted columns of its Jacobian row."""
        return self.cached("pattern", self._pattern)

    def _pattern(self):
        g = self.g
        out = []
        for r in self.rhs:
            cols = _bit_list(g.dbits(r, cg.STATE))
            rows = _bit_list(g.dbits(r, cg.LINROW)) if self.nlin else []
            if rows:
                cols = set(cols)
                for q in rows:
                    cols.update(self.linmap.cols(q))
                cols = sorted(cols)
            out.append(cols)
        return out

    def _diff_bits(self, n):
        g = self.g
        b = g.dbits(n, cg.STATE)
        return (b, g.dbits(n, cg.LINROW)) if self.nlin else (b, 0)

    def _diff_order(self):
        return [n for n in self.g.topo(self.rhs) if any(self._diff_bits(n))]

    def entry_work(self):
        """(W_ent, E_s): derivative nodes of the entry accumulation and
        state-dependent nodes of the right-hand side."""
        def work():
            order = self._diff_order()
            w = 0
            for n in order:
                b, q = self._diff_bits(n)
                w += bin(b).count("1") + bin(q).count("1")
            return w, len(order)
        return self.cached("entry_work", work)

    def entries(self):
        """dict (i, j) -> node of the nonzero entries of df/dx with the map
        rows held fixed."""
        return self.cached("entries", self._entries)[0]

    def lin_entries(self):
        """dict (i, r) -> node of df_i / d(C x)_r."""
        return self.cached("entries", self._entries)[1]

    def _entries(self):
        g = self.g
        n_x = self.n
        order = self._diff_order()
        roots = set(self.rhs)
        kids = {}
        refs = {}
        for n in order:
            ks = [a for a in g._dkids(n) if any(self._diff_bits(a))]
            kids[n] = ks
            for a in ks:
                refs[a] = refs.get(a, 0) + 1
        grad = {}
        for n in order:
            o = g.op[n]
            if o == cg.LEAF:
                at = g.attr[n]
                grad[n] = {at[1] if at[0] == cg.STATE else n_x + at[1]: g.ONE}
                continue
            acc = {}
            a = g.args[n]
            if o == cg.ADD:
                for c, t in zip(g.attr[n][1], a):
                    for j, d in grad.get(t, {}).items():
                        acc.setdefault(j, []).append(g.scale(c, d))
            elif o == cg.SELECT:
                cnd = a[0]
                for pos in (1, 2):
                    for j, d in grad.get(a[pos], {}).items():
                        br = (d, g.ZERO) if pos == 1 else (g.ZERO, d)
                        acc.setdefault(j, []).append(g.select(cnd, *br))
            else:
                for pos, part in self.ad.partials(n):
                    if g.is_zero(part):
                        continue
                    for j, d in grad.get(a[pos], {}).items():
                        acc.setdefault(j, []).append(g.mul(part, d))
            out = {}
            for j in sorted(acc):
                terms = acc[j]
                s = terms[0] if len(terms) == 1 else g.sum(terms)
                if not g.is_zero(s):
                    out[j] = s
            grad[n] = out
            for c in kids[n]:
                refs[c] -= 1
                if refs[c] == 0 and c not in roots:
                    grad.pop(c, None)
        a_ent, g_ent = {}, {}
        for i, r in enumerate(self.rhs):
            for j, d in grad.get(r, {}).items():
                if j < n_x:
                    a_ent[(i, j)] = d
                else:
                    g_ent[(i, j - n_x)] = d
        return a_ent, g_ent

    def colouring(self):
        """(colour per column or -1, number of colours) of a greedy
        distance-2 colouring of the column intersection graph."""
        def colour():
            pat = self.pattern()
            cols = [[] for _ in range(self.n)]
            for i, row in enumerate(pat):
                for j in row:
                    cols[j].append(i)
            order = sorted(range(self.n), key=lambda j: (-len(cols[j]), j))
            col = [-1] * self.n
            used = [set() for _ in range(self.n)]
            for j in order:
                if not cols[j]:
                    continue
                forbidden = set()
                for i in cols[j]:
                    forbidden |= used[i]
                c = 0
                while c in forbidden:
                    c += 1
                col[j] = c
                for i in cols[j]:
                    used[i].add(c)
            return col, max(col) + 1 if any(c >= 0 for c in col) else 0
        return self.cached("colouring", colour)

    def jacobian_strategy(self):
        """'entries' or 'colour'; CPPDE_JAC forces one."""
        forced = os.environ.get("CPPDE_JAC", "").strip().lower()
        if forced in ("entries", "colour", "color"):
            return "entries" if forced == "entries" else "colour"
        w, es = self.entry_work()
        if w <= max(4 * es, 50000):
            return "entries"
        if self.colouring()[1] <= self.n / 4:
            return "colour"
        return "entries"

    # -- stores of the contractions -----------------------------------------

    def contraction_stores(self, emit_jvp):
        """Stores of adjoint_terms, as [(name, kind, stores)].

        kind is 'x' (writes out[j] after out.assign(n, 0)), 'p' (adds
        sc*(...) to out[n_states + k]) or 'dot' (returns a scalar).
        """
        g = self.g
        sc = self.sc()
        n = self.n
        pairs = [("jac_t_vec", "dfdp_t_vec_axpy", self.rhs)]
        if emit_jvp:
            pairs.append(("jvp_x_t_vec", "jvp_p_t_vec_axpy", self.jv()))
            pairs.append(("dfdt_x_t_vec", "dfdt_p_t_vec_axpy", self.dfdt()))
        out = []
        for xname, pname, roots in pairs:
            xs, ps, ls = self.pullback(roots)
            stores = [(("vec", "out", j), xs[j], "=") for j in sorted(xs)]
            stores += [(lin_axpy(r), ls[r], "+=") for r in sorted(ls)]
            out.append((xname, "x", stores))
            out.append((pname, "p", [(("vec", "out", n + k), g.mul(sc, ps[k]), "+=")
                                     for k in sorted(ps)]))
        lam = self.lam()
        dot = g.sum([g.mul(d, lam[i]) for i, d in enumerate(self.dfdt())
                     if not g.is_zero(d)])
        out.append(("dfdt_dot", "dot", [(("ret",), dot, "=")]))
        return out


# ---------------------------------------------------------------------------
# Statements of ode_system and jacobian
# ---------------------------------------------------------------------------

_CLOCKED = (cg.F_STATE | cg.F_TIME | cg.F_FORCING | cg.F_FRATE | cg.F_VEC
            | cg.F_LINROW)


_AXPY_VEC = {"cpp": "cppde::axpy_row_dense(linmap_, {0}, %s, out);",
             "py": "R.axpy_row_dense(linmap_, {0}, %s, out)"}
_AXPY_MAT = {"cpp": "cppde::axpy_row_dense(linmap_, {0}, %s, J, {1});",
             "py": "R.axpy_row_dense(linmap_, {0}, %s, J, {1})"}
_AXPY_IDX = {"cpp": "cppde::axpy_row_idx(linmap_, {0}, %s, W.Ax, _gdst + {1});",
             "py": "R.axpy_row_idx(linmap_, {0}, %s, W.Ax, _gdst, {1})"}


def lin_axpy(r, i=None):
    """Call target adding expr * C[r, .] to out, or to row i of J."""
    if i is None:
        return ("call", _AXPY_VEC, (r,))
    return ("call", _AXPY_MAT, (r, i))


def lin_prelude(model, roots, T, src="x", out="_lin", decl=True):
    """Statements computing `out` = C src if a node of `roots` reads it."""
    g = model.g
    if not model.nlin:
        return []
    if out == "_lin":
        if not any(g.flags[r] & cg.F_LINROW for r in roots):
            return []
    elif not any(g.op[n] == cg.LEAF and g.attr[n][:2] == (cg.VEC, out)
                 for n in g.topo(roots) if g.flags[n] & cg.F_VEC):
        return []
    stmts = []
    if decl:
        stmts.append(("raw", {"cpp": "std::vector<%s> %s(%d);" % (T, out, model.nlin),
                              "py": "%s = [0.0] * %d" % (out, model.nlin)}))
    stmts.append(("raw", {"cpp": "cppde::apply(linmap_, %s, %s);" % (src, out),
                          "py": "R.apply(linmap_, %s, %s)" % (src, out)}))
    return stmts


def stores_prelude(model, stores, T):
    """lin_prelude for the rows of C x and C v read by `stores`."""
    roots = [n for _, n, _ in stores]
    return (lin_prelude(model, roots, T)
            + lin_prelude(model, roots, T, src="v", out="_linv"))


def ode_statements(model, T="double"):
    """Body of ode_system::operator(): dxdt[i] = f_i."""
    stores = [(("vec", "dxdt", i), r, "=") for i, r in enumerate(model.rhs)]
    return lin_prelude(model, model.rhs, T) + [block(model.g, stores, "_t")]


def block(g, stores, prefix, known=None):
    """Stores as one ('block', names, statements), with loops over classes
    of equal structure where cppde_struct.vector_block finds them."""
    loops = cs.vector_block(g, stores, prefix, known)
    if loops is not None:
        return ("block", dict(known or {}), loops)
    stmts, names = em.schedule(g, stores, prefix=prefix, known=known)
    return ("block", names, stmts)


def _seed_free_inputs(g, roots):
    """Non-atomic nodes without VEC leaves read by nodes with them."""
    out = []
    seen = set()
    for n in g.topo(roots):
        if not g.flags[n] & cg.F_VEC:
            continue
        for a in g.args[n]:
            if (a not in seen and not g.flags[a] & cg.F_VEC
                    and g.op[a] not in (cg.NUM, cg.NAMED, cg.BOOL, cg.LEAF)
                    and not g.is_bool(a)):
                seen.add(a)
                out.append(a)
    return out


def csc_layout(model, diagonal=True):
    """(rows, cols, ax, missing) of a sparse Jacobian: CSC-sorted pattern,
    with the missing diagonals if `diagonal`; ax[(i, j)] is the value index,
    missing the value indices of the added diagonals."""
    pat = model.pattern()
    pairs = [(i, j) for i, row in enumerate(pat) for j in row]
    have = {i for i, j in pairs if i == j}
    extra = [(i, i) for i in range(model.n) if i not in have] if diagonal else []
    allp = sorted(pairs + extra, key=lambda p: (p[1], p[0]))
    ax = {p: k for k, p in enumerate(allp)}
    return ([p[0] for p in allp], [p[1] for p in allp], ax,
            [ax[p] for p in extra])


_SUN_DENSE = {"cpp": "SM_ELEMENT_D(J, {0}, {1}) = %s;", "py": "J[{0}, {1}] = %s"}
_SUN_AXPY = {"cpp": "cppde::axpy_row_dense(linmap_, {0}, %s, _Jw, {1});",
             "py": "R.axpy_row_dense(linmap_, {0}, %s, _Jw, {1})"}
_SUN_AXPY_IDX = {"cpp": "cppde::axpy_row_idx(linmap_, {0}, %s, data, _gdst + {1});",
                 "py": "R.axpy_row_idx(linmap_, {0}, %s, data, _gdst, {1})"}


def jacobian_statements(model, sparse, strategy, T="double", cvode=False,
                        scoped=False):
    """Body of jacobian::operator(), writing -J into J or W and df/dt into dfdt,
    or with `cvode` of jac_fn, writing +J into the SUNMatrix. With a linear map
    J = A + G C; `scoped` binds the output tangents and opens an arena scope.
    """
    if cvode:
        return _cvode_jacobian(model, sparse, strategy)
    g = model.g
    n = model.n
    zero = {"cpp": "%s(0)" % T, "py": "0.0"}
    stmts = []
    if sparse:
        rows, cols, ax, missing = csc_layout(model)
        stmts.append(("raw", {"cpp": "const bool _init_consts = !W.pattern_built;",
                              "py": "_init_consts = not W.pattern_built"}))
        stmts.append(("if", {"cpp": "!W.pattern_built", "py": "not W.pattern_built"}, [
            ("table", "_rows", rows), ("table", "_cols", cols),
            ("raw", {"cpp": "W.build_pattern(%d, %d, _rows, _cols);" % (n, len(rows)),
                     "py": "W.build_pattern(%d, %d, _rows, _cols)" % (n, len(rows))})]))
        target = lambda i, j: ("vec", "W.Ax", ax[(i, j)])
    else:
        pairs = [(i, j) for i, row in enumerate(model.pattern()) for j in row]
        if 2 * len(pairs) > n * n:
            stmts.append(("raw", {"cpp": "J.set_zero();", "py": "J.set_zero()"}))
        elif pairs:
            stmts += [("table", "_dr", [p[0] for p in pairs]),
                      ("table", "_dc", [p[1] for p in pairs]),
                      ("for", "_k", str(len(pairs)),
                       [("assign", ("mat", "J", "_dr[_k]", "_dc[_k]"), zero)])]
        target = lambda i, j: ("mat", "J", i, j)
    if scoped:
        stmts += [("raw", {"cpp": "cppde::ad_traits::arm_outputs(%s, dfdt, x, params);"
                                  % ("W.Ax" if sparse else "J.data"), "py": None}),
                  ("raw", {"cpp": "cppde::dual_arena::scope _jac_arena_scope;",
                           "py": None})]
    dfdt = [(("vec", "dfdt", i), d, "=") for i, d in enumerate(model.dfdt())]

    if strategy == "entries":
        ent = model.entries()
        gent = model.lin_entries()
        gstores, gpos, gdst = [], set(), []
        for (i, r) in sorted(gent):
            cols_r = model.linmap.cols(r)
            if sparse:
                call = ("call", _AXPY_IDX, (r, len(gdst)))
                gdst += [ax[(i, j)] for j in cols_r]
            else:
                call = lin_axpy(r, i)
            gpos.update((i, j) for j in cols_r)
            gstores.append((call, g.neg(gent[(i, r)]), "+="))
        const, rest = [], []
        for (i, j) in sorted(ent):
            st = (target(i, j), g.neg(ent[(i, j)]), "=")
            if sparse and not (g.flags[ent[(i, j)]] & _CLOCKED) and (i, j) not in gpos:
                const.append(st)
            else:
                rest.append(st)
        if sparse and gstores:
            stmts.append(("table", "_gdst", gdst))
            zeros = sorted(ax[p] for p in gpos if p not in ent)
            if zeros:
                stmts += [("table", "_gz", zeros),
                          ("for", "_k", str(len(zeros)),
                           [("assign", ("vec", "W.Ax", "_gz[_k]"), zero)])]
        if const:
            stmts.append(("if", {"cpp": "_init_consts", "py": "_init_consts"},
                          [block(g, const, "_c")]))
        stores = rest + gstores + dfdt
        stmts += stores_prelude(model, stores, T)
        stmts.append(block(g, stores, "_t"))
    else:
        col, chi = model.colouring()
        cptr, cidx = [0], []
        for c in range(chi):
            cidx += [j for j in range(n) if col[j] == c]
            cptr.append(len(cidx))
        sptr, srow, sdst = [0], [], []
        pat = model.pattern()
        for c in range(chi):
            for i, row in enumerate(pat):
                for j in row:
                    if col[j] == c:
                        srow.append(i)
                        sdst.append(ax[(i, j)] if sparse else j)
            sptr.append(len(srow))
        jv = model.jv()
        tstores = [(("vec", "_jv", i), t, "=") for i, t in enumerate(jv)
                   if not g.is_zero(t)]
        # values that do not depend on the seeds are computed once, in _jo
        bound = _seed_free_inputs(g, [t for _, t, _ in tstores])
        sub = {b: g.vec("_jo", k) for k, b in enumerate(bound)}
        inner = block(g, [(tg, cs.raw_substitute(g, t, sub), op)
                          for tg, t, op in tstores], "_t")
        itables = [st for st in inner[2] if st[0] == "table"]
        inner = ("block", inner[1], [st for st in inner[2] if st[0] != "table"])
        outer = block(g, [(("vec", "_jo", k), b, "=") for k, b in enumerate(bound)]
                      + dfdt, "_o")
        write = (("vec", "W.Ax", "_sdst[_k]") if sparse
                 else ("mat", "J", "_srow[_k]", "_sdst[_k]"))
        colour = [("table", "_cptr", cptr), ("table", "_cidx", cidx),
                  ("table", "_sptr", sptr), ("table", "_srow", srow),
                  ("table", "_sdst", sdst),
                  ("raw", {"cpp": "std::vector<double> v(%d, 0.0);" % n,
                           "py": "v = [0.0] * %d" % n}),
                  ("raw", {"cpp": "std::vector<%s> _jv(%d);" % (T, n),
                           "py": "_jv = [0.0] * %d" % n})]
        colour += lin_prelude(model, [t for _, t, _ in tstores + dfdt], T)
        linv = lin_prelude(model, [t for _, t, _ in tstores], "double",
                           src="v", out="_linv", decl=False)
        if linv:
            colour.append(("raw", {"cpp": "std::vector<double> _linv(%d);" % model.nlin,
                                   "py": "_linv = [0.0] * %d" % model.nlin}))
        if bound:
            colour.append(("array", "_jo", len(bound)))
        colour += [outer] + itables
        colour.append(("for", "_c", str(chi), [
            ("loop", "_k", "_cptr[_c]", "_cptr[_c + 1]",
             [("assign", ("vec", "v", "_cidx[_k]"), "1.0")])]
            + linv
            + [inner]
            + [("loop", "_k", "_sptr[_c]", "_sptr[_c + 1]",
                [("assign", write, "-_jv[_srow[_k]]")]),
               ("loop", "_k", "_cptr[_c]", "_cptr[_c + 1]",
                [("assign", ("vec", "v", "_cidx[_k]"), "0.0")])]))
        stmts += colour
    if sparse:
        stmts += [("assign", ("vec", "W.Ax", k), zero) for k in missing]
    return stmts


def _cvode_jacobian(model, sparse, strategy):
    g = model.g
    n = model.n
    stmts = []
    if sparse:
        rows, cols, ax, _ = csc_layout(model, diagonal=False)
        colptr = [0] * (n + 1)
        for c in cols:
            colptr[c + 1] += 1
        for k in range(n):
            colptr[k + 1] += colptr[k]
        stmts += [
            ("raw", {"cpp": "sunindextype* indexptrs = SUNSparseMatrix_IndexPointers(J);",
                     "py": "indexptrs = J.indexptrs"}),
            ("raw", {"cpp": "sunindextype* indexvals = SUNSparseMatrix_IndexValues(J);",
                     "py": "indexvals = J.indexvals"}),
            ("raw", {"cpp": "sunrealtype* data = SUNSparseMatrix_Data(J);",
                     "py": "data = J.data"}),
            ("table", "_colptr", colptr), ("table", "_rowval", rows),
            ("for", "_k", str(n + 1),
             [("assign", ("vec", "indexptrs", "_k"), "_colptr[_k]")]),
            ("for", "_k", str(len(rows)),
             [("assign", ("vec", "indexvals", "_k"), "_rowval[_k]")])]
        target = lambda i, j: ("vec", "data", ax[(i, j)])
    else:
        target = lambda i, j: ("call", _SUN_DENSE, (i, j))
        if model.nlin:
            stmts.append(("raw", {
                "cpp": "struct { SUNMatrix M; sunrealtype& operator()(int i, int j) "
                       "{ return SM_ELEMENT_D(M, i, j); } } _Jw{J};",
                "py": "_Jw = J"}))
    if strategy == "entries":
        ent = model.entries()
        gent = model.lin_entries()
        gstores, gpos, gdst = [], set(), []
        for (i, r) in sorted(gent):
            cols_r = model.linmap.cols(r)
            if sparse:
                call = ("call", _SUN_AXPY_IDX, (r, len(gdst)))
                gdst += [ax[(i, j)] for j in cols_r]
            else:
                call = ("call", _SUN_AXPY, (r, i))
            gpos.update((i, j) for j in cols_r)
            gstores.append((call, gent[(i, r)], "+="))
        if sparse and gstores:
            stmts.append(("table", "_gdst", gdst))
            zeros = sorted(ax[p] for p in gpos if p not in ent)
            if zeros:
                stmts += [("table", "_gz", zeros),
                          ("for", "_k", str(len(zeros)),
                           [("assign", ("vec", "data", "_gz[_k]"), "0.0")])]
        stores = [(target(i, j), ent[(i, j)], "=") for (i, j) in sorted(ent)]
        stores += gstores
        stmts += stores_prelude(model, stores, "double")
        stmts.append(block(g, stores, "_t"))
        return stmts
    col, chi = model.colouring()
    cptr, cidx = [0], []
    for c in range(chi):
        cidx += [j for j in range(n) if col[j] == c]
        cptr.append(len(cidx))
    sptr, srow, sdst = [0], [], []
    pat = model.pattern()
    for c in range(chi):
        for i, row in enumerate(pat):
            for j in row:
                if col[j] == c:
                    srow.append(i)
                    sdst.append(ax[(i, j)] if sparse else j)
        sptr.append(len(srow))
    jv = model.jv()
    tstores = [(("vec", "_jv", i), t, "=") for i, t in enumerate(jv) if not g.is_zero(t)]
    bound = _seed_free_inputs(g, [t for _, t, _ in tstores])
    sub = {b: g.vec("_jo", k) for k, b in enumerate(bound)}
    inner = block(g, [(tg, cs.raw_substitute(g, t, sub), op) for tg, t, op in tstores], "_t")
    itables = [st for st in inner[2] if st[0] == "table"]
    inner = ("block", inner[1], [st for st in inner[2] if st[0] != "table"])
    outer = block(g, [(("vec", "_jo", k), b, "=") for k, b in enumerate(bound)], "_o")
    write = (("vec", "data", "_sdst[_k]") if sparse else
             ("var", {"cpp": "SM_ELEMENT_D(J, _srow[_k], _sdst[_k])",
                      "py": "J[_srow[_k], _sdst[_k]]"}))
    stmts += [("table", "_cptr", cptr), ("table", "_cidx", cidx),
              ("table", "_sptr", sptr), ("table", "_srow", srow),
              ("table", "_sdst", sdst),
              ("raw", {"cpp": "std::vector<double> v(%d, 0.0);" % n,
                       "py": "v = [0.0] * %d" % n}),
              ("raw", {"cpp": "std::vector<double> _jv(%d);" % n,
                       "py": "_jv = [0.0] * %d" % n})]
    stmts += lin_prelude(model, [t for _, t, _ in tstores], "double")
    linv = lin_prelude(model, [t for _, t, _ in tstores], "double",
                       src="v", out="_linv", decl=False)
    if linv:
        stmts.append(("raw", {"cpp": "std::vector<double> _linv(%d);" % model.nlin,
                              "py": "_linv = [0.0] * %d" % model.nlin}))
    if bound:
        stmts.append(("array", "_jo", len(bound)))
    stmts += [outer] + itables
    stmts.append(("for", "_c", str(chi), [
        ("loop", "_k", "_cptr[_c]", "_cptr[_c + 1]",
         [("assign", ("vec", "v", "_cidx[_k]"), "1.0")])]
        + linv + [inner]
        + [("loop", "_k", "_sptr[_c]", "_sptr[_c + 1]",
            [("assign", write, "_jv[_srow[_k]]")]),
           ("loop", "_k", "_cptr[_c]", "_cptr[_c + 1]",
            [("assign", ("vec", "v", "_cidx[_k]"), "0.0")])]))
    return stmts


# ---------------------------------------------------------------------------
# CVODE callbacks
# ---------------------------------------------------------------------------

def cvode_rhs_statements(model):
    """rhs_fn body: ydot_arr[i] = f_i."""
    stores = [(("vec", "ydot_arr", i), r, "=") for i, r in enumerate(model.rhs)]
    return lin_prelude(model, model.rhs, "double") + [block(model.g, stores, "_t")]


def _init_param_index(model):
    """state index -> index in params of the parameter `<state>_0`."""
    out = {}
    for name, i in model.init_names.items():
        if name in model.params:
            out[i] = model.params.index(name)
    return out


def cvode_sens_statements(model, sens_params):
    """sens_rhs1_fn body: ySdot = J yS + (df/dp) Phi'[params, iS], as one
    JVP. State j is seeded with yS_arr[j], parameter k with _Mp[k]; `<state>_0`
    is read as its parameter row."""
    g = model.g
    seeds = {g.state(j): g.vec("yS_arr", j) for j in range(model.n)}
    for k in sens_params:
        seeds[g.param(k)] = g.vec("_Mp", k)
    for i, k in _init_param_index(model).items():
        if k in sens_params:
            seeds[g.init(i)] = g.vec("_Mp", k)
    seeds.update(zip(model.linrows(), model.linv()))
    tan = model.ad.jvp(model.rhs, seeds)
    stores = [(("vec", "ySdot_arr", i), t, "=") for i, t in enumerate(tan)]
    pre = lin_prelude(model, model.rhs + tan, "double")
    pre += lin_prelude(model, tan, "double", src="yS_arr", out="_linv")
    return pre + [block(g, stores, "_t")]


def cvode_adjoint_statements(model):
    """(adj_rhs_fn body, adj_quad_fn body): lamdot = -J' lam and
    qdot = -(df/dp)' lam."""
    g = model.g
    adj = model.ad.vjp(model.rhs, model.lam(),
                       cg.F_STATE | cg.F_PARAM | cg.F_INIT | cg.F_LINROW)
    xs, ps, ls, inits = {}, {}, {}, {}
    for leaf, a in adj.items():
        at = g.attr[leaf]
        {cg.STATE: xs, cg.PARAM: ps, cg.LINROW: ls, cg.INIT: inits}[at[0]][at[1]] = a
    for i, k in _init_param_index(model).items():
        if i in inits:
            ps[k] = g.add(ps[k], inits[i]) if k in ps else inits[i]
    x_stores = [(("vec", "lamdot", j), g.neg(xs[j]) if j in xs else g.ZERO, "=")
                for j in range(model.n)]
    x_stores += [(("call", {"cpp": "cppde::axpy_row_dense(linmap_, {0}, %s, lamdot);",
                            "py": "R.axpy_row_dense(linmap_, {0}, %s, lamdot)"}, (r,)),
                  g.neg(ls[r]), "+=") for r in sorted(ls)]
    q_stores = [(("vec", "qdot", k), g.neg(ps[k]) if k in ps else g.ZERO, "=")
                for k in range(len(model.params))]
    x_body = stores_prelude(model, x_stores, "double") + [block(g, x_stores, "_t")]
    q_body = stores_prelude(model, q_stores, "double") + [block(g, q_stores, "_a")]
    return x_body, q_body


class CvodeEvent:
    """Event and root expressions of the CVODE backend as statement blocks.

    Derivatives in parameters are indexed by position in model.params; a
    `<state>_0` parameter stands for the initial value it reads.
    """

    def __init__(self, model):
        self.m = model
        self.init_param = _init_param_index(model)

    def value(self, node):
        """Statements returning node."""
        return [block(self.m.g, [(("ret",), node, "=")], "_e")]

    def cases_x(self, node):
        xs, _, _ = _grad(self.m, node, cg.F_STATE)
        return [(j, self.value(xs[j])) for j in sorted(xs)]

    def cases_p(self, node):
        _, ps, inits = _grad(self.m, node, cg.F_PARAM | cg.F_INIT)
        acc = dict(ps)
        for i, d in inits.items():
            k = self.init_param.get(i)
            if k is not None:
                acc[k] = self.m.g.add(acc[k], d) if k in acc else d
        return [(k, self.value(acc[k])) for k in sorted(acc)]

    def partial_t(self, node):
        g = self.m.g
        return self.value(self.m.ad.jvp([node], {g.time(): g.ONE})[0])


def cvode_lines(model, stmts, indent):
    """C++ lines of CVODE statements (double)."""
    return cpp_statements_raw(model, stmts, Scalar("double"), indent=indent)


# ---------------------------------------------------------------------------
# C++ rendering
# ---------------------------------------------------------------------------

class Scalar:
    """C++ scalar of a generated model.

    Args:
        name: type spelling, e.g. 'cppde::dual<double, 0>'.
        ad_level: derivative layers (0, 1, 2).
        arena: tangents live in cppde::dual_arena.
    """

    def __init__(self, name, ad_level=0, arena=False):
        self.name = str(name)
        self.level = int(ad_level)
        self.arena = bool(arena)
        self.style = "ad" if self.level > 0 else "double"

    def arena_scope(self):
        """Arena scope line for a first-order dual, else none."""
        if self.arena and self.level == 1:
            return ["    cppde::dual_arena::scope _rhs_arena_scope;"]
        return []


def cpp_statements(model, stores, scalar, prefix="_t", indent="    "):
    stmts, names = em.schedule(model.g, stores, prefix=prefix)
    pr = em.Printer(model.g, model.slot("cpp"), style=scalar.style,
                    names=names, ad_level=scalar.level)
    return em.render_cpp(stmts, pr, scalar.name, indent=indent)


def _struct_head(name, T):
    return [
        "struct %s {" % name,
        "  std::vector<%s> params;" % T,
        "  std::vector<const cppde::PchipForcing<%s>*> F;" % T,
        "",
        "  %s(const std::vector<%s>& p_," % (name, T),
        "  %s const std::vector<const cppde::PchipForcing<%s>*>& F_)"
        % (" " * len(name), T),
        "    : params(p_), F(F_) {}",
        "",
    ]


def cpp_adjoint_terms(model, scalar, emit_jvp):
    """Lines of `struct adjoint_terms` (see OdeModel.contraction_stores)."""
    T = scalar.name
    n = model.n
    lines = ["// Contractions of the reverse step."]
    lines += _struct_head("adjoint_terms", T)
    vec = "const std::vector<%s>&" % T
    for name, kind, stores in model.contraction_stores(emit_jvp):
        with_v = name.startswith("jvp_")
        args = ["%s x" % vec] + (["%s v" % vec] if with_v else [])
        args += ["%s lam" % vec, "const %s& t" % T]
        if kind == "x":
            args.append("std::vector<%s>& out" % T)
            head = "  void %s(" % name
        elif kind == "p":
            args += ["const %s& sc" % T, "%s* out" % T]
            head = "  void %s(" % name
        else:
            head = "  %s %s(" % (T, name)
        pad = " " * len(head)
        sig = head + (",\n" + pad).join(args) + ") const {"
        lines += sig.split("\n")
        used = "    (void)x; (void)t; (void)lam;" + (" (void)v;" if with_v else "")
        if kind == "p":
            used += " (void)sc;"
        lines.append(used)
        if kind == "x":
            lines.append("    out.assign(%du, %s(0.0));" % (n, T))
        if kind != "dot":
            lines += scalar.arena_scope()
        body = stores_prelude(model, stores, T)
        body.append(block(model.g, stores, "_%s" % _short(name)))
        lines += cpp_statements_raw(model, body, scalar)
        lines += ["  }", ""]
    lines += ["};", "// adjoint_terms writes into %d slots" % (n + len(model.params))]
    return lines


def cpp_ode_system(model, scalar):
    """Lines of `struct ode_system`."""
    T = scalar.name
    lines = ["// ODE system"] + _struct_head("ode_system", T)
    lines += ["  void operator()(const std::vector<%s>& x," % T,
              "                  std::vector<%s>& dxdt," % T,
              "                  const %s& t) {" % T,
              "    (void)x; (void)t;"]
    lines += scalar.arena_scope()
    lines += cpp_statements_raw(model, ode_statements(model, T), scalar)
    lines += ["  }", "};"]
    return lines


def cpp_jacobian(model, scalar, sparse, strategy):
    """Lines of `struct jacobian` writing -J (see jacobian_statements)."""
    T = scalar.name
    mat = ("cppde::csc_matrix<%s>& W" % T) if sparse else ("cppde::dense_matrix<%s>& J" % T)
    nnz = sum(len(r) for r in model.pattern())
    lines = ["// Jacobian (%s, %s, %d nonzeros), writes -J"
             % ("sparse" if sparse else "dense", strategy, nnz)]
    lines += _struct_head("jacobian", T)
    lines += ["  void operator()(const std::vector<%s>& x," % T,
              "                  %s," % mat,
              "                  const %s& t," % T,
              "                  std::vector<%s>& dfdt) {" % T,
              "    (void)x; (void)t;"]
    scoped = scalar.arena and scalar.level == 1
    lines += cpp_statements_raw(
        model, jacobian_statements(model, sparse, strategy, T, scoped=scoped), scalar)
    lines += ["  }", "};"]
    return lines


def cpp_noop_jacobian(scalar):
    """Jacobian stub of an explicit method."""
    T = scalar.name
    lines = ["// Jacobian stub (explicit method)"] + _struct_head("jacobian", T)
    lines += ["  void operator()(const std::vector<%s>& x," % T,
              "                  cppde::dense_matrix<%s>& J," % T,
              "                  const %s& t," % T,
              "                  std::vector<%s>& dfdt) {" % T,
              "    (void)x; (void)J; (void)t; (void)dfdt;",
              "  }", "};"]
    return lines


def cpp_statements_raw(model, stmts, scalar, indent="    "):
    pr = em.Printer(model.g, model.slot("cpp"), style=scalar.style,
                    ad_level=scalar.level)
    return em.render_cpp(stmts, pr, scalar.name, indent=indent)


def py_statements_raw(model, stmts, indent="    "):
    pr = em.Printer(model.g, model.slot("py"), style="py")
    return em.render_py(stmts, pr, indent=indent)


def _short(name):
    return "".join(w[0] for w in name.split("_")) + "_"


# ---------------------------------------------------------------------------
# Events and root functions
# ---------------------------------------------------------------------------

_METHODS = {"replace": "EventMethod::Replace", "add": "EventMethod::Add",
            "multiply": "EventMethod::Multiply"}


def valid_value(value):
    """False for None, NaN, booleans and the strings NA/NaN/None/true/false."""
    if value is None or isinstance(value, bool):
        return False
    if isinstance(value, (int, float)):
        return value == value
    s = str(value).strip().lower()
    return s not in ("", "none", "nan", "na", "true", "false")


def event_rows(events_df):
    """Rows of an events table as dicts; rows without `var` are dropped."""
    if events_df is None:
        return []
    d = events_df.to_dict("list") if hasattr(events_df, "to_dict") else dict(events_df)
    if not d:
        return []
    lens = [len(v) for v in d.values() if isinstance(v, (list, tuple))]
    n = max(lens) if lens else 1

    def get(key, i):
        v = d.get(key)
        if isinstance(v, (list, tuple)):
            return v[i] if i < len(v) else None
        return v

    rows = []
    for i in range(n):
        if get("var", i) is None:
            continue
        rows.append({k: get(k, i) for k in
                     ("var", "value", "time", "root", "method", "terminal",
                      "direction")} | {"index": i})
    return rows


def _grad(model, node, mask=cg.F_STATE | cg.F_PARAM | cg.F_INIT):
    """(state, param, init) gradient dicts of a scalar node."""
    g = model.g
    xs, ps, inits = {}, {}, {}
    for leaf, a in model.ad.vjp([node], [g.ONE], mask).items():
        k, idx = g.attr[leaf][0], g.attr[leaf][1]
        (xs if k == cg.STATE else ps if k == cg.PARAM else inits)[idx] = a
    return xs, ps, inits


class EventCode:
    """C++ of the events of one model in one scalar (see generate_event_code).

    Args:
        model: OdeModel (its right-hand side feeds g_dot and G_tt).
        scalar: Scalar.
        has_rhs: False emits no G_tt.
    """

    def __init__(self, model, scalar, has_rhs=True):
        self.m = model
        self.s = scalar
        self.has_rhs = has_rhs
        self.T = scalar.name
        self.V = "std::vector<%s>" % self.T
        self.lslot = model.slot("cpp", params="full_params")

    def _printer(self, slot=None, style=None, names=None):
        return em.Printer(self.m.g, slot or self.lslot, style=style or self.s.style,
                          names=names, ad_level=self.s.level)

    def body(self, stores, indent, prefix="_e", slot=None, style=None, scalar=None):
        stmts, names = em.schedule(self.m.g, stores, prefix=prefix)
        pr = self._printer(slot, style, names)
        return em.render_cpp(stmts, pr, scalar or self.T, indent=indent)

    def inline(self, node, slot=None, style=None):
        return self._printer(slot, style).expr(node)

    def value_lambda(self, node, comment, indent="    "):
        head = "%s[full_params, &F](const %s& x, const %s& t) -> %s {" % (
            indent, self.V, self.T, self.T)
        lines = [head, indent + "  (void)x; (void)t;"]
        lines += self.body([(("ret",), node, "=")], indent + "  ")
        lines.append("%s},  // %s" % (indent, comment))
        return lines

    # -- stores shared with the Python backend ----------------------------

    def dg_dx_stores(self, g_node):
        """out[j] = dg/dx_j for every state (zeros included)."""
        g = self.m.g
        xs, _, _ = _grad(self.m, g_node, cg.F_STATE)
        return [(("vec", "out", j), xs.get(j, g.ZERO), "=") for j in range(self.m.n)]

    def dg_dt_node(self, g_node):
        return self.m.ad.time_derivative([g_node])[0]

    def gtt_node(self, g_node):
        """Second total time derivative of g along the right-hand side."""
        m = self.m
        gdot = m.ad.lie([g_node], m.rhs_plain, forcings=False)[0]
        return m.ad.lie([gdot], m.rhs_plain, forcings=False)[0]

    def gdot_node(self, g_node):
        return self.m.ad.lie([g_node], self.m.rhs_plain, forcings=True)[0]

    def case_x_stores(self, node):
        if node is None:
            return []
        xs, _, _ = _grad(self.m, node, cg.F_STATE)
        return [(("vec", "out", j), xs[j], "=") for j in sorted(xs)]

    def case_p_stores(self, node):
        if node is None:
            return []
        m = self.m
        g = m.g
        sc = m.sc()
        _, ps, inits = _grad(m, node, cg.F_PARAM | cg.F_INIT)
        stores = [(("vec", "out", m.n + k), g.mul(sc, ps[k]), "+=") for k in sorted(ps)]
        return stores + [(("vec", "out", i), g.mul(sc, inits[i]), "+=")
                         for i in sorted(inits)]

    def case_t_stores(self, node):
        if node is None:
            return []
        m = self.m
        d = m.ad.time_derivative([node])[0]
        if m.g.is_zero(d):
            return []
        return [(("vec", "out", 0), m.g.mul(m.sc(), d), "+=")]

    # -- C++ ---------------------------------------------------------------

    def dg_dx(self, g_node, i, indent="    "):
        lines = ["%s// dg/dx for root event %d" % (indent, i),
                 "%s[full_params, &F](const %s& x, const %s& t, %s& out) {"
                 % (indent, self.V, self.T, self.V),
                 indent + "  (void)x; (void)t;"]
        lines += self.body(self.dg_dx_stores(g_node), indent + "  ")
        lines.append(indent + "},  // dg_dx")
        return lines

    def dg_dt(self, g_node, i, indent="    "):
        lines = ["%s// dg/dt for root event %d" % (indent, i)]
        return lines + self.value_lambda(self.dg_dt_node(g_node), "dg_dt", indent)

    def g_dot_dot(self, g_node, i, indent="    "):
        m = self.m
        g = m.g
        if m.forcings or not self.has_rhs:
            return [indent + "nullptr  // g_dot_dot (FD fallback)"]
        gtt = self.gtt_node(g_node)
        lvl = self.s.level
        peel = (lambda e: e) if lvl == 0 else (
            (lambda e: "(%s).val()" % e) if lvl == 1 else (lambda e: "(%s).val().val()" % e))
        local = {}
        decl = []
        for n in g.topo([gtt]):
            if g.op[n] != cg.LEAF:
                continue
            at = g.attr[n]
            k = at[0]
            if k == cg.STATE:
                name, src = "_gx%d" % at[1], "x[%d]" % at[1]
            elif k == cg.PARAM:
                name, src = "_gp%d" % at[1], "full_params[%d]" % (m.n + at[1])
            elif k == cg.INIT:
                name, src = "_gi%d" % at[1], "full_params[%d]" % at[1]
            elif k == cg.TIME:
                name, src = "_gt", "t"
            else:
                raise em.EmitError("G_tt reads an unsupported leaf")
            local[at] = name
            decl.append((name, src))
        lines = ["%s// G_tt for root event %d" % (indent, i),
                 "%s[full_params, &F](const %s& x, const %s& t) -> double {"
                 % (indent, self.V, self.T),
                 indent + "  (void)x; (void)t;"]
        for name, src in sorted(decl):
            lines.append("%s  double %s = %s;" % (indent, name, peel(src)))
        lines += self.body([(("ret",), gtt, "=")], indent + "  ", prefix="_g",
                           slot=lambda at: local[at], style="double",
                           scalar="double")
        lines.append(indent + "}  // g_dot_dot")
        return lines

    # -- forward event code ----------------------------------------------

    def forward(self, rows):
        m = self.m
        out = []
        for r in rows:
            i = r["index"]
            var = str(r["var"])
            if var not in m.states:
                raise ValueError("Event %d: unknown state variable '%s'" % (i, var))
            vidx = m.states.index(var)
            if not valid_value(r["value"]):
                raise ValueError("Event %d: 'value' is required but is NA/None" % i)
            h = m.parse(str(r["value"]), label="event %d value" % i)
            method = _METHODS.get(str(r["method"]).lower() if r["method"] is not None
                                  else "replace", "EventMethod::Replace")
            if valid_value(r["time"]):
                tnode = m.parse(str(r["time"]), label="event %d time" % i)
                out += ["  // Fixed event %d: %s at t = %s" % (i, var, _flat(r["time"])),
                        "  fixed_events.emplace_back(FixedEvent<%s, %s>{" % (self.V, self.T),
                        "    %s,  // time" % self.inline(tnode),
                        "    %d,    // state_index" % vidx]
                out += self.value_lambda(h, "value_func")
                out += ["    %s  // method" % method, "  });", ""]
            elif valid_value(r["root"]):
                gnode = m.parse(str(r["root"]), label="event %d root" % i)
                terminal = str(r["terminal"]).lower() == "true" if r["terminal"] else False
                try:
                    direction = int(r["direction"]) if r["direction"] is not None else 0
                except (TypeError, ValueError):
                    direction = 0
                out += ["  // Root event %d: %s when %s = 0" % (i, var, _flat(r["root"])),
                        "  root_events.push_back(RootEvent<%s, %s>{" % (self.V, self.T)]
                out += self.value_lambda(gnode, "func (root condition g)")
                out += ["    %d,  // state_index" % vidx]
                out += self.value_lambda(h, "value_func (h)")
                out += ["    %s,  // method" % method,
                        "    %s,     // terminal" % ("true" if terminal else "false"),
                        "    %d,    // direction" % direction]
                if terminal:
                    out += ["    nullptr,       // dg_dx (not needed for terminal)",
                            "    nullptr,       // dg_dt (not needed for terminal)",
                            "    nullptr        // g_dot_dot (not needed for terminal)"]
                else:
                    out += self.dg_dx(gnode, i)
                    out += self.dg_dt(gnode, i)
                    out += self.g_dot_dot(gnode, i)
                out += ["  });", ""]
            else:
                raise ValueError("Event %d: must specify either 'time' or 'root'" % i)
        return out

    # -- event_adjoint_terms ---------------------------------------------

    def _case(self, stores):
        return self.body(stores, "      ", slot=self.m.slot("cpp")) if stores else []

    def _case_x(self, node):
        return self._case(self.case_x_stores(node))

    def _case_p(self, node):
        return self._case(self.case_p_stores(node))

    def _case_t(self, node):
        return self._case(self.case_t_stores(node))

    def adjoint(self, rows):
        m = self.m
        T, V = self.T, self.V
        fx, fp, ftp, fht = [], [], [], []
        rx, rp, rgp, rgdx, rgdp, rht = [], [], [], [], [], []
        nf = nr = 0
        for r in rows:
            h = m.parse(str(r["value"])) if valid_value(r["value"]) else None
            if valid_value(r["time"]):
                tnode = m.parse(str(r["time"]))
                fx.append((nf, self._case_x(h)))
                fp.append((nf, self._case_p(h)))
                ftp.append((nf, self._case_p(tnode)))
                fht.append((nf, self._case_t(h)))
                nf += 1
            elif valid_value(r["root"]):
                gnode = m.parse(str(r["root"]))
                gdot = self.gdot_node(gnode)
                rx.append((nr, self._case_x(h)))
                rp.append((nr, self._case_p(h)))
                rgp.append((nr, self._case_p(gnode)))
                rgdx.append((nr, self._case_x(gdot)))
                rgdp.append((nr, self._case_p(gdot)))
                rht.append((nr, self._case_t(h)))
                nr += 1
        xargs = "int ev, const %s& x, const %s& t, %s& out" % (V, T, V)
        pargs = "int ev, const %s& x, const %s& t, const %s& sc, %s* out" % (V, T, T, T)
        xhead = ["    (void)x; (void)t;", "    out.assign(%du, %s(0.0));" % (m.n, T)]
        phead = ["    (void)x; (void)t;"]
        out = ["// Derivatives of the jump expressions, by event and kind."]
        out += _struct_head("event_adjoint_terms", T)
        out += _switch("fixed_dh_dx", fx, xargs, xhead)
        out += _switch("fixed_dh_dp_axpy", fp, pargs, phead)
        out += _switch("fixed_dtime_dp_axpy", ftp,
                       "int ev, const %s& sc, %s* out" % (T, T), [])
        out += _switch("root_dh_dx", rx, xargs, xhead)
        out += _switch("root_dh_dp_axpy", rp, pargs, phead)
        out += _switch("root_dg_dp_axpy", rgp, pargs, phead)
        out += _switch("root_gdot_dx", rgdx, xargs, xhead)
        out += _switch("root_gdot_dp_axpy", rgdp, pargs, phead)
        out += _switch("fixed_dh_dt_axpy", fht, pargs, phead)
        out += _switch("root_dh_dt_axpy", rht, pargs, phead)
        out += ["};", "// event_adjoint_terms writes into %d slots"
                % (m.n + len(m.params))]
        return out


def _flat(v):
    return " ".join(str(v).split())


def _switch(name, cases, args, head):
    lines = ["  void %s(%s) const {" % (name, args)] + head
    if not any(body for _, body in cases):
        return lines + ["    (void)ev;", "  }", ""]
    lines.append("    switch (ev) {")
    for idx, body in cases:
        if not body:
            continue
        lines.append("    case %d: {" % idx)
        lines += body
        lines.append("      break; }")
    return lines + ["    default: break;", "    }", "  }", ""]


def empty_event_adjoint_terms(n_states, n_params, T):
    """event_adjoint_terms of a model without events."""
    V = "std::vector<%s>" % T
    xargs = "int ev, const %s& x, const %s& t, %s& out" % (V, T, V)
    pargs = "int ev, const %s& x, const %s& t, const %s& sc, %s* out" % (V, T, T, T)
    xhead = ["    (void)x; (void)t; (void)ev;", "    out.assign(%du, %s(0.0));" % (n_states, T)]
    phead = ["    (void)x; (void)t; (void)ev; (void)sc; (void)out;"]
    out = ["// Derivatives of the jump expressions, by event and kind."]
    out += _struct_head("event_adjoint_terms", T)
    for name, args, head in (("fixed_dh_dx", xargs, xhead),
                             ("fixed_dh_dp_axpy", pargs, phead),
                             ("root_dh_dx", xargs, xhead),
                             ("root_dh_dp_axpy", pargs, phead),
                             ("root_dg_dp_axpy", pargs, phead),
                             ("root_gdot_dx", xargs, xhead),
                             ("root_gdot_dp_axpy", pargs, phead),
                             ("fixed_dh_dt_axpy", pargs, phead),
                             ("root_dh_dt_axpy", pargs, phead)):
        out += ["  void %s(%s) const {" % (name, args)] + head + ["  }", ""]
    out += ["  void fixed_dtime_dp_axpy(int ev, const %s& sc, %s* out) const {" % (T, T),
            "    (void)ev; (void)sc; (void)out;", "  }", ""]
    return out + ["};"]


def rootfunc_code(model, rootfunc, scalar):
    """Terminal root functions, or the steady-state check for 'equilibrate'."""
    T = scalar.name
    V = "std::vector<%s>" % T
    if rootfunc is None:
        return []
    if isinstance(rootfunc, str):
        if rootfunc.strip().lower() == "equilibrate":
            return ["", "  // --- Steady-state termination (rootfunc = 'equilibrate') ---",
                    "  auto ss_termination = make_steady_state_termination"
                    "<ode_system, %s, %s>(sys, root_tol);" % (V, T), ""]
        rootfunc = [rootfunc]
    ec = EventCode(model, scalar)
    lines = ["", "  // --- User-defined root function termination ---"]
    for i, text in enumerate(rootfunc):
        text = str(text).strip()
        if not text:
            continue
        node = model.parse(text, label="rootfunc %d" % i)
        lines += ["  // rootfunc[%d]: %s" % (i, _flat(text)),
                  "  root_events.push_back(RootEvent<%s, %s>{" % (V, T)]
        lines += ec.value_lambda(node, "func")
        lines += ["    0,            // state_index (ignored for terminal)",
                  "    [](const %s&, const %s&) { return %s(0.0); },  // value_func"
                  % (V, T, T),
                  "    EventMethod::Replace,  // method",
                  "    true,         // terminal = true",
                  "    0,            // direction = 0 (any crossing)",
                  "    nullptr,      // dg_dx (not needed for terminal)",
                  "    nullptr,      // dg_dt (not needed for terminal)",
                  "    nullptr       // g_dot_dot (not needed for terminal)",
                  "  });"]
    return lines + [""]


def fixed_event_times(model, rows):
    """Double expressions of the fixed-event times over the flat parameter
    vector `params`, or None (a root event, or a time reading a state, the
    clock or a forcing)."""
    out = []
    pr = em.Printer(model.g, model.slot("cpp"), style="double")
    for r in rows:
        if valid_value(r["root"]):
            return None
        if not valid_value(r["time"]):
            continue
        node = model.parse(str(r["time"]))
        if model.g.flags[node] & (cg.F_FORCING | cg.F_STATE | cg.F_TIME):
            return None
        out.append(pr.expr(node))
    return out or None


_MODELS = {}


def model_for(rhs, params, forcings):
    """OdeModel of (rhs, params, forcings), shared by consecutive calls."""
    linear = cs.linear_min_terms()
    key = (tuple((str(k), str(v)) for k, v in rhs.items()),
           tuple(as_list(params)), tuple(as_list(forcings)), linear)
    hit = _MODELS.get(key)
    if hit is None:
        if len(_MODELS) >= 2:
            _MODELS.pop(next(iter(_MODELS)))
        hit = _MODELS[key] = OdeModel(rhs, params, forcings, linear)
    return hit
