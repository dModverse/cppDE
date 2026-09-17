"""Structure of large right-hand sides.

E: long sums sum_j c_j * x_j with numeric c_j become rows of a constant map
C; the sum is replaced by a LINROW leaf whose value is (C x)[r]. Generated
code evaluates C x once per call and assembles Jacobian rows as G C.

D: statements of equal structure become loops over index tables
(vector_block).
"""

import math
import os
from fractions import Fraction

import cppde_emit as em
import cppde_graph as cg

# Fewest state terms sharing one factor that form a map row.
LIN_MIN_TERMS = 12

# Coefficients per line of the generated tables.
_PER_LINE = 8

# Fewest coefficients stored through a value pool.
_POOL_MIN = 256


def linear_min_terms(value=None):
    """Row threshold: `value`, else CPPDE_LINEAR, else LIN_MIN_TERMS; 0 disables."""
    if value is None:
        env = os.environ.get("CPPDE_LINEAR", "").strip()
        try:
            value = int(env) if env else LIN_MIN_TERMS
        except ValueError:
            value = LIN_MIN_TERMS
    return max(0, int(value))


class LinMap:
    """Constant rows over the state vector.

    Attributes:
        n_cols: number of states.
        rows: tuples ((j, c), ...) with increasing j; equal rows are stored once.
    """

    def __init__(self, n_cols):
        self.n_cols = n_cols
        self.rows = []
        self._index = {}

    def __len__(self):
        return len(self.rows)

    def add(self, row):
        """Index of `row`, appended if new."""
        key = tuple(row)
        r = self._index.get(key)
        if r is None:
            r = self._index[key] = len(self.rows)
            self.rows.append(key)
        return r

    def nnz(self):
        return sum(len(r) for r in self.rows)

    def cols(self, r):
        return [j for j, _ in self.rows[r]]

    def dense(self):
        """True if every row has a coefficient for every state."""
        return all(len(r) == self.n_cols for r in self.rows)

    def cpp(self, name="linmap_"):
        """Namespace-scope tables and the map object `name`; frequently
        repeated coefficients are stored once, in a pool."""
        vals = [em.cpp_literal(c) for row in self.rows for _, c in row]
        out = ["// Constant linear map: %d rows over %d states, %d coefficients."
               % (len(self.rows), self.n_cols, len(vals))]
        pool = {}
        for v in vals:
            pool.setdefault(v, len(pool))
        if len(vals) >= _POOL_MIN and 4 * len(pool) <= len(vals):
            idt = "unsigned short" if len(pool) <= 0xFFFF else "int"
            out += _array("const double", name + "pool_", list(pool))
            out += _array("const %s" % idt, name + "vid_", [str(pool[v]) for v in vals])
            out.append("const std::vector<double> %sval_ = cppde::linmap_values("
                       "%spool_, %svid_, %d);" % (name, name, name, len(vals)))
            val = name + "val_.data()"
        else:
            out += _array("const double", name + "val_", vals)
            val = name + "val_"
        if self.dense():
            out.append("const cppde::linmap_dense %s{%d, %d, %s};"
                       % (name, len(self.rows), self.n_cols, val))
            return out
        ptr = [0]
        for row in self.rows:
            ptr.append(ptr[-1] + len(row))
        idx = [str(j) for row in self.rows for j, _ in row]
        out += _array("const int", name + "ptr_", [str(p) for p in ptr])
        out += _array("const int", name + "idx_", idx)
        out.append("const cppde::linmap_csr %s{%d, %d, %sptr_, %sidx_, %s};"
                   % (name, len(self.rows), self.n_cols, name, name, val))
        return out


def _array(decl, name, items):
    if not items:
        return ["%s %s[1] = {0};" % (decl, name)]
    lines = ["%s %s[] = {" % (decl, name)]
    for k in range(0, len(items), _PER_LINE):
        lines.append("  " + ", ".join(items[k:k + _PER_LINE]) + ",")
    lines.append("};")
    return lines


def _state_term(g, t):
    """(j, factor) if t is x_j times a state-free factor, else None."""
    o = g.op[t]
    if o == cg.LEAF:
        at = g.attr[t]
        return (at[1], g.ONE) if at[0] == cg.STATE else None
    if o != cg.MUL:
        return None
    j = None
    others = []
    for a in g.args[t]:
        if g.op[a] == cg.LEAF and g.attr[a][0] == cg.STATE and j is None:
            j = g.attr[a][1]
        elif g.flags[a] & cg.F_STATE:
            return None
        else:
            others.append(a)
    if j is None:
        return None
    return j, (others[0] if len(others) == 1 else g.prod(others))


def _finite(c):
    return not isinstance(c, float) or math.isfinite(c)


def linear_rows(g, roots, n_states, min_terms):
    """Long linear sums below `roots` as map rows.

    A sum qualifies with at least `min_terms` terms c * x_j * f of one
    state-free factor f; that group becomes f * LINROW(r).

    Returns:
        (new roots, LinMap, dict LINROW leaf -> the sum it stands for), or
        (roots, None, {}) if nothing qualifies.
    """
    if not min_terms:
        return list(roots), None, {}
    lm = LinMap(n_states)
    mapping = {}
    expand = {}
    for node in g.topo(roots):
        if (g.op[node] != cg.ADD or len(g.args[node]) < min_terms
                or not g.flags[node] & cg.F_STATE):
            continue
        c0, cs = g.attr[node]
        groups = {}
        order = []
        rest = []
        for c, t in zip(cs, g.args[node]):
            hit = _state_term(g, t) if _finite(c) else None
            if hit is None:
                rest.append((c, t))
                continue
            j, f = hit
            if f not in groups:
                groups[f] = []
                order.append(f)
            groups[f].append((j, c, t))
        rows = []
        for f in order:
            items = groups[f]
            if len(items) < min_terms:
                rest += [(c, t) for _, c, t in items]
                continue
            items.sort(key=lambda it: it[0])
            r = lm.add((j, c) for j, c, _ in items)
            leaf = g.leaf(cg.LINROW, r)
            if leaf not in expand:
                expand[leaf] = g.add_terms([(c, g.state(j)) for j, c, _ in items])
            rows.append((1, g.mul(f, leaf)))
        if rows:
            mapping[node] = g.add_terms(rest + rows, c0)
    if not mapping:
        return list(roots), None, {}
    return g.substitute(roots, mapping), lm, expand


# ===========================================================================
# D: loops over statements of equal structure
# ===========================================================================

# Fewest instances of a class, or terms of a sum, emitted as a loop.
VEC_MIN_INSTANCES = 16

_ATOMIC = frozenset((cg.NUM, cg.NAMED, cg.BOOL, cg.LEAF))
_MULTI = ("m",)


def vector_min_instances(value=None):
    """Loop threshold: `value`, else CPPDE_VECTORISE, else VEC_MIN_INSTANCES;
    0 disables."""
    if value is None:
        env = os.environ.get("CPPDE_VECTORISE", "").strip()
        try:
            value = int(env) if env else VEC_MIN_INSTANCES
        except ValueError:
            value = VEC_MIN_INSTANCES
    return max(0, int(value))


def _target_split(t):
    """(shape, ints) of a store target; the ints vary between instances."""
    k = t[0]
    if k == "vec" and type(t[2]) is int:
        return ("vec", t[1]), (t[2],)
    if k == "mat" and type(t[2]) is int and type(t[3]) is int:
        return ("mat", t[1]), (t[2], t[3])
    if k == "call":
        args = tuple(t[2])
        if args and all(type(a) is int for a in args):
            return ("call", t[1]["cpp"], t[1]["py"], None), args
        return ("call", t[1]["cpp"], t[1]["py"], args), ()
    return t, ()


def _target_join(shape, idx):
    k = shape[0]
    if k == "call":
        args = tuple(idx) if shape[3] is None else shape[3]
        return ("call", {"cpp": shape[1], "py": shape[2]}, args)
    if not idx:
        return shape
    if k == "vec":
        return ("vec", shape[1], idx[0])
    return ("mat", shape[1], idx[0], idx[1])


class _Cone:
    """Key of one cone: its structure with leaves, value points and known
    nodes as slots numbered by first occurrence. A sum with at least `lim`
    terms of one structure becomes a reduction: those terms are keyed once,
    in a cone of their own each."""

    def __init__(self, an, reduce=True):
        self.an = an
        self.reduce = reduce
        self.slots = []
        self.slot_of = {}
        self.local = {}
        self.reds = []

    def ref(self, n, kind):
        s = self.slot_of.get(n)
        if s is None:
            s = self.slot_of[n] = len(self.slots)
            self.slots.append(n)
        return ("s", s, kind)

    def visit(self, n, top=False):
        an = self.an
        g = an.g
        o = g.op[n]
        if o == cg.NUM or o == cg.NAMED or o == cg.BOOL:
            return ("c", o, g.attr[n])
        if n in an.known:
            return self.ref(n, ("k",))
        if o == cg.LEAF:
            at = g.attr[n]
            k = at[0]
            if len(at) < 2 or k == cg.LOOPVAR or (k == cg.VEC and at[2] is None):
                return ("c", o, at)
            return self.ref(n, ("L", k, at[1] if k == cg.VEC else None))
        if n in an.value and not top:
            return self.ref(n, ("P", an.keyid[n]))
        hit = self.local.get(n)
        if hit is not None:
            return ("r", hit)
        self.local[n] = len(self.local)
        args = g.args[n]
        if self.reduce and o == cg.ADD and len(args) >= an.lim:
            tok = self._sum(n)
            if tok is not None:
                return tok
        return ("n", o, g.attr[n], tuple(self.visit(a) for a in args))

    def _sum(self, n):
        g = self.an.g
        c0, cs = g.attr[n]
        groups = {}
        order = []
        for c, t in zip(cs, g.args[n]):
            sub = _Cone(self.an, reduce=False)
            key = (sub.visit(t), c)
            if key not in groups:
                groups[key] = []
                order.append(key)
            groups[key].append((t, sub.slots))
        red = [k for k in order if len(groups[k]) >= self.an.lim]
        if not red:
            return None
        inner = set()
        for k in red:
            inner.update(t for t, _ in groups[k])
            self.reds.append((n, k[1], groups[k]))
        rest = tuple((c, self.visit(t)) for c, t in zip(cs, g.args[n]) if t not in inner)
        return ("A", c0, rest, tuple(k for k in red))


class _Analysis:
    """Value points, keys and classes of a list of stores.

    A value point is a non-boolean node used from more than one cone; every
    other node belongs to the cone of the one store or value point it feeds.
    An instance is a store ('s', index) or a value point ('v', node).
    """

    def __init__(self, g, stores, known, lim):
        self.g = g
        self.stores = stores
        self.known = known
        self.lim = lim
        roots = [n for _, n, _ in stores]
        self.order = [n for n in g.topo(roots) if n not in known]
        self.value = self._values(roots)
        self.keyid = {}
        self.slots = {}
        self.reds = {}
        self.classes = {}
        intern = {}
        for n in self.order:
            if n in self.value:
                self._add(("v", n), None, None, intern)
        for sid, (target, n, op) in enumerate(stores):
            self._add(("s", sid), n, (_target_split(target)[0], op), intern)

    def _values(self, roots):
        g = self.g
        users = {}
        for sid, n in enumerate(roots):
            users.setdefault(n, []).append(("s", sid))
        for m in self.order:
            for a in g.args[m]:
                users.setdefault(a, []).append(("n", m))
        owner = {}
        value = set()
        for n in reversed(self.order):
            if g.op[n] in _ATOMIC:
                continue
            own = None
            for kind, u in users.get(n, ()):
                if kind == "s":
                    o = ("s", u)
                elif u in value:
                    o = ("v", u)
                else:
                    o = owner[u]
                if own is None:
                    own = o
                elif o != own:
                    own = _MULTI
                if own == _MULTI:
                    break
            if own != _MULTI:
                owner[n] = own
            elif g.is_bool(n):
                owner[n] = _MULTI
            else:
                value.add(n)
                owner[n] = ("v", n)
        return value

    def _add(self, inst, root, head, intern):
        cone = _Cone(self)
        try:
            if inst[0] == "v":
                tok = cone.visit(inst[1], top=True)
            else:
                tok = head + (cone.visit(root),)
        except RecursionError:
            cone = _Cone(self)
            tok = ("u", inst)
        kid = intern.setdefault(tok, len(intern))
        if inst[0] == "v":
            self.keyid[inst[1]] = kid
        self.slots[inst] = cone.slots
        self.reds[inst] = cone.reds
        self.classes.setdefault(kid, []).append(inst)

    def is_value_class(self, kid):
        return self.classes[kid][0][0] == "v"


def _progression(values):
    """(a, b) with values[q] == a*q + b for all q, else None."""
    b = values[0]
    a = values[1] - b if len(values) > 1 else 0
    for q, v in enumerate(values):
        if v != a * q + b:
            return None
    return a, b


class _Tables:
    """Index tables of one block, shared by equal content."""

    def __init__(self, prefix):
        self.prefix = prefix
        self.names = {}

    def text(self, values, var):
        pr = _progression(values)
        if pr is not None:
            a, b = pr
            if a == 0:
                return str(b)
            t = var if abs(a) == 1 else "%d*%s" % (abs(a), var)
            if a < 0:
                return "(%d - %s)" % (b, t)
            return "(%s %s %d)" % (t, "+" if b > 0 else "-", abs(b)) if b else t
        key = tuple(values)
        name = self.names.get(key)
        if name is None:
            name = self.names[key] = "%si%d" % (self.prefix, len(self.names))
        return "%s[%s]" % (name, var)

    def statements(self):
        return [("table", name, list(vals)) for vals, name in self.names.items()]


def vector_block(g, stores, prefix, known=None, min_instances=None):
    """Statements computing `stores` with loops, or None if nothing loops.

    Statements whose cones have equal structure form a class; a class of at
    least `min_instances` (default: vector_min_instances()) runs as one loop
    that reads the differing leaves and shared values through index tables.
    A sum with that many terms of one structure runs as an inner loop.
    Shared values of a class live in the array `<prefix>v<class>`.
    """
    lim = vector_min_instances(min_instances)
    if not lim or not stores:
        return None
    known = dict(known or {})
    an = _Analysis(g, stores, known, lim)
    loops = {k for k, inst in an.classes.items()
             if len(inst) >= lim or an.reds[inst[0]]}
    if not loops:
        return None
    order = _unit_order(an, loops)
    if order is None:
        return None

    rank = {n: q for q, n in enumerate(an.order)}
    for kid, inst in an.classes.items():
        if an.is_value_class(kid) and not _self_reading(an, kid):
            inst.sort(key=lambda i: _slot_indices(an, i, rank))

    names = dict(known)
    arrays = {}
    pos = {}
    for kid, inst in an.classes.items():
        if not an.is_value_class(kid):
            continue
        if len(inst) == 1:
            names[inst[0][1]] = "%sv%d" % (prefix, kid)
            continue
        arrays[kid] = "%sv%d" % (prefix, kid)
        for q, (_, n) in enumerate(inst):
            pos[n] = q
            names[n] = "%s[%d]" % (arrays[kid], q)

    tables = _Tables(prefix)
    emit = _Emitter(an, names, arrays, pos, tables, prefix)
    body = []
    batch = []
    for kind, item in order:
        if kind == "loop":
            stmts = emit.loop(item)
            if stmts is not None:
                body += emit.batch(batch)
                batch = []
                body += stmts
                continue
            if any(an.reds[i] for i in an.classes[item]):
                return None
            batch += an.classes[item]
        else:
            batch.append(item)
    body += emit.batch(batch)
    return tables.statements() + body


def _self_reading(an, kid):
    return any(s in an.value and an.keyid[s] == kid
               for i in an.classes[kid] for s in an.slots[i])


def _slot_indices(an, inst, rank):
    g = an.g
    out = []
    for s in an.slots[inst]:
        if g.op[s] == cg.LEAF:
            at = g.attr[s]
            out.append(at[2] if at[0] == cg.VEC else at[1])
        else:
            out.append(rank.get(s, -1))
    return tuple(out)


def _unit_order(an, loops):
    """[('loop', class) or ('one', instance)] in a valid order, or None.

    Value classes follow their dependencies. Store classes follow in the
    order of their first store; for each target an '=' write must follow
    every earlier write and precede the later '+=' writes.
    """
    vkids = [k for k in an.classes if an.is_value_class(k)]
    deps = {k: set() for k in vkids}
    users = {k: [] for k in vkids}
    for k in vkids:
        inst = an.classes[k]
        rank = {i[1]: q for q, i in enumerate(inst)}
        for q, i in enumerate(inst):
            refs = list(an.slots[i])
            for _, _, terms in an.reds[i]:
                for _, ts in terms:
                    refs += ts
            for s in refs:
                if s not in an.value:
                    continue
                d = an.keyid[s]
                if d == k:
                    if k in loops and rank[s] >= q:
                        return None
                    continue
                if d not in deps[k]:
                    deps[k].add(d)
                    users[d].append(k)
    indeg = {k: len(deps[k]) for k in vkids}
    ready = [k for k in vkids if not indeg[k]]
    seq = []
    while ready:
        k = ready.pop()
        seq.append(k)
        for u in users[k]:
            indeg[u] -= 1
            if not indeg[u]:
                ready.append(u)
    if len(seq) != len(vkids):
        return None
    seq.sort(key=lambda k: an.order.index(an.classes[k][0][1])
             if len(vkids) < 64 else 0)
    skids = sorted((k for k in an.classes if not an.is_value_class(k)),
                   key=lambda k: an.classes[k][0][1])
    rank = {}
    for r, k in enumerate(skids):
        for i in an.classes[k]:
            rank[i[1]] = r
    top = {}
    last_set = {}
    for sid, (target, _, op) in enumerate(an.stores):
        key = _target_split(target)
        r = rank[sid]
        if op == "=":
            if top.get(key, -1) > r:
                return None
            last_set[key] = r
        elif last_set.get(key, -1) > r:
            return None
        top[key] = max(top.get(key, -1), r)
    order = []
    for k in _dep_sorted(seq, deps) + skids:
        if k in loops:
            order.append(("loop", k))
        else:
            order += [("one", i) for i in an.classes[k]]
    return order


def _dep_sorted(seq, deps):
    """`seq` reordered so that every class follows its dependencies, keeping
    the given order where possible."""
    done = set()
    out = []
    for k in seq:
        stack = [k]
        while stack:
            c = stack[-1]
            if c in done:
                stack.pop()
                continue
            todo = [d for d in deps[c] if d not in done]
            if todo:
                stack += todo
                continue
            done.add(c)
            out.append(c)
            stack.pop()
    return out


class _Emitter:
    def __init__(self, an, names, arrays, pos, tables, prefix):
        self.an = an
        self.g = an.g
        self.names = names
        self.arrays = arrays
        self.pos = pos
        self.tables = tables
        self.prefix = prefix
        self.counter = 0

    def batch(self, items):
        """Scalar statements of `items`, scheduled together."""
        if not items:
            return []
        an = self.an
        own = {i[1] for i in items if i[0] == "v"}
        known = {k: v for k, v in self.names.items() if k not in own}
        stores = []
        head = []
        scalars = {}
        for i in items:
            if i[0] == "s":
                stores.append(an.stores[i[1]])
                continue
            n = i[1]
            kid = an.keyid[n]
            arr = self.arrays.get(kid)
            if arr is None:
                scalars[n] = self.names[n]
                stores.append((("var", self.names[n]), n, "="))
                continue
            if self.pos[n] == 0:
                head.append(("array", arr, len(an.classes[kid])))
            stores.append((("vec", arr, self.pos[n]), n, "="))
        p = "%ss%d_" % (self.prefix, self.counter)
        self.counter += 1
        stmts, nm = em.schedule(self.g, stores, prefix=p, known=known)
        out = []
        for st in stmts:
            if st[0] == "decl" and st[2] in scalars:
                nm[st[2]] = scalars[st[2]]
                out.append(("decl", scalars[st[2]], st[2], st[3]))
            elif st[0] == "store" and st[1][0] == "var" and st[2] in scalars:
                if nm.get(st[2]) != scalars[st[2]]:
                    out.append(("decl", scalars[st[2]], st[2], False))
            else:
                out.append(st)
        return head + [("block", nm, out)]

    def _leaves(self, slot_rows, var):
        """Substitutes of the slots of the first row, or None."""
        g = self.g
        leaves = {}
        for k, s0 in enumerate(slot_rows[0]):
            refs = [row[k] for row in slot_rows]
            if all(r == s0 for r in refs):
                if s0 in self.names:
                    leaves[s0] = g.leaf(cg.LOOPVAR, "ref", self.names[s0])
                continue
            if g.op[s0] == cg.LEAF:
                at = g.attr[s0]
                if at[0] == cg.VEC:
                    idx = [g.attr[r][2] for r in refs]
                    base = (at[0], at[1])
                else:
                    idx = [g.attr[r][1] for r in refs]
                    base = (at[0],)
                text = self.tables.text(idx, var)
                leaves[s0] = g.leaf(cg.LOOPVAR, *(base + (text,)))
            elif s0 in self.pos:
                text = self.tables.text([self.pos[r] for r in refs], var)
                arr = self.arrays[self.an.keyid[s0]]
                leaves[s0] = g.leaf(cg.LOOPVAR, "ref", "%s[%s]" % (arr, text))
            else:
                return None
        return leaves

    def loop(self, kid):
        """A loop over the instances of class `kid`, or None if a slot
        cannot be read through a table."""
        g = self.g
        an = self.an
        inst = an.classes[kid]
        var = "_q"
        first = inst[0]
        leaves = self._leaves([an.slots[i] for i in inst], var)
        if leaves is None:
            return None
        pre = []
        sums = {}
        for rho, (add0, coef, terms0) in enumerate(an.reds[first]):
            rows = [an.reds[i][rho][2] for i in inst]
            counts = [len(r) for r in rows]
            ptr = [0]
            for c in counts:
                ptr.append(ptr[-1] + c)
            tl = self._leaves([ts for r in rows for _, ts in r], "_p")
            if tl is None:
                return None
            acc = "%sr%d_%d" % (self.prefix, kid, rho)
            term = raw_substitute(g, terms0[0][0], tl, top=False)
            if coef != 1:
                term = g._node(cg.ADD, (term,), (Fraction(0), (coef,)))
            st, nm = em.schedule(g, [(("var", acc), term, "+=")],
                                 prefix="%sp%d_%d_" % (self.prefix, kid, rho))
            pre += [("local", acc),
                    ("loop", "_p", self.tables.text(ptr[:-1], var),
                     self.tables.text(ptr[1:], var), [("block", nm, st)])]
            sums.setdefault(add0, []).append(
                ({t for t, _ in terms0}, g.leaf(cg.LOOPVAR, "ref", acc)))
        stmts = []
        if first[0] == "v":
            start = first[1]
            op = "="
            if kid in self.arrays:
                target = ("vec", self.arrays[kid], var)
                stmts.append(("array", self.arrays[kid], len(inst)))
            else:
                target = ("var", self.names[start])
                stmts.append(("local", self.names[start]))
            skel = raw_substitute(g, start, leaves, top=True, sums=sums)
        else:
            tgt, start, op = an.stores[first[1]]
            shape, _ = _target_split(tgt)
            cols = zip(*[_target_split(an.stores[i[1]][0])[1] for i in inst])
            target = _target_join(shape, tuple(self.tables.text(list(c), var)
                                               for c in cols))
            skel = raw_substitute(g, start, leaves, top=False, sums=sums)
        body, nm = em.schedule(g, [(target, skel, op)],
                               prefix="%sl%d_" % (self.prefix, kid))
        if len(inst) == 1:
            return stmts + pre + [("block", nm, body)]
        outer, inner = em.split_loop(g, body, cg.F_LOOPVAR)
        if outer:
            stmts.append(("block", nm, outer))
        stmts.append(("for", var, str(len(inst)), pre + [("block", nm, inner)]))
        return stmts


def raw_substitute(g, start, leaves, top=False, sums=None):
    """The cone of `start` with `leaves` substituted and argument order kept.

    `top` builds the start node even if it has a substitute; `sums` maps a
    sum node to [(terms, accumulator leaf)] that replace those terms.
    """
    memo = {}
    sums = sums or {}

    def build(n, first=False):
        if not first:
            hit = leaves.get(n)
            if hit is not None:
                return hit
        if g.op[n] in _ATOMIC:
            return n
        hit = memo.get(n)
        if hit is not None:
            return hit
        red = sums.get(n)
        if red is None:
            args = tuple(build(a) for a in g.args[n])
            r = n if args == g.args[n] else g._node(g.op[n], args, g.attr[n])
        else:
            inner = set().union(*(t for t, _ in red))
            c0, cs = g.attr[n]
            keep = [(c, build(a)) for c, a in zip(cs, g.args[n]) if a not in inner]
            keep += [(Fraction(1), acc) for _, acc in red]
            r = g._node(cg.ADD, tuple(a for _, a in keep),
                        (c0, tuple(c for c, _ in keep)))
        memo[n] = r
        return r

    return build(start, top)
