"""Scheduling and printing of graph nodes as C++ or Python statements.

`schedule` turns stores into statements. A node becomes a temporary if it is
used twice or its line exceeds MAX_INLINE; a once-used select stays inline.
Temporaries are ordered parameter-only, then time/forcing-dependent, then the
rest.

Statements are tuples:
  ('decl', name, node, is_bool)       a temporary
  ('store', target, node, op)         target op expr, op in '=', '+=', '-='
  ('comment', text)
  ('raw', {'cpp': text, 'py': text})  a backend's own line, either may be None
  ('if', {'cpp': cond, 'py': cond}, [statements])
  ('for', var, count, [statements])           var in [0, count)
  ('loop', var, lo, hi, [statements])         var in [lo, hi), bounds as text
  ('assign', target, value)                   value text, or {'cpp', 'py'}
  ('table', name, [int, ...])                 static const int table
  ('array', name, size)                       vector of the scalar type
  ('local', name)                             scalar variable, initially 0
  ('block', names, [statements])              statements printed with names
A target is ('vec', name, index), ('mat', name, i, j), ('var', name or
{'cpp', 'py'}), ('ret',) (return statement) or
('call', {'cpp': fmt, 'py': fmt}, args), the line fmt.format(*args) % expr.
Indices and args are ints or index text.
"""

import math
from fractions import Fraction

import cppde_graph as cg


class EmitError(ValueError):
    pass


# Largest inline tree size per line.
MAX_INLINE = 48

_ATOMIC = frozenset((cg.NUM, cg.NAMED, cg.BOOL, cg.LEAF))

# Functions whose C++ AD overloads take values only.
_SCALAR_ARGS = frozenset(("min", "max"))


# ---------------------------------------------------------------------------
# Scheduling
# ---------------------------------------------------------------------------

def _twice(g, n):
    """Children printed more than once in node n."""
    o = g.op[n]
    if o == cg.POW:
        b, e = g.args[n]
        if g.op[e] == cg.NUM and g.attr[e] == 2:
            return (b,)
    elif o == cg.CALL and g.attr[n] in ("Heaviside", "sign"):
        return (g.args[n][0],)
    return ()


def schedule(g, stores, prefix="_t", max_inline=MAX_INLINE, groups=True,
             known=None):
    """Statements for `stores`, a list of (target, node, op).

    Args:
        known: dict node -> name of values already in scope.

    Returns:
        (statements, dict node -> temporary name).
    """
    known = dict(known or {})
    roots = [n for _, n, _ in stores]
    order = [n for n in g.topo(roots) if n not in known]
    uses = {}
    for n in order:
        for c in g.args[n]:
            uses[c] = uses.get(c, 0) + 1
        for c in _twice(g, n):
            uses[c] = uses.get(c, 0) + 1
    for r in roots:
        uses[r] = uses.get(r, 0) + 1

    mat = {}
    cost = {}
    for n in order:
        o = g.op[n]
        if o in _ATOMIC:
            cost[n] = 1
            continue
        if uses.get(n, 0) >= 2 and not (o == cg.SELECT and uses[n] == 1):
            mat[n] = True
        if o == cg.CALL and g.attr[n] in _SCALAR_ARGS:
            # No expression-template overloads: arguments must be values.
            for a in g.args[n]:
                if g.op[a] not in _ATOMIC and a not in known:
                    mat[a] = True
                    cost[a] = 1
        c = 1
        kids = []
        for a in g.args[n]:
            if a in known or a in mat:
                c += 1
            else:
                c += cost[a]
                if g.op[a] not in _ATOMIC:
                    kids.append(a)
        if c > max_inline:
            kids.sort(key=lambda a: -cost[a])
            for a in kids:
                if c <= max_inline:
                    break
                if g.op[a] == cg.SELECT and uses.get(a, 0) == 1:
                    continue
                mat[a] = True
                c -= cost[a] - 1
        cost[n] = 1 if n in mat else c

    names = dict(known)
    decls = []
    counter = 0
    for n in order:
        if n not in mat:
            continue
        name = "%s%d" % (prefix, counter)
        counter += 1
        decls.append((_group(g, n) if groups else 0, len(decls), name, n))
    decls.sort()
    stmts = []
    for _, _, name, n in decls:
        stmts.append(("decl", name, n, g.is_bool(n)))
        names[n] = name
    for target, n, op in stores:
        stmts.append(("store", target, n, op))
    return stmts, names


_STATIC = cg.F_PARAM | cg.F_INIT
_CLOCK = _STATIC | cg.F_TIME | cg.F_FORCING | cg.F_FRATE


def _group(g, n):
    f = g.flags[n]
    if not (f & ~_STATIC):
        return 0
    if not (f & ~_CLOCK):
        return 1
    return 2


# ---------------------------------------------------------------------------
# Printing
# ---------------------------------------------------------------------------

# precedence: higher binds tighter
P_OR, P_AND, P_CMP, P_ADD, P_MUL, P_UNARY, P_ATOM = range(7)

# name -> (double spelling, AD spelling or None, template spelling)
CPP_CALLS = {
    "exp": ("std::exp", "cppde::exp", "exp"),
    "log": ("std::log", "cppde::log", "log"),
    "sin": ("std::sin", "cppde::sin", "sin"),
    "cos": ("std::cos", "cppde::cos", "cos"),
    "tan": ("std::tan", "cppde::tan", "tan"),
    "asin": ("std::asin", "cppde::asin", "asin"),
    "acos": ("std::acos", "cppde::acos", "acos"),
    "atan": ("std::atan", "cppde::atan", "atan"),
    "sinh": ("std::sinh", "cppde::sinh", "sinh"),
    "cosh": ("std::cosh", "cppde::cosh", "cosh"),
    "tanh": ("std::tanh", "cppde::tanh", "tanh"),
    "asinh": ("std::asinh", "cppde::asinh", "asinh"),
    "acosh": ("std::acosh", "cppde::acosh", "acosh"),
    "atanh": ("std::atanh", "cppde::atanh", "atanh"),
    "abs": ("std::fabs", "cppde::abs", "abs"),
    "min": ("std::min", "cppde::min", "min"),
    "max": ("std::max", "cppde::max", "max"),
    "erf": ("std::erf", None, None),
    "erfc": ("std::erfc", None, None),
    "gamma": ("std::tgamma", None, None),
    "loggamma": ("std::lgamma", None, None),
    "atan2": ("std::atan2", None, None),
    "beta": ("std::beta", None, None),
    "besselj": ("std::cyl_bessel_j", None, None),
    "bessely": ("std::cyl_neumann", None, None),
    "besseli": ("std::cyl_bessel_i", None, None),
    "besselk": ("std::cyl_bessel_k", None, None),
}

PY_CALLS = {name: "R." + name for name in (
    "exp", "log", "sin", "cos", "tan", "asin", "acos", "atan", "sinh", "cosh",
    "tanh", "asinh", "acosh", "atanh", "abs", "min", "max", "erf", "erfc",
    "gamma", "loggamma", "atan2", "floor", "ceiling", "sign", "Heaviside",
    "delta")}


def cpp_literal(v):
    if v == "nan":
        return "std::numeric_limits<double>::quiet_NaN()"
    if isinstance(v, Fraction):
        if v.denominator == 1 and abs(v.numerator) < (1 << 53):
            return "%d.0" % v.numerator
        v = float(v)
    if math.isinf(v):
        return ("std::numeric_limits<double>::infinity()" if v > 0
                else "-std::numeric_limits<double>::infinity()")
    if math.isnan(v):
        return "std::numeric_limits<double>::quiet_NaN()"
    s = repr(float(v))
    if "e" not in s and "." not in s:
        s += ".0"
    return s


def py_literal(v):
    if v == "nan":
        return "R.NAN"
    if isinstance(v, Fraction):
        if v.denominator == 1:
            return "%d.0" % v.numerator
        v = float(v)
    if math.isinf(v):
        return "R.INF" if v > 0 else "(-R.INF)"
    if math.isnan(v):
        return "R.NAN"
    return repr(float(v))


class Printer:
    """Expression printer.

    Args:
        slot: function leaf attr -> string.
        style: 'double', 'ad' or 'template' (C++), or 'py'.
        names: dict node -> temporary name.
        ad_level: derivative layers of the C++ scalar in template style.
    """

    def __init__(self, g, slot, style="double", names=None, ad_level=0):
        self.g = g
        self.slot = slot
        self.style = style
        self.py = style == "py"
        self.names = names if names is not None else {}
        self.ad_level = ad_level

    def with_names(self, names):
        """Copy printing the nodes in `names` as those names."""
        p = Printer(self.g, self.slot, self.style, names, self.ad_level)
        return p

    def lit(self, v):
        return py_literal(v) if self.py else cpp_literal(v)

    def expr(self, n, top=False):
        return self._p(n, top)[0]

    def _atom(self, n):
        name = self.names.get(n)
        if name is not None:
            return name, P_ATOM
        return None

    def _p(self, n, top=False):
        g = self.g
        if not top:
            hit = self._atom(n)
            if hit:
                return hit
        o = g.op[n]
        at = g.attr[n]
        if o == cg.NUM:
            s = self.lit(at)
            return s, (P_UNARY if s.startswith("-") else P_ATOM)
        if o == cg.NAMED:
            return self.lit(at[1]), P_ATOM
        if o == cg.BOOL:
            if self.py:
                return ("True" if at else "False"), P_ATOM
            return ("true" if at else "false"), P_ATOM
        if o == cg.LEAF:
            return self.slot(at), P_ATOM
        if o == cg.ADD:
            return self._add(n)
        if o == cg.MUL:
            return self._mul(g.args[n], Fraction(1))
        if o == cg.POW:
            return self._pow(n)
        if o == cg.CALL:
            return self._call(n), P_ATOM
        if o == cg.SELECT:
            c, a, b = g.args[n]
            f = "R.select" if self.py else "cppde::select"
            return "%s(%s, %s, %s)" % (f, self.expr(c), self.expr(a),
                                        self.expr(b)), P_ATOM
        if o == cg.CMP:
            a, b = g.args[n]
            return "%s %s %s" % (self._wrap(a, P_ADD), at,
                                 self._wrap(b, P_ADD)), P_CMP
        if o in (cg.AND, cg.OR):
            if self.py:
                sep = " and " if o == cg.AND else " or "
            else:
                sep = " && " if o == cg.AND else " || "
            prec = P_AND if o == cg.AND else P_OR
            return sep.join(self._wrap(a, prec + 1) for a in g.args[n]), prec
        if o == cg.NOT:
            a = g.args[n][0]
            if self.py:
                return "(not %s)" % self._wrap(a, P_ATOM), P_ATOM
            return "!%s" % self._wrap(a, P_ATOM), P_UNARY
        raise EmitError("cannot print " + cg.OP_NAMES[o])

    def _wrap(self, n, prec):
        s, p = self._p(n)
        return s if p >= prec else "(" + s + ")"

    def _add(self, n):
        g = self.g
        c0, cs = g.attr[n]
        parts = []
        for c, t in zip(cs, g.args[n]):
            neg = c < 0
            mag = -c if neg else c
            if g.op[t] == cg.MUL and not self._atom(t):
                body, _ = self._mul(g.args[t], mag)
            elif mag == 1:
                body = self._wrap(t, P_MUL)
            else:
                body = "%s*%s" % (self.lit(mag), self._wrap(t, P_MUL))
            parts.append((neg, body))
        if c0 != 0:
            parts.append((c0 < 0, self.lit(-c0 if c0 < 0 else c0)))
        s = ("-" if parts[0][0] else "") + parts[0][1]
        for neg, body in parts[1:]:
            s += (" - " if neg else " + ") + body
        if len(parts) == 1:
            return s, (P_UNARY if parts[0][0] else P_MUL)
        return s, P_ADD

    def _mul(self, factors, coef):
        g = self.g
        num, den = [], []
        for f in factors:
            if not self._atom(f) and g.op[f] == cg.POW:
                b, e = g.args[f]
                if g.op[e] == cg.NUM and g.attr[e] != "nan" and g.attr[e] < 0:
                    den.append(self._power(b, -g.attr[e]))
                    continue
            num.append(self._wrap(f, P_MUL))
        if coef != 1:
            num.insert(0, self.lit(coef))
        s = "*".join(num) if num else self.lit(1)
        if den:
            d = den[0] if len(den) == 1 else "(" + "*".join(den) + ")"
            if self.py:
                return "R.div(%s, %s)" % (s, d), P_ATOM
            s = "%s/%s" % (s, d)
        return s, P_MUL

    def _power(self, b, e):
        """b ** e for numeric e > 0, as a factor."""
        if e == 1:
            return self._wrap(b, P_ATOM)
        if e == 2:
            x = self._wrap(b, P_ATOM)
            return "(%s*%s)" % (x, x)
        if e == Fraction(1, 2):
            return "%s(%s)" % (self._fn("sqrt"), self.expr(b))
        return "%s(%s, %s)" % (self._fn("pow"), self.expr(b), self.lit(e))

    def _pow(self, n):
        g = self.g
        b, e = g.args[n]
        if g.op[e] == cg.NUM and g.attr[e] != "nan":
            v = g.attr[e]
            if v < 0:
                if self.py:
                    return "R.div(1.0, %s)" % self._power(b, -v), P_ATOM
                return "%s/%s" % (self.lit(1), self._power(b, -v)), P_MUL
            s = self._power(b, v)
            return s, P_ATOM
        return "%s(%s, %s)" % (self._fn("pow"), self.expr(b), self.expr(e)), P_ATOM

    def _fn(self, name):
        if self.py:
            return "R." + name
        if self.style == "double":
            return "std::" + name
        if self.style == "ad":
            return "cppde::" + name
        return name

    def _call(self, n):
        g = self.g
        name = g.attr[n]
        args = g.args[n]
        if self.py:
            f = PY_CALLS.get(name, "R.call_%s" % name)
            return "%s(%s)" % (f, ", ".join(self.expr(a) for a in args))
        if name == "Heaviside":
            x = self._wrap(args[0], P_ADD)
            return ("cppde::select(%s < 0.0, 0.0, cppde::select(%s == 0.0, 0.5, 1.0))"
                    % (x, x))
        if name == "sign":
            x = self._wrap(args[0], P_ADD)
            return ("cppde::select(%s > 0.0, 1.0, cppde::select(%s < 0.0, -1.0, 0.0))"
                    % (x, x))
        if name == "delta":
            x = self._wrap(args[0], P_ADD)
            return "cppde::select(%s == 0.0, 1.0, 0.0)" % x
        if name in ("floor", "ceiling"):
            fn = "std::floor" if name == "floor" else "std::ceil"
            if self.style == "double":
                return "%s(%s)" % (fn, self.expr(args[0]))
            # Value only; zero derivative.
            return "%s(cppde::value_of(%s))" % (fn, self.expr(args[0]))
        spell = CPP_CALLS.get(name)
        if spell is None:
            raise EmitError("function '%s' has no C++ equivalent" % name)
        if self.style == "double":
            f = spell[0]
        elif self.style == "ad":
            f = spell[1]
        else:
            f = spell[2]
        if f is None:
            if self.style == "template" and self.ad_level == 0:
                f = spell[0]
            else:
                raise EmitError(
                    "function '%s' has no overload for the AD scalar; it can "
                    "be used where the model is evaluated in plain double only"
                    % name)
        return "%s(%s)" % (f, ", ".join(self.expr(a) for a in args))


# ---------------------------------------------------------------------------
# Rendering statements
# ---------------------------------------------------------------------------

def cpp_target(t):
    k = t[0]
    if k == "vec":
        return "%s[%s]" % (t[1], t[2])
    if k == "mat":
        return "%s(%s,%s)" % (t[1], t[2], t[3])
    if k == "var":
        return _value(t[1], "cpp")
    raise EmitError("unknown target " + repr(t))


def py_target(t):
    k = t[0]
    if k == "vec":
        return "%s[%s]" % (t[1], t[2])
    if k == "mat":
        return "%s[%s, %s]" % (t[1], t[2], t[3])
    if k == "var":
        return _value(t[1], "py")
    raise EmitError("unknown target " + repr(t))


def split_loop(g, stmts, inner_mask):
    """(outer, inner): declarations free of `inner_mask` leaves move out."""
    outer, inner = [], []
    for st in stmts:
        if st[0] == "decl" and not (g.flags[st[2]] & inner_mask):
            outer.append(st)
        else:
            inner.append(st)
    return outer, inner


def _value(v, lang):
    return v[lang] if isinstance(v, dict) else v


def render_cpp(stmts, printer, scalar, indent="    "):
    lines = []
    for st in stmts:
        k = st[0]
        if k == "decl":
            _, name, n, is_bool = st
            ty = "bool" if is_bool else scalar
            lines.append("%sconst %s %s = %s;" % (indent, ty, name,
                                                  printer.expr(n, top=True)))
        elif k == "store":
            _, target, n, op = st
            if target[0] == "ret":
                lines.append("%sreturn %s;" % (indent, printer.expr(n)))
            elif target[0] == "call":
                lines.append(indent + target[1]["cpp"].format(*target[2])
                             % printer.expr(n))
            else:
                lines.append("%s%s %s %s;" % (indent, cpp_target(target), op,
                                              printer.expr(n)))
        elif k == "assign":
            _, target, v = st
            lines.append("%s%s = %s;" % (indent, cpp_target(target),
                                         _value(v, "cpp")))
        elif k == "table":
            _, name, values = st
            lines.append("%sstatic const int %s[] = {%s};"
                         % (indent, name, ",".join(str(int(x)) for x in values) or "0"))
        elif k == "array":
            lines.append("%sstd::vector<%s> %s(%d);" % (indent, scalar, st[1], st[2]))
        elif k == "local":
            lines.append("%s%s %s(0.0);" % (indent, scalar, st[1]))
        elif k == "comment":
            lines.append("%s// %s" % (indent, st[1]))
        elif k == "raw":
            if st[1].get("cpp") is not None:
                lines.append(indent + st[1]["cpp"])
        elif k == "if":
            lines.append("%sif (%s) {" % (indent, st[1]["cpp"]))
            lines += render_cpp(st[2], printer, scalar, indent + "  ")
            lines.append(indent + "}")
        elif k == "for":
            _, var, count, body = st
            lines.append("%sfor (int %s = 0; %s < %s; ++%s) {"
                         % (indent, var, var, count, var))
            lines += render_cpp(body, printer, scalar, indent + "  ")
            lines.append(indent + "}")
        elif k == "loop":
            _, var, lo, hi, body = st
            lines.append("%sfor (int %s = %s; %s < %s; ++%s) {"
                         % (indent, var, lo, var, hi, var))
            lines += render_cpp(body, printer, scalar, indent + "  ")
            lines.append(indent + "}")
        elif k == "block":
            lines += render_cpp(st[2], printer.with_names(st[1]), scalar, indent)
        else:
            raise EmitError("unknown statement " + k)
    return lines


def render_py(stmts, printer, indent="    "):
    lines = []
    for st in stmts:
        k = st[0]
        if k == "decl":
            _, name, n, _ = st
            lines.append("%s%s = %s" % (indent, name, printer.expr(n, top=True)))
        elif k == "store":
            _, target, n, op = st
            if target[0] == "ret":
                lines.append("%sreturn %s" % (indent, printer.expr(n)))
            elif target[0] == "call":
                lines.append(indent + target[1]["py"].format(*target[2])
                             % printer.expr(n))
            else:
                lines.append("%s%s %s %s" % (indent, py_target(target), op,
                                             printer.expr(n)))
        elif k == "assign":
            _, target, v = st
            lines.append("%s%s = %s" % (indent, py_target(target), _value(v, "py")))
        elif k == "table":
            _, name, values = st
            lines.append("%s%s = [%s]" % (indent, name,
                                          ", ".join(str(int(x)) for x in values)))
        elif k == "array":
            lines.append("%s%s = [0.0] * %d" % (indent, st[1], st[2]))
        elif k == "local":
            lines.append("%s%s = 0.0" % (indent, st[1]))
        elif k == "comment":
            lines.append("%s# %s" % (indent, st[1]))
        elif k == "raw":
            if st[1].get("py") is not None:
                lines.append(indent + st[1]["py"])
        elif k == "if":
            lines.append("%sif %s:" % (indent, st[1]["py"]))
            body = render_py(st[2], printer, indent + "    ")
            lines += body or [indent + "    pass"]
        elif k == "block":
            lines += render_py(st[2], printer.with_names(st[1]), indent)
        elif k in ("for", "loop"):
            if k == "for":
                _, var, count, body = st
                rng = "range(%s)" % count
            else:
                _, var, lo, hi, body = st
                rng = "range(%s, %s)" % (lo, hi)
            lines.append("%sfor %s in %s:" % (indent, var, rng))
            b = render_py(body, printer, indent + "    ")
            lines += b or [indent + "    pass"]
        else:
            raise EmitError("unknown statement " + k)
    return lines


# ---------------------------------------------------------------------------
# Slots of the ODE model
# ---------------------------------------------------------------------------

def model_slots(n_states, params_base=None, style="cpp", params="params",
                state="x", time="t", forcing="F", vectors=None):
    """Leaf speller for the flat ODE layout.

    State i -> x[i], initial value i -> params[i], parameter k ->
    params[params_base + k], forcing j -> (*F[j])(t), its rate ->
    F[j]->derivative(t), VEC (name, i) -> vectors.get(name, name)[i],
    map row r -> _lin[r]. An index may be index text. A LOOPVAR leaf
    (LOOPVAR, kind, ...) is the leaf (kind, ...) with index text, or with
    kind 'ref' the text itself.
    """
    base = n_states if params_base is None else params_base
    vectors = vectors or {}
    py = style == "py"

    def slot(at):
        k = at[0]
        if k == cg.LOOPVAR:
            return at[2] if at[1] == "ref" else slot(at[1:])
        if k == cg.STATE:
            return "%s[%s]" % (state, at[1])
        if k == cg.PARAM:
            if isinstance(at[1], str):
                return "%s[%d + %s]" % (params, base, at[1])
            return "%s[%d]" % (params, base + at[1])
        if k == cg.INIT:
            return "%s[%s]" % (params, at[1])
        if k == cg.TIME:
            return time
        if k == cg.FORCING:
            if py:
                return "%s[%s](%s)" % (forcing, at[1], time)
            return "(*%s[%s])(%s)" % (forcing, at[1], time)
        if k == cg.FRATE:
            if py:
                return "%s[%s].derivative(%s)" % (forcing, at[1], time)
            return "%s[%s]->derivative(%s)" % (forcing, at[1], time)
        if k == cg.LINROW:
            return "_lin[%s]" % at[1]
        if k == cg.VEC:
            name = vectors.get(at[1], at[1])
            if len(at) < 3 or at[2] is None:
                return name
            return "%s[%s]" % (name, at[2])
        raise EmitError("no slot for leaf " + cg.KIND_NAMES[k])

    return slot
