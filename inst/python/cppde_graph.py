"""Expression graph of the cppDE code generators: nodes, parser and
symbolic differentiation.

A node is an integer id into the parallel lists `Graph.op`, `Graph.args` and
`Graph.attr`. ADD holds c0 + sum(c_k * t_k) with attr (c0, (c_k, ...)); MUL a
coefficient-free product; POW (base, exponent). Numbers are Fraction or float.
Ids follow construction order and n-ary arguments are sorted by id, so the
graph does not depend on PYTHONHASHSEED.

`Parser` reads model strings with Python's `ast`; syntax or functions without
a rule fall back to SymPy. `AD` provides `jvp`, `vjp`, `time_derivative` and
`lie`.
"""

import ast
import math
import re
from fractions import Fraction

from cppsympy import (derivative_template, from_sympy, normalise_logic,
                      parse_error, safe_sympify)

# ---------------------------------------------------------------------------
# Operations
# ---------------------------------------------------------------------------

NUM, NAMED, BOOL, LEAF, ADD, MUL, POW, CALL, SELECT, CMP, AND, OR, NOT = range(13)
OP_NAMES = ("NUM", "NAMED", "BOOL", "LEAF", "ADD", "MUL", "POW", "CALL",
            "SELECT", "CMP", "AND", "OR", "NOT")

# Leaf kinds. A leaf's attr is a tuple whose first entry is the kind.
STATE, PARAM, INIT, TIME, FORCING, FRATE, VEC, LINROW, TABLE, LOOPVAR = range(10)
KIND_NAMES = ("STATE", "PARAM", "INIT", "TIME", "FORCING", "FRATE", "VEC",
              "LINROW", "TABLE", "LOOPVAR")

# Dependency flags, one bit per leaf kind.
F_STATE, F_PARAM, F_INIT, F_TIME, F_FORCING, F_FRATE, F_VEC, F_LINROW, \
    F_TABLE, F_LOOPVAR = (1 << k for k in range(10))

BOOL_OPS = frozenset((BOOL, CMP, AND, OR, NOT))
CMP_OPS = ("<", "<=", ">", ">=", "==", "!=")

# Nested sums up to this many terms are merged into the enclosing sum.
FLAT_LIMIT = 8

# Upper bound of the tree-size estimate.
_SIZE_CAP = 1 << 30

_INT_LIMIT = 1 << 53

# Largest denominator of a float read as a Fraction.
_DYADIC_LIMIT = 1 << 32


def as_number(v):
    """`v` as Fraction if rational or a float with a small power-of-two
    denominator, else as float.

    Equal numbers share a node whatever their type, so a float equal to a
    Fraction is stored as that Fraction.
    """
    t = type(v)
    if t is Fraction:
        return v
    if t is not float:
        if isinstance(v, Fraction):
            return v
        if isinstance(v, (bool, int)):
            return Fraction(int(v))
        v = float(v)
    if v.is_integer():
        return Fraction(int(v)) if abs(v) < _INT_LIMIT else v
    if math.isfinite(v):
        n, d = v.as_integer_ratio()
        if d <= _DYADIC_LIMIT:
            return Fraction(n, d)
    return v


def is_integer(v):
    return type(v) is Fraction and v.denominator == 1


def num_float(v):
    return float(v)


def _mulnum(a, b):
    if type(a) is Fraction and type(b) is Fraction:
        return a * b
    return as_number(float(a) * float(b))


def _addnum(a, b):
    if type(a) is Fraction and type(b) is Fraction:
        return a + b
    return as_number(float(a) + float(b))


def _pownum(b, e):
    """b**e folded, or None when the result is not a real number."""
    if isinstance(b, Fraction) and is_integer(e):
        n = int(e)
        if b == 0 and n < 0:
            return None
        if abs(n) > 4096:
            return as_number(float(b) ** n) if b != 0 else Fraction(0)
        return b ** n
    if isinstance(b, Fraction) and isinstance(e, Fraction):
        # Exact rational root.
        if b >= 0 and e.denominator <= 64:
            p, q = e.numerator, e.denominator
            rn = _iroot(b.numerator, q)
            rd = _iroot(b.denominator, q)
            if rn is not None and rd is not None:
                return _pownum(Fraction(rn, rd), Fraction(p))
    fb, fe = float(b), float(e)
    if fb < 0 and not float(fe).is_integer():
        return None
    try:
        r = fb ** fe
    except (OverflowError, ZeroDivisionError):
        return None
    if isinstance(r, complex):
        return None
    return as_number(r)


def _iroot(n, q):
    if n < 0:
        return None
    if n in (0, 1):
        return n
    r = round(n ** (1.0 / q))
    for c in (r - 1, r, r + 1):
        if c >= 0 and c ** q == n:
            return c
    return None


# ---------------------------------------------------------------------------
# Numeric functions
# ---------------------------------------------------------------------------

def _heaviside(x):
    return 0.0 if x < 0 else (0.5 if x == 0 else 1.0)


def _sign(x):
    return 1.0 if x > 0 else (-1.0 if x < 0 else 0.0)


def _delta(x):
    """1 at x == 0, else 0 (derivative of Heaviside)."""
    return 1.0 if x == 0 else 0.0


NUMERIC_FUNCS = {
    "exp": math.exp, "log": math.log, "sin": math.sin, "cos": math.cos,
    "tan": math.tan, "asin": math.asin, "acos": math.acos, "atan": math.atan,
    "sinh": math.sinh, "cosh": math.cosh, "tanh": math.tanh,
    "asinh": math.asinh, "acosh": math.acosh, "atanh": math.atanh,
    "abs": abs, "sign": _sign, "floor": math.floor, "ceiling": math.ceil,
    "Heaviside": _heaviside, "delta": _delta, "erf": math.erf,
    "erfc": math.erfc, "gamma": math.gamma, "loggamma": math.lgamma,
    "atan2": math.atan2, "min": min, "max": max, "_prod": lambda x: x,
}

# Functions with a built-in derivative rule.
ELEMENTARY = frozenset((
    "exp", "log", "sin", "cos", "tan", "asin", "acos", "atan", "sinh", "cosh",
    "tanh", "asinh", "acosh", "atanh", "abs", "sign", "floor", "ceiling",
    "Heaviside", "delta", "min", "max"))


class GraphError(ValueError):
    pass


class Graph:
    """Expression DAG.

    `flags[n]` is the bitmask of leaf kinds below node n, `size[n]` a capped
    tree size.
    """

    def __init__(self):
        self.op = []
        self.args = []
        self.attr = []
        self.flags = []
        self.size = []
        self._key = {}
        self._bits = {}
        self.ZERO = self.num(0)
        self.ONE = self.num(1)
        self.MINUS_ONE = self.num(-1)
        self.HALF = self.num(Fraction(1, 2))
        self.TRUE = self._node(BOOL, (), True)
        self.FALSE = self._node(BOOL, (), False)

    def __len__(self):
        return len(self.op)

    # -- raw construction -------------------------------------------------

    def _node(self, op, args, attr):
        key = (op, args, attr)
        nid = self._key.get(key)
        if nid is not None:
            return nid
        nid = len(self.op)
        self.op.append(op)
        self.args.append(args)
        self.attr.append(attr)
        if op == LEAF:
            fl = 1 << attr[0]
            sz = 1
        else:
            fl = 0
            sz = 1
            for a in args:
                fl |= self.flags[a]
                sz += self.size[a]
            if sz > _SIZE_CAP:
                sz = _SIZE_CAP
        self.flags.append(fl)
        self.size.append(sz)
        self._key[key] = nid
        return nid

    # -- queries ----------------------------------------------------------

    def is_num(self, n):
        return self.op[n] == NUM

    def value(self, n):
        return self.attr[n]

    def is_zero(self, n):
        return self.op[n] == NUM and self.attr[n] == 0

    def is_one(self, n):
        return self.op[n] == NUM and self.attr[n] == 1

    def is_bool(self, n):
        return self.op[n] in BOOL_OPS

    def kind(self, n):
        return self.attr[n][0] if self.op[n] == LEAF else None

    def depends(self, n, mask):
        return bool(self.flags[n] & mask)

    def bits(self, n, kind):
        """Bitset of the indices of `kind` leaves below node n."""
        return self._fold_bits(n, kind, ("s", kind), self._kids)

    def dbits(self, n, kind):
        """Bitset of `kind` leaves node n depends on differentiably.

        Conditions of selects, comparisons and logic, and the arguments of
        floor, ceiling and delta do not count.
        """
        return self._fold_bits(n, kind, ("d", kind), self._dkids)

    def _kids(self, m):
        return self.args[m]

    def _dkids(self, m):
        o = self.op[m]
        if o in BOOL_OPS:
            return ()
        if o == SELECT:
            return self.args[m][1:]
        if o == CALL and self.attr[m] in ("floor", "ceiling", "delta"):
            return ()
        return self.args[m]

    def _fold_bits(self, n, kind, key, kids):
        memo = self._bits.setdefault(key, {})
        hit = memo.get(n)
        if hit is not None:
            return hit
        want = 1 << kind
        flags = self.flags
        if not (flags[n] & want):
            return 0
        stack = [(n, iter(kids(n)))]
        while stack:
            m, it = stack[-1]
            pushed = False
            for a in it:
                if a in memo:
                    continue
                if not (flags[a] & want):
                    memo[a] = 0
                    continue
                if self.op[a] == LEAF:
                    memo[a] = 1 << self.attr[a][1] if self.attr[a][0] == kind else 0
                    continue
                stack.append((a, iter(kids(a))))
                pushed = True
                break
            if pushed:
                continue
            stack.pop()
            if m in memo:
                continue
            if self.op[m] == LEAF:
                memo[m] = 1 << self.attr[m][1] if self.attr[m][0] == kind else 0
                continue
            b = 0
            for a in kids(m):
                b |= memo.get(a, 0)
            memo[m] = b
        return memo[n]

    # -- leaves and constants ---------------------------------------------

    def num(self, v):
        v = as_number(v)
        if isinstance(v, float) and v != v:
            # A single NaN node.
            return self._node(NUM, (), "nan")
        return self._node(NUM, (), v)

    def named(self, name, value):
        return self._node(NAMED, (), (name, float(value)))

    def boolean(self, b):
        return self.TRUE if b else self.FALSE

    def leaf(self, kind, *idx):
        return self._node(LEAF, (), (kind,) + tuple(idx))

    def state(self, i):
        return self.leaf(STATE, i)

    def param(self, k):
        return self.leaf(PARAM, k)

    def init(self, i):
        return self.leaf(INIT, i)

    def time(self):
        return self.leaf(TIME)

    def forcing(self, j):
        return self.leaf(FORCING, j)

    def frate(self, j):
        return self.leaf(FRATE, j)

    def vec(self, name, i):
        return self.leaf(VEC, name, i)

    def numval(self, n):
        """Value of a NUM or NAMED node, else None."""
        o = self.op[n]
        if o == NUM:
            v = self.attr[n]
            return float("nan") if v == "nan" else v
        if o == NAMED:
            return self.attr[n][1]
        return None

    # -- sums -------------------------------------------------------------

    def add_terms(self, terms, const=0):
        """Node for const + sum(c * t for c, t in terms)."""
        const = as_number(const)
        acc = {}
        order = []
        for c, t in terms:
            c = as_number(c)
            if c == 0:
                continue
            o = self.op[t]
            if o == NUM:
                v = self.attr[t]
                if v == "nan":
                    const = float("nan")
                else:
                    const = _addnum(const, _mulnum(c, v))
                continue
            if o == ADD and len(self.args[t]) <= FLAT_LIMIT:
                c0, cs = self.attr[t]
                const = _addnum(const, _mulnum(c, c0))
                for ci, ti in zip(cs, self.args[t]):
                    self._acc(acc, order, ti, _mulnum(c, ci))
                continue
            self._acc(acc, order, t, c)
        items = [(t, acc[t]) for t in order if acc[t] != 0]
        if not items:
            return self.num(const)
        if len(items) == 1 and const == 0 and items[0][1] == 1:
            return items[0][0]
        items.sort()
        return self._node(ADD, tuple(t for t, _ in items),
                          (const, tuple(c for _, c in items)))

    @staticmethod
    def _acc(acc, order, t, c):
        if t in acc:
            acc[t] = _addnum(acc[t], c)
        else:
            acc[t] = c
            order.append(t)

    def add(self, *nodes):
        return self.add_terms([(1, n) for n in nodes])

    def sum(self, nodes):
        return self.add_terms([(1, n) for n in nodes])

    def sub(self, a, b):
        return self.add_terms([(1, a), (-1, b)])

    def neg(self, a):
        return self.add_terms([(-1, a)])

    def scale(self, c, a):
        return self.add_terms([(c, a)])

    def split_scale(self, n):
        """(c, t) with n == c * t and t free of a numeric factor."""
        if self.op[n] == ADD:
            c0, cs = self.attr[n]
            if c0 == 0 and len(cs) == 1:
                return cs[0], self.args[n][0]
        if self.op[n] == NUM:
            return self.attr[n], self.ONE
        return Fraction(1), n

    def add_const(self, n):
        return self.attr[n][0] if self.op[n] == ADD else 0

    def add_items(self, n):
        """(const, [(c, t), ...]) of node n read as a sum."""
        o = self.op[n]
        if o == ADD:
            c0, cs = self.attr[n]
            return c0, list(zip(cs, self.args[n]))
        if o == NUM:
            return self.attr[n], []
        return Fraction(0), [(Fraction(1), n)]

    # -- products ---------------------------------------------------------

    def exponent(self, e):
        """Exponent in mul_factors form: a number or ('n', node)."""
        if isinstance(e, tuple):
            return e
        if isinstance(e, (Fraction, float)):
            return as_number(e)
        # a node id
        if self.op[e] == NUM and self.attr[e] != "nan":
            return self.attr[e]
        return ("n", e)

    def _exp_node(self, e):
        return e[1] if isinstance(e, tuple) else self.num(e)

    def mul_factors(self, factors, _again=True):
        """Node for prod(b ** e for b, e in factors).

        `e` is a number (int, Fraction, float) or ('n', node); an int is a
        number, not a node id.
        """
        coef = Fraction(1)
        bases = {}
        order = []
        exp_terms = []
        stack = []
        for f, e in reversed(factors):
            if not isinstance(e, tuple):
                e = as_number(e)
            stack.append((f, e))
        while stack:
            f, e = stack.pop()
            if isinstance(e, tuple):
                en = e[1]
                if self.op[en] == NUM and self.attr[en] != "nan":
                    e = self.attr[en]
            if isinstance(e, tuple):
                if self.op[f] == NUM and self.attr[f] == 1:
                    continue
                self._acc_base(bases, order, f, e)
                continue
            if e == 0:
                continue
            if self._mul_numeric_exp(f, e, stack, exp_terms):
                continue
            if self.op[f] == NUM:
                v = self.attr[f]
                if v == "nan":
                    return self.num(float("nan"))
                r = _pownum(v, e)
                if r is None:
                    self._acc_base(bases, order, f, e)
                else:
                    coef = _mulnum(coef, r)
                continue
            self._acc_base(bases, order, f, e)
        out = []
        for b in order:
            e = bases[b]
            if not isinstance(e, tuple) and e == 0:
                continue
            out.append(self._pow_raw(b, self._exp_node(e)))
        if exp_terms:
            arg = self.add_terms(exp_terms)
            if not self.is_zero(arg):
                out.append(self.call("exp", arg))
        flat = []
        nested = False
        for n in out:
            o = self.op[n]
            if o == NUM:
                v = self.attr[n]
                if v == "nan":
                    return self.num(float("nan"))
                coef = _mulnum(coef, v)
            elif o == ADD:
                c, t = self.split_scale(n)
                if self.op[t] == ADD:
                    flat.append(n)
                else:
                    coef = _mulnum(coef, c)
                    if t != self.ONE:
                        flat.append(t)
                        nested = nested or self.op[t] == MUL
            else:
                flat.append(n)
                nested = nested or o == MUL
        if coef == 0:
            return self.ZERO
        if nested and _again:
            body = self.mul_factors([(n, 1) for n in flat], _again=False)
            return self.scale(coef, body) if coef != 1 else body
        if not flat:
            return self.num(coef)
        if len(flat) == 1:
            body = flat[0]
        else:
            body = self._node(MUL, tuple(sorted(flat)), None)
        if coef == 1:
            return body
        return self.add_terms([(coef, body)])

    def _mul_numeric_exp(self, f, e, stack, exp_terms):
        """Expand factor f under numeric exponent e; True if consumed."""
        o = self.op[f]
        if o == ADD:
            c0, cs = self.attr[f]
            if c0 == 0 and len(cs) == 1:
                c = cs[0]
                if is_integer(e) or c > 0:
                    r = _pownum(c, e)
                    if r is not None:
                        stack.append((self.args[f][0], e))
                        stack.append((self.num(r), Fraction(1)))
                        return True
            return False
        if o == MUL:
            if is_integer(e):
                for a in reversed(self.args[f]):
                    stack.append((a, e))
                return True
            return False
        if o == POW:
            if is_integer(e):
                b, pe = self.args[f]
                if self.op[pe] == NUM and self.attr[pe] != "nan":
                    stack.append((b, _mulnum(self.attr[pe], e)))
                else:
                    stack.append((b, ("n", self.scale(e, pe))))
                return True
            return False
        if o == CALL and self.attr[f] == "exp":
            exp_terms.append((e, self.args[f][0]))
            return True
        return False

    def _acc_base(self, bases, order, b, e):
        if b not in bases:
            bases[b] = e
            order.append(b)
            return
        prev = bases[b]
        if not isinstance(prev, tuple) and not isinstance(e, tuple):
            bases[b] = _addnum(prev, e)
            return
        s = self.add(self._exp_node(prev), self._exp_node(e))
        bases[b] = self.exponent(s)

    def mul(self, *nodes):
        if len(nodes) == 2:
            # number times leaf: the common term of a long sum
            a, b = nodes
            if self.op[a] == NUM:
                a, b = b, a
            if (self.op[b] == NUM and self.op[a] == LEAF
                    and self.attr[b] != "nan"):
                c = self.attr[b]
                return self.ZERO if c == 0 else self.add_terms([(c, a)])
        return self.mul_factors([(n, 1) for n in nodes])

    def prod(self, nodes):
        return self.mul_factors([(n, 1) for n in nodes])

    def div(self, a, b):
        return self.mul_factors([(a, 1), (b, -1)])

    def inv(self, a):
        return self.mul_factors([(a, -1)])

    def group(self, n):
        """n itself, or for a product an identity node that later products
        keep as one factor, so the product can be shared."""
        if self.op[n] != MUL:
            return n
        return self._node(CALL, (n,), "_prod")

    def pow_split(self, n):
        if self.op[n] == POW:
            return self.args[n]
        return n, self.ONE

    def mul_items(self, n):
        """(c, [(base, exponent node), ...]) of node n read as a product."""
        coef, t = self.split_scale(n)
        if t == self.ONE:
            return coef, []
        if self.op[t] == MUL:
            return coef, [self.pow_split(a) for a in self.args[t]]
        return coef, [self.pow_split(t)]

    # -- powers -----------------------------------------------------------

    def pow(self, b, e):
        """Node for b ** e; e is a node id, a Fraction or a float."""
        if isinstance(e, (Fraction, float)):
            return self.mul_factors([(b, e)])
        return self.mul_factors([(b, self.exponent(e))])

    def _pow_raw(self, b, e):
        if self.op[e] == NUM:
            v = self.attr[e]
            if v == 1:
                return b
            if v == 0:
                return self.ONE
        return self._node(POW, (b, e), None)

    def sqrt(self, a):
        return self.pow(a, Fraction(1, 2))

    # -- functions --------------------------------------------------------

    def call(self, name, *args):
        args = tuple(args)
        vals = [self.numval(a) for a in args]
        allnum = all(v is not None for v in vals)
        if name in ("min", "max"):
            return self._minmax(name, args)
        if allnum and name in NUMERIC_FUNCS:
            r = self._fold_call(name, vals)
            if r is not None:
                return r
        a0 = args[0] if args else None
        if name == "exp":
            if self.op[a0] == CALL and self.attr[a0] == "log":
                return self.args[a0][0]
        elif name == "log":
            if self.op[a0] == CALL and self.attr[a0] == "exp":
                return self.args[a0][0]
        elif name == "abs":
            o = self.op[a0]
            if o == CALL and self.attr[a0] in ("abs", "exp"):
                return a0
            c, t = self.split_scale(a0)
            if c != 1 and t != self.ONE:
                return self.scale(abs(c), self.call("abs", t))
        elif name in ("floor", "ceiling") and self.op[a0] == CALL \
                and self.attr[a0] in ("floor", "ceiling"):
            return a0
        return self._node(CALL, args, name)

    def _fold_call(self, name, vals):
        if name == "exp" and vals[0] == 0:
            return self.ONE
        if name == "log" and vals[0] == 1:
            return self.ZERO
        if name in ("sin", "tan", "asin", "atan", "sinh", "tanh", "asinh",
                    "atanh", "erf") and vals[0] == 0:
            return self.ZERO
        if name in ("cos", "cosh") and vals[0] == 0:
            return self.ONE
        if name in ("abs", "sign", "floor", "ceiling", "Heaviside", "delta"):
            v = vals[0]
            if isinstance(v, Fraction):
                if name == "abs":
                    return self.num(abs(v))
                if name == "sign":
                    return self.num((v > 0) - (v < 0))
                if name == "floor":
                    return self.num(math.floor(v))
                if name == "ceiling":
                    return self.num(math.ceil(v))
                if name == "Heaviside":
                    return self.num(Fraction(0) if v < 0 else
                                    (Fraction(1, 2) if v == 0 else Fraction(1)))
                return self.num(1 if v == 0 else 0)
        try:
            r = NUMERIC_FUNCS[name](*[float(v) for v in vals])
        except (ValueError, OverflowError, ZeroDivisionError):
            return None
        if isinstance(r, complex) or not math.isfinite(r):
            return None
        return self.num(float(r))

    def _minmax(self, name, args):
        flat = []
        for a in args:
            if self.op[a] == CALL and self.attr[a] == name:
                flat.extend(self.args[a])
            else:
                flat.append(a)
        nums = [self.numval(a) for a in flat]
        consts = [v for v in nums if v is not None]
        rest = []
        for a, v in zip(flat, nums):
            if v is None and a not in rest:
                rest.append(a)
        if consts:
            pick = min if name == "min" else max
            c = pick(consts, key=float)
            rest.append(self.num(c))
        if len(rest) == 1:
            return rest[0]
        # Binary nesting in argument order.
        out = rest[0]
        for a in rest[1:]:
            out = self._node(CALL, (out, a), name)
        return out

    # -- logic ------------------------------------------------------------

    def select(self, c, a, b):
        if c == self.TRUE:
            return a
        if c == self.FALSE:
            return b
        if a == b:
            return a
        return self._node(SELECT, (c, a, b), None)

    def cmp(self, op, a, b):
        va, vb = self.numval(a), self.numval(b)
        if va is not None and vb is not None:
            fa, fb = float(va), float(vb)
            r = {"<": fa < fb, "<=": fa <= fb, ">": fa > fb, ">=": fa >= fb,
                 "==": fa == fb, "!=": fa != fb}[op]
            return self.boolean(r)
        if a == b:
            return self.boolean(op in ("<=", ">=", "=="))
        return self._node(CMP, (a, b), op)

    def logic(self, op, args):
        flat = []
        absorb = self.FALSE if op == AND else self.TRUE
        neutral = self.TRUE if op == AND else self.FALSE
        for a in args:
            if self.op[a] == op:
                items = self.args[a]
            else:
                items = (a,)
            for x in items:
                if x == absorb:
                    return absorb
                if x == neutral or x in flat:
                    continue
                flat.append(x)
        if not flat:
            return neutral
        if len(flat) == 1:
            return flat[0]
        return self._node(op, tuple(sorted(flat)), None)

    def and_(self, *args):
        return self.logic(AND, args)

    def or_(self, *args):
        return self.logic(OR, args)

    def not_(self, a):
        if a == self.TRUE:
            return self.FALSE
        if a == self.FALSE:
            return self.TRUE
        if self.op[a] == NOT:
            return self.args[a][0]
        return self._node(NOT, (a,), None)

    def truth(self, n):
        """n as a condition: n itself if boolean, else n != 0."""
        if self.is_bool(n):
            return n
        return self.cmp("!=", n, self.ZERO)

    # -- traversal --------------------------------------------------------

    def topo(self, roots):
        """Nodes reachable from `roots`, children before parents."""
        seen = set()
        out = []
        for r in roots:
            if r in seen:
                continue
            stack = [(r, 0)]
            while stack:
                n, i = stack.pop()
                if i == 0:
                    if n in seen:
                        continue
                a = self.args[n]
                if i < len(a):
                    stack.append((n, i + 1))
                    c = a[i]
                    if c not in seen:
                        stack.append((c, 0))
                else:
                    if n not in seen:
                        seen.add(n)
                        out.append(n)
        return out

    def substitute(self, roots, mapping):
        """`roots` rebuilt with `mapping` (node -> node) applied."""
        memo = dict(mapping)
        for n in self.topo(roots):
            if n in memo:
                continue
            args = self.args[n]
            if not args:
                memo[n] = n
                continue
            new = tuple(memo[a] for a in args)
            if new == args:
                memo[n] = n
            else:
                memo[n] = self.rebuild(n, new)
        return [memo[r] for r in roots]

    def rebuild(self, n, args):
        """Operation of node n applied to `args`."""
        o = self.op[n]
        at = self.attr[n]
        if o == ADD:
            return self.add_terms(list(zip(at[1], args)), at[0])
        if o == MUL:
            return self.prod(args)
        if o == POW:
            return self.pow(args[0], args[1])
        if o == CALL:
            return self.call(at, *args)
        if o == SELECT:
            return self.select(*args)
        if o == CMP:
            return self.cmp(at, *args)
        if o in (AND, OR):
            return self.logic(o, args)
        if o == NOT:
            return self.not_(args[0])
        return self._node(o, args, at)

    def evaluate(self, roots, env):
        """Float values of `roots`; `env` maps leaf attr to float."""
        val = {}
        for n in self.topo(roots):
            val[n] = self._eval_node(n, val, env)
        return [val[r] for r in roots]

    def _eval_node(self, n, val, env):
        o = self.op[n]
        at = self.attr[n]
        a = self.args[n]
        if o == NUM:
            return float("nan") if at == "nan" else float(at)
        if o == NAMED:
            return at[1]
        if o == BOOL:
            return at
        if o == LEAF:
            return env[at]
        if o == ADD:
            s = float(at[0])
            for c, t in zip(at[1], a):
                s += float(c) * val[t]
            return s
        if o == MUL:
            p = 1.0
            for t in a:
                p *= val[t]
            return p
        if o == POW:
            b, e = val[a[0]], val[a[1]]
            try:
                r = b ** e
            except ZeroDivisionError:
                return math.inf
            except OverflowError:
                return math.inf
            return float("nan") if isinstance(r, complex) else r
        if o == CALL:
            args = [val[t] for t in a]
            try:
                return float(NUMERIC_FUNCS[at](*args))
            except (ValueError, OverflowError):
                return float("nan")
        if o == SELECT:
            return val[a[1]] if val[a[0]] else val[a[2]]
        if o == CMP:
            x, y = val[a[0]], val[a[1]]
            return {"<": x < y, "<=": x <= y, ">": x > y, ">=": x >= y,
                    "==": x == y, "!=": x != y}[at]
        if o == AND:
            return all(val[t] for t in a)
        if o == OR:
            return any(val[t] for t in a)
        if o == NOT:
            return not val[a[0]]
        raise GraphError("cannot evaluate " + OP_NAMES[o])

    def describe(self, n, names=None, depth=0):
        """Infix string of node n."""
        o = self.op[n]
        at = self.attr[n]
        a = self.args[n]
        d = lambda m: self.describe(m, names, depth + 1)
        if o == NUM:
            return str(at)
        if o == NAMED:
            return at[0]
        if o == BOOL:
            return "True" if at else "False"
        if o == LEAF:
            if names is not None and at in names:
                return names[at]
            return KIND_NAMES[at[0]].lower() + "".join("[%s]" % (x,) for x in at[1:])
        if o == ADD:
            parts = [str(at[0])] if at[0] != 0 else []
            parts += ["%s*%s" % (c, d(t)) for c, t in zip(at[1], a)]
            return "(" + " + ".join(parts) + ")"
        if o == MUL:
            return "(" + "*".join(d(t) for t in a) + ")"
        if o == POW:
            return "%s**%s" % (d(a[0]), d(a[1]))
        if o == CALL:
            return "%s(%s)" % (at, ", ".join(d(t) for t in a))
        if o == SELECT:
            return "select(%s, %s, %s)" % tuple(d(t) for t in a)
        if o == CMP:
            return "(%s %s %s)" % (d(a[0]), at, d(a[1]))
        if o in (AND, OR):
            return "(" + (" & " if o == AND else " | ").join(d(t) for t in a) + ")"
        if o == NOT:
            return "!" + d(a[0])
        return OP_NAMES[o]


# ===========================================================================
# Parser
# ===========================================================================

# Bare names read as numbers.
CONSTANTS = {
    'pi': math.pi,
    'E': math.e,
    'oo': math.inf,
    'euler_gamma': 0.5772156649015329,
}


class Fallback(Exception):
    """Raised when the expression needs the SymPy fallback."""


def _unary(name):
    return lambda g, a: g.call(name, _one(a, name))


def _one(a, name):
    if len(a) != 1:
        raise ValueError("%s() takes one argument, got %d" % (name, len(a)))
    return a[0]


def _log(g, a):
    if len(a) == 1:
        return g.call('log', a[0])
    if len(a) == 2:
        return g.div(g.call('log', a[0]), g.call('log', a[1]))
    raise ValueError("log() takes one or two arguments")


def _logb(base):
    return lambda g, a: g.div(g.call('log', _one(a, 'log')),
                              g.call('log', g.num(base)))


def _expb(base):
    return lambda g, a: g.call(
        'exp', g.mul(_one(a, 'exp'), g.call('log', g.num(base))))


def _recip(name):
    return lambda g, a: g.inv(g.call(name, _one(a, name)))


def _arc_recip(name):
    return lambda g, a: g.call(name, g.inv(_one(a, name)))


def _power(p):
    return lambda g, a: g.pow(_one(a, 'root'), p)


def _root(g, a):
    if len(a) != 2:
        raise Fallback()
    n = g.numval(a[1])
    if n is not None and isinstance(n, Fraction) and n != 0:
        return g.pow(a[0], 1 / n)
    return g.pow(a[0], g.inv(a[1]))


def _pow(g, a):
    if len(a) != 2:
        raise ValueError("pow() takes two arguments")
    return g.pow(a[0], a[1])


def _minmax(name):
    def f(g, a):
        if not a:
            raise ValueError("%s() needs an argument" % name)
        return g.call(name, *a)
    return f


def _round(g, a):
    return g.call('floor', g.add(_one(a, 'round'), g.HALF))


def _heaviside(g, a):
    if len(a) != 1:
        raise Fallback()
    return g.call('Heaviside', a[0])


def _atan2(g, a):
    if len(a) != 2:
        raise ValueError("atan2() takes two arguments")
    return g.call('atan2', a[0], a[1])


def _factorial(g, a):
    return g.call('gamma', g.add(_one(a, 'factorial'), g.ONE))


def _dirac(g, a):
    raise ValueError("DiracDelta has no value a model can evaluate")


def _logic(op):
    def f(g, a):
        return g.logic(op, [g.truth(x) for x in a])
    return f


def _not(g, a):
    return g.not_(g.truth(_one(a, 'Not')))


def _piecewise_flat(g, a):
    """piecewise(v1, c1, ..., otherwise) as nested selects. Without an
    otherwise branch SymPy decides whether the conditions are exhaustive."""
    if len(a) % 2 == 0:
        raise Fallback()
    out = a[-1]
    for i in range(len(a) - 3, -1, -2):
        out = g.select(g.truth(a[i + 1]), a[i], out)
    return out


FUNCS = {
    'exp': _unary('exp'), 'log': _log, 'ln': _unary('log'),
    'exp10': _expb(10), 'exp2': _expb(2),
    'log10': _logb(10), 'log2': _logb(2),
    'sin': _unary('sin'), 'cos': _unary('cos'), 'tan': _unary('tan'),
    'asin': _unary('asin'), 'acos': _unary('acos'), 'atan': _unary('atan'),
    'sinh': _unary('sinh'), 'cosh': _unary('cosh'), 'tanh': _unary('tanh'),
    'asinh': _unary('asinh'), 'acosh': _unary('acosh'),
    'atanh': _unary('atanh'),
    'cot': _recip('tan'), 'sec': _recip('cos'), 'csc': _recip('sin'),
    'coth': _recip('tanh'), 'sech': _recip('cosh'), 'csch': _recip('sinh'),
    'acot': _arc_recip('atan'), 'asec': _arc_recip('acos'),
    'acsc': _arc_recip('asin'), 'acoth': _arc_recip('atanh'),
    'asech': _arc_recip('acosh'), 'acsch': _arc_recip('asinh'),
    'atan2': _atan2,
    'sqrt': _power(Fraction(1, 2)), 'cbrt': _power(Fraction(1, 3)),
    'root': _root, 'pow': _pow,
    'abs': _unary('abs'), 'Abs': _unary('abs'), 'sign': _unary('sign'),
    'floor': _unary('floor'), 'ceiling': _unary('ceiling'), 'round': _round,
    'min': _minmax('min'), 'max': _minmax('max'),
    'Min': _minmax('min'), 'Max': _minmax('max'),
    'Heaviside': _heaviside, 'DiracDelta': _dirac,
    'erf': _unary('erf'), 'erfc': _unary('erfc'),
    'gamma': _unary('gamma'), 'loggamma': _unary('loggamma'),
    'factorial': _factorial,
    'And': _logic(AND), 'Or': _logic(OR), 'Not': _not,
    'piecewise': _piecewise_flat,
}

_AST_CMP = {ast.Lt: '<', ast.LtE: '<=', ast.Gt: '>', ast.GtE: '>=',
        ast.Eq: '==', ast.NotEq: '!='}

# Shortest source text of a cached subexpression.
_SUB_MIN = 64

# Chains with at least this many +/- operators are rewritten by _flat_sums.
_SUM_MIN = 64
_SUM_CALL = "__cppde_sum__"
_TOKEN = re.compile(r"\s*(?:(\d+\.?\d*(?:[eE][+-]?\d+)?|\.\d+(?:[eE][+-]?\d+)?)"
                    r"|([A-Za-z_]\w*)|(\*\*|//|<=|>=|==|!=|<<|>>|\S))")
# Tokens binding tighter than binary + and -.
_TIGHT = frozenset(("*", "/", "//", "%", "@", "**", "~", "+", "-", ".", "!"))
_KEYWORDS = frozenset(("and", "or", "not", "in", "is", "if", "else", "lambda"))


def _flat_sums(src):
    """`src` with long +/- chains as _SUM_CALL(+(t1), -(t2), ...), a flat
    call that Parser._walk reads as the same sum."""
    # frames: [items of the current segment, finished text, closing bracket]
    stack = [[[], [], None]]
    pos = 0
    n = len(src)
    while pos < n:
        m = _TOKEN.match(src, pos)
        if m is None or m.end() == pos:
            break
        pos = m.end()
        num, name, op = m.groups()
        frame = stack[-1]
        if num is not None or (name is not None and name not in _KEYWORDS):
            frame[0].append(("x", num if num is not None else name))
        elif op in ("(", "["):
            stack.append([[], [], ")" if op == "(" else "]"])
        elif op in (")", "]"):
            _flush(frame)
            stack.pop()
            if not stack or frame[2] != op:
                raise SyntaxError("unbalanced brackets")
            text = ("(" if op == ")" else "[") + " ".join(frame[1]) + op
            stack[-1][0].append(("x", text))
        elif op in _TIGHT:
            items = frame[0]
            binary = op in ("+", "-") and items and items[-1][0] == "x"
            items.append(("s" if binary else "o", op))
        else:
            _flush(frame)
            frame[1].append(op if op is not None else name)
    if src[pos:].strip() or len(stack) != 1:
        raise SyntaxError("cannot split the expression")
    _flush(stack[0])
    return " ".join(stack[0][1])


def _literal(n):
    """Value of a numeric literal or its negation, else None."""
    neg = type(n) is ast.UnaryOp and type(n.op) is ast.USub
    if neg:
        n = n.operand
    if type(n) is not ast.Constant or type(n.value) not in (int, float):
        return None
    v = as_number(n.value)
    return _addnum(Fraction(0), _mulnum(Fraction(-1), v)) if neg else v


def _is_sum_call(n):
    return (type(n) is ast.Call and isinstance(n.func, ast.Name)
            and n.func.id == _SUM_CALL)


def _flush(frame):
    items = frame[0]
    frame[0] = []
    if not items:
        return
    if sum(1 for k, _ in items if k == "s") < _SUM_MIN:
        frame[1].append(" ".join(t for _, t in items))
        return
    terms = []
    sign = "+"
    cur = []
    for k, t in items:
        if k == "s":
            terms.append("%s(%s)" % (sign, " ".join(cur)))
            sign, cur = t, []
        else:
            cur.append(t)
    terms.append("%s(%s)" % (sign, " ".join(cur)))
    frame[1].append("%s(%s)" % (_SUM_CALL, ", ".join(terms)))


class Parser:
    """Parser bound to one graph and one symbol table.

    Args:
        graph: Graph receiving the nodes.
        symbols: dict name -> node id.

    Results of `parse` are cached per string.
    """

    def __init__(self, graph, symbols):
        self.g = graph
        self.symbols = dict(symbols)
        self._cache = {}
        self._sub = {}
        self._src = None

    def parse(self, text, label=None):
        key = str(text).strip()
        hit = self._cache.get(key)
        if hit is not None:
            return hit
        node = self._parse(key, label)
        self._cache[key] = node
        return node

    def _parse(self, text, label):
        if text == "":
            raise parse_error(text, "empty expression", label=label)
        try:
            src = normalise_logic(text).replace('^', '**')
            try:
                tree = ast.parse(src, mode='eval')
            except RecursionError:
                src = _flat_sums(src)
                tree = ast.parse(src, mode='eval')
        except SyntaxError:
            return self._sympy(text, label)
        self._src = src if (src.isascii() and '\n' not in src) else None
        try:
            return self._walk(tree.body)
        except Fallback:
            return self._sympy(text, label)
        except (ValueError, TypeError, ZeroDivisionError) as e:
            raise parse_error(text, e, label=label) from None
        finally:
            self._src = None

    def _sympy(self, text, label):
        local = {}
        import sympy as sp
        for name in self.symbols:
            local[name] = sp.Symbol(name, real=True)
        expr = safe_sympify(text, local, label=label)
        try:
            return from_sympy(self.g, expr, self.symbols)
        except ValueError as e:
            raise parse_error(text, e, label=label) from None

    # -- the walk ---------------------------------------------------------

    def _walk(self, n):
        t = type(n)
        if (t is ast.BinOp and self._src is not None
                and n.end_col_offset - n.col_offset >= _SUB_MIN):
            # equal source text, equal node: repeated subexpressions are
            # walked once
            key = self._src[n.col_offset:n.end_col_offset]
            hit = self._sub.get(key)
            if hit is None:
                hit = self._sub[key] = self._walk_node(n, t)
            return hit
        return self._walk_node(n, t)

    def _walk_node(self, n, t):
        g = self.g
        if t is ast.BinOp and type(n.op) in (ast.Add, ast.Sub):
            terms = []
            self._terms(n, 1, terms)
            return g.add_terms(terms)
        if t is ast.BinOp:
            a = self._walk(n.left)
            b = self._walk(n.right)
            op = type(n.op)
            if op is ast.Mult:
                return g.mul(a, b)
            if op is ast.Div:
                if g.is_zero(b):
                    raise ValueError("division by zero")
                return g.div(a, b)
            if op is ast.Pow:
                return g.pow(a, b)
            if op is ast.BitAnd:
                return g.and_(g.truth(a), g.truth(b))
            if op is ast.BitOr:
                return g.or_(g.truth(a), g.truth(b))
            raise Fallback()
        if t is ast.UnaryOp:
            c = _literal(n)
            if c is not None:
                return g.num(c)
            a = self._walk(n.operand)
            op = type(n.op)
            if op is ast.USub:
                return g.neg(a)
            if op is ast.UAdd:
                return a
            if op in (ast.Not, ast.Invert):
                return g.not_(g.truth(a))
            raise Fallback()
        if t is ast.Constant:
            v = n.value
            if isinstance(v, bool):
                return g.boolean(v)
            if isinstance(v, int):
                return g.num(Fraction(v))
            if isinstance(v, float):
                return g.num(v)
            raise ValueError("unsupported literal %r" % (v,))
        if t is ast.Name:
            node = self.symbols.get(n.id)
            if node is not None:
                return node
            if n.id in CONSTANTS:
                return g.num(CONSTANTS[n.id])
            raise ValueError("unknown symbol '%s'" % n.id)
        if t is ast.Call:
            if not isinstance(n.func, ast.Name) or n.keywords:
                raise Fallback()
            name = n.func.id
            if name == 'Piecewise':
                return self._piecewise(n.args)
            if name == _SUM_CALL:
                return g.add_terms(self._sum_terms(n.args))
            f = FUNCS.get(name)
            if f is None:
                raise Fallback()
            return f(g, [self._walk(a) for a in n.args])
        if t is ast.Compare:
            left = self._walk(n.left)
            parts = []
            for op, right in zip(n.ops, n.comparators):
                sym = _AST_CMP.get(type(op))
                if sym is None:
                    raise Fallback()
                r = self._walk(right)
                parts.append(g.cmp(sym, left, r))
                left = r
            return parts[0] if len(parts) == 1 else g.and_(*parts)
        if t is ast.BoolOp:
            op = AND if isinstance(n.op, ast.And) else OR
            return g.logic(op, [g.truth(self._walk(v)) for v in n.values])
        raise Fallback()

    def _terms(self, n, sign, out):
        """Append (sign, node) for each operand of a +/- chain."""
        while type(n) is ast.BinOp and type(n.op) in (ast.Add, ast.Sub):
            s = -sign if type(n.op) is ast.Sub else sign
            out.append(self._term(n.right, s))
            n = n.left
        if _is_sum_call(n):
            out.extend((sign * s, t) for s, t in reversed(self._sum_terms(n.args)))
        else:
            out.append(self._term(n, sign))
        out.reverse()

    def _term(self, n, sign):
        """(coefficient, node) of one operand of a sum; a literal times a
        symbol is read without building the product."""
        if type(n) is ast.BinOp and type(n.op) is ast.Mult:
            a, b = n.left, n.right
            if type(a) is ast.Name:
                a, b = b, a
            if type(b) is ast.Name:
                c = _literal(a)
                node = self.symbols.get(b.id)
                if c is not None and node is not None:
                    return _mulnum(Fraction(sign), c), node
        return sign, self._walk(n)

    def _sum_terms(self, args):
        """Terms of a _SUM_CALL, a leading chain spliced in."""
        out = []
        rest = args
        head = args[0].operand
        if type(args[0].op) is ast.UAdd:
            if type(head) is ast.BinOp and type(head.op) in (ast.Add, ast.Sub):
                self._terms(head, 1, out)
                rest = args[1:]
            elif _is_sum_call(head):
                out = self._sum_terms(head.args)
                rest = args[1:]
        out += [self._term(a.operand, -1 if type(a.op) is ast.USub else 1)
                for a in rest]
        return out

    def _piecewise(self, args):
        g = self.g
        pairs = []
        for a in args:
            if not isinstance(a, ast.Tuple) or len(a.elts) != 2:
                raise Fallback()
            pairs.append((self._walk(a.elts[0]), g.truth(self._walk(a.elts[1]))))
        if not pairs:
            raise ValueError("Piecewise needs at least one branch")
        # Drop branches after an always-true condition.
        cut = len(pairs)
        for i, (_, c) in enumerate(pairs):
            if c == g.TRUE:
                cut = i + 1
                break
        pairs = pairs[:cut]
        if pairs[-1][1] != g.TRUE:
            # SymPy turns exhaustive conditions into a default branch.
            raise Fallback()
        out = pairs[-1][0]
        for v, c in reversed(pairs[:-1]):
            out = g.select(c, v, out)
        return out


# ===========================================================================
# Differentiation
# ===========================================================================

_NO_DERIV = frozenset((NUM, NAMED, BOOL, LEAF, CMP, AND,
                       OR, NOT))

# Functions with zero derivative everywhere they are defined.
_FLAT = frozenset(("floor", "ceiling", "delta"))


class AD:
    """Differentiation bound to one graph; local partials are cached."""

    def __init__(self, g):
        self.g = g
        self._partials = {}

    # -- local derivatives -------------------------------------------------

    def partials(self, n):
        """[(child position, d n / d child)] of node n."""
        hit = self._partials.get(n)
        if hit is None:
            hit = self._partials[n] = self._local(n)
        return hit

    def _local(self, n):
        g = self.g
        o = g.op[n]
        a = g.args[n]
        if o in _NO_DERIV:
            return []
        if o == ADD:
            return [(i, g.num(c)) for i, c in enumerate(g.attr[n][1])]
        if o == MUL:
            m = len(a)
            if m < 4:
                return [(i, g.prod([a[j] for j in range(m) if j != i]))
                        for i in range(m)]
            # prefix and suffix products, shared: linear, not quadratic, in m
            pre = [g.ONE] * m
            suf = [g.ONE] * m
            for i in range(1, m):
                pre[i] = g.group(g.mul(pre[i - 1], a[i - 1]))
            for i in range(m - 2, -1, -1):
                suf[i] = g.group(g.mul(a[i + 1], suf[i + 1]))
            return [(i, g.mul(pre[i], suf[i])) for i in range(m)]
        if o == POW:
            b, e = a
            out = []
            if g.op[e] == NUM:
                v = g.attr[e]
                out.append((0, g.scale(v, g.pow(b, _addnum(v, Fraction(-1))
                                                 if v != "nan" else float("nan")))))
            else:
                out.append((0, g.mul(e, g.pow(b, g.sub(e, g.ONE)))))
                out.append((1, g.mul(g.call("log", b), n)))
            return out
        if o == CALL:
            return self._call_partials(n)
        if o == SELECT:
            c = a[0]
            return [(1, g.select(c, g.ONE, g.ZERO)),
                    (2, g.select(c, g.ZERO, g.ONE))]
        raise GraphError("no derivative for " + OP_NAMES[o])

    def _call_partials(self, n):
        g = self.g
        name = g.attr[n]
        a = g.args[n]
        x = a[0] if a else None
        one = g.ONE
        if name == "_prod":
            return [(0, one)]
        if name == "exp":
            return [(0, n)]
        if name == "log":
            return [(0, g.inv(x))]
        if name == "sin":
            return [(0, g.call("cos", x))]
        if name == "cos":
            return [(0, g.neg(g.call("sin", x)))]
        if name == "tan":
            return [(0, g.add(g.pow(n, Fraction(2)), one))]
        if name == "asin":
            return [(0, g.pow(g.sub(one, g.pow(x, Fraction(2))), Fraction(-1, 2)))]
        if name == "acos":
            return [(0, g.neg(g.pow(g.sub(one, g.pow(x, Fraction(2))),
                                    Fraction(-1, 2))))]
        if name == "atan":
            return [(0, g.inv(g.add(g.pow(x, Fraction(2)), one)))]
        if name == "sinh":
            return [(0, g.call("cosh", x))]
        if name == "cosh":
            return [(0, g.call("sinh", x))]
        if name == "tanh":
            return [(0, g.sub(one, g.pow(n, Fraction(2))))]
        if name == "asinh":
            return [(0, g.pow(g.add(g.pow(x, Fraction(2)), one), Fraction(-1, 2)))]
        if name == "acosh":
            return [(0, g.mul_factors([(g.sub(x, one), Fraction(-1, 2)),
                                       (g.add(x, one), Fraction(-1, 2))]))]
        if name == "atanh":
            return [(0, g.inv(g.sub(one, g.pow(x, Fraction(2)))))]
        if name == "abs":
            return [(0, g.call("sign", x))]
        if name == "sign":
            return [(0, g.scale(2, g.call("delta", x)))]
        if name == "Heaviside":
            return [(0, g.call("delta", x))]
        if name in _FLAT:
            return [(0, g.ZERO)]
        if name == "min":
            y = a[1]
            return [(0, g.call("Heaviside", g.sub(y, x))),
                    (1, g.call("Heaviside", g.sub(x, y)))]
        if name == "max":
            y = a[1]
            return [(0, g.call("Heaviside", g.sub(x, y))),
                    (1, g.call("Heaviside", g.sub(y, x)))]
        if name in ("erf", "erfc"):
            c = 2.0 / (3.141592653589793 ** 0.5)
            d = g.scale(c, g.call("exp", g.neg(g.pow(x, Fraction(2)))))
            return [(0, d if name == "erf" else g.neg(d))]
        if name == "atan2":
            y, xx = a
            r2 = g.add(g.pow(xx, Fraction(2)), g.pow(y, Fraction(2)))
            return [(0, g.div(xx, r2)), (1, g.neg(g.div(y, r2)))]
        return self._template_partials(n)

    def _template_partials(self, n):
        g = self.g
        name = g.attr[n]
        a = g.args[n]
        out = []
        for i in range(len(a)):
            expr, znames = derivative_template(name, len(a), i)
            out.append((i, from_sympy(g, expr, dict(zip(znames, a)))))
        return out

    # -- forward ------------------------------------------------------------

    def jvp(self, roots, seeds):
        """Tangents of `roots` for `seeds`, a dict leaf node -> tangent node."""
        g = self.g
        mask = 0
        for leaf in seeds:
            mask |= g.flags[leaf]
        tan = {}
        zero = g.ZERO
        for n in g.topo(roots):
            if n in seeds:
                tan[n] = seeds[n]
                continue
            o = g.op[n]
            if o in _NO_DERIV or not (g.flags[n] & mask):
                tan[n] = zero
                continue
            a = g.args[n]
            if o == SELECT:
                ta, tb = tan[a[1]], tan[a[2]]
                tan[n] = zero if ta == zero and tb == zero else g.select(a[0], ta, tb)
                continue
            if o == ADD:
                tan[n] = g.add_terms([(c, tan[t]) for c, t in zip(g.attr[n][1], a)
                                      if tan[t] != zero])
                continue
            if all(tan[c] == zero for c in a):
                tan[n] = zero
                continue
            terms = []
            for i, d in self.partials(n):
                tc = tan[a[i]]
                if tc != zero and not g.is_zero(d):
                    # the tangent stays one factor: chained products do not grow
                    terms.append((1, g.mul(d, g.group(tc))))
            tan[n] = g.add_terms(terms)
        return [tan[r] for r in roots]

    def time_derivative(self, roots):
        """d/dt of `roots` with the forcing rates as chain terms."""
        g = self.g
        seeds = {g.time(): g.ONE}
        for n in g.topo(roots):
            if g.op[n] == LEAF and g.attr[n][0] == FORCING:
                seeds[n] = g.frate(g.attr[n][1])
        return self.jvp(roots, seeds)

    def lie(self, roots, rhs, forcings=True):
        """Total time derivative of `roots` along x' = rhs."""
        g = self.g
        seeds = {g.time(): g.ONE}
        for n in g.topo(roots):
            if g.op[n] != LEAF:
                continue
            k = g.attr[n][0]
            if k == STATE:
                seeds[n] = rhs[g.attr[n][1]]
            elif k == FORCING and forcings:
                seeds[n] = g.frate(g.attr[n][1])
        return self.jvp(roots, seeds)

    # -- reverse ------------------------------------------------------------

    def vjp(self, roots, cotangents, mask):
        """Adjoints of the leaves of kinds `mask` (F_* flags).

        Returns:
            dict leaf node -> adjoint node, for leaves with a nonzero adjoint.
        """
        g = self.g
        order = [n for n in g.topo(roots) if g.flags[n] & mask]
        pending = {}
        for r, c in zip(roots, cotangents):
            if not g.is_zero(c) and (g.flags[r] & mask):
                pending.setdefault(r, []).append(c)
        out = {}
        for n in reversed(order):
            terms = pending.pop(n, None)
            if not terms:
                continue
            adj = terms[0] if len(terms) == 1 else g.sum(terms)
            if g.is_zero(adj):
                continue
            o = g.op[n]
            if o == LEAF:
                out[n] = adj
                continue
            a = g.args[n]
            if o == ADD:
                for c, t in zip(g.attr[n][1], a):
                    if g.flags[t] & mask:
                        pending.setdefault(t, []).append(g.scale(c, adj))
                continue
            if o == SELECT:
                cnd = a[0]
                if g.flags[a[1]] & mask:
                    pending.setdefault(a[1], []).append(g.select(cnd, adj, g.ZERO))
                if g.flags[a[2]] & mask:
                    pending.setdefault(a[2], []).append(g.select(cnd, g.ZERO, adj))
                continue
            # the adjoint stays one factor: chained products do not grow
            adj = g.group(adj)
            for i, d in self.partials(n):
                t = a[i]
                if not (g.flags[t] & mask) or g.is_zero(d):
                    continue
                pending.setdefault(t, []).append(g.mul(d, adj))
        return out
