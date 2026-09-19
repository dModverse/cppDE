"""SymPy helpers of the cppDE code generators.

Parser fallback (`safe_sympify`, `from_sympy`), derivative templates of rare
functions, `normalise_logic` and `parse_error`.
"""

import ast
import keyword
import re
from functools import lru_cache
from tokenize import TokenError

import sympy as sp
from sympy.parsing.sympy_parser import (convert_xor, parse_expr,
                                        standard_transformations)


_MAX_EXPR_CHARS = 200

_LOGIC_RE = re.compile(r'&&|\|\||!(?!=)')
_LOGIC_WORDS = {'&&': ' and ', '||': ' or ', '!': ' not '}


class _BoolOpsToCalls(ast.NodeTransformer):
    """Rewrite the boolean operators as And / Or / Not calls."""

    def visit_BoolOp(self, node):
        self.generic_visit(node)
        name = 'And' if isinstance(node.op, ast.And) else 'Or'
        return ast.Call(func=ast.Name(id=name, ctx=ast.Load()),
                        args=node.values, keywords=[])

    def visit_UnaryOp(self, node):
        self.generic_visit(node)
        if not isinstance(node.op, ast.Not):
            return node
        return ast.Call(func=ast.Name(id='Not', ctx=ast.Load()),
                        args=[node.operand], keywords=[])


def normalise_logic(expr_str):
    """Rewrite `&&`, `||` and `!` as And / Or / Not calls, grouped by Python's
    parser; a textual swap to `&` and `|` would bind above the comparisons."""
    if not _LOGIC_RE.search(expr_str):
        return expr_str
    src = _logic_words(expr_str).strip()
    if src == expr_str.strip():
        return expr_str
    tree = _BoolOpsToCalls().visit(ast.parse(src, mode='eval'))
    return ast.unparse(ast.fix_missing_locations(tree))


def _logic_words(s):
    """`&&`, `||` and prefix `!` as and/or/not; a postfix `!` is kept."""
    out = []
    i = 0
    while i < len(s):
        two = s[i:i + 2]
        if two in ('&&', '||'):
            out.append(_LOGIC_WORDS[two])
            i += 2
            continue
        if s[i] == '!' and two != '!=':
            prev = s[:i].rstrip()[-1:]
            if not (prev.isalnum() or prev in ('_', ')', ']', '.')):
                out.append(_LOGIC_WORDS['!'])
                i += 1
                continue
        out.append(s[i])
        i += 1
    return ''.join(out)


def parse_error(expr_str, exc, label=None):
    """The exception to raise when parse_expr rejects an expression: one short
    line, cause dropped at the raise, since reticulate mis-indexes a long
    message and R sees a std::out_of_range."""
    what = "expression" if label is None else "expression '{}'".format(label)
    flat = " ".join(str(expr_str).split())
    if len(flat) > _MAX_EXPR_CHARS:
        flat = "{}... ({} characters)".format(flat[:_MAX_EXPR_CHARS], len(flat))
    reason = " ".join(str(exc).split())
    return ValueError("cannot parse {}: {} [{}]".format(what, flat, reason))


# ---------------------------------------------------------------------------
# Parsing a model string with SymPy
# ---------------------------------------------------------------------------

def sbml_piecewise(*args):
    """SBML `piecewise(v1, c1, ..., otherwise)` as sp.Piecewise."""
    pairs = [(args[i], args[i + 1]) for i in range(0, len(args) - 1, 2)]
    if len(args) % 2:
        pairs.append((args[-1], True))
    return sp.Piecewise(*pairs)


@lru_cache(maxsize=1)
def parse_dict():
    """Local dictionary for `parse_expr`.

    S, I, N, O, Q and C are plain symbols; exp10/exp2 are exp(x*log(b)).
    """
    return {
        'S': sp.Symbol('S'), 'I': sp.Symbol('I'), 'N': sp.Symbol('N'),
        'O': sp.Symbol('O'), 'Q': sp.Symbol('Q'), 'C': sp.Symbol('C'),
        'exp': sp.exp,
        'exp10': lambda x: sp.exp(x * sp.log(10)),
        'exp2': lambda x: sp.exp(x * sp.log(2)),
        'log': sp.log, 'ln': sp.log,
        'log10': lambda x: sp.log(x, 10),
        'log2': lambda x: sp.log(x, 2),
        'sin': sp.sin, 'cos': sp.cos, 'tan': sp.tan,
        'cot': sp.cot, 'sec': sp.sec, 'csc': sp.csc,
        'asin': sp.asin, 'acos': sp.acos, 'atan': sp.atan,
        'acot': sp.acot, 'asec': sp.asec, 'acsc': sp.acsc,
        'atan2': sp.atan2,
        'sinh': sp.sinh, 'cosh': sp.cosh, 'tanh': sp.tanh,
        'coth': sp.coth, 'sech': sp.sech, 'csch': sp.csch,
        'asinh': sp.asinh, 'acosh': sp.acosh, 'atanh': sp.atanh,
        'acoth': sp.acoth, 'asech': sp.asech, 'acsch': sp.acsch,
        'sqrt': sp.sqrt, 'cbrt': sp.cbrt, 'root': sp.root, 'pow': sp.Pow,
        'abs': sp.Abs, 'sign': sp.sign,
        'floor': sp.floor, 'ceiling': sp.ceiling,
        'round': lambda x: sp.floor(x + sp.Rational(1, 2)),
        'min': sp.Min, 'max': sp.Max,
        'factorial': sp.factorial, 'gamma': sp.gamma,
        'loggamma': sp.loggamma, 'digamma': sp.digamma,
        'polygamma': sp.polygamma, 'beta': sp.beta,
        'erf': sp.erf, 'erfc': sp.erfc, 'erfi': sp.erfi,
        'besselj': sp.besselj, 'bessely': sp.bessely,
        'besseli': sp.besseli, 'besselk': sp.besselk,
        'Heaviside': sp.Heaviside, 'DiracDelta': sp.DiracDelta,
        'And': sp.And, 'Or': sp.Or, 'Not': sp.Not,
        'KroneckerDelta': sp.KroneckerDelta, 'Piecewise': sp.Piecewise,
        'piecewise': sbml_piecewise,
        'pi': sp.pi, 'E': sp.E, 'oo': sp.oo, 'euler_gamma': sp.EulerGamma,
        're': sp.re, 'im': sp.im, 'conjugate': sp.conjugate, 'arg': sp.arg,
    }


IDENT_RE = re.compile(r'(?<![\.\w])[A-Za-z_][A-Za-z0-9_]*')
PY_RESERVED = frozenset(keyword.kwlist) | {'True', 'False', 'None'}
# Parse-dict names kept when used bare.
_CONSTANTS = frozenset(('pi', 'E', 'oo', 'euler_gamma', 'S', 'I', 'N', 'O',
                        'Q', 'C'))


def safe_sympify(expr_str, local_symbols=None, label=None):
    """Parse `expr_str` with SymPy; bare identifiers become real symbols.

    Raises:
        ValueError: from `parse_error` on a syntax error.
    """
    expr_str = normalise_logic(str(expr_str).strip())
    if expr_str == "0":
        return sp.Integer(0)
    local = dict(parse_dict())
    if local_symbols:
        local.update(local_symbols)
    for m in IDENT_RE.finditer(expr_str):
        name = m.group(0)
        if name in PY_RESERVED or (local_symbols and name in local_symbols):
            continue
        called = expr_str[m.end():].lstrip().startswith('(')
        if name not in local or (not called and name not in _CONSTANTS):
            local[name] = sp.Symbol(name, real=True)
    try:
        return parse_expr(expr_str, local_dict=local,
                          transformations=standard_transformations + (convert_xor,),
                          evaluate=True)
    except (SyntaxError, TokenError, TypeError) as e:
        raise parse_error(expr_str, e, label=label) from None


# ---------------------------------------------------------------------------
# SymPy to the expression graph
# ---------------------------------------------------------------------------

# SymPy function name -> graph function name, for the ones the graph knows.
_SYMPY_CALLS = {
    'exp': 'exp', 'log': 'log', 'sin': 'sin', 'cos': 'cos', 'tan': 'tan',
    'asin': 'asin', 'acos': 'acos', 'atan': 'atan', 'sinh': 'sinh',
    'cosh': 'cosh', 'tanh': 'tanh', 'asinh': 'asinh', 'acosh': 'acosh',
    'atanh': 'atanh', 'Abs': 'abs', 'sign': 'sign', 'floor': 'floor',
    'ceiling': 'ceiling', 'Heaviside': 'Heaviside', 'erf': 'erf',
    'erfc': 'erfc', 'gamma': 'gamma', 'loggamma': 'loggamma',
    'atan2': 'atan2', 'Min': 'min', 'Max': 'max',
}

_SYMPY_REL = {'StrictLessThan': '<', 'LessThan': '<=',
              'StrictGreaterThan': '>', 'GreaterThan': '>=',
              'Equality': '==', 'Unequality': '!='}


def from_sympy(g, expr, symbols):
    """Node for SymPy expression `expr` in graph `g`.

    Args:
        symbols: dict symbol name -> node id.

    Unknown functions become CALL nodes under their SymPy name.
    """
    from fractions import Fraction

    memo = {}

    def conv(e):
        hit = memo.get(e)
        if hit is not None:
            return hit
        r = _conv(e)
        memo[e] = r
        return r

    def _conv(e):
        if e is sp.true or e is sp.S.true:
            return g.TRUE
        if e is sp.false or e is sp.S.false:
            return g.FALSE
        if isinstance(e, sp.Symbol):
            node = symbols.get(e.name)
            if node is None:
                raise ValueError("unknown symbol '%s'" % e.name)
            return node
        if isinstance(e, sp.Integer):
            return g.num(Fraction(int(e)))
        if isinstance(e, sp.Rational):
            return g.num(Fraction(int(e.p), int(e.q)))
        if isinstance(e, sp.Float):
            return g.num(float(e))
        if e is sp.pi:
            return g.num(float(sp.pi))
        if e is sp.E:
            return g.num(float(sp.E))
        if e is sp.EulerGamma:
            return g.num(float(sp.EulerGamma))
        if e is sp.oo:
            return g.num(float('inf'))
        if e is sp.S.NegativeInfinity:
            return g.num(float('-inf'))
        if e is sp.nan:
            return g.num(float('nan'))
        if isinstance(e, sp.Add):
            return g.sum([conv(a) for a in e.args])
        if isinstance(e, sp.Mul):
            return g.prod([conv(a) for a in e.args])
        if isinstance(e, sp.Pow):
            b, x = e.args
            if isinstance(x, sp.Rational):
                return g.pow(conv(b), Fraction(int(x.p), int(x.q)))
            return g.pow(conv(b), conv(x))
        if isinstance(e, sp.Piecewise):
            out = None
            for val, cond in reversed(e.args):
                if out is None:
                    if cond is not sp.true:
                        raise ValueError(
                            "Piecewise needs an (expr, True) default branch: "
                            "without one the generated expression has no "
                            "value for some inputs.")
                    out = conv(val)
                else:
                    out = g.select(conv(cond), conv(val), out)
            return out
        cls = type(e).__name__
        if cls in _SYMPY_REL:
            return g.cmp(_SYMPY_REL[cls], conv(e.args[0]), conv(e.args[1]))
        if isinstance(e, sp.And):
            return g.and_(*[conv(a) for a in e.args])
        if isinstance(e, sp.Or):
            return g.or_(*[conv(a) for a in e.args])
        if isinstance(e, sp.Not):
            return g.not_(conv(e.args[0]))
        if isinstance(e, sp.DiracDelta):
            if len(e.args) > 1 and e.args[1] != 0:
                return g.ZERO
            return g.call('delta', conv(e.args[0]))
        if isinstance(e, sp.Heaviside):
            x = conv(e.args[0])
            if len(e.args) > 1 and e.args[1] != sp.Rational(1, 2):
                h0 = conv(e.args[1])
                return g.select(g.cmp('==', x, g.ZERO), h0, g.call('Heaviside', x))
            return g.call('Heaviside', x)
        if isinstance(e, sp.factorial):
            return g.call('gamma', g.add(conv(e.args[0]), g.ONE))
        if isinstance(e, sp.Function):
            name = _SYMPY_CALLS.get(cls)
            args = [conv(a) for a in e.args]
            if name is not None:
                if name == 'log' and len(args) == 2:
                    return g.div(g.call('log', args[0]), g.call('log', args[1]))
                return g.call(name, *args)
            return g.call(cls, *args)
        raise ValueError("cannot convert SymPy expression %s (%s)" % (e, cls))

    return conv(expr)


@lru_cache(maxsize=256)
def derivative_template(name, nargs, pos):
    """(d/dz_pos name(z_0, ..., z_{nargs-1}), symbol names), cached."""
    zs = sp.symbols('_z0:%d' % nargs, real=True)
    fn = getattr(sp, name, None)
    if fn is None:
        fn = sp.Function(name)
    d = sp.diff(fn(*zs), zs[pos])
    return d, tuple(z.name for z in zs)
