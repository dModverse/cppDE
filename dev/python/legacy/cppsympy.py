"""SymPy glue shared by the cppDE code generators.

Holds the C++ printer, the normalisation of the R and C logical operators, and
the parse-error wrapper. Piecewise is printed as cppde::select(cond, a, b)
rather than as a ternary: a ternary cannot type check when one branch is an
expression-template node and the other a literal, because neither converts to
the other, and cppde::select converts both branches to a common node instead.
"""

import ast
import re
from tokenize import TokenError

import sympy as sp
from sympy.printing.cxx import CXX17CodePrinter


class CppdePrinter(CXX17CodePrinter):
    """CXX17 printer with cppDE's piecewise form and model-symbol slots.

    `symbols` maps a model symbol name to the C++ it stands for (`x[0]`,
    `params[3]`, `(*F[1])(t)`). Printing the slot keeps user names out of the
    generated source, so no model symbol can collide with an identifier the
    generator emits around it, and none reaches SymPy's reserved-word suffix.
    """

    def __init__(self, settings=None, symbols=None):
        super().__init__(settings)
        self._symbols = {} if symbols is None else dict(symbols)

    def _print_Symbol(self, expr):
        slot = self._symbols.get(expr.name)
        if slot is None:
            return super()._print_Symbol(expr)
        return slot

    def _print_Piecewise(self, expr):
        if expr.args[-1].cond != True:  # noqa: E712  (sympy Boolean, not bool)
            raise ValueError(
                "Piecewise needs an (expr, True) default branch: without one "
                "the generated expression has no value for some inputs.")
        code = self._print(expr.args[-1].expr)
        for value, cond in reversed(expr.args[:-1]):
            code = "cppde::select({}, {}, {})".format(
                self._print(cond), self._print(value), code)
        return code


def is_boolean(expr):
    """True for a relation or a logical combination of them.

    CSE can lift a piecewise condition into its own temp, which then has to be
    declared bool rather than the model's numeric type.
    """
    return isinstance(expr, sp.logic.boolalg.Boolean)


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
    """Rewrite `&&`, `||` and `!` into a form SymPy's parser accepts.

    They bind below the comparisons while Python's `&` and `|` bind above, so
    a textual swap would regroup `a > b && c > d`. Python's own parser settles
    the grouping and the tree comes back as And / Or / Not calls.
    """
    if not _LOGIC_RE.search(expr_str):
        return expr_str
    # The word forms carry a leading space, which would read as an indent.
    src = _LOGIC_RE.sub(lambda m: _LOGIC_WORDS[m.group(0)], expr_str).strip()
    tree = _BoolOpsToCalls().visit(ast.parse(src, mode='eval'))
    return ast.unparse(ast.fix_missing_locations(tree))


def parse_error(expr_str, exc, label=None):
    """The exception to raise when parse_expr rejects an expression.

    The message stays short and on one line, and the chained cause is dropped
    at the raise: reticulate truncates a long message and then indexes it with
    an offset from the full one, which reaches R as a std::out_of_range.
    """
    what = "expression" if label is None else "expression '{}'".format(label)
    flat = " ".join(str(expr_str).split())
    if len(flat) > _MAX_EXPR_CHARS:
        flat = "{}... ({} characters)".format(flat[:_MAX_EXPR_CHARS], len(flat))
    reason = " ".join(str(exc).split())
    return ValueError("cannot parse {}: {} [{}]".format(what, flat, reason))
