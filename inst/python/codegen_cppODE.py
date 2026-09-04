"""
ODE and Jacobian C++ code generator for cppDE
==============================================

This module provides entry points for generating C++ code:

- generate_ode_cpp(...)
    Generates C++ code for the ODE right-hand side and its Jacobian,
    using SymPy for parsing and symbolic differentiation.

- generate_event_code(...)
    Generates C++ code for fixed-time and root-triggered events.
    For root events, analytical partial derivatives dg/dx and dg/dt
    of the root function g(x, t) are derived symbolically via SymPy
    and emitted as C++ lambdas.  These are used at runtime for the
    IFT-based saltation correction (see cppde_integrate_times.hpp):

      g_dot = sum_i(dg/dx_i * f_i) + dg/dt      (total time derivative of g)
      dt*   = -g / g_dot                          (event timing residual, AD quotient rule)
      x_star  = x + f * dt*                       (shift to event surface)
      x_after = event_map(x_star)                 (apply event action)
      x_final = x_after - f_after * dt*           (shift back to grid time)

    SFINAE dispatch at compile time:
      - double:  plain event action, no saltation needed
      - cppde::dual<double, N> (AD):  analytical saltation, correct first-order
        sensitivities
      - cppde::dual2nd<double, N> (AD2): analytical saltation, correct
        second-order via the dual quotient rule on dt* = -g / g_dot

- generate_rootfunc_code(...)
    Generates C++ code for root function based termination.
    Terminal root events do not require saltation gradients (no state
    modification), so dg_dx / dg_dt are set to nullptr.

Author: Simon Beyer
"""

import re
import keyword
import math
import numbers
import sympy as sp

def _sbml_piecewise(*args):
    """SBML's flat `piecewise(v1, c1, v2, c2, ..., otherwise)` as sp.Piecewise.

    libsbml's L3 formatter emits the branches as one flat argument list. An
    odd argument count means the last entry is the otherwise branch.
    """
    pairs = [(args[i], args[i + 1]) for i in range(0, len(args) - 1, 2)]
    if len(args) % 2:
        pairs.append((args[-1], True))
    return sp.Piecewise(*pairs)


from sympy.parsing.sympy_parser import (
    parse_expr,
    standard_transformations,
    convert_xor,
)
from cppsympy import (CppdePrinter, TokenError, is_boolean, normalise_logic,
                      parse_error)
from functools import lru_cache
from concurrent.futures import ThreadPoolExecutor
import os

# =====================================================================
# Safe parsing configuration
# =====================================================================

@lru_cache(maxsize=1)
def _get_safe_parse_dict():
    """
    Construct a local dictionary for SymPy parsing.
    Cached to avoid repeated dict creation.

    Explicitly overrides SymPy singletons (S, I, N, O, Q, C) that would
    otherwise shadow user symbols with the same name.
    exp10/exp2 are mapped to exp(x*log(10/2)) so that AD codegen
    never emits pow(10, x), which can be problematic with AD types.
    """
    return {
        # Override problematic SymPy singletons
        'S': sp.Symbol('S'),
        'I': sp.Symbol('I'),
        'N': sp.Symbol('N'),
        'O': sp.Symbol('O'),
        'Q': sp.Symbol('Q'),
        'C': sp.Symbol('C'),

        "exp": sp.exp,
        "exp10": lambda x: sp.exp(x * sp.log(10)),
        "exp2": lambda x: sp.exp(x * sp.log(2)),
        "log": sp.log,
        "ln": sp.log,
        "log10": lambda x: sp.log(x, 10),
        "log2": lambda x: sp.log(x, 2),
        "sin": sp.sin, "cos": sp.cos, "tan": sp.tan,
        "cot": sp.cot, "sec": sp.sec, "csc": sp.csc,
        "asin": sp.asin, "acos": sp.acos, "atan": sp.atan,
        "acot": sp.acot, "asec": sp.asec, "acsc": sp.acsc,
        "atan2": sp.atan2,
        "sinh": sp.sinh, "cosh": sp.cosh, "tanh": sp.tanh,
        "coth": sp.coth, "sech": sp.sech, "csch": sp.csch,
        "asinh": sp.asinh, "acosh": sp.acosh, "atanh": sp.atanh,
        "acoth": sp.acoth, "asech": sp.asech, "acsch": sp.acsch,
        "sqrt": sp.sqrt, "cbrt": sp.cbrt, "root": sp.root, "pow": sp.Pow,
        "abs": sp.Abs, "sign": sp.sign,
        "floor": sp.floor, "ceiling": sp.ceiling,
        "min": sp.Min, "max": sp.Max,
        "factorial": sp.factorial, "gamma": sp.gamma,
        "loggamma": sp.loggamma, "digamma": sp.digamma,
        "erf": sp.erf, "erfc": sp.erfc,
        "besselj": sp.besselj, "bessely": sp.bessely,
        "besseli": sp.besseli, "besselk": sp.besselk,
        "Heaviside": sp.Heaviside, "DiracDelta": sp.DiracDelta,
        "And": sp.And, "Or": sp.Or, "Not": sp.Not,
        "Piecewise": sp.Piecewise, "piecewise": _sbml_piecewise,
        "pi": sp.pi, "E": sp.E, "oo": sp.oo,
    }
_IDENT_RE = re.compile(r'(?<![\.\w])[A-Za-z_][A-Za-z0-9_]*')
_PY_RESERVED = frozenset(keyword.kwlist) | {'True', 'False', 'None'}

def _safe_sympify(expr_str, local_symbols=None):
    """Safely parse a string expression to SymPy."""
    expr_str = normalise_logic(str(expr_str).strip())
    if expr_str == "0":
        return sp.Integer(0)

    safe_local = dict(_get_safe_parse_dict())
    if local_symbols:
        safe_local.update(local_symbols)

    # Pre-declare any bare identifier as a Symbol so it shadows SymPy globals
    # like sp.beta / sp.zeta (FunctionClass) that would otherwise leak in via
    # parse_expr's default global_dict=sympy.__dict__ and break "10^beta".
    for name in _IDENT_RE.findall(expr_str):
        if name not in safe_local and name not in _PY_RESERVED:
            safe_local[name] = sp.Symbol(name, real=True)

    transformations = standard_transformations + (convert_xor,)
    try:
        return parse_expr(
            expr_str,
            local_dict=safe_local,
            transformations=transformations,
            evaluate=True,
        )
    except (SyntaxError, TokenError) as e:
        raise parse_error(expr_str, e) from None
def _ensure_double_literals(cpp_code):
    """Convert integer literals to double literals in C++ code."""
    sci_pattern = r'(\d+\.?\d*[eE][+-]?\d+)'
    
    sci_numbers = []
    def store_sci(match):
        sci_numbers.append(match.group(0))
        return f'__SCI_PLACEHOLDER_{len(sci_numbers)-1}__'
    
    temp = re.sub(sci_pattern, store_sci, cpp_code)
    int_pattern = r'(?<![a-zA-Z0-9_.\[])(\d+)(?![0-9.\]])'
    temp = re.sub(int_pattern, lambda m: m.group(1) + '.0', temp)
    
    for i, sci in enumerate(sci_numbers):
        temp = temp.replace(f'__SCI_PLACEHOLDER_{i}__', sci)
    
    return temp
# =====================================================================
# Cached symbol-to-slot mapping for _to_cpp
# =====================================================================

@lru_cache(maxsize=16)
def _get_printer(symbols):
    """Printer for one symbol-to-slot mapping, given as sorted name/slot pairs."""
    return CppdePrinter(symbols=dict(symbols))


class _ModelSymbols:
    """
    The C++ slot every state, param, forcing and time symbol prints as.

    The printer substitutes them, so no model symbol reaches the generated
    source: rewriting the finished source instead turns std::pow into
    params[3]::pow for a parameter named std.
    """
    __slots__ = ('printer', '_map')

    def __init__(self, states, params, n_states, forcings, use_initial_states):
        mapping = {}
        # forcings first (they may shadow param/state names)
        for i, f in enumerate(forcings):
            mapping[f] = f"(*F[{i}])(t)"
        # params
        for j, p in enumerate(params):
            mapping[p] = f"params[{n_states + j}]"
        # states
        if use_initial_states:
            for i, s in enumerate(states):
                mapping[s] = f"params[{i}]"
        else:
            for i, s in enumerate(states):
                mapping[s] = f"x[{i}]"
        # initial-value notation  state_0 -> params[i]
        for i, s in enumerate(states):
            mapping[f"{s}_0"] = f"params[{i}]"
        # time
        mapping["time"] = "t"

        self._map = mapping
        self.printer = _get_printer(tuple(sorted(mapping.items())))

    def slot(self, name):
        """The C++ a bare symbol prints as; an unmapped name stays as it is."""
        return self._map.get(name, name)
@lru_cache(maxsize=8)
def _get_symbols(states_tuple, params_tuple, n_states, forcings_tuple, use_initial_states):
    """Cached factory: the mapping is rebuilt only when the signature changes."""
    return _ModelSymbols(
        list(states_tuple), list(params_tuple), n_states,
        list(forcings_tuple), use_initial_states,
    )
# =====================================================================
# Math macro replacement map
# =====================================================================

_MATH_MACRO_MAP = {
    "M_E": "std::exp(1.0)",
    "M_LOG2E": "1.0 / std::log(2.0)",
    "M_LOG10E": "1.0 / std::log(10.0)",
    "M_LN2": "std::log(2.0)",
    "M_LN10": "std::log(10.0)",
    "M_PI": "std::acos(-1.0)",
    "M_PI_2": "(std::acos(-1.0) * 0.5)",
    "M_PI_4": "(std::acos(-1.0) * 0.25)",
    "M_1_PI": "(1.0 / std::acos(-1.0))",
    "M_2_PI": "(2.0 / std::acos(-1.0))",
    "M_2_SQRTPI": "(2.0 / std::sqrt(std::acos(-1.0)))",
    "M_SQRT2": "std::sqrt(2.0)",
    "M_SQRT1_2": "std::sqrt(0.5)",
}

# Precompiled regex for math macro replacement (single-pass)
_MATH_MACRO_PATTERN = re.compile(
    "|".join(re.escape(k) for k in sorted(_MATH_MACRO_MAP.keys(), key=len, reverse=True))
)

# Precompiled regex for std:: -> cppde:: math-function replacement
# (single-pass). The cppde namespace provides AD-aware overloads for the
# in-tree dual / dual2nd types via cppde_dual_math.hpp.
_AD_FN_PATTERN = re.compile(
    r'\bstd::(sin|cos|tan|asin|acos|atan|sinh|cosh|tanh|'
    r'asinh|acosh|atanh|exp|log|sqrt|pow|abs|min|max)\b'
)

_AD_PREFIX = "cppde"

# Precompiled whitespace collapse pattern
_WHITESPACE_PATTERN = re.compile(r"\s+")

# Precompiled regex for std::pow(expr, 2.0) to (expr)*(expr) optimization.
# Matches std::pow(ARG, 2.0) or std::pow(ARG, 2) where ARG may contain
# nested parentheses (up to 2 levels) or simple identifiers.
def _optimize_pow2(cpp_str):
    """Replace std::pow(expr, 2.0) with (expr)*(expr) for performance.
    Also handles cppde::pow for AD types.

    Uses parenthesis-counting instead of regex to avoid catastrophic
    backtracking on deeply nested expressions.
    """
    for prefix in ("std::pow(", "cppde::pow("):
        result = []
        i = 0
        plen = len(prefix)
        while i < len(cpp_str):
            if cpp_str[i:i+plen] == prefix:
                # Found pow(: find the comma separating args by counting parens
                depth = 1
                start = i + plen
                j = start
                comma_pos = -1
                while j < len(cpp_str) and depth > 0:
                    c = cpp_str[j]
                    if c == '(':
                        depth += 1
                    elif c == ')':
                        depth -= 1
                        if depth == 0:
                            break
                    elif c == ',' and depth == 1:
                        comma_pos = j
                    j += 1
                if comma_pos > 0 and depth == 0:
                    arg1 = cpp_str[start:comma_pos]
                    arg2 = cpp_str[comma_pos+1:j].strip()
                    if arg2 in ("2", "2.0"):
                        result.append(f"(({arg1})*({arg1}))")
                        i = j + 1
                        continue
                # No match or not pow2: emit original
                result.append(prefix)
                i += plen
            else:
                result.append(cpp_str[i])
                i += 1
        cpp_str = "".join(result)
    return cpp_str

def _negate_cpp_expr(expr_str):
    """Negate a C++ expression string for pre-negated Jacobian storage.
    Produces -(expr) for complex expressions, handles simple cases directly."""
    stripped = expr_str.strip()
    if stripped == '0.0' or stripped == '0':
        return stripped
    # -(expr) -> (expr), only if the parenthesis after the sign is the one
    # closed at the end.  In `-(a)/(b) + c` it is not, and only the first
    # term would flip.
    if stripped.startswith('-(') and stripped.endswith(')'):
        depth = 0
        close = -1
        for i, ch in enumerate(stripped[1:], start=1):
            if ch == '(':
                depth += 1
            elif ch == ')':
                depth -= 1
                if depth == 0:
                    close = i
                    break
        if close == len(stripped) - 1:
            return stripped[1:]  # -(expr) -> (expr)
    if stripped.startswith('-') and not stripped.startswith('-('):
        # -simple_expr -> simple_expr
        rest = stripped[1:]
        if re.match(r'^[a-zA-Z0-9_\[\]\.\*]+$', rest):
            return rest
    return f'-({stripped})'

# =====================================================================
# _to_cpp - Convert SymPy expression to C++ code
# =====================================================================

def _to_cpp(expr, states, params, n_states, num_type, forcings=None, use_initial_states=False):
    """Convert a SymPy expression to C++ code.

    The printer substitutes a slot for every model symbol, so the result holds
    no user name. Printer and slot map are cached per symbol set.
    """
    if forcings is None:
        forcings = []

    # Fast path: trivial expressions (skip expensive printer + regex chain)
    if expr is sp.S.Zero or expr == 0:
        return "0.0"
    if isinstance(expr, sp.Integer):
        return str(int(expr)) + ".0"
    if isinstance(expr, sp.Float):
        return repr(float(expr))

    symbols = _get_symbols(
        tuple(states), tuple(params), n_states,
        tuple(forcings), use_initial_states,
    )

    if isinstance(expr, sp.Symbol):
        return symbols.slot(expr.name)
    # Negative symbol: -x  (sp.Mul(-1, x))
    if (isinstance(expr, sp.Mul) and len(expr.args) == 2
            and expr.args[0] is sp.S.NegativeOne
            and isinstance(expr.args[1], sp.Symbol)):
        return "-" + symbols.slot(expr.args[1].name)

    # Full path: complex expressions via the slot-substituting printer
    cpp_code = symbols.printer.doprint(expr).replace("\n", " ")
    
    # Single-pass math macro replacement (precompiled regex)
    cpp_code = _MATH_MACRO_PATTERN.sub(lambda m: _MATH_MACRO_MAP[m.group(0)], cpp_code)
    
    # Single-pass std:: -> {prefix}:: replacement for AD types (precompiled regex)
    if num_type in ("AD", "AD2"):
        cpp_code = _AD_FN_PATTERN.sub(lambda m: f'{_AD_PREFIX}::{m.group(1)}', cpp_code)
    
    cpp_code = _WHITESPACE_PATTERN.sub("", cpp_code)
    cpp_code = _ensure_double_literals(cpp_code)
    cpp_code = _optimize_pow2(cpp_code)
    
    return cpp_code
# =====================================================================
# Main ODE generator
# =====================================================================

def generate_ode_cpp(
    rhs_dict,
    params_list,
    num_type="AD",
    fixed_states=None,
    fixed_params=None,
    forcings_list=None,
    sparse=None,
    skip_jacobian=False,
):
    """
    Generate C++ code for ODE system and Jacobian.

    Parameters:
    -----------
    rhs_dict : dict
        Dictionary mapping state names to RHS expressions
    params_list : list of str
        Parameter names
    num_type : str
        Numeric type ("AD", "AD2", "double")
    fixed_states : list of str, optional
        States to exclude from sensitivity
    fixed_params : list of str, optional
        Parameters to exclude from sensitivity
    forcings_list : list of str, optional
        Forcing function names
    
    Returns:
    --------
    dict with keys: "ode_code", "jac_code", "jac_matrix", "time_derivs",
                    "states", "params", "forcings"
    """
    # Normalize inputs
    if fixed_states is None:
        fixed_states = []
    if fixed_params is None:
        fixed_params = []
    if forcings_list is None:
        forcings_list = []
    if params_list is None:
        params_list = []

    if isinstance(params_list, str):
        params_list = [params_list]
    else:
        params_list = list(params_list)
    
    if isinstance(forcings_list, str):
        forcings_list = [forcings_list]
    else:
        forcings_list = list(forcings_list)
    
    if isinstance(fixed_states, str):
        fixed_states = [fixed_states]
    else:
        fixed_states = list(fixed_states)
    
    if isinstance(fixed_params, str):
        fixed_params = [fixed_params]
    else:
        fixed_params = list(fixed_params)

    states_list = list(rhs_dict.keys())
    odes_list = list(rhs_dict.values())
    n_states = len(states_list)

    # Define symbols
    states_syms = {n: sp.Symbol(n, real=True) for n in states_list}
    params_syms = {n: sp.Symbol(n, real=True) for n in params_list}
    forcing_syms = {n: sp.Symbol(n, real=True) for n in forcings_list}
    t = sp.Symbol("time", real=True)

    local_symbols = {}
    local_symbols.update(states_syms)
    local_symbols.update(params_syms)
    local_symbols.update(forcing_syms)
    local_symbols["time"] = t

    # =====================================================================
    # FAST PATH: Template-based structural deduplication
    # =====================================================================

    # Structurally repeated equations (MOL-discretised PDEs) fingerprint to a
    # canonical form, are parsed and differentiated once per distinct form, and
    # expand per instance: O(n_templates) sp.diff calls instead of O(n_states).
    dedup_result = None
    if n_states > 64 and not forcings_list and not skip_jacobian:
        dedup_result = _try_template_dedup(
            odes_list, states_list, params_list, n_states, num_type,
            forcings_list, t
        )

    if skip_jacobian:
        # EXPLICIT METHOD PATH: skip Jacobian entirely.
        # Only parse RHS + generate ODE code; emit a no-op Jacobian stub.
        exprs = [_safe_sympify(expr, local_symbols) for expr in odes_list]
        jac_matrix = None
        time_derivs = [sp.Integer(0)] * n_states

        ode_cpp_lines = _generate_ode_code_plain(
            exprs, states_list, params_list, n_states, num_type, forcings_list
        )
        jac_cpp_lines = _generate_noop_jacobian(n_states, num_type)

    elif dedup_result is not None:
        # Template dedup succeeded: use its results
        ode_cpp_lines = dedup_result["ode_cpp_lines"]
        jac_cpp_lines = dedup_result["jac_cpp_lines"]
        jac_matrix = dedup_result["jac_matrix"]
        time_derivs = dedup_result["time_derivs"]
    else:
        # STANDARD PATH: full symbolic parse + differentiate
        # Parse RHS expressions
        exprs = [_safe_sympify(expr, local_symbols) for expr in odes_list]

        states_syms_list = [states_syms[s] for s in states_list]
        states_syms_set = set(states_syms_list)

        # Compute Jacobian matrix: SPARSITY-AWARE + PARALLEL
        n_workers = min(os.cpu_count() or 4, n_states)
        use_parallel_jac = n_states > 8 and n_workers > 1

        if use_parallel_jac:
            jac_matrix = _compute_ode_jacobian_parallel(
                exprs, states_syms_list, states_syms_set, n_workers
            )
        else:
            jac_matrix = _compute_ode_jacobian_serial(
                exprs, states_syms_list, states_syms_set
            )

        # Compute time derivatives
        time_derivs = [sp.diff(expr, t) for expr in exprs]

        # Generate ODE + Jacobian code (plain, no CSE)
        ode_cpp_lines = _generate_ode_code_plain(
            exprs, states_list, params_list, n_states, num_type, forcings_list
        )
        jac_cpp_lines = _generate_jac_code_plain(
            jac_matrix, time_derivs, exprs, forcing_syms, forcings_list,
            states_list, params_list, n_states, num_type
        )

    # Sparsity analysis: decide dense vs sparse
    if skip_jacobian:
        pre_pattern = set()
    elif dedup_result is not None:
        pre_pattern = set(zip(dedup_result["jac_nnz_rows"], dedup_result["jac_nnz_cols"]))
    else:
        pre_pattern = _extract_jac_sparsity(jac_matrix, n_states)

    jac_nnz = len(pre_pattern)
    n2 = n_states * n_states
    jac_zeros_pct = 100.0 * (1.0 - jac_nnz / n2) if n2 else 0
    sorted_jpat = sorted(pre_pattern)

    # One threshold for both backends and both modes.
    use_sparse = decide_sparse(sparse, n_states, jac_nnz,
                               has_jacobian=not skip_jacobian)

    sparsity_stats = {
        'n': n_states,
        'jac_nnz': jac_nnz,
        'jac_zeros_pct': jac_zeros_pct,
        'jac_pattern': sorted_jpat,
    }

    # Sparse Jacobian stringification: triplet format (rows, cols, exprs).
    _zero = sp.Integer(0)
    _szero = sp.S.Zero
    jac_nnz_rows = []
    jac_nnz_cols = []
    jac_nnz_exprs = []

    states_syms_list = [states_syms[s] for s in states_list]
    states_syms_set = set(states_syms_list)

    if skip_jacobian:
        pass  # no-op Jacobian already generated; skip sparse stringification
    elif dedup_result is not None:
        jac_nnz_rows = dedup_result["jac_nnz_rows"]
        jac_nnz_cols = dedup_result["jac_nnz_cols"]
        jac_nnz_exprs = dedup_result["jac_nnz_exprs"]
    else:
        for i in range(n_states):
            free = exprs[i].free_symbols & states_syms_set
            for j, s in enumerate(states_syms_list):
                if s in free:
                    e = jac_matrix[i][j]
                    if e is not _zero and e is not _szero and e != 0:
                        jac_nnz_rows.append(i)
                        jac_nnz_cols.append(j)
                        jac_nnz_exprs.append(str(e))

    # --- Generate the Jacobian functor ---

    # Exactly one form is emitted, matching use_sparse: raw CSC with direct Ax[]
    # writes, or matrix<T> through J(i,j). With skip_jacobian the stub stands.
    if skip_jacobian:
        jac_code = "\n".join(jac_cpp_lines)
    elif use_sparse:
        jac_re = re.compile(r'^\s*J\((\d+),(\d+)\)\s*=\s*(.+);')
        num = num_type

        jac_nnz_count = len(jac_nnz_rows)

        # Determine which diagonal entries are missing from J's pattern.
        # The W matrix (I/gamma - J) needs all diagonals for the identity term.
        diag_in_pattern = set()
        for r, c in zip(jac_nnz_rows, jac_nnz_cols):
            if r == c:
                diag_in_pattern.add(r)
        missing_diags = sorted(set(range(n_states)) - diag_in_pattern)
        total_nnz = jac_nnz_count + len(missing_diags)

        # =================================================================
        # Sparse Jacobian via raw CSC (csc_matrix)
        # =================================================================

        # Entries are sorted into CSC order and given a linear Ax index. The first call
        # builds Ap/Ai from the static row/col arrays, every call writes W.Ax[k]
        # directly. Under template dedup a loop form uses a precomputed offset table.

        def _fmt_int_array(values, per_line=20):
            lines = []
            for i in range(0, len(values), per_line):
                chunk = values[i:i + per_line]
                lines.append("        " + ",".join(str(v) for v in chunk))
            return ",\n".join(lines)

        # --- Build CSC-sorted (row, col) to Ax index mapping ---
        # Merge real entries + missing diagonal zeros
        all_rows = list(jac_nnz_rows) + missing_diags
        all_cols = list(jac_nnz_cols) + missing_diags
        all_exprs = list(jac_nnz_exprs) + [None] * len(missing_diags)  # None = zero diag

        # Sort into CSC order: primary by col, secondary by row
        csc_order = sorted(range(total_nnz), key=lambda k: (all_cols[k], all_rows[k]))
        sorted_rows = [all_rows[k] for k in csc_order]
        sorted_cols = [all_cols[k] for k in csc_order]
        sorted_exprs = [all_exprs[k] for k in csc_order]

        # Build (row, col) to Ax index lookup for the loop-based path
        rc_to_ax = {}
        for ax_k, (r, c) in enumerate(zip(sorted_rows, sorted_cols)):
            rc_to_ax[(r, c)] = ax_k

        sparse_jac_lines = []
        sparse_jac_lines.append(f"// Sparse Jacobian: raw CSC, {n_states}x{n_states}, {total_nnz} nnz")
        sparse_jac_lines.append(f"// Stores NEGATED Jacobian (-J) for pre-negated W = (1/γh)I - J.")
        sparse_jac_lines.append(f"struct jacobian {{")
        sparse_jac_lines.append(f"  std::vector<{num}> params;")
        sparse_jac_lines.append(f"  std::vector<const cppde::PchipForcing<{num}>*> F;")
        sparse_jac_lines.append(f"")
        sparse_jac_lines.append(f"  jacobian(const std::vector<{num}>& p_,")
        sparse_jac_lines.append(f"           const std::vector<const cppde::PchipForcing<{num}>*>& F_)")
        sparse_jac_lines.append(f"    : params(p_), F(F_) {{}}")
        sparse_jac_lines.append(f"")
        sparse_jac_lines.append(f"  void operator()(const std::vector<{num}>& x,")
        sparse_jac_lines.append(f"                  cppde::csc_matrix<{num}>& W,")
        sparse_jac_lines.append(f"                  const {num}& t,")
        sparse_jac_lines.append(f"                  std::vector<{num}>& dfdt) {{")
        # No per-Jacobian arena scope: W entries and dfdt are not slab-bound, so their
        # tan_ pointers come from the arena, and the LU solver consumes them after this
        # call returns.

        # Gate constant-entry init on W being pattern-built, so a caller passing a fresh
        # csc_matrix gets the constants and later calls with the same W skip them.
        sparse_jac_lines.append(f"    const bool _init_consts = !W.pattern_built;")

        # First-call: build CSC pattern from static arrays
        rows_str = _fmt_int_array(sorted_rows)
        cols_str = _fmt_int_array(sorted_cols)
        sparse_jac_lines.append(f"    if (!W.pattern_built) {{")
        sparse_jac_lines.append(f"      static const int _rows[] = {{")
        sparse_jac_lines.append(rows_str)
        sparse_jac_lines.append(f"      }};")
        sparse_jac_lines.append(f"      static const int _cols[] = {{")
        sparse_jac_lines.append(cols_str)
        sparse_jac_lines.append(f"      }};")
        sparse_jac_lines.append(f"      W.build_pattern({n_states}, {total_nnz}, _rows, _cols);")
        sparse_jac_lines.append(f"    }}")

        use_loop = (dedup_result is not None and 'groups' in dedup_result)

        if use_loop:
            # =============================================================
            # LOOP-BASED: one loop per template group
            # =============================================================

            # The (row, dep[j]) to Ax offset is precomputed into a static table; the loop
            # body writes W.Ax[ax_offsets[c*n_jac+j]].
            groups = dedup_result['groups']
            tmpl_cpp = dedup_result['template_cpp']
            name_to_idx = dedup_result['name_to_state_idx']

            def _tmpl_to_loop_cpp(tmpl_str):
                result = _GENERIC_PATTERN.sub(lambda m: f'x[s[{m.group(1)}]]', tmpl_str)
                result = _WHITESPACE_PATTERN.sub(" ", result).strip()
                result = _ensure_double_literals(result)
                result = _optimize_pow2(result)
                return result

            # Emit shared data tables + Ax offset tables at function scope
            for tmpl_idx, (key, members) in enumerate(groups.items()):
                n_instances = len(members)
                n_deps = len(members[0][1])
                jac_positions = sorted(tmpl_cpp[key]['jac'].keys())
                n_jac = len(jac_positions)

                deps_data = []
                rows_data = []
                ax_offsets_data = []
                for expr_idx, dep_names in members:
                    rows_data.append(expr_idx)
                    dep_indices = [name_to_idx[dep_names[j]] for j in range(n_deps)]
                    deps_data.extend(dep_indices)
                    # Precompute Ax offset for each Jacobian entry of this instance
                    for j in jac_positions:
                        col_idx = dep_indices[j]
                        ax_k = rc_to_ax.get((expr_idx, col_idx), -1)
                        ax_offsets_data.append(ax_k)

                sparse_jac_lines.append(f"    // Template {tmpl_idx}: {n_instances} instances, {n_deps} deps, {n_jac} jac entries")
                sparse_jac_lines.append(f"    static const int t{tmpl_idx}_deps[] = {{")
                sparse_jac_lines.append(_fmt_int_array(deps_data))
                sparse_jac_lines.append(f"    }};")
                sparse_jac_lines.append(f"    static const int t{tmpl_idx}_rows[] = {{")
                sparse_jac_lines.append(_fmt_int_array(rows_data))
                sparse_jac_lines.append(f"    }};")
                sparse_jac_lines.append(f"    static const int t{tmpl_idx}_ax[] = {{")
                sparse_jac_lines.append(_fmt_int_array(ax_offsets_data))
                sparse_jac_lines.append(f"    }};")

            # Hot path: loops with direct W.Ax[] writes. Expressions are negated, because
            # W = (1/gamma*h)*I - J holds -J directly and factorize_W then needs no
            # O(nnz) negate-copy.

            # Constant entries, those not depending on x[], are written once in an init
            # block rather than on every call.
            for tmpl_idx, (key, members) in enumerate(groups.items()):
                n_instances = len(members)
                n_deps = len(members[0][1])
                jac_tmpl = tmpl_cpp[key]['jac']
                td_tmpl = tmpl_cpp[key]['time_deriv']
                jac_positions = sorted(jac_tmpl.keys())
                n_jac = len(jac_positions)

                # Classify entries: constant (no x[s[) vs state-dependent
                const_entries = []  # (local_j, negated_expr)
                state_entries = []  # (local_j, negated_expr)
                for local_j, j in enumerate(jac_positions):
                    loop_expr = _tmpl_to_loop_cpp(jac_tmpl[j])
                    neg_expr = _negate_cpp_expr(loop_expr)
                    if 'x[s[' in loop_expr:
                        state_entries.append((local_j, neg_expr))
                    else:
                        const_entries.append((local_j, neg_expr))

                # Constant entries: write once per fresh W (gated on _init_consts).
                if const_entries:
                    sparse_jac_lines.append(f"    // Template {tmpl_idx}: {len(const_entries)} constant + {len(state_entries)} state-dependent entries")
                    sparse_jac_lines.append(f"    if (_init_consts) {{")
                    sparse_jac_lines.append(f"      for (int c = 0; c < {n_instances}; ++c) {{")
                    for local_j, neg_expr in const_entries:
                        sparse_jac_lines.append(f"        W.Ax[t{tmpl_idx}_ax[c * {n_jac} + {local_j}]] = {neg_expr};")
                    sparse_jac_lines.append(f"      }}")
                    sparse_jac_lines.append(f"    }}")

                # State-dependent entries: always written
                sparse_jac_lines.append(f"    for (int c = 0; c < {n_instances}; ++c) {{")
                sparse_jac_lines.append(f"      const int* s = t{tmpl_idx}_deps + c * {n_deps};")
                for local_j, neg_expr in state_entries:
                    sparse_jac_lines.append(f"      W.Ax[t{tmpl_idx}_ax[c * {n_jac} + {local_j}]] = {neg_expr};")

                loop_dfdt = _tmpl_to_loop_cpp(td_tmpl)
                sparse_jac_lines.append(f"      dfdt[t{tmpl_idx}_rows[c]] = {loop_dfdt};")
                sparse_jac_lines.append(f"    }}")

        else:
            # =============================================================
            # PER-ENTRY: unrolled W.Ax[k] for small or irregular systems
            # =============================================================

            # Build (row, col) to expression mapping from the dense Jacobian code
            dense_jac_entries = {}
            dfdt_lines = []
            in_body = False
            brace_depth = 0
            other_lines = []  # CSE temporaries etc.
            for line in jac_cpp_lines:
                stripped = line.strip()
                if not in_body:
                    if 'dfdt)' in stripped and '{' in stripped:
                        in_body = True
                        brace_depth = 1
                    continue
                brace_depth += stripped.count('{') - stripped.count('}')
                if brace_depth <= 0:
                    break
                if 'set_zero' in stripped or 'zero_matrix' in stripped or '::Zero(' in stripped:
                    continue
                # Skip dense dirty-entry clearing code (static _dr/_dc arrays and J() zeroing loop)
                if 'static const int _dr[' in stripped or 'static const int _dc[' in stripped:
                    continue
                if '_dr[_k]' in stripped or '_dc[_k]' in stripped:
                    continue
                m = jac_re.match(line)
                if m:
                    row_idx, col_idx = int(m.group(1)), int(m.group(2))
                    expr_str = m.group(3)
                    dense_jac_entries[(row_idx, col_idx)] = expr_str
                elif 'dfdt[' in stripped:
                    dfdt_lines.append(line)
                else:
                    other_lines.append(line)

            # Emit CSE temporaries first
            for line in other_lines:
                sparse_jac_lines.append(line)

            # Ax writes in CSC order. dense_jac_entries are already negated by
            # _negate_cpp_expr, so negating again would flip the sign back to +J.
            for ax_k, (r, c, expr) in enumerate(zip(sorted_rows, sorted_cols, sorted_exprs)):
                if expr is not None:
                    # Real entry: look up the already-negated C++ expression
                    cpp_expr = dense_jac_entries.get((r, c), expr)
                    cpp_expr = _optimize_pow2(cpp_expr)
                    sparse_jac_lines.append(f"    W.Ax[{ax_k}] = {cpp_expr};")
                else:
                    # Missing diagonal: zero (will get identity term from factorize_W)
                    sparse_jac_lines.append(f"    W.Ax[{ax_k}] = {num}(0);")

            # Emit dfdt lines
            for line in dfdt_lines:
                sparse_jac_lines.append(line)

        sparse_jac_lines.append(f"  }}")
        sparse_jac_lines.append(f"}};")

        # Replace jac_code with the sparse version
        jac_code = "\n".join(sparse_jac_lines)
        jac_cpp_lines = sparse_jac_lines
    else:
        jac_code = "\n".join(jac_cpp_lines)

    # KLU auto-tuning: analyze pattern at codegen time
    klu_settings = None
    if use_sparse:
        klu_settings = analyze_klu_settings(n_states, jac_nnz_rows, jac_nnz_cols)

    return {
        "ode_code": "\n".join(ode_cpp_lines),
        "jac_code": jac_code,
        "jac_nnz_rows": jac_nnz_rows,
        "jac_nnz_cols": jac_nnz_cols,
        "jac_nnz_exprs": jac_nnz_exprs,
        "time_derivs": time_derivs if dedup_result is not None else [str(d) if d != 0 else "0" for d in time_derivs],
        "states": states_list,
        "params": params_list,
        "forcings": forcings_list,
        "use_sparse": use_sparse,
        "sparsity_stats": sparsity_stats,
        "klu_settings": klu_settings,
    }
# =====================================================================
# Parallel Jacobian computation for generate_ode_cpp
# =====================================================================

def _compute_ode_jac_row(expr, states_syms_list, states_syms_set):
    """Compute one row of the ODE Jacobian (for parallelization).
    Sparsity-aware: skips sp.diff when the state is absent from the expression."""
    free = expr.free_symbols & states_syms_set
    row = []
    for s in states_syms_list:
        if s in free:
            row.append(sp.diff(expr, s))
        else:
            row.append(sp.Integer(0))
    return row
def _compute_ode_jacobian_parallel(exprs, states_syms_list, states_syms_set, n_workers):
    """Compute ODE Jacobian in parallel across rows."""
    with ThreadPoolExecutor(max_workers=n_workers) as executor:
        futures = [
            executor.submit(_compute_ode_jac_row, expr, states_syms_list, states_syms_set)
            for expr in exprs
        ]
        # Preserve row order
        return [f.result() for f in futures]
def _compute_ode_jacobian_serial(exprs, states_syms_list, states_syms_set):
    """Compute ODE Jacobian serially."""
    jac_matrix = []
    for expr in exprs:
        free = expr.free_symbols & states_syms_set
        row = []
        for s in states_syms_list:
            if s in free:
                row.append(sp.diff(expr, s))
            else:
                row.append(sp.Integer(0))
        jac_matrix.append(row)
    return jac_matrix
def _generate_noop_jacobian(n_states, num_type):
    """Generate a no-op Jacobian struct for explicit methods (tsit5, adams).

    The struct satisfies the same interface as the real Jacobian functor
    so the generated C++ compiles, but the body is empty: explicit
    steppers never call the Jacobian."""
    return [
        "// No-op Jacobian (explicit method: Jacobian not needed at runtime)",
        "struct jacobian {",
        f"  std::vector<{num_type}> params;",
        f"  std::vector<const cppde::PchipForcing<{num_type}>*> F;",
        "",
        f"  jacobian(const std::vector<{num_type}>& p_,",
        f"           const std::vector<const cppde::PchipForcing<{num_type}>*>& F_)",
        "    : params(p_), F(F_) {}",
        "",
        f"  void operator()(const std::vector<{num_type}>& x,",
        f"                  cppde::dense_matrix<{num_type}>& J,",
        f"                  const {num_type}& t,",
        f"                  std::vector<{num_type}>& dfdt) {{",
        f"    // Explicit method: Jacobian never evaluated.",
        f"    (void)x; (void)J; (void)t; (void)dfdt;",
        "  }",
        "};",
    ]

# =====================================================================
# ODE and Jacobian code generation
# =====================================================================

# =====================================================================
# Common subexpression elimination
# =====================================================================

# sp.cse lifts a repeated subexpression into a num_type local, so the dual ET
# engine substitutes it by reference instead of rebuilding the tree at every
# occurrence. Skipped below 4 expressions, or when sympy finds nothing shared.

def _cse_temps(exprs, prefix='_cse_t'):
    """Apply sympy CSE; return (temps, simplified_exprs) or ([], exprs)."""
    if len(exprs) < 4:
        return [], list(exprs)
    syms = sp.numbered_symbols(prefix=prefix, cls=sp.Symbol)
    temps, simplified = sp.cse(list(exprs), symbols=syms, optimizations='basic')
    if len(temps) == 0:
        return [], simplified
    return temps, simplified


def _emit_cse_temps(temps, states_list, params_list, n_states, num_type, forcings_list):
    """Emit `const num_type _cse_tN = expr;` lines from sp.cse temps.

    Materialising into num_type (rather than `auto`) forces ET evaluation at
    the temp boundary so subsequent references are scalar dual loads instead
    of re-walks of the ET tree."""
    lines = []
    for sym, sub in temps:
        sub_cpp = _to_cpp(sub, states_list, params_list, n_states, num_type, forcings_list)
        # A piecewise condition can be lifted into its own temp, and that one
        # holds a truth value, not a model quantity.
        temp_type = "bool" if is_boolean(sub) else num_type
        lines.append(f"    const {temp_type} {sym.name} = {sub_cpp};")
    return lines


def _arena_scope_lines(num_type):
    """Per-RHS dual_arena::scope guard for the (non-nested) AD code path.

    Bounds the arena working set during a solve: every dual temporary
    (CSE locals, ET assignment buffers) bumps the thread-local arena;
    the scope rolls back to baseline when the RHS functor returns.

    Only safe for num_type == "AD" (single-level dual<double, N>): in that
    mode dxdt is slab-bound (is_dynamic_dual<dual<double, N>> = true), so
    write-to-slab uses the COPY-into-bound-buffer branch of move-assign,
    not the STEAL branch. Arena rollback then frees only CSE temps.

    NOT safe for num_type == "AD2" (nested dual<dual<double, N>, N>): the
    nested predicate is_dynamic_dual<dual<dual<double,N>,N>> is FALSE,
    so dxdt is NOT slab-bound. dxdt[i].tan_ starts at nullptr; an ET
    assignment from a temporary STEALS the rvalue's arena pointer. After
    scope rollback that pointer dangles, segfaulting the next read.
    Nested-AD therefore relies on the outer solveODE-level arena scope
    only: working-set growth is bounded by total RHS calls × per-RHS
    temps, but no per-call scope is safe.

    Non-AD (num_type == "double"): no arena involvement, no scope needed."""
    if num_type == "AD":
        return ["    cppde::dual_arena::scope _rhs_arena_scope;"]
    return []


def _generate_ode_code_plain(exprs, states_list, params_list, n_states, num_type, forcings_list):
    """Generate ODE system C++ code (with CSE)."""
    ode_cpp_lines = [
        "// ODE system",
        "struct ode_system {",
        f"  std::vector<{num_type}> params;",
        f"  std::vector<const cppde::PchipForcing<{num_type}>*> F;",
        "",
        f"  ode_system(const std::vector<{num_type}>& p_,",
        f"             const std::vector<const cppde::PchipForcing<{num_type}>*>& F_)",
        "    : params(p_), F(F_) {}",
        "",
        f"  void operator()(const std::vector<{num_type}>& x,",
        f"                  std::vector<{num_type}>& dxdt,",
        f"                  const {num_type}& t) {{",
    ]
    ode_cpp_lines += _arena_scope_lines(num_type)
    cse_temps, simplified = _cse_temps(exprs)
    ode_cpp_lines += _emit_cse_temps(
        cse_temps, states_list, params_list, n_states, num_type, forcings_list
    )
    for i, expr in enumerate(simplified):
        cpp = _to_cpp(expr, states_list, params_list, n_states, num_type, forcings_list)
        ode_cpp_lines.append(f"    dxdt[{i}] = {cpp};")
    ode_cpp_lines += ["  }", "};"]
    return ode_cpp_lines
def _generate_jac_code_plain(jac_matrix, time_derivs, exprs, forcing_syms, forcings_list,
                             states_list, params_list, n_states, num_type):
    """Generate Jacobian C++ code.

    Always produces DENSE signature (cppde::dense_matrix<T>& J, J(i,j) = ...).
    The caller in generate_ode_cpp converts to sparse if use_sparse is True.
    """
    _zero = sp.Integer(0)
    _szero = sp.S.Zero

    jac_cpp_lines = [
        "// Jacobian for stiff solver",
        "struct jacobian {",
        f"  std::vector<{num_type}> params;",
        f"  std::vector<const cppde::PchipForcing<{num_type}>*> F;",
        "",
        f"  jacobian(const std::vector<{num_type}>& p_,",
        f"           const std::vector<const cppde::PchipForcing<{num_type}>*>& F_)",
        "    : params(p_), F(F_) {}",
        "",
        f"  void operator()(const std::vector<{num_type}>& x,",
        f"                  cppde::dense_matrix<{num_type}>& J,",
        f"                  const {num_type}& t,",
        f"                  std::vector<{num_type}>& dfdt) {{",
        f"    J.set_zero();",
    ]
    # No per-Jacobian arena scope: J entries and dfdt are not slab-bound, so their
    # tan_ pointers come from the arena and the LU solver consumes them after this
    # call returns. The solveODE-level scope catches the growth at solve end.

    # Collect non-zero (i,j) positions first for dirty-index clearing
    jac_entries_plain = []
    for i in range(n_states):
        for j in range(n_states):
            e = jac_matrix[i][j]
            if e is not _zero and e is not _szero and e != 0:
                jac_entries_plain.append((i, j, e))

    # Replace set_zero with dirty-index clearing
    jac_cpp_lines[-1] = f"    // Clear only dirty entries from previous call (O(nnz) not O(n²))"
    n_dirty = len(jac_entries_plain)
    if n_dirty > 0:
        dirty_rows = [str(i) for i, j, e in jac_entries_plain]
        dirty_cols = [str(j) for i, j, e in jac_entries_plain]
        jac_cpp_lines.append(f"    static const int _dr[{n_dirty}] = {{{','.join(dirty_rows)}}};")
        jac_cpp_lines.append(f"    static const int _dc[{n_dirty}] = {{{','.join(dirty_cols)}}};")
        jac_cpp_lines.append(f"    for (int _k = 0; _k < {n_dirty}; ++_k) J(_dr[_k], _dc[_k]) = {num_type}(0);")

    # CSE across all non-zero Jacobian entries: in MM/Hill kinetics, the same
    # denominator (Km + x_j) often shows up in many df/dx_k entries: CSE lifts
    # those into _cse_t* temps materialised once per call.
    jac_exprs = [e for _, _, e in jac_entries_plain]
    jac_temps, jac_simplified = _cse_temps(jac_exprs, prefix='_cse_jt')
    jac_cpp_lines += _emit_cse_temps(
        jac_temps, states_list, params_list, n_states, num_type, []
    )

    # fill NEGATED entries
    for (i, j, _), e in zip(jac_entries_plain, jac_simplified):
        cpp = _to_cpp(e, states_list, params_list, n_states, num_type, [])
        neg_cpp = _negate_cpp_expr(cpp)
        jac_cpp_lines.append(f"    J({i},{j}) = {neg_cpp};")

    for i, expr in enumerate(exprs):
        cpp_code = _to_cpp(time_derivs[i], states_list, params_list, n_states, num_type, forcings_list)
        forcing_terms = []
        for j, fname in enumerate(forcings_list):
            df_dF = sp.diff(expr, forcing_syms[fname])
            if df_dF != 0:
                forcing_terms.append(
                    f"({_to_cpp(df_dF, states_list, params_list, n_states, num_type, forcings_list)})*F[{j}]->derivative(t)"
                )
        if forcing_terms:
            cpp_code = " + ".join([cpp_code] + forcing_terms) if cpp_code != "0" else " + ".join(forcing_terms)
        jac_cpp_lines.append(f"    dfdt[{i}] = {cpp_code};")
    jac_cpp_lines += ["  }", "};"]
    return jac_cpp_lines
# =====================================================================
# Template-based structural deduplication for large MOL systems
# =====================================================================

# Precompiled regex for generic placeholder substitution
_GENERIC_PATTERN = re.compile(r'_s(\d+)')
def _try_template_dedup(odes_list, states_list, params_list, n_states, num_type,
                         forcings_list, t_sym):
    """
    Attempt template-based structural deduplication for large ODE systems.

    For Method-of-Lines discretized PDEs and similar systems where many
    equations share the same algebraic structure (differing only in which
    state variables appear), this avoids redundant symbolic differentiation.

    Algorithm:
      1. Fingerprint: replace state names in each RHS string with positional
         placeholders (_s0, _s1, ...) to get a canonical form.
      2. Group expressions with identical canonical form.
      3. If dedup ratio (n_states / n_templates) >= 4, proceed.
      4. Parse & differentiate only the unique templates.
      5. Expand to all instances via fast string substitution.

    Returns None if dedup is not worthwhile (< 4× reduction).
    Otherwise returns a dict with all data needed by generate_ode_cpp.
    """
    states_set = set(states_list)

    # --- Step 1: String-based fingerprinting ---
    # Build regex to match any state name (longest first to avoid partial matches)
    sorted_names = sorted(states_list, key=len, reverse=True)
    name_pattern = re.compile(
        r'\b(' + '|'.join(re.escape(n) for n in sorted_names) + r')\b'
    )

    from collections import defaultdict
    groups = defaultdict(list)  # canonical_str -> [(expr_idx, [dep_state_names])]

    for i, expr_str in enumerate(odes_list):
        seen = {}
        def _replace_state(m, _seen=seen):
            sname = m.group(0)
            if sname in states_set:
                if sname not in _seen:
                    _seen[sname] = f'_s{len(_seen)}'
                return _seen[sname]
            return sname

        canonical = name_pattern.sub(_replace_state, expr_str)
        dep_names = sorted(seen.keys(), key=lambda s: seen[s])
        groups[canonical].append((i, dep_names))

    n_templates = len(groups)
    dedup_ratio = n_states / max(n_templates, 1)

    if dedup_ratio < 4:
        return None  # Not enough structural repetition

    # --- Step 2: Parse unique templates with generic symbols ---
    n_generic = max(len(members[0][1]) for members in groups.values())
    generic_syms = {f'_s{j}': sp.Symbol(f'_s{j}', real=True) for j in range(n_generic)}

    generic_local = dict(generic_syms)
    for p in params_list:
        generic_local[p] = sp.Symbol(p, real=True)
    generic_local['time'] = t_sym

    template_data = {}  # canonical -> {expr, jac_entries, time_deriv, n_deps}
    for key, members in groups.items():
        expr = _safe_sympify(key, generic_local)
        n_deps = len(members[0][1])

        # Compute Jacobian entries for this template
        free = expr.free_symbols
        jac_entries = {}  # dep_position -> (derivative_expr, generic_sym_name)
        for j in range(n_deps):
            gs = generic_syms[f'_s{j}']
            if gs in free:
                jac_entries[j] = sp.diff(expr, gs)

        # Time derivative
        time_deriv = sp.diff(expr, t_sym)

        template_data[key] = {
            'expr': expr,
            'jac_entries': jac_entries,
            'time_deriv': time_deriv,
            'n_deps': n_deps,
        }

    # --- Step 3: Convert templates to C++ strings ---
    # Params, states and time print as their slots; only the _sN placeholders
    # stay symbolic until the template is expanded per instance.
    tmpl_map = {n: f'x[{i}]' for i, n in enumerate(states_list)}
    tmpl_map.update({p: f'params[{n_states + i}]' for i, p in enumerate(params_list)})
    tmpl_map['time'] = 't'
    printer = _get_printer(tuple(sorted(tmpl_map.items())))

    def _template_to_cpp(sympy_expr):
        """Convert template expression to intermediate C++ (with _sN placeholders)."""
        if sympy_expr is sp.S.Zero or sympy_expr == 0:
            return "0.0"
        if isinstance(sympy_expr, sp.Integer):
            return str(int(sympy_expr)) + ".0"
        cpp = printer.doprint(sympy_expr).replace("\n", " ")
        cpp = _MATH_MACRO_PATTERN.sub(lambda m: _MATH_MACRO_MAP[m.group(0)], cpp)
        if num_type in ("AD", "AD2"):
            cpp = _AD_FN_PATTERN.sub(lambda m: f'{_AD_PREFIX}::{m.group(1)}', cpp)
        return cpp

    template_cpp = {}  # canonical -> {rhs_cpp, jac_cpp: {j: str}, time_deriv_cpp}
    for key, tdata in template_data.items():
        rhs_cpp = _template_to_cpp(tdata['expr'])
        jac_cpp = {j: _template_to_cpp(d) for j, d in tdata['jac_entries'].items()}
        td_cpp = _template_to_cpp(tdata['time_deriv'])
        template_cpp[key] = {'rhs': rhs_cpp, 'jac': jac_cpp, 'time_deriv': td_cpp}

    # --- Step 4: Expand the templates to their instances ---
    name_to_state_idx = {n: i for i, n in enumerate(states_list)}

    def _expand_template(template_str, dep_names):
        """Expand a template C++ string by substituting _sN by the concrete states."""
        cpp = _GENERIC_PATTERN.sub(
            lambda m: f'x[{name_to_state_idx[dep_names[int(m.group(1))]]}]', template_str)
        cpp = _WHITESPACE_PATTERN.sub("", cpp)
        cpp = _ensure_double_literals(cpp)
        return cpp

    # --- Step 5: Generate ODE code ---
    ode_cpp_lines = [
        "// ODE system",
        "struct ode_system {",
        f"  std::vector<{num_type}> params;",
        f"  std::vector<const cppde::PchipForcing<{num_type}>*> F;",
        "",
        f"  ode_system(const std::vector<{num_type}>& p_,",
        f"             const std::vector<const cppde::PchipForcing<{num_type}>*>& F_)",
        "    : params(p_), F(F_) {}",
        "",
        f"  void operator()(const std::vector<{num_type}>& x,",
        f"                  std::vector<{num_type}>& dxdt,",
        f"                  const {num_type}& t) {{",
    ]
    ode_cpp_lines += _arena_scope_lines(num_type)
    for key, members in groups.items():
        rhs_tmpl = template_cpp[key]['rhs']
        for expr_idx, dep_names in members:
            cpp = _expand_template(rhs_tmpl, dep_names)
            ode_cpp_lines.append(f"    dxdt[{expr_idx}] = {cpp};")
    ode_cpp_lines += ["  }", "};"]

    # --- Step 6: Generate Jacobian code ---
    jac_cpp_lines = [
        "// Jacobian for stiff solver",
        "struct jacobian {",
        f"  std::vector<{num_type}> params;",
        f"  std::vector<const cppde::PchipForcing<{num_type}>*> F;",
        "",
        f"  jacobian(const std::vector<{num_type}>& p_,",
        f"           const std::vector<const cppde::PchipForcing<{num_type}>*>& F_)",
        "    : params(p_), F(F_) {}",
        "",
        f"  void operator()(const std::vector<{num_type}>& x,",
        f"                  cppde::dense_matrix<{num_type}>& J,",
        f"                  const {num_type}& t,",
        f"                  std::vector<{num_type}>& dfdt) {{",
    ]
    # No per-Jacobian arena scope (see _generate_jac_code_plain comment).

    jac_nnz_rows = []
    jac_nnz_cols = []
    jac_nnz_exprs = []

    # Build symbolic jac_matrix for sparsity pattern code generation
    # Store symbolic entry for sparsity analysis
    _zero_sym = sp.Integer(0)
    jac_matrix = [[_zero_sym] * n_states for _ in range(n_states)]

    # Collect all (row, col, cpp_expr) entries first
    jac_entry_list = []
    for key, members in groups.items():
        jac_tmpl = template_cpp[key]['jac']
        jac_sym = template_data[key]['jac_entries']
        for expr_idx, dep_names in members:
            for j, jac_str in jac_tmpl.items():
                actual_state = dep_names[j]
                col_idx = name_to_state_idx[actual_state]
                cpp = _expand_template(jac_str, dep_names)
                jac_entry_list.append((expr_idx, col_idx, cpp))

                jac_nnz_rows.append(expr_idx)
                jac_nnz_cols.append(col_idx)
                jac_nnz_exprs.append(str(jac_sym[j]))

                # Store symbolic entry for sparsity analysis
                jac_matrix[expr_idx][col_idx] = jac_sym[j]

    # Emit dirty-index clearing (O(nnz) not O(n²))
    n_dirty = len(jac_entry_list)
    jac_cpp_lines.append(f"    // Clear only dirty entries from previous call (O(nnz) not O(n²))")
    if n_dirty > 0:
        dirty_rows = [str(r) for r, c, e in jac_entry_list]
        dirty_cols = [str(c) for r, c, e in jac_entry_list]
        jac_cpp_lines.append(f"    static const int _dr[{n_dirty}] = {{{','.join(dirty_rows)}}};")
        jac_cpp_lines.append(f"    static const int _dc[{n_dirty}] = {{{','.join(dirty_cols)}}};")
        jac_cpp_lines.append(f"    for (int _k = 0; _k < {n_dirty}; ++_k) J(_dr[_k], _dc[_k]) = {num_type}(0);")

    # Emit NEGATED Jacobian assignments
    for expr_idx, col_idx, cpp in jac_entry_list:
        neg_cpp = _negate_cpp_expr(cpp)
        jac_cpp_lines.append(f"    J({expr_idx},{col_idx}) = {neg_cpp};")

    # dfdt
    for key, members in groups.items():
        td_tmpl = template_cpp[key]['time_deriv']
        for expr_idx, dep_names in members:
            cpp = _expand_template(td_tmpl, dep_names)
            jac_cpp_lines.append(f"    dfdt[{expr_idx}] = {cpp};")

    jac_cpp_lines += ["  }", "};"]

    # Time derivatives as strings (for return value)
    time_derivs = []
    for key, members in groups.items():
        td = template_data[key]['time_deriv']
        for expr_idx, dep_names in members:
            time_derivs.append((expr_idx, td))
    time_derivs.sort(key=lambda x: x[0])
    time_derivs_ordered = [str(td) if td != 0 else "0" for _, td in time_derivs]

    return {
        "ode_cpp_lines": ode_cpp_lines,
        "jac_cpp_lines": jac_cpp_lines,
        "jac_matrix": jac_matrix,
        "time_derivs": time_derivs_ordered,
        "jac_nnz_rows": jac_nnz_rows,
        "jac_nnz_cols": jac_nnz_cols,
        "jac_nnz_exprs": jac_nnz_exprs,
        # Extra data for loop-based sparse codegen
        "groups": dict(groups),              # canonical -> [(expr_idx, dep_names)]
        "template_cpp": template_cpp,        # canonical -> {rhs, jac: {j: str}, time_deriv}
        "name_to_state_idx": name_to_state_idx,  # state_name -> int
    }
# =====================================================================
# Forcing initialization code generation
# =====================================================================

def generate_forcing_init_code(n_forcings, num_type="AD"):
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
# Event code generation with analytical saltation gradients
# =====================================================================

# For a root event g(x, t), dg/dx_i and dg/dt are differentiated at build time
# and emitted as lambdas, so the runtime gets the IFT saltation correction
# analytically. Terminal and steady-state events pass nullptr and skip it.

def _get_list_value(dict_or_df, key, index, n_events):
    """Extract a value from a dict-of-lists."""
    if key not in dict_or_df:
        return None
    value = dict_or_df[key]
    if isinstance(value, (list, tuple)):
        if index < len(value):
            return value[index]
        return None
    return value
def _is_valid_value(value):
    """Check whether a value is meaningful (not NA/None/boolean)."""
    if value is None:
        return False
    if isinstance(value, bool):
        return False
    if isinstance(value, numbers.Number):
        try:
            if value != value:  # NaN check
                return False
        except Exception:
            pass
        return True
    str_val = str(value).lower().strip()
    if str_val in {"none", "nan", "na", ""}:
        return False
    if str_val in {"true", "false"}:
        return False
    return True
def _parse_value_or_expression(value, local_symbols, states_list, params_list, 
                                n_states, num_type, forcings_list=None):
    """Interpret a value as numeric literal or symbolic expression."""
    if forcings_list is None:
        forcings_list = []
    
    if not _is_valid_value(value):
        return None

    value_str = str(value).strip()

    # Try numeric literal first
    try:
        float(value_str)
        return value_str
    except (ValueError, TypeError):
        pass

    # Parse as symbolic expression
    try:
        value_expr = _safe_sympify(value_str, local_symbols)
        value_code = _to_cpp(value_expr, states_list, params_list, n_states, 
                            num_type, forcings=forcings_list)
        return str(value_code)
    except Exception as e:
        raise ValueError(f"Failed to parse expression '{value_str}': {e}")
def _generate_root_gradient_lambdas(root_expr, states_list, params_list,
                                     n_states, num_type, forcings_list,
                                     local_symbols, event_idx,
                                     rhs_exprs=None):
    """
    Symbolically differentiate root function g(x, t) and emit C++ lambdas
    for dg/dx (vector), dg/dt (scalar), and G_tt (scalar double).

    G_tt = d(g_dot)/dt is the total second time derivative of g along the
    ODE trajectory, needed for the second-order IFT correction of dt*.
    It is computed by substituting dx_i/dt = f_i into the total derivative:

      G_tt = f^T · H_g · f + grad_g^T · (J_f · f + df/dt)
           + 2·(d²g/dxdt)^T · f + d²g/dt²

    When rhs_exprs is None (ODE RHS not available), G_tt is set to nullptr
    and the runtime falls back to scalar finite difference.

    Parameters
    ----------
    root_expr : sp.Expr
        The parsed SymPy expression for g(x, t).
    states_list : list of str
        State variable names.
    params_list : list of str
        Parameter names.
    n_states : int
        Number of state variables.
    num_type : str
        Numeric type ("AD", "AD2", "double").
    forcings_list : list of str
        Forcing function names.
    local_symbols : dict
        Symbol table for SymPy parsing.
    event_idx : int
        Event index (for comments).
    rhs_exprs : list of sp.Expr, optional
        Parsed SymPy expressions for the ODE RHS (f_0, f_1, ...).
        If provided, G_tt is computed analytically.

    Returns
    -------
    dg_dx_lines : list of str
        C++ lambda lines for dg_dx.
    dg_dt_lines : list of str
        C++ lambda lines for dg_dt.
    g_dot_dot_lines : list of str
        C++ lambda lines for g_dot_dot (G_tt), or ["    nullptr  // g_dot_dot"].
    """
    state_type = f"std::vector<{num_type}>"
    t = local_symbols["time"]

    # --- dg/dx_i for each state ---
    states_syms = [local_symbols[s] for s in states_list]
    dg_dx_exprs = []
    free = root_expr.free_symbols
    for s_sym in states_syms:
        if s_sym in free:
            dg_dx_exprs.append(sp.diff(root_expr, s_sym))
        else:
            dg_dx_exprs.append(sp.Integer(0))

    # --- dg/dt ---
    dg_dt_expr = sp.diff(root_expr, t)

    # sp.diff cannot differentiate through the forcing interpolation, so for a
    # forcing appearing in g the chain rule term (dg/dF_k) * F[k]->derivative(t)
    # is added to the explicit dg/dt.
    forcing_chain_terms = []
    for k, fname in enumerate(forcings_list):
        f_sym = local_symbols[fname]
        if f_sym in free:
            dg_df = sp.diff(root_expr, f_sym)
            if dg_df != 0:
                dg_df_cpp = _to_cpp(dg_df, states_list, params_list, n_states,
                                    num_type, forcings_list)
                dg_df_cpp = str(dg_df_cpp).replace("params[", "full_params[")
                forcing_chain_terms.append(
                    f"({dg_df_cpp})*(*F[{k}]).derivative(t)"
                )

    # --- Emit dg_dx lambda ---
    # All n_states partials go into the output vector, zero entries explicitly:
    # out may be uninitialised.
    dg_dx_lines = []
    dg_dx_lines.append(f"    // dg/dx for root event {event_idx} (analytical, codegen)")
    dg_dx_lines.append(f"    [full_params, &F](const {state_type}& x, const {num_type}& t, {state_type}& out) {{")

    for j, dexpr in enumerate(dg_dx_exprs):
        cpp = _to_cpp(dexpr, states_list, params_list, n_states,
                      num_type, forcings_list)
        cpp = str(cpp).replace("params[", "full_params[")
        dg_dx_lines.append(f"      out[{j}] = {cpp};")

    dg_dx_lines.append(f"    }},  // dg_dx")

    # --- Emit dg_dt lambda ---
    dg_dt_lines = []
    dg_dt_lines.append(f"    // dg/dt for root event {event_idx} (analytical, codegen)")
    dg_dt_lines.append(f"    [full_params, &F](const {state_type}& x, const {num_type}& t) -> {num_type} {{")

    dg_dt_cpp = _to_cpp(dg_dt_expr, states_list, params_list, n_states,
                        num_type, forcings_list)
    dg_dt_cpp = str(dg_dt_cpp).replace("params[", "full_params[")

    if forcing_chain_terms:
        # Combine explicit dg/dt with forcing chain rule terms
        all_terms = [dg_dt_cpp] + forcing_chain_terms if dg_dt_cpp != "0.0" else forcing_chain_terms
        dg_dt_lines.append(f"      return {' + '.join(all_terms)};")
    else:
        dg_dt_lines.append(f"      return {dg_dt_cpp};")

    dg_dt_lines.append(f"    }},  // dg_dt")

    # --- Compute and emit G_tt = d(g_dot)/dt along trajectory ---

    # G_tt is the total time derivative of g_dot = sum(dg/dx_i * f_i) + dg/dt along
    # the trajectory, taken symbolically with dx_i/dt substituted by f_i.
    g_dot_dot_lines = []
    if rhs_exprs is not None and not forcings_list:
        # Build g_dot symbolically: sum(dg/dx_i * f_i) + dg/dt
        g_dot_sym = dg_dt_expr
        for i, s_sym in enumerate(states_syms):
            g_dot_sym += dg_dx_exprs[i] * rhs_exprs[i]

        # Total time derivative of g_dot: sum_j(dg_dot/dx_j * f_j) + dg_dot/dt
        G_tt_sym = sp.diff(g_dot_sym, t)
        g_dot_free = g_dot_sym.free_symbols
        for j, s_sym in enumerate(states_syms):
            if s_sym in g_dot_free:
                G_tt_sym += sp.diff(g_dot_sym, s_sym) * rhs_exprs[j]

        G_tt_sym = sp.powsimp(G_tt_sym)

        # Emit as a double-returning lambda: x and full_params are AD types but G_tt
        # needs scalar values only, so the scalars go into local doubles first.
        g_dot_dot_lines.append(f"    // G_tt = d(g_dot)/dt for root event {event_idx} (analytical, codegen)")
        g_dot_dot_lines.append(f"    [full_params, &F](const {state_type}& x, const {num_type}& t) -> double {{")

        # Locals are named by position: no model symbol reaches the source, and
        # none can shadow the lambda parameters (x, t, full_params).

        # .val() extracts the scalar from an AD type, const-correct: none for double,
        # one level for dual, two for dual2nd.
        if num_type == "AD2":
            xtr = lambda expr: f"({expr}).val().val()"
        elif num_type == "AD":
            xtr = lambda expr: f"({expr}).val()"
        else:
            xtr = lambda expr: expr

        gtt_map = {}
        for j, sname in enumerate(states_list):
            g_dot_dot_lines.append(f"      double _gx{j} = {xtr(f'x[{j}]')};")
            gtt_map[sname] = f"_gx{j}"
        for j, pname in enumerate(params_list):
            g_dot_dot_lines.append(f"      double _gp{j} = {xtr(f'full_params[{n_states + j}]')};")
            gtt_map[pname] = f"_gp{j}"
        g_dot_dot_lines.append(f"      double _gt = {xtr('t')};")
        gtt_map["time"] = "_gt"

        G_tt_cpp = _get_printer(tuple(sorted(gtt_map.items()))).doprint(G_tt_sym)
        # Replace non-standard math macros (single-pass precompiled regex)
        G_tt_cpp = _MATH_MACRO_PATTERN.sub(lambda m: _MATH_MACRO_MAP[m.group(0)], G_tt_cpp)

        g_dot_dot_lines.append(f"      return {G_tt_cpp};")
        g_dot_dot_lines.append(f"    }}  // g_dot_dot")
    else:
        # Forcings present or RHS not available: fall back to FD at runtime
        g_dot_dot_lines.append(f"    nullptr  // g_dot_dot (FD fallback)")

    return dg_dx_lines, dg_dt_lines, g_dot_dot_lines
def fixed_event_time_exprs(events_df, states_list, params_list, n_states,
                           forcings_list=None):
    """C expressions for the fixed-event times, in plain double.

    An event whose time is not among the requested output times adds its own
    row to the result, so the batch entry has to know the event times before
    the solve if it wants to size the output exactly. The expressions address
    the same flat [states, params] vector the solver is handed, so evaluating
    them up front costs nothing.

    Returns one expression per fixed-time event, or None when the model has
    none, has a root event, or an expression reaches for a forcing (which is
    not available before the solve).
    """
    import sympy as sp

    if isinstance(states_list, str):
        states_list = [states_list]
    else:
        states_list = list(states_list)
    if isinstance(params_list, str):
        params_list = [params_list]
    else:
        params_list = list(params_list)
    forcings_list = [] if forcings_list is None else (
        [forcings_list] if isinstance(forcings_list, str) else list(forcings_list))

    if events_df is None or len(events_df) == 0:
        return None
    events_dict = events_df.to_dict("list") if hasattr(events_df, "to_dict") else events_df

    local_symbols = {name: sp.Symbol(name, real=True)
                     for name in list(states_list) + list(params_list) + list(forcings_list)}
    local_symbols["time"] = sp.Symbol("time", real=True)

    list_lengths = [len(v) for v in events_dict.values() if isinstance(v, (list, tuple))]
    n_events = max(list_lengths) if list_lengths else 1

    out = []
    for i in range(n_events):
        if _get_list_value(events_dict, "var", i, n_events) is None:
            continue
        time_raw = _get_list_value(events_dict, "time", i, n_events)
        root_raw = _get_list_value(events_dict, "root", i, n_events)
        if _is_valid_value(root_raw):
            return None                      # root event: the grid is dynamic
        if not _is_valid_value(time_raw):
            continue
        code = _parse_value_or_expression(
            time_raw, local_symbols, states_list, params_list, n_states,
            "double", forcings_list)
        if code is None:
            return None
        for f in forcings_list:
            if f in str(code):
                return None                  # needs a forcing: not known up front
        out.append(str(code))
    return out or None


def generate_event_code(events_df, states_list, params_list, n_states,
                        num_type="AD", forcings_list=None, rhs_dict=None):
    """
    Generate C++ initialization lines for fixed-time and root events.

    For each root event with root function g(x, t), this function:
      1. Parses g to a SymPy expression
      2. Computes dg/dx_i (i = 0..n_states-1) and dg/dt symbolically
      3. Emits C++ lambdas for the RootEvent's dg_dx and dg_dt members
         (used by the analytical IFT-based saltation correction at runtime)
      4. If rhs_dict is provided, computes G_tt = d(g_dot)/dt analytically
         for the second-order IFT correction (otherwise FD fallback at runtime)

    Terminal root events get nullptr for dg_dx / dg_dt / g_dot_dot since they
    don't modify state and don't need saltation correction.

    Fixed-time events don't need dg_dx / dg_dt (the timing residual
    dt_corr = t_event - scalar(t_event) is known directly from the
    AD type of t_event).

    Parameters
    ----------
    rhs_dict : dict, optional
        Dictionary mapping state names to RHS expression strings.
        When provided, enables analytical G_tt computation for the
        second-order IFT correction.
    """
    if forcings_list is None:
        forcings_list = []
    if states_list is None:
        states_list = []
    if params_list is None:
        params_list = []
    
    # Normalize inputs
    if isinstance(states_list, str):
        states_list = [states_list]
    else:
        states_list = list(states_list)
    if isinstance(params_list, str):
        params_list = [params_list]
    else:
        params_list = list(params_list)
    if isinstance(forcings_list, str):
        forcings_list = [forcings_list]
    else:
        forcings_list = list(forcings_list)
    
    if events_df is None or len(events_df) == 0:
        return []

    if hasattr(events_df, "to_dict"):
        events_dict = events_df.to_dict("list")
    else:
        events_dict = events_df

    # Create symbols
    states_syms = {name: sp.Symbol(name, real=True) for name in states_list}
    params_syms = {name: sp.Symbol(name, real=True) for name in params_list}
    t = sp.Symbol("time", real=True)

    local_symbols = {}
    local_symbols.update(states_syms)
    local_symbols.update(params_syms)
    for name in forcings_list:
        local_symbols[name] = sp.Symbol(name, real=True)
    local_symbols["time"] = t

    event_lines = []
    state_type = f"std::vector<{num_type}>"

    # Parse ODE RHS expressions for analytical G_tt computation
    rhs_exprs_parsed = None
    if rhs_dict is not None:
        try:
            rhs_exprs_parsed = [
                _safe_sympify(str(rhs_dict[s]), local_symbols) for s in states_list
            ]
        except Exception:
            rhs_exprs_parsed = None  # fall back to FD if parsing fails

    # Determine number of events
    list_lengths = []
    for v in events_dict.values():
        if isinstance(v, (list, tuple)):
            list_lengths.append(len(v))
    n_events = max(list_lengths) if list_lengths else 1

    for i in range(n_events):
        var_raw = _get_list_value(events_dict, "var", i, n_events)
        if var_raw is None:
            continue

        var_name = str(var_raw)
        if var_name not in states_list:
            raise ValueError(f"Event {i}: unknown state variable '{var_name}'")
        var_idx = states_list.index(var_name)

        value_raw = _get_list_value(events_dict, "value", i, n_events)
        value_code = _parse_value_or_expression(
            value_raw, local_symbols, states_list, params_list, n_states, 
            num_type, forcings_list
        )
        if value_code is None:
            raise ValueError(f"Event {i}: 'value' is required but is NA/None")
        value_code = str(value_code).replace("params[", "full_params[")

        method_raw = _get_list_value(events_dict, "method", i, n_events)
        method = str(method_raw).lower() if method_raw is not None else "replace"
        method_map = {
            "replace": "EventMethod::Replace",
            "add": "EventMethod::Add",
            "multiply": "EventMethod::Multiply",
        }
        method_code = method_map.get(method, "EventMethod::Replace")

        time_raw = _get_list_value(events_dict, "time", i, n_events)
        time_code = _parse_value_or_expression(
            time_raw, local_symbols, states_list, params_list, n_states, 
            num_type, forcings_list
        )

        root_raw = _get_list_value(events_dict, "root", i, n_events)
        root_code = _parse_value_or_expression(
            root_raw, local_symbols, states_list, params_list, n_states, 
            num_type, forcings_list
        )

        if time_code is not None:
            # ============================================================
            # Fixed-time event
            # ============================================================

            # No dg_dx / dg_dt: the saltation correction uses dt_corr = t_event -
            # scalar(t_event), built from the AD type of the event time.
            time_code = str(time_code).replace("params[", "full_params[")
            
            event_lines.append(f"  // Fixed event {i}: {var_name} at t = {time_raw}")
            event_lines.append(f"  fixed_events.emplace_back(FixedEvent<{state_type}, {num_type}>{{")
            event_lines.append(f"    {time_code},  // time")
            event_lines.append(f"    {var_idx},    // state_index")
            event_lines.append(f"    [full_params, &F](const {state_type}& x, const {num_type}& t) -> {num_type} {{")
            event_lines.append(f"      return {value_code};")
            event_lines.append(f"    }},  // value_func")
            event_lines.append(f"    {method_code}  // method")
            event_lines.append(f"  }});")
            event_lines.append("")
            
        elif root_code is not None:
            # ============================================================
            # Root-finding event
            # ============================================================

            # g(x, t) is parsed symbolically for the dg/dx and dg/dt that the IFT
            # saltation correction needs at runtime.
            root_code = str(root_code).replace("params[", "full_params[")
            
            # Get optional parameters
            terminal_raw = _get_list_value(events_dict, "terminal", i, n_events)
            terminal = "true" if terminal_raw and str(terminal_raw).lower() == "true" else "false"
            
            direction_raw = _get_list_value(events_dict, "direction", i, n_events)
            try:
                direction = int(direction_raw) if direction_raw is not None else 0
            except (ValueError, TypeError):
                direction = 0

            # --- Analytical gradients of g for saltation correction ---

            # A non-terminal event modifies state and needs the full correction, so dg/dx
            # and dg/dt come from the SymPy expression of g. A terminal event only stops
            # the integration and gets nullptr.
            is_terminal = (terminal == "true")

            event_lines.append(f"  // Root event {i}: {var_name} when {root_raw} = 0")
            event_lines.append(f"  root_events.push_back(RootEvent<{state_type}, {num_type}>{{")
            event_lines.append(f"    [full_params, &F](const {state_type}& x, const {num_type}& t) -> {num_type} {{")
            event_lines.append(f"      return {root_code};")
            event_lines.append(f"    }},  // func (root condition g)")
            event_lines.append(f"    {var_idx},  // state_index")
            event_lines.append(f"    [full_params, &F](const {state_type}& x, const {num_type}& t) -> {num_type} {{")
            event_lines.append(f"      return {value_code};")
            event_lines.append(f"    }},  // value_func (h)")
            event_lines.append(f"    {method_code},  // method")
            event_lines.append(f"    {terminal},     // terminal")
            event_lines.append(f"    {direction},    // direction")

            if is_terminal:
                # Terminal event: no state modification, no saltation needed
                event_lines.append(f"    nullptr,       // dg_dx (not needed for terminal)")
                event_lines.append(f"    nullptr,       // dg_dt (not needed for terminal)")
                event_lines.append(f"    nullptr        // g_dot_dot (not needed for terminal)")
            else:
                # Non-terminal: generate analytical dg/dx, dg/dt, g_dot_dot lambdas
                root_raw_str = str(root_raw).strip()
                root_sympy = _safe_sympify(root_raw_str, local_symbols)

                dg_dx_lines, dg_dt_lines, g_dot_dot_lines = _generate_root_gradient_lambdas(
                    root_sympy, states_list, params_list, n_states,
                    num_type, forcings_list, local_symbols, i,
                    rhs_exprs=rhs_exprs_parsed
                )
                event_lines.extend(dg_dx_lines)
                event_lines.extend(dg_dt_lines)
                event_lines.extend(g_dot_dot_lines)

            event_lines.append(f"  }});")
            event_lines.append("")
            
        else:
            raise ValueError(f"Event {i}: must specify either 'time' or 'root'")

    return event_lines
# =====================================================================
# Root function code generation
# =====================================================================

def generate_rootfunc_code(rootfunc, states_list, params_list, n_states,
                           num_type="AD", forcings_list=None):
    """
    Generate C++ code for root function based termination.
    
    Handles two cases:
    1. rootfunc = "equilibrate": steady-state detection
    2. rootfunc = list of expressions: stop when any crosses zero
    """
    if forcings_list is None:
        forcings_list = []
    
    if rootfunc is None:
        return []
    
    # Handle single string vs list
    if isinstance(rootfunc, str):
        if rootfunc.strip().lower() == "equilibrate":
            return _generate_equilibrate_code(num_type)
        else:
            rootfunc = [rootfunc]
    
    if not isinstance(rootfunc, (list, tuple)):
        raise ValueError(f"rootfunc must be 'equilibrate' or a list of expressions, got: {type(rootfunc)}")
    
    return _generate_user_rootfunc_code(
        rootfunc, states_list, params_list, n_states, num_type, forcings_list
    )
def _generate_equilibrate_code(num_type):
    """
    Generate C++ code for steady-state detection (direct threshold check).

    Checks after each step whether max(|dxdt|) across all AD levels
    is below root_tol.  Uses the termination callback mechanism :
    no root-finding or sign-change detection needed.
    """
    state_type = f"std::vector<{num_type}>"

    return [
        "",
        "  // --- Steady-state termination (rootfunc = 'equilibrate') ---",
        f"  auto ss_termination = make_steady_state_termination<ode_system, {state_type}, {num_type}>(sys, root_tol);",
        ""
    ]
def _generate_user_rootfunc_code(rootfunc_list, states_list, params_list, 
                                  n_states, num_type, forcings_list):
    """
    Generate C++ code for user-defined root expressions (terminal events).

    Terminal root events stop integration when g(x, t) crosses zero.
    They don't modify state, so no saltation correction is needed :
    dg_dx / dg_dt are set to nullptr.
    """
    states_syms = {name: sp.Symbol(name, real=True) for name in states_list}
    params_syms = {name: sp.Symbol(name, real=True) for name in params_list}
    t = sp.Symbol("time", real=True)
    
    local_symbols = {}
    local_symbols.update(states_syms)
    local_symbols.update(params_syms)
    for name in forcings_list:
        local_symbols[name] = sp.Symbol(name, real=True)
    local_symbols["time"] = t
    
    state_type = f"std::vector<{num_type}>"
    lines = [
        "",
        "  // --- User-defined root function termination ---"
    ]
    
    for i, expr_str in enumerate(rootfunc_list):
        expr_str = str(expr_str).strip()
        if not expr_str:
            continue
        
        try:
            expr = _safe_sympify(expr_str, local_symbols)
            root_code = _to_cpp(expr, states_list, params_list, n_states, 
                               num_type, forcings_list)
            root_code = str(root_code).replace("params[", "full_params[")
        except Exception as e:
            raise ValueError(f"Failed to parse rootfunc expression '{expr_str}': {e}")
        
        lines.extend([
            f"  // rootfunc[{i}]: {expr_str}",
            f"  root_events.push_back(RootEvent<{state_type}, {num_type}>{{",
            f"    [full_params, &F](const {state_type}& x, const {num_type}& t) -> {num_type} {{",
            f"      return {root_code};",
            f"    }},  // func",
            f"    0,            // state_index (ignored for terminal)",
            f"    [](const {state_type}&, const {num_type}&) {{ return {num_type}(0.0); }},  // value_func",
            f"    EventMethod::Replace,  // method",
            f"    true,         // terminal = true",
            f"    0,            // direction = 0 (any crossing)",
            f"    nullptr,      // dg_dx (not needed for terminal)",
            f"    nullptr,      // dg_dt (not needed for terminal)",
            f"    nullptr       // g_dot_dot (not needed for terminal)",
            f"  }});",
        ])
    
    lines.append("")
    return lines
# =====================================================================
# Sparse LU pattern code generation
# =====================================================================

def _extract_jac_sparsity(jac_matrix, n):
    """Extract non-zero positions from a symbolic Jacobian matrix."""
    pattern = set()
    _zero = sp.Integer(0)
    _szero = sp.S.Zero
    for i in range(n):
        row = jac_matrix[i]
        for j in range(n):
            expr = row[j]
            if expr is not _zero and expr is not _szero and expr != 0:
                pattern.add((i, j))
    return pattern


def decide_sparse(sparse, n_states, jac_nnz, has_jacobian=True):
    """
    Dense vs sparse linear solver, shared by the cppDE and CVODE backends.

    `sparse` arrives already normalised by R (None = auto, True/False = pinned);
    the KLU-availability downgrade happens there, above this call.

    Auto-detection selects sparse from 8 states unless the structural Jacobian
    is denser than 0.4.  Both bounds are calibrated on the benchmark suite
    (benchmarks/run-benchmarks.R --sparse-sweep).

    `has_jacobian` is False for explicit methods, whose nnz count is 0; the
    auto path stays dense there.
    """
    if sparse is not None:
        return bool(sparse)
    if not has_jacobian or n_states <= 0:
        return False
    density = jac_nnz / float(n_states * n_states)
    return (n_states >= 8) and (density <= 0.4)


def analyze_klu_settings(n, jac_nnz_rows, jac_nnz_cols):
    """
    Analyze the Jacobian sparsity pattern to determine optimal KLU settings.

    Returns a dict with:
      - use_btf (bool):  True if BTF decomposition is beneficial
      - ordering (int):  0=AMD, 1=COLAMD

    BTF decision:
      Find strongly connected components of the directed graph defined
      by the Jacobian pattern.  If nblocks > 1, BTF can decompose the
      problem into smaller independent blocks, so BTF pays.
      If nblocks == 1 (strongly connected), BTF is pure overhead and stays off.

    Ordering decision:
      For PDE-like patterns (uniform row degree, wide bandwidth),
      AMD typically produces less fill-in than COLAMD.
      For irregular patterns (high degree variance, hub nodes),
      COLAMD often wins.
      Heuristic: if the coefficient of variation of row degrees > 0.5, use COLAMD.
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
