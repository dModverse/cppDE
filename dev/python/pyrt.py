"""Runtime of the Python backend of cppde_emit (imported as R).

Functions accept float or (nested) Dual and return IEEE results instead of
raising.
"""

import math

INF = math.inf
NAN = math.nan


def _val(x):
    while isinstance(x, Dual):
        x = x.v
    return x


class Dual:
    """Forward-mode dual v + sum(d[k] * eps_k); v, d[k] float or Dual."""

    __slots__ = ("v", "d")

    def __init__(self, v, d):
        self.v = v
        self.d = tuple(d)

    @staticmethod
    def lift(x, n):
        if isinstance(x, Dual):
            return x
        z = 0.0 * x if not isinstance(x, (int, float)) else 0.0
        return Dual(x, (z,) * n)

    def __repr__(self):
        return "Dual(%r, %r)" % (self.v, self.d)

    def _n(self):
        return len(self.d)

    # arithmetic
    def __add__(self, o):
        if isinstance(o, Dual):
            return Dual(self.v + o.v, [a + b for a, b in zip(self.d, o.d)])
        return Dual(self.v + o, self.d)

    __radd__ = __add__

    def __sub__(self, o):
        if isinstance(o, Dual):
            return Dual(self.v - o.v, [a - b for a, b in zip(self.d, o.d)])
        return Dual(self.v - o, self.d)

    def __rsub__(self, o):
        return Dual(o - self.v, [-a for a in self.d])

    def __neg__(self):
        return Dual(-self.v, [-a for a in self.d])

    def __pos__(self):
        return self

    def __mul__(self, o):
        if isinstance(o, Dual):
            return Dual(self.v * o.v,
                        [a * o.v + self.v * b for a, b in zip(self.d, o.d)])
        return Dual(self.v * o, [a * o for a in self.d])

    __rmul__ = __mul__

    def __truediv__(self, o):
        return div(self, o)

    def __rtruediv__(self, o):
        return div(o, self)

    def __pow__(self, o):
        return pow(self, o)

    def __rpow__(self, o):
        return pow(o, self)

    # comparisons use the innermost value
    def __lt__(self, o):
        return _val(self) < _val(o)

    def __le__(self, o):
        return _val(self) <= _val(o)

    def __gt__(self, o):
        return _val(self) > _val(o)

    def __ge__(self, o):
        return _val(self) >= _val(o)

    def __eq__(self, o):
        return _val(self) == _val(o)

    def __ne__(self, o):
        return _val(self) != _val(o)

    def __bool__(self):
        return bool(_val(self))

    __hash__ = None


def _chain(x, f, df):
    """f(x) for Dual x, with df the derivative of f."""
    fv = f(x.v)
    dv = df(x.v)
    return Dual(fv, [dv * a for a in x.d])


def _safe(fn, x):
    try:
        r = fn(x)
    except (ValueError, ZeroDivisionError):
        return NAN
    except OverflowError:
        return INF
    return r


def div(a, b):
    if isinstance(a, Dual) or isinstance(b, Dual):
        if not isinstance(b, Dual):
            return Dual(div(a.v, b), [div(x, b) for x in a.d])
        inv = div(1.0, b.v)
        if not isinstance(a, Dual):
            return Dual(a * inv, [-(a * inv * inv) * x for x in b.d])
        q = a.v * inv
        return Dual(q, [(x - q * y) * inv for x, y in zip(a.d, b.d)])
    try:
        return a / b
    except ZeroDivisionError:
        if a != a or a == 0:
            return NAN
        return math.copysign(INF, a) * math.copysign(1.0, b)


def pow(a, b):
    if isinstance(b, Dual):
        # a**b = exp(b*log(a))
        return exp(b * log(a))
    if isinstance(a, Dual):
        fv = pow(a.v, b)
        dv = b * pow(a.v, b - 1.0) if b != 0 else 0.0
        return Dual(fv, [dv * x for x in a.d])
    if a == 0 and b < 0:
        # IEEE pole: -inf only for -0.0 and an odd integer exponent
        odd = float(b).is_integer() and int(b) % 2 == 1
        return math.copysign(INF, a) if odd else INF
    try:
        r = math.pow(a, b)
    except ValueError:
        return NAN
    except OverflowError:
        return INF
    return r


def sqrt(x):
    if isinstance(x, Dual):
        s = sqrt(x.v)
        return Dual(s, [div(a, 2.0 * s) for a in x.d])
    return NAN if x < 0 else math.sqrt(x)


def _unary(fn, dfn):
    def f(x):
        if isinstance(x, Dual):
            return _chain(x, f, dfn)
        return _safe(fn, x)
    return f


def _log(x):
    if x == 0:
        return -INF
    if x != x or x < 0:
        return NAN
    return INF if x == INF else math.log(x)


def _atanh(x):
    if x == 1:
        return INF
    if x == -1:
        return -INF
    return _safe(math.atanh, x)


exp = _unary(math.exp, lambda v: exp(v))
log = _unary(_log, lambda v: div(1.0, v))
sin = _unary(math.sin, lambda v: cos(v))
cos = _unary(math.cos, lambda v: -sin(v))
tan = _unary(math.tan, lambda v: 1.0 + tan(v) * tan(v))
asin = _unary(math.asin, lambda v: div(1.0, sqrt(1.0 - v * v)))
acos = _unary(math.acos, lambda v: -div(1.0, sqrt(1.0 - v * v)))
atan = _unary(math.atan, lambda v: div(1.0, 1.0 + v * v))
sinh = _unary(math.sinh, lambda v: cosh(v))
cosh = _unary(math.cosh, lambda v: sinh(v))
tanh = _unary(math.tanh, lambda v: 1.0 - tanh(v) * tanh(v))
asinh = _unary(math.asinh, lambda v: div(1.0, sqrt(v * v + 1.0)))
acosh = _unary(math.acosh, lambda v: div(1.0, sqrt(v * v - 1.0)))
atanh = _unary(_atanh, lambda v: div(1.0, 1.0 - v * v))
erf = _unary(math.erf,
             lambda v: (2.0 / math.sqrt(math.pi)) * exp(-(v * v)))
erfc = _unary(math.erfc,
              lambda v: -(2.0 / math.sqrt(math.pi)) * exp(-(v * v)))


def abs(x):
    if isinstance(x, Dual):
        return -x if _val(x) < 0 else x
    return math.fabs(x)


def _const(fn):
    def f(x):
        return fn(_val(x))
    return f


floor = _const(lambda v: float(math.floor(v)) if math.isfinite(v) else v)
ceiling = _const(lambda v: float(math.ceil(v)) if math.isfinite(v) else v)


def Heaviside(x):
    v = _val(x)
    return 0.0 if v < 0 else (0.5 if v == 0 else 1.0)


def sign(x):
    v = _val(x)
    return 1.0 if v > 0 else (-1.0 if v < 0 else 0.0)


def delta(x):
    return 1.0 if _val(x) == 0 else 0.0


def select(c, a, b):
    return a if c else b


def min(a, b):
    return b if _val(b) < _val(a) else a


def max(a, b):
    return b if _val(a) < _val(b) else a


def gamma(x):
    if isinstance(x, Dual):
        raise TypeError("gamma has no AD rule in the harness")
    return _safe(math.gamma, x)


def loggamma(x):
    if isinstance(x, Dual):
        raise TypeError("loggamma has no AD rule in the harness")
    return _safe(math.lgamma, x)


def atan2(y, x):
    if isinstance(x, Dual) or isinstance(y, Dual):
        r2 = x * x + y * y
        n = len(x.d) if isinstance(x, Dual) else len(y.d)
        xv, yv = (x.v if isinstance(x, Dual) else x), (y.v if isinstance(y, Dual) else y)
        dx = x.d if isinstance(x, Dual) else (0.0,) * n
        dy = y.d if isinstance(y, Dual) else (0.0,) * n
        rv = r2.v if isinstance(r2, Dual) else r2
        return Dual(atan2(yv, xv), [div(xv * b - yv * a, rv) for a, b in zip(dx, dy)])
    return math.atan2(y, x)


# ---------------------------------------------------------------------------
# Seeding and reading
# ---------------------------------------------------------------------------

def seed(values, dirs):
    """Duals with values[i] and tangents dirs[i]."""
    return [Dual(v, d) for v, d in zip(values, dirs)]


def identity_seed(values):
    n = len(values)
    return [Dual(v, [1.0 if k == i else 0.0 for k in range(n)])
            for i, v in enumerate(values)]


def nested_identity(values):
    """Nested duals, both layers seeded with the identity."""
    n = len(values)
    out = []
    for i, v in enumerate(values):
        inner = Dual(v, [1.0 if k == i else 0.0 for k in range(n)])
        tan = [Dual(1.0 if k == i else 0.0, (0.0,) * n) for k in range(n)]
        out.append(Dual(inner, tan))
    return out


def value(x):
    return _val(x)


def tangent(x, k):
    if not isinstance(x, Dual):
        return 0.0
    return _val(x.d[k]) if k < len(x.d) else 0.0


def hessian(x, k, m):
    """Second derivative (k, m) of a nested dual."""
    if not isinstance(x, Dual):
        return 0.0
    tk = x.d[k]
    if not isinstance(tk, Dual):
        return 0.0
    return _val(tk.d[m])


# ---------------------------------------------------------------------------
# Stand-ins for the C++ runtime
# ---------------------------------------------------------------------------

def row_dot(m, r, x):
    """sum_j C[r, j] * x[j] of a cppde_struct.LinMap."""
    s = 0.0
    for j, c in m.rows[r]:
        s = s + float(c) * x[j]
    return s


def apply(m, x, y):
    for r in range(len(m.rows)):
        y[r] = row_dot(m, r, x)


def axpy_row_dense(m, r, a, y, i=None):
    """y[j] += a * C[r, j], or y[i, j] with a row i."""
    for j, c in m.rows[r]:
        if i is None:
            y[j] = y[j] + a * float(c)
        else:
            y[i, j] = y[i, j] + a * float(c)


def axpy_row_idx(m, r, a, ax, dst, off):
    for k, (_, c) in enumerate(m.rows[r]):
        ax[dst[off + k]] = ax[dst[off + k]] + a * float(c)


def apply_t_add(m, w, y):
    for r in range(len(m.rows)):
        axpy_row_dense(m, r, w[r], y)


class Forcing:
    """PchipForcing stand-in: f(t) and derivative(t)."""

    def __init__(self, f, df):
        self.f = f
        self.df = df

    def __call__(self, t):
        return self.f(t)

    def derivative(self, t):
        return self.df(t)


class DenseMatrix:
    """dense_matrix<T> stand-in with [i, j] access."""

    def __init__(self, n, zero=0.0):
        self.n = n
        self.data = [[zero] * n for _ in range(n)]

    def __getitem__(self, ij):
        return self.data[ij[0]][ij[1]]

    def __setitem__(self, ij, v):
        self.data[ij[0]][ij[1]] = v

    def set_zero(self):
        self.data = [[0.0] * self.n for _ in range(self.n)]


class SunSparse:
    """SUNSparseMatrix stand-in (CSC arrays written by the callback)."""

    def __init__(self, n, nnz):
        self.n = n
        self.indexptrs = [0] * (n + 1)
        self.indexvals = [0] * nnz
        self.data = [0.0] * nnz

    def dense(self):
        m = [[0.0] * self.n for _ in range(self.n)]
        for j in range(self.n):
            for p in range(self.indexptrs[j], self.indexptrs[j + 1]):
                m[self.indexvals[p]][j] = self.data[p]
        return m


class CscMatrix:
    """csc_matrix<T> stand-in."""

    def __init__(self):
        self.pattern_built = False
        self.n = 0
        self.Ap = []
        self.Ai = []
        self.Ax = []

    def build_pattern(self, n, nnz, rows, cols):
        self.n = n
        counts = [0] * (n + 1)
        for c in cols:
            counts[c + 1] += 1
        for j in range(n):
            counts[j + 1] += counts[j]
        self.Ap = counts
        cursor = list(counts[:n])
        self.Ai = [0] * nnz
        for k in range(nnz):
            pos = cursor[cols[k]]
            cursor[cols[k]] += 1
            self.Ai[pos] = rows[k]
        self.Ax = [0.0] * nnz
        self.pattern_built = True

    def dense(self):
        m = [[0.0] * self.n for _ in range(self.n)]
        for j in range(self.n):
            for p in range(self.Ap[j], self.Ap[j + 1]):
                m[self.Ai[p]][j] = self.Ax[p]
        return m
