"""Generator wrappers that log each call as a JSON line.

Environment: CPPDE_TAP_REAL (directory of the real generators), CPPDE_TAP_LOG
(log directory, one file per process).
"""

import importlib.util
import json
import os
import sys

REAL = os.environ["CPPDE_TAP_REAL"]
LOG = os.environ["CPPDE_TAP_LOG"]

if REAL not in sys.path:
    sys.path.insert(0, REAL)


def load_real(name):
    spec = importlib.util.spec_from_file_location(
        "_real_" + name, os.path.join(REAL, name + ".py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _plain(x):
    if hasattr(x, "to_dict"):
        try:
            return {k: [_plain(v) for v in vs]
                    for k, vs in x.to_dict("list").items()}
        except TypeError:
            pass
    if isinstance(x, dict):
        return {str(k): _plain(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_plain(v) for v in x]
    if isinstance(x, (str, int, bool)) or x is None:
        return x
    if isinstance(x, float):
        return x if x == x else None
    try:
        return float(x)
    except (TypeError, ValueError):
        return str(x)


def record(module, func, kwargs, args=()):
    try:
        os.makedirs(LOG, exist_ok=True)
        path = os.path.join(LOG, "calls-%d.jsonl" % os.getpid())
        with open(path, "a", encoding="utf-8") as f:
            f.write(json.dumps({"module": module, "func": func,
                                "args": _plain(list(args)),
                                "kwargs": _plain(kwargs)}) + "\n")
    except Exception:
        pass


def wrap(module, real, names):
    def make(name):
        fn = getattr(real, name)

        def wrapped(*args, **kwargs):
            record(module, name, kwargs, args)
            return fn(*args, **kwargs)
        wrapped.__name__ = name
        return wrapped
    return {name: make(name) for name in names}
