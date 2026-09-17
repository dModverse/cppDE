import _tap

_real = _tap.load_real("codegen_cvode")
globals().update(_tap.wrap("codegen_cvode", _real, ("generate_cvode_cpp",)))


def __getattr__(name):
    return getattr(_real, name)
