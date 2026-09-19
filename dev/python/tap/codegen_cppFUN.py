import _tap

_real = _tap.load_real("codegen_cppFUN")
globals().update(_tap.wrap("codegen_cppFUN", _real, ("generate_fun_cpp",)))


def __getattr__(name):
    return getattr(_real, name)
