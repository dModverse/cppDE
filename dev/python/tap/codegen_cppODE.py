import _tap

_real = _tap.load_real("codegen_cppODE")
globals().update(_tap.wrap("codegen_cppODE", _real, (
    "generate_ode_cpp", "generate_event_code", "generate_rootfunc_code",
    "fixed_event_time_exprs")))


def __getattr__(name):
    return getattr(_real, name)
