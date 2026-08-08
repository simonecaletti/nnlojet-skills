#!/usr/bin/env python3
"""One IR-limit vocabulary for all probe generators.

The mode families are the ones genuine_modes.py (write-spike-test)
defines: ss (single soft), sco (single collinear), ds (double soft),
tc (triple collinear), sc (soft + collinear), dc (double collinear).
A limit is named as

    {"family": "sc", "unresolved": [i, j, k], "spectators": [l, m]}

with the sc convention "first index soft, remaining pair collinear"
and dc "two pairs in order". This module turns such a limit into the
Fortran statements that drive phase space there (librambo get_*
generators; mass conventions of the check-program suite), so specs can
name limits instead of hand-writing `limit_call` — which remains the
escape hatch for anything not covered.

Self-contained; imported by gen_probe.py and antenna_probe.py.
"""

FAMILIES = ("ss", "sco", "ds", "tc", "sc", "dc")

# em conventions: soft-type drivers take the recoil-system mass
# em = sqrts*sqrt(1-x); collinear-type the cluster mass em = sqrts*x.
_SOFT = "em1=sqrts_proc*dsqrt(1d0-xs)"
_COLL = "em1=sqrts_proc*xs"


def _call(name, *args):
    alist = ",".join(str(a) for a in args if a != "")
    return f"call {name}(sqrts_proc,{alist})"


def limit_lines(family, npar, unresolved, spectators):
    """-> list of Fortran statements (unindented) driving the limit."""
    u = list(unresolved)
    S = ",".join(str(x) for x in spectators)
    if family == "ss":
        return [_SOFT, _call(f"get_ss{npar}", u[0], "em1", S)]
    if family == "sco":
        return [_COLL, _call(f"get_sco{npar}", u[0], u[1], "em1", S)]
    if family == "ds":
        return [_SOFT, _call(f"get_ds{npar}", u[0], u[1], "em1", S)]
    if family == "tc":
        return [_COLL, _call(f"get_tc{npar}", u[0], u[1], u[2], "em1", S)]
    if family == "sc":  # u[0] soft, u[1]||u[2]
        return [_SOFT, "em2=sqrts_proc*xs",
                _call(f"get_sc{npar}", u[0], u[1], u[2], "em1", "em2", S)]
    if family == "dc":  # u[0]||u[1], u[2]||u[3]
        return [_COLL, "em2=sqrts_proc*xs",
                _call(f"get_dc{npar}", u[0], u[1], "em1",
                      u[2], u[3], "em2", S)]
    raise ValueError(f"unknown family {family!r} (want one of {FAMILIES})")


def s_expr(family, npar, unresolved):
    """Dimensionless vanishing-scale expression for the pole scan:
    the invariant for collinear-type families, the scan parameter x
    for soft-type ones (no single invariant vanishes there)."""
    u = list(unresolved)
    if family == "sco":
        return f"abs(s{npar}({u[0]},{u[1]}))/sqrts_proc**2"
    if family == "tc":
        return (f"(abs(s{npar}({u[0]},{u[1]}))+abs(s{npar}({u[0]},{u[2]}))"
                f"+abs(s{npar}({u[1]},{u[2]})))/sqrts_proc**2")
    if family == "sc":
        return f"abs(s{npar}({u[1]},{u[2]}))/sqrts_proc**2"
    if family == "dc":
        return f"abs(s{npar}({u[0]},{u[1]}))/sqrts_proc**2"
    return "xs"  # ss, ds


def mode_name(family, unresolved):
    """Same naming as genuine_modes.classify."""
    u = list(unresolved)
    if family == "ss":
        return f"{u[0]} soft"
    if family == "sco":
        return f"{u[0]}||{u[1]}"
    if family == "ds":
        return f"{u[0]},{u[1]} soft"
    if family == "tc":
        return f"{u[0]}||{u[1]}||{u[2]}"
    if family == "sc":
        return f"{u[0]} soft + {u[1]}||{u[2]}"
    if family == "dc":
        return f"{u[0]}||{u[1]} + {u[2]}||{u[3]}"
    raise ValueError(family)


def from_spec(limit, npar):
    """Spec dict -> (lines, s_expr, name). Validates arity."""
    fam = limit["family"]
    u = limit["unresolved"]
    need = {"ss": 1, "sco": 2, "ds": 2, "tc": 3, "sc": 3, "dc": 4}[fam]
    if len(u) != need:
        raise ValueError(f"{fam} needs {need} unresolved indices, "
                         f"got {len(u)}")
    sp = limit.get("spectators", [])
    return (limit_lines(fam, npar, u, sp),
            s_expr(fam, npar, u),
            mode_name(fam, u))


def selftest():
    """Structural checks only — no physics."""
    ls, se, nm = from_spec({"family": "sc", "unresolved": [5, 3, 4],
                            "spectators": [6, 7]}, 7)
    assert ls[-1] == "call get_sc7(sqrts_proc,5,3,4,em1,em2,6,7)", ls
    assert "em2=" in ls[1] and "dsqrt" in ls[0]
    assert se == "abs(s7(3,4))/sqrts_proc**2"
    assert nm == "5 soft + 3||4"
    ls, se, nm = from_spec({"family": "dc", "unresolved": [3, 4, 5, 6],
                            "spectators": []}, 6)
    assert ls[-1] == "call get_dc6(sqrts_proc,3,4,em1,5,6,em2)", ls
    ls, se, nm = from_spec({"family": "ss", "unresolved": [4],
                            "spectators": [3, 5, 6, 7]}, 7)
    assert ls == ["em1=sqrts_proc*dsqrt(1d0-xs)",
                  "call get_ss7(sqrts_proc,4,em1,3,5,6,7)"]
    assert se == "xs" and nm == "4 soft"
    try:
        from_spec({"family": "tc", "unresolved": [3, 4]}, 7)
        raise AssertionError("arity not checked")
    except ValueError:
        pass
    print("irlimits selftest OK")


if __name__ == "__main__":
    selftest()
