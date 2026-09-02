#!/usr/bin/env python3
"""Per-pole COEFFICIENT balance between a term's X40 lines and the
iterated counterterms that cancel their spurious single poles.

pole_ledger.py checks that a spurious single pole has SOME negative
iterated partner -- existence.  It says so itself: "a block touching the
right invariants with the wrong coefficients or mapping still passes".
This script is the arithmetic the SKILL then asks you to do by hand:

    "for EACH pole, the coefficients of the iterated lines sharing that
     pole must sum against the X40's coefficient on that pole.  Verify
     the per-pole sums directly in your file."

For every invariant that is a spurious single of at least one X40 line it
prints the signed coefficient sums

    x40 = sum of (sign*coeff) over X40 lines singular there
    itr = sum of (sign*coeff) over iterated lines (>=2 antennae) singular
          there

and their total.  Usage:

  pairing_balance.py TERM.map --spec TERM.spec.json [--strict] [--quiet]
  pairing_balance.py --selftest

MEASURED SCOPE -- what this catches and what it does NOT.
  `x40 + itr == 0` is NOT an invariant.  One iterated line carries
  several poles and cancels different X40 halves in different kinematic
  regions, so a correct term shows nonzero totals.  Calibrated on one
  validated 5-parton epem term in four configurations, three of them
  with a known spike-test verdict:

    configuration                        sum|total|   spike test
    validated                                12       21/21 modes
    wrong RADIATOR arguments in a
      counterterm family                     12       15/21 modes
    counterterm family missing               20       ~0.4 everywhere
    two rival families composed together     16       broken

  Read that table before relying on this script:
   * it SEPARATES count/coefficient errors -- a missing family, a
     duplicated family -- which is exactly the class pole_ledger passes,
     because existence holds and only the arithmetic is wrong;
   * it is BLIND to argument errors.  Rows 1 and 2 are identical: same
     lines, same coefficients, different antenna arguments.  Only a run
     finds that (run-spike-test, then fit_lines.py).
   * there is NO absolute pass mark.  The signal is the DIFFERENCE from
     a term known to be correct, which is what --baseline compares.

  Hard errors, independent of any baseline:
   * DANGLING  -- x40 != 0 and itr == 0: nothing can cancel that pole
     (pole_ledger reports this too; kept here so one command suffices).
   * SAME-SIGN -- x40 and itr have the same sign: a counterterm that
     ADDS to the pole it is supposed to remove.
   * ORPHAN    -- itr != 0 and x40 == 0 on a single-unresolved pole.
  --strict additionally errors on any |total| > --tol.  Use it only when
  a sibling family is known to balance exactly -- not as a general gate.
"""
import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from _map_parser import parse_map, norm                      # noqa: E402
from _colour_algebra import make_valid_cluster               # noqa: E402
import pole_ledger as PL                                     # noqa: E402


def build_lines(mapfile, spec, sheet):
    """LineCheck objects, same construction as run_ledger (checks off)."""
    flavours = spec["flavours"]
    fn_name, fn_args, terms = parse_map(mapfile)
    unknown = [p for p in flavours if p not in fn_args]
    if unknown:
        raise SystemExit(f"ERROR: spec flavours {unknown} are not FN "
                         f"arguments {fn_args}")
    sink = {"errors": [], "warnings": [], "info": []}
    out = []
    for aN, term in terms:
        lc = PL.LineCheck(aN, term, flavours, sheet, sink,
                          subleading=(str(spec.get("colour", "")).lower()
                                      == "subleading"))
        if lc.antennae:
            out.append(lc)
    return out, flavours


def weight(lc):
    return lc.term["sign"] * lc.term["coeff"]


def x40_singles(lc, valid_cluster):
    """Spurious single poles of a single-X40 line (ledger's rule)."""
    prim4 = [(t, a, e) for (t, a, e) in lc.primary() if len(a) == 4]
    if not prim4 or len(lc.antennae) != 1:
        return None, []
    token, args, entry = prim4[0]
    meas = entry.get("measured", {})
    singles = []
    for key in meas.get("sco", {}):
        p, q = (int(v) for v in key.split(","))
        pair = frozenset([norm(args[p - 1]), norm(args[q - 1])])
        if all(isinstance(x, str) for x in pair) and valid_cluster(pair):
            singles.append(("sco", pair))
    for bslot in meas.get("ss", {}):
        leaf = norm(args[int(bslot) - 1])
        if isinstance(leaf, str) and lc.flavours.get(leaf) == "g":
            singles.append(("ss", frozenset([leaf])))
    return token, singles


def polename(kind, inv):
    if kind == "sco":
        return "s(" + ",".join(sorted(inv)) + ")"
    if kind == "ss":
        return "soft " + next(iter(inv))
    return f"{kind}(" + ",".join(sorted(inv)) + ")"


def analyse(lines, flavours):
    valid_cluster = make_valid_cluster(flavours)
    x40, iterated = [], []
    for lc in lines:
        tok, singles = x40_singles(lc, valid_cluster)
        if tok:
            x40.append((lc, tok, singles))
        elif len(lc.antennae) >= 2:
            iterated.append((lc, lc.poles()))

    poles = {}
    for lc, tok, singles in x40:
        for p in singles:
            poles.setdefault(p, {"x40": [], "itr": []})["x40"].append(lc)
    # SINGLE-unresolved poles only.  ds/tc poles of an X40 are its
    # legitimate job, not a spurious single, and are not collected on the
    # x40 side; letting them in on the iterated side alone manufactures
    # ORPHANs on a correct term (measured on a validated term).
    for lc, pl in iterated:
        for p in (q for q in pl if q[0] in ("sco", "ss")):
            if p in poles:
                poles[p]["itr"].append(lc)
            else:
                poles.setdefault(p, {"x40": [], "itr": []})["itr"].append(lc)
    return poles


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mapfile")
    ap.add_argument("--spec", required=True)
    ap.add_argument("--datasheet", default=PL.DEFAULT_DATASHEET)
    ap.add_argument("--tol", type=float, default=1e-6)
    ap.add_argument("--strict", action="store_true")
    ap.add_argument("--quiet", action="store_true",
                    help="print only poles that are flagged")
    ap.add_argument("--baseline", metavar="GOOD.map",
                    help="a term known correct; report per-pole DIFFERENCES "
                         "against it, which is the only reading with an "
                         "absolute meaning (see MEASURED SCOPE)")
    a = ap.parse_args()

    spec = json.load(open(a.spec))
    sheet = json.load(open(a.datasheet))
    lines, flavours = build_lines(a.mapfile, spec, sheet)
    poles = analyse(lines, flavours)
    if not poles:
        print("no X40 line with a spurious single pole, and no iterated "
              "line: nothing to balance")
        return 0

    errs, warns = [], []
    rows = []
    for p in sorted(poles, key=lambda k: (k[0], sorted(k[1]))):
        d = poles[p]
        sx = sum(weight(lc) for lc in d["x40"])
        si = sum(weight(lc) for lc in d["itr"])
        tot = sx + si
        tag = ""
        nm = polename(*p)
        if abs(sx) > a.tol and abs(si) <= a.tol:
            tag = "DANGLING"
            errs.append(f"{nm}: X40 weight {sx:+g} with no iterated "
                        f"counterterm — nothing can cancel it")
        elif abs(sx) > a.tol and sx * si > 0:
            tag = "SAME-SIGN"
            errs.append(f"{nm}: X40 {sx:+g} and iterated {si:+g} have the "
                        f"same sign — the counterterm ADDS to the pole")
        elif abs(sx) <= a.tol and abs(si) > a.tol:
            tag = "ORPHAN"
            errs.append(f"{nm}: iterated weight {si:+g} but no X40 carries "
                        f"this pole — over-subtracts a limit nothing "
                        f"restores")
        elif abs(tot) > a.tol:
            tag = "unbalanced"
            (errs if a.strict else warns).append(
                f"{nm}: x40 {sx:+g} + iterated {si:+g} = {tot:+g}")
        rows.append((nm, sx, si, tot, tag,
                     [lc.aN for lc in d["x40"]],
                     [lc.aN for lc in d["itr"]]))

    print(f"{'pole':<18}{'x40':>7}{'iter':>7}{'total':>8}  "
          f"{'flag':<11}lines")
    for nm, sx, si, tot, tag, ax, ai in rows:
        if a.quiet and not tag:
            continue
        lab = ("X40 " + ",".join("a%d" % n for n in sorted(ax))
               if ax else "X40 -")
        lab += "  |  " + ("itr " + ",".join("a%d" % n for n in sorted(ai))
                          if ai else "itr -")
        print(f"{nm:<18}{sx:>+7g}{si:>+7g}{tot:>+8g}  {tag:<11}{lab}")

    if a.baseline:
        blines, bfl = build_lines(a.baseline, spec, sheet)
        bp = analyse(blines, bfl)
        base = {}
        for p, d in bp.items():
            base[p] = (sum(weight(l) for l in d["x40"])
                       + sum(weight(l) for l in d["itr"]))
        print(f"\n  vs baseline {os.path.basename(a.baseline)}")
        diffs = 0
        for p in sorted(set(poles) | set(base),
                        key=lambda k: (k[0], sorted(k[1]))):
            here = poles.get(p)
            t = 0.0 if here is None else (
                sum(weight(l) for l in here["x40"])
                + sum(weight(l) for l in here["itr"]))
            b = base.get(p, 0.0)
            if abs(t - b) > a.tol:
                diffs += 1
                print(f"    {polename(*p):<18} this {t:+g}  baseline "
                      f"{b:+g}   delta {t - b:+g}")
                errs.append(f"{polename(*p)}: differs from baseline by "
                            f"{t - b:+g}")
        if not diffs:
            print("    identical per-pole sums — no count/coefficient "
                  "error relative to the baseline\n    (argument errors "
                  "are invisible here: see MEASURED SCOPE)")

    print()
    for e in errs:
        print("ERROR   " + e)
    for w in warns:
        print("warning " + w)
    print(f"\n{len(errs)} error(s), {len(warns)} warning(s)"
          + ("" if a.strict else
             "   [--strict turns imbalance into an error]"))
    return 1 if errs else 0


def selftest():
    """Structural only — a synthetic datasheet, no physics answers."""
    import tempfile
    sheet = {
        "XX30": {"arity": 3, "species": ["Q", "g", "Q"],
                 "measured": {"sco": {"1,2": {}, "2,3": {}},
                              "ss": {"2": {}}}},
        "YY40": {"arity": 4, "species": ["Q", "g", "g", "Q"],
                 "measured": {"sco": {"2,3": {}}, "ss": {}}},
    }
    spec = {"flavours": {"a": "q1", "b": "g", "c": "g", "d": "qb1"}}

    def run(body, extra=()):
        d = tempfile.mkdtemp()
        mp = os.path.join(d, "T.map")
        open(mp, "w").write("FN:=T(a,b,c,d):\nXX:=\n" + body + "\n:\n")
        lines, fl = build_lines(mp, spec, sheet)
        return analyse(lines, fl)

    # one X40 (+1) with pole s(b,c); one iterated line (-1) sharing it
    p = run("+YY40(a,b,c,d)*M(a,d)*a1\n"
            "-XX30(a,b,c)*XX30(a,c,d)*M(a,d)*a2\n")
    key = ("sco", frozenset(["b", "c"]))
    assert key in p, list(p)
    sx = sum(weight(l) for l in p[key]["x40"])
    si = sum(weight(l) for l in p[key]["itr"])
    assert sx == 1 and si == -1 and abs(sx + si) < 1e-9, (sx, si)

    # halve the counterterm -> unbalanced but not flagged as an error kind
    p = run("+YY40(a,b,c,d)*M(a,d)*a1\n"
            "-1/2*XX30(a,b,c)*XX30(a,c,d)*M(a,d)*a2\n")
    sx = sum(weight(l) for l in p[key]["x40"])
    si = sum(weight(l) for l in p[key]["itr"])
    assert abs(sx + si - 0.5) < 1e-9, (sx, si)

    # same-sign counterterm must be detectable
    p = run("+YY40(a,b,c,d)*M(a,d)*a1\n"
            "+XX30(a,b,c)*XX30(a,c,d)*M(a,d)*a2\n")
    sx = sum(weight(l) for l in p[key]["x40"])
    si = sum(weight(l) for l in p[key]["itr"])
    assert sx > 0 and si > 0, (sx, si)

    # X40 alone -> dangling
    p = run("+YY40(a,b,c,d)*M(a,d)*a1\n")
    assert p[key]["itr"] == []
    assert polename("ss", frozenset(["b"])) == "soft b"
    print("pairing_balance selftest OK")


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        selftest()
    else:
        sys.exit(main())
