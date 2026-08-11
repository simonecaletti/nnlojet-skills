#!/usr/bin/env python3
"""Audit a written .map against the rung-0 block prediction.

STRUCTURAL ONLY: it checks which '# block:'-marked blocks exist and how
many terms each carries — MISSING blocks, EXTRA blocks, and line-count
mismatches. It cannot check antenna letters, arguments or mappings:
those are measured (probe-me-ir-structure), and the invariant-level
bookkeeping is pole_ledger.py's job. Order of operations:
audit_blocks.py (counts) -> pole_ledger.py (invariants) -> generate ->
spike test (numbers).

Counts are per-ordering upper bounds, so a colour-summed term may
legitimately differ by a factor equal to the number of orderings —
investigate any other discrepancy before generating.

Usage:
  python audit_blocks.py spec.json TERM.map [--modes modes.json]
  python audit_blocks.py --selftest

Spec and --modes are exactly predict_blocks.py's (one spec file drives
predict / emit / audit alike); prediction warnings go to stderr.
Exit status: 1 on any MISSING/COUNT/EXTRA finding, else 0.
"""
import argparse
import re
import sys

from _block_structure import LABEL, MARKER
from predict_blocks import BLOCK_ORDER, predict, predict_from_files


def audit(pred, text):
    """Compare a written .map against the prediction."""
    actual, cur = {}, None
    for ln in text.splitlines():
        m = MARKER.match(ln)
        if m:
            cur = m.group(1)
            actual.setdefault(cur, 0)
        elif cur is not None:
            actual[cur] += len(LABEL.findall(ln))
    if not actual:
        return ["no '# block:' markers found — the file is not composed from "
                "a block master, so the audit cannot run."]
    stems = {k.split("_")[0] for k in actual}
    rep = []
    for b in BLOCK_ORDER:
        if b not in pred["counts"]:
            continue
        n = pred["counts"][b]
        if b not in stems:
            rep.append(f"MISSING  {b}: predicted {n} line(s), block absent")
            continue
        got = sum(v for k, v in actual.items() if k.split("_")[0] == b)
        if got != n:
            rep.append(f"COUNT    {b}: predicted {n}, file has {got} "
                       "(per-ordering counts are upper bounds — a colour-"
                       "summed term may legitimately differ)")
        else:
            rep.append(f"ok       {b}: {n}")
    # The documented S,a anti-pattern: averaging over BOTH radiators and BOTH
    # reduced-ME orderings with 1/2 factors.  It reproduces the collinear limit
    # (so the single-unresolved modes read 1.000000 and look like a pass) and
    # matches no X40 — every downstream pairing then fails for reasons that
    # look like they belong to the block being edited.  Detectable statically.
    cur, tot, half = None, {}, {}
    for ln in text.splitlines():
        m = MARKER.match(ln)
        if m:
            cur = m.group(1)
            continue
        if cur and cur.split("_")[0].startswith("Sa"):
            n = len(LABEL.findall(ln))
            if n:
                tot[cur] = tot.get(cur, 0) + n
                if re.search(r"(?<![\w/])1/2\s*\*", ln):
                    half[cur] = half.get(cur, 0) + n
    for b in sorted(tot):
        if tot[b] >= 4 and half.get(b, 0) == tot[b]:
            rep.append(
                f"RADIATOR {b}: {tot[b]} lines, ALL at coefficient 1/2 — the "
                "'average over radiators and orderings with 1/2 factors' "
                "anti-pattern (SKILL.md): it reproduces the collinear limit "
                "and matches no X40.  Sweep the radiator axis (G30-radiator, "
                "coefficient 1, one line per ordering vs E30-radiator) BEFORE "
                "building S,b1/S,b2 on it.")
    for k in sorted(actual):
        if k != "_default" and k.split("_")[0] not in pred["counts"]:
            rep.append(f"EXTRA    {k}: present but not predicted — either a "
                       "second colour ordering, or a spurious block")
    return rep


def selftest():
    """Audit invariants only (prediction invariants live in
    predict_blocks.py's selftest, skeleton invariants in
    emit_skeleton.py's)."""
    from emit_skeleton import emit_markers
    spec = {"chain": [1, 3, 4, 5, 2],
            "flavours": {"1": "q1", "3": "g", "4": "g", "5": "g", "2": "qb1"},
            "initial": [1, 2],
            "unresolvable": [3, 4, 5]}
    p = predict(spec)

    # the emitted skeleton audits clean against its own prediction
    skel = emit_markers(p)
    assert all(r.startswith("ok") for r in audit(p, skel)), audit(p, skel)

    # the audit detects under-count, missing blocks and extra blocks
    rep = audit(p, "# block: Sa\n+X*a1\n")
    assert any(r.startswith("COUNT") and " Sa:" in r for r in rep), rep
    assert any(r.startswith("MISSING") for r in rep), rep
    rep = audit(p, "# block: Sa\n+X*a1\n# block: Zz\n+Y*a2\n")
    assert any(r.startswith("EXTRA") for r in rep), rep
    assert audit(p, "no markers here")[0].startswith("no '# block:'")
    print("audit_blocks selftest OK")


def main():
    if "--selftest" in sys.argv:
        selftest()
        return
    ap = argparse.ArgumentParser()
    ap.add_argument("spec")
    ap.add_argument("mapfile", metavar="TERM.map")
    ap.add_argument("--modes", help="genuine_modes.py --json output")
    args = ap.parse_args()
    pred = predict_from_files(args.spec, args.modes)
    rep = audit(pred, open(args.mapfile).read())
    for r in rep:
        print(r)
    # a gate that cannot fail is not a gate
    raise SystemExit(0 if all(r.startswith("ok") for r in rep) else 1)


if __name__ == "__main__":
    main()
