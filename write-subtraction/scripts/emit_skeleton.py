#!/usr/bin/env python3
"""Emit a composable MASTER SKELETON from the rung-0 block prediction.

Writes one placeholder term per predicted line, grouped under the
'# block:' markers compose_blocks.py composes from, aN labels gap-free.
Replace each TODO_* with antenna(args)*reducedME(args)*JETnm(args), then
compose the subset under test.

Usage:
  python emit_skeleton.py spec.json [--modes modes.json] [-o master.map]
  python emit_skeleton.py --selftest

Spec and --modes are exactly predict_blocks.py's (one spec file drives
predict / emit / audit alike); prediction warnings go to stderr.
"""
import argparse
import json
import sys

from _block_structure import LABEL
from predict_blocks import BLOCK_ORDER, predict, predict_from_files


def secondary_pair(spec):
    """True when the channel carries two distinct quark flavours, i.e. a
    SECONDARY q qb pair.  Then 'which parton emits it' — the Born gluon
    (G30FF, coefficient 1, one line per reduced-ME ordering) or the quark
    pair (E30FF) — is a free but LOAD-BEARING choice: it fixes which X40
    family S,b1/S,b2 must use.  It is emitted as an axis so it cannot be
    hard-coded by accident."""
    tags = set()
    for f in (spec or {}).get("flavours", {}).values():
        if f.startswith("qb"):
            tags.add(f[2:])
        elif f.startswith("q"):
            tags.add(f[1:])
    return len(tags) >= 2


def emit_markers(pred, fn="FN_PLACEHOLDER", spec=None):
    """A compose_blocks.py-composable skeleton: one placeholder term per
    predicted line, grouped under '# block:' markers, aN gap-free. The
    label is placed BEFORE the trailing comment so it survives as a
    label, not as commented-out text."""
    split_sa = secondary_pair(spec)
    body, i, axis_g, axis_e = [], 0, [], []

    def emit(name, sel, fam=None, coeff=None):
        nonlocal i
        body.append(f"# block: {name}")
        for ln in sel:
            i += 1
            stem = ((fam or ln["family_hint"]).replace("*", "x")
                    .replace(" ", "").replace("-", "_"))
            tag = "_".join(map(str, ln["unresolved"]))
            body.append(f"{ln['sign']}TODO_{stem}_{tag}*a{i}"
                        f"    # rad={ln['radiators']} {ln['config']} "
                        f"coeff~{coeff or ln['coeff_hint']}")

    for b in BLOCK_ORDER:
        sel = [l for l in pred["lines"] if l["block"] == b]
        if not sel:
            continue
        if split_sa and b.startswith("Sa"):
            emit(f"{b}_g30", sel, fam="G30FF", coeff="1")
            emit(f"{b}_e30", sel, fam="E30FF", coeff="1")
            axis_g.append(f"{b}_g30")
            axis_e.append(f"{b}_e30")
        else:
            emit(b, sel)

    blocks = ",".join((f"{b}_g30" if (split_sa and b.startswith("Sa")) else b)
                      for b in BLOCK_ORDER if b in pred["counts"])
    axis = ([f"# axis: radiator = g30:{','.join(axis_g)} "
             f"| e30:{','.join(axis_e)}",
             "# !! The radiator choice above is an AXIS, not a default. Both "
             "options reproduce the single-unresolved limits exactly, so a "
             "mode reading 1.000000 does NOT confirm it; what it fixes is "
             "which X40 family S,b1/S,b2 can pair with. SWEEP IT."]
            if axis_g else [])
    return "\n".join(
        ["# PREDICTED SKELETON — placeholders, NOT a valid .map.",
         "# Replace each TODO_* with antenna(args)*reducedME(args)*JETnm(args),",
         f"# then: compose_blocks.py compose <this> --blocks {blocks} -o TERM.map"]
        + axis
        + [f"# !! {w}" for w in pred["warnings"]]
        + [f"FN:={fn}:", "XX:="] + body + [":"]) + "\n"


def selftest():
    """Skeleton invariants only (prediction invariants live in
    predict_blocks.py's selftest; the round-trip against the audit lives
    in audit_blocks.py's)."""
    spec = {"chain": [1, 3, 4, 5, 2],
            "flavours": {"1": "q1", "3": "g", "4": "g", "5": "g", "2": "qb1"},
            "initial": [1, 2],
            "unresolvable": [3, 4, 5]}
    p = predict(spec)
    skel = emit_markers(p)
    # emitted skeleton is gap-free
    labels = [int(x) for x in LABEL.findall(skel)]
    assert labels == list(range(1, len(p["lines"]) + 1)), labels
    # the label must not be swallowed by the trailing comment
    for ln in skel.splitlines():
        if "*a" in ln and "#" in ln:
            assert ln.index("*a") < ln.index("#"), ln
    # every predicted block appears as a marker, in canonical order
    markers = [ln.split(":", 1)[1].strip() for ln in skel.splitlines()
               if ln.startswith("# block:")]
    assert markers == [b for b in BLOCK_ORDER if b in p["counts"]], markers
    # prediction warnings are carried into the skeleton, visibly
    lin = predict({"chain": [1, 2, 3, 4, 5], "cyclic": False,
                   "flavours": {str(i): "g" for i in [1, 2, 3, 4, 5]},
                   "unresolvable": [1, 5]})
    assert all(f"# !! {w}" in emit_markers(lin) for w in lin["warnings"])
    # a channel with a SECONDARY quark pair must emit the radiator axis,
    # with both options present as real blocks
    two = {"chain": [1, 3, 4, 5, 2],
           "flavours": {"1": "q1", "3": "g", "4": "q2", "5": "qb2",
                        "2": "qb1"},
           "initial": [], "unresolvable": [3, 4, 5]}
    sk2 = emit_markers(predict(two), spec=two)
    assert "# axis: radiator = g30:" in sk2, sk2
    assert "_g30" in sk2 and "_e30" in sk2, sk2
    lab2 = [int(x) for x in LABEL.findall(sk2)]
    assert lab2 == list(range(1, len(lab2) + 1)), lab2
    # a single-flavour channel must NOT invent the axis
    assert "# axis: radiator" not in emit_markers(p, spec=spec)
    print("emit_skeleton selftest OK")


def main():
    if "--selftest" in sys.argv:
        selftest()
        return
    ap = argparse.ArgumentParser()
    ap.add_argument("spec")
    ap.add_argument("--modes", help="genuine_modes.py --json output")
    ap.add_argument("-o", "--output")
    args = ap.parse_args()
    pred = predict_from_files(args.spec, args.modes)
    with open(args.spec) as fh:
        spec = json.load(fh)
    out = emit_markers(pred, spec=spec)
    if args.output:
        with open(args.output, "w") as fh:
            fh.write(out)
        print(f"wrote {args.output}: {len(pred['lines'])} placeholder line(s)")
    else:
        sys.stdout.write(out)


if __name__ == "__main__":
    main()
