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
import sys

from _block_structure import LABEL
from predict_blocks import BLOCK_ORDER, predict, predict_from_files


def emit_markers(pred, fn="FN_PLACEHOLDER"):
    """A compose_blocks.py-composable skeleton: one placeholder term per
    predicted line, grouped under '# block:' markers, aN gap-free. The
    label is placed BEFORE the trailing comment so it survives as a
    label, not as commented-out text."""
    body, i = [], 0
    for b in BLOCK_ORDER:
        sel = [l for l in pred["lines"] if l["block"] == b]
        if not sel:
            continue
        body.append(f"# block: {b}")
        for ln in sel:
            i += 1
            stem = (ln["family_hint"].replace("*", "x").replace(" ", "")
                    .replace("-", "_"))
            tag = "_".join(map(str, ln["unresolved"]))
            body.append(f"{ln['sign']}TODO_{stem}_{tag}*a{i}"
                        f"    # rad={ln['radiators']} {ln['config']} "
                        f"coeff~{ln['coeff_hint']}")
    blocks = ",".join(b for b in BLOCK_ORDER if b in pred["counts"])
    return "\n".join(
        ["# PREDICTED SKELETON — placeholders, NOT a valid .map.",
         "# Replace each TODO_* with antenna(args)*reducedME(args)*JETnm(args),",
         f"# then: compose_blocks.py compose <this> --blocks {blocks} -o TERM.map"]
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
    out = emit_markers(pred)
    if args.output:
        with open(args.output, "w") as fh:
            fh.write(out)
        print(f"wrote {args.output}: {len(pred['lines'])} placeholder line(s)")
    else:
        sys.stdout.write(out)


if __name__ == "__main__":
    main()
