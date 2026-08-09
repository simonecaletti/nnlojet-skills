#!/usr/bin/env python3
"""Block composer for NNLOJET maple subtraction files (.map).

Represents a .map as named blocks of term lines, composes a file from a
chosen subset of blocks, and renumbers the aN labels automatically so
they stay gap-free (the generator aborts on gaps). This turns "does
this block belong here?" into a single short test instead of a
hand-renumbering session.

Master-file convention (the '# block:' / '# axis:' format itself lives
in _block_structure.py, shared with audit_blocks.py and
emit_skeleton.py — emit_skeleton.py writes a starting master for you):
inside the XX:= body, group term lines under maple comment markers

    # block: <name>

Everything before the first marker belongs to block "_default".
Everything outside the XX:= body (FN line, header comments, the
terminating ':') is structural and always preserved.

VARIANT AXES. Toggling whole blocks explores only the line set you
already wrote. The choices INSIDE a line — which pair member a radiator
absorbs, which antenna of a species, which coefficient — are just as
free, and they are NOT independent: flipping an S,a cluster forces the
matching flip in every S,b2 line built on it. Declare such a choice as
an AXIS whose options name whole alternative block sets, so co-dependent
lines can only move together:

    # axis: absorb = near:Sa_f1n,Sb2_f1Gn,Sb2_f1mn | far:Sa_f1f,Sb2_f1Gf,Sb2_f1mf
    # axis: sb1    = E40a:Sb1_f1,Sb1_f2 | full:Sb1_f1F,Sb1_f2F | none:

An option with an empty block list ("none:") means "omit this axis".
`enumerate` then emits the cross-product as ready-to-compose block
lists — the input scan_blocks.py (run-spike-test) sweeps unattended.

Usage:
  compose_blocks.py list      master.map
  compose_blocks.py axes      master.map
  compose_blocks.py enumerate master.map --fixed Sa_g --axes absorb,sb1
  compose_blocks.py compose   master.map --blocks Sa,Sb1,Sb2 -o out.map
  compose_blocks.py compose   master.map --select absorb=near,sb1=E40a \
                              --fixed Sa_g -o out.map
  compose_blocks.py renumber  file.map  -o out.map     # only fix aN gaps
  compose_blocks.py --selftest

The composer is purely structural: it never invents, reorders, or
edits physics content — it selects lines you wrote and renumbers
labels.
"""
import argparse
import itertools
import sys

from _block_structure import (LABEL, split_map, blocks_of, renumber, compose,
                              axes_of, check_axes, expand)


def cmd_axes(text):
    ax = check_axes(text)
    if not ax:
        print("  (no '# axis:' declarations in this master)")
        return
    total = 1
    for name, opts in ax.items():
        total *= len(opts)
        print(f"  {name:12s} {len(opts)} option(s)")
        for opt, blks in opts:
            print(f"      {opt:10s} -> {','.join(blks) if blks else '(none)'}")
    print(f"  cross-product: {total} combination(s)")


def cmd_enumerate(text, fixed, want):
    ax = check_axes(text)
    names = want or list(ax)
    for n in names:
        if n not in ax:
            raise SystemExit(f"ERROR: unknown axis '{n}'")
    for combo in itertools.product(*[[o for o, _ in ax[n]] for n in names]):
        sel = dict(zip(names, combo))
        label = "|".join(f"{k}={v}" for k, v in sel.items())
        print(f"{label}\t{','.join(expand(text, fixed, sel))}")


def cmd_list(text):
    _, body, _ = split_map(text)
    for n, ls in blocks_of(body):
        nt = sum(len(LABEL.findall(ln)) for ln in ls)
        print(f"  {n:24s} {nt:3d} term(s), {len(ls)} line(s)")


def selftest():
    """Structural self-test — encodes no physics. Checks gap-free
    renumbering, subset selection, idempotence, header preservation."""
    master = "\n".join([
        "# header comment",
        "FN:=DUMMY(1,2,3,4):",
        "XX:=",
        "# block: Sa",
        "+ANT1(3,4,5)*RED1([3,4],[4,5])*JET22([3,4],[4,5])*a1",
        "+ANT2(5,4,3)*RED2([5,4],[4,3])*JET22([5,4],[4,3])*a2",
        "# block: Sb1",
        "+ANT4(3,4,5,6)*RED3([3,4,5],[6,5,4])*JET22([3,4,5],[6,5,4])*a3",
        "# block: Sb2",
        "-ANT1(3,4,5)*ANT1([3,4],[4,5],6)*RED3(...)*JET22(...)*a4",
        "-1/2*ANT2(5,4,3)*ANT1([5,4],[4,3],6)*RED3(...)*JET22(...)*a5",
        ":",
    ]) + "\n"
    # full composition is idempotent and gap-free
    full, n = compose(master, ["Sa", "Sb1", "Sb2"])
    assert n == 5
    labels = [int(x) for x in LABEL.findall(full)]
    assert labels == list(range(1, 6)), labels
    again, n2 = compose(full, ["Sa", "Sb1", "Sb2"])
    assert again == full and n2 == n, "not idempotent"
    # subset drops the right lines and renumbers without gaps
    sub, n = compose(master, ["Sa", "Sb2"])
    assert n == 4
    assert "ANT4" not in sub, "unselected block leaked"
    labels = [int(x) for x in LABEL.findall(sub)]
    assert labels == [1, 2, 3, 4], labels
    assert "FN:=DUMMY(1,2,3,4):" in sub, "header lost"
    assert sub.rstrip().endswith(":"), "terminator lost"
    # unknown block refused
    try:
        compose(master, ["nope"])
        raise AssertionError("unknown block accepted")
    except SystemExit:
        pass
    # zero-term composition refused
    try:
        compose(master, [])
        raise AssertionError("empty composition accepted")
    except SystemExit:
        pass

    # --- variant axes (structural only: no physics encoded) ---------
    vm = "\n".join([
        "# axis: pick = A:Sa | B:Sb1 | none:",
        "# axis: ct   = on:Sb2 | off:",
        "FN:=DUMMY(1,2,3,4):",
        "XX:=",
        "# block: Sa",
        "+ANT1(3,4,5)*RED1([3,4],[4,5])*JET22([3,4],[4,5])*a1",
        "# block: Sb1",
        "+ANT4(3,4,5,6)*RED3([3,4,5],[6,5,4])*JET22([3,4,5],[6,5,4])*a2",
        "# block: Sb2",
        "-ANT1(3,4,5)*ANT1([3,4],[4,5],6)*RED3(...)*JET22(...)*a3",
        ":",
    ]) + "\n"
    ax = check_axes(vm)
    assert list(ax) == ["pick", "ct"], ax
    assert len(ax["pick"]) == 3 and len(ax["ct"]) == 2
    # expansion honours fixed blocks, dedupes, and drops empty options
    assert expand(vm, [], {"pick": "A", "ct": "on"}) == ["Sa", "Sb2"]
    assert expand(vm, ["Sa"], {"pick": "A"}) == ["Sa"], "duplicate leaked"
    assert expand(vm, [], {"pick": "none", "ct": "off"}) == []
    # a composed selection equals the equivalent explicit --blocks
    a, _ = compose(vm, expand(vm, [], {"pick": "B", "ct": "on"}))
    b, _ = compose(vm, ["Sb1", "Sb2"])
    assert a == b, "axis expansion != explicit block list"
    # unknown axis / option refused
    for bad in ({"nope": "A"}, {"pick": "Z"}):
        try:
            expand(vm, [], bad)
            raise AssertionError(f"bad selection accepted: {bad}")
        except SystemExit:
            pass
    # an axis naming a block that does not exist is refused
    try:
        check_axes(vm.replace("A:Sa |", "A:Ghost |"))
        raise AssertionError("axis with unknown block accepted")
    except SystemExit:
        pass
    # a single-option axis is refused (nothing to vary)
    try:
        axes_of("# axis: solo = only:Sa\n")
        raise AssertionError("single-option axis accepted")
    except SystemExit:
        pass
    print("compose_blocks selftest OK")


def main():
    if "--selftest" in sys.argv:
        selftest()
        return
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd",
                    choices=["list", "axes", "enumerate", "compose",
                             "renumber"])
    ap.add_argument("mapfile")
    ap.add_argument("--blocks", default="")
    ap.add_argument("--fixed", default="",
                    help="blocks always included (axes add to these)")
    ap.add_argument("--axes", default="",
                    help="enumerate: axes to cross-product (default: all)")
    ap.add_argument("--select", default="",
                    help="compose: axis=option,axis=option")
    ap.add_argument("-o", "--output")
    args = ap.parse_args()
    text = open(args.mapfile).read()
    fixed = [s for s in args.fixed.split(",") if s]
    if args.cmd == "list":
        cmd_list(text)
        return
    if args.cmd == "axes":
        cmd_axes(text)
        return
    if args.cmd == "enumerate":
        cmd_enumerate(text, fixed,
                      [s for s in args.axes.split(",") if s])
        return
    if args.cmd == "compose":
        if args.select:
            sel = {}
            for item in args.select.split(","):
                if "=" not in item:
                    raise SystemExit(
                        f"ERROR: --select item '{item}' is not axis=option")
                k, v = item.split("=", 1)
                sel[k.strip()] = v.strip()
            blocks = expand(text, fixed, sel)
        else:
            blocks = fixed + [s for s in args.blocks.split(",")
                              if s and s not in fixed]
        out, n = compose(text, blocks)
    else:  # renumber
        head, body, tail = split_map(text)
        new_body, n = renumber(body)
        out = "\n".join(head + new_body + tail) + "\n"
    if args.output:
        with open(args.output, "w") as fh:
            fh.write(out)
        print(f"wrote {args.output}: {n} term(s), labels a1..a{n}")
    else:
        sys.stdout.write(out)


if __name__ == "__main__":
    main()
