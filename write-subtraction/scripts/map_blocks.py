#!/usr/bin/env python3
"""Block composer for NNLOJET maple subtraction files (.map).

Represents a .map as named blocks of term lines, composes a file from a
chosen subset of blocks, and renumbers the aN labels automatically so
they stay gap-free (the generator aborts on gaps). This turns "does
this block belong here?" into a single short test instead of a
hand-renumbering session.

Master-file convention: inside the XX:= body, group term lines under
maple comment markers

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
  map_blocks.py list      master.map
  map_blocks.py axes      master.map
  map_blocks.py enumerate master.map --fixed Sa_g --axes absorb,sb1
  map_blocks.py compose   master.map --blocks Sa,Sb1,Sb2 -o out.map
  map_blocks.py compose   master.map --select absorb=near,sb1=E40a \
                          --fixed Sa_g -o out.map
  map_blocks.py renumber  file.map  -o out.map     # only fix aN gaps
  map_blocks.py --selftest

The composer is purely structural: it never invents, reorders, or
edits physics content — it selects lines you wrote and renumbers
labels.
"""
import argparse
import itertools
import re
import sys

MARKER = re.compile(r"^\s*#\s*block:\s*(\S+)")
AXIS = re.compile(r"^\s*#\s*axis:\s*(\w+)\s*=\s*(.+?)\s*$")
LABEL = re.compile(r"\*a(\d+)\b")


def axes_of(text):
    """-> ordered dict {axis: [(option, [blocks]), ...]} from '# axis:' lines."""
    out = {}
    for ln in text.splitlines():
        m = AXIS.match(ln)
        if not m:
            continue
        name, rhs = m.group(1), m.group(2)
        opts = []
        for chunk in rhs.split("|"):
            chunk = chunk.strip()
            if not chunk:
                continue
            if ":" not in chunk:
                raise SystemExit(
                    f"ERROR: axis '{name}' option '{chunk}' is not <opt>:<blocks>")
            opt, blks = chunk.split(":", 1)
            opts.append((opt.strip(),
                         [b.strip() for b in blks.split(",") if b.strip()]))
        if len(opts) < 2:
            raise SystemExit(f"ERROR: axis '{name}' needs >= 2 options")
        if name in out:
            raise SystemExit(f"ERROR: axis '{name}' declared twice")
        out[name] = opts
    return out


def check_axes(text):
    """Every block an axis names must exist. Returns the axis table."""
    _, body, _ = split_map(text)
    known = {n for n, _ in blocks_of(body)}
    ax = axes_of(text)
    for name, opts in ax.items():
        for opt, blks in opts:
            missing = [b for b in blks if b not in known]
            if missing:
                raise SystemExit(
                    f"ERROR: axis '{name}' option '{opt}' names unknown "
                    f"block(s): {','.join(missing)}")
    return ax


def expand(text, fixed, selection):
    """fixed: [blocks]; selection: {axis: option} -> ordered block list."""
    ax = check_axes(text)
    out = list(fixed)
    for name, opt in selection.items():
        if name not in ax:
            raise SystemExit(f"ERROR: unknown axis '{name}'")
        match = [b for o, b in ax[name] if o == opt]
        if not match:
            have = ",".join(o for o, _ in ax[name])
            raise SystemExit(
                f"ERROR: axis '{name}' has no option '{opt}' (have: {have})")
        for b in match[0]:
            if b not in out:
                out.append(b)
    return out


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


def split_map(text):
    """-> (head_lines, body_lines, tail_lines). Body = inside XX:= ... :"""
    lines = text.splitlines()
    try:
        istart = next(i for i, ln in enumerate(lines)
                      if ln.strip().replace(" ", "").startswith("XX:="))
    except StopIteration:
        raise SystemExit("ERROR: no 'XX:=' line found — not a term .map?")
    # terminating line: last line that is exactly ':' (possibly spaces)
    iend = None
    for i in range(len(lines) - 1, istart, -1):
        if lines[i].strip() == ":":
            iend = i
            break
    if iend is None:
        raise SystemExit("ERROR: no terminating ':' line after XX:=")
    return lines[:istart + 1], lines[istart + 1:iend], lines[iend:]


def blocks_of(body):
    """-> ordered list of (name, [lines])."""
    out = []
    cur_name, cur = "_default", []
    for ln in body:
        m = MARKER.match(ln)
        if m:
            out.append((cur_name, cur))
            cur_name, cur = m.group(1), [ln]
        else:
            cur.append(ln)
    out.append((cur_name, cur))
    # drop an empty _default
    return [(n, ls) for n, ls in out if n != "_default" or
            any(s.strip() for s in ls)]


def renumber(lines):
    """Renumber every *aN occurrence sequentially, gap-free from 1."""
    counter = [0]

    def sub(_m):
        counter[0] += 1
        return f"*a{counter[0]}"

    return [LABEL.sub(sub, ln) for ln in lines], counter[0]


def compose(text, selected):
    head, body, tail = split_map(text)
    blks = blocks_of(body)
    names = [n for n, _ in blks]
    unknown = [s for s in selected if s not in names]
    if unknown:
        raise SystemExit(f"ERROR: unknown block(s) {unknown}; "
                         f"available: {names}")
    kept = []
    for n, ls in blks:
        if n == "_default" or n in selected:
            kept.extend(ls)
    new_body, nterms = renumber(kept)
    if nterms == 0:
        raise SystemExit("ERROR: composition selects zero terms")
    return "\n".join(head + new_body + tail) + "\n", nterms


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
    print("map_blocks selftest OK")


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
