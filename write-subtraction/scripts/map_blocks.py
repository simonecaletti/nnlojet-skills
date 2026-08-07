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

Usage:
  map_blocks.py list    master.map
  map_blocks.py compose master.map --blocks Sa,Sb1,Sb2 -o out.map
  map_blocks.py renumber file.map  -o out.map     # only fix aN gaps
  map_blocks.py --selftest

The composer is purely structural: it never invents, reorders, or
edits physics content — it selects lines you wrote and renumbers
labels.
"""
import argparse
import re
import sys

MARKER = re.compile(r"^\s*#\s*block:\s*(\S+)")
LABEL = re.compile(r"\*a(\d+)\b")


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
    print("map_blocks selftest OK")


def main():
    if "--selftest" in sys.argv:
        selftest()
        return
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["list", "compose", "renumber"])
    ap.add_argument("mapfile")
    ap.add_argument("--blocks", default="")
    ap.add_argument("-o", "--output")
    args = ap.parse_args()
    text = open(args.mapfile).read()
    if args.cmd == "list":
        cmd_list(text)
        return
    if args.cmd == "compose":
        sel = [s for s in args.blocks.split(",") if s]
        out, n = compose(text, sel)
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
