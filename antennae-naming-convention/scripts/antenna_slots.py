#!/usr/bin/env python3
"""Antenna slot-order audit: which entry points can be called
positionally from a .map, and which need a wrapper first.

The maple cluster rule is positional: X40(a,b,c,d) clusters
[a,b,c],[d,c,b]; X30(a,b,c) clusters [a,b],[b,c] (getpmapIK in
maple/getpmap.map). The Fortran entry points declare their momentum
slots with canonical names i1..i9 — but NOT all declare them in
ascending positional order (reversed and fully permuted orders coexist
in the tree). A .map line written positionally against a permuted
declaration evaluates the right function with the WRONG momentum map,
silently.

This script is the cheap FIRST lookup (before the pole scan, not
instead of it):
  - name grammar  -> species                (antennae-naming-convention)
  - THIS SCRIPT   -> slot plumbing          (declaration vs positional call)
  - pole scan     -> convention             (which slots radiate; probe-me-ir-structure)

Usage:
  antenna_slots.py [--root <repo root>] [--name <antenna substring>] [--all]
  antenna_slots.py --selftest

Output:
  - table of entry points grouped by slot-declaration order;
  - a FLAG on every entry point whose slot order is not ascending —
    i.e. a positional maple call clusters different legs than the
    dummy names imply — with a ready wrapper recipe;
  - cross-check against maple/notation.map: tokens in
    ant30set/ant31set/ant40set with no Fortran entry point, and
    antenna-like entry points not registered in any ant*fortset.

A flag means "check before positional use", not "wrong": some species'
canonical chains are legitimately non-consecutive. The pole scan
decides; the wrapper makes the positional call unambiguous either way.
"""
import argparse
import os
import re
import sys

FUNC = re.compile(
    r"^\s{0,6}(?:double\s+precision\s+)?function\s+(\w+)\s*\(([^)]*)\)",
    re.IGNORECASE | re.MULTILINE)
SLOT = re.compile(r"^i[1-9]$")
SETDEF = re.compile(r"^\s*(\w+set)\s*:=\s*(.*?):\s*$",
                    re.MULTILINE | re.DOTALL)
ANTENNA_LIKE = re.compile(r"(30|31|40)")

SRC_DIRS = ["src/X30", "src/X31", "src/X40",
            "src/X30int/FF", "src/X30int/IF",
            "src/X30int/FI", "src/X30int/II", "src/X30int"]


# ---------------- Fortran side ----------------

def parse_decls(text):
    """-> [(name, [args])] for every function declaration in text."""
    out = []
    for m in FUNC.finditer(text):
        args = [a.strip() for a in m.group(2).split(",") if a.strip()]
        out.append((m.group(1), args))
    return out


def slot_info(args):
    """-> (slots, extras, ascending) — slots are the i1..i9 args in
    positional order; ascending = positional order matches name order."""
    slots = [a for a in args if SLOT.match(a)]
    extras = [a for a in args if not SLOT.match(a)]
    nums = [int(a[1:]) for a in slots]
    return slots, extras, nums == sorted(nums)


def wrapper_recipe(name, args):
    """One-line wrapper putting the positional order in ascending slot
    order, plus its two registrations."""
    slots, extras, _ = slot_info(args)
    canon = sorted(slots, key=lambda a: int(a[1:]))
    head = ",".join(canon + extras)
    body = ",".join(args)
    fam = "X40" if "40" in name else ("X31" if "31" in name else "X30")
    return (
        f"      function {name}w({head})\n"
        f"      implicit double precision (a-h,o-z)\n"
        f"      {name}w = {name}({body})\n"
        f"      return\n"
        f"      end\n"
        f"c  register: (1) add the wrapper's file to the {fam} source list\n"
        f"c      in NNLOJET.mk;  (2) add the token to the matching ant*set\n"
        f"c      in maple/notation.map AND '<token>={name}w' to its\n"
        f"c      ant*fortset, so makefort* resolves it.\n")


# ---------------- notation.map side ----------------

def parse_notation(text):
    """Parse maple set assignments; resolve unions.
    -> (members, fortmap): members[name] = set of tokens;
       fortmap = {token: fortran_name} from the *fortset definitions."""
    raw = {}
    # capture 'name := { ... }:' including multiline braces, and unions
    for m in re.finditer(r"(\w+)\s*:=\s*([^:]*?):", text, re.DOTALL):
        name, rhs = m.group(1), m.group(2).strip()
        if not name.endswith("set"):
            continue
        raw[name] = rhs

    def resolve(name, seen=None):
        seen = seen or set()
        if name in seen or name not in raw:
            return set()
        seen.add(name)
        rhs = raw[name]
        if "{" in rhs:
            inner = rhs[rhs.index("{") + 1:rhs.rindex("}")]
            return {t.strip() for t in inner.replace("\n", ",").split(",")
                    if t.strip()}
        # union chain
        out = set()
        for part in re.split(r"\bunion\b", rhs):
            out |= resolve(part.strip(), seen)
        return out

    members = {n: resolve(n) for n in raw}
    fortmap = {}
    for n in raw:
        if not n.endswith("fortset"):
            continue
        for entry in members[n]:
            if "=" in entry:
                tok, fort = entry.split("=", 1)
                fortmap[tok.strip()] = fort.strip()
    return members, fortmap


# ---------------- audit ----------------

def audit(root, name_filter=None, show_all=False, out=print):
    decls = []
    for d in SRC_DIRS:
        p = os.path.join(root, d)
        if not os.path.isdir(p):
            continue
        for f in sorted(os.listdir(p)):
            if not f.endswith(".f"):
                continue
            try:
                text = open(os.path.join(p, f), errors="replace").read()
            except OSError:
                continue
            for name, args in parse_decls(text):
                decls.append((name, args, os.path.join(d, f)))
    if name_filter:
        decls = [d for d in decls if name_filter.lower() in d[0].lower()]

    groups = {}
    flagged = []
    for name, args, path in decls:
        slots, extras, asc = slot_info(args)
        if not slots:
            continue
        key = ",".join(slots)
        groups.setdefault(key, []).append(name)
        if not asc:
            flagged.append((name, args, path))

    out("=== slot-declaration orders (positional) ===")
    for key in sorted(groups, key=lambda k: (len(k.split(",")), k)):
        names = groups[key]
        mark = "" if slot_info(key.split(","))[2] else "   <-- PERMUTED"
        shown = names if (show_all or len(names) <= 8) \
            else names[:8] + [f"... +{len(names) - 8}"]
        out(f"  ({key}){mark}: {', '.join(shown)}")

    out(f"\n=== flagged: positional call clusters differently than the "
        f"dummy names imply ({len(flagged)}) ===")
    for name, args, path in flagged:
        out(f"\n  {name}({','.join(args)})   [{path}]")
        out("  wrapper recipe:")
        for ln in wrapper_recipe(name, args).splitlines():
            out("    " + ln)

    # cross-check notation.map
    npath = os.path.join(root, "maple", "notation.map")
    if os.path.isfile(npath):
        members, fortmap = parse_notation(open(npath,
                                               errors="replace").read())
        tokens = set()
        for s in ("ant30set", "ant31set", "ant40set"):
            tokens |= members.get(s, set())
        declnames = {n for n, _, _ in decls}
        missing_fort = sorted(
            t for t in tokens
            if fortmap.get(t, f"Full{t}") not in declnames
            and t not in declnames)
        registered = set(fortmap.values()) | set(fortmap.keys())
        unregistered = sorted(
            n for n in declnames
            if ANTENNA_LIKE.search(n) and n not in registered)
        out(f"\n=== notation.map tokens with no Fortran entry point "
            f"({len(missing_fort)}) ===")
        out("  " + (", ".join(missing_fort) or "none"))
        out(f"\n=== antenna-like entry points not registered in any "
            f"ant*fortset ({len(unregistered)}) ===")
        out("  " + (", ".join(unregistered) or "none"))
    else:
        out("\n(notation.map not found — cross-check skipped)")


# ---------------- selftest ----------------

def selftest():
    """Synthetic sources and notation text — encodes no real antenna's
    order. Checks declaration parsing, slot extraction, the ascending
    test, wrapper emission and the notation cross-check."""
    src = (
        "      function AAA40(i1,i2,i3,i4,ipset)\n"
        "      end\n"
        "      function BBB40(i1,i3,i4,i2,ipset)\n"
        "      end\n"
        "      double precision function CCC31(i3,i2,i1,ipset,renscale)\n"
        "      end\n"
        "      function helper(x,y)\n"
        "      end\n")
    decls = parse_decls(src)
    assert [d[0] for d in decls] == ["AAA40", "BBB40", "CCC31", "helper"]
    s, e, asc = slot_info(decls[0][1])
    assert s == ["i1", "i2", "i3", "i4"] and e == ["ipset"] and asc
    s, e, asc = slot_info(decls[1][1])
    assert not asc and s == ["i1", "i3", "i4", "i2"]
    s, e, asc = slot_info(decls[2][1])
    assert not asc and e == ["ipset", "renscale"]
    assert slot_info(decls[3][1])[0] == []          # no slots -> skipped

    w = wrapper_recipe("BBB40", decls[1][1])
    assert "function BBB40w(i1,i2,i3,i4,ipset)" in w
    assert "BBB40w = BBB40(i1,i3,i4,i2,ipset)" in w
    w = wrapper_recipe("CCC31", decls[2][1])
    assert "function CCC31w(i1,i2,i3,ipset,renscale)" in w
    assert "CCC31w = CCC31(i3,i2,i1,ipset,renscale)" in w

    notation = (
        "antXXFFset:={TOKA, TOKB}:\n"
        "antXXIFset:={TOKC}:\n"
        "ant40set:= antXXFFset union antXXIFset:\n"
        "ant30set:={}:\n"
        "ant31set:={}:\n"
        "XXfortset:={\nTOKA=AAA40,\nTOKB=BBB40w\n}:\n"
        "ant40fortset:=XXfortset:\n")
    members, fortmap = parse_notation(notation)
    assert members["ant40set"] == {"TOKA", "TOKB", "TOKC"}, members
    assert fortmap == {"TOKA": "AAA40", "TOKB": "BBB40w"}
    # TOKC has no fortmap entry and no FullTOKC declaration -> missing;
    # CCC31 is antenna-like and unregistered
    lines = []
    import tempfile
    with tempfile.TemporaryDirectory() as root:
        os.makedirs(os.path.join(root, "src", "X40"))
        os.makedirs(os.path.join(root, "maple"))
        open(os.path.join(root, "src", "X40", "t.f"), "w").write(src)
        open(os.path.join(root, "maple", "notation.map"), "w") \
            .write(notation)
        audit(root, out=lines.append)
    txt = "\n".join(lines)
    assert "PERMUTED" in txt
    assert "BBB40w = BBB40(i1,i3,i4,i2,ipset)" in txt
    assert "TOKC" in txt.split("no Fortran entry point")[1]
    assert "CCC31" in txt.split("not registered")[1]
    print("antenna_slots selftest OK")


def main():
    if "--selftest" in sys.argv:
        selftest()
        return
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".")
    ap.add_argument("--name", help="filter by antenna-name substring")
    ap.add_argument("--all", action="store_true",
                    help="do not abbreviate group listings")
    args = ap.parse_args()
    audit(args.root, args.name, args.all)


if __name__ == "__main__":
    main()
