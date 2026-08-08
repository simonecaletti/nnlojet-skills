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
    dummy names imply — with an ARGUMENT-PERMUTATION recipe;
  - cross-check against maple/notation.map: tokens in
    ant30set/ant31set/ant40set with no Fortran entry point, and
    antenna-like entry points not registered in any ant*fortset;
  - a scan for UNREGISTERED FORWARDING WRAPPERS in src/X30|X31|X40
    (one-line functions renaming another antenna that no ant*fortset
    resolves): their headers carry unverified claims and rot silently
    — do not trust them, and do not call them from a .map.

A flag means "check before positional use", not "wrong": some species'
canonical chains are legitimately non-consecutive.

THE FIX IS A PERMUTED CALL, NEVER A WRAPPER (write-subtraction's
argument-alignment rule): permute the .map ARGUMENTS so the generic
positional cluster rule lands on the legs the antenna actually
clusters — e.g. the two halves of a Full X40 are used as
X40a(A,B,C,D) together with X40b(A,D,C,B).  A wrapper restoring
ascending slot order gives the RIGHT function value with the WRONG
clusters, and it looks correct because the value-level checks agree
with each other.  Which permutation is right is a physics question:
read it off the measured pole graph / verified split identity in the
antenna DATASHEET (antenna_datasheet.py show <name>), or run the pole
scan if the antenna is not in the datasheet yet.
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


def call_recipe(name, args):
    """Permuted-ARGUMENT call recipe for a permuted declaration: how a
    positional .map call realises the named-slot assignment, plus the
    datasheet pointer that decides whether that is the physically
    right permutation.  NO wrapper is emitted: a wrapper restoring
    ascending order gives the right VALUE with the WRONG clusters
    under the generic positional cluster rule (write-subtraction)."""
    slots, extras, _ = slot_info(args)
    canon = sorted(slots, key=lambda a: int(a[1:]))
    letters = "ABCDEFGHI"[:len(slots)]
    named = {s: letters[i] for i, s in enumerate(canon)}
    positional = "".join(named[s] for s in slots)
    return (
        f"  named legs ({','.join(canon)}) = ({','.join(letters)})  "
        f"are realised by the positional call\n"
        f"      {name}({','.join(positional)}"
        + ("," + ",".join(extras) if extras else "") + ")\n"
        f"  under which the generic cluster rule produces the "
        f"clusters of THAT leg\n"
        f"  assignment. Whether it is the right permutation for your "
        f"line is decided\n"
        f"  by the measured pole graph / split identity: "
        f"antenna_datasheet.py show {name}\n"
        f"  (pole scan via probe-me-ir-structure if not in the "
        f"datasheet).\n"
        f"  Do NOT write a slot-reordering wrapper: right value, "
        f"wrong clusters.\n")


# ---------------- notation.map side ----------------

def parse_notation(text):
    """Parse maple set assignments; resolve unions.
    -> (members, fortmap): members[name] = set of tokens;
       fortmap = {token: fortran_name} from the *fortset definitions."""
    # brace-balanced capture up to the first ':' at depth 0 — a colon
    # INSIDE a braced set (the LaTeX comment sets have them) must not
    # truncate the rhs, which used to crash resolve() on '}'-less text
    raw = {}
    for m in re.finditer(r"(\w+set)\s*:=", text):
        name = m.group(1)
        depth, k = 0, m.end()
        while k < len(text):
            ch = text[k]
            if ch in "{[(":
                depth += 1
            elif ch in "}])":
                depth -= 1
            elif ch == ":" and depth == 0:
                break
            k += 1
        raw[name] = text[m.end():k].strip()

    def resolve(name, seen=None):
        seen = seen or set()
        if name in seen or name not in raw:
            return set()
        seen.add(name)
        rhs = raw[name]
        for op, cl in (("{", "}"), ("[", "]")):
            if op in rhs:
                if cl not in rhs:
                    return set()      # malformed rhs — skip, don't crash
                inner = rhs[rhs.index(op) + 1:rhs.rindex(cl)]
                out = set()
                depth = 0
                cur = ""
                for ch in inner.replace("\n", ","):
                    if ch == "(":
                        depth += 1
                    elif ch == ")":
                        depth -= 1
                    if ch == "," and depth == 0:
                        if cur.strip():
                            out.add(cur.strip())
                        cur = ""
                    else:
                        cur += ch
                if cur.strip():
                    out.add(cur.strip())
                # expand op(NAME) indirections through the named list
                final = set()
                for t in out:
                    m2 = re.match(r"op\(\s*(\w+)\s*\)$", t)
                    if m2:
                        final |= resolve(m2.group(1), seen)
                    else:
                        final.add(t)
                return final
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


def forwarders(text):
    """Detect one-line forwarding wrappers in a source file:
    functions whose only executable statement is
    NAME = OTHER(...).  -> [(name, other)]"""
    out = []
    chunks = re.split(r"(?=^\s{0,6}(?:double\s+precision\s+)?"
                      r"function\s+\w+)", text, flags=re.I | re.M)
    for ch in chunks:
        m = FUNC.match(ch)
        if not m:
            continue
        name = m.group(1)
        body = []
        for ln in ch.splitlines()[1:]:
            s = ln.strip().lower()
            if not s or ln[:1] in "cC*!" or s.startswith("!"):
                continue
            if s.startswith(("implicit", "return", "end", "use ")):
                continue
            body.append(ln.strip())
        if len(body) == 1:
            m2 = re.match(rf"{re.escape(name)}\s*=\s*(\w+)\s*\(",
                          body[0], re.I)
            if m2 and m2.group(1).lower() != name.lower():
                out.append((name, m2.group(1)))
    return out


# ---------------- audit ----------------

def audit(root, name_filter=None, show_all=False, out=print):
    decls = []
    fwds = []
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
            for name, other in forwarders(text):
                fwds.append((name, other, os.path.join(d, f)))
    alldecls = decls          # the notation cross-check needs them ALL
    if name_filter:
        decls = [d for d in decls if name_filter.lower() in d[0].lower()]
        fwds = [d for d in fwds if name_filter.lower() in d[0].lower()]

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
        for ln in call_recipe(name, args).splitlines():
            out("  " + ln)

    # cross-check notation.map
    npath = os.path.join(root, "maple", "notation.map")
    if os.path.isfile(npath):
        members, fortmap = parse_notation(open(npath,
                                               errors="replace").read())
        tokens = set()
        for s in ("ant30set", "ant31set", "ant40set"):
            tokens |= members.get(s, set())
        declnames = {n for n, _, _ in alldecls}
        if name_filter:
            tokens = {t for t in tokens
                      if name_filter.lower() in t.lower()}
        missing_fort = sorted(
            t for t in tokens
            if fortmap.get(t, f"Full{t}") not in declnames
            and t not in declnames)
        registered = set(fortmap.values()) | set(fortmap.keys())
        unregistered = sorted(
            n for n in declnames
            if ANTENNA_LIKE.search(n) and n not in registered
            and (not name_filter or name_filter.lower() in n.lower()))
        out(f"\n=== notation.map tokens with no Fortran entry point "
            f"({len(missing_fort)}) ===")
        out("  " + (", ".join(missing_fort) or "none"))
        out(f"\n=== antenna-like entry points not registered in any "
            f"ant*fortset ({len(unregistered)}) ===")
        out("  " + (", ".join(unregistered) or "none"))
        stale = [(n, o, p) for n, o, p in fwds if n not in registered]
        out(f"\n=== UNREGISTERED forwarding wrappers "
            f"({len(stale)}) ===")
        if stale:
            for n, o, p in stale:
                out(f"  {n} -> {o}   [{p}]")
            out("  ^ local renamings no ant*fortset resolves: a .map "
                "cannot call them, their")
            out("    header claims are unverified and rot silently — "
                "measure the TARGET with")
            out("    antenna_datasheet.py / the pole scan instead of "
                "trusting the comment.")
        else:
            out("  none")
    else:
        out("\n(notation.map not found — cross-check skipped)")


# ---------------- selftest ----------------

def selftest():
    """Synthetic sources and notation text — encodes no real antenna's
    order. Checks declaration parsing, slot extraction, the ascending
    test, the permuted-call recipe, the forwarding-wrapper scan and
    the notation cross-check (including a colon inside a braced set)."""
    src = (
        "      function AAA40(i1,i2,i3,i4,ipset)\n"
        "      end\n"
        "      function BBB40(i1,i3,i4,i2,ipset)\n"
        "      end\n"
        "      double precision function CCC31(i3,i2,i1,ipset,renscale)\n"
        "      end\n"
        "      function helper(x,y)\n"
        "      end\n"
        "      function BBB40bb(i1,i2,i3,i4,ipset)\n"
        "      implicit double precision (a-h,o-z)\n"
        "      BBB40bb = BBB40(i1,i2,i4,i3,ipset)\n"
        "      return\n"
        "      end\n")
    decls = parse_decls(src)
    assert [d[0] for d in decls] == \
        ["AAA40", "BBB40", "CCC31", "helper", "BBB40bb"]
    s, e, asc = slot_info(decls[0][1])
    assert s == ["i1", "i2", "i3", "i4"] and e == ["ipset"] and asc
    s, e, asc = slot_info(decls[1][1])
    assert not asc and s == ["i1", "i3", "i4", "i2"]
    s, e, asc = slot_info(decls[2][1])
    assert not asc and e == ["ipset", "renscale"]
    assert slot_info(decls[3][1])[0] == []          # no slots -> skipped

    # permuted-call recipe: declared (i1,i3,i4,i2) -> named
    # (i1,i2,i3,i4)=(A,B,C,D) realised by positional (A,C,D,B)
    w = call_recipe("BBB40", decls[1][1])
    assert "BBB40(A,C,D,B,ipset)" in w
    assert "Do NOT write a slot-reordering wrapper" in w
    w = call_recipe("CCC31", decls[2][1])
    assert "CCC31(C,B,A,ipset,renscale)" in w

    # forwarding-wrapper detection
    fw = forwarders(src)
    assert fw == [("BBB40bb", "BBB40")], fw

    notation = (
        "antXXFFset:={TOKA, TOKB}:\n"
        "antXXIFset:={TOKC}:\n"
        "ant40set:= antXXFFset union antXXIFset:\n"
        "ant30set:={}:\n"
        "ant31set:={}:\n"
        "XXfortset:={\nTOKA=AAA40,\nTOKB=BBB40\n}:\n"
        "ant40fortset:=XXfortset:\n"
        # a colon INSIDE a braced set must not truncate the capture —
        # this input crashed the pre-fix parser
        "commentset:={\"TOKA\"=\"Eqs.~(5.27): see paper\"}:\n")
    members, fortmap = parse_notation(notation)
    assert members["ant40set"] == {"TOKA", "TOKB", "TOKC"}, members
    assert fortmap == {"TOKA": "AAA40", "TOKB": "BBB40"}
    # TOKC has no fortmap entry and no FullTOKC declaration -> missing;
    # CCC31 is antenna-like and unregistered; BBB40bb is an
    # unregistered forwarding wrapper
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
    assert "BBB40(A,C,D,B,ipset)" in txt
    assert "TOKC" in txt.split("no Fortran entry point")[1]
    assert "CCC31" in txt.split("not registered")[1]
    assert "BBB40bb -> BBB40" in txt.split("forwarding wrappers")[1]
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
