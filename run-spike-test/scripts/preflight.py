#!/usr/bin/env python3
"""Preflight for a spike test: validate the TEST PROGRAM and the build
chain BEFORE believing any ratio it prints.

Two failure classes are silent, look exactly like a broken .map, and
between them account for the most expensive false diagnoses:

 1. STALE BUILD. The chain is .map -> (maple) auto*.f -> .o -> binary,
    plus check*.f -> .o -> binary. A break anywhere leaves an executable
    that tests something other than what you are reading. Measured: a
    check*.f edited three minutes AFTER its binary was linked made a
    CORRECT term read ~45% outliers with median 0.9988 on two
    triple-collinear modes; rebuilt, the same .map gave 0/757 and
    1.000000.

 2. INCOMPLETE AZIMUTHAL AVERAGING. An antenna is spin-averaged, so a
    GLUON-parent collinear cluster only agrees with the matrix element
    after averaging over the cluster's azimuth. The check program must
    rotate EVERY leg of such a cluster (`nrot=N`, `mrot(1..N)` = the
    cluster). Rotating a subset averages nothing at all, and the symptom
    -- ~45% outliers on precisely the gluon-parent collinear modes -- is
    the single most physically plausible-looking artefact the harness
    produces, because those are the modes whose failure you would
    believe.

    QUARK-parent clusters need no rotation; flagging them would be a
    false positive. So the audit computes each cluster's parent by net
    flavour from the program's own channel banner and only requires
    rotation where the parent is a gluon.

Usage:
  preflight.py test/process/epem check5to3 [--map-dir maple/process/epem]
               [--src-dir src/process/epem] [--root .]

Exit 1 on a hard finding (stale chain, or a gluon-parent cluster that is
unrotated / partly rotated).
"""

import argparse
import glob
import json
import os
import re
import sys

CASE_RE = re.compile(r"^\s*case\s*\(\s*(\d+)\s*\)", re.I)
STITLE_RE = re.compile(r"stitle\s*=\s*'([^']*)'", re.I)
NROT_RE = re.compile(r"^\s*nrot([12])\s*=\s*(\d+)", re.I)
MROT_RE = re.compile(r"^\s*mrot([12])\s*\(\s*(\d+)\s*\)\s*=\s*(\d+)", re.I)


def parse_banner(text):
    """index -> parton name, from the program's own usage banner.

    Looks for `NAME(n)` tokens on a line describing the channel, e.g.
    '  1: ep(1) em(2) -> qb1(3) g(4) Q(5) Qb(6) q1(7)  [C1g0Z]'.
    """
    for line in text.split("\n"):
        if "->" not in line:
            continue
        rhs = line.split("->", 1)[1]
        hits = re.findall(r"([A-Za-z][A-Za-z0-9]*)\s*\(\s*(\d+)\s*\)", rhs)
        if len(hits) >= 3:
            return {int(n): nm for nm, n in hits}
    return {}


def parse_spec(root, proc):
    """index -> parton name from maple/process/<proc>/*.spec.json.

    Authoritative and format-independent: the generated banner does not
    always carry indices (`-> qb g Q Qb q`), so parsing prose is a
    fallback, not the primary source.
    """
    for path in sorted(glob.glob(os.path.join(
            root, "maple", "process", proc, "*.spec.json"))):
        try:
            d = json.load(open(path))
        except (OSError, ValueError):
            continue
        tab = d.get("fs_partons") or d.get("partons") or {}
        if tab:
            try:
                return {int(k): v for k, v in tab.items()}
            except (TypeError, ValueError):
                continue
    return {}


def flavour(name):
    """parton name -> ('g',None) | ('q',tag) | ('qb',tag).

    Convention in the banners: q/Q/q1/q2 = quark, anything ending in 'b'
    before the tag (qb, Qb, qb1) = antiquark, g = gluon. The TAG groups
    same-flavour pairs: q1/qb1 share a tag, Q/Qb share another.
    """
    n = name.strip()
    if n.lower() == "g":
        return ("g", None)
    m = re.match(r"^([A-Za-z]+?)(b?)(\d*)$", n)
    if not m:
        return (None, None)
    stem, bar, tag = m.group(1), m.group(2), m.group(3)
    # 'Qb' / 'qb1' -> antiquark; 'Q' / 'q1' -> quark. Tag = digits if
    # present, else the stem's case ('q' vs 'Q' distinguishes the pairs
    # in banners that use Q/Qb for the secondary pair).
    kind = "qb" if bar else "q"
    key = tag if tag else stem.lower() + ("U" if stem[0].isupper() else "l")
    return (kind, key)


def parent(idxs, banner):
    """Net flavour of a collinear cluster: 'g', 'q', or None (illegal)."""
    net = {}
    for i in idxs:
        k, tag = flavour(banner.get(i, ""))
        if k is None:
            return None
        if k == "g":
            continue
        net[tag] = net.get(tag, 0) + (1 if k == "q" else -1)
    left = {t: v for t, v in net.items() if v}
    if not left:
        return "g"
    if len(left) == 1 and abs(list(left.values())[0]) == 1:
        return "q"
    return None


def clusters_of(title):
    """Mode title -> list of collinear clusters (lists of indices).

    '3||4||7' -> [[3,4,7]];  '3||4 + 5||6' -> [[3,4],[5,6]];
    '4 soft + 3||7' -> [[3,7]];  '5,6 soft' / '4 soft' -> [].
    """
    body = re.sub(r"\[(GENUINE|dead)\]", "", title).strip()
    out = []
    for part in body.split("+"):
        part = part.strip()
        if "||" not in part:
            continue                       # a soft factor, not a cluster
        part = re.sub(r"\bsoft\b", "", part).strip()
        legs = [p.strip() for p in part.split("||")]
        if all(p.isdigit() for p in legs) and len(legs) >= 2:
            out.append([int(p) for p in legs])
    return out


def parse_modes(text):
    """case(N) -> {title, rot: {1: set, 2: set}, nrot: {1: n, 2: n}}."""
    modes, cur = {}, None
    seen_case = False
    for line in text.split("\n"):
        if line[:1] in ("c", "C", "*", "!"):
            continue
        m = CASE_RE.match(line)
        if m:
            seen_case = True
            cur = modes.setdefault(int(m.group(1)),
                                   {"title": "", "rot": {1: set(), 2: set()},
                                    "nrot": {1: 0, 2: 0}})
            continue
        if not seen_case or cur is None:
            continue
        m = STITLE_RE.search(line)
        if m:
            cur["title"] = m.group(1)
        m = NROT_RE.match(line)
        if m:
            cur["nrot"][int(m.group(1))] = int(m.group(2))
        m = MROT_RE.match(line)
        if m:
            cur["rot"][int(m.group(1))].add(int(m.group(3)))
    return modes


def audit_rotation(path, banner=None):
    text = open(path, errors="replace").read()
    banner = banner or parse_banner(text)
    if not banner:
        return None                        # UNCHECKED, not "clean"
    modes = parse_modes(text)
    if not any(m["title"] for m in modes.values()):
        return None                        # no stitle metadata to audit
    print("  partons: " + " ".join(f"{i}={n}"
                                   for i, n in sorted(banner.items())))
    findings = []
    for n in sorted(modes):
        md = modes[n]
        if "dead" in md["title"]:
            continue
        rotsets = [s for s in md["rot"].values() if s]
        declared = {1: md["nrot"][1], 2: md["nrot"][2]}
        for cl in clusters_of(md["title"]):
            if parent(cl, banner) != "g":
                continue                    # quark parent: no averaging
            match = [s for s in rotsets if s == set(cl)]
            if match:
                # legs listed; check the COUNT agrees with the set
                k = [i for i, s in md["rot"].items() if s == set(cl)][0]
                if declared[k] != len(cl):
                    findings.append(
                        (n, md["title"], cl,
                         f"nrot{k}={declared[k]} but the cluster has "
                         f"{len(cl)} legs — partial rotation"))
                continue
            partial = [s for s in rotsets if s & set(cl)]
            if partial:
                findings.append(
                    (n, md["title"], cl,
                     f"rotated legs {sorted(partial[0])} != cluster "
                     f"{cl} — partial rotation, averages nothing"))
            else:
                findings.append(
                    (n, md["title"], cl,
                     "gluon-parent cluster is NOT rotated at all"))
    return findings


def newer(a, b):
    return os.path.exists(a) and os.path.exists(b) \
        and os.path.getmtime(a) > os.path.getmtime(b)


def audit_freshness(root, testdir, target, srcdir, mapdir):
    binary = os.path.join(root, testdir, target)
    src = os.path.join(root, testdir, target + ".f")
    bad = []
    if not os.path.exists(binary):
        print(f"  binary {testdir}/{target} does not exist — build it")
        return [("build", "binary missing")]
    if newer(src, binary):
        bad.append(("check source", f"{target}.f is NEWER than the binary "
                                    f"— the run tests the OLD program"))
    if srcdir and os.path.isdir(os.path.join(root, srcdir)):
        for f in sorted(os.listdir(os.path.join(root, srcdir))):
            if not f.startswith("auto") or not f.endswith(".f"):
                continue
            p = os.path.join(root, srcdir, f)
            if newer(p, binary):
                bad.append(("generated term",
                            f"{srcdir}/{f} is NEWER than the binary"))
    if mapdir and os.path.isdir(os.path.join(root, mapdir)):
        for f in sorted(os.listdir(os.path.join(root, mapdir))):
            if not f.endswith(".map") or f.startswith("auto"):
                continue
            mp = os.path.join(root, mapdir, f)
            gen = os.path.join(root, srcdir or "",
                               "auto" + f[:-4] + ".f")
            if os.path.exists(gen) and newer(mp, gen):
                bad.append(("regeneration pending",
                            f"{mapdir}/{f} is NEWER than its auto*.f — "
                            f"run maple makefort<RR|RV> first"))
    return bad


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("testdir")
    ap.add_argument("target")
    ap.add_argument("--root", default=".")
    ap.add_argument("--src-dir", default="")
    ap.add_argument("--map-dir", default="")
    a = ap.parse_args()

    proc = os.path.basename(a.testdir.rstrip("/"))
    srcdir = a.src_dir or f"src/process/{proc}"
    mapdir = a.map_dir or f"maple/process/{proc}"
    src = os.path.join(a.root, a.testdir, a.target + ".f")

    print(f"=== build freshness: {a.testdir}/{a.target} ===")
    stale = audit_freshness(a.root, a.testdir, a.target, srcdir, mapdir)
    for kind, msg in stale:
        print(f"  STALE [{kind}] {msg}")
    if not stale:
        print("  ok — binary is newer than its check source, generated "
              "terms and .map files")

    print(f"\n=== azimuthal-rotation audit: {a.target}.f ===")
    if not os.path.exists(src):
        print(f"  ! {src} not found")
        rot = []
    else:
        rot = audit_rotation(src, parse_spec(a.root, proc))
        if rot is None:
            print("  UNCHECKED — no channel banner / no stitle metadata "
                  "in this program, so the rotation rule cannot be "
                  "verified mechanically. This is NOT a pass: check by "
                  "hand that every gluon-parent collinear cluster is "
                  "fully rotated, or regenerate the program with "
                  "gen_spike_test.py (write-spike-test).")
            rot = []
        else:
            for n, title, cl, msg in rot:
                print(f"  UNAVERAGED mode {n:>3} '{title}': cluster {cl} "
                      f"-> gluon parent; {msg}")
            if not rot:
                print("  ok — every gluon-parent collinear cluster of "
                      "every genuine mode is fully rotated")

    if stale or rot:
        print("\nFIX THESE BEFORE READING ANY RATIO. Both classes produce "
              "~45% outliers with median ~1 on exactly the modes whose "
              "failure looks most physical, and are indistinguishable "
              "from a broken .map.")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
