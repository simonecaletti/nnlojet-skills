#!/usr/bin/env python3
"""Per-line weight attribution for a subtraction term under spike test.

Answers "which .map lines are responsible for this mode" in ONE run,
instead of one regenerate+rebuild cycle per hypothesis.

Prerequisite: the term's auto*.f was generated with the opt-in dump
(`maple makefortRR -Diprocess=<N> -Dwtdebug=1`, same for RV — default
off, byte-identical output without the flag) and the check program runs
with environment variable NNLOJET_WTDEBUG=1 and OMP_NUM_THREADS=1.
Each evaluation then prints one line per .map term:

    WTDBG <function> <i> <jpass> <wt>     (RR; RV omits jpass)

Usage:
  wt_attribute.py --map TERM.map [--fn NAME] --x 1e-7 1e-8 1e-9 \\
      (--log run.out | --cmd "./check4to2 1 12 12 200 1")
  wt_attribute.py --selftest

  --map     the .map source; aN labels are gap-free by construction, so
            dump index i maps 1:1 onto the i-th *aN term in file order.
  --fn      filter dumps to one generated function (needed when several
            subtraction functions run in the same execution).
  --x       the x values of the scan, shallow to deep, in the order the
            program visits them.
  --marker  'after' (default): the program prints its per-x summary
            line AFTER that x's points (gen_spike_test.py layout);
            'before': header printed before the points.
  --cmd     run this command (NNLOJET_WTDEBUG=1, OMP_NUM_THREADS=1 are
            set for you); --log parses an existing capture instead.

Output, per .map line: pass rate, median |wt| per x, fitted exponent
alpha (|wt| ~ x^-alpha), magnitude relative to the largest line at the
deepest x, and the source term text. Lines whose exponent reaches the
mode's family maximum are the ones active in this limit; a line with an
anomalous exponent or magnitude is the suspect. Run this BEFORE block
bisection (see run-spike-test SKILL.md).
"""
import argparse
import math
import re
import statistics
import subprocess
import sys

LABEL = re.compile(r"\*a(\d+)\b")
XMARK = re.compile(r"(?<![A-Za-z])x\s*=\s*([0-9][0-9.]*[eEdD][+-]?\d+)")


# ---------------- .map side ----------------

def map_terms(text):
    """-> {label_int: term_text} parsed from the XX:= body."""
    try:
        start = text.index("XX:=")
    except ValueError:
        raise SystemExit("ERROR: no 'XX:=' in map file")
    body = text[start:]
    terms = {}
    prev_end = 0
    for m in LABEL.finditer(body):
        raw = body[prev_end:m.end()]
        # strip the XX:= opener and comment-only lines from the chunk
        lines = [ln.strip() for ln in raw.splitlines()]
        lines = [ln for ln in lines
                 if ln and not ln.startswith("#") and ln != "XX:="]
        terms[int(m.group(1))] = " ".join(lines)
        prev_end = m.end()
    if not terms:
        raise SystemExit("ERROR: no *aN labels found in map file")
    return terms


# ---------------- dump side ----------------

def parse_stream(lines, fn=None):
    """-> (records, marker_positions). records = list of (pos, {i: (jpass, wt)});
    marker_positions = list of (pos, xvalue)."""
    records, markers = [], []
    cur = {}
    cur_start = 0
    pos = 0
    for ln in lines:
        pos += 1
        if ln.startswith("WTDBG"):
            tok = ln.split()
            if len(tok) < 3:
                continue
            name = tok[1]
            if fn and name != fn:
                continue
            i = int(tok[2])
            if len(tok) >= 5:
                jp, wt = int(tok[3]), float(tok[4])
            else:
                jp, wt = 1, float(tok[3])
            if i in cur:            # new evaluation record begins
                records.append((cur_start, cur))
                cur = {}
            if not cur:             # position of the record's FIRST line
                cur_start = pos
            cur[i] = (jp, wt)
        else:
            m = XMARK.search(ln)
            if m:
                markers.append((pos, float(m.group(1).replace("d", "e")
                                           .replace("D", "e"))))
    if cur:
        records.append((cur_start, cur))
    return records, markers


def assign_x(records, markers, xs, marker="after"):
    """-> {x: [record dicts]} using marker positions; falls back to a
    single unassigned group when no markers are found."""
    if not markers:
        return {None: [r for _, r in records]}
    groups = {x: [] for x in xs}
    for rpos, rec in records:
        if marker == "after":
            # record belongs to the first marker AT OR AFTER it
            cand = [(mp, mx) for mp, mx in markers if mp >= rpos]
            mx = cand[0][1] if cand else markers[-1][1]
        else:
            cand = [(mp, mx) for mp, mx in markers if mp <= rpos]
            mx = cand[-1][1] if cand else markers[0][1]
        # snap to the nearest requested x (marker prints are rounded)
        best = min(xs, key=lambda x: abs(math.log10(x) - math.log10(mx)))
        groups[best].append(rec)
    return groups


# ---------------- analysis ----------------

def attribute(terms, groups, xs):
    """-> rows: (label, passrate, {x: med|wt|}, alpha, relmag)."""
    labels = sorted(terms)
    rows = []
    usable_xs = [x for x in xs if groups.get(x)]
    deep = usable_xs[-1] if usable_xs else None
    # largest median |wt| at the deepest x, for relative magnitude
    ref = 0.0
    meds = {}
    for lab in labels:
        meds[lab] = {}
        npass = ntot = 0
        for x in usable_xs:
            vals = []
            for rec in groups[x]:
                if lab in rec:
                    jp, wt = rec[lab]
                    ntot += 1
                    npass += jp
                    vals.append(abs(wt))
            meds[lab][x] = statistics.median(vals) if vals else 0.0
        if deep is not None and meds[lab].get(deep, 0.0) > ref:
            ref = meds[lab][deep]
        passrate = npass / ntot if ntot else 0.0
        rows.append([lab, passrate])
    for row in rows:
        lab = row[0]
        # fitted exponent across consecutive usable x
        alphas = []
        for x1, x2 in zip(usable_xs, usable_xs[1:]):
            m1, m2 = meds[lab][x1], meds[lab][x2]
            if m1 > 0 and m2 > 0:
                alphas.append(-(math.log10(m2) - math.log10(m1))
                              / (math.log10(x2) - math.log10(x1)))
        alpha = statistics.mean(alphas) if alphas else float("nan")
        relmag = (meds[lab][deep] / ref) if (deep and ref > 0) else 0.0
        row.extend([meds[lab], alpha, relmag])
    return rows


def report(rows, terms, xs):
    usable = [x for x in xs if any(r[2].get(x, 0) > 0 for r in rows)]
    print(f"{'aN':>4s} {'pass':>5s} {'alpha':>7s} {'rel|wt|':>9s}  "
          + "  ".join(f"med@{x:.0e}" for x in usable))
    for lab, passrate, meds, alpha, relmag in rows:
        medstr = "  ".join(f"{meds.get(x, 0.0):9.2e}" for x in usable)
        astr = f"{alpha:7.2f}" if alpha == alpha else "    n/a"
        print(f" a{lab:<3d} {passrate:5.2f} {astr} {relmag:9.2e}  {medstr}")
        src = terms[lab]
        print(f"      {src[:100]}{'...' if len(src) > 100 else ''}")


def selftest():
    """Synthetic dump + synthetic .map — no physics encoded. Checks the
    label mapping, x association (both marker layouts) and the fitted
    exponents on constructed power laws."""
    maptext = "\n".join([
        "FN:=DUMMY(1,2,3):",
        "XX:=",
        "# block: one",
        "+F1(3,4,5)*R1([3,4],[4,5])*JET(1)*a1",
        "+F2(3,4,5)*R2([3,4],[4,5])*JET(1)*a2",
        "-F3(3,4,5)*R3([3,4],[4,5])*JET(1)*a3",
        ":",
    ])
    terms = map_terms(maptext)
    assert sorted(terms) == [1, 2, 3] and "F2(" in terms[2]
    xs = [1e-7, 1e-8, 1e-9]

    def dump(x, n=4):
        # wt1 ~ x^-2, wt2 ~ x^-1, wt3 ~ const
        out = []
        for _ in range(n):
            out.append(f"WTDBG DUMFN     1  1 {x ** -2:23.15e}")
            out.append(f"WTDBG DUMFN     2  1 {x ** -1:23.15e}")
            out.append(f"WTDBG DUMFN     3  0 {7.0:23.15e}")
        return out

    for marker in ("after", "before"):
        lines = []
        for x in xs:
            if marker == "before":
                lines.append(f"  x= {x:9.1e} header")
            lines += dump(x)
            if marker == "after":
                lines.append(f"  x= {x:9.1e}  n= 4 summary")
        records, markers = parse_stream(lines, fn="DUMFN")
        assert len(records) == 12, len(records)
        assert len(markers) == 3
        groups = assign_x(records, markers, xs, marker)
        assert all(len(groups[x]) == 4 for x in xs), \
            {x: len(g) for x, g in groups.items()}
        rows = attribute(terms, groups, xs)
        got = {lab: alpha for lab, _, _, alpha, _ in rows}
        assert abs(got[1] - 2.0) < 1e-6, got
        assert abs(got[2] - 1.0) < 1e-6, got
        assert abs(got[3] - 0.0) < 1e-6, got
        pr = {lab: p for lab, p, *_ in rows}
        assert pr[3] == 0.0 and pr[1] == 1.0
    # fn filter drops foreign functions
    records, _ = parse_stream(["WTDBG OTHER 1 1 1.0e0"], fn="DUMFN")
    assert not records
    print("wt_attribute selftest OK")


def repo_root(start=None):
    """Nearest ancestor holding NNLOJET.mk (the maple/ generators live
    next to it). None if we are not inside a tree."""
    import os
    d = os.path.abspath(start or os.getcwd())
    while True:
        if os.path.exists(os.path.join(d, "NNLOJET.mk")):
            return d
        parent = os.path.dirname(d)
        if parent == d:
            return None
        d = parent


def check_generator_support():
    """Fail BEFORE a regenerate+rebuild cycle if this tree's maple
    generators have no -Dwtdebug hook. The dump is opt-in at generation
    time, so without the hook there is nothing to attribute and the run
    is wasted — which is exactly the cost this script exists to avoid.
    """
    import os
    root = repo_root()
    if root is None:
        return                      # not in a tree: let the run decide
    gens = [os.path.join(root, "maple", g)
            for g in ("makefortRR", "makefortRV", "makeformVV")]
    present = [g for g in gens if os.path.exists(g)]
    if not present:
        return
    for g in present:
        try:
            if "wtdebug" in open(g, errors="replace").read():
                return              # at least one generator supports it
        except OSError:
            return
    raise SystemExit(
        "ERROR: none of this tree's maple generators supports -Dwtdebug\n"
        "       (checked: " + ", ".join(os.path.relpath(g, root)
                                        for g in present) + ")\n"
        "  The per-line dump is emitted by the GENERATOR, so no run of the\n"
        "  check program can produce WTDBG lines here — regenerating and\n"
        "  rebuilding first would be wasted.\n"
        "  Use instead: block bisection with compose_blocks.py compose\n"
        "  (write-subtraction) + scan_blocks.py, which needs no dump.\n"
        "  To enable attribution in this tree, add the wtdebug branch to\n"
        "  maple/makefortRR (emit `if(wtdebug=1)` write statements of the\n"
        "  form  WTDBG <fn> <i> <jpass> <wt>  next to each wt(i) fill).")


def main():
    if "--selftest" in sys.argv:
        selftest()
        return
    ap = argparse.ArgumentParser()
    ap.add_argument("--map", required=True)
    ap.add_argument("--fn")
    ap.add_argument("--x", nargs="+", type=float, required=True)
    ap.add_argument("--marker", choices=["after", "before"],
                    default="after")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--log")
    g.add_argument("--cmd")
    args = ap.parse_args()

    terms = map_terms(open(args.map).read())
    if args.log:
        lines = open(args.log).read().splitlines()
    else:
        import os
        check_generator_support()
        env = dict(os.environ, NNLOJET_WTDEBUG="1", OMP_NUM_THREADS="1")
        proc = subprocess.run(args.cmd, shell=True, env=env,
                              capture_output=True, text=True)
        lines = proc.stdout.splitlines()
        if not any(ln.startswith("WTDBG") for ln in lines):
            raise SystemExit(
                "ERROR: no WTDBG lines in output. Was the auto*.f "
                "generated with -Dwtdebug=1 (autogen-subtraction)?")
    records, markers = parse_stream(lines, fn=args.fn)
    groups = assign_x(records, markers, args.x, args.marker)
    if None in groups:
        print("WARNING: no x markers found in output — aggregating all "
              "evaluations, exponents unavailable")
        groups = {args.x[-1]: groups[None]}
    rows = attribute(terms, groups, args.x)
    report(rows, terms, args.x)


if __name__ == "__main__":
    main()
