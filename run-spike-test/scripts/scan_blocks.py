#!/usr/bin/env python3
"""Unattended sweep over subtraction-term configurations.

Composes a .map from each candidate block set (write-subtraction's
compose_blocks.py), regenerates + rebuilds, runs the spike test, scores the
GENUINE modes, and prints a ranked table. Turns "is this block set
better?" from a manual regenerate-build-run-read cycle into one command.

Why it exists: reasoning narrows the candidate space but does not empty
it, and the remaining choices (which S,b1, which iterated family, which
absorption order, which coefficient) are cheap to test and expensive to
argue about. A sweep also PROVES a negative — if no configuration in the
declared space scores better, the missing piece is a line form that is
not in the master at all, which is a different (and much more useful)
conclusion than "I tried some things".

Usage
  # explicit block sets, one per line, "label<TAB>blk,blk,..." on stdin
  compose_blocks.py enumerate master.map --fixed Sa --axes absorb,sb1 \\
    | scan_blocks.py --master master.map --out <TERM>.map \\
        --regen "bash .../regen_rebuild.sh -n 13 -l RR -s src/process/P -t test/process/P -m checkNtoM" \\
        --run   "./checkNtoM 1 1 80 60 2" --run-cwd test/process/P \\
        --genuine 4,8,13,14 --results scan.txt

  scan_blocks.py --selftest

Options
  --genuine     comma list of GENUINE mode numbers, or
  --modes-spec  a genuine_modes.py spec (write-spike-test) to derive them
  --tol         |median-1| tolerance for "exact" (default 1e-4); a mode
                counts as exact only if it ALSO has zero outliers
  --near        |median-1| tolerance for "near" (default 0.10)
  --keep        directory for the per-combination raw outputs

Two harness lessons are baked in, because both cost a debugging round:
  * the loop is synchronous — never poll with `pgrep -f <this script>`,
    which matches the poller's own command line and never terminates;
  * ranking parses the score fields as integers — never `sort` on a
    line whose first field is "[12]".
"""
import argparse
import os
import re
import subprocess
import sys

MODE = re.compile(r"^\s*mode\s+(\d+):")
XLINE = re.compile(r"x=\s*(\S+).*?med=\s*(\S+).*?outl=\s*(\d+)")


def parse_run(text):
    """-> {mode: (median, outliers)} using the FIRST x line of each mode."""
    out, cur = {}, None
    for ln in text.splitlines():
        m = MODE.match(ln)
        if m:
            cur = int(m.group(1))
            continue
        m = XLINE.search(ln)
        if m and cur is not None and cur not in out:
            try:
                out[cur] = (float(m.group(2)), int(m.group(3)))
            except ValueError:
                out[cur] = (float("nan"), 10 ** 9)
    return out


def score(res, genuine, tol, near_tol):
    """-> (exact, near, bad) lists of mode numbers."""
    exact, near, bad = [], [], []
    for k in sorted(genuine):
        if k not in res:
            bad.append(k)
            continue
        md, outl = res[k]
        if md == md and abs(md - 1.0) < tol and outl == 0:
            exact.append(k)
        elif md == md and abs(md - 1.0) < near_tol:
            near.append(k)
        else:
            bad.append(k)
    return exact, near, bad


def genuine_from_spec(path):
    """Derive the GENUINE mode numbers via write-spike-test's classifier."""
    here = os.path.dirname(os.path.abspath(__file__))
    gm = os.path.normpath(os.path.join(
        here, "..", "..", "write-spike-test", "scripts", "genuine_modes.py"))
    if not os.path.exists(gm):
        raise SystemExit(f"ERROR: genuine_modes.py not found at {gm}")
    out = subprocess.run([sys.executable, gm, path],
                         capture_output=True, text=True).stdout
    modes = []
    for i, ln in enumerate(out.splitlines(), 0):
        if "[GENUINE]" in ln:
            modes.append(i)
    if not modes:
        raise SystemExit("ERROR: genuine_modes.py reported no GENUINE modes; "
                         "pass --genuine explicitly")
    raise SystemExit(
        "ERROR: mode NUMBERS must match your check program's ordering; "
        "pass --genuine explicitly (genuine_modes.py prints families, not "
        "the program's numbering)")


def run(cmd, cwd=None):
    p = subprocess.run(cmd, shell=True, cwd=cwd, capture_output=True,
                       text=True)
    return p.returncode, p.stdout + p.stderr


def selftest():
    """Structural only — encodes no physics and runs no build."""
    sample = "\n".join([
        " mode 4: a,b soft [GENUINE] ",
        "   x= 1.0E-07 n= 10 nan= 0 max|ME|= 1.0E+05 med=    1.000000 "
        "min= 1.00E+00 max= 1.00E+00 outl=   0",
        "   x= 1.0E-08 n= 10 nan= 0 max|ME|= 1.0E+06 med=    9.000000 "
        "min= 1.00E+00 max= 1.00E+00 outl=  10",
        " mode 5: c||d [dead] ",
        "   x= 1.0E-07 n= 10 nan= 0 max|ME|= 1.0E-09 med=    0.500000 "
        "min= 1.00E-01 max= 9.00E-01 outl=   9",
        " mode 6: e||f [GENUINE] ",
        "   x= 1.0E-07 n= 10 nan= 0 max|ME|= 1.0E+05 med=    1.030000 "
        "min= 1.00E+00 max= 1.10E+00 outl=   4",
        " mode 7: g soft [GENUINE] ",
        "   x= 1.0E-07 n= 10 nan= 0 max|ME|= 1.0E+05 med=    1.700000 "
        "min= 1.00E+00 max= 3.00E+00 outl=  10",
    ])
    res = parse_run(sample)
    assert res[4] == (1.0, 0), res
    assert res[6][1] == 4 and res[7][0] == 1.7
    # only the FIRST x line of a mode is used
    assert res[4][0] == 1.0, "later x line overwrote the first"
    ex, ne, bad = score(res, [4, 6, 7], 1e-4, 0.10)
    assert ex == [4] and ne == [6] and bad == [7], (ex, ne, bad)
    # a genuine mode absent from the output counts as bad, never exact
    ex, ne, bad = score(res, [4, 99], 1e-4, 0.10)
    assert ex == [4] and bad == [99], (ex, bad)
    # zero outliers is required even when the median is perfect
    res2 = {1: (1.0, 3)}
    ex, ne, bad = score(res2, [1], 1e-4, 0.10)
    assert ex == [] and ne == [1], (ex, ne)
    # ranking is numeric, not lexicographic ("10" must beat "9")
    rows = [("a", 9, 0), ("b", 10, 0)]
    rows.sort(key=lambda r: (-r[1], -r[2]))
    assert rows[0][0] == "b", "ranking is lexicographic"
    print("scan_blocks selftest OK")


def main():
    if "--selftest" in sys.argv:
        selftest()
        return
    ap = argparse.ArgumentParser()
    ap.add_argument("--master", required=True)
    ap.add_argument("--out", required=True,
                    help="path the composed .map is written to")
    ap.add_argument("--regen", required=True)
    ap.add_argument("--run", required=True)
    ap.add_argument("--run-cwd", default=None)
    ap.add_argument("--regen-cwd", default=None)
    ap.add_argument("--genuine", default="")
    ap.add_argument("--modes-spec", default="")
    ap.add_argument("--tol", type=float, default=1e-4)
    ap.add_argument("--near", type=float, default=0.10)
    ap.add_argument("--results", default="scan_results.txt")
    ap.add_argument("--keep", default=None)
    args = ap.parse_args()

    if args.genuine:
        genuine = [int(s) for s in args.genuine.split(",") if s.strip()]
    elif args.modes_spec:
        genuine = genuine_from_spec(args.modes_spec)
    else:
        raise SystemExit("ERROR: give --genuine (or --modes-spec)")

    here = os.path.dirname(os.path.abspath(__file__))
    mb = os.path.normpath(os.path.join(
        here, "..", "..", "write-subtraction", "scripts", "compose_blocks.py"))
    combos = []
    for ln in sys.stdin:
        ln = ln.rstrip("\n")
        if not ln.strip():
            continue
        label, _, blocks = ln.partition("\t")
        combos.append((label.strip(), (blocks or label).strip()))
    if not combos:
        raise SystemExit("ERROR: no combinations on stdin")

    if args.keep:
        os.makedirs(args.keep, exist_ok=True)
    rows = []
    with open(args.results, "w") as fh:
        for i, (label, blocks) in enumerate(combos, 1):
            tag = f"[{i}/{len(combos)}] {label}"
            rc, _ = run(f"{sys.executable} {mb} compose {args.master} "
                        f"--blocks '{blocks}' -o {args.out}")
            if rc:
                fh.write(f"COMPOSE-FAIL  {label}  ||  {blocks}\n"); fh.flush()
                print(f"{tag}  COMPOSE-FAIL", flush=True); continue
            rc, log = run(args.regen, cwd=args.regen_cwd)
            if rc:
                fh.write(f"REGEN-FAIL    {label}  ||  {blocks}\n"); fh.flush()
                print(f"{tag}  REGEN-FAIL", flush=True); continue
            rc, out = run(args.run, cwd=args.run_cwd)
            if args.keep:
                with open(os.path.join(args.keep, f"run_{i}.txt"), "w") as g:
                    g.write(out)
            ex, ne, bad = score(parse_run(out), genuine, args.tol, args.near)
            rows.append((label, blocks, len(ex), len(ne), len(bad)))
            fh.write(f"EXACT {len(ex):>3}/{len(genuine)}  NEAR {len(ne):>3}  "
                     f"BAD {len(bad):>3}  ||  {label}  ||  {blocks}\n")
            fh.flush()
            print(f"{tag}  exact={len(ex)} near={len(ne)} bad={len(bad)}",
                  flush=True)

    rows.sort(key=lambda r: (-r[2], -r[3]))
    print(f"\n=== ranked ({len(rows)} scored, {len(genuine)} genuine modes) ===")
    print(f"{'EXACT':>5} {'NEAR':>5} {'BAD':>4}   configuration")
    for label, blocks, e, n, b in rows:
        print(f"{e:>5} {n:>5} {b:>4}   {label or blocks}")
    if rows and rows[0][2] < len(genuine):
        print("\nNo configuration reached all genuine modes. Within the "
              "declared space this is a NEGATIVE RESULT: the missing piece "
              "is a line form absent from the master, not a subset of it. "
              "Next step is measurement (probe-me-ir-structure), not more "
              "searching.")


if __name__ == "__main__":
    main()
