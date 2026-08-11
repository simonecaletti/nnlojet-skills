#!/usr/bin/env python3
"""Solve a subtraction term's line coefficients against the matrix element,
and decide whether the term is CLOSABLE by recoefficienting or needs a
structurally new block.

This is the question a block sweep cannot answer. scan_blocks.py searches
the configurations you WROTE; when none of them scores, the useful next
statement is not "I tried some things" but "no assignment of coefficients
to these lines satisfies the failing and the passing modes at once, so the
missing object is not in their span". That is a rank statement about the
design matrix, and this script computes it.

Per accepted phase-space point the spike test evaluates

    ME  =  s * sum_i  c_i * wt_i    (c_i = 1 for the term as written)

with an overall SIGN s fixed by the generator: makefortRR emits
`wt(i)=bino(ix,partons,-wt(i)*facnorm)`, so the dumped wt_i are MINUS the
subtraction contributions and s = -1.  The sign is auto-detected per run
(median of ME / sum wt) and reported, because getting it wrong silently
re-anchors the whole fit: "as written" lands at -1 instead of +1, every
residual against the current term is meaningless, and the verdict inverts.

wt_attribute.py reports each wt_i one mode at a time; here the same dumps
are fed to a least-squares solve for the c_i, over any set of modes at
once.

Prerequisites
  * the term's auto*.f generated with the per-line dump
    (`maple makefortRR -Diprocess=<N> -Dwtdebug=1`; wtdebug_install.py, or
    the recipe wt_attribute.py prints when the branch is absent);
  * the check program emitting its full ME per point as
        WTDBG MEFULL 0 0 <wt1>
    which gen_spike_test.py writes under NNLOJET_WTDEBUG=1;
  * the run made with NNLOJET_WTDEBUG=1 and OMP_NUM_THREADS=1.

Usage
  fit_lines.py --map TERM.map --fn NAME --log run.out \\
      --modes 4,8,14,15,16,20 --hold 67,71,74,75,76,77,78
  fit_lines.py --selftest

  --modes   modes the fit is asked to satisfy (default: every mode found)
  --hold    modes that currently PASS and must not be broken.  With --hold
            the verdict below is available; without it only coefficients
            are reported.
  --anchor  ridge target, 'current' (=1, default) or 'zero'.  Degenerate
            directions are pinned to the anchor instead of running away to
            +-1e4, which is what an unregularised solve does here and what
            makes its output unreadable.
  --lambda  ridge strength (default 1e-2).
  --deepest N   use only the N deepest x blocks per mode (default 2).

Reading the verdict
  The test is the RATIO of the held-mode residual under the fitted
  coefficients to the same residual under the term as written -- not an
  absolute cut.  Gluon-parent collinear modes have a per-point azimuthal
  floor (only the rotation average converges), so a correct term can sit
  at ~1e-1 there; both numbers carry that floor and the ratio divides it out.

  CLOSABLE     one coefficient set satisfies --modes and --hold together.
               The printed rationals are the term you should write.
  NOT CLOSABLE the failing modes can only be fixed by coefficients that
               break the held ones (or cannot be fixed at all).  Stop
               tuning: derive the missing block (write-subtraction's S,c
               construction, or a residue fit for the missing X40).

Coefficients are judged by MEMBERSHIP in the antenna family
{2, 1, 2/3, 1/2, 1/3, 1/4, 1/9} (antennae-naming-convention), never by
decimal proximity: 0.667 is 2/3 and fine, 0.70 is nothing.

Rotated modes.  A collinear mode with a gluon-parent cluster evaluates the
term n_azim times per point.  gen_spike_test.py prints MEFULL immediately
after the UNROTATED evaluation, so the stream per point reads

    [block]  MEFULL  [rot 1] [rot 2] ... [next block] MEFULL ...

and each MEFULL is paired with the block current at that moment; the
rotation blocks that follow are overwritten by the next point's block and
never enter the fit.  Pairing a rotated block with that ME silently
fabricates constraint rows -- the failure mode this parser exists to
avoid, and the reason the selftest emits rotations in the real order.
"""
import argparse
import re
import sys
from fractions import Fraction

LABEL = re.compile(r"\*a(\d+)\b")
DUMP = re.compile(r"WTDBG\s+(\S+)\s+(\d+)\s+(-?\d+)\s+([-+.\dEe]+)")
MELINE = re.compile(r"WTDBG\s+MEFULL\s+0\s+0\s+([-+.\dEe]+)")
XHDR = re.compile(r"x=\s*[\d.eE+-]+\s+n=")
MODEHDR = re.compile(r"^\s*mode\s+(\d+):")
FAMILY = [Fraction(2), Fraction(1), Fraction(2, 3), Fraction(1, 2),
          Fraction(1, 3), Fraction(1, 4), Fraction(1, 9)]


def map_terms(text):
    """-> {label:int -> term text}, same convention as wt_attribute.py."""
    try:
        body = text[text.index("XX:="):]
    except ValueError:
        raise SystemExit("ERROR: no 'XX:=' in map file")
    terms, prev = {}, 0
    for m in LABEL.finditer(body):
        chunk = [ln.strip() for ln in body[prev:m.end()].splitlines()]
        chunk = [ln for ln in chunk
                 if ln and not ln.startswith("#") and ln != "XX:="]
        terms[int(m.group(1))] = " ".join(chunk)
        prev = m.end()
    if not terms:
        raise SystemExit("ERROR: no *aN labels found in map file")
    return terms


def parse_log(lines, fn=None):
    """-> {mode: [ [ (me, {i: wt}) ... ] per x block ]}."""
    out, mode, blocks = {}, None, None
    cur, prev = {}, None
    for ln in lines:
        m = MODEHDR.match(ln)
        if m:
            mode = int(m.group(1))
            blocks = out.setdefault(mode, [[]])
            cur, prev = {}, None
            continue
        m = DUMP.match(ln.strip())
        if m and m.group(1) != "MEFULL":
            if fn and m.group(1) != fn:
                continue
            i = int(m.group(2))
            if i in cur:                 # a rotation block began
                prev, cur = dict(cur), {}
            cur[i] = float(m.group(4))
            continue
        m = MELINE.search(ln)
        if m:
            blk = dict(cur) if cur else prev
            if blk and blocks is not None:
                blocks[-1].append((float(m.group(1)), blk))
            cur, prev = {}, None
            continue
        if XHDR.search(ln) and blocks is not None:
            blocks.append([])
            cur, prev = {}, None
    return {k: [b for b in v if b] for k, v in out.items() if any(v)}


def detect_sign(points):
    """s such that ME ~ s * sum(wt); -1 for makefortRR output."""
    rs = []
    for me, w in points:
        t = sum(w.values())
        if t:
            rs.append(me / t)
    if not rs:
        return 1.0, float("nan")
    rs.sort()
    med = rs[len(rs) // 2]
    return (-1.0 if med < 0 else 1.0), med


def rows(points, cols, sign=1.0):
    """Design rows scaled per point so no single deep point dominates."""
    A, b = [], []
    for me, w in points:
        r = [sign * w.get(c, 0.0) for c in cols]
        if not any(r):
            continue
        s = max(max(abs(x) for x in r), abs(me)) or 1.0
        A.append([x / s for x in r])
        b.append(me / s)
    return A, b


def solve(A, b, lam, anchor):
    """Ridge solve without numpy: (A'A + lam I) c = A'b + lam*anchor."""
    n = len(A[0])
    ata = [[sum(A[k][i] * A[k][j] for k in range(len(A))) for j in range(n)]
           for i in range(n)]
    atb = [sum(A[k][i] * b[k] for k in range(len(A))) for i in range(n)]
    for i in range(n):
        ata[i][i] += lam
        atb[i] += lam * anchor
    # Gaussian elimination with partial pivoting
    M = [row[:] + [atb[i]] for i, row in enumerate(ata)]
    for i in range(n):
        p = max(range(i, n), key=lambda r: abs(M[r][i]))
        if abs(M[p][i]) < 1e-300:
            raise SystemExit("ERROR: singular normal equations even with a "
                             "ridge; raise --lambda")
        M[i], M[p] = M[p], M[i]
        for r in range(n):
            if r == i:
                continue
            f = M[r][i] / M[i][i]
            for c in range(i, n + 1):
                M[r][c] -= f * M[i][c]
    return [M[i][n] / M[i][i] for i in range(n)]


def residual(A, b, c):
    num = sum((sum(A[k][i] * c[i] for i in range(len(c))) - b[k]) ** 2
              for k in range(len(A)))
    den = sum(x * x for x in b) or 1.0
    return (num / den) ** 0.5


def snap(v):
    """Nearest family member if within tolerance, else None."""
    for f in FAMILY:
        for s in (1, -1):
            if abs(v - s * float(f)) < 0.03:
                return s * f
    if abs(v) < 0.03:
        return Fraction(0)
    return None


def report(cols, c, terms, quiet):
    print(f"  {'aN':>4}  {'fitted':>9}  {'family':>7}   line")
    for k, v in zip(cols, c):
        f = snap(v)
        tag = str(f) if f is not None else "--"
        flag = "" if f is not None else "   <-- outside the family"
        txt = terms.get(k, "")
        if quiet and abs(v - 1.0) < 0.03:
            continue
        print(f"  a{k:<3}  {v:+9.4f}  {tag:>7}{flag}")
        if not quiet:
            print(f"        {txt[:96]}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--map", required=True)
    ap.add_argument("--fn")
    ap.add_argument("--log", required=True)
    ap.add_argument("--modes")
    ap.add_argument("--hold")
    ap.add_argument("--anchor", choices=["current", "zero"], default="current")
    ap.add_argument("--lambda", dest="lam", type=float, default=1e-2)
    ap.add_argument("--deepest", type=int, default=2)
    ap.add_argument("--quiet", action="store_true",
                    help="print only lines whose fitted value differs from 1")
    a = ap.parse_args()

    terms = map_terms(open(a.map).read())
    cols = sorted(terms)
    data = parse_log(open(a.log, errors="replace"), a.fn)
    if not data:
        raise SystemExit(
            "ERROR: no WTDBG points parsed.\n"
            "  * generated with -Dwtdebug=1?   (wtdebug_install.py --check)\n"
            "  * check program emitting 'WTDBG MEFULL'?  (gen_spike_test.py)\n"
            "  * run with NNLOJET_WTDEBUG=1 and OMP_NUM_THREADS=1?")

    def pts(modes):
        out = []
        for m in modes:
            for blk in data.get(m, [])[-a.deepest:]:
                out += blk
        return out

    want = ([int(x) for x in a.modes.split(",")] if a.modes
            else sorted(data))
    hold = [int(x) for x in a.hold.split(",")] if a.hold else []
    anchor = 1.0 if a.anchor == "current" else 0.0

    sign, med = detect_sign(pts(sorted(data)))
    print(f"sign convention: ME = {sign:+.0f} * sum(wt)   "
          f"(median ME/sum wt = {med:.4f})")
    if abs(abs(med) - 1.0) > 0.25:
        print("  WARNING: median far from +-1 -- the term is not close to "
              "the ME on\n  these modes, or dumps and ME are misaligned; "
              "treat the fit with care.")
    A, b = rows(pts(want), cols, sign)
    if not A:
        raise SystemExit(f"ERROR: no points for modes {want}")
    print(f"modes fitted {want}   ({len(A)} points, {len(cols)} lines, "
          f"ridge {a.lam} at {anchor:g})")
    c = solve(A, b, a.lam, anchor)
    print(f"  residual on fitted modes: {residual(A, b, c):.2e}")
    report(cols, c, terms, a.quiet)

    if not hold:
        print("\nno --hold given: verdict unavailable (pass the modes that "
              "currently PASS to get one)")
        return
    Ah, bh = rows(pts(hold), cols, sign)
    if not Ah:
        raise SystemExit(f"ERROR: no points for held modes {hold}")
    rh = residual(Ah, bh, c)
    r1 = residual(Ah, bh, [1.0] * len(cols))
    ratio = rh / max(r1, 1e-12)
    print(f"\n  held modes {hold}")
    print(f"    residual with the fitted coefficients : {rh:.2e}")
    print(f"    residual with the term as written     : {r1:.2e}")
    print(f"    ratio fitted/as-written               : {ratio:.1f}")
    # The verdict compares the two, never an absolute threshold.  On a
    # gluon-parent collinear mode a single UNROTATED point deviates at
    # O(1) in azimuth -- only the rotation average converges -- so even a
    # perfect term has a per-point residual floor of ~1e-1 there.  An
    # absolute cut would call a correct term NOT CLOSABLE; the ratio is
    # blind to the floor because both numbers carry it.
    if ratio < 3.0:
        print("\nVERDICT: CLOSABLE — one coefficient set satisfies both "
              "sets.\n  Write the rationals above; anything outside the "
              "family means the\n  basis is still incomplete.")
    else:
        print("\nVERDICT: NOT CLOSABLE by recoefficienting.\n"
              "  The failing modes are only fixed by coefficients that "
              "break the held\n  ones, so their deficit is not in the span "
              "of these lines.  Derive the\n  missing block (write-"
              "subtraction: S,c construction / residue fit for a\n  missing "
              "X40) instead of scanning coefficients.")


def selftest():
    """Structural only: synthetic dumps with known coefficients."""
    terms = map_terms("XX:=\n+A(1)*a1\n+B(2)*a2\n+C(3)*a3\n:\n")
    assert sorted(terms) == [1, 2, 3], terms

    # ME = 1*w1 + 0.5*w2 + 1*w3.  Real stream order: block, MEFULL, rotations.
    truth = [1.0, 0.5, 1.0]
    lines = ["  mode   9: x [GENUINE]"]
    for n in range(1, 9):
        w = [n * 1.0, 0.1 * n * n + 1.0, 1.0 / n + 0.5]
        me = sum(t * x for t, x in zip(truth, w))
        for i, x in enumerate(w, 1):
            lines.append(f"WTDBG T {i} 1  {x:.15E}")
        lines.append(f"WTDBG MEFULL 0 0  {me:.15E}")
        for _ in range(2):                     # rotations: must never be used
            for i in range(1, 4):
                lines.append(f"WTDBG T {i} 1  {9.9e9:.15E}")
    d = parse_log(lines, "T")
    assert list(d) == [9], d
    assert len(d[9][0]) == 8, len(d[9][0])
    for me, blk in d[9][0]:
        assert max(abs(v) for v in blk.values()) < 1e6, \
            "rotation block leaked into the fit"

    A, b = rows(d[9][0], [1, 2, 3])
    c = solve(A, b, 1e-9, 1.0)
    assert all(abs(x - t) < 1e-4 for x, t in zip(c, truth)), c
    assert residual(A, b, c) < 1e-6
    assert snap(0.6667) == Fraction(2, 3) and snap(0.70) is None
    assert snap(-0.5) == Fraction(-1, 2)
    print("fit_lines selftest OK")


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        selftest()
    else:
        main()
