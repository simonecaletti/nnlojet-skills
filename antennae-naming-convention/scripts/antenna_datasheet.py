#!/usr/bin/env python3
"""Antenna datasheet: a MEASURED, re-verifiable operative description of
each antenna function — the fields a subtraction-term author needs and
that source headers do not reliably state.

Per antenna the datasheet records:
  - positional slot species (q/g per argument slot)  [declared, letter table]
  - the pole graph, MEASURED: which 2-parton invariants are singular
    and with what power; which legs can go soft and WHICH DIPOLE the
    soft eikonal sits on (fitted, with coefficient); which triple
    invariants / double-soft pairs are singular
  - poles it does NOT have (implicit: absent from the lists)
  - for 4-slot antennae: whether the generic cluster rule
    X40(A,B,C,D) -> [A,B,C],[D,C,B] supports every pole
    (`requires_split`): a pole on the ENDPOINT pair (slots 1,4) cannot
    be carried by either cluster, so the Full composite must be split
    into halves with different mappings
  - for Full composites: the halves and the MEASURED calling identity,
    e.g. Full(A,B,C,D) = a(A,B,C,D) + b(A,D,C,B)  — numerically
    verified, so a stale header comment cannot mislead
  - the J21/J22 integrated-dipole coefficient lines mentioning the
    antenna (parsed from maple/form/common/J2[12].map)
  - registration: maple token, Fortran entry point, slot declaration

Everything measured carries provenance (date, test dir, x scan, npt)
and `verify` re-measures and diffs — the datasheet cannot silently rot
the way a hand-written header (e.g. an unregistered wrapper's comment)
can.

Usage:
  antenna_datasheet.py measure NAME [NAME...] --testdir test/process/epem
                       [--root .] [--npt 40] [--xs 1e-5,1e-6,1e-7]
                       [--data <json>] [--keep]
  antenna_datasheet.py static  [--root .] [--data <json>]
  antenna_datasheet.py show    NAME [--data <json>]
  antenna_datasheet.py verify  NAME --testdir ... [--data <json>]
  antenna_datasheet.py render  [--data <json>]          # markdown
  antenna_datasheet.py --selftest

NAME is a maple token (A30FF, E40, E40a, d30FF, ...); the Fortran entry
point is resolved through notation.map's ant*fortset (fallback:
Full<NAME>). The probe needs a BUILT spike-test object directory
(--testdir) of any process with 7-parton phase space; the antenna is
evaluated on legs 3.. of that phase space, so the datasheet is
process-independent.

Interpretation thresholds (recorded raw, so the reader can re-judge):
  sco pole power   = round(alpha), alpha = -dlog med|A| / dlog x
  soft-capable leg = ss alpha > 1.5   (eikonal alpha = 2)
  singular tc/ds   = alpha > 1.5      (raw alpha kept alongside)
Soft dipole fit: per soft leg, A is least-squares fitted onto the
eikonals 2 s_ac/(s_ab s_bc) of every candidate radiator pair (a,c);
stable rational coefficients across the x scan name the dipole(s).
"""
import argparse
import datetime
import itertools
import json
import os
import re
import subprocess
import sys
import tempfile

# positional slot species per antenna token: 'q' = quark-type (quark or
# antiquark; primes = secondary pair), 'g' = gluon.  This is
# antenna-level physics (the parent-process letter table of
# antennae-naming-convention), NOT process-specific.  Positional means
# "as the Fortran entry point orders its arguments" — for permuted
# declarations (FullC40 declares (i1,i3,i4,i2)) the species follow the
# POSITIONAL order.  None = unknown, measure/read source first.
SPECIES = {
    "A30FF": ["q", "g", "qb"],
    "D30FF": ["q", "g", "g"],
    "d30FF": ["q", "g", "g"],
    "E30FF": ["q", "qp", "qbp"],
    "F30FF": ["g", "g", "g"],
    "f30FF": ["g", "g", "g"],
    "G30FF": ["g", "q", "qb"],
    "A40":   ["q", "g", "g", "qb"],
    "At40":  ["q", "g", "g", "qb"],
    "B40":   ["q", "qp", "qbp", "qb"],
    "C40":   ["q", "q", "qb", "qb"],     # FullC40 declares (i1,i3,i4,i2)
    "D40":   ["q", "g", "g", "g"],
    "E40":   ["q", "qp", "qbp", "g"],
    "E40a":  ["q", "qp", "qbp", "g"],
    # E40b DECLARES (i1,i5,i4,i3): positional species follow the
    # declaration, not the named ascending order — the exact trap the
    # slot audit warns about, here applied to species.
    "E40b":  ["q", "g", "qbp", "qp"],
    "Et40":  ["q", "qp", "qbp", "g"],
    "F40":   ["g", "g", "g", "g"],
    "G40":   ["g", "q", "qb", "g"],
    "Gt40":  ["g", "q", "qb", "g"],
    "H40":   ["g", "q", "qb", "g"],
}

# Slot tag constraints: which slot groups must hold a same-flavour
# quark/antiquark PAIR (opposite kind, equal tag) and which must all
# carry one tag.  'Q' in SPECIES means quark-KIND (quark or antiquark —
# charge conjugation is free); 'g' means gluon.
# Only the UNRESOLVED-PAIR slots are constrained: a g->q qb splitting
# needs a same-flavour particle/antiparticle pair.  Radiator endpoints
# (A30 slots 1,3; A40/B40 slots 1,4) are NOT: colour-ordered
# multi-quark MEs legitimately hang the antenna between quarks of
# different flavours (the antenna is flavour-blind; the reduced ME
# carries the flavours).
SLOT_PAIRS = {
    "E30FF": [(2, 3)],
    "G30FF": [(2, 3)],
    "B40":   [(2, 3)],
    "C40":   [(2, 3)],
    "E40":   [(2, 3)],
    "E40a":  [(2, 3)],
    "E40b":  [(3, 4)],   # declared (i1,i5,i4,i3): the pair sits at 3,4
    "Et40":  [(2, 3)],
    "G40":   [(2, 3)],
    "Gt40":  [(2, 3)],
    "H40":   [(2, 3)],
}

# Reduced-cluster flavour semantics: what parton each generic cluster
# BECOMES.  ("Q", n) = quark kind+tag inherited from the arg in slot n;
# ("g", None) = gluon.  NAIVE NET-FLAVOUR ARITHMETIC IS WRONG for the
# E/G families: E30FF(l,i,j) clusters [l,i],[i,j] where [l,i] carries
# l's flavour (the pair (i,j) is the unresolved system), and
# G30FF(a,b,c)'s [a,b] stays a GLUON although net(g+q)=q.  This table
# is the correct bookkeeping, per antenna, cluster1 = [1,2](3-slot) /
# [1,2,3](4-slot), cluster2 = [2,3] / [4,3,2].
REDUCED = {
    "A30FF": [("Q", 1), ("Q", 3)],
    "D30FF": [("Q", 1), ("g", None)],
    "d30FF": [("Q", 1), ("g", None)],
    "E30FF": [("Q", 1), ("g", None)],
    "F30FF": [("g", None), ("g", None)],
    "f30FF": [("g", None), ("g", None)],
    "G30FF": [("g", None), ("g", None)],
    "A40":   [("Q", 1), ("Q", 4)],
    "At40":  [("Q", 1), ("Q", 4)],
    "B40":   [("Q", 1), ("Q", 4)],
    "C40":   [("Q", 1), ("Q", 4)],
    "D40":   [("Q", 1), ("g", None)],
    "E40":   [("Q", 1), ("g", None)],
    "E40a":  [("Q", 1), ("g", None)],
    "E40b":  [("Q", 1), ("g", None)],
    "Et40":  [("Q", 1), ("g", None)],
    "F40":   [("g", None), ("g", None)],
    "F40a":  [("g", None), ("g", None)],
    "F40b":  [("g", None), ("g", None)],
    "G40":   [("g", None), ("g", None)],
    "Gt40":  [("g", None), ("g", None)],
    "H40":   [("g", None), ("g", None)],
}

# Full composites whose halves are registered separately: candidate
# split identities to VERIFY numerically (never trust a header).
# Second-half permutations tried: X40 -> (A,B,C,D) and (A,D,C,B);
# X30 double -> (A,C,B).
SPLIT_CANDIDATES = {
    "E40":   ("E40a", "E40b"),
    "D30FF": ("d30FF", "d30FF"),
    "F30FF": ("f30FF", "f30FF"),
    "F40":   ("F40a", "F40b"),
}

SOFT_ALPHA_MIN = 1.5
POWER_TOL = 0.30


# ---------------- repo lookups ----------------

def repo_fortmap(root):
    """token -> fortran entry name, from notation.map ant*fortset."""
    path = os.path.join(root, "maple", "notation.map")
    text = open(path, errors="replace").read()
    fortmap = {}
    for m in re.finditer(r"(\w+fortset)\s*:=\s*\{(.*?)\}\s*:", text,
                         re.DOTALL):
        for entry in m.group(2).replace("\n", ",").split(","):
            if "=" in entry:
                tok, fort = entry.split("=", 1)
                fortmap[tok.strip()] = fort.strip()
    return fortmap


def repo_decls(root):
    """fortran name -> (args list, path) for src/X30|X31|X40."""
    func = re.compile(r"^\s{0,6}(?:double\s+precision\s+)?function\s+"
                      r"(\w+)\s*\(([^)]*)\)", re.I | re.M)
    out = {}
    for d in ("src/X30", "src/X31", "src/X40"):
        p = os.path.join(root, d)
        if not os.path.isdir(p):
            continue
        for f in sorted(os.listdir(p)):
            if not f.endswith(".f"):
                continue
            text = open(os.path.join(p, f), errors="replace").read()
            for m in func.finditer(text):
                args = [a.strip() for a in m.group(2).split(",")
                        if a.strip()]
                out[m.group(1)] = (args, os.path.join(d, f))
    return out


def resolve(token, fortmap, decls):
    """token -> (fortran name, n slots)."""
    fort = fortmap.get(token)
    if fort is None or fort not in decls:
        for cand in (token, "Full" + token):
            if cand in decls:
                fort = cand
                break
    if fort is None or fort not in decls:
        raise KeyError(f"no Fortran entry point for {token!r}")
    args, _ = decls[fort]
    nslot = len([a for a in args if re.match(r"^i[1-9]$", a)])
    return fort, nslot


def j2_coefficients(root, token):
    """Lines of J21.map / J22.map mentioning cal<token> (FF block name
    convention: A30FF -> calA30FF, E40 -> calE40FF)."""
    cal = "cal" + token + ("FF" if not token.endswith(("FF", "IF",
                                                       "FI", "II"))
                           else "")
    hits = []
    for fn in ("J21.map", "J22.map"):
        path = os.path.join(root, "maple", "form", "common", fn)
        if not os.path.isfile(path):
            continue
        text = open(path, errors="replace").read()
        # entries are 'NAME(args) = rhs,' possibly multiline
        for m in re.finditer(r"(J2\w+)\([^)]*\)\s*=\s*(.*?),\s*\n\s*\n",
                             text + "\n\n", re.DOTALL):
            rhs = " ".join(m.group(2).split())
            if re.search(r"(?<![\w])" + re.escape(cal) + r"(?![\w])", rhs):
                hits.append(f"{m.group(1)} = {rhs}")
    return hits


# ---------------- probe emission ----------------

PROG_HEAD = """      program dsprobe
c     generated by antenna_datasheet.py -- do not edit
      use Process_mod
      use KinData_mod
      use Scales_mod
      use Mapping_mod
      implicit real*8(a-h,o-y)
      implicit complex*16(z)
      common/plotmode/iplot
      common/koppel/nf,nc,b0
      parameter (nxs={nxs}, npt={npt})
      dimension xslist(nxs)
      dimension v0(100000)
      dimension bmat(3,3),cvec(3),wsol(3)
      data xslist / {xsdata} /
      call init_proc("{process}")
      call init_map()
      call setSqrts_proc(1000d0)
      iplot=2
      nf=5
      nc=3
      call setScales()
      call init_kin(3,10)
"""

PROG_TAIL = """      call destroy_kin()
      stop
      end program dsprobe

      subroutine pmed(v,n,res)
      implicit real*8(a-h,o-z)
      dimension v(n), w(100000)
      res=0d0
      if (n.lt.1) return
      do i=1,n
        w(i)=v(i)
      end do
      do i=2,n
        t=w(i)
        j=i-1
        do while (j.ge.1)
          if (w(j).le.t) exit
          w(j+1)=w(j)
          j=j-1
        end do
        w(j+1)=t
      end do
      if (mod(n,2).eq.1) then
        res=w((n+1)/2)
      else
        res=0.5d0*(w(n/2)+w(n/2+1))
      end if
      return
      end

      subroutine gsolve(b,c,w,n,ok)
      implicit real*8(a-h,o-z)
      dimension b(3,3),c(3),w(3)
      integer ok
      ok=1
      do k=1,n
        ip=k
        do i=k+1,n
          if (abs(b(i,k)).gt.abs(b(ip,k))) ip=i
        end do
        if (abs(b(ip,k)).lt.1d-280) then
          ok=0
          return
        end if
        if (ip.ne.k) then
          do j=1,n
            t=b(k,j)
            b(k,j)=b(ip,j)
            b(ip,j)=t
          end do
          t=c(k)
          c(k)=c(ip)
          c(ip)=t
        end if
        do i=k+1,n
          f=b(i,k)/b(k,k)
          do j=k,n
            b(i,j)=b(i,j)-f*b(k,j)
          end do
          c(i)=c(i)-f*c(k)
        end do
      end do
      do i=n,1,-1
        t=c(i)
        do j=i+1,n
          t=t-b(i,j)*w(j)
        end do
        w(i)=t/b(i,i)
      end do
      return
      end
"""

COLL_MODE = """
      write(*,'(a,a)') '@mode ', '{tag}'
      do ix=1,nxs
        xs=xslist(ix)
        nacc=0
        do ii=1,npt
{gen}
          a = {target}
          call rotp7({ri},{rj})
          a2 = {target}
          if (a.eq.a .and. a2.eq.a2) then
            a=0.5d0*(a+a2)
          endif
          if (a.ne.a) cycle
          nacc=nacc+1
          v0(nacc)=abs(a)
        end do
        call pmed(v0,nacc,r0)
        write(*,'(a,es10.2,i6,es15.6)') '@x ', xs, nacc, r0
      end do
"""

SOFT_MODE = """
      write(*,'(a,a)') '@mode ', '{tag}'
      write(*,'(a,a)') '@cands ', '{cands}'
      do ix=1,nxs
        xs=xslist(ix)
        nacc=0
        do i=1,3
          cvec(i)=0d0
          do j=1,3
            bmat(i,j)=0d0
          end do
        end do
        do ii=1,npt
          em=sqrts_proc*dsqrt(1d0-xs)
          call get_ss7(sqrts_proc,{b},em,{rest})
          a = {target}
          if (a.ne.a) cycle
{eiks}
          wgt=1d0/(a*a+1d-290)
          do i=1,{nc}
            cvec(i)=cvec(i)+wgt*ee{arr}(i)*a
            do j=1,{nc}
              bmat(i,j)=bmat(i,j)+wgt*ee{arr}(i)*ee{arr}(j)
            end do
          end do
          nacc=nacc+1
          v0(nacc)=abs(a)
        end do
        call pmed(v0,nacc,r0)
        call gsolve(bmat,cvec,wsol,{nc},iok)
        write(*,'(a,es10.2,i6,es15.6)') '@x ', xs, nacc, r0
        if (iok.eq.1) then
          write(*,'(a,3es15.6)') '@dip ', (wsol(i),i=1,{nc})
        else
          write(*,'(a)') '@dip degenerate'
        end if
      end do
"""

SPLIT_BLOCK = """
      write(*,'(a,a)') '@split ', '{label}'
      do ii=1,40
        call get_gen7(sqrts_proc,3,4,5,6,7)
        f = {full}
        g = {sum}
        if (f.ne.f .or. g.ne.g) cycle
        write(*,'(a,es15.6)') '@sp ', abs(f-g)/(abs(f)+1d-290)
      end do
"""


def emit_probe(token, fort, nslot, process, npt, xs, splits):
    """Fortran source for the full datasheet measurement of one antenna."""
    legs = [3, 4, 5, 6][:nslot]
    allp = [3, 4, 5, 6, 7]
    target = f"{fort}({','.join(map(str, legs))},7)"
    xsdata = ", ".join(f"{x:.1e}".replace("e", "d") for x in xs)
    out = PROG_HEAD.format(nxs=len(xs), npt=npt, xsdata=xsdata,
                           process=process)

    def rest_of(used):
        return [p for p in allp if p not in used]

    # single collinear
    for i, j in itertools.combinations(legs, 2):
        r = rest_of([i, j])
        gen = (f"          em=sqrts_proc*dsqrt(xs)\n"
               f"          call get_sco7(sqrts_proc,{i},{j},em,"
               f"{r[0]},{r[1]},{r[2]})")
        out += COLL_MODE.format(tag=f"sco {i} {j}", gen=gen,
                                target=target, ri=i, rj=j)
    # single soft + dipole fit
    for b in legs:
        others = [p for p in legs if p != b]
        cands = list(itertools.combinations(others, 2))
        r = rest_of([b])
        eiks = ""
        for k, (a, c) in enumerate(cands, 1):
            eiks += (f"          ee1({k})=2d0*kin(7)%s({a},{c})/"
                     f"(kin(7)%s({a},{b})*kin(7)%s({c},{b}))\n")
        out += SOFT_MODE.format(
            tag=f"ss {b}", b=b,
            cands=" ".join(f"{a},{c}" for a, c in cands),
            rest=",".join(map(str, r)), eiks=eiks.rstrip("\n"),
            nc=len(cands), arr="1", target=target)
    # replace the ee1 array name by reusing v0? simpler: declare ee1
    out = out.replace("      dimension v0(100000)",
                      "      dimension v0(100000)\n"
                      "      dimension ee1(3)")
    if nslot == 4:
        # triple collinear
        for tri in itertools.combinations(legs, 3):
            r = rest_of(list(tri))
            i, j, k = tri
            gen = (f"          em=sqrts_proc*dsqrt(xs)\n"
                   f"          call get_tc7(sqrts_proc,{i},{j},{k},em,"
                   f"{r[0]},{r[1]})")
            out += COLL_MODE.format(tag=f"tc {i} {j} {k}", gen=gen,
                                    target=target, ri=i, rj=j)
        # double soft
        for i, j in itertools.combinations(legs, 2):
            r = rest_of([i, j])
            gen = (f"          em=sqrts_proc*dsqrt(1d0-xs)\n"
                   f"          call get_ds7(sqrts_proc,{i},{j},em,"
                   f"{r[0]},{r[1]},{r[2]})")
            out += COLL_MODE.format(tag=f"ds {i} {j}", gen=gen,
                                    target=target, ri=i, rj=j)
    # split identities
    for label, full_call, sum_call in splits:
        out += SPLIT_BLOCK.format(label=label, full=full_call,
                                  sum=sum_call)
    out += PROG_TAIL
    return out


def split_probe_calls(token, fortmap, decls, nslot):
    """Candidate split identities to test for a Full composite."""
    if token not in SPLIT_CANDIDATES:
        return []
    ha, hb = SPLIT_CANDIDATES[token]
    try:
        fa, _ = resolve(ha, fortmap, decls)
        fb, _ = resolve(hb, fortmap, decls)
        ff, _ = resolve(token, fortmap, decls)
    except KeyError:
        return []
    legs = [3, 4, 5, 6][:nslot]
    A = ",".join(map(str, legs))
    full = f"{ff}({A},7)"
    out = []
    if nslot == 4:
        p, q, r, s = legs
        perms = {"ABCD": (p, q, r, s), "ADCB": (p, s, r, q)}
        for pname, perm in perms.items():
            B = ",".join(map(str, perm))
            out.append((f"{ha}(ABCD)+{hb}({pname})", full,
                        f"{fa}({A},7)+{fb}({B},7)"))
    else:
        p, q, r = legs
        out.append((f"{ha}(ABC)+{hb}(ACB)", full,
                    f"{fa}({A},7)+{fb}({p},{r},{q},7)"))
        out.append((f"{ha}(ABC)+{hb}(BCA)+{hb}(CAB)", full,
                    f"{fa}({A},7)+{fb}({q},{r},{p},7)"
                    f"+{fb}({r},{p},{q},7)"))
    return out


# ---------------- run + parse ----------------

def run_probe(src, testdir, keep=False):
    objdir = os.path.join(testdir, "obj")
    mods = os.path.join(objdir, "mod")
    with tempfile.TemporaryDirectory() as td:
        fsrc = os.path.join(td, "dsprobe.f")
        open(fsrc, "w").write(src)
        subprocess.run(
            ["gfortran", "-ffixed-line-length-none", "-c", fsrc,
             "-J", mods, "-I", mods, "-o", os.path.join(td, "dsprobe.o")],
            check=True, capture_output=True, text=True)
        objs = [os.path.join(objdir, f) for f in sorted(os.listdir(objdir))
                if f.endswith(".o") and not f.startswith("check")
                and not f.startswith("probe")]
        subprocess.run(
            ["gfortran", "-o", os.path.join(td, "dsprobe"),
             os.path.join(td, "dsprobe.o")] + objs,
            check=True, capture_output=True, text=True)
        env = dict(os.environ, OMP_NUM_THREADS="1")
        res = subprocess.run([os.path.join(td, "dsprobe")], env=env,
                             capture_output=True, text=True, timeout=1800)
        if keep:
            kept = os.path.join(testdir, "dsprobe_kept.f")
            open(kept, "w").write(src)
        return res.stdout


def parse_probe(text, xs):
    """probe stdout -> {mode_tag: {"med": [...], "n": [...],
    "cands": [...], "dip": [[...]...], "split": med_rel}}."""
    out = {}
    cur = None
    for line in text.splitlines():
        line = line.strip()
        if line.startswith("@mode "):
            cur = line[6:].strip()
            out[cur] = {"med": [], "n": [], "dip": [], "cands": None}
        elif line.startswith("@cands ") and cur:
            out[cur]["cands"] = line[7:].split()
        elif line.startswith("@x ") and cur:
            parts = line[3:].split()
            out[cur]["n"].append(int(parts[1]))
            out[cur]["med"].append(float(parts[2]))
        elif line.startswith("@dip ") and cur:
            rest = line[5:].split()
            if rest == ["degenerate"]:
                out[cur]["dip"].append(None)
            else:
                out[cur]["dip"].append([float(v) for v in rest])
        elif line.startswith("@split "):
            cur = "split " + line[7:].strip()
            out[cur] = {"sp": []}
        elif line.startswith("@sp ") and cur and cur.startswith("split"):
            out[cur]["sp"].append(float(line[4:]))
    return out


def alpha_of(meds, xs):
    """log-log slope of med|A| vs x (positive = singular)."""
    pts = [(x, m) for x, m in zip(xs, meds) if m > 0]
    if len(pts) < 2:
        return None
    import math
    sl = []
    for (x1, m1), (x2, m2) in zip(pts, pts[1:]):
        sl.append(-(math.log10(m2) - math.log10(m1))
                  / (math.log10(x2) - math.log10(x1)))
    sl.sort()
    return sl[len(sl) // 2]


def snap(v, tol=0.01):
    """Snap a fitted coefficient to the nearest small rational."""
    from fractions import Fraction
    for fr in (Fraction(0), Fraction(1), Fraction(2), Fraction(1, 2),
               Fraction(1, 3), Fraction(2, 3), Fraction(1, 4),
               Fraction(1, 9), Fraction(-1), Fraction(-1, 2)):
        if abs(v - float(fr)) < tol * max(1.0, abs(float(fr))):
            return str(fr)
    return None


def digest(parsed, xs, nslot):
    """parsed probe output -> datasheet 'measured' dict."""
    legs = [3, 4, 5, 6][:nslot]
    slot_of = {p: i + 1 for i, p in enumerate(legs)}
    meas = {"sco": {}, "ss": {}, "tc": {}, "ds": {}, "splits": {}}
    for tag, d in parsed.items():
        parts = tag.split()
        if parts[0] == "sco":
            i, j = int(parts[1]), int(parts[2])
            al = alpha_of(d["med"], xs)
            if al is None or al < 0.5:
                continue
            key = f"{slot_of[i]},{slot_of[j]}"
            meas["sco"][key] = {"alpha": round(al, 2),
                                "power": int(round(al))}
        elif parts[0] == "ss":
            b = int(parts[1])
            al = alpha_of(d["med"], xs)
            if al is None or al < SOFT_ALPHA_MIN:
                continue
            entry = {"alpha": round(al, 2), "dipoles": {}}
            if nslot == 4:
                entry["note"] = ("soft limit = eik x reduced X30, not "
                                 "eik x const: fit the residue with "
                                 "probe-me-ir-structure mode 3")
            if d["cands"] and d["dip"] and nslot == 3:
                last = next((v for v in reversed(d["dip"])
                             if v is not None), None)
                if last:
                    for cand, coeff in zip(d["cands"], last):
                        r = snap(coeff)
                        if r not in (None, "0"):
                            a, c = (int(v) for v in cand.split(","))
                            entry["dipoles"][
                                f"{slot_of[a]},{slot_of[c]}"] = r
                        elif r is None and abs(coeff) > 0.02:
                            a, c = (int(v) for v in cand.split(","))
                            entry["dipoles"][
                                f"{slot_of[a]},{slot_of[c]}"] = \
                                f"~{coeff:.3f} (not rational)"
            meas["ss"][str(slot_of[b])] = entry
        elif parts[0] == "tc":
            i, j, k = (int(v) for v in parts[1:4])
            al = alpha_of(d["med"], xs)
            if al is None or al < SOFT_ALPHA_MIN:
                continue
            key = ",".join(str(slot_of[v]) for v in (i, j, k))
            meas["tc"][key] = {"alpha": round(al, 2)}
        elif parts[0] == "ds":
            i, j = int(parts[1]), int(parts[2])
            al = alpha_of(d["med"], xs)
            if al is None or al < SOFT_ALPHA_MIN:
                continue
            key = f"{slot_of[i]},{slot_of[j]}"
            meas["ds"][key] = {"alpha": round(al, 2)}
        elif parts[0] == "split":
            sp = sorted(d["sp"])
            if sp:
                med = sp[len(sp) // 2]
                meas["splits"][" ".join(parts[1:])] = {
                    "median_rel_residual": f"{med:.2e}",
                    "holds": bool(med < 1e-9)}
    return meas


def derive_split_flag(entry):
    """requires_split for 4-slot antennae: any sco pole not supported
    by the generic clusters [1,2,3],[4,3,2]."""
    if entry.get("arity") != 4:
        entry["requires_split"] = False
        return
    supported = {frozenset(p) for p in
                 [(1, 2), (1, 3), (2, 3), (4, 3), (4, 2)]}
    bad = []
    meas = entry.get("measured", {})
    for key in meas.get("sco", {}):
        pair = frozenset(int(v) for v in key.split(","))
        if pair not in supported:
            bad.append(key)
    for b, d in meas.get("ss", {}).items():
        for dip in d.get("dipoles", {}):
            pass  # soft dipoles need no cluster support of their own
    entry["requires_split"] = bool(bad)
    entry["unsupported_poles"] = bad


# ---------------- commands ----------------

def default_data_path():
    return os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "..", "assets", "antenna_datasheet.json")


def load_data(path):
    if os.path.isfile(path):
        return json.load(open(path))
    return {}


def save_data(path, data):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    json.dump(data, open(path, "w"), indent=1, sort_keys=True)


def cmd_measure(args):
    root = os.path.abspath(args.root)
    fortmap = repo_fortmap(root)
    decls = repo_decls(root)
    data = load_data(args.data)
    xs = [float(v) for v in args.xs.split(",")]
    for token in args.names:
        fort, nslot = resolve(token, fortmap, decls)
        splits = split_probe_calls(token, fortmap, decls, nslot)
        src = emit_probe(token, fort, nslot, args.process, args.npt,
                         xs, splits)
        print(f"[{token}] fortran={fort} slots={nslot} ... ",
              end="", flush=True)
        try:
            outtxt = run_probe(src, os.path.join(root, args.testdir),
                               keep=args.keep)
        except subprocess.CalledProcessError as e:
            print("FAILED to build probe")
            sys.stderr.write((e.stderr or "")[-2000:])
            continue
        parsed = parse_probe(outtxt, xs)
        meas = digest(parsed, xs, nslot)
        entry = data.setdefault(token, {})
        entry.update({
            "fortran": fort, "arity": nslot,
            "slots_declared": ",".join(decls[fort][0]),
            "source": decls[fort][1],
            "species": SPECIES.get(token),
            "measured": meas,
            "provenance": {
                "date": datetime.date.today().isoformat(),
                "testdir": args.testdir, "process": args.process,
                "npt": args.npt, "xs": xs},
        })
        derive_split_flag(entry)
        if token in SPLIT_CANDIDATES:
            ha, hb = SPLIT_CANDIDATES[token]
            holds = [k for k, v in meas["splits"].items() if v["holds"]]
            entry["halves"] = [ha, hb]
            entry["split_identity"] = holds[0] if holds else \
                "NONE OF THE TESTED IDENTITIES HOLDS"
        save_data(args.data, data)
        npol = (len(meas['sco']), len(meas['ss']), len(meas['tc']),
                len(meas['ds']))
        print(f"sco={npol[0]} ss={npol[1]} tc={npol[2]} ds={npol[3]}"
              + (f"  split: {entry.get('split_identity')}"
                 if token in SPLIT_CANDIDATES else "")
              + ("  REQUIRES SPLIT" if entry.get("requires_split")
                 else ""))


def cmd_static(args):
    root = os.path.abspath(args.root)
    data = load_data(args.data)
    for token, entry in data.items():
        entry["j2_coefficients"] = j2_coefficients(root, token)
        # the tables are authoritative for declared data — always refresh
        entry["species"] = SPECIES.get(token, entry.get("species"))
        entry["slot_pairs"] = SLOT_PAIRS.get(token, [])
        entry["reduced"] = REDUCED.get(token)
    save_data(args.data, data)
    print(f"enriched {len(data)} entries "
          f"(J21/J22 + species + reduced-cluster semantics)")


def fmt_entry(token, e):
    out = [f"=== {token} ===",
           f"  fortran   : {e.get('fortran')}"
           f"({e.get('slots_declared')})   [{e.get('source')}]",
           f"  species   : {e.get('species')}"]
    m = e.get("measured", {})
    sco = m.get("sco", {})
    out.append("  collinear poles (slot pairs): "
               + (", ".join(f"s({k}) power {v['power']}"
                            for k, v in sorted(sco.items())) or "none"))
    ss = m.get("ss", {})
    for b, d in sorted(ss.items()):
        dips = ", ".join(f"({k}) coeff {v}"
                         for k, v in d["dipoles"].items()) or "?"
        out.append(f"  soft leg {b}: alpha {d['alpha']}, eikonal on"
                   f" dipole {dips}")
    if not ss:
        out.append("  soft legs : none")
    if e.get("arity") == 4:
        out.append("  tc singular triples: "
                   + (", ".join(sorted(m.get("tc", {}))) or "none"))
        out.append("  ds singular pairs  : "
                   + (", ".join(sorted(m.get("ds", {}))) or "none"))
        out.append(f"  requires_split: {e.get('requires_split')}"
                   + (f"  (unsupported poles: "
                      f"{e.get('unsupported_poles')})"
                      if e.get("requires_split") else ""))
    if e.get("halves"):
        out.append(f"  halves: {e['halves']}   verified identity: "
                   f"{e.get('split_identity')}")
    for line in e.get("j2_coefficients", []):
        out.append(f"  J2: {line}")
    pr = e.get("provenance", {})
    out.append(f"  measured {pr.get('date')} in {pr.get('testdir')}"
               f" (npt={pr.get('npt')}, xs={pr.get('xs')})")
    return "\n".join(out)


def cmd_show(args):
    data = load_data(args.data)
    for token in args.names:
        if token not in data:
            print(f"{token}: not in datasheet "
                  f"(run 'measure {token}' first)")
            continue
        print(fmt_entry(token, data[token]))
        print()


def cmd_render(args):
    data = load_data(args.data)
    print("# Antenna datasheet (generated — do not edit by hand)\n")
    for token in sorted(data):
        print("```")
        print(fmt_entry(token, data[token]))
        print("```\n")


def cmd_verify(args):
    """Re-measure and diff against the stored entry."""
    data = load_data(args.data)
    stored = {t: json.loads(json.dumps(data.get(t))) for t in args.names}
    cmd_measure(args)
    fresh = load_data(args.data)
    bad = 0
    for t in args.names:
        keys = ("measured", "requires_split", "split_identity")
        diffs = [k for k in keys
                 if json.dumps(stored[t].get(k) if stored[t] else None,
                               sort_keys=True)
                 != json.dumps(fresh[t].get(k), sort_keys=True)]
        if stored[t] is None:
            print(f"{t}: no stored entry — now measured")
        elif diffs:
            bad += 1
            print(f"{t}: CHANGED fields {diffs} — the datasheet was "
                  f"stale (or the antenna source changed)")
        else:
            print(f"{t}: verified, unchanged")
    sys.exit(1 if bad else 0)


# ---------------- selftest ----------------

def selftest():
    """Parser + derivation logic on synthetic output. Encodes NO real
    antenna's pole graph — fake names, fake numbers."""
    xs = [1e-5, 1e-6, 1e-7]
    fake = "\n".join([
        "@mode sco 3 4",
        "@x  1.0E-05    40   1.0E+05",
        "@x  1.0E-06    40   1.0E+06",
        "@x  1.0E-07    40   1.0E+07",
        "@mode sco 3 6",       # endpoint pair of a 4-slot antenna
        "@x  1.0E-05    40   2.0E+05",
        "@x  1.0E-06    40   2.0E+06",
        "@x  1.0E-07    40   2.0E+07",
        "@mode sco 4 5",
        "@x  1.0E-05    40   3.0E+00",   # not singular
        "@x  1.0E-06    40   3.0E+00",
        "@x  1.0E-07    40   3.0E+00",
        "@mode ss 4",
        "@cands 3,5 3,6 5,6",
        "@x  1.0E-05    40   1.0E+10",
        "@dip  1.0000001 0.0000001 0.4999999",
        "@x  1.0E-06    40   1.0E+12",
        "@dip  0.9999999 -0.0000001 0.5000001",
        "@x  1.0E-07    40   1.0E+14",
        "@dip  1.0000000 0.0000000 0.5000000",
        "@mode tc 3 4 5",
        "@x  1.0E-05    40   1.0E+10",
        "@x  1.0E-06    40   1.0E+12",
        "@x  1.0E-07    40   1.0E+14",
        "@mode ds 4 5",
        "@x  1.0E-05    40   1.0E+02",   # dead
        "@x  1.0E-06    40   1.0E+02",
        "@x  1.0E-07    40   1.0E+02",
        "@split halfA(ABCD)+halfB(ADCB)",
        "@sp  1.2E-13", "@sp  0.8E-13", "@sp  2.0E-13",
        "@split halfA(ABCD)+halfB(ABCD)",
        "@sp  0.3E+00", "@sp  0.4E+00", "@sp  0.5E+00",
    ])
    parsed = parse_probe(fake, xs)
    meas = digest(parsed, xs, 4)
    assert meas["sco"]["1,2"]["power"] == 1
    assert meas["sco"]["1,4"]["power"] == 1          # endpoint pair
    assert "2,3" not in meas["sco"]                  # non-singular
    assert meas["ss"]["2"]["alpha"] == 2.0
    # arity 4: dipole fit deferred to the residue fitter, note recorded
    assert meas["ss"]["2"]["dipoles"] == {}
    assert "mode 3" in meas["ss"]["2"]["note"]
    # arity 3: the dipole fit IS recorded, coefficients snapped
    fake3 = "\n".join([
        "@mode ss 4",
        "@cands 3,5",
        "@x  1.0E-05    40   1.0E+10",
        "@dip  0.9999999",
        "@x  1.0E-06    40   1.0E+12",
        "@dip  1.0000001",
        "@x  1.0E-07    40   1.0E+14",
        "@dip  1.0000000",
    ])
    meas3 = digest(parse_probe(fake3, xs), xs, 3)
    assert meas3["ss"]["2"]["dipoles"] == {"1,3": "1"}
    assert "1,2,3" in meas["tc"]
    assert meas["ds"] == {}                          # dead pair dropped
    assert meas["splits"]["halfA(ABCD)+halfB(ADCB)"]["holds"]
    assert not meas["splits"]["halfA(ABCD)+halfB(ABCD)"]["holds"]

    entry = {"arity": 4, "measured": meas}
    derive_split_flag(entry)
    assert entry["requires_split"] and \
        entry["unsupported_poles"] == ["1,4"]
    # supported-only pole set -> no split needed
    entry2 = {"arity": 4, "measured":
              {"sco": {"1,2": {}, "2,3": {}, "3,4": {}}, "ss": {}}}
    derive_split_flag(entry2)
    assert not entry2["requires_split"]

    # coefficient snapping
    assert snap(0.667) == "2/3" and snap(0.5001) == "1/2"
    assert snap(0.7) is None and snap(-0.0000003) == "0"
    print("antenna_datasheet selftest OK")


def main():
    if "--selftest" in sys.argv:
        selftest()
        return
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)

    def common(p, names=True):
        if names:
            p.add_argument("names", nargs="+")
        p.add_argument("--data", default=default_data_path())
        p.add_argument("--root", default=".")

    p = sub.add_parser("measure")
    common(p)
    p.add_argument("--testdir", required=True)
    p.add_argument("--process", default="epem")
    p.add_argument("--npt", type=int, default=40)
    p.add_argument("--xs", default="1e-5,1e-6,1e-7")
    p.add_argument("--keep", action="store_true")
    p.set_defaults(func=cmd_measure)

    p = sub.add_parser("static")
    common(p, names=False)
    p.set_defaults(func=cmd_static)

    p = sub.add_parser("show")
    common(p)
    p.set_defaults(func=cmd_show)

    p = sub.add_parser("render")
    common(p, names=False)
    p.set_defaults(func=cmd_render)

    p = sub.add_parser("verify")
    common(p)
    p.add_argument("--testdir", required=True)
    p.add_argument("--process", default="epem")
    p.add_argument("--npt", type=int, default=40)
    p.add_argument("--xs", default="1e-5,1e-6,1e-7")
    p.add_argument("--keep", action="store_true")
    p.set_defaults(func=cmd_verify)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
