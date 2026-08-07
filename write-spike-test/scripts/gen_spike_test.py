#!/usr/bin/env python3
"""Generate a complete NNLOJET spike-test program (check<N>to<M>[loop].f)
from a short JSON spec.

Usage:  python gen_spike_test.py spec.json > check5to3.f
        python gen_spike_test.py --selftest

The mode table for n final-state partons is a mechanical function of n
(single soft, single collinear, and for RR: double soft, triple
collinear, soft-collinear, double collinear); this generator builds it
programmatically, marks every mode GENUINE/dead via the reduced-Born
rule (genuine_modes.py), applies azimuthal averaging on collinear modes
with a gluon parent, and reports per (mode, x) a summary line —
n / nan / max|ME| / median / min / max / outlier count — not raw event
dumps. Hand-writing 80 case blocks invites silent indexing errors;
generating them does not.

Spec fields:
  process      : init_proc name, e.g. "epem"
  sqrts        : e.g. 1000.0
  init_kin     : [npartons, arg2] as in existing check programs
  npar         : kinematics-set index = total momenta (6..8)
  nborn        : Born-level final-state parton count (for the name)
  contribution : "R" | "RV" | "RR" (R/RV: single-unresolved families
                 only; RV additionally names the file ...loop.f)
  fs_partons   : {"<momentum index>": "<flavour>"} — flavours as in
                 genuine_modes.py ("g", "q<tag>", "qb<tag>")
  born         : list of legal Born flavour lists (genuine_modes.py)
  decl_lines   : verbatim DECLARATION lines for the main program
  setup_lines  : verbatim EXECUTABLE init lines (must leave iplot!=0)
  cuts_call    : verbatim, e.g. "call ecuts_epem(1,7,ipass)"
  channels     : [{"itype": 1, "comment": "...", "me": "<full ME expr>",
                   "sub": "<subtraction expr>"}]  — identical momentum
                 arguments in me and sub
  test_use_lines / test_decl_lines : optional, for the test()/tests()
                 functions (default: same modules as the main program)
  xs_soft      : x values for soft-driven modes  (default 1e-5..1e-7)
  xs_coll      : x values for collinear modes    (default 1e-5..1e-7)

Build with the makefile pattern in the write-spike-test skill
(-ffixed-line-length-none is assumed).
"""
import itertools
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from genuine_modes import classify, cluster_parent  # noqa: E402

RR_FAMS = ["ds", "tc", "sc", "dc", "ss", "sco"]
R_FAMS = ["ss", "sco"]


def enumerate_modes(fs, fams):
    """[(family, name, unresolved, spectators, rot1, rot2)] — names must
    match genuine_modes.classify exactly."""
    idx = sorted(int(i) for i in fs)
    fsset = set(idx)
    out = []

    def g_parent(pair):
        return cluster_parent([fs[pair[0]], fs[pair[1]]]) == "g"

    for fam in fams:
        if fam == "ds":
            for i, j in itertools.combinations(idx, 2):
                out.append((fam, f"{i},{j} soft", (i, j),
                            tuple(sorted(fsset - {i, j})), None, None))
        elif fam == "tc":
            for i, j, k in itertools.combinations(idx, 3):
                out.append((fam, f"{i}||{j}||{k}", (i, j, k),
                            tuple(sorted(fsset - {i, j, k})),
                            (i, j) if cluster_parent(
                                [fs[i], fs[j], fs[k]]) == "g" else None,
                            None))
        elif fam == "sc":
            for i in idx:
                for j, k in itertools.combinations(
                        [x for x in idx if x != i], 2):
                    out.append((fam, f"{i} soft + {j}||{k}", (i, j, k),
                                tuple(sorted(fsset - {i, j, k})),
                                (j, k) if g_parent((j, k)) else None, None))
        elif fam == "dc":
            for (i, j), (k, l) in itertools.combinations(
                    itertools.combinations(idx, 2), 2):
                if {i, j} & {k, l}:
                    continue
                out.append((fam, f"{i}||{j} + {k}||{l}", (i, j, k, l),
                            tuple(sorted(fsset - {i, j, k, l})),
                            (i, j) if g_parent((i, j)) else None,
                            (k, l) if g_parent((k, l)) else None))
        elif fam == "ss":
            for i in idx:
                out.append((fam, f"{i} soft", (i,),
                            tuple(sorted(fsset - {i})), None, None))
        elif fam == "sco":
            for i, j in itertools.combinations(idx, 2):
                out.append((fam, f"{i}||{j}", (i, j),
                            tuple(sorted(fsset - {i, j})),
                            (i, j) if g_parent((i, j)) else None, None))
    return out


def gen_call(fam, npar, u, sp):
    def call(name, *args):
        # drop empty tails (e.g. dc with zero spectators at npar=6)
        alist = ",".join(str(a) for a in args if a != "")
        return f"call {name}(sqrts_proc,{alist})"

    S = ",".join(str(x) for x in sp)
    if fam == "ss":
        return ["em1=sqrts_proc*dsqrt(1d0-xs)",
                call(f"get_ss{npar}", u[0], "em1", S)]
    if fam == "sco":
        return ["em1=sqrts_proc*xs",
                call(f"get_sco{npar}", u[0], u[1], "em1", S)]
    if fam == "ds":
        return ["em1=sqrts_proc*dsqrt(1d0-xs)",
                call(f"get_ds{npar}", u[0], u[1], "em1", S)]
    if fam == "tc":
        return ["em1=sqrts_proc*xs",
                call(f"get_tc{npar}", u[0], u[1], u[2], "em1", S)]
    if fam == "sc":  # i soft, j||k
        return ["em1=sqrts_proc*dsqrt(1d0-xs)",
                "em2=sqrts_proc*xs",
                call(f"get_sc{npar}", u[0], u[1], u[2], "em1", "em2", S)]
    if fam == "dc":  # i||j , k||l
        return ["em1=sqrts_proc*xs",
                "em2=sqrts_proc*xs",
                call(f"get_dc{npar}", u[0], u[1], "em1",
                     u[2], u[3], "em2", S)]
    raise ValueError(fam)


def f77(lines, indent="      "):
    return "\n".join(indent + ln if ln.strip() else ln for ln in lines)


def fmt_x(x):
    return f"{x:.1e}".replace("e", "d")


def generate(spec):
    fs = {int(k): v for k, v in spec["fs_partons"].items()}
    npar = spec["npar"]
    nfs = len(fs)
    contrib = spec["contribution"]
    fams = RR_FAMS if contrib == "RR" else R_FAMS
    modes = enumerate_modes(fs, fams)
    cls = {(m["family"], m["name"]): m
           for m in classify(fs, spec["born"], fams)}
    xs_soft = spec.get("xs_soft", [1e-5, 1e-6, 1e-7])
    xs_coll = spec.get("xs_coll", [1e-5, 1e-6, 1e-7])
    chans = spec["channels"]
    prog = f"check{nfs}to{spec['nborn']}" + \
        ("loop" if contrib == "RV" else "")

    L = []
    L.append(f"      program {prog}")
    L.append("c     generated by gen_spike_test.py"
             " -- see write-spike-test skill")
    L.append(f"c     process {spec['process']}, contribution {contrib},"
             f" {len(modes)} modes, {len(chans)} channel(s)")
    L.append("      use Process_mod")
    L.append("      use KinData_mod")
    L.append("      use Scales_mod")
    L.append("      use Mapping_mod")
    L.append("      implicit real*8(a-h,o-y)")
    L.append("      implicit complex*16(z)")
    L.append("      common/plotmode/iplot")
    L.append(f"      parameter (nmodes={len(modes)})")
    L.append("      dimension x(3)")
    L.append("      dimension vrat(100000)")
    L.append("      character*64 stitle")
    L.append("      character*32 arg")
    L.extend("      " + ln for ln in spec.get("decl_lines", []))
    L.append("")
    L.append("c     CLI: ./" + prog + " ITYPE [MODELO MODEHI] [IPOINT] [ILOW]")
    L.append("      nargs = iargc()")
    L.append("      if (nargs.lt.1) then")
    L.append(f"        write(*,*) 'usage: ./{prog} ITYPE"
             " [MODELO MODEHI] [IPOINT] [ILOW]'")
    L.append("        write(*,*) 'channels:'")
    for c in chans:
        L.append(f"        write(*,*) '  {c['itype']}: "
                 f"{c.get('comment', '')}'")
    L.append("        stop")
    L.append("      end if")
    L.append("      call getarg(1,arg)")
    L.append("      read(arg,*) itype")
    L.append("      modelo=1")
    L.append("      modehi=nmodes")
    L.append("      ipoint=100")
    L.append("      ilow=1")
    L.append("      if (nargs.ge.3) then")
    L.append("        call getarg(2,arg)")
    L.append("        read(arg,*) modelo")
    L.append("        call getarg(3,arg)")
    L.append("        read(arg,*) modehi")
    L.append("      end if")
    L.append("      if (nargs.ge.4) then")
    L.append("        call getarg(4,arg)")
    L.append("        read(arg,*) ipoint")
    L.append("      end if")
    L.append("      if (nargs.ge.5) then")
    L.append("        call getarg(5,arg)")
    L.append("        read(arg,*) ilow")
    L.append("      end if")
    L.append("      if (ipoint.gt.100000) ipoint=100000")
    L.append("      if (modelo.lt.1) modelo=1")
    L.append("      if (modehi.gt.nmodes) modehi=nmodes")
    L.append("")
    L.append(f'      call init_proc("{spec["process"]}")')
    L.append("      call init_map()")
    L.append(f"      call setSqrts_proc({spec['sqrts']}d0)")
    L.append("      iplot=2")
    L.append("      call setScales()")
    L.append(f"      call init_kin({spec['init_kin'][0]},"
             f"{spec['init_kin'][1]})")
    L.extend("      " + ln for ln in spec.get("setup_lines", []))
    L.append("")
    L.append("      do mode=modelo,modehi")
    L.append("c       -- mode metadata (title, x set, rotations) --")
    L.append("        mroti=0")
    L.append("        mrotj=0")
    L.append("        mrt2i=0")
    L.append("        mrt2j=0")
    L.append("        select case(mode)")
    for m, (fam, name, u, sp, rot1, rot2) in enumerate(modes, 1):
        c = cls[(fam, name)]
        tag = "GENUINE" if c["genuine"] else "dead"
        xset = list(xs_soft if fam in ("ss", "ds", "sc") else xs_coll)
        xset = (xset + [xset[-1]] * 3)[:3]     # always fill x(1..3)
        L.append(f"        case({m})")
        L.append(f"          stitle='{name} [{tag}]'")
        L.append(f"c         {c['reason']}")
        for k, xv in enumerate(xset, 1):
            L.append(f"          x({k})={fmt_x(xv)}")
        if rot1:
            L.append(f"          mroti={rot1[0]}")
            L.append(f"          mrotj={rot1[1]}")
        if rot2:
            L.append(f"          mrt2i={rot2[0]}")
            L.append(f"          mrt2j={rot2[1]}")
    L.append("        end select")
    L.append("        write(*,*)")
    L.append("        write(*,'(a,i3,2a)') 'mode ', mode, ': ', stitle")
    L.append("")
    L.append("        do ix=ilow,3")
    L.append("          xs=x(ix)")
    L.append("          nacc=0")
    L.append("          nnan=0")
    L.append("          wmax=0d0")
    L.append("          do ii=1,ipoint")
    L.append("            select case(mode)")
    for m, (fam, name, u, sp, rot1, rot2) in enumerate(modes, 1):
        L.append(f"            case({m})")
        for ln in gen_call(fam, npar, u, sp):
            L.append("              " + ln)
    L.append("            end select")
    L.append("            ipass=1")
    L.append("            " + spec["cuts_call"])
    L.append("            if (ipass.ne.1) cycle")
    L.append("            wt1=test(itype)")
    L.append("            wt2=tests(itype)")
    L.append("c           azimuthal average (pi/2 rotation about the")
    L.append("c           collinear axis) for gluon-parent clusters")
    L.append("            if (mroti.ne.0) then")
    L.append(f"              call rotp{npar}(mroti,mrotj)")
    L.append("              wt1=0.5d0*(wt1+test(itype))")
    L.append("              wt2=0.5d0*(wt2+tests(itype))")
    L.append("            end if")
    L.append("            if (mrt2i.ne.0) then")
    L.append(f"              call rotp{npar}(mrt2i,mrt2j)")
    L.append("              wt1=0.5d0*(wt1+test(itype))")
    L.append("              wt2=0.5d0*(wt2+tests(itype))")
    L.append("            end if")
    L.append("            if (wt1.ne.wt1 .or. wt2.ne.wt2) then")
    L.append("              nnan=nnan+1")
    L.append("              cycle")
    L.append("            end if")
    L.append("            if (abs(wt1).gt.wmax) wmax=abs(wt1)")
    L.append("            if (abs(wt2).lt.1d-300) cycle")
    L.append("            nacc=nacc+1")
    L.append("            vrat(nacc)=wt1/wt2")
    L.append("          end do")
    L.append("          if (nacc.lt.3) then")
    L.append("            write(*,'(a,es9.1,a)') '  x=', xs,")
    L.append("     .        '  too few accepted points'")
    L.append("            cycle")
    L.append("          end if")
    L.append("          call pmedian(vrat,nacc,rmed)")
    L.append("          rmin=vrat(1)")
    L.append("          rmax=vrat(1)")
    L.append("          nout=0")
    L.append("          do k=1,nacc")
    L.append("            if (vrat(k).lt.rmin) rmin=vrat(k)")
    L.append("            if (vrat(k).gt.rmax) rmax=vrat(k)")
    L.append("            if (abs(vrat(k)-rmed).gt.5d-2) nout=nout+1")
    L.append("          end do")
    L.append("          write(*,'(a,es9.1,a,i5,a,i4,a,es10.2,a,f12.6,")
    L.append("     .      a,es10.2,a,es10.2,a,i4)')")
    L.append("     .      '  x=',xs,' n=',nacc,' nan=',nnan,")
    L.append("     .      ' max|ME|=',wmax,' med=',rmed,")
    L.append("     .      ' min=',rmin,' max=',rmax,' outl=',nout")
    L.append("        end do")
    L.append("      end do")
    L.append("")
    L.append("      call destroy_kin()")
    L.append("      stop")
    L.append(f"      end program {prog}")

    # test()/tests() channel functions
    use_lines = spec.get("test_use_lines",
                         ["use Process_mod", "use KinData_mod"])
    tdecl = spec.get("test_decl_lines", [])
    for fname, key in (("test", "me"), ("tests", "sub")):
        L.append("")
        L.append("*" * 72)
        L.append(f"*     {fname}(itype): "
                 + ("full ME" if key == "me" else "subtraction")
                 + " per channel")
        L.append("*" * 72)
        L.append(f"      function {fname}(itype)")
        L.extend("      " + ln for ln in use_lines)
        L.append("      implicit real*8(a-h,o-y)")
        L.append("      implicit complex*16(z)")
        L.extend("      " + ln for ln in tdecl)
        L.append("      select case(itype)")
        for c in chans:
            L.append(f"      case({c['itype']})")
            if c.get("comment"):
                L.append(f"c       {c['comment']}")
            L.append(f"        {fname} = {c[key]}")
        L.append("      case default")
        L.append(f"        write(*,*) 'unknown channel ', itype")
        L.append("        stop")
        L.append("      end select")
        L.append("      return")
        L.append("      end")

    L.append("""
************************************************************************
*     median of v(1:n) (insertion sort on a copy)
************************************************************************
      subroutine pmedian(v,n,res)
      implicit real*8(a-h,o-z)
      dimension v(n), w(100000)
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
      end""")
    return "\n".join(L) + "\n"


def selftest():
    """Structural/combinatoric self-test — no expected physics answer.
    Verifies the mechanical mode count (e.g. 80 for 5 partons at RR)
    and the well-formedness of the emitted program."""
    from math import comb
    fs5 = {"3": "q1", "4": "qb1", "5": "g", "6": "q2", "7": "qb2"}
    spec = {"process": "p", "sqrts": 1000.0, "init_kin": [5, 10],
            "npar": 7, "nborn": 3, "contribution": "RR",
            "fs_partons": fs5, "born": [["q1", "qb1", "g"]],
            "decl_lines": [], "setup_lines": ["continue"],
            "cuts_call": "call dummy_cuts(1,7,ipass)",
            "channels": [{"itype": 1, "comment": "chan one",
                          "me": "dummy_me(3,4,5,6,7)",
                          "sub": "dummy_sub(3,4,5,6,7)"}]}
    n = 5
    modes = enumerate_modes({int(k): v for k, v in fs5.items()}, RR_FAMS)
    cnt = {}
    for fam, *_ in modes:
        cnt[fam] = cnt.get(fam, 0) + 1
    assert cnt["ss"] == comb(n, 1) and cnt["sco"] == comb(n, 2)
    assert cnt["ds"] == comb(n, 2) and cnt["tc"] == comb(n, 3)
    assert cnt["sc"] == n * comb(n - 1, 2)
    assert cnt["dc"] == comb(n, 2) * comb(n - 2, 2) // 2
    assert len(modes) == 80, len(modes)   # the n=5 identity
    src = generate(spec)
    nchan = len(spec["channels"])
    assert src.count("case(") >= 2 * len(modes) + 2 * nchan, \
        "missing case blocks"
    assert src.count(f"call get_ds7") == cnt["ds"]
    assert src.count(f"call get_tc7") == cnt["tc"]
    assert src.count("function test(") == 1
    assert src.count("function tests(") == 1
    assert "end program check5to3" in src
    for fname in ("test", "tests"):
        assert f"{fname} = dummy_" in src
    # R contribution: single-unresolved families only
    spec_r = dict(spec, contribution="R", nborn=4)
    modes_r = enumerate_modes({int(k): v for k, v in fs5.items()}, R_FAMS)
    assert len(modes_r) == comb(n, 1) + comb(n, 2)
    src_r = generate(spec_r)
    assert "call get_ds7" not in src_r and "call get_tc7" not in src_r
    assert "end program check5to4" in src_r
    # 4-parton RR at npar=6 exercises the zero-spectator dc path
    fs4 = {"3": "q1", "4": "qb1", "5": "g", "6": "g"}
    spec6 = dict(spec, npar=6, init_kin=[4, 10], fs_partons=fs4,
                 cuts_call="call dummy_cuts(1,6,ipass)",
                 channels=[{"itype": 1, "me": "dummy_me(3,4,5,6)",
                            "sub": "dummy_sub(3,4,5,6)"}])
    src6 = generate(spec6)
    assert ",)" not in src6, "trailing comma in a generated call"
    assert "call get_dc6" in src6
    # every emitted line: fixed-form comment or starts at col 7+ or label
    for ln in src.splitlines():
        if ln and not ln.startswith(("c", "*", " ")):
            raise AssertionError(f"bad fixed-form line: {ln!r}")
    print("gen_spike_test selftest OK")


def main():
    if "--selftest" in sys.argv:
        selftest()
        return
    spec = json.load(open(sys.argv[1]))
    sys.stdout.write(generate(spec))


if __name__ == "__main__":
    main()
