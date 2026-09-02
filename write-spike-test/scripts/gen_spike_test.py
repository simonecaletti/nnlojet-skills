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
  xs_coll      : x values for collinear modes    (default 1e-4..1e-6)
                 — the pass criterion is a plateau at 1 ACROSS the
                 scan, not one deep evaluation (run-spike-test).
                 The emitted calls set a cluster MASS em = sqrts*x, so a
                 collinear mode probes s_ij/s ~ x**2: the defaults stay
                 above the double-precision floor, and a deeper override
                 is warned about rather than silently producing roundoff
                 that reads exactly like a missing counterterm.
  n_azim       : orientations in the azimuthal average of a gluon-parent
                 cluster (default 4; 2 is NOT enough — see below)

The emitted program also prints, per accepted point and gated on
NNLOJET_WTDEBUG=1, its full ME as

    WTDBG MEFULL 0 0 <wt1>

which is what fit_lines.py (run-spike-test) solves against the per-line
WTDBG dumps of the subtraction term.  It is printed BEFORE the azimuthal
rotations so that it pairs with the unrotated evaluation; moving it after
them silently mis-pairs every rotated mode's constraint rows.

Azimuthal averaging.  An antenna is spin-averaged, so it reproduces a
gluon-parent collinear limit only after averaging over the cluster's
azimuth.  Two things matter and both were wrong in the first version of
this generator:
  * the WHOLE cluster must rotate about its parent axis.  Rotating two
    legs of a triple-collinear cluster averages nothing: on a 5-parton
    epem term the tc modes 3||4||7 and 4||5||6 sat at ~50% outliers
    with the 2-leg rotation and went to 1.000000 with ZERO outliers once
    all three legs moved together.
  * `n_azim` orientations kill harmonics up to cos((n_azim-1)phi).  4 is
    enough for single- and triple-collinear clusters (exact above), but
    soft+collinear modes carry higher harmonics: same term, mode
    "4 soft + 5||6" measured 40% / 18% / 10% outliers at n_azim = 4 / 8 /
    16, median -> 1.000000.  A residual that SHRINKS as n_azim grows is
    the average, not the subtraction term; one that does not is real.

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
                # azimuthal average of a gluon-parent cluster must rotate
                # the WHOLE cluster about the parent direction; rotating
                # only the first two legs leaves the correlation in.
                out.append((fam, f"{i}||{j}||{k}", (i, j, k),
                            tuple(sorted(fsset - {i, j, k})),
                            (i, j, k) if cluster_parent(
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


#: Below this invariant ratio a double-precision evaluation of the ME is
#: dominated by roundoff, so the spike ratio stops measuring the physics.
#: Calibrated on a 5-parton epem term (250 pts/mode): s_ij/s = 1e-14
#: still gives median 0.999997 with ZERO outliers on the g->qqbar
#: collinear modes, while 1e-16 gives 13/100 outliers and 1e-18 collapses
#: the median to 0.65.  So the floor sits between 1e-14 and 1e-16.
PRECISION_FLOOR = 1e-15


def warn_precision(xs_soft, xs_coll, fams):
    """Refuse to let a scan silently run past double precision.

    Collinear-driven families set em = sqrts*x, so the invariant ratio is
    x**2.  A scan whose deepest point is below ~1e-13 produces ratios that
    drift away from 1 for numerical reasons only — and the failure is
    indistinguishable by eye from a missing counterterm, so it burns a
    debugging cycle every time.
    """
    for fam, xs in (("xs_coll (sco/tc/dc)", xs_coll),
                    ("xs_soft, sc's collinear leg", xs_soft)):
        if fam.startswith("xs_soft") and "sc" not in fams:
            continue            # only sc uses xs_soft for a collinear leg
        deep = [xv for xv in xs if xv ** 2 < PRECISION_FLOOR]
        if not deep:
            continue
        print(f"WARNING: {fam}: x="
              + ", ".join(f"{xv:.0e}" for xv in deep)
              + f" probe s_ij/s down to {min(deep) ** 2:.0e}, below the "
              f"{PRECISION_FLOOR:.0e} double-precision floor — ratios "
              f"there are roundoff, not physics. A genuine mode that is "
              f"exact at the shallow x and drifts to ~0.6-0.8 only at "
              f"the deep one is this, not a missing counterterm.",
              file=sys.stderr)


def generate(spec):
    fs = {int(k): v for k, v in spec["fs_partons"].items()}
    npar = spec["npar"]
    nfs = len(fs)
    contrib = spec["contribution"]
    fams = RR_FAMS if contrib == "RR" else R_FAMS
    modes = enumerate_modes(fs, fams)
    cls = {(m["family"], m["name"]): m
           for m in classify(fs, spec["born"], fams)}
    # default x scans: the typical working range 1e-7..1e-9 (deeper —
    # down to 1e-10 — can work but is not automatically better; the
    # object is a PLATEAU at 1 across the scan, and instabilities at
    # the deepest x are expected and reportable, not silently
    # droppable). Override per family with xs_soft / xs_coll.
    # `sc` sets its COLLINEAR leg from xs_soft (em2 = sqrts*x), so the
    # soft default has to respect the same floor as xs_coll; ss/ds alone
    # would tolerate 1e-9 (their scale is linear in x, not quadratic).
    xs_soft = spec.get("xs_soft", [1e-5, 1e-6, 1e-7])
    # COLLINEAR default is NOT the same as the soft one.  The generated
    # calls set the cluster mass as em = sqrts*x, so the invariant ratio
    # probed is s_ij/s ~ x**2: x=1e-8 already means 1e-16, i.e. machine
    # epsilon, and the ratio degrades for pure roundoff while looking
    # exactly like a subtraction error (a genuine mode drifting to 0.6-0.8
    # only at the deepest x, exact at the shallowest).  Note this differs
    # from the hand-written legacy check programs, which use
    # em = sqrts*sqrt(x) and therefore probe s_ij/s ~ x.
    xs_coll = spec.get("xs_coll", [1e-4, 1e-5, 1e-6])
    warn_precision(xs_soft, xs_coll, fams)
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
    L.append("      dimension mrot1(4), mrot2(4)")
    L.append(f"      parameter (n_azim={spec.get('n_azim', 4)})")
    L.append("      parameter (pi=3.141592653589793238d0)")
    L.append("      parameter (dazim=2d0*pi/dble(n_azim))")
    L.append("      character*64 stitle")
    L.append("      character*32 arg")
    L.append("      character*8 wtdbgenv")
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
    L.append("      call get_environment_variable('NNLOJET_WTDEBUG',"
             "wtdbgenv)")
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
    L.append("        nrot1=0")
    L.append("        nrot2=0")
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
        for nr, rot in ((1, rot1), (2, rot2)):
            if not rot:
                continue
            L.append(f"          nrot{nr}={len(rot)}")
            for p, leg in enumerate(rot, 1):
                L.append(f"          mrot{nr}({p})={leg}")
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
    # per-point ME for fit_lines.py (run-spike-test).  Printed here, BEFORE
    # the azimuthal rotations below, so it pairs with the UNROTATED
    # evaluation -- fit_lines relies on that ordering.  Runtime-gated, so a
    # normal run is byte-identical.
    L.append("            if (wtdbgenv.eq.'1') write(*,'(A,E23.15)')")
    L.append("     .        'WTDBG MEFULL 0 0 ', wt1")
    L.append("c           azimuthal average for gluon-parent clusters:")
    L.append("c           4 orientations (0, pi/2, pi, 3pi/2) about the")
    L.append("c           CLUSTER axis, rotating every leg of the cluster.")
    L.append("c           A 2-point average kills cos(2phi) only; the")
    L.append("c           triple-collinear splitting functions carry")
    L.append("c           higher harmonics, and rotating a subset of the")
    L.append("c           cluster does not average anything at all.")
    L.append("c           The 4th call restores the original orientation,")
    L.append("c           so cluster 2's average starts from the same point.")
    for nr in (1, 2):
        L.append(f"            if (nrot{nr}.gt.0) then")
        L.append("              s1=wt1")
        L.append("              s2=wt2")
        L.append(f"              do irot=1,n_azim-1")
        L.append(f"                call rotcl{npar}(mrot{nr},nrot{nr},dazim)")
        L.append("                s1=s1+test(itype)")
        L.append("                s2=s2+tests(itype)")
        L.append("              end do")
        L.append(f"              call rotcl{npar}(mrot{nr},nrot{nr},dazim)")
        L.append("              wt1=s1/dble(n_azim)")
        L.append("              wt2=s2/dble(n_azim)")
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

    L.append(f"""
************************************************************************
*     rotate the n momenta listed in idx by pi/2 about their COMMON
*     axis (the cluster's parent direction).  librambo's rotp{npar} does
*     this for exactly two legs; a collinear cluster of three needs all
*     three moved together or the azimuthal average is not an average.
************************************************************************
      subroutine rotcl{npar}(idx,n,theta)
      use KinData_mod
      implicit real*8(a-h,o-z)
      dimension idx(n), ax(3), pd(3,{npar})
      ct=cos(theta)
      st=sin(theta)
      omct=1d0-ct
      ax(1)=0d0
      ax(2)=0d0
      ax(3)=0d0
      do m=1,n
        do ii=1,3
          ax(ii)=ax(ii)+p{npar}(ii,idx(m))
        end do
      end do
      aa=sqrt(ax(1)**2+ax(2)**2+ax(3)**2)
      if (aa.le.0d0) return
      ux=ax(1)/aa
      uy=ax(2)/aa
      uz=ax(3)/aa
      do m=1,n
        j=idx(m)
        pd(1,m)=(ct+ux*ux*omct)*p{npar}(1,j)
     .         +(ux*uy*omct-uz*st)*p{npar}(2,j)
     .         +(ux*uz*omct+uy*st)*p{npar}(3,j)
        pd(2,m)=(ux*uy*omct+uz*st)*p{npar}(1,j)
     .         +(ct+uy*uy*omct)*p{npar}(2,j)
     .         +(uy*uz*omct-ux*st)*p{npar}(3,j)
        pd(3,m)=(ux*uz*omct-uy*st)*p{npar}(1,j)
     .         +(uy*uz*omct+ux*st)*p{npar}(2,j)
     .         +(ct+uz*uz*omct)*p{npar}(3,j)
      end do
      do m=1,n
        j=idx(m)
        do ii=1,3
          p{npar}(ii,j)=pd(ii,m)
        end do
      end do
      call fillS_kin({npar})
      return
      end""")
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
