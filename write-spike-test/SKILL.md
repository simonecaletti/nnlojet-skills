---
name: write-spike-test
description: >
  Create a new NNLOJET spike test (check*.f program + makefile under
  test/process/<PROC>/) for a process or contribution that has none. Use
  when run-spike-test finds no test for the requested process/contribution,
  or when the user asks to set up spike tests for a new process, a new
  channel set, or a new contribution (R, RR, RV).
---

# Writing a new spike test

Never start from scratch: copy the nearest existing test and adapt.

- **epemZH2bb** is the modern template (local `obj/`, self-`mkdir` output
  dir, auto `-fallow-argument-mismatch`): lepton-initiated, FF limits
  only.
- **DISWp** is the reference for hadron-initiated processes with IF/II
  limits (its `check5to3.f` enumerates 66 modes: double soft, triple
  collinear FS/IS, soft+collinear, double collinear combinations) and for
  the factored layout `test5to3<PROC>.f` / `tests5to3<PROC>.f` /
  `stype5to3<PROC>.f`.

## Files

`test/process/<PROC>/`: `makefile`, `makedepend.sh` (copy), and one
program per contribution, named `check<N>to<M>[loop].f` (N = partons in
the full ME, M = partons in the born-level kinematics; `loop` = RV).

## Makefile

Pattern (epemZH2bb style):

```make
BASE = $(PWD)/../../..
-include $(BASE)/NNLOJET.mk
FC = gfortran
FFLAGS = -O0 ... -ffixed-line-length-none -J $(MODS_DIR) -I$(INCPATH)
VPATH += $(BASE)/driver/core $(BASE)/src/process/<DIR> ...
LIBFILES = $(SPIKECORE) $(CORE) $(CORETEST) $(ANTENNA) $(ANTENNAINT) \
           $(MATRIX_<PROC>) $(EMPTY_MC)
LIBFILES += auto<term1>.f
LIBFILES += auto<term2>.f
```

Key macros from `NNLOJET.mk`: `SPIKECORE` (kinematics/mapping core),
`CORETEST` (rambo generators + `writespikeRR.f`/`writespikeRV.f`),
`MATRIX_<PROC>` (the full MEs), and `EMPTY_MC = null.f null_hfrag.f90` —
stubs for `bino`/`getqcdnorm` so the SAME `auto*.f` used in the MC runs
standalone in the test (no duplicated subtraction sources). List every
subtraction `auto*.f` under test in `LIBFILES`. Use a local `obj/`, not
the main build's `$(BASE)/obj`.

## Harness facts (each rediscovery costs a build cycle)

- 7-particle limit generators live in `src/rambo/librambo.f` and mirror
  the 6-particle ones argument-for-argument: `get_ss7 get_sco7 get_ds7
  get_tc7 get_sc7 get_dc7`, plus `rotp7(i,j)` (π/2 rotation about the
  collinear axis, for azimuthal averaging).
- Init sequence: `init_proc`, `init_map`, `setSqrts_proc`,
  `setScales`, `init_kin(nPartons,10)`; cuts via
  `ecuts_epem(1,N,ipass)`-style calls; and `common/plotmode/iplot`
  MUST be nonzero — otherwise `ecuts_*` in `src/core/null.f` prints
  "incorrect version ... used for production" and `stop`s before any
  output.
- Mode counts, 5 final-state partons: 10 double soft + 10 triple
  collinear + 30 soft-collinear + 15 double collinear + 5 single soft
  + 10 single collinear = 80. GENERATE the mode blocks with a short
  script — hand-writing 80 `case` blocks invites typos (this is how
  the current `epem/check5to3.f` was produced).

## Check program structure

From `check3to2.f`/`check4to2.f` (epemZH2bb), keep this skeleton:

1. **Give it a CLI**: `./check ITYPE [MODELO MODEHI] [IPOINT] [ILOW]`
   via `getarg`, with defaults and a usage line when ITYPE is absent.
   Hard-coded `do itype=n,n` / `do mode=1,65` loops mean every
   narrowing of a test costs an edit+recompile+relink cycle — the
   argument-driven rebuild made the debug loop several times faster.
2. `iplot` flag: `2` = print per limit (default for Claude runs), `1`
   = histogram + gnuplot output (user-driven; create the output dir
   with `call system("mkdir -p "//sdir)`). Print SUMMARY STATISTICS
   per mode — `n / max|ME| / median / mean / min / max` — not raw
   event dumps; the max|ME| column is what separates genuine from
   junk modes at a glance (run-spike-test).
3. Mode loop — one `mode` per infrared limit, each defining:
   - `stitle`/`sname` (human-readable limit name — run-spike-test reports
     these, make them precise: `'"6 soft"'`, `'"5||6 collinear"'`);
   - window `xmin`, `xmax` around 1 and approach values `x(1..3)`
     (typically `1d-5, 1d-6, 1d-7`; softer, e.g. `1d-4..1d-6`, for
     delicate limits);
   - the constrained phase-space call: `get_ss<n>` for soft
     (`emij = sqrts*sqrt(1d0-x(ips))`), `get_sco<n>` for collinear
     (`emij = sqrts*x(ips)`), passing the spectator indices;
   - azimuthal rotations for g→gg / g→qqb collinear limits (see
     `doc/process/VFH/texfiles/spikesAndRotation.tex`); pure q||g limits
     don't need them.
4. Jet-function call: the process's `MYJET` from `maple/iprocess.map`
   (e.g. `call ecuts_epem_vh(0,npar,ipass)`).
5. `test(itype)` / `tests(itype)` functions: `case(itype)` pairs
   selecting full ME and subtraction with IDENTICAL momentum arguments;
   number channels sequentially and comment each case with the partonic
   channel.
6. Per point: `wt1 = test(itype)`, `wt2 = tests(itype)`,
   `rat = wt1/wt2`; print for `iplot=2`, histogram via
   `writespikeRR`/`writespikeRV` for `iplot=1`; count NaNs and
   out-of-window points.

**Mode coverage is the correctness-critical part**: enumerate EVERY
unresolved limit of the contribution — R: each soft gluon, each collinear
pair (FS and, for hadronic processes, IS); RR additionally: double soft,
each triple-collinear cluster, soft+collinear, and double-collinear
(FF+FF, IF+FF, ...) combinations; RV: same single-unresolved set as R.
Use DISWp's `check5to3.f` mode grouping as the checklist model.

## After writing

Build with `make -j8`, then hand over to **run-spike-test**. A brand-new
term also needs the **autogen-subtraction** checklist (NNLOJET.mk +
this makefile + the case entries). `hg add` the new files (no commit
unless asked).
