---
name: run-spike-test
description: >
  Build, run, and interpret NNLOJET spike tests (test/process/<PROC>/check*.f)
  that validate antenna-subtraction terms against the full matrix element in
  infrared limits. Use whenever the user asks to spike-test a channel, check
  a subtraction term, verify infrared cancellation, or after a subtraction
  term was edited/registered (write-subtraction / autogen-subtraction
  skills). Runs with iplot=2 only; iplot=1 plot production is user-driven.
---

# Running spike tests

A spike test evaluates, for phase-space points driven into an infrared
limit, the ratio R = fullME / subtraction. Pass criterion: **the ratio
approaches 1 asymptotically as the infrared parameter x decreases
(x(1) > x(2) > x(3)), for many points in the phase space, for every
GENUINE limit of the channel** — the mode lists also probe
configurations that are not limits at all (see classification below).

## Locate the test

- `test/process/<PROC>/` — note the mapping to maple-level processes is
  not 1:1: maple `epemZH` → test `epemZH2bb` / `epemZH2gg` / `epemZH2gaga`.
- Executable naming tells the contribution: `check<N>to<M>.f` with small
  N−M gap = R (NLO), `check<N>to<M>loop.f` = RV, larger N−M gap = RR
  (e.g. epemZH2bb: `check3to2` = R, `check3to2loop` = RV, `check4to2` =
  RR; DISWp: `check4to3`, `check4to3loop`, `check5to3`). Confirm from the
  header comment.
- The channel numbering is process-dependent: read the `case(itype)`
  blocks in the `test`/`tests` functions (in `check*.f` itself or in
  factored `test*<PROC>.f` / `tests*<PROC>.f` files) to see which
  (ME, subtraction) pair each CHANNEL integer selects.
- The limits ("modes") are hard-coded in the source: each `mode` sets a
  name (`'6 soft'`, `'5||6 collinear'`, ...), the window `xmin/xmax`, the
  approach values `x(1..3)`, and the phase-space generator call
  (`get_ss*` soft, `get_sco*` collinear).
- No test exists for this process/contribution → **write-spike-test**
  skill.
- VV/U has no unresolved radiation and CANNOT be spike-tested — its
  numerical validation is the pole check (**run-pole-check** skill),
  which also covers the ε-poles of V and RV.

## Build

```bash
cd test/process/<PROC>
make check4to2 -j8
```

Traps:
- Some older makefiles (DISWp style) share `OBJDIR = $(BASE)/obj` with the
  main build — stale-`.mod` conflicts if the driver was built with other
  flags, and `make clean` there wipes the main build's objects. Newer ones
  (epemZH style) use a local `obj/`.
- gfortran ≥ 10 needs `-fallow-argument-mismatch`; newer makefiles add it
  automatically, older ones don't — add to FFLAGS if compilation fails on
  argument mismatch.
- Some makefiles (DISWp) require `lhapdf-config` on PATH.
- Parallel `make -jN` can die with "Fatal Error: Cannot open module file
  'dis_mod.mod'" — a module-dependency race, not a real error. Re-run make
  (may take two passes) or use `-j1`.

## Run — iplot=2 only

**Determine the argument convention FIRST — it is not uniform.** Of the
159 check*.f programs: 113 take NO argument (the channel is the
hard-coded `do itype=<n>,<n>` loop — selecting a channel means editing
that line and REBUILDING; e.g. `epem/check5to3.f` ~line 79); 36 take one
argument, CHANNEL (e.g. `DIS/check5to3.f`, `epemZH2bb/check4to2.f`);
5 take two, IORDER CHANNEL (the loop/RV programs, e.g.
`epem/check4to3loop.f`); 5 are special-purpose, not channel-driven
(`GGJfc/check_njet.f`, `Z/check4to2.f`, ...). The convention varies
WITHIN a process dir (`epem/check5to3.f` none, `epem/check4to3loop.f`
two) — it cannot be inferred from the process or the filename. Detect
it: run the binary bare — every argument-taking program prints a usage
line and stops, harmlessly; or grep the source for `getarg` (no hits =
hard-coded; then grep `do itype=` for the line to edit).

`iplot` is hard-coded in the check source. Ensure it is set to `2`
(prints a few points per limit to stdout) before building:

```fortran
      iplot = 2           ! print events in each limit
```

Then (one-argument programs; adapt per the convention above):

```bash
OMP_NUM_THREADS=1 ./check4to2 <CHANNEL>
```

**iplot=2 defaults often cannot show the pass criterion.** E.g.
`epem/check5to3.f` sets `ipoint=5, ilow=3, iup=3` for iplot=2 — only
the DEEPEST x, 5 points: no trend visible. Cheap fix, still entirely
within iplot=2: temporarily set `ilow=1` (prints all three x values —
also needed for the scaling-exponent classification below) and raise
`ipoint` (~100; 100 points × 3 x values ran in ~3 s), rebuild, run —
then RESTORE both values. This is NOT iplot=1 and carries none of its
cost. Also check the MODE-loop bound: e.g. `check5to3.f` has
`do mode=1,65`, silently skipping modes 66–80 — ALL the single-soft and
single-collinear limits; raise it to cover them.

Do NOT run with `iplot=1` on your own initiative: it histograms ~1000
accepted points per (mode, x) over all modes — hours per channel — and
its output is gnuplot files meant for human inspection. If the user
explicitly asks for the plots, set `iplot=1`, rebuild, launch it for
them, and let the user run gnuplot and read the plots (older programs
also need the output dir created first, e.g. `mkdir 4to2` — some create
it themselves, some don't).

## Interpret

For each mode, stdout shows per point: `wt1` (full ME), `wt2`
(subtraction), `rat`. Check across the three x values:

- **Pass**: `rat` → 1 for (essentially) all printed points, and the
  spread around 1 shrinks as x decreases. Must hold for every mode.
- **Fail — report precisely**: state the CHANNEL and the failing mode(s)
  by name (e.g. "channel 5, limit 5||6 collinear"). That is the actionable
  output: the failing limit identifies which subtraction lines are
  suspect (limit → antenna mapping: antennae-naming-convention skill).
- Symptom guide (apply ONLY to modes classified GENUINE):
  - ratio → constant ≠ 1 in a limit → wrong normalisation / colour factor
    on the covering line;
  - ratio diverges or → 0 in a GENUINE limit → missing or wrong antenna
    for that limit, or wrong mapped arguments (in a dead mode, rat → 0
    is meaningless — see classification);
  - ratio oscillates around 1 without narrowing in a g→gg or g→qqb
    collinear limit → azimuthal-rotation issue (see
    `doc/process/VFH/texfiles/spikesAndRotation.tex`), not necessarily a
    wrong .map;
  - NaN counters → broken momentum mapping.

### Classify every mode BEFORE reading ratios

(Empirically validated on 240 channel×mode combinations, epem
C1g0/Ct1g0/B3g0 — 240/240 with the rule below.)

1. **Decode the channel species**: match the `test(itype)` call's
   argument order against the ME's particle content in
   `driver/maple/<PROC>.map` (e.g. C1g0Zepem(i4,i5,i6,i7,i3,...) vs
   `[qb,g,Q,Qb,q,...]` → 5=gluon, 3,4 and 6,7 = qq̄ pairs).
2. **Reduced-Born rule** — a mode is GENUINE iff, after (a) replacing
   every collinear cluster by its parent parton, where a cluster is
   valid iff its NET FLAVOUR is a single parton (q∥g, g∥g, same-flavour
   q∥q̄ pass; cross-flavour q∥Q fails), and (b) deleting every soft
   parton — gluons, or a same-flavour q q̄ PAIR going soft together —
   the remaining state is a legal Born for the process. Two traps this
   rule fixes:
   - composites are NOT decomposed pairwise: {q,q̄,Q} is a genuine
     triple-collinear cluster (parent Q) although its q∥Q pair is dead;
   - all-components-genuine is NOT sufficient: in B3g0 (q̄ggg q), q̄∥q
     is species-valid but collapsing it leaves ggg — no Born → dead;
     likewise a soft same-flavour q q̄ pair is genuine ONLY if deleting
     it leaves a Born (C1g0 "3,4" genuine, B3g0 "3,7" not).
   - Do NOT use colour adjacency in the argument list: check programs
     sum over colour orderings, so non-adjacent pairs are fully
     connected (B3g0 3/5 is genuine).
3. **Verdict per mode from the scaling exponent, not |wt1| or the
   median.** With `ilow=1` (three x values) fit
   p = −Δlog10|wt1|/Δlog10 x: it comes out integer-quantized. Genuine
   modes reach the family maximum (double soft 4, soft-collinear 3,
   triple/double collinear 2, single soft 2, single collinear 1);
   p ≤ p_genuine−1 = sub-singular or pointless, with no channel
   dependence. Quick triage: SORT the modes of each family by max|ME|
   — genuine and junk separate at a glance (gaps of 3–15 decades
   observed). Within a family only: never compare |wt1| across
   families, their x-sets differ (cross-family ranges overlap by ~9
   decades).
4. **median(rat)→1 is NOT sufficient**: sub-singular modes routinely
   read 1.0000 with tiny spread because a lower, correctly-subtracted
   limit dominates — a false pass that tests nothing. Conversely
   rat→0 or O(10) noise on a DEAD mode is normal (the subtraction's
   mapping just doesn't apply there), not a failure.
5. **Use the MEDIAN across points, not the max deviation**, for the
   genuine modes: azimuthal-correlation terms leave a few-% point
   spread at every x while the median sits at 1.0000 — that is a pass;
   single outlier points (even −0.04 among 100 points at median
   1.0000) do not fail a mode.
6. **Run a sibling channel whose .map is untouched as a control**
   (epem: Ct1g0 next to C1g0). Artefacts common to both are
   harness/kinematics, not your term.

## On failure

Report the failing channel + limit(s). If the user wants a fix, invoke
the **write-subtraction** skill (the failing limit pinpoints the lines),
then **autogen-subtraction**, then rerun here. Iterate until all limits
pass.
