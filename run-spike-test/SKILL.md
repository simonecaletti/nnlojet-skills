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
(x(1) > x(2) > x(3)), for many points in the phase space, and for ALL
limits of the channel.**

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

## Run — iplot=2 only

`iplot` is hard-coded in the check source. Ensure it is set to `2`
(prints ~10 points per limit to stdout) before building:

```fortran
      iplot = 2           ! print 10 events in each limit
```

Then:

```bash
OMP_NUM_THREADS=1 ./check4to2 <CHANNEL>
```

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
  suspect.
- Symptom guide:
  - ratio → constant ≠ 1 in a limit → wrong normalisation / colour factor
    on the covering line;
  - ratio diverges or → 0 in a limit → missing or wrong antenna for that
    limit, or wrong mapped arguments;
  - ratio oscillates around 1 without narrowing in a g→gg or g→qqb
    collinear limit → azimuthal-rotation issue (see
    `doc/process/VFH/texfiles/spikesAndRotation.tex`), not necessarily a
    wrong .map;
  - NaN counters → broken momentum mapping.

## On failure

Report the failing channel + limit(s). If the user wants a fix, invoke
the **write-subtraction** skill (the failing limit pinpoints the lines),
then **autogen-subtraction**, then rerun here. Iterate until all limits
pass.
