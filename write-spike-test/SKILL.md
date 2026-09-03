---
name: write-spike-test
description: >
  Create a new NNLOJET spike test (check*.f program + makefile under
  test/process/<PROC>/) for a process or contribution that has none, or
  classify which of its modes are genuine limits. Use when run-spike-test
  finds no test for the requested process/contribution, when the user asks
  to set up spike tests for a new process, a new channel set, or a new
  contribution (R, RR, RV), or to decide which modes of a channel are
  genuine vs dead (reduced-Born rule).
---

# Writing a new spike test

**Generate, don't transcribe.** The mode table for n final-state
partons is a mechanical function of n, and which modes are genuine is a
mechanical function of the flavour content (reduced-Born rule).
Hand-writing 80 case blocks — or hand-classifying 80 modes — is pure
transcription and an easy place to introduce a silent indexing error.

```bash
# which modes are genuine limits, which are dead (reduced-Born rule):
python .claude/skills/write-spike-test/scripts/genuine_modes.py spec.json

# the complete check program (mode table, azimuthal averaging on
# gluon-parent collinear modes, per-mode summary statistics, CLI):
python .claude/skills/write-spike-test/scripts/gen_spike_test.py spec.json > check5to3.f

# both have structural self-tests that encode no physics answer:
python .../genuine_modes.py --selftest && python .../gen_spike_test.py --selftest
```

Spec inputs (documented in each script's docstring): the process init
block copied from any existing check program of the process family, the
final-state partons with flavours (`g`, `q<tag>`, `qb<tag>`), the legal
Born flavour sets, optional colour-structure flags, and the channel list —
(ME, subtraction) entry-point pairs with IDENTICAL momentum arguments.
Set `colour: subleading` for `Bt`/`Bty` matrix elements and antennae, and
`interference: true` for `D`/`Dy`; the defaults are leading colour and no
interference, preserving existing specs. `genuine_modes.py` is also the
per-mode checklist the whole validation hangs on: run it before reading
any spike-test output (run-spike-test applies the same rule).

**Second consumer**: `predict_blocks.py` (write-subtraction, rung 0)
reads `genuine_modes.py --json` output to decide which legs admit an
unresolved limit — without it, it falls back to a rule that over-counts.
It reads a SUPERSET of this spec, so one file can drive both, but note
the two use different keys: this script reads `partons` (final state
only), `born`, and optional `colour`/`interference`; `predict_blocks.py`
reads `chain` (the full colour
ordering, initial legs included) and `flavours`. The flavour information
therefore appears twice and is NOT cross-validated — keep them in sync
by hand, or keep two files.

```bash
python .claude/skills/write-spike-test/scripts/genuine_modes.py spec.json --json > modes.json
python .claude/skills/write-subtraction/scripts/predict_blocks.py spec.json --modes modes.json
```

What the generated program does for you:
- **CLI** `./check ITYPE [MODELO MODEHI] [IPOINT] [ILOW]` with a usage
  line — no edit+recompile cycle to narrow a run;
- mode table built programmatically for the contribution (R/RV: single
  soft + single collinear; RR additionally: double soft, triple
  collinear, soft-collinear, double collinear), each mode titled with
  its GENUINE/dead classification;
- azimuthal averaging (`rotp<n>` pi/2 rotation) on collinear modes
  whose cluster parent is a gluon (see
  `doc/process/VFH/texfiles/spikesAndRotation.tex` for the delicate
  cases);
- per (mode, x) a **summary line** — `n / nan / max|ME| / median /
  min / max / outlier count` — not raw event dumps; the median is what
  makes a 240-line run readable at a glance, and the max|ME| column
  separates genuine from junk modes (run-spike-test).

Copying the nearest existing test remains the fallback for layouts the
generator does not cover (factored `test*<PROC>.f` files, exotic
initial states): epemZH2bb is the modern lepton-initiated template,
DISWp the reference for hadron-initiated processes with IF/II limits.

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
the main build's `$(BASE)/obj`. The generated program needs
`-ffixed-line-length-none` (already in the pattern above).

## Harness facts (each rediscovery costs a build cycle)

- Limit generators live in `src/rambo/librambo.f`, one family per
  multiplicity with parallel argument conventions (unresolved indices,
  then the driven invariant mass, then spectators): `get_ss<n>
  get_sco<n> get_ds<n> get_tc<n> get_sc<n> get_dc<n>` for n = 6,7,8
  (n = 5: `get_ss5`/`get_sco5` only), plus `rotp<n>(i,j)` (pi/2
  rotation about the collinear axis, for azimuthal averaging). Mass
  conventions: soft-type drivers take `em = sqrts*sqrt(1-x)` (mass of
  the recoil system), collinear-type `em = sqrts*x`.
- Init sequence: `init_proc`, `init_map`, `setSqrts_proc`,
  `setScales`, `init_kin(nPartons,10)`; cuts via
  `ecuts_epem(1,N,ipass)`-style calls; and `common/plotmode/iplot`
  MUST be nonzero — otherwise `ecuts_*` in `src/core/null.f` prints
  "incorrect version ... used for production" and `stop`s before any
  output.
- Mode counts are pure combinatorics (the generator's self-test checks
  them): for n final-state partons at RR, C(n,2) double soft + C(n,3)
  triple collinear + n·C(n−1,2) soft-collinear + 3·C(n,4) double
  collinear + n single soft + C(n,2) single collinear — n=5 gives
  10+10+30+15+5+10 = 80.

## x scan defaults and the chat report

The generated program defaults to x = `1e-7, 1e-8, 1e-9` per family
(override with `xs_soft`/`xs_coll`; down to `1e-10` can work but is
not automatically better). The pass criterion is a **plateau at 1
across the scan**, never a single deep evaluation; instability at the
deepest x is expected and must be REPORTED (naming that x), not
silently dropped — see run-spike-test. The per-mode summary lines the
program prints (median/min/max/outliers per x) are exactly the raw
material for the **limit table that must be reported in chat** after
every run on a new or edited term: one row per genuine limit,
including passes, ordered by family; dead modes in one line
(run-spike-test documents the table format).

## After writing

Build with `make -j8`, then hand over to **run-spike-test**. A brand-new
term also needs the **autogen-subtraction** checklist (NNLOJET.mk +
this makefile + the case entries). `hg add` the new files (no commit
unless asked).
