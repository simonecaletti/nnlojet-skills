---
name: write-runcard
description: >
  Compose, modify, or debug an NNLOJET runcard (.run file). Use whenever the
  user asks to write, create, adapt, or fix a runcard, set up a run for a
  process, add histograms/selectors/cuts/scales to a run, or when NNLOJET
  fails to parse a runcard. Also use to translate a physics request
  ("NNLO ZJ at 13 TeV with pT cuts...") into a complete runcard.
---

# Writing an NNLOJET runcard

Parser: `driver/core/IOHelper.f90`. Comments are `!`; keys are
order-independent within a block; Fortran `d0` literals accepted.

## Start from a template, not from scratch

Selection order for a base card:

1. `examples/<PROC>/<arxiv#>/*.run` — paper benchmarks (with reference
   results in sibling `*result*` dirs); `examples/<PROC>/test/*.run` —
   lighter variants.
2. `test/regression_tests/runcards/<PROC>/*.run` — small, fast skeletons.
3. `driver/example.run` — annotated grammar reference (ZJ); the SELECTORS
   syntax comments live here.
4. The hepforge examples page (https://nnlojet.hepforge.org/examples.html)
   mirrors `examples/` — only needed if a process is absent locally (site is
   behind an anti-bot wall; prefer the repo).

Legacy positional cards (first line `NNLOJET_RUNCARD`, line-numbered
values): convert with `driver/bin/patch_runcard.sh <card>` first. Its awk
reads the card by absolute line numbers, so any deviation silently corrupts
the output — always diff the result against the original (`.bak`).

## Translating the user's physics request

When the user describes observables, cuts, scales, or precision beyond just
the process name, map them to runcard sections:

- **Process + collider + energy** → `PROCESS` block (`sqrts`, `jet =
  antikt[0.4]`, `decay_type`, ...). Check the process exists:
  `ls src/process/` or `driver/process/` (do NOT trust `-listprocs`, it is
  hard-coded text, not build-introspection).
- **"histogram of X" / "distribution in X"** → `HISTOGRAMS` line. First
  verify the observable name with `./NNLOJET -listobs <PROC>`
  (build-accurate: calls init_proc + init_obs; if it dies with a
  `libLHAPDF.so` loader error, see the get-lhapdf-lib-path skill). If the
  observable does not exist, tell the user and offer the `add-observable`
  skill — do not invent names.
- **"cut on X" / fiducial cuts** → `SELECTORS` line. Selectors use the SAME
  registered observable names as histograms and scales (single registry,
  case-insensitive). Same `-listobs` check, same add-observable fallback.
- **"scale choice mT/HT/fixed"** → `SCALES` block; dynamic scales are
  observables too.
- **Order (LO/NLO/NNLO) or specific contributions** → `CHANNELS` wildcards;
  specific partonic channels → numeric IDs (see CHANNELS below).
- **Precision/statistics** → `warmup`/`production` events and iterations in
  `RUN`.
- Anything the user did not specify: keep the template's values and flag
  the ones that matter (pdf, tcut, seeds) rather than silently inheriting.

## Block reference

Seven blocks. Required: PROCESS, RUN, SELECTORS, SCALES, CHANNELS.
HISTOGRAMS required unless warmup-only. PARAMETERS optional.
Full key lists: grep `extract_option` in `driver/core/IOHelper.f90`.

### PROCESS <NAME> ... END_PROCESS

Common keys: `collider` (`LHC|PP`, `TEVATRON|PPBAR`, `DIS|EP`, `PE`,
`LEP|EPEM`, `PA|LHCASY`), `sqrts` (or `Ebeam1`+`Ebeam2`), `jet`
(`none|antikt|akt|ca|cam|kt|flavkt|flavantikt|jade|durham`, radius in
`[...]`), `jet_exclusive`, `jet_recomb` (default `v4`), `decay_type`,
`V_NC` (`ZGAMMA|GAMMA|Z`), `L_NC`; DIS: `psdis`, `dis_frame`
(`FIXT|BREIT|LAB`), `dis_eventshapes`, `eprc`.

### RUN <ID> ... END_RUN

Common keys: `pdf = NAME[member]` (also `pdf1`/`pdf2` for asymmetric
beams), `tcut`, `iseed`, `warmup = N[niter]`, `production = N[niter]`,
`iplot`, `point_check`, `pole_check`, `multi_channel`, `angular_average`,
`cache_kinematics`, `scale_coefficients`, `reset_vegas_grid`.

### PARAMETERS ... END_PARAMETERS  (optional)

`MASS[X] = value`, `WIDTH[X] = value`.

### SELECTORS ... END_SELECTORS

Verbs: `select|accept`, `reject`, `select_if_valid`, `reject_if_valid`.
Prefer single-line `select obs min = a max = b` over separate
select/reject lines. Lines combine with logical AND; wrap alternatives in
`OR ... END_OR`.

### HISTOGRAMS ... END_HISTOGRAMS

Per line: `obs nbins min max` — all three mandatory, no defaults from the
code. Options: explicit edges `obs [b0,b1,...,bn]`, `binning = log`,
`> filename` alias — MANDATORY when the same observable appears twice,
otherwise output files clash — `fac=`, `mu0=`, `grid=`, `output_type`,
`cumulant = -1|0|+1`, `run_group=`, nested `HISTOGRAM_SELECTORS ...
END_HISTOGRAM_SELECTORS`, `COMPOSITE ... END_COMPOSITE`. `cross` is the
cross-section pseudo-observable (no binning args).

### SCALES ... END_SCALES

`muf = [fac *] {value|obsname}` + `mur = ...` per line; first line is the
central scale. A product of two observables is rejected by the parser.

### CHANNELS [region = a|b|all] ... END_CHANNELS

Wildcards: `ALL LO NLO NNLO V R VV RV RR`, or numeric IDs from
`driver/process/<PROC>/selectchannel<PROC>.f`. For channel breakdowns by
initial state / matrix element, generate ID lists with
`python/gen_channels/` (`./gen-channels <PROC> -rc -IP1 u -IP2 g ...`;
must run from its own directory — cwd-relative config path).

## Validate and smoke-test

1. `cd driver && ./NNLOJET -listobs <PROC>` — every observable used in
   HISTOGRAMS, SELECTORS, HISTOGRAM_SELECTORS, and SCALES must appear.
2. Short smoke run in a scratch directory (running in `driver/` litters it
   with `.dat`/`.log`/grid files):

   ```bash
   mkdir -p /tmp/nnlojet_smoke && cd /tmp/nnlojet_smoke
   <path>/driver/NNLOJET -run card.run   # with small warmup, e.g. warmup = 1000[2]
   ```

3. Check the `.log` and that expected `.dat` files appear, named
   `<PROC>.<RUNID>.<PART>.<obs>.s<seed>.dat`.

An immediate `error while loading shared libraries: libLHAPDF.so` from
either command is an environment problem, not a runcard problem — see
the get-lhapdf-lib-path skill.

## Gotchas

- **tcut is baked into the warmup grid filename**
  (`<PROC>.<RUNID>.y1.00E-07.LO`). Changing `tcut` between warmup and
  production silently discards the grid — production restarts cold with no
  error.
- Missing `END_<BLOCK>` → "invalid runcard: can't find closing ... block!".
- "parse_hist_io: missing/invalid `bins = ...` option" → a histogram line
  lacks `nbins`/`min`/`max`.
- "only warmup run => skipping HISTOGRAMS parsing" → informational, not an
  error: HISTOGRAMS is ignored for warmup-only cards.
- Unknown observable name → print + hard stop from `getIdFromName_obs`;
  fix the name or add the observable (add-observable skill).
- `-iseed` / `-imember` CLI flags override the runcard values — useful for
  seed scans without editing the card.
