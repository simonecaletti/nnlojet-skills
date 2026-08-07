---
name: probe-me-ir-structure
description: >
  Measure a matrix element's infrared structure NUMERICALLY: which colour
  dipoles a soft gluon connects, with what coefficient; whether a
  permutation is an f1<->f2 flavour symmetry; an antenna's flavour
  content; whether a crossed .map's reduced-ME argument order survived
  the crossing. Use whenever reading the colour flow out of amplitude
  source (gen*amp, colour-ordering sums, interference terms) stalls or
  would take more than a few minutes — fitting ME/redME onto an antenna
  basis over a few hundred limit points answers it in minutes with exact
  rational coefficients. First question of every S,a block. NOT a
  subtraction-term validator (that is run-spike-test): this probes the
  MATRIX ELEMENT.
---

# Probing a matrix element's IR structure

Fit, over phase-space points driven into an IR limit,

```
ME / redME  =  sum_d  c_d * X30(a, k, b)        (soft mode)
ME          =  sum_d  c_d * [X * redME_d]       (collinear mode)
```

by least squares. The coefficients c_d ARE the colour structure: which
dipoles exist, with what weight. Verified result (C1g0Z, soft gluon,
400 points, x=1e-10): `(3,5) -> 1.00000000`, `(6,7) -> 1.00000000`,
all four cross-pair candidates `0.00000000` — in minutes, after hours
of amplitude-source reading had produced no conclusion.

## Confidence — state of validation (n = 1)

- **Soft limits: high confidence, verified once** (one ME, one process).
  Expect EXACT RATIONALS. Non-rational output means the basis is
  incomplete or a reduced-ME argument order is wrong — that is itself
  the diagnostic.
- **Collinear limits: UNTESTED — expect degraded results.** The method
  needs a different basis construction (below), not just more points.
  Treat 3 stable digits as a good collinear result, not 8.
- Not a validator for a subtraction term — that is run-spike-test.

## Usage

```bash
python .claude/skills/probe-me-ir-structure/scripts/gen_probe.py spec.json > probe_gen.f
cd test/process/<PROC>
gfortran -c probe_gen.f -J obj/mod -I obj/mod            # seconds — reuse
gfortran -o probe_gen probe_gen.o $(ls obj/*.o | grep -v <checkprog>.o)  # spike-test objects
OMP_NUM_THREADS=1 ./probe_gen
```

The spec is ~15 lines (see `scripts/example_spec.json`); a hand-editable
template with identical structure is `assets/probe_template.f` for
one-offs. Both derive from the session-verified probe; the GENERATOR
itself is untested end-to-end — read its emitted Fortran once before
trusting a new spec shape.

## Boilerplate the probe must contain (each omission is a dead program)

- `common/plotmode/iplot` with `iplot != 0` — otherwise `ecuts_epem` in
  `src/core/null.f` prints "incorrect version of ecuts_epem used for
  production" and `stop`s before any output.
- The reduced ME returns **NaN** unless `call unset_map()` then
  `call set_map(N, N-1, (/antenna/), (/pmap/), ipass)` ran first — it
  reads `kin(N-1)` via `s6`/`zA6`, which `get_*7` never fills. Reduced
  momentum indices are `j1..j9` from `Mapping_mod`.
- Process init: `init_proc`, `init_map`, `setSqrts_proc`, `setScales`,
  `init_kin`, jet-parameter and flavour commons (`nf1/nf2/nfB1`, ...)
  as in the check programs — pass them via the spec's `setup_lines`.
- Per-point normalisation of the design matrix (each point scaled by
  1/sum(b^2)) so no single deep point dominates.

## Soft vs collinear — the part that actually differs

**Soft** (eikonal factorisation is exact at leading power and
spin-independent):
- all dipoles share the SAME reduced ME → divide it out, fit the ratio;
- one `set_map` serves every basis term (leading-power reduced
  kinematics are mapping-insensitive);
- basis functions have distinct poles `1/(s_ak s_kb)` → well-conditioned,
  and conditioning IMPROVES as x -> 0. Go deep (x ~ 1e-10).

**Collinear** — two independent problems, both fatal if ignored:
1. *Spin correlations.* Quark parent (q->qg): scalar factorisation holds
   directly. Gluon parent (g->qqb, g->gg): the splitting function
   carries a k_perp^mu k_perp^nu term that vanishes only after
   azimuthal averaging — use `rotp7(i,j)` (pi/2 rotation about the
   collinear axis) and average, exactly as the spike tests do.
2. *Basis degeneracy.* Every antenna containing the collinear pair has
   the SAME 1/s_ij pole with the SAME residue: as x -> 0 the design
   matrix collapses toward rank 1 and only the SUM of coefficients is
   recoverable. The inverse of the soft case. Fix: the basis element
   must be `antenna x reduced ME`, UN-divided — the antennae are
   degenerate but the reduced MEs (different colour orderings of
   B2g0Z, f1B* vs f2B*) are not, and that lifts the degeneracy. Work at
   moderate x (~1e-6), scan several x, judge by COEFFICIENT STABILITY
   across x rather than one deep evaluation.

## Solver requirements (the original probe.f got this wrong)

Never clamp a small pivot (probe.f's `|pivot|<1e-30 -> 1e-30` silently
fabricates coefficients on a rank-deficient system — and a collinear
basis IS rank-deficient). The generator's solver instead: pivots
relative to the matrix scale, REFUSES to print coefficients when
degenerate (reports which basis directions are undetermined), and
always prints the fit residual. No residual printed = do not trust the
run.

## Questions this answers (all arose in one session)

- which colour dipoles, what coefficient (the S,a block question);
- is this permutation an f1<->f2 symmetry? (measure it, don't copy);
- coefficient discrepancies between processes (the 1/2 vs 1/4 E30FF
  question) — measure both sides;
- an antenna's flavour content when source comments are silent;
- whether a crossed .map's reduced-ME argument order survived
  un-crossing (fit with both orders: the wrong one is not rational).

## Prior art / non-overlap

`test/process/epem/check_ME_*` and `*_OL_*` programs are ME-vs-OpenLoops
validation (verified: `check_ME_epemRR.f` header "test of OL and NNLOJET
matrix elements") — a different job. The original one-off
`test/process/epem/probe.f` is superseded by this skill's assets and can
be removed from the test directory.

Related: write-subtraction (uses this right after pinning the
reduced-ME set), run-spike-test (validator, downstream),
antennae-naming-convention (basis-function names).
