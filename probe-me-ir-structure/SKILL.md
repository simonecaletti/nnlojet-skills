---
name: probe-me-ir-structure
description: >
  Measure the infrared structure of a MATRIX ELEMENT or an ANTENNA
  FUNCTION numerically. Three modes: (1) dipole/colour fits — which
  colour dipoles a soft gluon connects, with what coefficient; flavour
  symmetries; whether a crossed .map's reduced-ME argument order
  survived the crossing; (2) POLE SCAN of an antenna — which invariants
  it is singular in and with what power (its pole graph), which
  argument slots are radiators vs unresolved, which dipole its soft
  limit sits on; (3) RESIDUE FITTER — what a higher antenna (X40)
  reduces to on each of its boundaries, fitted onto a basis of
  lower-antenna × reduced-ME products (the input for iterated
  counterterm construction). Use whenever reading colour flow or an
  antenna's convention out of source/papers stalls, and ALWAYS before
  using an antenna you have not used before. NOT a subtraction-term
  validator (that is run-spike-test).
---

# Probing the IR structure of a matrix element or antenna function

One shared machinery — phase-space drivers into IR limits, a
least-squares fit with degeneracy refusal — behind three questions:

1. **Dipole fit** (`scripts/gen_probe.py`): fit `ME/redME` (soft) or
   `ME` (collinear) onto an antenna basis → exact rational
   coefficients = the colour structure.
2. **Pole scan** (`scripts/antenna_probe.py`): drive an antenna into
   EVERY limit of the phase space → its pole graph.
3. **Residue fit** (`scripts/gen_probe.py`, `target` = an antenna):
   fit an X40's residue in a limit onto lower-antenna(×reduced-ME)
   products → what it reduces to on that boundary.

## The measure-before-you-use rule

**Before using an antenna you have not used before, run the pole
scan.** Do not infer its convention from the letter, the paper
equation number, or a sibling `.map`. The letter fixes the SPECIES
(antennae-naming-convention); which slots are hard radiators, which
are unresolved, and which dipole the soft limit actually sits on are a
separate CONVENTION that most `src/X40/*.f` headers do not state — and
the two are independent. A `Full` composite additionally divides its
singularities between sub-antennae in a way the composite's own name
does not reveal: scan each half too (the halves cover DIFFERENT
limits; a half with fewer poles is not "safer" — it is incomplete).

## Mode 2 — pole scan (antenna_probe.py)

```bash
python .claude/skills/probe-me-ir-structure/scripts/antenna_probe.py spec.json > polescan_gen.f
# compile/link exactly like a probe (below); run with OMP_NUM_THREADS=1
```

Spec: process init block (as in the check programs), `npar`,
`fs_partons`, and `target` = the antenna call with its kinematics
index, e.g. `"FullB40(3,4,5,6,7)"`. The program visits every single
soft, 2-parton collinear, double soft and 3-parton collinear limit and
prints per (limit, x) the **median** of `|s^p × antenna|` for
p = 0, 1, 2 (median, not mean — single deep points dominate a mean).

Reading it — the pole graph:
- **Collinear modes**: pole power = smallest p with median`|s^p A|`
  stable as x decreases (independent of the harness's x convention).
- **Soft modes**: read the printed exponent
  `alpha = −dlog10(median|A|)/dlog10(x)`.
- Radiators/unresolved/colour chain follow mechanically: the interior
  of the chain is what the antenna is singular in; endpoints are the
  radiators. The program prints a per-mode consistency check (slope
  agreement across consecutive x pairs) — trust nothing flagged
  INCONSISTENT.

## Mode 1 — dipole fit (gen_probe.py)

Fit, over phase-space points driven into an IR limit,

```
ME / redME  =  sum_d  c_d * X30(a, k, b)        (soft mode)
ME          =  sum_d  c_d * [X * redME_d]       (collinear mode)
```

by least squares. The coefficients c_d ARE the colour structure.
Verified on a 5-parton RR channel (soft gluon, 400 points, x=1e-10):
the genuine dipoles came out at exactly `1.00000000` and every
spurious candidate at `0.00000000` — in minutes, after hours of
amplitude-source reading had produced no conclusion.

## Mode 3 — residue fit (gen_probe.py, antenna target)

Set `target` to the higher antenna (any Fortran expression is legal —
`full_me` is a legacy synonym), the limit driver to the boundary of
interest, and the basis to candidate `lower-antenna × reduced-ME`
products (basis entries may omit `map`/`redme` for pure
antenna-on-antenna fits). The fitted coefficients tell you what the
X40 reduces to there — precisely the input needed to construct the
matching iterated counterterm (write-subtraction), and the measured
check that an SS combination reproduces a soft-eikonal difference
before it is written into a `.map`.

**Degeneracy caveat**: symmetry-related orientations of the same
antenna are degenerate in a soft limit. Start from a 2–3 element basis
and grow it, rather than listing every candidate and getting a
rank-deficiency refusal.

## Confidence — state of validation

- **Soft limits: high confidence.** Expect EXACT RATIONALS.
  Non-rational output means the basis is incomplete or a reduced-ME
  argument order is wrong — that is itself the diagnostic.
- **X40-level fits in double-unresolved limits: VALIDATED.** They
  behave exactly like the soft-gluon case — exact rational
  coefficients, residuals at the 1e-14 level.
- **Collinear limits**: need the basis construction below (basis =
  antenna × reduced ME, un-divided); judge by coefficient stability
  across an x scan, treat 3 stable digits as a good result.
- Not a validator for a subtraction term — that is run-spike-test.

## Usage

```bash
python .claude/skills/probe-me-ir-structure/scripts/gen_probe.py spec.json > probe_gen.f
cd test/process/<PROC>
gfortran -c probe_gen.f -J obj/mod -I obj/mod            # seconds — reuse
gfortran -o probe_gen probe_gen.o $(ls obj/*.o | grep -v <checkprog>.o)  # spike-test objects
OMP_NUM_THREADS=1 ./probe_gen
```

The spec is ~15 lines (see `scripts/example_spec.json`); a
hand-editable template with identical structure is
`assets/probe_template.f` for one-offs. Both scripts have a
`--selftest` that checks emitted-code structure WITHOUT encoding any
physics answer (`python gen_probe.py --selftest`); run it after any
edit to the generators.

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
  (A bare antenna pole scan needs no jet cuts; `cuts_call` is optional
  there.)
- Per-point normalisation of the design matrix (each point scaled by
  1/sum(b^2)) so no single deep point dominates.

## Soft vs collinear — the part that actually differs

**Soft** (eikonal factorisation is exact at leading power and
spin-independent):
- all dipoles share the SAME reduced ME → divide it out, fit the ratio;
- one `set_map` serves every basis term (leading-power reduced
  kinematics are mapping-insensitive);
- basis functions have distinct poles `1/(s_ak s_kb)` → well-conditioned,
  and conditioning IMPROVES as x -> 0. Go deep (x ~ 1e-10). But mind
  the degeneracy caveat above when the basis contains symmetry-related
  orientations of one antenna.

**Collinear** — two independent problems, both fatal if ignored:
1. *Spin correlations.* Quark parent (q->qg): scalar factorisation holds
   directly. Gluon parent (g->qqb, g->gg): the splitting function
   carries a k_perp^mu k_perp^nu term that vanishes only after
   azimuthal averaging — use `rotp<n>(i,j)` (pi/2 rotation about the
   collinear axis) and average, exactly as the spike tests do
   (antenna_probe.py does this by default in collinear modes).
2. *Basis degeneracy.* Every antenna containing the collinear pair has
   the SAME 1/s_ij pole with the SAME residue: as x -> 0 the design
   matrix collapses toward rank 1 and only the SUM of coefficients is
   recoverable. The inverse of the soft case. Fix: the basis element
   must be `antenna x reduced ME`, UN-divided — the antennae are
   degenerate but the reduced MEs (different colour orderings, f1* vs
   f2* sectors) are not, and that lifts the degeneracy. Work at
   moderate x (~1e-6), scan several x, judge by COEFFICIENT STABILITY
   across x rather than one deep evaluation.

## Solver requirements

Never clamp a small pivot (a `|pivot|<1e-30 -> 1e-30` clamp silently
fabricates coefficients on a rank-deficient system — and a collinear
basis IS rank-deficient). The generator's solver instead: pivots
relative to the matrix scale, REFUSES to print coefficients when
degenerate (reports that only coefficient sums are determined), and
always prints the fit residual. No residual printed = do not trust the
run. These two behaviours — residual report and rank-deficiency
refusal — are exactly what the `--selftest` protects.

## Questions this answers

- which colour dipoles, what coefficient (the S,a block question);
- which invariants an unfamiliar antenna is singular in, with what
  power; which slots radiate (every new-antenna placement question);
- what an X40 reduces to on each boundary (the S,b2/S,c construction
  input — see write-subtraction);
- whether the eikonal a term produces in a radiator-unresolved limit
  lands on a dipole the full ME actually has (the argument-alignment
  cross-check of write-subtraction);
- is this permutation an f1<->f2 symmetry? (measure it, don't copy);
- coefficient discrepancies between processes — measure both sides;
- an antenna's flavour content when source comments are silent;
- whether a crossed .map's reduced-ME argument order survived
  un-crossing (fit with both orders: the wrong one is not rational).

## Prior art / non-overlap

`test/process/epem/check_ME_*` and `*_OL_*` programs are ME-vs-OpenLoops
validation — a different job. Any leftover one-off `probe.f` in a test
directory is superseded by this skill's generators.

Related: write-subtraction (pole graphs and residue fits feed the
argument-alignment and counterterm-construction procedures there),
run-spike-test (validator, downstream; its mode→block diagnosis reads
the pole graphs produced here), antennae-naming-convention
(basis-function names; species-vs-convention distinction).
