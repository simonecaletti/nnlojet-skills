---
name: probe-me-ir-structure
description: >
  Measure the infrared structure of a MATRIX ELEMENT or an ANTENNA
  FUNCTION numerically. Four modes: (1) dipole/colour fits — which
  colour dipoles a soft gluon connects, with what coefficient; flavour
  symmetries; whether a crossed .map's reduced-ME argument order
  survived the crossing; (2) POLE SCAN of an antenna — which invariants
  it is singular in and with what power (its pole graph), which
  argument slots are radiators vs unresolved, which dipole its soft
  limit sits on; (3) RESIDUE FITTER — what a higher antenna (X40)
  reduces to on each of its boundaries, fitted onto a basis of
  lower-antenna × reduced-ME products (the input for iterated
  counterterm construction); (4) Z-PROFILE — the collinear counterpart
  of (3) for a bare-antenna target, where a least-squares fit is
  structurally rank-deficient. Use whenever reading colour flow or an
  antenna's convention out of source/papers stalls, and ALWAYS before
  using an antenna you have not used before. NOT a subtraction-term
  validator (that is run-spike-test).
---

# Probing the IR structure of a matrix element or antenna function

One shared machinery — phase-space drivers into IR limits, a
least-squares fit with degeneracy refusal, a z-binned profiler —
behind four questions,
each feeding a specific construction in write-subtraction:

1. **Dipole fit** (`scripts/gen_probe.py`): fit `ME/redME` (soft) or
   `ME` (collinear) onto an antenna basis → exact rational
   coefficients = the colour structure. **Feeds: the S,a block's
   dipole content.**
2. **Pole scan** (`scripts/antenna_probe.py`): drive an antenna into
   EVERY limit of the phase space → its pole graph. **Feeds: the
   antenna's slot/radiator convention** (after `antenna_slots.py`
   settled the cheaper slot-plumbing question —
   antennae-naming-convention).
3. **Residue fit** (`scripts/gen_probe.py`, `target` = an antenna;
   runnable example: `assets/example_spec_residue.json`): fit an
   X40's residue in a limit onto lower-antenna(×reduced-ME) products
   → what it reduces to on that boundary. **Feeds: the S,b2 iterated
   counterterms and the S,c soft differences** — including deciding
   between the two writings of an iterated counterterm (fit both
   orderings; the wrong one is not rational).
4. **z-profile** (`scripts/residue_profile.py`): the collinear
   counterpart of (3), for a BARE ANTENNA target where (3) is
   rank-deficient by construction. **Feeds: the same S,b2/S,c
   constructions, and specifically the decision "wrong coefficient vs
   wrong object vs needs a sub-antenna split"**, which a fitted number
   cannot express.

All three name their limits in ONE mode vocabulary
(`scripts/irlimits.py` — the same `ss/sco/ds/tc/sc/dc` + indices that
`genuine_modes.py` in write-spike-test uses): a spec says
`"limit": {"family": "sc", "unresolved": [ISOFT,IA,IB], ...}` instead
of hand-writing the generator call; verbatim `limit_call` remains the
escape hatch.

## Arrive with an expected answer

**Predict first, then measure.** Before any fit, run rung 0 of
write-subtraction:

```bash
python .claude/skills/write-subtraction/scripts/predict_blocks.py spec.json
```

It costs seconds and tells you
which block the line belongs to, its sign, and the rational family its
coefficient must live in. A fit with an expected answer is a hypothesis
test; a fit without one silently absorbs structural errors as fitted
numbers.

**The coefficient family.** Antennae are normalised colour-ordered MEs
of a parent process, so fitted coefficients come from a small discrete
set — observed across `maple/form/common/J21.map` and `J22.map`:
`{2, 1, 2/3, 1/2, 1/3, 1/4, 1/9}`. It is **not** all unit fractions
(`2*calC40FF`, `-2/3*calG30FF*calF30FF`). Derivation and the per-antenna
table: antennae-naming-convention. Read the three outcomes as:

| fit returns | meaning | action |
|---|---|---|
| the predicted rational | confirmed | proceed |
| a different family member (`1/2` for `1`) | symmetry factor, or a Full-vs-split question | local fix to the line |
| something outside the set, or rank-deficiency refusal | basis incomplete, reduced-ME argument order wrong, **or the block structure is wrong** | do NOT write the number down; re-audit the structure first |

The third row is what a bare measurement cannot distinguish. Judge
membership, not decimals — `0.667` is `2/3` and fine; collinear fits are
only good to ~3 digits (below), so a decimal threshold would reject real
coefficients.

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

**Where mode 3 structurally cannot be used**: a BARE ANTENNA target in
a COLLINEAR limit. Every basis entry then carries the same `1/s_ij`
pole, the design matrix collapses to rank 1, and only the coefficient
SUM is determined. The documented cure (basis = antenna × reduced ME)
lifts the degeneracy through the reduced MEs, so it applies only when
the target is a matrix element; with a bare antenna there is nothing to
lift it. Symptom: coefficients that drift with x and are not rational,
while the residual still looks small. Use mode 4 instead — do not
"fix" it by shrinking the basis to one entry and reading the number,
which silently discards the question you were asking.

## Mode 4 — z-differential residue profile (residue_profile.py)

The collinear counterpart of mode 3: do not solve, PROFILE. Bin points
by the collinear momentum fraction z (measured from invariants,
`z = s(i,r)/(s(i,r)+s(j,r))` with r a hard spectator) and report
`median(target/candidate)` per z-bin with its spread, one candidate at
a time — no linear system, so no degeneracy.

```bash
python .claude/skills/probe-me-ir-structure/scripts/residue_profile.py spec.json > profile_gen.f
# compile/link exactly like a probe (below); run with OMP_NUM_THREADS=1
```

Reading the profile — this is the whole point of the mode:

| profile | conclusion |
|---|---|
| flat in z, rational constant | candidate is right; the constant is its coefficient |
| flat in z, irrational constant | right species, wrong normalisation or a missing partner term |
| **varies with z** | wrong species, OR the candidate mixes soft and collinear pieces that the target weights differently — the signature that a **sub-antenna split** is needed, not a coefficient |
| flat but wide per-bin spread | the bulk cancels and a TAIL does not. Two causes, distinguish before concluding: (a) gluon parent + wrong spin structure (azimuthal terms survive the built-in `rotp` average); (b) jet-cut acceptance differing between target and candidate, which zeroes one of them on part of the sample — check by re-profiling with `cuts_call` removed. A ~10% spread with a quark parent is (b) far more often than (a) |

The middle two rows are the ones a least-squares fit cannot report at
all: it returns a number either way.

**Profile COMPLETE decompositions, not single lines.** Alternative
complete decompositions of the same block (e.g. gluon-radiator
`G30 x redME` vs the quark-radiator pair `1/2(E30+E30) x redME`) are
EACH flat at 1.000000 in the limit they cover — a per-line profile
cannot rank them, and the tightest single-line spread is not evidence
(observed: it selected the leading-colour arrangement for a term whose
reference implementation uses the quark-radiator pair). Put each full
candidate decomposition in ONE candidate via `"terms"` (several
expr/map pieces summed), and decide between decompositions by what
they force DOWNSTREAM: the radiator choice fixes which X40 family the
S,b1/S,b2 blocks must use (write-subtraction). Spec fields mirror the other
generators (`coll` names the collinear pair, `zref` the hard reference,
`candidates` a list of independent expressions).

## Confidence — state of validation

- **Soft limits: high confidence.** Expect EXACT RATIONALS, and
  specifically members of the family above. Non-rational output means
  the basis is incomplete, a reduced-ME argument order is wrong, or the
  block structure is wrong — that is itself the diagnostic.
- **X40-level fits in double-unresolved limits: VALIDATED.** They
  behave exactly like the soft-gluon case — exact rational
  coefficients, residuals at the 1e-14 level.
- **Collinear limits, ME target**: need the basis construction below
  (basis = antenna × reduced ME, un-divided); judge by coefficient
  stability across an x scan, treat 3 stable digits as a good result.
- **Collinear limits, ANTENNA target**: mode 3 does not apply at all
  (rank-1 design matrix). Use mode 4 and read the z-profile.
- Not a validator for a subtraction term — that is run-spike-test.

## Usage

```bash
python .claude/skills/probe-me-ir-structure/scripts/gen_probe.py spec.json > probe_gen.f
cd test/process/<PROC>
gfortran -c probe_gen.f -J obj/mod -I obj/mod            # seconds — reuse
gfortran -o probe_gen probe_gen.o $(ls obj/*.o | grep -v <checkprog>.o)  # spike-test objects
OMP_NUM_THREADS=1 ./probe_gen
```

The spec is ~15 lines (`scripts/example_spec.json` for mode 1,
`assets/example_spec_residue.json` for mode 3 — the latter also shows
the soft-while-collinear `sc` boundary driver); a hand-editable
template with identical structure is `assets/probe_template.f` for
one-offs. Both scripts have a
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
  and conditioning IMPROVES as x -> 0. Deep x (down to ~1e-10) is
  allowed here — but the validated object is still STABILITY of the
  coefficients across the scan (typical range 1e-7...1e-9, deeper as a
  check), not one deep evaluation; if the deepest x goes numerically
  unstable, report the x at which it sets in rather than dropping it
  silently. Mind the degeneracy caveat above when the basis contains
  symmetry-related orientations of one antenna.

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
