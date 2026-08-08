---
name: write-subtraction
description: >
  Write, edit, or fix an antenna-subtraction term at the maple level
  (maple/process/<DIR>/*.map) in NNLOJET — including building one from
  scratch when no analogous term exists in any neighbouring process:
  deriving the block skeleton from colour connection before any run,
  aligning an unfamiliar antenna's arguments with the cluster rule,
  constructing the S,b2 / S,c / S,d blocks from measured pole graphs, and
  iterating block hypotheses with the composer script. Also answers "why
  does this block exist / why this sign / how many lines". Use whenever the
  user asks to write, modify, or debug a subtraction term, fix a channel
  that fails a spike test, or add missing infrared limits to a
  subtraction. This skill covers ONLY the maple .map file; generating
  Fortran and registering it is the autogen-subtraction skill, and
  validation is the run-spike-test skill.
---

# Writing NNLOJET subtraction terms (maple level)

A subtraction term must mimic every infrared (unresolved) limit of the
corresponding channel of the full matrix element: each single-unresolved
limit for R, plus double-unresolved and one-loop×single limits for RR/RV.
A spike-test failure in a specific limit means the term line(s) covering
that limit are wrong or missing.

## Cost ladder — cheap lookups first

Before re-deriving anything by measurement or build cycle, walk this
ladder from the top; each rung answers in seconds what the rung below
answers in minutes-to-hours:

0. **Derive the block skeleton** — which blocks exist, how many lines
   each carries, their signs, and the rational family the coefficients
   live in are all FIXED by colour connection; they do not need
   measuring. `predict_blocks.py` (below) emits them in seconds. Do
   this FIRST: everything on the rungs below then has an expected
   answer, and a disagreement localises to the block STRUCTURE instead
   of being absorbed as a fitted number.
1. **Name lookups** — me-naming-convention (ME grammar, crossing
   relatives), antennae-naming-convention (species, layer,
   FF/IF/FI/II correspondences); list-channels / list-processes for
   what exists.
2. **Static scans** — `antenna_slots.py` (slot plumbing,
   antennae-naming-convention), `genuine_modes.py` (which modes are
   real, write-spike-test), precedent grep over generated `auto*.f`
   (section below).
3. **One-run measurements** — per-line attribution
   (`wt_attribute.py`, run-spike-test), pole scan, dipole fit
   (probe-me-ir-structure).
4. **Multi-run measurements** — residue fits (probe-me-ir-structure
   mode 3).
5. **Build cycles** — the block composer below plus
   autogen-subtraction's regenerate+rebuild wrapper. Never hand-roll
   the regenerate–restore–rebuild sequence; the wrapper is that
   sequence.

## Where things live

- `maple/process/<DIR>/*.map` — one dir per maple-level process.
- `maple/iprocess.map` — authoritative process list: iprocess number,
  directory, jet function (`MYJET`, e.g. `ecuts_epem_vh`, `ecuts_dis`),
  subprocess index ranges, extra Fortran modules. NOTE: mapping to
  test/driver processes is not 1:1 (maple `epemZH` = iprocess 320/330/340 →
  test-level `epemZH2bb` / `epemZH2gg` / `epemZH2gaga`).
- `maple/notation.map` — THE token dictionary (antenna sets, soft sets,
  integrated dipoles, ME multiplicity sets). Read it before inventing any
  token. Rendered: `maple/notation.pdf`.
- Review aids: `maple/process/<DIR>/view<PROC><RR|RV|VV>.tex` and generated
  `auto*.tex` — LaTeX rendering of every term.

## File naming

`<initial state><ME layer><n>g<l><PROC><Type>.map`, e.g. `qqbBt2g0ZS.map`.
Type suffix = contribution: `SNLO` = R (NLO real), `TNLO` = V (NLO),
`S` = RR, `T` = RV, `U` = VV. Decay processes use `x` as separator:
`B0g0ZHepemxBy1g0HSNLO.map` (production ME × decay ME).

## Anatomy

```maple
# Differential R subtraction term for
# e+(1) e-(2) -> q(3) qb(4) (H -> b(i) g(j) bb(k)).

FN:=B0g0ZHepemxBy1g0HSNLO(1,2,3,4,i,j,k):

XX:=
+A30FF(i,j,k)*B0g0ZHepemxBy0g0H(1,2,3,4,[i,j],[j,k])*JET22([i,j],[j,k])*a1
:
```

Every line of `XX` is:

```
±  antenna(...) × reducedME(mapped args) × JETnm(...) × aN
```

- **Antennae** (from `notation.map`; full grammar + limit coverage:
  antennae-naming-convention skill): tree 3-parton `ant30set` split by
  configuration — FF (`A30FF, D30FF, d30FF, E30FF, F30FF, ...`), IF
  (`gA30IF, qA30IF, qE30IF, ...`), FI, II (`qqA30II, ...`); 4-parton
  `ant40set` (`A40, At40, B40, D40, E40, Et40, ...`) for double-unresolved;
  one-loop `ant31set` (`A31FF, ...`) in RV only. Soft eikonal sums:
  `SSset = {SFF, SIF, SFI, SII}`. Integrated dipoles `J21*`/`J22*`
  (`J21QGFF`, ...) appear in RV/V terms.
- **Reduced MEs** are the Fortran-level ME names (`B0g0Z`,
  `B0g0ZHepemxBy0g0H`, ...; naming grammar and crossing rules:
  me-naming-convention skill). Mapped (clustered) momenta are bracketed:
  `[i,j]` = single mapping, `[[i,j],[j,k]]` = iterated mapping.
- **Jet functions** `JETnm(args)`: `n` = resolved jets, `m` = final-state
  partons in the reduced kinematics (`JET00 ... JET33`). Args must be the
  correct mapped momenta of the reduced ME, but their ORDER is free — the
  momentum list in the generated `set_map` call is canonical (i1, i2,
  then the cluster representatives in cluster order, then the remaining
  spectators in ascending index order) and is derived from the antenna
  cluster, not from the JET arguments (verifiable in any generated
  `auto*.f`: the emitted set_map list routinely differs from the JET
  argument order of its `.map` line). Do not agonise over JET ordering.
- **Which mapped momenta a cluster produces**: a 3-parton antenna
  cluster (a,b,c) yields `[a,b]` and `[b,c]`; a 4-parton cluster
  (a,b,c,d) yields `[a,b,c]` and `[d,c,b]`. Bracket contents are
  order-insensitive up to reversal (`[k,i,j]` ≡ `[j,i,k]`). Use these to
  write the reduced-ME and JET arguments. Source of truth: `getpmapIK`
  in `maple/getpmap.map` (~line 146).
- **`aN`** is a sequential line label: `a1, a2, a3, ...` with NO gaps (the
  checker aborts on gaps). Each `aN` becomes one `wt(N)` slot and one
  independent mapping/jet block in the generated Fortran — numbering is
  load-bearing; renumber consistently when inserting/removing lines.

Optional: `colflag:=true:` (groups several reduced MEs under one antenna),
`XX:=expand( ... ):` in some RV files.

**`colflag` is load-bearing, not cosmetic** — 474 `.map` files set it
(concentrated in `4jet/`, `WJJJ/`, `WJJ/`, `3jet/`), and it changes what
the generator emits in three ways (`makefortRR` ~lines 119–122, 327–340,
349–375, 396–400; `makefortRV` similarly). Default is `false`, reset per
subtraction term inside the loop:

1. it **skips `expand()`** when splitting on the `aN` labels, preserving
   the factorised `ant30*(ME_1 + ME_2 + ...)` shape instead of
   distributing it into separate terms;
2. it collects a **LIST** of matrix elements (`matM0:=[...]`) rather
   than a single one;
3. **the flavour mapping is taken from `matM0[1]` — the FIRST element of
   that list only** (`FLAVxx:=subs(pmap,matM0[1])`). This is the trap:
   reorder the MEs inside the bracket and you silently change the
   emitted `set_flav_perm`, with no error and no `.map`-level symptom.
   If a colflag term fails in a way that looks like a wrong flavour
   sector, check the FIRST ME in the group before anything else.

**Comment discipline — a `.map` is source, not a notebook:**

- The file HEADER carries only what is needed to READ the file: the
  process line, the momentum/flavour assignment of the `FN` arguments,
  and the reduced-ME conventions. A few lines.
- Inside `XX:=`, ONE comment line per BLOCK of subtraction lines,
  stating the block's origin (which limit family or colour structure
  it comes from) — the same `# block: <name>` markers the composer
  uses. No per-line commentary.
- NO status, validation results, diagnosis, or to-do lists in the
  `.map`. Those belong in the chat summary (run-spike-test's limit
  table) and, if they must persist, in a separate notes file — never
  in the term. A session that has just debugged a term must REMOVE its
  scratch commentary before finishing.

## Rung 0 — derive the block skeleton before measuring anything

Colour connection fixes the STRUCTURE of an RR term completely: which
blocks exist, how many lines each carries, every sign, and the rational
family every coefficient belongs to. None of that needs a run. Derive
it first, then let the measurement stack adjudicate.

```bash
# predict_blocks reads a SUPERSET of genuine_modes' spec, so one file can
# drive both — but they read different keys (partons/born vs
# chain/flavours) and the flavour info is NOT cross-validated between
# them. Keep them in sync by hand, or keep two files.
python .claude/skills/write-spike-test/scripts/genuine_modes.py spec.json --json > modes.json
python .claude/skills/write-subtraction/scripts/predict_blocks.py spec.json --modes modes.json
python .claude/skills/write-subtraction/scripts/predict_blocks.py spec.json --emit-markers > master.map
python .claude/skills/write-subtraction/scripts/predict_blocks.py spec.json --audit TERM.map
python .claude/skills/write-subtraction/scripts/predict_blocks.py --selftest
```

It prints, per unresolved pair, the colour-connection class and the
block it implies; per line, the sign, the FF/IF/FI/II configuration, the
antenna family hint and the expected coefficient. `--emit-markers`
writes a `# block:`-marked skeleton that feeds straight into
`map_blocks.py`; `--audit` compares a written `.map` back against the
prediction and reports missing, extra and mis-counted blocks.

**The one fact it encodes** (hep-ph/0505111 §2.3.1): `dσ^{S,a}` alone
**vanishes** in genuinely colour-connected double-unresolved limits and
yields **exactly twice** the matrix element in the other three classes.
Hence S,b1 is `+` (S,a left a hole) while S,c and S,d are `−` (S,a
over-filled it). **If a line seems to want the opposite sign, the
CLASSIFICATION is wrong, not the sign.** Derivation:
`references/block-counting.md` §1.

Two consequences worth carrying:

- **Counts are fixed, not discovered**: two S,b2 lines per X40 (one per
  single-unresolved limit of that X40), two S,c lines per
  almost-connected pair (both strong orderings), one S,d line per
  **unordered** disjoint pair — no product of two antenna
  configurations may appear twice.
- **Coefficients come from a small discrete rational set**, because an
  antenna is a normalised colour-ordered ME of a parent process. The
  observed set across `maple/form/common/J21.map` and `J22.map` is
  `{2, 1, 2/3, 1/2, 1/3, 1/4, 1/9}` — note it is **not** all unit
  fractions (`2*calC40FF` and `-2/3*calG30FF*calF30FF` are both live in
  `J22.map`). Use it as a MEMBERSHIP test on a fitted coefficient, never
  as a prediction of which member applies: a fit landing outside the set
  means the basis is incomplete, an argument order is wrong, or the
  structure is wrong — not that the fitted number should be written
  down. Derivation and the per-antenna table:
  `references/block-counting.md` §5.

**Scope, so it is not over-trusted.** The prediction is per COLOUR
ORDERING (adjacency is the right criterion for block assignment inside
one ordered line, and the wrong one for deciding which limits are
genuine — the two uses are opposite, see step 1 below); counts are
per-ordering upper bounds. It does not predict antenna letters, slot
conventions, which X40 to use, or reduced-ME argument orders. Those are
measured, as before.

Full derivation, the colour-neighbouring sub-class, the half-eikonal
rule and the free static checks: `references/block-counting.md` (load
when a block's existence, sign or count is in question, or when a
prediction and a measurement disagree).

## Building the lines: from the ME's infrared limits

**First, pin the reduced-ME set.** Look up your full ME in
`maple/FLAVlist<iprocessname>.map` (iprocessname from
`maple/iprocess.map`; the file is GENERATED by
`maple makeflavlist -Diprocess=<N>` and consumed by all three code
generators — see autogen-subtraction). Each entry is
`[[incoming], FullME, [reduced MEs]]`. Read the structure off it
before writing any line: both `f1X*` AND `f2X*` present ⇒ the term
needs two flavour sectors; plain vs `t`-modified reduced MEs ⇒ the
colour level of the reduction (compare any pair of neighbouring
leading/subleading entries in the file). If the entry is missing,
derive it: apply the reduced-Born rule below to every limit; the union
of the resulting reduced MEs is the entry. Honest scope: the set fixes
WHICH reduced MEs appear — not how many lines, not the colour
ordering, not the argument order.

**Then MEASURE the colour structure** — do not read it off amplitude
source (colour-ordering sums plus interference defeat source reading):
the probe-me-ir-structure skill fits ME/redME onto a dipole basis and
returns exact rational coefficients in minutes. It answers which
dipoles the S,a block needs, f1↔f2 symmetry questions, and coefficient
normalisations, before any line is written.

**Measurement ADJUDICATES the rung-0 prediction; it does not replace
it.** Come to every fit with an expected answer — the block, the sign
and the coefficient family from rung 0 — and treat the fit as a
hypothesis test rather than an open question. The three outcomes are
diagnostic, and the third is the one that is otherwise missed:

- fit returns the predicted rational → confirmed, proceed;
- fit returns a different member of the family (`1/2` where `1` was
  predicted) → a symmetry factor or a Full-vs-split question; local,
  fix the line;
- fit lands **outside** the family, or refuses on rank-deficiency → do
  NOT write the fitted number. Either the basis is incomplete, a
  reduced-ME argument order is wrong, **or the block structure itself is
  wrong** — and only a prediction to compare against makes that third
  cause visible at all. Re-run `predict_blocks.py --audit` on the file
  before touching coefficients.

Judge membership, not decimals: `0.667` is `2/3`, a live coefficient,
while `0.7` is nothing. Collinear fits are only good to ~3 digits
(probe-me-ir-structure), so decide against the rational set, never
against a decimal threshold.

The line list is DERIVED from the full ME, not invented:

1. **Enumerate the limits of the full ME — by the reduced-Born rule**
   (empirically validated, 240/240 mode classifications on epem
   C1g0/Ct1g0/B3g0). A configuration is a genuine limit iff, after
   replacing every collinear cluster by its parent parton (a cluster is
   valid iff its NET FLAVOUR is a single parton: q∥g, g∥g, same-flavour
   q∥q̄; cross-flavour q∥Q is dead) and deleting every soft parton
   (gluons; or a same-flavour q q̄ PAIR going soft together), what
   remains is a legal Born state of the process. Consequences:
   - a SINGLE quark never goes soft; g→qq̄ splittings have a collinear
     limit but no soft one;
   - a same-flavour qq̄ pair has a double-soft limit ONLY if deleting
     it leaves a Born (in C1g0 both qq̄ pairs qualify; in B3g0 deleting
     q̄q leaves ggg — no Born, dead);
   - composites are NOT decomposed pairwise: {q,q̄,Q} is a genuine
     triple-collinear cluster (parent Q) although its q∥Q pair is
     dead, and conversely two individually-genuine pairs can combine
     into a dead double limit if collapsing both breaks the Born;
   - do not enumerate from colour adjacency in the argument list —
     colour-summed MEs connect non-adjacent partons. Adjacency only
     guides which radiator pair an individual (colour-ordered) LINE
     uses; for those, read the analogous existing `.map`.
   Write no line for a non-limit — a spurious subtraction is as wrong
   as a missing one (it survives integration un-cancelled). Note the
   spike test will NOT flag it directly: on dead modes ratio→0 or
   noise is normal harness behaviour either way; the layer check and
   integration are where a spurious term shows.
2. **One antenna per limit CLUSTER, not per limit.** Pick the antenna
   whose two hard radiators are the neighbours of the unresolved
   parton(s) and whose species match (antennae-naming-convention);
   arguments are (hard, unresolved, hard) — middle leg unresolved. One
   `A30FF(a,b,c)` line covers b-soft AND both a∥b, b∥c collinear at
   once; do not write separate lines per limit. Use sub-antennae
   (`d30`, `f30`) where a full antenna would double-count limits shared
   between overlapping clusters. **Before using an antenna you have
   not used before, run the pole scan** (probe-me-ir-structure): the
   letter fixes the species but NOT the slot convention (which
   arguments are radiators, which unresolved, which dipole the soft
   limit sits on) — do not infer the convention from the letter, the
   paper, or a sibling `.map`. Then follow the alignment procedure
   below.
3. **Build the reduced ME** from the cluster rules above: remove the
   unresolved parton, substitute the mapped momenta for the radiators;
   the reduced ME is the (n−1)-parton amplitude of the resulting
   flavour content (find its name via me-naming-convention). JET
   arguments = the reduced final-state momenta (order free).
4. **RR assembly by block** (the dσ^S decomposition of Fig. 3,
   arXiv:1301.4693). For each PAIR of unresolved partons, its COLOUR
   CONNECTION decides the block:
   - **S,a** — `X30 × M_{n-1}` for each single-unresolved cluster;
   - **S,b1** — COLOUR-CONNECTED pair (the two unresolved partons share
     a radiator between them): one `X40 × M_{n-2}`;
   - **S,b2** — MINUS the iterated `X30 × X30 × M_{n-2}` overlap
     between S,a and S,b1, with ±1/2 symmetry factors where clusters
     are symmetric. **Construct it from measurement, not from a
     shortcut**: the set of iterated counterterms is indexed by
     **(S,a line, X40 line) pairs that share a singular invariant**.
     Take the pole graphs of both (probe-me-ir-structure pole scan);
     if they share no pole, there is no overlap and no counterterm.
     For each pair that does share one, the counterterm is the product
     of the two antennae with the second evaluated on the first's
     mapped momenta, and its coefficient is fixed by requiring
     cancellation on the shared boundary — verify with the residue
     fitter rather than by pattern-matching a neighbouring file. The
     COUNT of counterterms is a property of the pole graphs, NOT of
     the term count in the reduced ME's NLO subtraction. (The rule
     "minus the S,a antenna times the NLO subtraction of its reduced
     ME" agrees with this only when the reduced ME has a single
     unresolved sector; with several sectors it generates counterterms
     for overlaps that do not exist, which destroys limits that were
     previously exact.)

     **Reconciling this with rung 0.** These are two different counts,
     not a disagreement. The rule REJECTED above is "one counterterm per
     term in the reduced ME's NLO subtraction" — that remains wrong.
     What rung 0 fixes is a different quantity: TWO per X40 per colour
     ordering, one for each single-unresolved limit **of that X40**
     (eq. 2.17 subtracts `X30_{ijk}X30_{IKl}` *and* `X30_{jkl}X30_{iJL}`)
     — and that is precisely what the pole graphs must reproduce. Use
     rung 0's number as the per-ordering expectation and the pole graphs
     to resolve the multiplicity across sectors and orderings; if the
     two differ by more than the number of orderings, one is wrong and
     it is worth finding out which before building.

     **The colour-neighbouring case, where pole-sharing under-counts.**
     Inside the colour-connected class there is a sub-case (0505111
     §2.3.1): two neighbouring pairs going collinear independently —
     one pair inside the antenna, the other formed by the remaining
     antenna momentum and its colour-connected neighbour. There the
     `X30×X30` products are NOT vanishing; each equals the double
     unresolved limit, and the line has a SECOND job beyond cancelling
     the X40's spurious single: it cancels S,a's DOUBLED double limit.
     Pole-sharing sees only the first job. `predict_blocks.py` flags
     these pairs `[COLOUR-NEIGHBOURING]`. Empirical signature of
     missing it: all single-unresolved modes pass, the colour-connected
     double limits pass, and a DOUBLE-COLLINEAR mode whose two clusters
     are chain-adjacent sits near 2 (or near 0.5) rather than 1.

     **The two writings of an iterated counterterm are NOT
     equivalent.** An overlap between two S,a lines can be written
     with either antenna first, the second evaluated on the first's
     mapped momenta — and the two writings are numerically different
     objects; the wrong one over-subtracts by an O(1) factor, not a
     small residue. The choice is determined, not free: **the correct
     writing is the one that is regular in every single-unresolved
     limit that the other antenna alone already reproduces exactly.**
     Decide it by measurement, in one run per candidate, without
     composing a new term: test both writings with per-line
     attribution (run-spike-test's `wt_attribute.py`) in the affected
     single-unresolved modes, or fit both as basis entries in a
     residue fit (probe-me-ir-structure mode 3) — the wrong one is
     not rational. Empirical signature: correct ordering →
     `1.000000` with zero outliers in the affected single-unresolved
     modes; wrong ordering → a large-magnitude, often sign-flipped
     ratio in exactly those modes, while the DOUBLE-unresolved modes
     may still look plausible (they do not discriminate).
   - **S,c** — ALMOST-COLOUR-CONNECTED pair (separated by one hard
     radiator): a large-angle soft correction, written as SS-difference
     blocks multiplying `X30 × M` lines, shape
     `(SFF(..)+SFF(..)−SFF(..)−SFF(..))*X30(..)*M`. **Construction,
     not just shape**: S,c is the difference between the exact eikonal
     of the TRUE ME dipole and the eikonal the iterated X30×X30
     product actually produces in the soft limit. Both are measurable:
     fit the soft residue of the X40 and of the iterated product
     separately (residue fitter) and take the difference; then confirm
     numerically that the candidate SS combination reproduces that
     difference BEFORE writing it into the `.map`.

     **Two structural rules that decide most of it before measuring.**
     (i) *Half-eikonal / sub-antenna rule* (0505111 eq. 2.24): S,c uses
     a SUB-antenna `x30`, not a Full composite. The sub-antenna
     contains only the `m∥l` collinear limit and NOT `l∥K`, and yields
     **half** the soft eikonal in the `l`-soft limit. The reason is the
     usable criterion: **`K` is a MAPPED momentum, and the matrix
     element has no collinear limit with a mapped momentum** — a Full
     antenna would be singular where the ME is finite, a spurious
     singularity no counterterm can repair. Generalised: *a mapped
     momentum must never appear as a collinear partner in a subtraction
     line.* That test picks the half; the pole graph then confirms what
     it covers. (ii) *The eikonal's hard legs are free*: in `S_abc =
     2·s_ac/(s_ab·s_bc)` the hard momenta `a, c` **need not be the
     antenna's own radiators** (0710.0346) — which is why the block can
     mix momentum sets at all, and why `SS1` takes a `(jpset, ipset)`
     pair rather than one kinematics index (wiring below).
     Count: **two S,c lines per almost-connected pair per colour
     ordering**, one per strong ordering (eq. 2.24 writes both).

     **S,c wiring at the Fortran level** (none of it inferable from
     `notation.map`; a reduced tree may contain no surviving example
     to pattern-match, so the accepted syntax and its emitted Fortran
     are carried as assets in this skill:
     `assets/sc_block_skeleton.map` + `assets/sc_block_emitted.f`):
     the generator routes a term through the soft path when its
     antenna content is SS-set functions times at most one X30
     (makefortRR branches "sum SS * MM0" / "X30 * sum SS * ML0",
     which toggle `insoft` and take the soft leg from the SS middle
     argument), emitting a `wtsoft` accumulation of
     `SS1(j1,i3,j2,jpset,ipset)` calls — radiators j1,j2 on the
     MAPPED set, soft leg i3 on the ORIGINAL set, invariants from the
     `s{ipset}on{jpset}` cross-set commons filled by the `fillson*`
     calls of the emitted set_map chain, which the caller must
     guarantee ran. `SS(i1,i3,i2,ipset)` is the unmapped-radiators
     variant: it reads the soft momentum from `common /soft/
     psoft(4)`, filled by `makesoft(i1,ipset)` (src/map/libmap.f) —
     which currently has NO caller in the tree, so a hand-written
     probe using `SS` must call it itself. Verify a new S,c block BY
     REGENERATING, not by assertion: find any process with an SS
     block (`grep -rl 'SFF(' maple/process --include='*S.map'`), run
     `makeRRcheck`/`makefortRR` for it, and compare the emitted
     `wtsoft` block against the asset's shape and your own.
   - **S,d** — COLOUR-UNCONNECTED pair (two disjoint dipoles): a plain
     product `X30 × X30 × M_{n-2}`, one line per pair of DISJOINT
     clusters that are each singular (pole graphs again) — nothing
     genuinely new at NNLO, but the lines must be there. The sum is
     over **UNORDERED** pairs: the paper's restriction is explicit
     (0505111 eq. 2.26, *"such that no product of two antenna
     configurations appears twice"*), so a double sum over both
     clusters double-counts every line. **Diagnostic
     signature of a missing S,d**: the double-collinear (and
     double-soft-adjacent) modes of DISTINCT radiator pairs fail while
     all single limits and the colour-connected double limits pass,
     and the failing ratio ME/subtraction sits ABOVE 1 and grows as x
     decreases (under-subtraction).

   **The minimal buildable/testable unit is S,a + S,b1 + S,b2 for one
   flavour sector — never S,a alone.** The X40 carries a spurious
   single-unresolved pole and the iterated X30×X30 counterterm cancels
   it: they only work as a pair, and the iterated terms are the SAME
   order as S,a in a single-unresolved limit, not subleading. Verified
   failure signatures of partial builds: S,a alone → soft-gluon mode
   reads ~0.4 (X30×M lines still contain the gluon in their reduced
   ME); S,a + iterated only → WORSE, negative O(1) ratios. All three
   together → 1.000000 first try. Consistency check on the pairing —
   **per pole of the X40, not per X40**: a single X40 can have one
   pole cancelled by one family of iterated terms and another pole by
   a different family, so a global "X40 coefficients sum against
   iterated coefficients" rule is neither necessary nor sufficient.
   Enumerate the X40's poles (mode-2 pole scan,
   probe-me-ir-structure); for EACH pole, the coefficients of the
   iterated lines sharing that pole must sum against the X40's
   coefficient on that pole. Verify the per-pole sums directly in
   your file rather than trusting any single reference file to exist.

   **Choosing WHICH X40**: sum what the iterated X30×X30 counterterms
   give in the collapsing limit, then pick the X40 whose MEASURED
   limit (pole scan + residue fit; `maple/notation.pdf` for
   orientation) reproduces that sum. Do NOT pick by colour-adjacency
   intuition — it over-includes: two colour-connected quarks do not by
   themselves imply a B40-type antenna if the iterated sum matches a
   different X40's limit.
   RV (dσ^T): **T,a** = minus the integrated counterparts (`J21`) of
   the R-term's antennae × `M_n` (cancels the RV poles); **T,b** =
   `X30 × M^{1-loop}` plus the `(X31 + X30·J21)` closures; **T,c** =
   the integrated S,c. Every S block reappears integrated in T or U —
   the arrow structure that run-layer-check verifies.
5. **Completeness check before generating**: every limit passing the
   reduced-Born rule must be covered by at least one line, AND every
   reduced ME in the FLAVlist entry must appear in the file — a
   FLAVlist ME with no corresponding line is a missing sector (the
   exact failure mode of forgetting an f2 block). Add the structural
   half of the check, which costs nothing:

   ```bash
   python .claude/skills/write-subtraction/scripts/predict_blocks.py \
          spec.json --modes modes.json --audit <TERM>.map
   ```

   It reports blocks predicted but absent (`MISSING`), present but not
   predicted (`EXTRA` — either a second colour ordering or a spurious
   block), and line-count mismatches. Counts are per-ordering upper
   bounds, so a colour-summed term may legitimately differ by a factor
   equal to the number of orderings — investigate any other
   discrepancy before generating. The check program's
   mode list (`stitle`) is a superset checklist — apply the rule to
   classify each mode (run-spike-test carries the same rule plus the
   empirical verdict via the |wt1| scaling exponent), then demand
   coverage only for the genuine ones. `makeRRcheck`'s
   `autoRRX40/M0/SS.map` split shows your classification back to you.

## Aligning an antenna's argument list with the cluster rule

The cluster rule `X40(a,b,c,d) → [a,b,c],[d,c,b]` presumes a colour
chain `a–b–c–d`. Whether a given Fortran antenna actually realises
that chain under those arguments is a SEPARATE convention, documented
nowhere — and it can differ between a `Full` composite and its own
sub-antennae. Writing a `.map` line therefore couples two independent
conventions; get the pairing wrong and the term is singular in limits
the ME is finite in, and every downstream counterterm attempt fails
for reasons that look like physics errors. The generic procedure
(applies to D40, G40, H40, and every IF/FI/II variant):

0. **Run the slot audit first** (`antenna_slots.py`,
   antennae-naming-convention): if the entry point is flagged — its
   Fortran dummy arguments are not declared in ascending positional
   order — write the emitted wrapper BEFORE anything else, and write
   the `.map` against the wrapper. This is a static scan; it costs
   seconds and would otherwise cost build cycles.
1. **Measure the antenna's pole graph** (probe-me-ir-structure pole
   scan).
2. **Read the chain off the pole graph**: the endpoints are the
   radiators, the interior legs the unresolved pair.
3. **Assign physical partons to slots so that BOTH hold**: every
   cluster implied by the cluster rule has single-parton net flavour,
   AND the physical colour chain matches the measured chain.
4. **Cross-check by measurement** (mandatory — this is the step that
   is easy to skip and expensive to skip): in the limit where a
   RADIATOR becomes unresolved, the antenna must reproduce the eikonal
   of a dipole that actually exists in the full matrix element (fit it
   with probe-me-ir-structure). If it lands on a dipole the ME does
   not have, the assignment is wrong — no counterterm can repair it.

**Splitting a `Full` composite on all-final clusters** (the companion
rule): the halves carry different momentum mappings AND cover
DIFFERENT limits — choosing one half silently drops a whole class of
limits while looking locally correct. Measure each half's pole graph
and confirm their UNION covers every limit the block is responsible
for. A half with fewer poles is not "safer" — it is incomplete.

Three rules that let you PREDICT the split instead of only verifying
it — apply them first, then measure to confirm:

1. **One soft limit per sub-antenna; a `g→gg` collinear shared between
   two.** That is the split rule, and it follows from what an antenna
   is (a normalised colour-ordered ME of a parent process, split into as
   many configurations as the parent has radiator-pair choices) — so it
   applies to an antenna you have never seen. Provenance and the
   coefficient table: antennae-naming-convention.
2. **No mapped momentum as a collinear partner** (the general form of
   the S,c half-eikonal rule above): the half whose collinear limits all
   involve genuine, unmapped legs is the correct one.
3. **Angular averaging must close within a SINGLE phase-space
   mapping.** 0710.0346 §3.4: *"this average has to take place within
   each phase space mapping."* This is a DESIGN constraint on the
   split, not just a property of the spike-test harness: a split that
   sends one collinear limit through two different mappings is NOT
   rescued by averaging afterwards, however well the union of pole
   graphs looks. The authors state they checked this explicitly for the
   decompositions of `E40` and `D40`. Symptom: a `g→gg` or `g→qq̄`
   collinear mode whose spread does not narrow as x decreases even
   though `rotp` averaging is on (run-spike-test reads this as an
   "azimuthal-rotation issue" — if the term splits a composite, suspect
   the split, not the harness).

Related: angular averaging is not always SUFFICIENT either. In some
colour factors it leaves uncancelled 1/ε poles, and that failure is the
historical origin of the large-angle soft (S,c) blocks — see
`references/block-counting.md` §4.

## Finding a precedent for an antenna

To see how an antenna has been used before, grep for it across
`maple/process/`, `src/process/` and `driver/process/` — INCLUDING the
generated `auto*.f` and library files. Generated Fortran survives when
`.map` sources do not (resets, never-written terms), and a generated
call site still shows the argument order, the mapping (`set_map`
list), the coefficient, and the counterterm it pairs with:

```bash
grep -rl 'FullG40' maple/process/ src/process/ driver/process/
grep -B5 -A10 'FullG40' src/process/*/auto*.f   # call site + set_map + wt()
```

## Iterating structural hypotheses: the block composer

Debugging RR structure means testing "does this block belong here?"
many times, and each test costs a regenerate–compile–link cycle; the
`aN` labels must stay gap-free, so hand-inserting or deleting lines
means renumbering everything below. Use the block composer instead:

```bash
# master file: your .map with '# block: <name>' comment markers
# (predict_blocks.py --emit-markers writes a starting master for you)
python .claude/skills/write-subtraction/scripts/map_blocks.py list  master.map
python .claude/skills/write-subtraction/scripts/map_blocks.py compose master.map \
       --blocks Sa,Sb1,Sb2_f1 -o <TERM>.map      # renumbers aN gap-free
python .claude/skills/write-subtraction/scripts/map_blocks.py --selftest
```

Keep ALL candidate lines in the master, group them under block
markers, compose the subset under test, then run the one-command
regenerate+rebuild wrapper (autogen-subtraction) and the spike test.
This turns a structural hypothesis into a single short test — and
makes the block bisection of run-spike-test practical.

## Deriving a term by crossing an existing one

Most "new" terms are not new: the same full ME usually exists in
another process with legs crossed to the initial state, and its `.map`
is the most efficient starting point. Find sources by grepping
`maple/process/*/` for the ME stem with crossing prefixes/suffixes
(e.g. a `C1g0Z*` term has one-leg-crossed relatives in `DIS/` and
two-leg-crossed ones in `WJ/`, `ZJJ/`); a one-leg-crossed source's
final-state sector typically maps ~1:1 onto one flavour sector of the
all-final term under a relabelling of the parton letters.

**A crossing is a DERIVATION, not a transcription.** Re-derive each of
these explicitly — every one is a silent, build-cycle-eating trap:

1. **Full vs split antennae.** An initial-state cluster is a pure
   momentum rescaling, written `[1]`/`[[1]]` — unambiguous — so a Full
   composite antenna is safe on a crossed leg. Un-crossed to all-final
   kinematics the two halves carry DIFFERENT momentum mappings: Full
   `E40` must become `E40a`+`E40b`, Full `D30FF` → `d30FF`+`d30FF`
   (one-leg-crossed DIS-type files show both patterns side by side:
   Full on `[1]` clusters, split halves on all-final clusters). After
   splitting, apply the union rule above: measure each half's pole
   graph and confirm together they cover every limit of the block.
   Symptom of transcribing: one dipole fails BOTH its single-collinear
   limits while its mirror dipole reads 1.0000.
2. **Iterated-counterterm cluster ordering** is not crossing-invariant:
   an X30 argument order that is unambiguous when one leg is the
   crossed `[1]` must be re-derived for the all-final version —
   clustered momentum FIRST, mirroring the corresponding X40-half
   quark cluster. The transcribed order leaves a 1.2–3.4× residual
   visible only in a couple of single-collinear modes. When in doubt,
   measure both orders with probe-me-ir-structure: the wrong one is
   not rational.
3. **Missing sectors.** The crossed file is INCOMPLETE by
   construction: sectors whose partons sit in the initial state there
   simply do not exist in it (WJ/ZJJ have no f2 sector; DIS only part
   of it) — and the file cannot tell you what is missing. Enumerate
   the full process's limits independently (FLAVlist + reduced-Born
   rule above) and write the absent blocks yourself.
4. **Antenna configuration renames**: FF↔IF↔FI↔II per leg
   (`qA30IF`↔`A30FF`, `qpE30IF`↔`E30FF`, `qD30IF`↔`D30FF`, and the
   conversion variants `ga30IFgtoq`, `gd30IF`, ...) — correspondence
   table in antennae-naming-convention.
5. **JET arity** follows the final-state parton count
   (DIS `JET22/JET23` → epem `JET33/JET34`).

Then validate as usual; the untouched-sibling control channel in
run-spike-test is what separates crossing artefacts from harness
effects.

(A fully worked instance of this whole workflow belongs in
`references/worked-example-C1g0ZepemS.md` — WHEN PRESENT. It contains
EXERCISE ANSWERS for the reset-C1g0ZepemS clean-room rebuild, so it is
routinely quarantined out of the repo; do not load it unless the user
explicitly asks, and withhold it when that rebuild is run as a test. If
the file is absent, that is the quarantine, not a broken link.)

(`references/block-counting.md` carries the scheme-level derivation —
the factor-of-two counting argument, the colour-neighbouring sub-class,
the half-eikonal rule, the coefficient family, and the static checks the
papers give for free. It contains no process-specific answers, so it is
safe to load during a clean-room rebuild; load it when a block's
existence, sign or count is in question, or when a prediction and a
measurement disagree.)

## Structural patterns per contribution type

(X30/X31/X40 here and above are generic CLASSES — X is a wildcard for
the concrete antenna letter A/D/E/F/G...; no antenna is literally named
X30. See antennae-naming-convention.)

- **`*SNLO` (R)**: lines of `X30 × M_{n-1} × JET × aN`, one per
  single-unresolved limit.
- **`*S` (RR)**: (a) `X30 × M_{n-1}` single-unresolved blocks;
  (b) `X40 × M_{n-2}` genuine double-unresolved;
  (c) `− X30 × X30 × M_{n-2}` iterated (removes double counting between a
  and b — typically negative sign);
  (d) `SS × M` soft eikonal blocks where needed.
- **`*T` (RV)**: (a) `− J21 × M_n × JET` integrated-dipole lines;
  (b) `X30 × M^{1-loop}_{n-1}` (tree antenna × loop ME);
  (c) `( X31 + X30 × J21 ) × M_{n-1}` (loop antenna + integrated Sb2).

The RR consistency split (`autoRRX40.map`, `autoRRM0.map`, `autoRRSS.map`,
generated by `maple makeRRcheck -Diprocess=N`) encodes which lines belong
to which class — regenerate and inspect it when debugging structure.

## Fixing a term after a failed spike test

The failing mode reported by run-spike-test names the limit (e.g.
"5||6 collinear", "6 soft", "triple collinear 567"). **Re-run rung 0
first** — `predict_blocks.py --audit <TERM>.map` costs seconds and, if
a whole block is missing or an extra one is present, tells you so
before any per-line reasoning. Then use the mode→block
diagnosis in run-spike-test: only blocks whose antennae have a
pole in that mode's invariants can be responsible, the ratio's
CHARACTER (stable offset vs wide spread vs tails) narrows the error
class, and block bisection with the composer above plus the
regenerate+rebuild wrapper localises it. Then check the suspect lines:
correct antenna type for the parton species (quark/gluon,
FF/IF/FI/II), correct mapped arguments in the reduced ME and JET
function, correct relative sign of the iterated `X30×X30` lines, no
missing limit entirely. When an analogous term exists in a
neighbouring process, compare against it; when it does not, use the
precedent-grep procedure above — and fall back on measurement, never
on guessing.

## Next step

After editing the `.map`, use the **autogen-subtraction** skill to
generate the Fortran and hook it into the spike-test build, then
**run-spike-test** to validate. Loop back here if a limit still fails.
