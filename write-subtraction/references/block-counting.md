# Why the blocks are what they are — the counting argument

Load this when a block's EXISTENCE, SIGN or COUNT is in question, or when
a prediction from `predict_blocks.py` disagrees with a measurement. It is
not needed to write a term whose structure is already settled.

Source: hep-ph/0505111 §2.3 (eqs. 2.16–2.28), 1301.4693 §3.1 (eqs. 3.9–3.13),
0710.0346 §3 for the large-angle soft terms. Every claim below is a
statement about the SCHEME, not about a slot convention — conventions are
still measured (probe-me-ir-structure).

---

## 1. The one fact everything follows from

Start from `dσ^{S,a}` alone: the NLO subtraction term with the jet
function demoted from `J_n^{(n+1)}` to `J_n^{(n)}`, so an extra parton is
allowed to go unresolved while still forming n jets. Ask what it does in
each double-unresolved region. The paper's answer (0505111 §2.3.1):

> **`dσ^{S,a}` yields TWICE the (m+2)-parton matrix element in all
> colour-unconnected, almost-colour-unconnected and colour-neighbouring
> double unresolved limits, while VANISHING in all genuinely
> colour-connected double unresolved limits.**

The factor of two is not an accident and does not need measuring. A
double-unresolved limit requires both `p_j` (the antenna's own unresolved
leg) and one other momentum `p_o` in the reduced matrix element to become
unresolved. **The roles of `p_j` and `p_o` can be interchanged**, and both
assignments are present in the sum — two identical terms.

Everything else is bookkeeping against that sentence:

| region | what S,a does | therefore |
|---|---|---|
| genuinely colour-connected | vanishes | need a positive block: **S,b1** (`+X40`) |
| colour-neighbouring | doubles | the `X30×X30` products are non-zero and cancel it: **S,b2** does double duty |
| almost-colour-connected | doubles | need a negative block: **S,c** (`−`) |
| colour-unconnected | doubles | need a negative block: **S,d** (`−`) |

So the signs are derived, not conventional: **S,b1 is `+` because S,a left
a hole; S,c and S,d are `−` because S,a over-filled it.** If you find
yourself wanting a `+` on an S,d line, the classification is wrong, not
the sign.

---

## 2. The four classes, and the fifth thing hiding in one of them

Classification is of the PAIR of unresolved partons, with respect to a
single colour ordering (chain distance `d` between them):

- **`d = 1` colour-connected** — the two sit between one pair of hard
  radiators. One `X40` (S,b1), plus **two** `X30×X30` counterterms
  (S,b2), one per single-unresolved limit of that `X40`
  (0505111 eq. 2.17 subtracts `X⁰_{ijk}X⁰_{IKl}` *and* `X⁰_{jkl}X⁰_{iJL}`).
- **`d = 2` almost-colour-connected** — colour-disconnected but sharing
  one radiator. **Two** lines, one per strong ordering (eq. 2.24 writes
  both `X⁰_{ijk}x⁰_{mlK}` and `X⁰_{klm}x⁰_{ijK}`), overall `−`.
- **`d ≥ 3` colour-unconnected** — disjoint dipoles. **One** line per
  **unordered** pair, overall `−`. The paper's restriction is explicit:
  *"such that no product of two antenna configurations appears twice."*
  A double sum here double-counts.
- **colour-neighbouring** — NOT a fourth class of pair. It is a property
  of a **pair of CLUSTERS**, arising inside `d = 1`: the two unresolved
  partons cluster OUTWARD onto their respective radiators, so two
  chain-adjacent clusters go collinear independently. Concretely, on a
  chain `… i, j, k, l …` with `j, k` unresolved, it is the configuration
  `(i∥j)` and `(k∥l)`. It exists only when `i` and `l` are distinct and
  **both** outward clusters are flavour-valid (net flavour a single
  parton) — so a `d = 1` pair does not automatically have one. Where it
  does exist, the `X40` correctly approximates the ME **and** the
  `X30×X30` products are non-vanishing, each equal to the
  double-unresolved limit. `predict_blocks.py` applies exactly this test
  and reports the cluster pair it found.

**Why colour-neighbouring matters operationally.** The skill's working
rule for S,b2 — *"indexed by (S,a line, X40 line) pairs that share a
singular invariant"* — is a proxy for "cancel the X40's spurious single
limits". In a colour-neighbouring configuration the same line has a
**second job**: cancelling S,a's doubled double limit. Pole-sharing sees
only the first job, so the proxy can under-count exactly there.

*Signature of getting this wrong:* every single-unresolved mode passes and
the colour-connected double limits pass, but a **double-collinear** mode
whose two clusters are chain-adjacent sits at ~2 (S,a's doubling left
uncancelled) or ~0.5. `predict_blocks.py` flags these pairs
`[COLOUR-NEIGHBOURING]`.

---

## 3. Sub-antennae in S,c: the half-eikonal rule

S,c does not use full antennae. 0505111 eq. (2.24) uses `x⁰`, a
**sub-antenna**, defined by two properties:

- it contains **only** the `m∥l` collinear limit, **not** the `l∥K` one;
- in the soft limit of `l` it yields **half** the soft eikonal factor.

The reason is the criterion you actually need when choosing which half to
use: **`K` is a MAPPED (composite hard) momentum, and the matrix element
has no collinear limit with a mapped momentum.** A full antenna would be
singular where the ME is finite — a spurious singularity that no
counterterm can repair.

Generalised, and worth stating as a rule: **a mapped momentum must never
appear as a collinear partner in a subtraction line.** That single test
decides most Full-vs-split questions before any measurement, and it is
the missing half of the skill's "measure each half's pole graph" step —
the pole graph tells you what a half covers; this tells you what it is
*allowed* to cover.

The eikonal itself carries one more freedom, easy to miss (0710.0346):

> *"The hard momenta a, c do not need to be equal to the hard momenta i, k
> in the antenna phase space — **they can be arbitrary on-shell
> momenta**."*

with `S_abc = 2 s_ac / (s_ab s_bc)` (implemented as `SS1` in
`src/X30/SS1.f`, which reads its invariants from the cross-set
`s{ipset}on{jpset}` commons — mapped and unmapped legs mixed in one
call). That freedom is what generates the six-term difference structure
`S_{mapped,mapped} − S_{half-mapped} − S_{unmapped} + …`: the block is a
difference between the eikonal the TRUE ME dipole has and the eikonal the
iterated `X30×X30` product actually produces, and those two live on
different momentum sets.

---

## 4. Where the large-angle soft terms came from

Worth knowing because it prevents a wrong mental model.

hep-ph/0505111 has **no** large-angle soft term — its decomposition is
exactly `S = a + b + c + d` and no eikonal object appears. They were added
in the **revised** 0710.0346, for a specific reason:

> *"in the N² and N⁰ colour factor, the **angular averaging is not
> sufficient** to cancel the 1/ε poles in the four-parton one-loop
> subtraction terms."*

i.e. the subtraction terms **themselves** introduce spurious large-angle
soft limits, which angular averaging does not remove. 1301.4693 then
folds them into S,c (its §3.1 defines S,c as almost-colour-connected
*"and including large angle soft radiation"*), which is the grouping the
repo uses. Appendix B of that paper is devoted to the integrated IF
large-angle soft term as a distinct ingredient; in-tree that is
`src/X30int/SSintIF.f`, documented in `doc/LAST/SSintIF.tex`
("LAST" = large-angle soft term).

**Consequence for debugging:** an unresolved 1/ε in a colour factor where
the angular average "should" have worked is not necessarily a broken
antenna — it may be a missing S,c/LAST block.

---

## 5. Coefficients are unit fractions, and you can predict which

Antennae are **normalised colour-ordered matrix elements of a parent
process** — not functions built by imposing limits. 0505111 §1:

> *"If normalised appropriately, these full four-parton tree-level and
> three-parton one-loop matrix elements can be interpreted as antenna
> functions at NNLO."*

Parents: `γ* → qq̄ + partons` (quark–antiquark), `χ̃ → g̃g + partons` via
the Haber–Wyler effective Lagrangian (quark–gluon), `H → gg + partons` via
`L = −(λ/4) H F²` (gluon–gluon).

The normalisation is then **1/(number of distinct antenna configurations
the parent contains)** — hep-ph/0502110 states it for the gluonic case:
the `H→ggg` matrix element *"contains three different antenna
configurations … the effect of the symmetrisation over the three gluons
is that these three antenna configurations are averaged over."*

This reproduces the FF block of `maple/form/common/J21.map` exactly (the
quark-loop / `N_F` family carries an infix **`h`**, not an `_NF` suffix —
grep `J21hQGFF`):

| antenna | parent | configurations | J21 entry |
|---|---|---|---|
| `A30` | `γ* → qq̄g` | 1 (quark endpoints fixed) | `J21QQFF = calA30FF` |
| `D30` | `χ̃ → g̃gg` | 2 (either gluon outermost) | `J21QGFF = 1/2*calD30FF` |
| `E30` | `χ̃ → g̃q'q̄'` | 2 | `J21hQGFF = 1/2*calE30FF` |
| `F30` | `H → ggg` | 3 (cyclic) | `J21GGFF = 1/3*calF30FF` |
| `G30` | `H → gqq̄` | 1 | `J21hGGFF = calG30FF` |

and, at `J22` level, `1/4·calF40FF` and `−1/9·calF30FF²` for `H → gggg`
(4 configurations; the product of two `1/3`s).

**The rule does NOT extend to every X40.** `J22.map`'s FF block contains
`calA40FF` 1, `1/2*calD40FF`, `calE40FF` 1, `calG40FF` 1, `1/2*calH40FF`,
`calB40FF` 1, **`2*calC40FF`** and **`−2/3*calG30FF*calF30FF`** — so the
observed set is not all unit fractions. Read `J22.map` for the case at
hand; do not extrapolate.

**Use.** A residue or dipole fit should return a member of the observed
set `{2, 1, 2/3, 1/2, 1/3, 1/4, 1/9}`. A fit landing **outside** it is
telling you the basis is incomplete, an argument order is wrong, or the
block structure is wrong — not that the fitted number belongs in the
`.map`. Judge membership, not decimals: `0.667` is `2/3` and is fine,
and collinear fits are only good to ~3 digits anyway. That turns an
open-ended fit into a hypothesis test.

**Same rule predicts the Full-composite split.** hep-ph/0502110:

> *"the matrix element has to be split into three individual antenna
> configurations. **Each individual antenna configuration contains only
> one soft limit.** Each collinear `g→gg` is **split between the two
> antenna configurations** appropriate to the two final-state gluons
> involved in the splitting."*

So: one soft limit per sub-antenna, and a `g→gg` collinear shared between
two of them. That predicts the split for an antenna you have never seen —
then measure the halves' pole graphs to confirm.

X40 coefficients are less well evidenced: only `A40 → 1` and `F40 → 1/4`
are directly readable in `maple/form/common/J22.map`. Read that file
rather than extrapolating.

---

## 6. Constraints the papers give for free

Each is a static check on a written `.map`, cheaper than any run:

1. **`1/S_n` is per channel, not per antenna.** It appears in every master
   formula attached to the FULL final-state multiplicity (`S_{m+2}` for
   RR, `S_{m+1}` for RV). Do not fold it into an antenna coefficient.
2. **S,d sums unordered pairs.** No product of two antenna configurations
   may appear twice.
3. **S,b2 has two lines per X40 PER COLOUR ORDERING**, one per
   single-unresolved limit of that X40. The per-ordering qualifier is
   load-bearing: a file covering several orderings or flavour sectors
   carries a multiple of that, which is what the pole-graph indexing in
   the SKILL resolves. What the count is NOT is "one per term in the
   reduced ME's NLO subtraction" — a different quantity, and wrong. The
   per-pole consistency check is likewise **per pole, not per X40** (a
   single X40 can have one pole cancelled by one family and another by a
   different family).
4. **S,c has two lines per almost-connected pair per ordering** — both
   strong orderings of that pair.
5. **No mapped momentum as a collinear partner** (§3).
6. **The iterated integral does not factorise.** 0505111, after eq. 2.25:
   the second antenna integral *"will pick up ε-dependent factors from the
   first integral (both integrals are fully independent only in four
   dimensions) … the analytic integration will **not** yield the product of
   two independent integrated NLO antenna functions."* Relevant the moment
   you reason about what an S,b2 line becomes one layer up.
7. **Angular averaging must close within a single phase-space mapping**
   (0710.0346 §3.4): *"this average has to take place within each phase
   space mapping."* This is a DESIGN constraint on how you may split a
   composite, not only a property of the spike-test harness — a split that
   sends one collinear limit through two different mappings is not
   rescued by averaging afterwards.

---

## 7. What this does NOT replace

The counting argument fixes **structure**: which blocks, how many lines,
what sign, which rational family the coefficient lives in. It says nothing
about:

- which argument slots of a Fortran antenna are radiators
  (`antenna_slots.py`, then the pole scan);
- which concrete `X40` reproduces a given iterated sum (residue fit — the
  paper's classification does **not** pick the letter for you);
- which of the two writings of an iterated counterterm is correct
  (measure; the wrong one is not rational);
- reduced-ME names and argument orders (`FLAVlist`, me-naming-convention,
  and measurement).

Predict, then measure. The value of predicting first is that the
measurement now has an expected answer, and a disagreement localises to
the structure instead of being absorbed as a fitted number.
