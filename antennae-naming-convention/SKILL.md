---
name: antennae-naming-convention
description: >
  Decode or construct NNLOJET antenna-function names (A30FF, d30FF, qA30IF,
  qqA30II, At40, A31FF, J21QGFF, ...), map each antenna to the infrared
  limits it subtracts, and explain the NNLO antenna-scheme structure: the
  dsigma^S/T/U blocks (S,a...S,d; T,a...T,c; U,A/B/C), colour-connected vs
  almost- vs un-connected pairs, and the J2 integrated dipoles / Catani I
  operators. Use for questions like "what is d30FF", "which antenna covers
  the 5||6 collinear limit", "what is J21 / J2^(1)", "what does S,b2 mean",
  "colour-connected vs unconnected", "what cancels the VV poles", or when
  interpreting spike-test failures and layer-check residues. The name
  grammar fixes an antenna's SPECIES only; slot PLUMBING is audited
  mechanically (antenna_slots.py — declaration order vs the positional
  cluster rule, with permuted-call recipes) and slot CONVENTIONS are
  MEASURED and CACHED in the antenna datasheet
  (antenna_datasheet.py — pole graphs, soft dipoles, split identities,
  reduced-cluster flavour semantics), falling back on the pole scan
  (probe-me-ir-structure) for antennae not yet in it. Read-only
  reference for the repo — modify nothing there.
---

# NNLOJET antenna naming and limit coverage

Authoritative token dictionary: `maple/notation.map` (sets `ant30set`,
`ant31set`, `ant40set`, `SSset`; rendered as `maple/notation.pdf`).
Fortran implementations: `src/X30` (tree 3-parton), `src/X31` (one-loop
3-parton), `src/X40` (tree 4-parton), `src/X30int/` (integrated
antennae — a FLAT directory of 14 generated files, configuration is a
filename SUFFIX not a subdirectory: `autoX30FFint.f`, `autoX30IFint.f`,
`autoX30FIint.f`, `autoX30IIint.f`, `autoJ21{FF,IF,FI,II}.f`,
`autoP0IF.f`, `autoP0FI.f`).

**Species vs convention — the one-line rule: before USING an antenna
you have not used before, look it up in the DATASHEET; if absent, run
the pole scan (probe-me-ir-structure).** This skill's letter grammar
fixes the SPECIES (which parton kinds, which limit families). Which
argument slots are the hard radiators, which are unresolved, how a
`Full` composite divides its singularities between its sub-antennae,
and which dipole its soft limit actually sits on are a separate
CONVENTION, and the two are independent. Most `src/X40/*.f` headers
state only a paper equation number, and `maple/notation.map` gives the
token, not the convention — do not infer it from the letter, the
paper, or a sibling `.map`; measure it (or read the cached
measurement).

## The antenna datasheet — measured operative description, cached

`assets/antenna_datasheet.json` (generated, never hand-edited) carries
per antenna the fields a term author needs and source headers do not
reliably state — each with provenance, re-verifiable at any time:

- MEASURED pole graph: singular 2-parton invariants with powers; which
  legs go soft and (3-parton antennae) which dipole the eikonal sits
  on, with snapped rational coefficient; singular triple-collinear /
  double-soft configurations. Absence from a list = the antenna does
  NOT have that pole (e.g. G30FF/E30FF carry no soft limit and no
  gluon–quark collinear — the poles an antenna does NOT have are as
  load-bearing as the ones it has).
- `requires_split`: whether the generic positional cluster rule
  `X40(A,B,C,D) -> [A,B,C],[D,C,B]` supports every measured pole. An
  ENDPOINT-pair pole s(1,4) cannot sit in either cluster, so a Full
  composite carrying one (FullE40's quark||gluon pole, FullD40,
  FullG40, FullF40) MUST be split into halves with their own mappings.
- for splittable composites: the halves and the NUMERICALLY VERIFIED
  calling identity (e.g. `FullE40(A,B,C,D) = E40a(A,B,C,D) +
  E40b(A,D,C,B)`; `FullF30FF(A,B,C) = f30FF(ABC)+f30FF(BCA)+
  f30FF(CAB)`) — a stale header comment cannot mislead, because the
  identity is re-checked on random phase-space points.
- declared positional slot species (quark-kind vs gluon), same-flavour
  pair constraints on the UNRESOLVED slots, and the reduced-cluster
  flavour semantics (what each cluster BECOMES — naive net-flavour
  arithmetic is WRONG for the E/G families and this field is the
  correct bookkeeping).
- the J21/J22 integrated-dipole coefficient lines mentioning the
  antenna, parsed from `maple/form/common/J2[12].map`.
- 4-parton soft limits factor as eik x reduced X30, not eik x const —
  the datasheet says so per leg and defers to the residue fitter
  (probe-me-ir-structure mode 3).

```bash
python .claude/skills/antennae-naming-convention/scripts/antenna_datasheet.py \
    show E40                                  # read the cached entry
python .../antenna_datasheet.py measure <NAME> --testdir test/process/<P>
python .../antenna_datasheet.py verify  <NAME> --testdir ...   # re-measure+diff
python .../antenna_datasheet.py static  --root .   # refresh declared fields
python .../antenna_datasheet.py residue <NAME> --limit sco:2,3 \
       --reduces-to "..." --method "..." [--confidence measured|inferred]
```

`measure` needs any BUILT spike-test object dir with 7-parton phase
space; the antenna is evaluated on its own legs, so entries are
process-independent. The datasheet feeds `pole_ledger.py`
(write-subtraction) — the static bookkeeping check of a whole `.map`.

### The operative contract — read this before placing an antenna

`show` prints, under `--- operative contract ---`, the four facts you
actually need to write a `.map` line. All of them are DERIVED from the
measured entry, so a newly measured antenna gets its contract free:

- **`call`** — the letter frame (A,B,C,D) and the declared slot order,
  so a permuted entry point cannot be miscalled;
- **`cluster`** — what each cluster the generic rule produces *becomes*:
  `[A,D,C] -> quark, flavour of A` means A absorbed D. **This is the
  legality test, not net-flavour arithmetic on the cluster members.**
  A cluster like `[qb, g, Q]` looks illegal by flavour counting and is
  perfectly legal here; rejecting an antenna on that reasoning is a
  known, expensive mistake;
- **`UNSUPPORTED` / `does NOT`** — poles in the graph that this
  antenna's clusters cannot represent, naming the partner half whose
  cluster can. That sentence is the difference between "E40b is
  unusable" and "s(A,D) is E40b's job";
- **`residue`** — what it reduces to on a boundary, which is the input
  to the matching iterated counterterm. Record one with `residue`
  whenever you establish it (a residue fit is `measured`; a cancellation
  inferred from a spike-test mode is `inferred` — say which, and the
  note field carries the evidence).

Halves always come as a pair; `show` says so explicitly and
`pole_ledger.py` turns it into an error with the calls to add.

**Slot plumbing — the cheap first lookup (run BEFORE the pole scan,
not instead of it).** The maple cluster rule is positional, but the
Fortran entry points do not all declare their `i1..i9` slots in
ascending positional order — reversed and fully permuted declarations
coexist in the tree, and a positional `.map` call against a permuted
one evaluates the right function with the wrong momentum map,
silently. Audit mechanically:

```bash
python .claude/skills/antennae-naming-convention/scripts/antenna_slots.py \
    --root <repo root> [--name <antenna>]        # --selftest available
```

It groups every entry point under `src/X30|X31|X40|X30int` by
slot-declaration order, FLAGS the permuted ones with an
ARGUMENT-PERMUTATION recipe (how a positional call realises the named
legs, e.g. `E40b(A,D,C,B)` for named `(i1,i3,i4,i5)=(A,B,C,D)`),
cross-checks `ant30set`/`ant31set`/`ant40set` against the Fortran in
both directions (tokens with no entry point; entry points never
registered), and scans for UNREGISTERED FORWARDING WRAPPERS — local
one-line renamings no `ant*fortset` resolves (e.g. an `E40bb.f`):
their header claims are unverified and rot silently, so measure the
TARGET instead of trusting the comment. **The script never emits a
slot-reordering wrapper**: a wrapper restoring ascending order gives
the right VALUE with the WRONG clusters under the generic positional
cluster rule (write-subtraction's argument-alignment rule) — permute
the `.map` ARGUMENTS instead, and let the datasheet's measured pole
graph / split identity decide which permutation is physically right.
Three different questions, in increasing cost order: the NAME fixes
the species, THIS SCRIPT fixes the slot plumbing, the DATASHEET /
POLE SCAN fixes the convention. A flag means "check before positional
use", not "wrong" — some species' canonical chains are legitimately
non-consecutive; the measurement decides.

## The name grammar

```
[crossed partons] <Letter> [t|tt] <30|31|40> <FF|IF|FI|II> [_g|_q|_frag|GtoQ]
```

**"X" is a WILDCARD, not an antenna.** The generic classes X30 / X31 /
X40 used throughout (papers, these skills, and the repo's directory
names `src/X30`, `src/X31`, `src/X40`, `src/X30int`,
`autoRRX40.map`, `calX30`, ...) mean "any tree 3-parton / one-loop
3-parton / tree 4-parton antenna". No antenna is literally named X30 —
concrete tokens always carry a letter: `A30FF`, `d30FF`, `E30FF`,
`B40`, `A31FF`, ...

- **`30`/`31`/`40`** = tree 3-parton (single unresolved) / one-loop
  3-parton (RV layer) / tree 4-parton (double unresolved).
- **Configuration** `FF|IF|FI|II` = where the two hard radiators sit
  (Final-Final, Initial-Final, ...). The leading lowercase parton
  letters name WHICH partons are crossed into the initial state:
  `gA30IF` (gluon initial), `qA30IF`, `qqA30II`, `qpG30IF`, ...
- **Lowercase letter** (`d30`, `f30`, `a30`) = sub-antenna; the
  uppercase antenna is the symmetrised sum, e.g.
  `FullD30FF = d30FF(i1,i2,i3) + d30FF(i1,i3,i2)` (src/X30/FullD30FF.f).
- **`t`** = tilde: subleading-colour contribution (same meaning as for
  matrix elements): `At40`, `Et40`, `Gt40`.
- **Suffixes**: `_g`/`_q` = fragmentation variants, `GtoQ` = g→q
  conversion variants (`ga30IFGtoQ`).

## Crossing between configurations (FF ↔ IF ↔ FI ↔ II)

Crossing radiators to the initial state maps an antenna to its
IF/FI/II variant; the crossed partons are named by the leading
lowercase letters. Observed correspondences (from `notation.map` and
usage, e.g. `DIS/qC1g0ZDISS.map`):

| FF | IF (one crossed) | II (two crossed) |
|---|---|---|
| `A30FF` | `qA30IF`, `gA30IF` | `qqA30II`, `qgA30II`, `gqA30II` |
| `D30FF` / `d30FF` | `qD30IF`, `gD30IF`, `gd30IF` | `qgD30II`, `ggD30II`, ... |
| `E30FF` | `qE30IF`, `qpE30IF` | `qpqpE30II`, `qqpE30II`, ... |
| `G30FF` | `qpG30IF`, `gG30IF` | `qqG30II`, `gqG30II`, ... |
| `F30FF` | `F30IF`, `gf30IF` | `F30II`, `ggF30II` |

plus conversion variants (`ga30IFgtoq`, `gd30IFGtoQ`, ...). FI mirrors
IF. These are DISTINCT functions with distinct integrated counterparts
— never substitute configurations.

**Inferring an FF antenna's flavour content from its crossed names.**
When source comments are silent, the crossed variants in
`maple/notation.map` spell the legs out: e.g. `calgE40IF_q`,
`calqpE40IF_q`, `calqbpE40IF_g` together imply E40 carries a gluon
leg, a quark leg, and a primed (secondary) quark pair — each crossed
name tells you which species was moved to the initial state. This
deduction generalises to any antenna. For a definitive numerical
check of flavour content, use probe-me-ir-structure.

**Cluster rule source of truth**: the mapping
`X40(A,B,C,D) → [A,B,C], [D,C,B]` (and the 3-parton
`(a,b,c) → [a,b],[b,c]`) is implemented in `getpmapIK` in
`maple/getpmap.map` (~line 146) — verify there, not in prose.

**The Full-composite rule for crossed legs**: an initial-state cluster
is a pure momentum rescaling, written `[1]`/`[[1]]` in the maple files
— unambiguous. Therefore a Full composite antenna (`E40`, `D30FF`) is
safe on a line whose cluster involves a crossed leg, but MUST be split
into its halves (`E40a`+`E40b`, `d30FF`+`d30FF`) on all-final-state
clusters, because the halves carry different momentum mappings AND
cover DIFFERENT limits: choosing one half silently drops a whole class
of limits while looking locally correct. Measure each half's pole
graph (probe-me-ir-structure) and confirm the union covers every limit
the block is responsible for — a half with fewer poles is not "safer",
it is incomplete. This is the central trap when deriving a term by
crossing an existing one — procedure and symptoms in
write-subtraction.

## Where the letters come from — and what that predicts

Antennae are not functions built by imposing limits. They are
**normalised colour-ordered matrix elements of a parent process**
(hep-ph/0505111 §1: *"If normalised appropriately, these full
four-parton tree-level and three-parton one-loop matrix elements can be
interpreted as antenna functions at NNLO"*), and the letter records
which parent:

| letters | parent process | hard radiators |
|---|---|---|
| `A`, `B`, `C` | `γ* → qq̄ + partons` | quark–antiquark |
| `D`, `E` | `χ̃ → g̃g + partons` (Haber–Wyler effective Lagrangian) | quark–gluon |
| `F`, `G`, `H` | `H → gg + partons` via `L = −(λ/4) H F²` | gluon–gluon |

Two things follow that are otherwise only discoverable by measurement:

**1. The normalisation is 1/(number of antenna configurations the parent
contains).** hep-ph/0502110, on the gluonic case: the `H→ggg` matrix
element *"contains three different antenna configurations … the effect
of the symmetrisation over the three gluons is that these three antenna
configurations are averaged over."* This reproduces the FF block of
`maple/form/common/J21.map` exactly (note the quark-loop/`N_F` family
carries an infix **`h`** — grep `J21hQGFF`, not `J21_NF`):

| antenna | parent | configs | J21 entry |
|---|---|---|---|
| `A30` | `γ*→qq̄g` | 1 (quark endpoints fixed) | `J21QQFF = calA30FF` |
| `D30` | `χ̃→g̃gg` | 2 (either gluon outermost) | `J21QGFF = 1/2*calD30FF` |
| `E30` | `χ̃→g̃q'q̄'` | 2 | `J21hQGFF = 1/2*calE30FF` |
| `F30` | `H→ggg` | 3 (cyclic) | `J21GGFF = 1/3*calF30FF` |
| `G30` | `H→gqq̄` | 1 | `J21hGGFF = calG30FF` |

**Use**: a dipole or residue fit should land in the small rational set
observed across `J21.map` and `J22.map` — `{2, 1, 2/3, 1/2, 1/3, 1/4,
1/9}`. It is **not** all unit fractions: `J22tQQFF` contains
`+2*calC40FF` and `J22hGGFF` contains `-2/3*calG30FF*calF30FF`. A
non-member is a diagnostic, not a number to write down.

**The X40 coefficients are not derivable from this rule** — read
`maple/form/common/J22.map` for the case at hand. Its FF block alone
gives `calA40FF` 1, `1/2*calD40FF`, `calE40FF` 1, `calG40FF` 1,
`1/2*calH40FF`, `calB40FF` 1, `2*calC40FF`, `1/4*calF40FF`,
`1/2*calAt40FF`, `1/2*calEt40FF`, `1/2*calGt40FF`.

**2. The Full-composite split is predictable.** hep-ph/0502110: *"the
matrix element has to be split into three individual antenna
configurations. **Each individual antenna configuration contains only
one soft limit.** Each collinear `g→gg` is **split between the two
antenna configurations** appropriate to the two final-state gluons
involved in the splitting."* So a soft limit belongs to exactly one
half, and a `g→gg` collinear is shared by two. That predicts the split
for an antenna you have never seen — then measure the halves' pole
graphs to confirm (the measure-before-you-use rule stands; what changes
is that you now arrive with an expected answer).

Neither point overrides the species-vs-convention rule above: provenance
fixes normalisation and limit content, **not** which argument slot is
which.

## Letter ↔ radiators ↔ limits subtracted

The letter encodes the hard-radiator pair and the unresolved
parton(s) — hence WHICH infrared limits the antenna covers. For the
3-parton (single-unresolved) antennae:

| antenna | radiators + unresolved | limits covered |
|---|---|---|
| `A30(q,g,qb)` | quark-antiquark, gluon unresolved | soft g; q∥g both sides |
| `D30(q,g,g)` (sub: `d30`) | quark-gluon | soft g; q∥g; g∥g |
| `E30(q,q',qb')` | quark line, secondary pair | g→q'q̄' collinear only (no soft) |
| `F30(g,g,g)` (sub: `f30`) | gluon-gluon | soft g; g∥g |
| `G30(g,q,qb)` | gluon line, secondary pair | g→qq̄ collinear only |

Verified in source: `FullA30FF` is the eikonal
`s12/s23 + s23/s12 + 2 s13 s123/(s12 s23)` with the middle argument
soft (comment: "i.e. for i2 soft" — the ARGUMENT ORDER encodes which
leg is unresolved).

The 4-parton `X40` letters follow the same radiator logic with two
unresolved partons (`A40` = qq̄ + gg, `D40` = qg + gg, `B40`/`C40` =
secondary quark pairs on a quark line, `F40`/`G40`/`H40` = gluonic /
multi-pair; `At40` etc. subleading colour). They cover the
double-unresolved limits: double soft, triple collinear,
soft+collinear, double collinear. Note the E30/G30 "no soft" rule is a
single-unresolved statement — a same-flavour qq̄ PAIR does have a
double-soft limit at X40 level (B40-type), provided deleting the pair
leaves a legal Born (reduced-Born rule, see run-spike-test). For the
exact per-letter content check `maple/notation.pdf` or the header of
the `src/X40` file — and when assigning a specific double limit, do
not guess: measure the antenna's pole graph
(probe-me-ir-structure pole scan) before placing it.

`X31` (one-loop 3-parton, `A31FF`, ...) cover the same single limits
as their X30 partner at one loop; they appear only in RV (`*T`) terms.

**Debugging use**: a spike test failing in limit L on channel C →
the suspect subtraction lines are those whose antenna covers L with
the unresolved parton(s) of L among its arguments (see
run-spike-test / write-subtraction).

## The NNLO scheme map (Fig. 3 of arXiv:1301.4693)

The whole NNLO subtraction is three layers of blocks; every unintegrated
block reappears integrated one layer up (the figure's arrows):

- **dσ^S (RR)**: S,a = `X30·M_{n+1}` (single-unresolved); S,b1 =
  `X40·M_n` (colour-connected double-unresolved); S,b2 =
  `−X30·X30·M_n` (iterated overlap removal); S,c = large-angle soft
  SS-difference blocks (almost-colour-connected pairs); S,d = disjoint
  `X30×X30` products (colour-unconnected pairs).
- **dσ^T (RV)**: T,a = `−J2^(1)·M_{n+1}` (cancels the RV poles); T,b =
  `X30·M^(1)` and `(X31 + X30·J21)·M_n` (RV's own unresolved limits);
  T,c = the integrated S,c.
- **dσ^U (VV)**: U,A = `−J2^(1)·M^(1)`; U,B = `−½ J2^(1)⊗J2^(1)·M`;
  U,C = `−J2^(2)·M` — pure integrated dipoles cancelling the two-loop
  poles.

**Colour connection of an RR pair decides its block** (the criterion,
hep-ph/0505111 §2.3): colour-connected (radiators shared between the two
unresolved partons) → one X40 (S,b1) plus two iterated X30×X30
counterterms (S,b2); almost-colour-connected (separated by a single hard
radiator) → two lines, one per strong ordering, carrying the SS soft
correction (S,c); colour-unconnected (disjoint dipoles) → plain product,
**unordered** pairs only (S,d). At subleading colour the classification
applies per colour structure — the tilde antennae implement the 1/nc
tower (full-colour dijets: arXiv:1310.3993).

**The signs are derived, not conventional.** `dσ^{S,a}` alone vanishes
in genuinely colour-connected double limits and yields **exactly twice**
the matrix element in the other three classes (the double limit needs
both the antenna's unresolved leg and one other momentum, and their
roles interchange). So S,b1 is `+` because S,a left a hole, and S,c/S,d
are `−` because S,a over-filled it. Full argument, plus the
**colour-neighbouring** sub-case hiding inside the colour-connected
class — where the X30×X30 products do NOT vanish and pole-sharing
under-counts — in write-subtraction's `references/block-counting.md`.
This classification is emitted mechanically by
`.claude/skills/write-subtraction/scripts/predict_blocks.py`.

**Provenance caveat on S,c.** hep-ph/0505111 has NO large-angle soft
term; its decomposition is exactly `S = a+b+c+d` with no eikonal object.
The SS blocks were introduced in the revised arXiv:0710.0346 because
*"the angular averaging is not sufficient to cancel the 1/ε poles"* in
some colour factors — i.e. the subtraction terms themselves introduce
spurious large-angle soft limits. arXiv:1301.4693 then folded them into
S,c, which is the grouping the repo uses. Practical consequence: an
uncancelled 1/ε where the angular average "should" have worked may be a
missing S,c/LAST block rather than a broken antenna.

**Integrated objects**: `SSset = {SFF, SIF, SFI, SII}` are MAPLE tokens,
and there is no 1:1 Fortran file per token — the tree has three
unintegrated soft functions (`src/X30/SS.f`, `SS1.f`, `SSII.f`) and
three integrated ones (`src/X30int/SSint.f`, `SSintIF.f`, `SSintFI.f`;
no `SSintII.f`). `SS1(j1,i3,j2,jpset,ipset)` is the eikonal
`2·s_{j1j2}/(s_{i3j1}·s_{i3j2})` with the radiators on the `jpset`
momenta and the soft leg on the `ipset` ones — mapped and unmapped legs
mixed in one call, which is the point (see write-subtraction's S,c
section). Documented in `doc/LAST/SSintIF.tex` ("LAST" = large-angle
soft term).

Every unintegrated antenna has an integrated counterpart in the virtual
layers: X30 ↔ `J21*` integrated dipoles (`J21QGFF`, `J21GQFI`, ...,
defined in `maple/form/common/J21.map`, Fortran in
`src/X30int/autoJ21*.f`), X40/X31 ↔ `J22*` (defined in
`maple/form/common/J22.map`); assembled into the `J2^(ℓ)` operators,
which absorb the mass factorisation for initial states and are related
to Catani's I operators.

**Mass factorisation lives INSIDE the J21/J22 definitions**, not in a
separate module (there is none in `driver/core/`). The live token names
are `gamma1qq(z1)`, `gamma1gg(z1)`, `gamma2gq(z1)`, ... plus `P0set` /
`P1set` — declared in `maple/notation.map` (`gamma1set`, `gamma2set`,
~lines 338–408). Grep for `gamma1`, not for `Gamma`/`MF`/`massfact`,
which do not exist. The absorption pattern is visible in
`maple/form/common/J21.map`: the kernel attaches to the dipole touching
that initial leg, with a factor matching the antenna's own coefficient,
and a **minus sign** for identity-changing dipoles; FF entries and
N_F-type dipoles carry no kernel at all:

```
J21QQIF(1,3)     = calqA30IF(1,3)      - gamma1qq(z1)
J21GQIF(1,3)     = calgD30gqIF(1,3)    - 1/2*gamma1gg(z1)
J21QQgtoqIF(1,3) = -1/2*calgA30IF(1,3) - Sgtoq*gamma1qg(z1)
J21QQFF(3,4)     = calA30FF(3,4)                      # FF: no kernel
```

Fortran side: `src/X30int/autoP0IF.f`, `autoP0FI.f`.

**There is no J22 Fortran.** No `src/X40int/`, no `src/X31int/`, no
`*40*int` symbols anywhere in `src/`. J22 exists only in the FORM layer
(`maple/form/common/J22.map` → `autoJ22.frm`, `autoJ22.h`, `doJ22`;
catalogue in `doc/J22/J22_catalog.tex`) and is consumed by the pole
check. Do not go looking for it under `src/`. This bookkeeping
— crossings and symmetry factors included — is exactly the Fig.-3
arrow structure that run-layer-check verifies per process; its failure
residues are printed in this language (leftover `calX30`/`J21` symbols
name the missing or mis-factored integrated antenna).

## Fortran level

Generated subtraction code calls the `Full<antenna>` wrappers with one
extra trailing argument — the kinematics-set index:
`FullA30FF(i5,i6,i7,7)` evaluates on `kin(7)` (7-parton phase space).
Unlike matrix elements, crossed antennae (`FF` vs `qA30IF` vs
`qqA30II`) are DISTINCT functions with distinct integrated
counterparts — not one-line aliases; never substitute one
configuration for another.

## Lookup recipe

1. Classify the limit or radiator pair → letter (table above);
   configuration from where the radiators sit (FF/IF/FI/II).
2. Confirm the token exists in `maple/notation.map` (`ant30set` /
   `ant31set` / `ant40set` — grep the exact name; makeproc and
   makefortRR key on exact names).
3. Implementation: `src/X30|X31|X40/<Full...>.f`; integrated partner:
   `src/X30int/autoX30<config>int.f` / `autoJ21<config>.f` (flat dir,
   suffix not subdir); symbolic `cal*` names:
   `test/layer_check/include/X30.inc` etc.
4. Usage precedent: grep the exact name across `maple/process/`,
   `src/process/` and `driver/process/` INCLUDING generated `auto*.f`
   — generated Fortran survives when `.map` sources do not, and a
   generated call site shows argument order, mapping, coefficient and
   the counterterm it pairs with (procedure in write-subtraction).
5. Before first use in a new term: `antenna_datasheet.py show <name>`
   (cached measured convention); if the antenna is not in the
   datasheet, `measure` it, or pole-scan directly
   (probe-me-ir-structure) — see the species-vs-convention rule at the
   top.

Related: subtraction-term structure (which lines carry which antenna)
→ write-subtraction; ME names → me-naming-convention.

## References

- hep-ph/0505111 — Gehrmann-De Ridder, Gehrmann, Glover: "Antenna
  Subtraction at NNLO" — defines the antenna functions and the
  colour-connection classes.
- 0711.4711 — e+e- event shapes at NNLO — first full-scale
  application (the epem processes in this repo).
- 1301.4693 — Currie, Glover, Wells: "Infrared Structure at NNLO" —
  the J2^(ℓ) integrated-dipole formulation; Fig. 3 is the scheme map
  above.
- 1310.3993 — full-colour NNLO gluonic dijets — subleading-colour /
  tilde-antenna structures in practice.
