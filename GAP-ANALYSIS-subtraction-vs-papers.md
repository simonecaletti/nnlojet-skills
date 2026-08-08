# Gap analysis: `nnlojet-skills` subtraction logic vs. the antenna-subtraction papers

Compared against:

- **hep-ph/0505111** — Gehrmann-De Ridder, Gehrmann, Glover, *Antenna Subtraction at NNLO*, JHEP 0509:056
- **1301.4693** — Currie, Glover, Wells, *Infrared Structure at NNLO Using Antenna Subtraction*, JHEP 04(2013)066

Skills reviewed: `write-subtraction`, `antennae-naming-convention`, `probe-me-ir-structure`,
`autogen-subtraction`, `run-spike-test`, `run-pole-check`, `run-layer-check`,
`write-spike-test`, `me-naming-convention`.
Cross-checked against the NNLOJET tree at `/home/scaletti/hepsoftware/NNLOJET`.

> **Retrieval caveat.** Neither paper could be fetched in full (the fetch tool truncates
> at ~70k characters). Verified verbatim: 0505111 §1–§2 complete (eqs. 2.1–2.38) + Table 1;
> 1301.4693 §1–§3.1.3 complete (eqs. 1.1–3.13) + Tables 4/5/6 + full ToC. Statements about
> 0505111 §§4–11 (incl. §8.5 *Angular terms*) and 1301.4693 §3.1.4 onward (S,c formula, S,d,
> T,a–T,c, U,A–U,C, the qq̄→gg worked example) are drawn from the companion papers
> **0710.0346** (e⁺e⁻→3j, which introduced the large-angle soft terms), **1310.3993**
> (gluonic dijets, explicit J₂⁽²⁾), **hep-ph/0502110** (gg antennae) and are flagged as such.
> The gaps listed below that depend on those sections are marked ⚠ and should be re-verified
> against a local PDF.

---

## Verdict in one paragraph

The skills are a **very good engineering guide to one quadrant of the problem**: constructing
and debugging the *unintegrated real* subtraction term (`*S.map`, `*SNLO.map`) for a given
channel. Within that quadrant they are in places *better* than the papers — the papers do not
tell you that `Full` composites must be split on all-final clusters, that the two writings of an
iterated counterterm are inequivalent, or that antenna slot conventions are undocumented and must
be measured. Those are real, hard-won operational facts that no paper contains.

The holes are of three kinds:

1. **A whole half of the scheme is missing.** The T (RV) and U (VV) `.map` files get three
   bullet lines and one bullet line respectively — yet the tree holds **928 `*T.map`** and
   **959 `*U.map`** files against 892 `*S.map`. Everything about *how to build* them is absent.
2. **Mass factorisation and the initial-state machinery are nearly absent** — which is the
   entire subject of 1301.4693.
3. **The derivational logic is replaced by "measure it".** That is pragmatic and often
   correct, but several of the papers' *counting arguments* answer in a line questions the
   skills currently answer with a build cycle, and generalise where a measurement does not.

---

## Tier 1 — structural holes

### 1.1 Mass factorisation is essentially undocumented

**Papers.** 1301.4693 eqs. (3.1)–(3.4) make the MF counterterms co-equal with RR/RV/VV:

```
dσ̂_NNLO = ∫_{n+2} dσ̂^RR + ∫_{n+1}( dσ̂^RV + dσ̂^{MF,1} ) + ∫_n ( dσ̂^VV + dσ̂^{MF,2} )   (3.1)

dσ̂^{MF,1} = − ∫(dz1/z1)(dz2/z2) (αs N/2π) C̄(ε) Γ⁽¹⁾_{ij;kl}(z1,z2) ( dσ̂^R_NLO − dσ̂^S_NLO )   (3.2)

dσ̂^{MF,2} = − ∫(dz1/z1)(dz2/z2) { (αs N/2π)² C̄(ε)² Γ⁽²⁾_{ij;kl} dσ̂_LO
                                 + (αs N/2π) C̄(ε) Γ⁽¹⁾_{ij;kl} ( dσ̂^V_NLO − dσ̂^T_NLO ) }   (3.3)

Γ⁽²⁾_{ij;kl}(z1,z2) = δ(1−z2)δ_lj Γ⁽²⁾_ki(z1) + δ(1−z1)δ_ki Γ⁽²⁾_lj(z2)
                     + Γ⁽¹⁾_ki(z1) Γ⁽¹⁾_lj(z2)                                          (3.4)
```

The paper's organising trick is that the kernels are **not** a separate additive block: they are
absorbed pairwise into the integrated dipoles (its Tables 5 and 6), so that only genuine
final-state poles survive in the virtual subtraction term. The absorption rule is specific —
the kernel attaches to the dipole touching that initial leg, with multiplicity ½ or 1 matching
the antenna's symmetry factor, `δ(1−x)` on the spectator fraction, a **minus sign** for
identity-changing antennae, and **no kernel at all** for N_F-type dipoles.

**Skills.** One subordinate clause in `antennae-naming-convention`: *"assembled into the J₂⁽ℓ⁾
operators, which absorb the mass factorisation for initial states"*. Nothing else. No mention
that:

- the live maple names are `gamma1qq(z1)`, `gamma1gg(z1)`, `gamma2gq(z1)`, plus `P0set`/`P1set`
  (declared in `maple/notation.map` ~lines 338–408) — **not** `Gamma`, `MF`, `massfact`;
- they live inside `maple/form/common/J21.map` and `maple/form/common/J22.map`;
- there is **no** mass-factorisation module in `driver/core/` — this is entirely a maple/FORM
  layer concept, with Fortran only in `src/X30int/autoP0IF.f`, `autoP0FI.f`, `autoJ21{FF,FI,IF,II}.f`.

In-tree confirmation of the sign/factor pattern the papers describe:

```
J21QQIF(1,3)      = calqA30IF(1,3)      - gamma1qq(z1)
J21GQIF(1,3)      = calgD30gqIF(1,3)    - 1/2*gamma1gg(z1)
J21QQgtoqIF(1,3)  = -1/2*calgA30IF(1,3) - Sgtoq*gamma1qg(z1)
J21GGqtogIF(1,3)  = -calqpG30IF(1,3)    - Sqtog*gamma1gq(z1)
J21QQFF(3,4)      = calA30FF(3,4)                          # FF: no kernel
```

**Impact.** Anyone asked to add or debug an initial-state process cannot find the MF terms from
the skills at all — not even the vocabulary to grep for. `run-layer-check` and `run-pole-check`
will *report* the failure, and `write-subtraction` will send them to the `.map` file, but the
skill set contains no statement that these objects exist.

### 1.2 Building the T (RV) and U (VV) terms has no procedure

**Skills.** `write-subtraction` has one paragraph on RV inside step 4 (`T,a = −J21×M_n`,
`T,b = X30×M^{1-loop}` + `(X31+X30·J21)`, `T,c` = integrated S,c) and a three-bullet
"structural patterns" list. `antennae-naming-convention` gives one line each for `U,A`/`U,B`/`U,C`.
`autogen-subtraction` covers the *mechanics* of `makefortRV`/`makeformVV` but explicitly not the
content. There is no analogue of the RR block-by-block construction procedure, no
completeness check, no failure-mode catalogue, and no worked structure.

**Reality in the tree.** The `*U.map` files are structurally *different from anything the skills
describe*. `maple/process/H/qgA2g2HU.map` contains **three** top-level objects, not one:

```
FN  := qgA2g2HU(1,2,H):
VV  := ... ;  VV := -VV*Dd([1-z1])*Dd([1-z2]):     # two-loop pole structure, c_i labels
XX  := ... a1 ... a5 ;                              # the subtraction term
VVD := 0:
```

with tokens the skills never mention: `calqgG40II`, `calqgG40nadjII`, `calqgG31II` (integrated
X40/X31), `gamma1gq(z1)`, `gamma2gq(z1)` (MF kernels), `b0*e^(-1)` (renormalisation), `Sqtog`/
`Sgtoq` (flavour-change tags), `QQ(s12)` (scale), and arguments in `z1,z2` which `makeformVV`
rewrites to `x1,x2` at line 154. `makeformVV` maps `calXXX → XXXint` (lines 236–254).

The `*T.map` files likewise vary in ways the skill does not flag: `Z/qgBt2g1ZT.map` uses
`XX:=expand(...)` and `J21…` symbols; `H/A4g1HXT.map` uses neither, writing the integrated
antennae longhand as `calF30FF`, `calgF30IF`, `calggF30II`, and carries a `#goto/#from/#combine
with/#cancel with` bookkeeping convention (tracking cancellation against `A5g0HXS` and
`A3g2HXU`) that is nowhere documented.

**Impact.** This is the single largest hole. `*U.map` is the *biggest* term family in the repo
and has zero construction guidance.

### 1.3 The `b0` / renormalisation-scale block is missing

**Papers.** 0505111 eq. (2.34): `X¹₃` is renormalised at `s_ijk` while the one-loop reduced ME is
at `μ²`, so one must substitute

```
X¹_{ijk} → X¹_{ijk} + (β0/ε)·((4π)^ε e^{−εγ}/8π²)·X⁰_{ijk}·[ (s_ijk)^{−ε} − (μ²)^{−ε} ]
```

and *"these terms will in general be kept apart in the construction of the colour-ordered
subtraction terms, since they all share a common colour structure β0."* 1301.4693's ToC names
**T,b3 = "One-loop renormalisation subtraction term"** as a block in its own right.

**Skills.** The RV pattern list has `T,a / T,b / T,c` — three blocks, no b3, no β0 anywhere in
the skill set.

**In the tree** it is unmissable:

```
J22QQFF(3,4) = calA40FF + calA31FF - 1/2*calA30FF*calA30FF + b0/e*calA30FF
J22GGFF(3,4) = 1/4*calF40FF + 1/3*calF31FF - 1/9*calF30FF*calF30FF + b0/e/3*calF30FF
```

### 1.4 Convolution ⊗ vs. product ×

**Papers.** J₂⁽²⁾ contains a genuine **convolution** over momentum fractions, `−[X⁰₃ ⊗ X⁰₃]`,
for the IF/II configurations (1310.3993 eqs. 2.20–2.22; 1211.2710 devotes an appendix to these
convolutions). Only in FF does it degenerate to a product.

Related and more subtle — 0505111, immediately after eq. (2.25): in the iterated (almost-colour-
unconnected) case *"some care has to be taken in the second dipole integral … which will pick up
**ε-dependent factors from the first integral** (both integrals are fully independent only in
four dimensions) … the analytic integration will **not** yield the product of two independent
integrated NLO antenna functions."*

**Skills.** No `⊗` anywhere. `antennae-naming-convention` states the layer-check arrow structure
as *"every unintegrated antenna has an integrated counterpart in the virtual layers"* with no
caveat — which is exactly the naive reading the paper warns against.

---

## Tier 2 — conceptual holes inside the RR logic that *is* covered

### 2.1 The "factor of two" counting argument — why S,c and S,d exist at all

**Paper (0505111 §2.3.1, verbatim).** *"(2.16) yields **twice** the (m+2)-parton matrix element
in all colour-unconnected, almost colour-unconnected and colour-neighbouring double unresolved
limits, while **vanishing** in all genuinely colour-connected double unresolved limits."*
The factor of two is derived: the double-unresolved limit needs both `p_j` (the antenna's
unresolved leg) and one other `p_o` to go unresolved, and *the roles of `p_j` and `p_o` can be
interchanged*, giving two identical terms.

That single sentence explains, without measurement, (i) why S,c and S,d carry an **overall minus
sign**, (ii) why there is no S,d-like term for colour-connected pairs, and (iii) why S,a alone
over-subtracts by exactly a factor 2 in those regions.

**Skills.** State the blocks and their signs as facts ("typically negative sign") and route every
"is this right?" question to a residue fit. The counting argument is absent. This is the clearest
case where a paragraph of theory replaces a build cycle.

### 2.2 The colour-neighbouring sub-class is missing from the classification

**Paper.** 0505111 identifies a **fourth** configuration inside the colour-connected class:
two neighbouring pairs going collinear independently, one pair inside the antenna, the other
formed by the remaining antenna momentum and its colour-connected neighbour. Here `X⁰₄`
approximates the ME correctly *and* the `X⁰₃X⁰₃` products are non-vanishing, each equal to the
double-unresolved limit — they are there to cancel the doubled terms of S,a, not (only) the
spurious singles of `X⁰₄`. *"Each configuration of this type is contained precisely twice in (2.16)."*

**Skills.** A clean trichotomy: colour-connected / almost-colour-connected / unconnected. The
S,b2 construction rule is *"the set of iterated counterterms is indexed by (S,a line, X40 line)
pairs that share a singular invariant"* — an operational proxy that is exactly where the
colour-neighbouring case is most likely to mis-count, because there the counterterm's job is
dual (cancel `X⁰₄`'s spurious single **and** cancel S,a's doubled double) and the pole-sharing
criterion sees only the first job.

### 2.3 The sub-antenna `x⁰` in S,c — the half-eikonal rule

**Paper (0505111 eq. 2.24 + text).** *"`x⁰_{mlK}` denotes a **sub-antenna**, containing only the
collinear limit of m with l, but **not** the collinear limit of l with K; in the soft limit of l,
this sub-antenna yields **half** the soft eikonal factor."* The reason is stated: `K` is a
*mapped* (composite hard) momentum, so a collinear `l∥K` singularity would be spurious — the ME
has no such limit.

**Skills.** `write-subtraction` says *"use sub-antennae (`d30`, `f30`) where a full antenna would
double-count limits shared between overlapping clusters"* and gives the Full-vs-split rule for
crossing, with the instruction to measure each half's pole graph and check the union. The
**half-eikonal** property and the **"a mapped momentum must not be a collinear partner"**
criterion are both absent — and that criterion is precisely what tells you *which* half to use,
which the skill currently answers with two measurement runs.

### 2.4 The eikonal's hard legs need not be the antenna's radiators

**Paper (0710.0346 ⚠).** *"Those soft factors are associated with an antenna phase space mapping
(i,j,k)→(I,K). **The hard momenta a, c do not need to be equal to the hard momenta i, k in the
antenna phase space — they can be arbitrary on-shell momenta.**"* With `S_abc = 2 s_ac/(s_ab s_bc)`.

**Skills.** `write-subtraction`'s S,c section is excellent on the *Fortran wiring* (radiators
`j1,j2` on the MAPPED set, soft leg `i3` on the ORIGINAL set, `s{ipset}on{jpset}` cross-set
commons filled by `fillson*`) — better than any paper. But it never states this as a **physics
freedom**, and that freedom is what generates the six-term difference structure
`S_{mapped,mapped} − S_{half-mapped} − S_{unmapped} + …` in the first place. The skill says
"construct it by measuring the difference of soft residues"; the paper says what the difference
*is made of*.

Confirmed in-tree: `src/X30/SS1.f` computes `abs(2*s12/s13/s23)` with `(jpset, ipset)` selecting
which momentum sets the invariants come from — exactly the mixed-set eikonal.

### 2.5 The origin of large-angle soft terms is mis-attributed ⚠

**Paper history.** hep-ph/0505111 has **no** large-angle soft term — its decomposition is exactly
`S = a+b+c+d` (eq. 2.28) with no eikonal object anywhere. They were introduced in the **revised**
0710.0346, for a specific reason: *"in the N² and N⁰ colour factor, the **angular averaging is not
sufficient** to cancel the 1/ε poles in the four-parton one-loop subtraction terms"* — i.e. the
subtraction terms themselves introduce spurious large-angle soft limits.

**Skills.** `antennae-naming-convention` attributes the colour-connection criterion to
hep-ph/0505111 and lists S,c as *"large-angle soft SS-difference blocks"* in the same breath —
conflating the 2005 four-fold decomposition with the later five-fold one, and severing the causal
link between the angular-average failure and the existence of the SS blocks.

1301.4693 does fold large-angle soft into S,c (its §3.1 lists S,c as *"two almost colour-connected
unresolved partons … **and including large angle soft radiation**"*), so the skills' grouping
matches the 2013 paper — but the *reason* is still missing, and Appendix B of 1301.4693 is
devoted to the integrated IF large-angle soft term as a distinct ingredient.

### 2.6 Angular correlations are treated as a test artefact, not a design constraint

**Paper (0710.0346 §3.4 ⚠).** *"The angular terms … average to zero after integration over the
antenna phase space. To ensure numerical stability and reliability, **this average has to take
place within each phase space mapping**."* The authors state they verified this explicitly for
the decompositions of `E⁰₄` and `D⁰₄`.

**Skills.** `run-spike-test`: *"ratio oscillates around 1 without narrowing in a g→gg or g→qq̄
collinear limit → azimuthal-rotation issue, not necessarily a wrong .map."* `probe-me-ir-structure`
and `write-spike-test` correctly apply `rotp<n>` π/2 averaging. All of this is downstream —
it is about *measuring*.

Missing: the **design constraint**. "The average must happen within a single phase-space mapping"
forbids splitting a collinear limit across two mappings unless each averages correctly alone.
That rule belongs in `write-subtraction`'s *"Splitting a `Full` composite on all-final clusters"*
section, which currently says only "confirm the union covers every limit" — union of *poles* is
not the same test as *per-mapping angular closure*.

Also missing: that angular averaging is **not always sufficient**, and that its failure is the
origin of the large-angle soft blocks (see 2.5).

### 2.7 Normalisation factors and the `s_ik` vs `s_ijk` argument

**Paper (0505111).** Integrated antennae carry `(8π²(4π)^{−ε} e^{εγ})^k`, **one power per
unresolved parton** — squared for `X⁰₄`, single for `X¹₃` (eqs. 2.11 / 2.23 / 2.35). And in
`dσ^{VS,1,a}` (eq. 2.31) the integrated antenna is evaluated at the **two-particle invariant**
`X⁰_{ijk}(s_ik)`, not `s_ijk`, because the 3→2 map identifies `p_i,p_k` with `p̃_I,p̃_K`; the
paper warns the resulting integrals *"differ from the standard tree-level three-parton antenna
integrals by normalisation factors."*

**Skills.** Neither appears. This is precisely the class of error `run-pole-check` diagnoses as
*"ratio a constant ≠ 1 = wrong factor or crossing in the integrated-dipole assembly"* — the skill
knows the symptom but not the cause.

### 2.8 Antenna provenance explains the split rule the skills present as an empirical trap

**Papers.** Antennae are normalised colour-ordered matrix elements of specific parent processes:
`γ*→qq̄+partons` (quark–antiquark), `χ̃→g̃g+partons` via the Haber–Wyler effective Lagrangian
(quark–gluon), `H→gg+partons` via `L = −(λ/4)H F²` (gluon–gluon). hep-ph/0502110 spells out the
consequence: *"the matrix element has to be split into three individual antenna configurations.
**Each individual antenna configuration contains only one soft limit.** Each collinear g→gg is
**split between the two antenna configurations** appropriate to the two final-state gluons."*

That sentence *is* the derivation of the skills' Full-composite splitting rule, and it also
explains the numerical coefficients visible in the repo (`1/3*calF30FF`, `1/9*calF30FF²`,
`1/4*calF40FF` vs. `1*calA30FF`) as symmetry/averaging factors of the parent process.

**Skills.** Present the split rule as *"the central trap when deriving a term by crossing"* with
the instruction to measure. Correct, but non-generalising: an agent that knows *why* can predict
the split for an antenna it has never seen; one that only knows the trap must measure every time.

---

## Tier 3 — smaller items and factual drift

| # | Item | Papers say | Skills say |
|---|---|---|---|
| 3.1 | Identical-particle factors | `1/S_n` appears explicitly in **every** master formula (0505111 eqs. 2.1, 2.15–2.17, 2.24, 2.26, 2.31, 2.33, 2.36), attached to the **full** multiplicity of the channel, not the antenna | mentioned only as "±1/2 symmetry factors where clusters are symmetric" in S,b2 |
| 3.2 | S,d sum restriction | *"such that **no product of two antenna configurations appears twice**"* — an unordered-pair sum | "one line per pair of DISJOINT clusters that are each singular" — the no-double-count rule is implicit at best |
| 3.3 | `dσ^{VS,1,c}` sum restriction | `p_i,p_k` may coincide with `p_n,p_p,p̃_N,p̃_P` but **not** with the unresolved `p_o` (eq. 2.36) | absent |
| 3.4 | Mapping requirements | Explicit 4-item list: momentum conservation, on-shellness, correct reduction in exact limits, **no spurious singularities**; and separate maps for analytic integration (tripole, eq. 2.20) vs. numerical implementation | cluster rule only (`getpmapIK`) — correct and verified, but the *requirements* the map satisfies are not stated |
| 3.5 | `∫dσ^S` splits in two | 1301.4693 eq. (3.8): the same `dσ^S` splits into an `∫₁` piece landing in RV and an `∫₂` piece landing in VV | `run-layer-check` describes "every block reappears one layer up" — the two-way split is not stated |
| 3.6 | No generic formula | 0505111: *"the precise nature of cancellations … **differs considerably among the different colour structures. Therefore, no generic formula can be stated.**"* | skills imply a uniform arrow structure |
| 3.7 | `colflag` | — | one parenthetical ("groups several reduced MEs under one antenna"). **474 `.map` files set it.** `makefortRR` (lines 327–400) does three things with it: skips `expand()`, collects a **list** of MEs, and uses **only `matM0[1]`** for the flavour mapping. That last is a live trap: wrong first element ⇒ silently wrong `set_flav_perm`. |

### Factual drift to fix (verified against the tree)

- `antennae-naming-convention` and `write-subtraction`: *"Implementations: `src/X30int/<config>/`"*
  and *"`FF`/`IF`/`FI`/`II` subdirs of `X30int`"`. **`src/X30int/` is flat** — 14 files,
  configuration is a filename suffix (`autoX30FFint.f`, `autoJ21IF.f`, …). The subdir claim also
  appears in the repo's `CLAUDE.md`.
- `antennae-naming-convention`: *"`SSset = {SFF, SIF, SFI, SII}` — soft eikonal functions
  (`src/X30/SS*.f`)"* implies one file per token. The tree has **three** files:
  `SS.f`, `SS1.f`, `SSII.f` (no `SSFF.f`/`SSIF.f`/`SSFI.f`); integrated:
  `SSint.f`, `SSintIF.f`, `SSintFI.f` (no `SSintII.f`).
- `antennae-naming-convention`: *"X40/X31 ↔ `J22*`; … Implementations: `src/X30int/<config>/`"*.
  **There is no J22 Fortran at all** — no `src/X40int/`, no `src/X31int/`, no `*40*int` symbols
  in `src/`. J22 exists only in the FORM layer (`maple/form/common/J22.map`, `autoJ22.frm`,
  `autoJ22.h`, `doJ22`) and is consumed by the pole check. `run-pole-check` gets this right
  ("the `J21`/`J22` lines of the `*T`/`*TNLO`/`*U` maple files"); the naming skill does not.
- Undocumented in-tree resources the skills should point at:
  `doc/J22/J22_catalog.tex`, `doc/FF-antennae-decomposition/decompG.tex`,
  `doc/LAST/SSintIF.tex` (LAST = large-angle soft term), `doc/Flavour Changing/Flavour Changing J22.tex`,
  `doc/makeformVV/makeformVV.tex`, `doc/process/*/texfiles/subtractionTerms/{LO,NLO,NNLO}.tex`,
  `doc/*/MFterms/VVMF*.txt`. Only `doc/process/VFH/texfiles/spikesAndRotation.tex` is currently
  referenced (and it does exist).

---

## What the skills have that the papers do not

Worth recording, so none of it is lost in a revision:

- Antenna **slot conventions** are undocumented in source and must be measured; `antenna_slots.py`
  finds non-ascending Fortran declarations that make a positional `.map` call silently wrong.
- The **two writings of an iterated counterterm are numerically different**, the correct one is
  the one regular in the single-unresolved limits the other antenna already reproduces, and the
  wrong one gives O(1) sign-flipped ratios in exactly those modes while double-unresolved modes
  still look fine.
- The **reduced-Born rule** for classifying genuine vs. dead limits (240/240 validated), including
  the traps: composites not decomposed pairwise, all-components-genuine insufficient, colour
  adjacency irrelevant because check programs sum orderings.
- **S,a + S,b1 + S,b2 is the minimal testable unit**; partial builds fail with specific
  signatures (S,a alone → ~0.4 on soft-gluon; S,a+iterated → negative O(1)).
- Per-pole (not per-X40) consistency of iterated coefficients.
- `FLAVlist<proc>.map` as the authoring-time reduced-ME lookup, and that it is itself generated.
- The whole measurement stack: pole scan, dipole fit, residue fit, per-line `WTDBG` attribution,
  block composer, `regen_rebuild.sh`.
- Spike-test epistemics: plateau not depth; median not max; sub-singular modes give false passes;
  dead modes give meaningless ratios; the untouched-sibling control channel.

---

## Suggested remediation, in priority order

1. **New skill `write-integrated-subtraction`** (or split `write-subtraction` into `-real` /
   `-virtual`): the `*T.map` and `*U.map` construction procedure. Must cover the three-object
   `VV`/`XX`/`VVD` layout, `cal*` naming and the `makeformVV` `cal→int` rewrite, `z1/z2` vs
   `x1/x2`, `Sqtog`/`Sgtoq`, `QQ(s)`, and the `#goto/#combine with/#cancel with` bookkeeping.
2. **New skill or major section on mass factorisation**: `gamma1*`/`gamma2*`/`P0set`/`P1set`,
   `J21.map`/`J22.map` as the place they live, the absorption pattern (½ factors, `δ(1−x)` on the
   spectator, minus for identity-changing, none for N_F), the `b0/ε` renormalisation term, and
   `⊗` vs `×`.
3. **Add the derivations to `write-subtraction` step 4**: the factor-of-two counting argument, the
   colour-neighbouring sub-class, the half-eikonal sub-antenna rule, the "mapped momenta must not
   be collinear partners" criterion, and the arbitrary-hard-legs property of the eikonal. Keep the
   measurement procedures — the theory tells you what to expect, the measurement confirms it.
4. **Promote the angular-average constraint** from `run-spike-test` (diagnostic) into
   `write-subtraction`'s Full-splitting section (design rule), and state the
   angular-average-insufficiency → large-angle-soft causal chain.
5. **Fix the factual drift** in §3 above (flat `X30int`, three `SS*.f` files, no J22 Fortran) and
   add pointers to `doc/J22/`, `doc/LAST/`, `doc/FF-antennae-decomposition/`.
6. **Expand the `colflag` note** in `write-subtraction` to include the `matM0[1]` flavour-mapping
   trap.
