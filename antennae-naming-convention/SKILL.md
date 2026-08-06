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
  interpreting spike-test failures and layer-check residues. Read-only
  reference — modify nothing.
---

# NNLOJET antenna naming and limit coverage

Authoritative token dictionary: `maple/notation.map` (sets `ant30set`,
`ant31set`, `ant40set`, `SSset`; rendered as `maple/notation.pdf`).
Fortran implementations: `src/X30` (tree 3-parton), `src/X31` (one-loop
3-parton), `src/X40` (tree 4-parton), `src/X30int/{FF,IF,FI,II}`
(integrated antennae, split by configuration).

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
the `src/X40` file — do not guess when assigning a specific double
limit.

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
hep-ph/0505111): colour-connected (radiators shared between the two
unresolved partons) → one X40 (S,b1); almost-colour-connected
(separated by a single hard radiator) → iterated X30×X30 + SS soft
correction (S,b2+S,c); colour-unconnected (disjoint dipoles) → plain
product (S,d). At subleading colour the classification applies per
colour structure — the tilde antennae implement the 1/nc tower
(full-colour dijets: arXiv:1310.3993).

**Integrated objects**: `SSset = {SFF, SIF, SFI, SII}` — soft eikonal
functions (`src/X30/SS*.f`). Every unintegrated antenna has an
integrated counterpart in the virtual layers: X30 ↔ `J21*` integrated
dipoles (`J21QGFF`, `J21GQFI`, ..., `maple/form/common/J21.map`),
X40/X31 ↔ `J22*`; assembled into the `J2^(ℓ)` operators, which absorb
the mass factorisation for initial states and are related to Catani's
I operators. Implementations: `src/X30int/<config>/`. This bookkeeping
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
   `src/X30int/<config>/`; symbolic `cal*` names:
   `test/layer_check/include/X30.inc` etc.

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
