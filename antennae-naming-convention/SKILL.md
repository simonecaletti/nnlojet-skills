---
name: antennae-naming-convention
description: >
  Decode or construct NNLOJET antenna-function names (A30FF, d30FF, qA30IF,
  qqA30II, At40, A31FF, J21QGFF, ...) and map each antenna to the infrared
  limits it subtracts. Use for questions like "what is d30FF", "which
  antenna covers the 5||6 collinear limit", "what does the IF suffix mean",
  "what is J21", or when interpreting spike-test failures and layer-check
  residues. Read-only reference — modify nothing.
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
soft+collinear, double collinear. For the exact per-letter content
check `maple/notation.pdf` or the header of the `src/X40` file — do
not guess when assigning a specific double limit.

`X31` (one-loop 3-parton, `A31FF`, ...) cover the same single limits
as their X30 partner at one loop; they appear only in RV (`*T`) terms.

**Debugging use**: a spike test failing in limit L on channel C →
the suspect subtraction lines are those whose antenna covers L with
the unresolved parton(s) of L among its arguments (see
run-spike-test / write-subtraction).

## Soft functions and integrated counterparts

- `SSset = {SFF, SIF, SFI, SII}` — soft eikonal functions, summed
  inside RR/RV terms (`src/X30/SS*.f`).
- **Every unintegrated antenna has an integrated counterpart** added
  back in the virtual layers: X30 ↔ `J21*` integrated dipoles
  (`J21QGFF`, `J21GQFI`, ..., definitions in
  `maple/form/common/J21.map`), X40/X31 ↔ `J22*`. Integrated antennae
  live in `src/X30int/<config>/`. This X30↔J21 bookkeeping — with
  crossings and symmetry factors — is exactly what run-layer-check
  verifies; its failure residues are printed in this language
  (leftover `calX30`/`J21` symbols name the missing or mis-factored
  integrated antenna).

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
