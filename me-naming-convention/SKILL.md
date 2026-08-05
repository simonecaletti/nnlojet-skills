---
name: me-naming-convention
description: >
  Decode or construct NNLOJET matrix-element names (B1g0Z, Bt2g1Zepem,
  sC0g1ZDIS, ...) and explain how one matrix element serves several
  processes via crossing. Use for questions like "what is Bt2g1Z", "which
  ME do I need for qqb to Z+2g at one loop", "can epem reuse the Z matrix
  elements", "why does B2g0Zepem exist next to B2g0Z", or whenever another
  skill needs to identify or look up a matrix element. Read-only
  reference — modify nothing.
---

# NNLOJET matrix-element naming and crossing

## The name grammar

```
[s|f1|f2|Full] <Letter> [t|tt|ttt|h|hh|th|ht] <n> g <l> <PROCtag> [a|b|c] [_OL]
```

- **Letter — quark-line structure**: `A` 0 quark pairs (pure gluon /
  gluon-fusion Higgs), `B` 1 pair, `C` 2 pairs distinct flavour, `D` 2
  pairs identical, `E` 3 pairs distinct, `F` 3 pairs one identical, `G`
  3 pairs all identical. Identical-flavour types (D, F, G) carry
  interference terms and live in makeproc's `identset`.
- **Modifiers**: `t` = tilde, colour-subleading (one power of 1/nc²);
  stacked `tt`/`ttt` = further suppressed; `h` = hat, closed fermion
  loop (one power of nf); `hh` = nf²; `th`/`ht` mixed. These correlate
  1:1 with channel colour factors (`Bt → -1/nc**2`, `Bh → nf/nc`, ...;
  table in add-process-to-driver).
- **`<n>` = number of gluons, `<l>` = number of loops**: `B1g0Z` = 1
  quark pair, 1 gluon, tree; `B1g2Z` = 2 loops; `A2g1H` = gg→H at 1 loop.
- **Prefixes**: `s` = symmetrised over the qqb pair
  (`sB1g0Z = 0.5*(B1g0Z(i1,i2,i3,..) + B1g0Z(i3,i2,i1,..))`), `Full` =
  colour-dressed wrapper, `f1`/`f2` = wrapper inheriting the flavour of
  the first/second quark line.
- **`<PROCtag>`**: the process family (`Z`, `W`, `H`, `GG`, none for
  pure jets) plus a crossing suffix (`Zepem`, `ZDIS`, `A` for photon) —
  see crossing below. `x` infix = production×decay product
  (`B0g0ZHepemxBy0g0H`); `y` marks decay-side MEs.
- **Trailing `a`/`b`** = W-emission region (up-type vs down-type line),
  `c` = further variant; `X`/`Y` = distinct colour orderings (jet maps).
  `_OL` = OpenLoops amplitude.

Perturbative bookkeeping from a known Born: V = same particle content,
l+1 loops; R = +1 parton, same loops; RR = +2 partons; RV = +1 parton
+1 loop; VV = +2 loops.

## Crossing: one implementation, many processes

MEs are written in the **all-outgoing convention**. A process fixes
which legs are crossed into the initial state; the crossing is entirely
in the momentum arguments — the SAME Fortran function serves every
crossing. Visible in generated code (`selectchannelepem.f`):

```fortran
chan_str = 'ep em to db g d'
me_str   = 'B1g0Z(3,4,5,1,2)'    ! e+e- annihilation: leptons in slots 1,2
```

versus the hadronic `Z` process calling `B1g0Z` with partons crossed
in. Argument position = leg in the all-outgoing amplitude; crossed legs
are charge-conjugated (makeproc's `con` substitution does this when
generating channels).

**The crossing suffixes are trivial aliases, not reimplementations.**
`src/process/epem/EPEMwrapper.f` and `src/process/DIS/DISwrapper.f`
consist of one-line wrappers, including loop amplitudes:

```fortran
      B2g0Zepem = B2g0Z(i1,i2,i3,i4,i5,i6)
      B1g2ZDIS  = B1g2Z(i1,i2,i3,i4,i5,rs2)
      B3g0A     = B3g0Z(i1,i2,i3,i4,i5,i6,i7)   ! photon = Z tag reused
```

They exist so each process family has distinct function names in the
flat link namespace (and for makeproc's bookkeeping); any
crossing-dependent analytic continuation or sign convention is handled
inside the shared implementation or absorbed in flags (e.g. the `idis`
flag compensating a sign in the DIS subtraction). Consequences:

- **Colour factors are invariant under crossing** — the epem/ZJ/DIS
  channel lists carry identical factors; only ME suffixes and lepton
  symbols differ.
- Before assuming an ME exists for your crossing, check the wrapper
  file (`<PROC>wrapper.f` / `EPEMwrapper.f` / `DISwrapper.f`) and grep
  `src/process/<DIR>/` — some variants are wrapped, some only exist
  symmetrised (`sB2g1ZDIS` with no bare `B2g1ZDIS`).
- Never treat suffixed and unsuffixed names as interchangeable in maple
  files or driver code — makeproc keys on the exact name.

## Looking up an ME

1. Resolve the process family folder via
   `add-process-to-driver/databases/proc_registry.yml` (no process is
   named `DY` — that is `Z`'s src/maple folder).
2. Construct the candidate name from the grammar: Letter from the
   quark-pair structure, count gluons and loops, append the family tag.
3. Grep `src/process/<DIR>/` for the bare amplitude file (`B1g0Z.f`),
   the wrapper (`*wrapper.f`), and the subtraction users
   (`auto<ME><S|T|U|SNLO|TNLO>.f`).
4. Where each family's amplitudes live and which processes share them:
   the registry's `src` column (e.g. `Z`, `ZJ`, `epem*`, `DIS*` all
   draw on `src/process/DY` + `src/process/Z`).

Related: channel colour-factor table and `identset`/region-a/b traps →
add-process-to-driver; reduced MEs inside subtraction terms (bracketed
mapped momenta) → write-subtraction.
