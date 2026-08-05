---
name: add-process-to-driver
description: >
  Write the process-definition maple file driver/maple/<PROC>.map that
  declares the channel structure (LO, R, V, RR, RV, VV: matrix elements,
  particle content, colour factors) of an NNLOJET process. Use when the user
  asks to add a process to the driver, write a <PROC>.map, or set up the
  channel/matrix-element structure for a new or crossed process. Default
  scope is writing the .map file only; makeproc registration and generation
  run only if the user asks.
---

# Adding a process to the driver (the <PROC>.map file)

`driver/maple/<PROC>.map` is Maple source, `read` verbatim by
`driver/maple/makeproc`: `:=` assignment, `:` terminator, `#` comments,
`**` powers. The filename base MUST equal the process name used at
registration (`read cat(myname,`.map`)`).

Always start from the closest existing map (`epem.map` for e+e−, `ZJ.map`
/ `WpJ.map` for hadronic, `DIS.map` for DIS, `H2.map` for decays) and
study 2–3 neighbours before writing.

## Process name ≠ folder name — resolve FIRST

The runcard/driver process name and the maple/src folder names differ for
several core processes. Resolve every name against
`databases/proc_registry.yml` (next to this skill) before touching
anything. The classic traps:

| runcard/driver name | maple/process | src/process |
|---|---|---|
| `Z` (Drell-Yan) | `DY` | `DY` |
| `ZJ` | `Z` | `Z` |
| `GJ` | `G` | `G` |
| `HJ` | `H` | `H` |
| `Wp`/`Wm` | `W` | `W` |

There is NO process named `DY` — "the DY process" means the `Z` process
whose subtraction terms live in `maple/process/DY`. (The registry also
records both iprocess numbers per process; `test/layer_check/databases/
proc_folder_db.yml` is the layer_check tool's own copy of the folder
mapping — keep the two in sync.)

## Skeleton (epem.map as model)

```maple
Ofac:="ave*g2qcd*(4d0*pi*amz)**2*nc*(nc**2-1)/(2d0*nc)*8d0";
parset:={u,ub,d,db,q,qb,Q,Qb,g}:
psymset:=[g,q,qb]:

dress1:={{q=u,qb=ub,nqqb=nup},{q=d,qb=db,nqqb=ndown}}:
dress2:={{Q=u,Qb=ub,nQQb=nup},{Q=d,Qb=db,nQQb=ndown}}:

LO:=[
[B1g0Z,[qb,g,q,em,ep],1]
]:
R1:=[
[B2g0Zepem,[qb,g,g,q,em,ep],1],
[Bt2g0Zepem,[qb,gt,gt,q,em,ep],-1/nc**2],
[C0g0Zepem,[qb,Q,Qb,q,em,ep],1/nc],
[D0g0Zepem,[qb,q,qb,q,em,ep],-1/nc**2]
]:
V1:=[
[B1g1Zepem,[qb,g,q,em,ep],1],
[Bt1g1Zepem,[qb,g,q,em,ep],-1/nc**2],
[Bh1g1Zepem,[qb,g,q,em,ep],nf/nc]
]:
RR:=[ ... ]:   RV:=[ ... ]:   VV:=[ ... ]:

XX:=[LO,R1,V1,RR,RV,VV]:
```

- `Ofac` — Fortran string, emitted verbatim as `facB=` in
  `qcdnorm<PROC>.f`: spin/colour average `ave`, coupling powers, EW
  constants, overall nc factors. makeproc rewrites `g2qcd` →
  `(4d0*pi*alphas_scl)` and `amz`/`amzw` → `alpha_param()`.
- `parset` — near-vestigial (one use in the hadron-fragmentation block);
  copy `{u,ub,d,db,q,qb,Q,Qb,g}` unless a neighbour differs.
- `psymset` — final-state identical-particle symmetrisation: every
  species listed is summed over its permutations by makeproc and divided
  by n!. Standard `[g,q,qb]`. **Therefore NEVER hand-write 1/2 or 1/4
  symmetry factors into production colour factors — that double-counts.**
  (Explicit symmetry factors exist only as the 4th element of decay-map
  entries, cf. the notation header of `H2bb.map`.)
- `dress1/2/3` — flavour dressing, ONE list per independent quark line
  (`q/qb` → dress1, `Q/Qb` → dress2, `R/Rb` → dress3); makeproc takes the
  Cartesian product to generate concrete flavour channels; `nqqb=nup` is
  the generation-multiplicity counter. W-coupled lines dress
  flavour-changing (see `WpJ.map`, plus file-level `qFlavChange`).
- `XX:=[LO,R1,V1,RR,RV,VV]:` — copy VERBATIM. It is identical across
  essentially all NNLO maps; use `[]` for layers not implemented
  (`XX:=[LO,R1,V1,[],[],[]]:` for NLO-only, cf. `DIS_h.map`).

## Channel entries: [ME, particle content, colour factor]

**Particle content** — Maple list, all-outgoing convention, and its order
IS the Fortran argument order of the ME. Hard constraint: crossable
partons FIRST, colourless/decay particles LAST (makeproc permutes only
the leading block); in EW maps the photon comes before the leptons.
Symbols: `q/qb` first quark line, `Q/Qb` second, `R/Rb` third; `gt` =
photon-like (abelian) gluon of subleading-colour amplitudes; `em,ep` /
`lm,lp` / `nu` leptons; `ph1,ph2` Higgs-EFT gluon legs; `G` crossable
photon; explicit `u/ub/d/db` only where dressing cannot produce the
pairing (region a/b, see below).

**Finding the matrix elements** — nomenclature
`[s|f1|f2|Full]<Letter>[t|tt|h|hh|th]<n>g<l><PROCtag>[a|b|c][_OL]`:

- Letter = quark-line structure: `A` 0 pairs (pure gluon), `B` 1 pair,
  `C` 2 distinct pairs, `D` 2 identical, `E`/`F`/`G` 3 pairs
  (distinct / one identical / all identical).
- `t` = tilde, colour-subleading; `h` = hat, closed fermion loop (nf);
  combinations stack (`tt`, `th`, `hh`).
- `<n>` = number of gluons, `<l>` = number of loops. `s` prefix =
  symmetrised over the qqb pair; `x` infix = production×decay
  (`B0g0ZHepemxBy0g0H`).
- Given the Born ME: **V-layer MEs = same particle content, l+1 loops;
  R-layer MEs = one extra parton, same loops.** RR = +2 partons, RV =
  +1 parton +1 loop, VV = +2 loops.
- Check existence: grep `src/process/<DIR>/` for the amplitude files
  (bare MEs like `B1g0Z.f`, plus `auto<ME><S|T|U|SNLO|TNLO>.f`); crossing
  variants carry suffixes (`Zepem`, `ZDIS`, none). Colour factors do NOT
  change under crossing (epem/ZJ/DIS lists are identical) — only ME
  suffixes and lepton symbols do.

**Colour factors** — Maple expressions in `nc`, `nf`. Patterns correlate
with the name modifiers:

| ME type | typical factor |
|---|---|
| leading (B, colour in Ofac) | `1` |
| `t` (tilde) | `-1/nc**2` |
| `tt` | `+1/nc**4` |
| `h` (nf loop) | `nf/nc` |
| `hh` | `nf**2/nc**2` |
| `th`/`ht` | `-nf/nc**3` |
| C-type (2 distinct pairs) | `1/nc` |
| D-type (identical pairs) | `-1/nc**2` |

One order lower in alpha_s the `(nc**2-1)` may sit in the factors instead
of `Ofac` (compare `Z.map` vs `ZJ.map`). Higgs maps use the magic Wilson
symbols `C1,C2,C2h` (makeproc changes control flow on them — reserved).

**DISCLAIMER — always include when delivering the file**: colour factors
are not trivially guessable; the table above gives the conventional
pattern, but relative signs and `(nc**2-1)` placement are
convention-dependent per process family. Present the full channel lists
to the user and ask them to double-check every colour factor against a
neighbouring process before use.

## Gotchas

- **Region a/b** (W/VH processes, RR and beyond): 4-quark amplitudes with
  the W on the up-type (`a`) vs down-type (`b`) line must be written with
  explicit `u/ub/d/db` flavours, and the ME names must be added to
  makeproc's four hard-coded a/b name sets (~lines 1283, 1286, 1849,
  1859) — otherwise channels are SILENTLY mis-crossed.
- Any new identical-flavour ME (D/F/G-type) must be in makeproc's
  `identset` (~lines 156–174).
- **`imapprocess` mis-numbering is SILENT**: setting it to a number from
  the wrong table (makeproc vs iprocess.map), or leaving it stale while
  `mymapdir` points somewhere real, makes the generated check scripts
  read a *different process's* subtraction terms with no error. Live
  example in this tree: Z_EW carries `imapprocess:=50`, which now
  resolves to DYtest's iprocess.map entry. Cross-check against
  `databases/proc_registry.yml` and grep BOTH tables for the number
  before picking it.
- Optional 4th element of an entry = Maple set of options: `{OL_qcd=n}` /
  `{OL_ew=n}` (mandatory for `*_OL` amplitudes), `{EW_real=n}`,
  `{EW_virt=n}`, `{CC=n}` (per-amplitude qFlavChange).
- Decays: either `applydecay` (`H2.map`: read prod map, `decay:=...`,
  `XX:=applydecay(...)`) or `decaylist:=[H2bb]:` — the decaylist route
  only handles a SINGLE entry correctly (see
  `epemZHhad_prodXdecay_plan.md` at repo root).
- File-level flags (`flag_EW`, `qFlavChange`, `nulllist`, `select_FS`,
  ...) are reset by makeproc before each read — set only what differs
  from defaults.

## Registration and generation — ONLY IF THE USER ASKS

Default deliverable is the `.map` file. Then ASK the user whether to
register and run makeproc. If yes:

There are TWO independent process-numbering tables — do not confuse them:

- `driver/maple/makeproc` numbers (header comment + `elif(iprocess = N)`
  dispatch): select which process makeproc generates.
- `maple/iprocess.map` numbers (`getiprocess()`): a SEPARATE table giving
  the maple subtraction `directory`, `iprocessname` (selects
  `maple/FLAVlist<name>.map`), and the RR/RV/VV channel ranges.

The link between them is `imapprocess`: it is an index into
**`maple/iprocess.map`**, NOT into makeproc. It is consumed only by
`makecheck(...)` (`makeproc:3110`), which bakes
`maple makeRRcheck -Diprocess=<imapprocess>` lines into the generated
`autocheck*` scripts (`mymapdir` supplies their `read ../<dir>auto*.map`
paths). Worked example: Drell-Yan is `iprocess = 1` (`myname:="Z"`) in
makeproc but `imapprocess = 5` in iprocess.map — unrelated numbers.
Both numberings per process are recorded in `databases/proc_registry.yml`.

1. In `driver/maple/makeproc`: add the header comment line
   (`# iprocess = N:  "<PROC>":`) and the dispatch branch:

   ```maple
   elif(iprocess = N)then
     myname:="<PROC>":
     mynameflag:="HADRON":        # or "DIS" or "epem"
     mydir:="../process/<PROC>/":
     mymapdir:="":                # or ../../maple/process/<DIR>/ for checks
     imapprocess:=0:              # iprocess.map number; 0 = no check-script link
     MYJET:=`ecuts_vj`:           # process cuts routine, cf. neighbours
   ```

2. Set the loop bound `for iprocess from N to N do` (~line 183).
3. `cd driver/maple && maple makeproc` (needs the maple binary; run from
   that directory — relative paths). Output into `driver/process/<PROC>/`:
   `sig{B,R,V,RR,RV,VV}*.f`, `sig{SNLO,TNLO,S,T,U}*.f`, `lumi*.f`,
   `selectchannel*.f`, `qcdnorm*.f`, `array<PROC>.dat`, `*.config`,
   dummy subtraction stubs.
4. `autoAddFortran.py PROC.map` can automate step 1–2 plus the driver
   plumbing (`driver/Makefile`, `NNLOJET.f`, `initialiseproc.f`) — but it
   assumes hadron-collider defaults (`mynameflag:="HADRON"`,
   `MYJET:=ecuts_vj`), misses G-type ME names in its regex, and has an
   `imappprocess` typo in its template — review everything it writes.

Verification of the generated channels (`array<PROC>.dat`,
`selectchannel<PROC>.f`) is done entirely by the user — do not attempt
it. Repo is Mercurial: `hg add` new files, never commit unless asked.
