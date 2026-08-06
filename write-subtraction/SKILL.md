---
name: write-subtraction
description: >
  Write, edit, or fix an antenna-subtraction term at the maple level
  (maple/process/<DIR>/*.map) in NNLOJET. Use whenever the user asks to write,
  modify, or debug a subtraction term, fix a channel that fails a spike test,
  or add missing infrared limits to a subtraction. This skill covers ONLY the
  maple .map file; generating Fortran and registering it is the
  autogen-subtraction skill, and validation is the run-spike-test skill.
---

# Writing NNLOJET subtraction terms (maple level)

A subtraction term must mimic every infrared (unresolved) limit of the
corresponding channel of the full matrix element: each single-unresolved
limit for R, plus double-unresolved and one-loop×single limits for RR/RV.
A spike-test failure in a specific limit means the term line(s) covering
that limit are wrong or missing.

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
  cluster, not from the JET arguments (`Ct1g0ZepemS.map` line a2: the
  JET34 order differs from the emitted set_map list). Do not agonise
  over JET ordering.
- **Which mapped momenta a cluster produces**: a 3-parton antenna
  cluster (a,b,c) yields `[a,b]` and `[b,c]`; a 4-parton cluster
  (a,b,c,d) yields `[a,b,c]` and `[d,c,b]`. Bracket contents are
  order-insensitive up to reversal (`[k,i,j]` ≡ `[j,i,k]`). Use these to
  write the reduced-ME and JET arguments.
- **`aN`** is a sequential line label: `a1, a2, a3, ...` with NO gaps (the
  checker aborts on gaps). Each `aN` becomes one `wt(N)` slot and one
  independent mapping/jet block in the generated Fortran — numbering is
  load-bearing; renumber consistently when inserting/removing lines.

Optional: `colflag:=true:` (groups several reduced MEs under one antenna),
`XX:=expand( ... ):` in some RV files.

## Building the lines: from the ME's infrared limits

The line list is DERIVED from the full ME, not invented:

1. **Enumerate the limits of the full ME.** Single-unresolved (R/SNLO):
   every gluon soft, plus each collinear splitting (q∥g, g∥g, g→qq̄)
   between colour-connected partons. Colour connection is the organising
   principle: for colour-ordered B-type amplitudes it is adjacency in
   the argument list; for multi-quark-line (C/D/E/F/G) and
   subleading-colour structures read the connections off the analogous
   existing `.map` — do not guess. RR adds the double-unresolved set:
   double soft, triple collinear, soft+collinear, double collinear
   between disjoint pairs.
2. **One antenna per limit CLUSTER, not per limit.** Pick the antenna
   whose two hard radiators are the neighbours of the unresolved
   parton(s) and whose species match (antennae-naming-convention);
   arguments are (hard, unresolved, hard) — middle leg unresolved. One
   `A30FF(a,b,c)` line covers b-soft AND both a∥b, b∥c collinear at
   once; do not write separate lines per limit. Use sub-antennae
   (`d30`, `f30`) where a full antenna would double-count limits shared
   between overlapping clusters.
3. **Build the reduced ME** from the cluster rules above: remove the
   unresolved parton, substitute the mapped momenta for the radiators;
   the reduced ME is the (n−1)-parton amplitude of the resulting
   flavour content (find its name via me-naming-convention). JET
   arguments = the reduced final-state momenta (order free).
4. **RR assembly order** (worked example: `Ct1g0ZepemS.map`):
   (a) `X30 × M_{n-1}` for each single-unresolved cluster (its a1–a8);
   (b) `X40 × M_{n-2}` for each colour-connected double-unresolved
   cluster (a9, `B40`);
   (c) MINUS the iterated `X30 × X30 × M_{n-2}` overlap between (a) and
   (b), with ±1/2 symmetry factors where clusters are symmetric
   (a10–a21);
   (d) wide-angle soft corrections as SS-difference blocks multiplying
   `X30 × M` lines: `(SFF(..)+SFF(..)−SFF(..)−SFF(..))*E30FF(..)*M`
   (a22, a27).
   RV: minus the integrated counterparts (`J21`) of the R-term's
   antennae × `M_n`, plus `X30 × M^{1-loop}`, plus the
   `(X31 + X30·J21)` closures.
5. **Completeness check before generating**: every mode of the
   process's check program (its `stitle` list is the limit checklist)
   must be covered by at least one line; `makeRRcheck`'s
   `autoRRX40/M0/SS.map` split shows your classification back to you.

## Structural patterns per contribution type

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
"5||6 collinear", "6 soft", "triple collinear 567"). Locate the `XX` lines
whose antenna covers that limit (antenna arguments contain the unresolved
parton(s)) and check: correct antenna type for the parton species
(quark/gluon, FF/IF/FI/II), correct mapped arguments in the reduced ME and
JET function, correct relative sign of the iterated `X30×X30` lines, no
missing limit entirely. Compare against the analogous term in a
neighbouring process (`maple/process/` is full of worked examples).

## Next step

After editing the `.map`, use the **autogen-subtraction** skill to
generate the Fortran and hook it into the spike-test build, then
**run-spike-test** to validate. Loop back here if a limit still fails.
