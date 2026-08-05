---
name: list-processes
description: >
  Orient inside the NNLOJET repository at the process level: list the
  driver-level and maple-level process names, resolve the correspondence
  between them, explain naming suffixes, and locate the directories of a
  process across all layers. Use for questions like "which processes
  exist", "what is the maple folder for ZJ", "what does the _h suffix
  mean", "where is the code for DIS". Read-only skill — modify nothing.
---

# NNLOJET process taxonomy and correspondences

## The four layers of a process

| layer | location | content |
|---|---|---|
| maple sources | `maple/process/<DIR>/` | subtraction terms (.map) |
| Fortran MEs + subtractions | `src/process/<DIR>/` | amplitudes, auto*.f |
| driver glue | `driver/process/<PROC>/` | sig*.f, selectchannel, qcdnorm |
| spike tests | `test/process/<PROC>/` | check*.f |

Directory names are NOT the same across layers — resolve correspondences
before answering.

## Driver-level names

- `ls driver/process/` — the compiled process dirs.
- User-facing runcard names are ALIASED in
  `driver/core/Process.f90::init_proc` (e.g. `epL* → DIS*`,
  `epN* → DISWm*`, `epNb* → DISWp*`) — a runcard name may not match any
  directory; grep init_proc for the mapping.
- `./NNLOJET -listprocs` prints a HARD-CODED, hand-maintained list — do
  not trust it as build introspection; prefer the directories +
  `Process.f90`.

## Maple-level names and the correspondence

- `ls maple/process/` — maple dirs.
- **`maple/iprocess.map` is authoritative**: iprocess number → directory,
  jet function (`MYJET`), subprocess index ranges, extra modules.
- `driver/maple/makeproc` has its OWN numbering for driver generation
  (header comment block, e.g. `112: epemZH2bb`).
- maple↔driver folder correspondence for mismatching names:
  `test/layer_check/databases/proc_folder_db.yml` (e.g. `Z: {maple: DY,
  driver: Z}`, `ZJ: {maple: Z, driver: ZJ}`, `GJ: {maple: G, driver: G}`).
- Not 1:1: one maple dir can serve several driver/test processes — maple
  `epemZH` (iprocess 320/330/340) → `epemZH2bb` / `epemZH2gg` /
  `epemZH2gaga`; one src dir can serve several drivers similarly.

## Suffix conventions (consistent across the whole tree)

- `_h` = hadronic, i.e. process WITH FRAGMENTATION (identified hadrons);
- `_EW` = electroweak corrections;
- `J`, `JJ`, `JJJ`, `JJJJ` = additional jet multiplicity;
- `unsym` = unsymmetrised variant;
- `fc` = flavour-changing (e.g. `GGJfc`, `2jetfc`);
- `prod` = production-only map (decays applied on top via `applydecay`);
- decay-specified variants: `H2/H3/H4` (2/3/4-particle Higgs decay),
  `Hto2p...`, `epemZH2bb` etc.;
- `_frag` = photon/hadron fragmentation variants;
- `_OL` (ME names) = OpenLoops amplitudes.

## Related lookups

- Citations for a process: `driver/bin/nnlojet-references.py <PROC>`.
- Which observables a process has: list-observables skill.
- Channel content: list-channels skill. Test coverage:
  list-test-coverage skill.
