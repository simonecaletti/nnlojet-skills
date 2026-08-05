---
name: list-test-coverage
description: >
  Report what validation exists for a given NNLOJET process: spike tests,
  layer-check status, regression-test runcards, and example/paper
  runcards. Use for questions like "is ZJ spike-tested", "does DIS pass
  the layer check", "are there regression tests for epem", "which paper
  runcards exist for H". Read-only skill — modify nothing.
---

# Test coverage of an NNLOJET process

Check the four validation layers; report per layer what exists and,
where recorded, its status. Names may differ across layers — resolve via
the list-processes skill first (e.g. maple epemZH → test epemZH2bb/2gg).

## 1. Spike tests (numerical, unintegrated subtractions vs ME)

`test/process/<PROC>/` — existence of the dir + which `check<N>to<M>.f`
programs it holds tells which contributions are testable (`loop` suffix
= RV; small N−M gap = R, large = RR). No dir → no spike tests
(write-spike-test skill can create them). To run: run-spike-test skill.

## 2. Layer check (symbolic, integrated-subtraction linking)

`test/layer_check/layer_checked.md` — the status board (per process:
pass at NLO / NNLO, ⏳ NLO-only, ❌ known failures). READ ONLY — this
file is manually maintained; never modify it. Whether a process is
runnable at all: `test/layer_check/databases/proc_list.yml` +
`proc_folder_db.yml`. To run: run-layer-check skill.

## 3. Regression tests (end-to-end integrated output)

`test/regression_tests/runcards/<PROC>/` — one subdir per process,
`.run` files inside; runcards containing an uncommented `point_check`
token also run in fast mode (`tests.py -pc`). Validation reference data:
`test/regression_tests/VALIDATION_OUTPUT/<PROC>/<runcard>.run.tar.gz`
(existence = validated reference exists; TEST_DETAILS.txt inside records
the hg revision it was generated at).

## 4. Example / paper runcards

`examples/<PROC>/<arxiv#>/` — published-paper benchmark runcards with
sibling `*result*` dirs (`data/*.dat`, `Histograms/*.pdf`) as reference
output; `examples/<PROC>/test/` — lighter test cards.

## Extras

- Per-process standalone Fortran checks: `test/process/<PROC>/` beyond
  spike programs (some dirs hold ME cross-checks, e.g. against MadGraph
  — see local READMEs).
- Citations for the process: `driver/bin/nnlojet-references.py <PROC>`.
- Core-infrastructure tests (not per-process): `test/core/`,
  `test/antenna/`, `test/flav/`, `test/hfrag/`, `test/layer_check/`,
  `test/pynnlojet/`.
