---
name: run-layer-check
description: >
  Run the NNLOJET layer check (test/layer_check) for a process and report
  pass/fail. Use whenever the user asks to run a layer check, verify that
  integrated subtraction terms are correctly linked in the virtual layers,
  or check the consistency of selectchannel/qcdnorm/bridge files against the
  maple subtraction sources. Symbolic Python+FORM tool — no compilation, no
  NNLOJET binary, no runcard involved.
---

# Running the NNLOJET layer check

The layer check verifies symbolically that every unintegrated antenna in
the real layers (SNLO, S) has its correctly-integrated counterpart in the
virtual layers (TNLO, T, U) with the right crossing, symmetry, and colour
factors — including the finite flavour-changing dipoles that pole checks
cannot see. In scheme language it proves, per process, the arrow
structure of Fig. 3 of arXiv:1301.4693: each dσ^S block reappearing in
dσ^T/dσ^U as `J21`/`J22` integrated dipoles (see
antennae-naming-convention). Complementary to spike tests
(run-spike-test), which validate the unintegrated terms numerically, and
to pole checks (run-pole-check), which validate the integrated terms
numerically and pointwise.

## Requirements

- Python >= 3.11, with `pyyaml`, `tabulate`, `argcomplete`
  (`pip install -r requirements.txt`; argcomplete is imported
  unconditionally, so it is required).
- The FORM binary `form` on PATH. Maple is NOT needed.
- cwd MUST be `test/layer_check` — all paths in the tool are relative.

## Step 1 — check the databases

`--proc_name` takes the **driver-level** process name (as in
`driver/process/<P>/selectchannel<P>.f`). If the maple-level and
driver-level folder names differ, the process MUST have an entry in
`databases/proc_folder_db.yml` before running:

```yaml
Z:   {maple: DY, driver: Z}
ZJ:  {maple: Z,  driver: ZJ}
GJ:  {maple: G,  driver: G}
```

Absent an entry, both default to the process name. Verify the two folders
exist (`maple/process/<maple>/`, `driver/process/<driver>/`) and add the
entry if missing. `databases/proc_list.yml` only drives `--all`; add the
process there too if it should be part of the full sweep.

## Step 2 — run

```bash
cd test/layer_check
./layer_check.py --proc_name <P> --skip        # NLO + NNLO
./layer_check.py --proc_name <P> --skip --nlo_only
```

`--skip` is mandatory when run non-interactively — without it the script
blocks on a y/n prompt. Runtime is seconds to ~a minute per process.
Other flags (use only if the user asks): `--channel "u,ub;d,db"`
(restrict luminosities; partial sums need not vanish unless the channel
set is closed), `--debug` (channel tables), `--chan_num` (tag factors
with channel numbers), `--form` (reuse existing `src/<P>_layers.inc`),
`--all` (full sweep over proc_list.yml).

## Step 3 — report

The check is automatic: the FORM sums `sumNLO` / `sumNNLO` must be
literally `0`.

- Output "Layer check passed at NLO/NNLO" → **report to the user that the
  test passed** (per order).
- Otherwise → **report that it did not pass** (and at which order). The
  script prints the non-cancelling residue; do not attempt fixes on your
  own — leave interpretation to the user unless asked. (The residue
  symbols — `calX30`, `J21*`, ... — decode via the
  antennae-naming-convention skill.)
- `RuntimeError ... Check layer_check.log` → FORM itself failed; report
  and point at `layer_check.log`.
- `ValueError: "Channels and QCD factors must be in equal number"` →
  `selectchannel<P>.f` and `qcdnorm<P>.f` are out of sync — this is a
  real inconsistency in the process files, not a tool problem; report it
  as such.

## layer_checked.md — hands off

A run without `--channel` rewrites the tracked status table
`layer_checked.md` in place. This file is updated **manually by the
developers only**: after every run, restore it with

```bash
hg revert layer_checked.md
```

and never edit it. (`src/` output is hg-ignored; `include/*.inc` are
tracked — do not modify them.)
