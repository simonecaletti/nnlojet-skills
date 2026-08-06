---
name: register-subtraction
description: >
  Generate the Fortran for an NNLOJET subtraction term from its maple .map
  file (makefortRR / makefortRV / makeformVV) and register it so the spike
  test compiles and links it. Use after writing or editing a .map file
  (write-subtraction skill), or when the user asks to regenerate auto*.f
  files or hook a subtraction term into a spike test. Scope is
  spike-test-only: full-process/driver registration is a separate concern.
---

# Generating and registering a subtraction term (spike-test scope)

Chain position: write-subtraction → **register-subtraction** →
run-spike-test.

## Scope

IN: `maple/` generators → `src/process/<DIR>/auto*.f` → `NNLOJET.mk` +
test makefile + check program.
OUT (do NOT touch): `driver/maple/` (`makeproc`, `updateproc.sh`) and
`driver/process/<PROC>/sig*.f` — those serve the full matrix elements and
full-process registration, handled by a separate skill.

## Step 1 — consistency check (fast, do this first)

```bash
cd maple                              # MUST run from maple/: relative paths
maple makeRRcheck -Diprocess=<N>      # RR; makeRVcheck for RV
```

`<N>` comes from `maple/iprocess.map` (e.g. epemZH2bb = 320, DY = 5).
Catches `aXX` numbering gaps ("ERROR: Gap in line numberings") and
regenerates the structure split `autoRRX40.map` / `autoRRM0.map` /
`autoRRSS.map` in the process dir.

## Step 2 — generate the Fortran

```bash
cd maple
maple makefortRR -Diprocess=<N>       # all *S / *SNLO files of the process
maple makefortRV -Diprocess=<N>       # all *T / *TNLO files
```

- Output: `src/process/<DIR>/auto<name>.f` (dated header
  `generated using makefortRR on ...`). The generator loops over ALL
  terms of the process, not just the edited one — and the current
  generator emits different boilerplate from what is checked in
  (`double precision` instead of `dimension` for the partons/facnorm
  arrays; `set_flav_perm`/`unset_flav_perm` moved inside the
  `if(ipass.eq.1)` block), so EVERY auto*.f of the process shows as
  modified even when only one changed semantically (all ten epem files
  rewrote for a one-file edit). **Back up `src/process/<DIR>/` before
  generating; afterwards restore the files you did not intend to touch**,
  so the working tree shows only your real change.
- Verify the edit changed only what you intended: diff the PHYSICS lines
  against the backup — extract the lines matching
  `set_map|wt\(|getqcdnorm|set_flav_perm|\* FF:` from old and new, strip
  leading whitespace, and diff those. That isolates the semantic content
  from the boilerplate drift above.
- Watch stdout for `ERROR: <file>.map > term #N > invalid ME argument` or
  `left-over list` — the primary symptom of a malformed .map; go back to
  write-subtraction.
- **VV is two-stage and needs FORM**: `maple makeformVV -Diprocess=<N>`
  writes `auto*.frm` into the process dir plus a shell script
  `makeformVV-update`; running that script executes `form auto*.frm`,
  moves the resulting `auto*.f` to `src/process/<DIR>/`, and `hg add`s it.
- In the generated Fortran: `aN` → `wt(N)`, antenna → `FullA30FF(...)`,
  `JETnm` → the process jet-function call (e.g. `ecuts_epem_vh`), mappings
  → `call set_map(...)`, reduced-ME args → `j1..jN`. Useful for reading,
  never for editing.

**No maple (or no FORM) available**: report it to the user and stop the
generation step — NEVER hand-edit `auto*.f` to mimic what the generator
would do. Only already-checked-in `auto*.f` can be spike-tested then.

## Step 3 — register (only needed for NEW terms)

For a regenerated existing term, nothing to register — go to
run-spike-test. For a new `auto<name>.f`, three touchpoints:

1. `NNLOJET.mk` (repo root): append the filename to the process's
   `SUBTRACTION_<PROC>` variable (e.g. `SUBTRACTION_VH += auto....f` for
   epemZH). VPATH is flat — the filename must be globally unique.
2. `test/process/<PROC>/makefile`: append to the subtraction block,
   `LIBFILES += auto<name>.f`.
3. The check program (`check<N>to<M>[loop].f`, or the factored-out
   `test*<PROC>.f` / `tests*<PROC>.f` files): add a `case(<itype>)` pair —
   full ME in `test(itype)`, subtraction in `tests(itype)` — with matching
   momentum arguments. Pick the next free channel number.

## Traps

- `driver/process/<PROC>/` may contain stale COPIES of `auto*.f`; only
  `src/process/<DIR>/` is on the build VPATH. Ignore/flag the copies,
  never edit them.
- The generators hardcode `read ./notation.map` and
  `fortdir := ../src/...` — running from anywhere but `maple/` fails.
- Mercurial repo: `hg add src/process/<DIR>/auto<name>.f` for new files
  (makeformVV-update already does it for VV). Do not commit unless asked.

## Next step

Use the **run-spike-test** skill. If a limit fails there, loop back to
write-subtraction → here → run-spike-test.
