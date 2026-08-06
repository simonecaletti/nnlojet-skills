---
name: clone-process
description: >
  Derive a new NNLOJET process from an existing one instead of authoring
  from scratch. Use when the user asks to clone, copy, duplicate, or derive
  a process, set up a process "as a copy of" another, create a sandbox or
  variant of an existing process, or add a crossing of an existing process.
  Default deliverable is the copied .map + registration; generation
  (makeproc) and building happen only if the user asks.
---

# Cloning a process

Scope: PROCESS-level cloning — folders, maps, registration numbers.
Deriving a single subtraction TERM by crossing another process's
`.map` (un-crossing legs, antenna renames, Full→split) is a different
job: write-subtraction, section "Deriving a term by crossing an
existing one".

Prior art in-repo: `maple/process/GG/shell/GGfromZ.sh` (generates GG
subtraction terms from Z's), and the `DYtest` clone of `Z` in this tree
(makeproc iprocess 127, iprocess.map 50) — use DYtest as the worked
example. Colour factors are unchanged under crossing; only ME suffixes
and lepton symbols differ — the crossing suffixes are one-line wrappers
around a single shared implementation (see the me-naming-convention
skill and add-process-to-driver).

## 1. Resolve the source process FIRST

Names and folders differ (no process is named `DY` — that is `Z`'s
maple/src folder). Look the source up in
`../add-process-to-driver/databases/proc_registry.yml`: both iprocess
numbers and the maple / driver / src folders.

## 2. Decide the sharing boundary — ASK the user

Three independent axes:

- matrix elements (`src/process/<DIR>`)
- subtraction terms (`maple/process/<DIR>`)
- generated driver files (`driver/process/<PROC>`)

Sharing MEs while copying subtraction terms is a normal choice for a
sandbox (DYtest does exactly this: own maple dir, shares `DY` MEs and
`FLAVlistZ.map`). Do not assume — ask which axes the clone owns.

## 3. Copy

- `driver/maple/<SRC>.map` → `driver/maple/<NEW>.map`. The filename base
  MUST equal the `myname` you register (makeproc does
  `read cat(myname,`.map`)`, `makeproc:1112`).
- If subtraction terms are copied:
  `maple/process/<SRCDIR>/` → `maple/process/<NEW>/`.

## 4. Register in driver/maple/makeproc

Header comment line (`# iprocess = N:  "<NEW>":`) + `elif` dispatch
block copied from the source branch with `myname`, `mydir`, `mymapdir`
repointed; keep `MYJET` and `mynameflag` from the source. Pick an unused
N — and grep BOTH numbering tables for it first: makeproc numbers and
`maple/iprocess.map` numbers are independent, and a collision is silent
(Z_EW's stale `imapprocess:=50` vs DYtest's iprocess.map entry 50 is the
in-tree cautionary example; see add-process-to-driver).

## 5. iprocess.map entry (only if the clone owns a subtraction dir)

Add to `maple/iprocess.map`: a new number, the header comment line, and
a `getiprocess()` branch copying the source's channel ranges with
`directory` repointed to the new folder. Keep `iprocessname` pointing at
the SOURCE's `FLAVlist<name>.map` when matrix elements are shared (no
new FLAVlist needed — DYtest keeps `iprocessname:=Z`). Set the makeproc
branch's `imapprocess` to this new number.

## 6. Verify without generating (cheap, always do)

- iprocess.map branch resolves and points at the right directory:

  ```bash
  cd maple
  printf 'read `./iprocess.map`:\nprint(getiprocess(<N>,1,RR)):\n' > /tmp/chk.mpl
  maple -q /tmp/chk.mpl
  ```

- makeproc edit parses, without side effects: copy makeproc, make the
  loop range empty (`for iprocess from N+1 to N do`), run
  `maple <copy>` from `driver/maple/`; exit 0 with only pre-existing
  warnings means the edit parses. Delete the temp copy afterwards.

## 7. Hand off

Generation and building only if the user asks:
**add-process-to-driver** (set the loop bound, run `maple makeproc`) →
**link-to-driver** (auto*.f stubs + build/core registration). If the
verification ever invokes `./NNLOJET` and it fails with a `libLHAPDF.so`
loader error, see get-lhapdf-lib-path.

Mercurial: `hg add` every new file (`driver/maple/<NEW>.map`,
`maple/process/<NEW>/`); do not commit unless asked.
