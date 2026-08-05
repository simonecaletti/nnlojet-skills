---
name: fast-rebuild
description: Use when rebuilding NNLOJET after editing existing Fortran/C/C++ source files (not after adding/removing files or changing includes) and a full dependency-checked build would be too slow. Skips makedepend.sh dependency regeneration for a much faster incremental compile.
---

# Fast rebuild (skip dependency regeneration)

NNLOJET's `driver/makefile` regenerates per-file `.d` dependency files via `makedepend.sh` on
every build (the `depend` target, included via `-include $(DEPS)`). This scan is slow across the
whole tree and is unnecessary when you've only edited the *bodies* of existing source files
without changing `include`/`use`/module dependencies.

## Usage

```
cd driver
make skipdepend=true
```

No `-j` — this build must run serially (see below). (`make mode=skipdepend` is equivalent —
`driver/makefile` maps `mode=skipdepend` to `skipdepend=true` internally.)

This builds the `NNLOJET` target without regenerating dependency files. Note: with
`skipdepend=true` the `DEPS` variable is left **empty** (`driver/makefile:602-608` — the `else`
branch that populates it is skipped), so `-include $(DEPS)` at `driver/makefile:740` includes
nothing and the existing `.d` files in `driver/deps/` are never read. Object selection is by
`.o` timestamp ONLY; recorded dependencies play no part.

## Why not `-j`

Empty `DEPS` means there are NO inter-object ordering constraints at all, so a parallel make
compiles modules concurrently with their consumers. Recognisable symptom: you add a symbol to a
module (say a new function in `EvalFuncs.f90`) and bind it in `Observables.f90`, and the
parallel build dies with

```
Error: Symbol '<new_symbol>' at (1) has no IMPLICIT type
```

because the consumer compiled against the stale `.mod`. Serial `make skipdepend=true` avoids
the race.

## When NOT to use this

Do a normal `make -j` (no `skipdepend`) if you:
- added or removed a source file
- changed which modules/includes a file depends on (new `use`/`include` statements)
- added or renamed a public module entity (function, variable, type) consumed by another file
- are unsure why a build is behaving unexpectedly (stale deps can mask real breakage)

If a fast rebuild produces confusing link errors or missing symbols, fall back to a full
`make -j` (or `make clean && make -j`) before debugging further.
