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
make skipdepend=true -j
```

(`make mode=skipdepend -j` is equivalent — `driver/makefile` maps `mode=skipdepend` to
`skipdepend=true` internally.)

This builds the `NNLOJET` target without regenerating dependency files, compiling only the
objects whose source actually changed (per existing `.o` timestamps/deps).

## When NOT to use this

Do a normal `make -j` (no `skipdepend`) if you:
- added or removed a source file
- changed which modules/includes a file depends on (new `use`/`include` statements)
- are unsure why a build is behaving unexpectedly (stale deps can mask real breakage)

If a fast rebuild produces confusing link errors or missing symbols, fall back to a full
`make -j` (or `make clean && make -j`) before debugging further.
