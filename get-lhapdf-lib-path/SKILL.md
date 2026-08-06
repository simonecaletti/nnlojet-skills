---
name: get-lhapdf-lib-path
description: >
  Fix the NNLOJET startup failure "error while loading shared libraries:
  libLHAPDF.so: cannot open shared object file". Use whenever ./NNLOJET
  (or any tool linking LHAPDF) dies with a loader error before printing
  any output, or before running the binary in a fresh shell.
---

# LHAPDF library path

```bash
export LD_LIBRARY_PATH="$(lhapdf-config --libdir):$LD_LIBRARY_PATH"
```

Needed once per shell before any `./NNLOJET` invocation. The symptom is
a loader error *before any NNLOJET output appears*:

```
./NNLOJET: error while loading shared libraries: libLHAPDF.so:
cannot open shared object file: No such file or directory
```

This is an environment problem, not a build failure — the binary is
fine. Confirm with `ldd driver/NNLOJET | grep -i lhapdf` (should resolve
after the export).

No general wrapper handles this: the only in-repo occurrence is
`driver/singlerun.sh`, which sets `LD_LIBRARY_PATH` from its 3rd
argument for cluster jobs — direct `./NNLOJET` runs need the export
above.
