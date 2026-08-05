---
name: list-observables
description: >
  List which observables are available for a given NNLOJET process, or
  find which processes provide a given observable. Use for questions like
  "what observables does ZJ have", "can I cut on ptj3 in DIS", "is there
  an mll observable for W processes". Read-only skill — modify nothing.
---

# Listing NNLOJET observables

One registry serves everything: an observable name is simultaneously
usable in runcard HISTOGRAMS, SELECTORS, HISTOGRAM_SELECTORS, and SCALES
(mur/muf), case-insensitively.

## Preferred: the binary

```bash
cd driver && ./NNLOJET -listobs <PROCESS>
```

Build-accurate (calls init_proc + init_obs, then prints every bound
observable with its description). Use the USER-FACING process name (the
one a runcard would use — aliases are resolved by init_proc). Note it
initialises the jet algorithm as `jade` for EPEM* processes, `antikt`
otherwise.

## Without a binary: read the source

`driver/core/Observables.f90::init_obs()`:

1. Find the `select case (name_proc)` branch containing the process
   (main select starts ~line 152; names are UPPERCASE; one branch covers
   many related processes, e.g. `case ("Z", "ZJ", "ZJ_H", ...)`).
2. Add the unconditional "default observables" block (~line 3054):
   njets, npartons, jet kinematics (ptj1..., yj1..., energy-ordered
   `*_energy` variants), etc. — available to every process.
3. Process-family extras appended after the main select: identified
   hadrons, P2B, epemZH block (~line 3501), forced DIS set (~line 3544:
   q2, x, y, W2, ... always active for DIS).
4. Skip entries with `digest=.false.` — internal, hidden from -listobs.

Each `bind_obs("name", "description", ...)` line = one observable;
`min_njets=`/`min_npar=` show its validity requirements (e.g. ptj2
needs ≥ 2 jets).

## Reverse lookup (which processes have observable X)

grep `bind_obs("X"` in `driver/core/Observables.f90` and report the
enclosing `case` branches, plus whether it sits in the default block
(= all processes).

## If the observable does not exist

Say so and point to the add-observable skill — never invent names: an
unknown observable in a runcard is a hard stop at parse time.
