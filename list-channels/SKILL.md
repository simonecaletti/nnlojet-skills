---
name: list-channels
description: >
  List the partonic channels of an NNLOJET process — the numeric channel
  IDs used in runcard CHANNELS blocks, their partonic content, matrix
  elements, and perturbative layer. Use for questions like "which channels
  does ZJ have", "what is channel 5 of epem", "which channel numbers are
  the qg-initiated NLO ones", "what goes in the CHANNELS block".
  Read-only skill — modify nothing.
---

# Listing NNLOJET partonic channels

## Sources, in order of preference

1. **`driver/process/<PROC>/selectchannel<PROC>.f`** — the ground truth.
   The `getchannel<PROC>` subroutine maps each channel ID to its strings:

   ```fortran
   chan_str = 'ep em to db g d     '
   me_str   = 'B1g0Z(3,4,5,1,2)    '
   sub_str  = 'epemB2g0ZepemSNLO(...)'   ! empty for LO/V-type
   ```

   Layer headers (`c -- real`, `c -- virtual`, `c -- double real`, ...)
   partition the IDs into LO / R / V / RR / RV / VV.
2. **`driver/process/<PROC>/array<PROC>.dat`** — the full generated
   channel table (header line documents the record layout: iproc,
   initial state, final state, ME, colour/symmetry factors, PDG codes).
3. **`gen-channels`** for filtered queries:

   ```bash
   cd python/gen_channels     # MUST run from here (cwd-relative config)
   ./gen_channels.py <PROC> -rc -IP1 u -IP2 g        # runcard ID list
   ./gen_channels.py <PROC> -t                        # org-mode table
   ```

   Selectors: `-IP1`/`-IP2` (initial partons, multiple values = OR),
   `-ME` (matrix element), rejectors prefixed `r` (`-rNF 2`). Python 3.
   See `python/gen_channels/EXAMPLES.md`.

## Runcard CHANNELS block

Wildcards (no IDs needed): `ALL`, `LO`, `NLO`, `NNLO`, `V`, `R`, `VV`,
`RV`, `RR`. Numeric IDs from the sources above for specific partonic
channels; optional `region = a|b|all` for the RR phase-space split.

## Notes

- Channel IDs are process-specific and NOT stable across processes;
  never reuse numbers between processes.
- Driver-level process name required (aliases resolved in
  `Process.f90::init_proc` — see list-processes skill).
- The channel decomposition mirrors the `driver/maple/<PROC>.map` channel
  lists ([ME, particle content, colour factor] per layer) — for the
  generating definition, see that file (add-process-to-driver skill).
