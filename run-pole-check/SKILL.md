---
name: run-pole-check
description: >
  Run and interpret the NNLOJET numerical pole checks: verify that the
  explicit infrared epsilon-poles of the virtual contributions (V, RV, VV)
  are cancelled pointwise by their subtraction terms (TNLO, T, U). Use when
  the user asks to check pole cancellation, validate a virtual or two-loop
  contribution, run pole_check, or verify the U/VV layer — which has no
  unresolved radiation and therefore CANNOT be spike-tested; this is its
  only numerical validation.
---

# Running pole checks

The virtual layers carry explicit 1/ε poles; the subtraction adds them
back as integrated dipoles (`J2^(1)` for V/RV, `J2^(1)`⊗`J2^(1)` and
`J2^(2)` for VV — the T,a and U blocks of Fig. 3, arXiv:1301.4693; see
antennae-naming-convention). Pole cancellation is POINTWISE: at any
single phase-space point, matrix and subtraction must agree at every
pole order. No infrared limit is involved — this is the key contrast
with spike tests.

Validation triangle: run-spike-test = unintegrated terms, numerically,
in IR limits; run-layer-check = integrated bookkeeping, symbolically;
**run-pole-check = integrated terms, numerically, pointwise**. Their
blind spots differ: the layer check sees the finite flavour-changing
dipoles that have no pole; the pole check sees numerical/crossing
errors that symbolic cancellation cannot.

## Route 1 — runcard `pole_check` (full-channel sweep, primary)

1. Template: a `*_point_chk.run` runcard for the process
   (`test/regression_tests/runcards/<PROC>/`), which already has the
   `point_check` flow enabled — the pole check runs inside it
   (`sig.f:202-207`: `polecheck` is called when `pole_check_run` and
   the virtual weights are nonzero).
2. In the RUN block set:

   ```
   pole_check = .true.
   ```

   (boolean; parsed at `IOHelper.f90:1842` → `pole_check_run`).
3. `CHANNELS`: select the virtual layers — `V`, `RV`, or `VV` (or
   numeric IDs). Non-virtual channels are skipped by the check.
4. Run in a scratch directory (`driver/` litters otherwise); if the
   binary dies on a `libLHAPDF.so` loader error, see get-lhapdf-lib-path.
   `OMP_NUM_THREADS=1` for readable ordered output.

## Output and pass criterion (`sig.f` polecheck, ~:300-360)

Per channel, each singled out in turn:

- 1-loop layers (V, RV): `Pole check ieorder: -2` then `-1`, each with
  `matrix:`, `subtraction:`, `ratio:`.
- 2-loop layer (VV): `Pole check ipole: -4` … `-1` (for ipole ≥ −2 the
  one-loop insertions are set to the same ε order; `common /order2l/`).

**PASS** = `ratio: 1.000000...` to numerical precision at EVERY printed
order and channel — pointwise, no trend to inspect. `no pole` (both
zero) is fine at orders where the channel has none. **FAIL** signatures:
matrix printed with no subtraction (or vice versa) = missing
counterpart at that order; ratio a constant ≠ 1 = wrong factor or
crossing in the integrated-dipole assembly — run-layer-check names the
offending `J21`/`J22` symbolically. Precision-level wobble (1e-10) is
numerics, not physics.

## Route 2 — test programs (single-term isolation)

The two-argument spike programs take the ε order directly:

```bash
./check4to3loop IORDER CHANNEL     # IORDER: 0 = finite, -1, -2 = poles
```

(the 5 two-arg programs of the run-spike-test census, e.g.
`epem/check4to3loop.f`, which checks that ME and subtraction agree for
"pole and finite term" also inside IR limits; plus dedicated programs
like `GGJfc/check_poles_V.f`, `check_poles_T_a.f`). Use this route to
isolate ONE term/channel after Route 1 flags it; build/argument/iplot
mechanics as in run-spike-test.

## On failure

Report channel + layer + ε order + ratio. The fix target is usually the
integrated side: the `J21`/`J22` lines of the `*T`/`*TNLO`/`*U` maple
files or their crossing factors (write-subtraction → autogen-subtraction
→ rerun), cross-checked symbolically with run-layer-check. Repo is
Mercurial; revert scratch edits, commit nothing unless asked.
