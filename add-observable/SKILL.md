---
name: add-observable
description: >
  Add a new observable to NNLOJET so it can be used in runcard HISTOGRAMS,
  SELECTORS, HISTOGRAM_SELECTORS, and SCALES blocks. Use whenever the user asks
  to add, register, define, or expose an observable, a new histogram variable,
  a new selector/cut variable, or a new dynamic-scale variable — adding a
  selector IS adding an observable in NNLOJET; there is no separate registry.
---

# Adding an observable to NNLOJET

One `bind_obs` call makes a name simultaneously available to histograms,
selectors, histogram selectors, dynamic scales (`mur`/`muf`), `mu0=` and
`fac=`. Lookup is case-insensitive via `getIdFromName_obs`.

## Files you touch

1. `driver/core/EvalFuncs.f90` (module `EvalFuncs_mod`) — the eval function.
   For DIS-specific quantities use `driver/core/EvalDIS.f90` instead (see DIS
   note below).
2. `driver/core/Observables.f90` — the `bind_obs` registration inside
   `init_obs()`.

epemZH processes register in the extra `case ("EPEMZH2BB", ...)` block after
the main select, ~line 3501 of `Observables.f90`.

## Step 1 — write the eval function

Two common interfaces (defined in `driver/core/EvalObs.f90`):

```fortran
! simple: eval=
real(kind=wp) function eval_pt_V(npar)
  integer, intent(in) :: npar
  real(kind=wp), dimension(4) :: pV
  pV(1:4) = kin(npar)%p(1:4, npar - 1) + kin(npar)%p(1:4, npar)
  eval_pt_V = v1_pt(pV)
end function eval_pt_V

! parametrised: eval_mem= with imem (e.g. jet index)
real(kind=wp) function eval_pt_jet(npar, ijet)
  integer, intent(in) :: npar, ijet
  eval_pt_jet = v1_pt(kin(npar)%pjets(:, ijet))
end function eval_pt_jet
```

Kinematics come from `kin(npar)` (`KinData_mod`); helpers like `v1_pt`,
`v1_y`, `v2_delta_R` are in `ObsHelper_mod`. A *negative* `imem` is just a
plain parameter passed as `abs(imem)`; a *positive* `imem` means the function
memoizes internally and must reset its cache when called with `imem=0` — use
negative unless you implement caching.

## Step 2 — register in init_obs()

Find the correct `select case (name_proc)` branch in `Observables.f90`
(main select starts ~line 152; `name_proc` is uppercase). Process-agnostic
observables go in the unconditional "default observables" block (~line 3054).

```fortran
call bind_obs("mll", "invariant mass of Z", &
              EvalObs_t(eval=eval_minv_V), ifac=1)
call bind_obs("ptlm", "transverse momentum of l-", &
              EvalObs_t(eval_mem=eval_pt_lep, imem=-1), ifac=2)
call bind_obs("abs_yj1", "|yj1|", &
              EvalObs_t(eval_mem=eval_y_jet, imem=-1, decorator=abs_deco, min_njets=1))
```

`EvalObs_t` fields you may set:
- `eval=` / `eval_mem=`+`imem=` — the function pointer (exactly one).
- `min_njets=`, `min_npar=`, `min_nphotons=`, `min_flav_njets=` — validity
  guards checked *before* evaluation; e.g. `min_njets=2` for `ptj2`. Without
  the guard the eval can read garbage momenta on low-multiplicity events.
- `decorator=` — `abs_deco`, `sqrt_deco`, `sqrt_abs_deco`, `sqr_deco`
  (from `ObsHelper_mod`), applied after eval.
- `is_integer=.true.` — for counters (njets-like); makes selector
  comparisons tolerance-based.
- `needs_flav=.true.` — auto-enables flavour tagging when the observable is
  requested.
- `digest=.false.` — hide from `-listobs` (internal observables only).

`ifac` (optional bind_obs arg) feeds the min-shat phase-space optimisation in
`Selectors.f90`: `ifac=1` for mass/energy-like lower cuts (`m > cut` ⇒
`shat > cut²`), `ifac=2` for pt-like (`shat > (2·cut)²`). **When unsure,
omit it** — a wrong value tightens `smin` incorrectly, which is a physics
bug, not a crash.

Name limit 32 chars, description 224. `bind_obs` does no duplicate checking
— grep the branch first to avoid registering a clashing name.

## Step 3 — rebuild and verify

```bash
cd driver && make skipdepend=true -j     # OK: no new use statements/files
./NNLOJET -listobs <PROCESS>             # new name must appear with its desc
```

`-listobs` is build-accurate (it calls init_proc + init_obs); no extra step
is needed for the observable to show up. If you created a *new* .f90 file
instead, add it to `MODULES` in NNLOJET.mk and run a full `make -j`.

Then smoke-test with a runcard using the observable in a HISTOGRAM line
(`nbins`, `min`, `max` are mandatory — binning is runcard-only, bind_obs sets
no defaults) and/or a SELECTORS line. Run in a scratch dir so `driver/`
doesn't collect stray `.dat`/grid files. If the observable changes physics
selection for existing runcards, run the relevant regression tests
(`test/regression_tests/tests.py -rc ...`).

## Special cases

- **Energy-sorted jet observables**: `kin(npar)%pjets` is pt-ordered by
  default. For energy ordering, do not re-sort pjets; use the helper
  `get_jet_index_by_energy(npar, Nth)` (EvalFuncs.f90 ~line 3687) inside the
  eval function, following `eval_ejN_energy` / `eval_yjN_energy`. Bind one
  observable per rank with `imem=-N` and `min_njets=N`, naming convention
  `<obs>jN_energy`:

  ```fortran
  call bind_obs("ej2_energy", "energy of sub-leading jet (ordered in energy)", &
                EvalObs_t(eval_mem=eval_ejN_energy, imem=-2, min_njets=2))
  ```

  A new energy-sorted quantity only needs one new eval function reusing the
  same index helper.
- **DIS**: most DIS observables are getters over module state precomputed in
  `DIS.f90` (`eval_Q2_dis = value_q2_dis`), not `kin(npar)` calculations. Add
  the `value_*_dis` computation in `DIS.f90`, the getter in `EvalDIS.f90`
  (`EvalDIS_mod` is use'd *inside* init_obs, not at module scope). If
  internal code needs the observable by id, mirror the force-activation
  pattern at Observables.f90 ~line 3544 (`getIdFromName_obs(...,
  no_error=.true.)` into an `id_*_dis` variable in `DIS_mod`).
- **Composite** observables: `eval_comp=` + `is_composite=.true.`, 4-arg
  interface — called without `isub` returns the number of sub-observables;
  with `isub` returns that value and may set `wgt`/`valid`. Template:
  `EvalFuncs.f90` ~line 542. Needs `COMPOSITE ... END_COMPOSITE` in the
  runcard.
- **Manual** observables (`is_manual=.true.`): no eval; value pushed via
  `setManual_obs` (see p2b_qt, ~line 3478).
- `cross` and `COMPOSITE` are histogram-level special cases in
  `Histograms.f90`, not observables — never a model to copy.
- `cumulant` is a per-histogram runcard option, unrelated to registration.
