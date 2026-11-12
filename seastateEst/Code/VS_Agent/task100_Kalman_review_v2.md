# 3-DOF Kalman follow-up (v2)

## Direct answers

- **Keep mooring enabled in both runs for the comparison you described.** The current `main.fst` and `main_lin.fst` both have `CompMooring=3` and reference the same `MoorDyn.dat`. Mooring contributes the low-frequency horizontal restoring force and therefore affects the long-period surge mode. Disabling it is a useful separate ablation, not a fix for the current estimator; it changes the physical mode you want to estimate. Keep the mooring setup, initial state, and operating condition aligned.
- **Using OpenFAST velocity truth for the radiation convolution is acceptable for an offline validation stage.** Label it as an oracle input, then remove it in a later sensor-only stage.
- **The full-record result says the heave flatline is not fixed by the candidate Q/R values.** On all 12,001 samples (300 s), the candidate's heave estimate had standard deviation about `0.00105 m`, compared with `0.0309 m` for OpenFAST; its heave RMSE was `0.0405 m`. Changing the continuous force-bias intensity multiplier from `0` through `1e12` left estimate standard deviation in the narrow range `0.0010–0.0015 m`. That is evidence to stop tuning that scalar and check excitation/model observability.

## Important configuration clarification

The current files confirm the intended DFT reference: `HydroDyn.dat` has `ExctnMod=1` (DFT) and `RdtnMod=1` (convolution). `HydroDyn_lin.dat`, used by `main_lin.fst` to produce the `.lin` data, has `ExctnMod=2` and `RdtnMod=2`. This is consistent with the stated goal of validating against DFT output while extracting mechanical coefficients from a linearization, but it means the `.lin` operating point at 299 s came from a different excitation model/run than the DFT `main.out`. The first review's statement that those two settings were reversed was incorrect.

`main_lin.fst` linearizes at `t=299 s` with `CalcSteady=False`; that operating point is not the initial condition of the DFT record. Its active heave position is about `-0.4314 m`, while the DFT record starts at `-0.44 m`. More importantly, its heave acceleration at the linearization point is about `0.00245 m/s^2`, while the DFT run's heave acceleration has mean near zero. The v1 adjustment to subtract the linearization acceleration from the entire DFT record barely changed the full-run result, but the operating point should still be aligned before treating the model as validated. Prefer an equilibrium/quiet-water linearization or a point from a matching run, rather than a late-time point from the separate state-space-excitation run.

## Why Q/R tuning cannot yet solve it

The estimator uses an excitation state-space fit (`marin_semi.ssexctn`) while the reference simulation computes excitation by DFT. Since the user has independently found that their excitation forces differ significantly, the model's heave forcing is not the forcing that generated the measured heave acceleration. A large acceleration innovation may then be explained by excitation/wave states or force-bias states instead of a coherent heave displacement trajectory. Q/R can trade off those explanations, but cannot make the excitation model correct.

Also, surge/heave positions are not directly measured. A low-frequency acceleration measurement does not uniquely determine displacement without sufficiently accurate dynamics and forcing. The long-period surge response is particularly sensitive to mooring restoring and force bias; keep the MoorDyn stiffness consistent.

## Suggested validation sequence

1. **Freeze the physical case:** same MoorDyn file/flags, DOFs, hydrostatic data, initial state, and wave realization. Keep DFT as the reference if that is the target.
2. **Use a matching operating point:** linearize a quiet-water equilibrium or a matching DFT run. Avoid applying the `t=299 s` state and acceleration from the separate state-space-excitation run to the full DFT history.
3. **Validate the mechanical model before the estimator:** use OpenFAST truth only as diagnostic data to check whether the extracted `M`, `K`, `C`, and radiation-memory force predict heave/surge acceleration for known states and known excitation. Do not tune a filter to compensate for a failed forward-model check.
4. **Isolate heave:** use a heave-only case with aligned operating point and known excitation force. First compare predicted and OpenFAST heave acceleration; then add the heave accelerometer measurement to the KF. This distinguishes bad `Yx`/sign/state ordering from excitation mismatch and weak observability.
5. **Replace the excitation model deliberately:** for the oracle validation, feed the reference DFT excitation force to the structural estimator as a known input and confirm that the heave trajectory is recovered. Then replace this oracle with the excitation estimate intended for the proof of concept and quantify the loss. Do not call the oracle-force result the sensor-only estimator.
6. **Tune after these checks:** set accelerometer/inclinometer `R` from the chosen idealized measurement-noise assumption, scale `P0` by state units, and tune only the remaining model-discrepancy covariance using full-record innovations and trajectory errors.

## Reproducible diagnostic

`_temp/diagnose_full_record_v1.py` compares its vectorized radiation convolution against `FRad.Conv_Rdtn` at several points and runs the candidate through the full 300 s record. The plot is `_temp/surge_heave_full_record_v1.png`. This full-record test supersedes the earlier 30 s-only score: the first 30 s looked more promising, but the complete record does not show usable heave tracking.
