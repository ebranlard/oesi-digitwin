# 3-DOF Kalman filter review (v1)

## Bottom line

This is not only a Q/R tuning problem. The original filter mixes absolute OpenFAST output coordinates with perturbation states from the `.lin` model, assigns an acceleration-style variance to a pitch-position measurement, and uses the OpenFAST velocity truth to calculate the radiation-memory input. There is also a model/run configuration mismatch to resolve before the covariances can be called final.

The versioned candidate script applies the first two corrections, uses an exact discrete process covariance, and includes a provisional covariance candidate. It does not make the translational position estimates reliable enough to call the 3-DOF estimator complete.

## Findings

1. **Operating-point offset is a fundamental state error.** `main_lin.1.lin` contains nonzero operating-point values (`x`): active positions are approximately `[0.0004994, -0.4314, 9.284e-8]` m/rad, while the initial `main.out` positions are `[0, -0.44, 0]`. The linearized matrices describe deviations around that operating point. The original script sets the heave state to the absolute `-0.44` m and compares absolute pitch to a perturbation-state output. With a heave restoring coefficient near `11.4 s^-2`, that absolute initialization introduces a several-m/s² false restoring acceleration. The candidate subtracts the `.lin` operating point from position, velocity, and acceleration channels; it adds the position offset back for plotting. Note the `.lin` was requested at 299 s while the estimator starts its output record at 0 s.

2. **Pitch `R` is in the wrong scale in the original.** Pitch is measured as position in radians, but its original variance is `1e-3 rad^2` (`sigma = 0.0316 rad`, or 1.81 degrees). The simulated pitch signal standard deviation is only `0.00091 rad` (0.052 degrees), so this setting effectively discards the pitch measurement. The candidate uses `1e-7 rad^2` as a provisional simulation value (`sigma = 0.000316 rad`, 0.018 degrees); a real sensor's specified noise should replace it.

3. **The translational `R` was overconfident against this model's innovations.** Original `R_accel = 1e-5 (m/s^2)^2`; with the operating point corrected, the short-window innovations were about `0.0143` and `0.0107 m/s^2` for surge and heave. With the candidate process noise and other corrections, setting `R_accel = 1e-4 (m/s^2)^2` gave mean NIS `3.41` for 3 measurements over the first 30 s, close to the expected order of 3. Treat this as a simulation-calibrated starting point, not sensor calibration.

4. **`Q = Qc * dt` is not the main failure.** Exact Van Loan discretization differed from that approximation by only 0.53% in relative Frobenius norm for this model and timestep. The candidate nevertheless uses exact discretization. The wave filter in [Method.pdf](../Method.pdf) assumes unit white-noise spectral density, whereas the original sets its intensity to `0.2^2`; the candidate uses 1.0 to match that assumption. The force-bias terms are continuous force-noise spectral intensities, not acceleration variances. Their original magnitudes (`1e-2` etc.) are not physically interpretable as written; the candidate's `1e12` multiplier is only a provisional fit to this short simulation segment.

5. **The test setup leaks unmeasured truth into the filter.** At every step, radiation memory force is computed from `velocity_history`, i.e. OpenFAST's true velocities (all DOFs), even though the proposed measurement set contains only translational acceleration and pitch position. This is useful for a model-only offline check, but is not a causal implementation for those sensors. A deployable estimator needs velocity estimates and a causal memory-force update, or should include radiation states in the estimator.

6. **The linearization and truth runs are not using identical HydroDyn modes.** The script loads `main_lin.1.lin` from `main_lin.fst`/`HydroDyn_lin.dat`, whose current settings use `ExctnMod=1` and `RdtnMod=1`, but validates against `main.out` from `main.fst`/`HydroDyn.dat`, whose current settings use `ExctnMod=2` and `RdtnMod=2`. Those are DFT/convolution versus state-space modes. Align the runs or demonstrate equivalence before attributing residuals to Q/R.

7. **Absolute surge/heave position remains weakly constrained.** The measurement set has no direct surge or heave position/velocity observation; accelerometers alone are poor references for low-frequency displacement, especially when excitation and force-bias states are also uncertain. In the first 30 s candidate test, surge/heave position RMSE remained `0.116 m` / `0.0374 m`, versus truth standard deviations `0.0706 m` / `0.0261 m`. Pitch RMSE was `0.000393 rad` versus `0.000699 rad` truth standard deviation. This is a structural observability/measurement limitation, not something to solve by continuing to tune Q/R blindly.

8. **`P0 = 0.01 I` is not unit-consistent.** The state vector mixes metres, velocities, radiation-excitation states, wave states, and forces. Initial covariance should be defined per state family from initialization uncertainty; the current uniform diagonal remains a provisional placeholder in the candidate.

## Versioned work and test

- Candidate: `task100_usingKalman_v1.py`
- Reproducible short-window diagnostics: `_temp/diagnose_kalman_v1.py`
- Review of the first 1,201 samples (`dt = 0.025 s`, 30 s), using the existing true-velocity convolution as an offline input:
  - Candidate settings/state convention: `[0.116 m, 0.0374 m, 0.000393 rad]`, mean NIS `3.41`.
- This is a short-segment benchmark, not a full-record optimum or a sensor-noise calibration. The translational position error and truth-velocity dependency remain open issues.

## Recommended next steps

1. Make the linearization and output run use the same HydroDyn modes and operating point.
2. Decide which real sensors provide surge/heave low-frequency position or velocity; add those measurement rows if available.
3. Replace the truth-velocity radiation convolution with a causal estimate-based calculation or radiation-state model.
4. Set `R` from sensor specifications and estimate state-family `P0` / force-bias spectral intensities from aligned residuals over the full record.
