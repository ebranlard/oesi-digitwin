import pathlib
import sys

import numpy as np
from scipy.signal import fftconvolve


CODE_DIR = pathlib.Path(r"C:\Jonah\Uni\PhD\OESI\seastateEst\Code")
VERSION_DIR = CODE_DIR / "VS_Agent"
sys.path[:0] = [str(CODE_DIR), r"C:\Jonah\Uni\PhD\_libs\welib"]

source_path = VERSION_DIR / "task100_usingKalman_v3.py"
source = source_path.read_text(encoding="utf-8")
source = source.replace(
    "import notebook_setup # type: ignore\nnotebook_setup.enable_autoreload()\n", ""
)
setup_end = source.index("# --- Initialization & Execution Loop ---")
ns = {"__name__": "__main__"}
exec(compile(source[:setup_end], str(source_path), "exec"), ns)

from welib.kalman import kalman


time = ns["time"]
dt = ns["dt"]
velocity = ns["velocity_history"]
kernel_df = ns["dfHDsum"]
max_kernel = min(len(kernel_df), int(60 / dt))
kernel = kernel_df.iloc[:max_kernel, 2:38].to_numpy(dtype=float).reshape(
    -1, 6, 6
)

memory_force = np.zeros((6, len(time)))
for output_dof in range(6):
    for input_dof in range(6):
        memory_force[output_dof] -= fftconvolve(
            velocity[input_dof],
            kernel[:, output_dof, input_dof],
            mode="full",
        )[: len(time)]
memory_force *= dt
memory_force += 0.5 * dt * np.einsum("ij,jt->it", kernel[0], velocity)
for i in range(1, len(time)):
    oldest_lag = min(i, max_kernel - 1)
    memory_force[:, i] += (
        0.5 * dt * kernel[oldest_lag] @ velocity[:, i - oldest_lag]
    )
memory_force[:, 0] = -0.5 * dt * kernel[0] @ velocity[:, 0]

for i in [1, 100, 1200, max_kernel - 1, max_kernel, len(time) - 1]:
    reference = ns["FRad"].Conv_Rdtn(
        kernel_df, velocity[:, : i + 1], RdtnTMax=60
    ).flatten()
    np.testing.assert_allclose(
        memory_force[:, i], reference, rtol=2e-10, atol=2e-7
    )

active = ns["active_dof"]
memory_force = memory_force[active].T
memory_force_deviation = memory_force - memory_force[0]
excitation_force = ns["F_exc_history"]
excitation_force_deviation = excitation_force - excitation_force[0]

q_true = ns["displacement_history"][active].T
qdot_true = ns["velocity_history"][active].T
accel_true = ns["accel_history"][active].T
q_deviation = q_true - q_true[0]
qdot_deviation = qdot_true - qdot_true[0]
accel_deviation = accel_true - accel_true[0]

accel_model = (
    -q_deviation @ ns["inv_M_K"].T
    - qdot_deviation @ ns["inv_M_C"].T
    + (memory_force_deviation + excitation_force_deviation) @ ns["M_inv"].T
)
accel_residual = accel_deviation - accel_model
restoring_accel = -q_deviation @ ns["inv_M_K"].T
damping_accel = -qdot_deviation @ ns["inv_M_C"].T
memory_accel = memory_force_deviation @ ns["M_inv"].T
excitation_accel = excitation_force_deviation @ ns["M_inv"].T
heave_target_accel = (
    accel_deviation[:, 1]
    - damping_accel[:, 1]
    - memory_accel[:, 1]
    - excitation_accel[:, 1]
)
heave_k_fit = -np.dot(q_deviation[:, 1], heave_target_accel) / np.dot(
    q_deviation[:, 1], q_deviation[:, 1]
)

measurements = np.column_stack(
    [
        accel_deviation[:, 0],
        accel_deviation[:, 1],
        q_deviation[:, 2],
    ]
)
inputs = np.hstack((memory_force_deviation, excitation_force_deviation))

n_states = ns["A"].shape[0]
z = np.zeros(n_states)
P = np.eye(n_states) * 1e-2
P[-ns["n_dofs"] :, -ns["n_dofs"] :] = 0.0
estimated_states = np.zeros((len(time), n_states))
innovations = np.zeros((len(time), ns["n_dofs"]))
heave_idx = [name for name, enabled in ns["dof"].items() if enabled].index(
    "Heave"
)
heave_state_indices = (
    heave_idx,
    ns["n_dofs"] + heave_idx,
    2 * ns["n_dofs"] + heave_idx,
)
heave_gain = np.zeros((len(time), 3))
heave_correction = np.zeros((len(time), 3))

for i in range(1, len(time)):
    z_pred = ns["Xxd"] @ z + ns["Xud"] @ inputs[i]
    y_pred = ns["Yx"] @ z_pred + ns["Yu"] @ inputs[i]
    innovations[i] = measurements[i] - y_pred
    z, P, gain = kalman.EstimateKFTimeStep(
        inputs[i],
        measurements[i],
        z,
        ns["Xxd"],
        ns["Xud"],
        ns["Yx"],
        ns["Yu"],
        P,
        ns["Q"],
        ns["R"],
    )
    estimated_states[i] = z
    heave_gain[i] = gain[list(heave_state_indices), heave_idx]
    heave_correction[i] = heave_gain[i] * innovations[i, heave_idx]

estimated_q = (
    estimated_states[:, : ns["n_dofs"]] + ns["q_reference"]
)
heave_error = estimated_q[:, heave_idx] - q_true[:, heave_idx]
print("Samples/duration/dt:", len(time), time[-1], dt)
print("Initial heave reference [m]:", ns["q_reference"][heave_idx])
print("True heave std [m]:", np.std(q_true[:, heave_idx]))
print("Estimated heave std [m]:", np.std(estimated_q[:, heave_idx]))
print("Heave position RMSE [m]:", np.sqrt(np.mean(heave_error**2)))
print(
    "Acceleration model residual RMS [active DOFs]:",
    np.sqrt(np.mean(accel_residual**2, axis=0)),
)
print(
    "Acceleration deviation RMS [active DOFs]:",
    np.sqrt(np.mean(accel_deviation**2, axis=0)),
)
print("Restoring acceleration RMS:", np.sqrt(np.mean(restoring_accel**2, axis=0)))
print("Damping acceleration RMS:", np.sqrt(np.mean(damping_accel**2, axis=0)))
print("Radiation acceleration RMS:", np.sqrt(np.mean(memory_accel**2, axis=0)))
print("Excitation acceleration RMS:", np.sqrt(np.mean(excitation_accel**2, axis=0)))
print(
    "Heave restoring coefficient from .lin [1/s^2] and trajectory-fit [1/s^2]:",
    ns["inv_M_K"][heave_idx, heave_idx],
    heave_k_fit,
)
print(
    "Heave innovation RMS [m/s^2]:",
    np.sqrt(np.mean(innovations[:, heave_idx] ** 2)),
)
print(
    "Heave measurement gain RMS [q, qdot, bias]:",
    np.sqrt(np.mean(heave_gain**2, axis=0)),
)
print(
    "Heave measurement correction RMS [q, qdot, bias]:",
    np.sqrt(np.mean(heave_correction**2, axis=0)),
)
print("All state estimates finite:", np.isfinite(estimated_states).all())

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(11, 4))
ax.plot(time, q_true[:, heave_idx], color="black", label="OpenFAST")
ax.plot(time, estimated_q[:, heave_idx], color="crimson", label="KF")
ax.set_xlabel("Time [s]")
ax.set_ylabel("Heave [m]")
ax.grid(True, linestyle=":")
ax.legend()
fig.tight_layout()
plot_path = VERSION_DIR / "_temp" / "heave_incremental_validation_v1.png"
fig.savefig(plot_path, dpi=140)
print("Heave plot saved to", plot_path)
