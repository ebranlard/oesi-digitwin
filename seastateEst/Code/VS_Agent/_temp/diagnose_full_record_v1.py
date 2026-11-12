import contextlib
import io
import pathlib
import sys

import numpy as np
from scipy.signal import fftconvolve


CODE_DIR = pathlib.Path(r"C:\Jonah\Uni\PhD\OESI\seastateEst\Code")
VERSION_DIR = CODE_DIR / "VS_Agent"
OUTPUT_DIR = VERSION_DIR / "_temp"
sys.path[:0] = [str(CODE_DIR), r"C:\Jonah\Uni\PhD\_libs\welib"]

source_path = VERSION_DIR / "task100_usingKalman_v1.py"
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
velocity = ns["velocity_perturbation_history"]
kernel_df = ns["dfHDsum"]
max_kernel = min(len(kernel_df), int(60 / dt))
kernel = kernel_df.iloc[:max_kernel, 2:38].to_numpy(dtype=float).reshape(
    -1, 6, 6
)

# Compute the same truncated trapezoidal radiation convolution as FRad.Conv_Rdtn,
# but all at once so the full OpenFAST record can be checked efficiently.
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
    if not np.allclose(memory_force[:, i], reference, rtol=2e-10, atol=2e-7):
        print(
            "Convolution mismatch",
            i,
            "fast",
            memory_force[:, i],
            "reference",
            reference,
        )
    np.testing.assert_allclose(
        memory_force[:, i], reference, rtol=2e-10, atol=2e-7
    )

active = ns["active_dof"]
memory_force = memory_force[active]
measurements = np.column_stack(
    [
        ns["accel_history"][0],
        ns["accel_history"][2],
        ns["displacement_perturbation_history"][4],
    ]
)

z = np.zeros(ns["A"].shape[0])
active_names = [name for name, enabled in ns["dof"].items() if enabled]
z[active_names.index("Pitch")] = (
    ns["displacement_perturbation_history"][4, 0]
)
P = np.eye(len(z)) * 1e-2
estimated_states = np.zeros((len(time), len(z)))
estimated_states[0] = z
innovation_values = np.zeros((len(time), 3))

for i in range(1, len(time)):
    u = memory_force[:, i]
    y = measurements[i]
    z_pred = ns["Xxd"] @ z + ns["Xud"] @ u
    innovation_values[i] = y - (ns["Yx"] @ z_pred + ns["Yu"] @ u)
    z, P, _ = kalman.EstimateKFTimeStep(
        u, y, z, ns["Xxd"], ns["Xud"], ns["Yx"], ns["Yu"], P, ns["Q"], ns["R"]
    )
    estimated_states[i] = z

estimated_q = estimated_states[:, : ns["n_dofs"]].copy()
estimated_q += ns["q_operating_point"]
true_q = ns["displacement_history"][active].T
error = estimated_q - true_q
print(f"Full record: {len(time)} samples, duration={time[-1]:.1f}s, dt={dt}")
print("Position RMSE [surge, heave, pitch] =", np.sqrt(np.mean(error**2, axis=0)))
print("Position truth std =", np.std(true_q, axis=0))
print("Position estimate std =", np.std(estimated_q, axis=0))
print("Position estimate range =", np.ptp(estimated_q, axis=0))
print("Innovation RMS =", np.sqrt(np.mean(innovation_values[1:] ** 2, axis=0)))
print("All estimates finite =", np.isfinite(estimated_states).all())
print("Heave acceleration truth std =", np.std(ns["accel_history"][2]))
print("Acceleration reference offset =", ns["qddot_operating_point"])

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, axes = plt.subplots(2, 1, figsize=(12, 7), sharex=True)
axes[0].plot(time, true_q[:, 0], label="OpenFAST", color="black", linewidth=1)
axes[0].plot(time, estimated_q[:, 0], label="KF", color="crimson", linewidth=1)
axes[0].set_ylabel("Surge [m]")
axes[0].grid(True, linestyle=":")
axes[0].legend()
axes[1].plot(time, true_q[:, 1], label="OpenFAST", color="black", linewidth=1)
axes[1].plot(time, estimated_q[:, 1], label="KF", color="crimson", linewidth=1)
axes[1].set_ylabel("Heave [m]")
axes[1].set_xlabel("Time [s]")
axes[1].grid(True, linestyle=":")
axes[1].legend()
fig.tight_layout()
plot_path = OUTPUT_DIR / "surge_heave_full_record_v1.png"
fig.savefig(plot_path, dpi=140)
print("Plot saved to", plot_path)


def exact_qd(A, Qc, step):
    n = A.shape[0]
    block = np.block([[A, Qc], [np.zeros_like(A), -A.T]]) * step
    result = ns["expm"](block)
    Ad = result[:n, :n]
    Qd = result[:n, n:] @ Ad.T
    return 0.5 * (Qd + Qd.T)


base_force_intensity = np.array(
    [
        ns["var_f_surge"],
        ns["var_f_heave"],
        ns["var_f_pitch"],
    ]
)
for scale in [0.0, 1e6, 1e8, 1e10, 1e12]:
    W_test = np.diag(np.r_[ns["var_u_eta"], scale * base_force_intensity])
    Qc_test = ns["B_w"] @ W_test @ ns["B_w"].T
    Q_test = exact_qd(ns["A"], Qc_test, dt)
    z_test = np.zeros(ns["A"].shape[0])
    z_test[active_names.index("Pitch")] = (
        ns["displacement_perturbation_history"][4, 0]
    )
    P_test = np.eye(len(z_test)) * 1e-2
    q_test = np.zeros((len(time), 3))
    for i in range(1, len(time)):
        u = memory_force[:, i]
        z_test, P_test, _ = kalman.EstimateKFTimeStep(
            u,
            measurements[i],
            z_test,
            ns["Xxd"],
            ns["Xud"],
            ns["Yx"],
            ns["Yu"],
            P_test,
            Q_test,
            ns["R"],
        )
        q_test[i] = z_test[:3]
    q_test += ns["q_operating_point"]
    err_test = q_test - true_q
    print(
        f"bias scale {scale:.0e}: RMSE={np.sqrt(np.mean(err_test**2, axis=0))}; "
        f"estimate std={np.std(q_test, axis=0)}"
    )

raw_acceleration_measurements = measurements.copy()
raw_acceleration_measurements[1:, :2] += ns["qddot_operating_point"][
    [0, 1]
]
z_raw = np.zeros(ns["A"].shape[0])
z_raw[active_names.index("Pitch")] = (
    ns["displacement_perturbation_history"][4, 0]
)
P_raw = np.eye(len(z_raw)) * 1e-2
q_raw = np.zeros((len(time), 3))
for i in range(1, len(time)):
    z_raw, P_raw, _ = kalman.EstimateKFTimeStep(
        memory_force[:, i],
        raw_acceleration_measurements[i],
        z_raw,
        ns["Xxd"],
        ns["Xud"],
        ns["Yx"],
        ns["Yu"],
        P_raw,
        ns["Q"],
        ns["R"],
    )
    q_raw[i] = z_raw[:3]
q_raw += ns["q_operating_point"]
print(
    "Without subtracting linearization acceleration:",
    "RMSE=",
    np.sqrt(np.mean((q_raw - true_q) ** 2, axis=0)),
    "estimate std=",
    np.std(q_raw, axis=0),
)
