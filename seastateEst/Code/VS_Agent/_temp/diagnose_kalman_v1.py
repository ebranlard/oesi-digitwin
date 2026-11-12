import contextlib
import io
import pathlib
import sys

import numpy as np
from scipy.linalg import expm


CODE_DIR = pathlib.Path(r"C:\Jonah\Uni\PhD\OESI\seastateEst\Code")
VERSION_DIR = CODE_DIR / "VS_Agent"
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
from openfast_toolbox.io import FASTLinearizationFile


def exact_process_covariance(A, Qc, dt):
    n = A.shape[0]
    van_loan = np.block(
        [[A, Qc], [np.zeros_like(A), -A.T]]
    ) * dt
    result = expm(van_loan)
    Ad = result[:n, :n]
    Qd = result[:n, n:] @ Ad.T
    return Ad, 0.5 * (Qd + Qd.T)


def run_filter(label, Qd, R):
    n_steps = min(1201, len(ns["time"]))
    z = np.zeros(ns["A"].shape[0])
    active_names = [name for name, enabled in ns["dof"].items() if enabled]
    if "Pitch" in active_names:
        pitch_index = active_names.index("Pitch")
        z[pitch_index] = ns["displacement_history"][4, 0] - q_operating_point[2]
    P = np.eye(len(z)) * 1e-2
    q_true = (
        ns["displacement_history"][ns["active_dof"], :n_steps].T
        - q_operating_point
    )
    estimates = np.zeros_like(q_true)
    estimates[0] = z[:3]
    innovations = []
    nis_values = []

    for i in range(1, n_steps):
        u = memory_force[:, i]
        y = measurements[i]
        z_pred = ns["Xxd"] @ z + ns["Xud"] @ u
        P_pred = ns["Xxd"] @ P @ ns["Xxd"].T + Qd
        innovation = y - (ns["Yx"] @ z_pred + ns["Yu"] @ u)
        S = ns["Yx"] @ P_pred @ ns["Yx"].T + R
        nis_values.append(float(innovation @ np.linalg.solve(S, innovation)))
        innovations.append(innovation)
        z, P, _ = kalman.EstimateKFTimeStep(
            u, y, z, ns["Xxd"], ns["Xud"], ns["Yx"], ns["Yu"], P, Qd, R
        )
        estimates[i] = z[:3]

    innovations = np.asarray(innovations)
    error = estimates[1:] - q_true[1:]
    print(
        f"{label}: q_RMSE={np.sqrt(np.mean(error**2, axis=0))}; "
        f"truth_std={np.std(q_true[1:], axis=0)}; "
        f"innovation_RMS={np.sqrt(np.mean(innovations**2, axis=0))}; "
        f"mean_NIS={np.mean(nis_values):.4g}; "
        f"finite={np.isfinite(estimates).all()}"
    )


print(f"Samples={len(ns['time'])}, dt={ns['dt']}, states={ns['A'].shape[0]}")
print(
    "Initial true q [surge, heave, pitch rad]=",
    ns["displacement_history"][ns["active_dof"], 0],
)
print(
    "Measured signal std [surge acc, heave acc, pitch rad]=",
    np.std(
        np.vstack(
            [
                ns["accel_history"][0],
                ns["accel_history"][2],
                ns["displacement_history"][4],
            ]
        ),
        axis=1,
    ),
)
print("Configured measurement sigma=", np.sqrt(np.diag(ns["R"])))

Ad_exact, Q_exact = exact_process_covariance(
    ns["A"], ns["Q_continuous"], ns["dt"]
)
print(
    "Current Q=Qc*dt relative error vs Van Loan Qd=",
    np.linalg.norm(ns["Q"] - Q_exact) / np.linalg.norm(Q_exact),
)
print(
    "Exact Ad consistency error=",
    np.linalg.norm(Ad_exact - ns["Xxd"]) / np.linalg.norm(ns["Xxd"]),
)

n_steps = min(1201, len(ns["time"]))
linearization = FASTLinearizationFile(ns["filepath_KC_lin"])
state_descriptions = linearization["x_info"]["Description"]
state_operating_point = np.asarray(linearization["x"])
dof_keywords = ["surge", "sway", "heave", "roll", "pitch", "yaw"]
active_keywords = [
    keyword
    for keyword, active in zip(dof_keywords, ns["active_dof"])
    if active
]
position_operating_point = []
velocity_operating_point = []
for keyword in active_keywords:
    position_index = next(
        i
        for i, description in enumerate(state_descriptions)
        if description.strip().startswith("ED Platform")
        and keyword in description.lower()
    )
    velocity_index = next(
        i
        for i, description in enumerate(state_descriptions)
        if description.strip().startswith("ED First time derivative")
        and keyword in description.lower()
    )
    position_operating_point.append(state_operating_point[position_index])
    velocity_operating_point.append(state_operating_point[velocity_index])
q_operating_point = np.asarray(position_operating_point)
velocity_operating_point = np.asarray(velocity_operating_point)
velocity_perturbation_history = ns["velocity_history"].copy()
velocity_perturbation_history[ns["active_dof"]] -= velocity_operating_point[:, None]
print("Linearization q operating point=", q_operating_point)
print("Linearization qdot operating point=", velocity_operating_point)
print(
    "Initial perturbation q=",
    ns["displacement_history"][ns["active_dof"], 0] - q_operating_point,
)

memory_force = np.zeros((ns["n_dofs"], n_steps))
with contextlib.redirect_stdout(io.StringIO()):
    for i in range(1, n_steps):
        memory_force[:, i] = ns["FRad"].Conv_Rdtn(
            ns["dfHDsum"],
            velocity_perturbation_history[:, : i + 1],
            RdtnTMax=60,
        ).flatten()[ns["active_dof"]]

measurements = np.column_stack(
    [
        ns["accel_history"][0, :n_steps],
        ns["accel_history"][2, :n_steps],
        ns["displacement_history"][4, :n_steps] - q_operating_point[2],
    ]
)

run_filter("current settings", ns["Q"], ns["R"])
run_filter("exact Q, current W/R", Q_exact, ns["R"])

R_accel_x0_1 = ns["R"].copy()
R_accel_x0_1[:2, :2] *= 0.1
run_filter("exact Q, acceleration R x0.1", Q_exact, R_accel_x0_1)

R_accel_x10 = ns["R"].copy()
R_accel_x10[:2, :2] *= 10.0
run_filter("exact Q, acceleration R x10", Q_exact, R_accel_x10)

W_unit_eta = ns["W"].copy()
W_unit_eta[0, 0] = 1.0
Qc_unit_eta = ns["B_w"] @ W_unit_eta @ ns["B_w"].T
_, Q_unit_eta = exact_process_covariance(ns["A"], Qc_unit_eta, ns["dt"])
run_filter("exact Q, unit wave-noise intensity", Q_unit_eta, ns["R"])

R_tighter_pitch = ns["R"].copy()
R_tighter_pitch[2, 2] = 1e-7
run_filter("exact Q/unit wave, pitch R=1e-7 rad^2", Q_unit_eta, R_tighter_pitch)

W_large_bias = W_unit_eta.copy()
W_large_bias[1:, 1:] *= 1e-2
Qc_large_bias = ns["B_w"] @ W_large_bias @ ns["B_w"].T
_, Q_large_bias = exact_process_covariance(ns["A"], Qc_large_bias, ns["dt"])
run_filter(
    "exact Q/unit wave, bias spectral intensity x1e10",
    Q_large_bias,
    R_tighter_pitch,
)
run_filter(
    "bias intensity x1e10, acceleration R x10",
    Q_large_bias,
    R_accel_x10,
)

W_larger_bias = W_unit_eta.copy()
W_larger_bias[1:, 1:] *= 1e2
Qc_larger_bias = ns["B_w"] @ W_larger_bias @ ns["B_w"].T
_, Q_larger_bias = exact_process_covariance(ns["A"], Qc_larger_bias, ns["dt"])
run_filter(
    "exact Q/unit wave, bias spectral intensity x1e14",
    Q_larger_bias,
    R_tighter_pitch,
)
run_filter(
    "bias intensity x1e14, acceleration R x10",
    Q_larger_bias,
    R_accel_x10,
)
