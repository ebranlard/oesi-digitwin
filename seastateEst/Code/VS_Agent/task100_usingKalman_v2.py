#%%
# Enable Autoreload
import notebook_setup # type: ignore
notebook_setup.enable_autoreload()

# General packages
import numpy as np
import matplotlib.pyplot as plt
import os

# OpenFAST toolbox
from openfast_toolbox.io import FASTOutputFile # type: ignore
from openfast_toolbox.io import FASTInputFile # type: ignore

# Welib
from welib.kalman import kalman

# Seastate Estimator
from seastateEst import FRad
from seastateEst import FExc
from seastateEst import linData
from seastateEst import kalman_helper
#%%

# =============================================================================
# --- Parameters defining the system and simulation ---
# =============================================================================
# Set the DOFs
dof = {
    "Surge": True,
    "Sway": False,
    "Heave": True,
    "Roll": False,
    "Pitch": True,
    "Yaw": False,
}
active_dof = list(dof.values())
n_dofs = sum(dof.values())

# --- Structure mass
# filepath_mass = r"C:\Jonah\Uni\PhD\OESI\seastateEst\OpenFAST_cases\IEA-15_Full_LINEAR_Falling\main.1.lin" #TODO fix Ainf to account for active_dof, M already does this
# M_s, __ = linData.rigidBodyMass(filepath_mass)
# # Added mass at infinite frequency
# filepath_Ainf = r"C:\Jonah\Uni\PhD\OESI\seastateEst\OpenFAST_cases\IEA-15_Full\HydroData\IEA-15-240-RWT-UMaineSemi.1"
# A_inf = FRad.Ainf(filepath_Ainf, WtrDens=1025)
# # Total mass
# M = M_s + A_inf
# # Filter for active dofs
# M = M[np.ix_(active_dof, active_dof)]
# # Compute inverse mass
# M_inv = np.linalg.inv(M)
M, M_inv = linData.rigidBodyMass(r'C:\Jonah\Uni\PhD\OESI\seastateEst\OpenFAST_cases\IEA-15_PlatformOnly_v5.0\main_lin.1.lin')

# --- Stiffness (K) and damping (C)
# filepath_KC = r'C:\Jonah\Uni\PhD\OESI\seastateEst\OpenFAST_cases\IEA-15_PlatformOnly_LINEAR\No_waves' #TODO chose from what OpenFAST simulation I really want K and C
# filepath_KC_lin = str(os.path.join(filepath_KC, 'main.1.lin'))
filepath_KC = r'C:\Jonah\Uni\PhD\OESI\seastateEst\OpenFAST_cases\IEA-15_PlatformOnly_v5.0'
filepath_KC_lin = str(os.path.join(filepath_KC, 'main_lin.1.lin'))

inv_M_K = linData.invMK(filepath_KC_lin, dof)   # M^(-1)*K
inv_M_C = linData.invMC(filepath_KC_lin, dof)   # M^(-1)*C
# Sanity check if same DoFs are enabled as in this python file as in the OpenFAST ED file
filepath_KC_ED = os.path.join(filepath_KC, 'ED.dat')
ED_file = FASTInputFile(filepath_KC_ED)
assert dof["Surge"] == ED_file["PtfmSgDOF"], "Surge DOF mismatch!"
assert dof["Sway"]  == ED_file["PtfmSwDOF"], "Sway DOF mismatch!"
assert dof["Heave"] == ED_file["PtfmHvDOF"], "Heave DOF mismatch!"
assert dof["Roll"]  == ED_file["PtfmRDOF"],  "Roll DOF mismatch!"
assert dof["Pitch"] == ED_file["PtfmPDOF"],  "Pitch DOF mismatch!"
assert dof["Yaw"]   == ED_file["PtfmYDOF"],  "Yaw DOF mismatch!"

# --- Load measurements from OpenFAST
filepath_OF = r'C:\Jonah\Uni\PhD\OESI\seastateEst\OpenFAST_cases\IEA-15_PlatformOnly_v5.0'
dfOut = FASTOutputFile(os.path.join(filepath_OF, 'main.out')).toDataFrame()
dfHDsum = FRad.RdtnKrnl(os.path.join(filepath_OF, 'main.HD.sum'))

# Velocity and acceleration channels in dfOut
displacement_channels = of_displs_full = ['PtfmSurge_[m]', 'PtfmSway_[m]', 'PtfmHeave_[m]', 'PtfmRoll_[deg]', 'PtfmPitch_[deg]', 'PtfmYaw_[deg]']
velocity_channels = ['PtfmTVxi_[m/s]', 'PtfmTVyi_[m/s]', 'PtfmTVzi_[m/s]', 'PtfmRVxi_[deg/s]', 'PtfmRVyi_[deg/s]', 'PtfmRVzi_[deg/s]']
accel_channels = ['PtfmTAxi_[m/s^2]', 'PtfmTAyi_[m/s^2]', 'PtfmTAzi_[m/s^2]', 'PtfmRAxi_[deg/s^2]', 'PtfmRAyi_[deg/s^2]', 'PtfmRAzi_[deg/s^2]']

# Extract values from dfOut
displacement_history = dfOut[displacement_channels].values.T
displacement_history[3:6, :] = np.radians(displacement_history[3:6, :]) # Convert Roll, Pitch and Yaw displacements from deg to rad
velocity_history = dfOut[velocity_channels].values.T
velocity_history[3:6, :] = np.radians(velocity_history[3:6, :]) # Convert Roll, Pitch and Yaw velocities from deg/s to rad/s
accel_history = dfOut[accel_channels].values.T
accel_history[3:6, :] = np.radians(accel_history[3:6, :])
time = dfOut['Time_[s]'].to_numpy()

# --- Precompute excitation force from the OpenFAST wave-elevation history
seastate_file = FASTInputFile(os.path.join(filepath_OF, 'SeaState.dat'))
filepath_WAMIT_excitation = os.path.join(filepath_OF, 'HydroData', 'marin_semi.3')
K_exc, omega_exc = FExc.ExctnKrnl(
    filepath_WAMIT_excitation,
    beta=seastate_file["WaveDir"],
)
eta_t = dfOut['Wave1Elev_[m]'].to_numpy()
eta_omega, omega_fft = FExc.complex_WaveElev(eta_t, time)
dt = time[1] - time[0]
F_exc_history_full = FExc.IFFT_Exctn(
    K_exc, omega_exc, eta_omega, omega_fft, dt
)
if F_exc_history_full.shape != (len(time), 6):
    raise ValueError(
        f"Expected excitation history shape {(len(time), 6)}, "
        f"got {F_exc_history_full.shape}"
    )
F_exc_history = F_exc_history_full[:, active_dof]

# Subtract mean from heave
# displacement_history[2, :] -= np.mean(displacement_history[2, :])
# accel_history[2, :] -= np.mean(accel_history[2, :])

# =============================================================================
# --- Build System matrices ---
# =============================================================================

A = np.block([
    [np.zeros((n_dofs, n_dofs)), np.eye(n_dofs), np.zeros((n_dofs, n_dofs))],
    [-inv_M_K, -inv_M_C, M_inv],
    [np.zeros((n_dofs, n_dofs)), np.zeros((n_dofs, n_dofs)), np.zeros((n_dofs, n_dofs))]
])

B = np.block([
    [np.zeros((n_dofs, 2 * n_dofs))],
    [M_inv, M_inv],
    [np.zeros((n_dofs, 2 * n_dofs))]
])

# --- Error levels of the measurements
# R_diags = [1e-5, 1e-5, 1e-5, 1e-3, 1e-3, 1e-3]
# R_diags = np.array(R_diags)[active_dof]         # Filter for active DoFs
# R = np.diag(R_diags)
sigma_accel = 1e-2
sigma_pitch = np.sqrt(1e-7)
R_diags = [
    sigma_accel**2,
    sigma_accel**2,
    sigma_accel**2,
    1e-3,
    sigma_pitch**2,
    1e-3,
]
R_diags = np.array(R_diags)[active_dof]         # Filter for active DoFs
R = np.diag(R_diags)

# --- Error levels of the model
# Noise mapping matrix
B_w = np.block([
    [np.zeros((n_dofs, n_dofs))],
    [np.zeros((n_dofs, n_dofs))],
    [np.eye(n_dofs)]
])

# Continuous white-noise spectral intensities (not discrete per-sample variances).
var_f_surge = 1e-2
var_f_sway = 1e-4
var_f_heave = 1e-2        
var_f_roll = 1e-2
var_f_pitch = 1e-1        
var_f_yaw = 1e-2
bias_noise_scale = 1e12   # Provisional 30 s simulation calibration; not sensor-independent

# Construct the diagonal elements for W and apply active_dof filter
W_diags_full = [var_f_surge, var_f_sway, var_f_heave, var_f_roll, var_f_pitch, var_f_yaw]
W_diags = list(np.array(W_diags_full)[active_dof] * bias_noise_scale)
W = np.diag(W_diags)

# Project into full state space process noise
Q_continuous = B_w @ W @ B_w.T

# =============================================================================
# --- Discretization & Measurement Mapping ---
# =============================================================================

# 1. Compute time step from OpenFAST time array
dt = time[1] - time[0]

# 2. Discretize the state matrices
# Xxd, Xud = kalman.KFDiscretize(A, B, dt, method='exponential')
Xxd, Xud = kalman_helper.discretize_exact(A, B, dt)

# 3. Scale process noise to discrete-time
Q = Q_continuous * dt

# 4. Both force inputs enter acceleration through M_inv.
Yx, Yu_one_force = kalman_helper.build_mixed_measurement_matrices(
    A, B[:, :n_dofs], dof
)
Yu = np.hstack((Yu_one_force, Yu_one_force))
# =============================================================================
# --- Initialization & Execution Loop ---
# =============================================================================

# Initialize states and state error covariance matrix
n_states = A.shape[0]
z0 = np.zeros(n_states)          # Initial guess for all states (zeros)
z0[1] = -0.44                   # Initial heave displacement
P0 = np.eye(n_states) * 1e-2     # Initial state uncertainty baseline

# Pre-allocate array to store estimated states over time
estimated_states_history = np.zeros((len(time), n_states))
estimated_states_history[0, :] = z0

innovation_history = np.zeros((len(time), n_dofs))  #TODO: This line can be deleted, just for debugging
innovation_nis_per_dof = np.full(len(time), np.nan)
gain_heave_measurement = np.zeros((len(time), n_states)) #TODO: This line can be deleted, just for debugging

# Execution loop over the OpenFAST time-series
for i in range(1, len(time)):
    
    # --- Calculate Radiation Memory Force (F_mem)
    current_velocity_slice = velocity_history[:, :i+1]

    # Compute Fmem using only the historical data available up to step n
    F_mem = FRad.Conv_Rdtn(dfHDsum, current_velocity_slice, RdtnTMax=60)[active_dof]
    
    # 2. Extract inputs and measurements for current time step n
    u1 = np.concatenate((F_mem.flatten(), F_exc_history[i]))
    # y1 = accel_history[:, i][active_dof]    # Measured accelerations (shape: n_dofs,)

    # Dynamically Assemble Mixed Vector y1 for current step i
    y1_elements = []
    for idx in range(6):
        if active_dof[idx]:
            # If it is an active translational DOF (Surge=0, Heave=2), extract ACCELERATION
            if idx < 3:
                y1_elements.append(accel_history[idx, i])
            # If it is an active rotational DOF (Pitch=4), extract INCLINATION POSITION
            else:
                y1_elements.append(displacement_history[idx, i])
                
    y1 = np.array(y1_elements)            # Mixed measurement vector (shape: 3,)

    # Predicted state and measurement before the KF correction
    z_pred = Xxd @ z0 + Xud @ u1
    y_pred = Yx @ z_pred + Yu @ u1

    innovation = y1 - y_pred #TODO: This line can be deleted, just for debugging
    P_pred = Xxd @ P0 @ Xxd.T + Q
    innovation_covariance = Yx @ P_pred @ Yx.T + R
    innovation_nis_per_dof[i] = (
        innovation @ np.linalg.solve(innovation_covariance, innovation)
    ) / n_dofs
    
    # 3. Process time-step through the Kalman Filter engine
    z0, P0, Kk = kalman.EstimateKFTimeStep(u1, y1, z0, Xxd, Xud, Yx, Yu, P0, Q, R)

    innovation_history[i, :] = innovation #TODO: This line can be deleted, just for debugging
    gain_heave_measurement[i, :] = Kk[:, 1] #TODO: This line can be deleted, just for debugging
    
    # 4. Store the updated states
    estimated_states_history[i, :] = z0

print("Kalman Filter execution complete.")

# Print a tab-separated summary row for the comparison document.
estimated_q_for_metrics = estimated_states_history[:, :n_dofs]
true_q_for_metrics = displacement_history[[0, 2, 4], :].T
estimated_q_increment = estimated_q_for_metrics - estimated_q_for_metrics[0]
true_q_increment = true_q_for_metrics - true_q_for_metrics[0]
position_rmse = np.sqrt(np.mean((estimated_q_increment - true_q_increment) ** 2, axis=0))
position_correlation = np.array([
    np.corrcoef(estimated_q_increment[:, j], true_q_increment[:, j])[0, 1]
    if np.std(estimated_q_increment[:, j]) > 0 and np.std(true_q_increment[:, j]) > 0
    else np.nan
    for j in range(n_dofs)
])
valid_innovations = innovation_history[1:]
innovation_mean = np.mean(valid_innovations, axis=0)
innovation_rms = np.sqrt(np.mean(valid_innovations ** 2, axis=0))
force_bias_rms = np.sqrt(np.mean(estimated_states_history[:, -n_dofs:] ** 2, axis=0))

def _format_triplet(values):
    return "/".join(f"{value:.4g}" for value in values)

print("\nCOPY-PASTE SUMMARY (grouped values are Surge/Heave/Pitch):")
print("Run\tQ/R settings\tPosition RMSE [m/m/rad]\tPosition correlation\tInnovation mean [m/s^2,m/s^2,rad]\tInnovation RMS [m/s^2,m/s^2,rad]\tNIS/DOF\tForce-bias RMS [N/N/Nm]\tWave RMSE [m]")
print("\t".join([
    "Prescribed F_exc (v2)",
    f"sig_acc/heave/pitch={sigma_accel:.4g}/{sigma_accel:.4g}/{sigma_pitch:.4g};bias_scale={bias_noise_scale:.4g};var_f={var_f_surge:.4g}/{var_f_heave:.4g}/{var_f_pitch:.4g}",
    _format_triplet(position_rmse),
    _format_triplet(position_correlation),
    _format_triplet(innovation_mean),
    _format_triplet(innovation_rms),
    f"{np.nanmean(innovation_nis_per_dof):.4g}",
    _format_triplet(force_bias_rms),
    "N/A",
]))

# %%
# =============================================================================
# --- Plotting 3-DOF Platform Displacements (Surge, Heave, Pitch) ---
# =============================================================================

# Define labels, units, and explicit names matching your active 3 DOFs
dof_names = ["Surge", "Heave", "Pitch"]
units = ["m", "m", "rad"]
n_active_dofs = len(dof_names)

# Extract estimated positions (First 3 columns of your state history matrix)
estimated_q = estimated_states_history[:, 0:n_active_dofs]

# Extract corresponding raw true channels from OpenFAST for validation
of_displ_channels = ['PtfmSurge_[m]', 'PtfmHeave_[m]', 'PtfmPitch_[deg]']

if all(ch in dfOut.columns for ch in of_displ_channels):
    # If you filtered dfOut earlier, make sure it matches the full time axis.
    # Otherwise, extract the values directly:
    of_true_q = dfOut[of_displ_channels].values
    of_true_q[:, 2] = np.radians(of_true_q[:, 2]) # Convert Pitch from deg to rad
    plot_comparison = True
else:
    plot_comparison = False

# Construct a clean 3x1 vertical stack layout
fig, axes = plt.subplots(n_active_dofs, 1, figsize=(10, 8), sharex=True)

for idx in range(n_active_dofs):
    ax = axes[idx]
    
    # Plot Kalman Filter State Estimates
    ax.plot(time, estimated_q[:, idx], label="KF Estimate", color="crimson", linewidth=1.5)
    
    # Overlay OpenFAST baseline truths if available
    if plot_comparison:
        ax.plot(time, of_true_q[:, idx], label="OpenFAST Truth", color="black", linestyle="--", alpha=0.7)
        
    ax.set_title(f"Platform {dof_names[idx]}")
    ax.set_ylabel(f"Displacement [{units[idx]}]")
    ax.grid(True, linestyle=":", alpha=0.6)
    # ax.set_xlim([50, 300])
    
    if idx == 0: 
        ax.legend(loc="upper right")

# Set shared X-axis label only on the bottom-most subplot (Pitch)
axes[-1].set_xlabel("Time [s]")

plt.tight_layout()
plt.show()

#%%
# plt.figure()
# plt.plot(time, estimated_states_history[:, 41], label="KF Estimate", color="crimson", linewidth=1.5)

# %%
