#%%
# Enable Autoreload
import notebook_setup # type: ignore
notebook_setup.enable_autoreload()

# General packages
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.linalg import expm

# OpenFAST toolbox
from openfast_toolbox.io import FASTOutputFile # type: ignore
from openfast_toolbox.io import FASTInputFile # type: ignore
from openfast_toolbox.io import FASTLinearizationFile # type: ignore

# Welib
from welib.kalman import kalman

# Seastate Estimator
from seastateEst import FRad
from seastateEst import linData
from seastateEst import read_ssFiles
from seastateEst import wave_filter
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

# The linearized state is a perturbation about the operating point in the .lin file.
linearization = FASTLinearizationFile(filepath_KC_lin)
linear_state_descriptions = linearization["x_info"]["Description"]
linear_state_operating_point = np.asarray(linearization["x"])
linear_state_derivative = np.asarray(linearization["xdot"])
dof_keywords = ["surge", "sway", "heave", "roll", "pitch", "yaw"]
position_state_indices = []
velocity_state_indices = []
for keyword, is_active in zip(dof_keywords, active_dof):
    if not is_active:
        continue
    position_state_indices.append(
        next(
            i
            for i, description in enumerate(linear_state_descriptions)
            if description.strip().startswith("ED Platform")
            and keyword in description.lower()
        )
    )
    velocity_state_indices.append(
        next(
            i
            for i, description in enumerate(linear_state_descriptions)
            if description.strip().startswith("ED First time derivative")
            and keyword in description.lower()
        )
    )

q_operating_point = linear_state_operating_point[position_state_indices]
qdot_operating_point = linear_state_operating_point[velocity_state_indices]
qddot_operating_point = linear_state_derivative[velocity_state_indices]

# Sanity check if same DoFs are enabled as in this python file as in the OpenFAST ED file
filepath_KC_ED = os.path.join(filepath_KC, 'ED.dat')
ED_file = FASTInputFile(filepath_KC_ED)
assert dof["Surge"] == ED_file["PtfmSgDOF"], "Surge DOF mismatch!"
assert dof["Sway"]  == ED_file["PtfmSwDOF"], "Sway DOF mismatch!"
assert dof["Heave"] == ED_file["PtfmHvDOF"], "Heave DOF mismatch!"
assert dof["Roll"]  == ED_file["PtfmRDOF"],  "Roll DOF mismatch!"
assert dof["Pitch"] == ED_file["PtfmPDOF"],  "Pitch DOF mismatch!"
assert dof["Yaw"]   == ED_file["PtfmYDOF"],  "Yaw DOF mismatch!"

# --- SS_EXCTN matrices
# filepath_ssexctn = r'C:\Jonah\Uni\PhD\OESI\seastateEst\OpenFAST_cases\Parametric_study\Lock_Rad\HydroData\marin_semi_Fit0.2_2_2_2_2_2_4.ssexctn'
filepath_ssexctn = r'C:\Jonah\Uni\PhD\OESI\seastateEst\OpenFAST_cases\IEA-15_PlatformOnly_v5.0\HydroData\marin_semi.ssexctn'
A_exc, B_exc, C_exc, n_exc, m_exc = read_ssFiles.read_ssexctn(filepath_ssexctn)

# Reduce Exctn matrices depening on active dofs
A_exc, B_exc, C_exc, n_exc = read_ssFiles.reduceExctn(dof, A_exc, B_exc, C_exc, n_exc, m_exc)

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

# Express measurements and radiation-convolution velocities as perturbations.
displacement_perturbation_history = displacement_history.copy()
displacement_perturbation_history[active_dof] -= q_operating_point[:, None]
velocity_perturbation_history = velocity_history.copy()
velocity_perturbation_history[active_dof] -= qdot_operating_point[:, None]
accel_history[active_dof] -= qddot_operating_point[:, None]

# --- Wave shaping filter
seastate_file = FASTInputFile(os.path.join(filepath_OF, 'SeaState.dat')) # Match the parameters to the OF simulation
hs = seastate_file["WaveHs"]                # Significant wave height
tp = seastate_file["WaveTp"]                # Peak spectral period
zeta = 0.1

wave_Shapefilter = dict()
wave_Shapefilter['S_eta(f)'], wave_Shapefilter['k_w'], wave_Shapefilter['omega_0'] = wave_filter.waveFilter2Ord(zeta, hs, tp)
A_eta = np.block([[0, 1], [-wave_Shapefilter['omega_0']**2, -2 * zeta * wave_Shapefilter['omega_0']]])
B_eta = np.array([[0], [wave_Shapefilter['k_w']]])
C_eta = np.array([[1, 0]])
n_eta = A_eta.shape[0]

# =============================================================================
# --- Build System matrices ---
# =============================================================================

A = np.block([
    [np.zeros((n_dofs, n_dofs)), np.eye(n_dofs),             np.zeros((n_dofs, n_exc)), np.zeros((n_dofs, n_eta)), np.zeros((n_dofs, n_dofs))],
    [-inv_M_K,                   -inv_M_C,                   M_inv @ C_exc,             np.zeros((n_dofs, n_eta)), M_inv                     ],
    [np.zeros((n_exc, n_dofs)),  np.zeros((n_exc, n_dofs)),  A_exc,                     B_exc @ C_eta,             np.zeros((n_exc, n_dofs)) ],
    [np.zeros((n_eta, n_dofs)),  np.zeros((n_eta, n_dofs)),  np.zeros((n_eta, n_exc)),  A_eta,                     np.zeros((n_eta, n_dofs)) ],
    [np.zeros((n_dofs, n_dofs)), np.zeros((n_dofs, n_dofs)), np.zeros((n_dofs, n_exc)), np.zeros((n_dofs, n_eta)), np.zeros((n_dofs, n_dofs))]
])

B = np.block([
    [np.zeros((n_dofs, n_dofs))],
    [M_inv],
    [np.zeros((n_exc, n_dofs))],
    [np.zeros((n_eta, n_dofs))],
    [np.zeros((n_dofs, n_dofs))]
])

# --- Measurement variance (accelerations in m/s^2; pitch position in rad)
# Provisional values from the first 30 s OpenFAST consistency test; replace the
# measurement variances with sensor specifications for deployment.
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
    [np.zeros((n_dofs, 1)), np.zeros((n_dofs, n_dofs))],
    [np.zeros((n_dofs, 1)), np.zeros((n_dofs, n_dofs))],
    [np.zeros((n_exc, 1)),  np.zeros((n_exc, n_dofs)) ],
    [B_eta,                 np.zeros((n_eta, n_dofs)) ],
    [np.zeros((n_dofs, 1)), np.eye(n_dofs)            ]
])

# Continuous white-noise spectral intensities (not discrete per-sample variances).
var_u_eta = 1.0          # Unit white-noise intensity assumed by the wave-spectrum fit
var_f_surge = 1e-2
var_f_sway = 1e-4
var_f_heave = 1e-2        
var_f_roll = 1e-2
var_f_pitch = 1e-1        
var_f_yaw = 1e-2
bias_noise_scale = 1e12   # Provisional 30 s simulation calibration; not sensor-independent

# Construct the diagonal elements for W and apply active_dof filter
W_diags_full = [var_u_eta, var_f_surge, var_f_sway, var_f_heave, var_f_roll, var_f_pitch, var_f_yaw] 
W_diags = [W_diags_full[0]] + list(np.array(W_diags_full[1:])[active_dof] * bias_noise_scale)
W = np.diag(W_diags)
# W_diags = [var_u_eta] + [var_f_trans]*3 + [var_f_rot]*3
# W_diags = [W_diags[0]] + list(np.array(W_diags[1:])[active_dof])    # Filter for active DoFs
# W = np.diag(W_diags)

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

# 3. Discretize the continuous process-noise covariance using Van Loan's method.
n_states = A.shape[0]
van_loan = np.block([
    [A, Q_continuous],
    [np.zeros_like(A), -A.T],
]) * dt
van_loan_exp = expm(van_loan)
Ad_q = van_loan_exp[:n_states, :n_states]
Q = van_loan_exp[:n_states, n_states:] @ Ad_q.T
Q = 0.5 * (Q + Q.T)

# 4. Construct Measurement Matrices (Yx and Yu)
Yx, Yu = kalman_helper.build_mixed_measurement_matrices(A, B, dof)
# Yx = A[n_dofs : 2 * n_dofs, :]
# Yu = B[n_dofs : 2 * n_dofs, :]
# Yx = np.block([
#     [-inv_M_K, -inv_M_C, M_inv @ C_exc, np.zeros((n_dofs, n_eta)), np.zeros((n_dofs, n_dofs))]
# ])

# Yu = M_inv  # Direct throughput mapping of F_mem on physical acceleration


# =============================================================================
# --- Initialization & Execution Loop ---
# =============================================================================

# Initialize states and state error covariance matrix
n_states = A.shape[0]
z0 = np.zeros(n_states)          # Initial guess for all states (zeros)
active_names = [name for name, is_active in dof.items() if is_active]
if "Pitch" in active_names:
    pitch_state_index = active_names.index("Pitch")
    z0[pitch_state_index] = displacement_perturbation_history[4, 0]
P0 = np.eye(n_states) * 1e-2     # Initial state uncertainty baseline

# Pre-allocate array to store estimated states over time
estimated_states_history = np.zeros((len(time), n_states))
estimated_states_history[0, :] = z0

innovation_history = np.zeros((len(time), n_dofs))  #TODO: This line can be deleted, just for debugging
gain_heave_measurement = np.zeros((len(time), n_states)) #TODO: This line can be deleted, just for debugging

# Execution loop over the OpenFAST time-series
for i in range(1, len(time)):
    
    # --- Calculate Radiation Memory Force (F_mem)
    current_velocity_slice = velocity_perturbation_history[:, :i+1]

    # Compute Fmem using only the historical data available up to step n
    F_mem = FRad.Conv_Rdtn(dfHDsum, current_velocity_slice, RdtnTMax=60)[active_dof]
    
    # 2. Extract inputs and measurements for current time step n
    u1 = F_mem.flatten()                  # Input vector (shape: n_dofs,)
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
                y1_elements.append(displacement_perturbation_history[idx, i])
                
    y1 = np.array(y1_elements)            # Mixed measurement vector (shape: 3,)

    # Predicted state and measurement before the KF correction
    z_pred = Xxd @ z0 + Xud @ u1
    y_pred = Yx @ z_pred + Yu @ u1

    innovation = y1 - y_pred #TODO: This line can be deleted, just for debugging
    
    # 3. Process time-step through the Kalman Filter engine
    z0, P0, Kk = kalman.EstimateKFTimeStep(u1, y1, z0, Xxd, Xud, Yx, Yu, P0, Q, R)

    innovation_history[i, :] = innovation #TODO: This line can be deleted, just for debugging
    gain_heave_measurement[i, :] = Kk[:, 1] #TODO: This line can be deleted, just for debugging
    
    # 4. Store the updated states
    estimated_states_history[i, :] = z0

print("Kalman Filter execution complete.")

# %%
# =============================================================================
# --- Plotting 3-DOF Platform Displacements (Surge, Heave, Pitch) ---
# =============================================================================

# Define labels, units, and explicit names matching your active 3 DOFs
dof_names = ["Surge", "Heave", "Pitch"]
units = ["m", "m", "rad"]
n_active_dofs = len(dof_names)

# Extract estimated positions (First 3 columns of your state history matrix)
estimated_q = estimated_states_history[:, 0:n_active_dofs] + q_operating_point

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

# =============================================================================
# --- Plotting Wave Elevation (Sea-State Tracking) ---
# =============================================================================
estimated_wave_elevation = estimated_states_history[:, 2 * n_dofs + n_exc]

# Match the styling, sizing, and colors of your structural figures
plt.figure(figsize=(10, 3))
plt.plot(time, estimated_wave_elevation, label="KF Estimate", color="crimson", linewidth=1.5)
plt.plot(time, dfOut['Wave1Elev_[m]'].to_numpy(), label="OpenFAST Truth", color="black", linestyle="--", alpha=0.7)

# Match structural plot labeling and grid style
plt.title("Wave Elevation")
plt.ylabel("Displacement [m]")
plt.xlabel("Time [s]")
plt.grid(True, linestyle=":", alpha=0.6)
plt.legend(loc="upper right")

plt.tight_layout()
plt.show()

#%%
# plt.figure()
# plt.plot(time, estimated_states_history[:, 41], label="KF Estimate", color="crimson", linewidth=1.5)

# %%
