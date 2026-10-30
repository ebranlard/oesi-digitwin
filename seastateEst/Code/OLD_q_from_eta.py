""" 
Get q from a prescribed eta with a time integration scheme
"""
# Enable Autoreload
import notebook_setup # type: ignore
notebook_setup.enable_autoreload()

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# OpenFAST toolbox
from openfast_toolbox.io import FASTOutputFile # type: ignore
from openfast_toolbox.io import FASTLinearizationFile # type: ignore
# Welib
from welib.system.mech_system import MechSystem
from welib.system.statespacelinear import LinearStateSpace
from welib.tools.signal_analysis import correlated_signal, step
from scipy.interpolate import interp1d
# Seastate Estimator
from seastateEst import linData
from seastateEst import read_ssFiles
from seastateEst import readWAMIT


# --- Parameters defining the system and simulation
# Set the DOFs
dof = {
    "Surge": False,
    "Sway": False,
    "Heave": True,
    "Roll": False,
    "Pitch": False,
    "Yaw": False,
}
n_dofs = sum(dof.values())

# Structure mass
filepath_mass = r"C:\Jonah\Uni\PhD\OESI\OpenFAST_cases\IEA-15_Full_LINEAR_Falling\main.1.lin"
M_s, __ = linData.rigidBodyMass(filepath_mass)
# Added mass at infinite frequency
filepath_Ainf = r"C:\Jonah\Uni\PhD\OESI\OpenFAST_cases\IEA-15_Full\HydroData\IEA-15-240-RWT-UMaineSemi.1"
A_inf = readWAMIT.Ainf(filepath_Ainf)
# Total mass
M = M_s + A_inf
active_dof = list(dof.values()) # Use this and next line to make sure M has correct dimensions
M = M[np.ix_(active_dof, active_dof)]
M_inv = np.linalg.inv(M)

# Stiffness (K) and damping (C)
filepath_KC = r'C:\Jonah\Uni\PhD\OESI\OpenFAST_cases\IEA-15_PlatformOnly_LINEAR\No_waves\main.1.lin' #TODO this linearization is not good, K and C values are off
inv_M_K = linData.invMK(filepath_KC, dof)   # M^(-1)*K
inv_M_C = linData.invMC(filepath_KC, dof)   # M^(-1)*C

# SS and SS_EXCTN matrices
A_r, B_r, C_r, n_r = read_ssFiles.read_ss(r'C:\Jonah\Uni\PhD\OESI\OpenFAST_cases\IEA-15_PlatformOnly_LINEAR\No_waves\HydroData\IEA-15-240-RWT-UMaineSemi.ss')
A_exc, B_exc, C_exc, n_exc = read_ssFiles.read_ssexctn(r'C:\Jonah\Uni\PhD\OESI\OpenFAST_cases\IEA-15_PlatformOnly_LINEAR\No_waves\HydroData\IEA-15-240-RWT-UMaineSemi.ssexctn')

# Load wave elevation time series from OpenFAST
# df = FASTOutputFile(r'C:\Jonah\Uni\PhD\OESI\OpenFAST_cases\IEA-15_PlatformOnly_LINEAR\main.out').toDataFrame() # No waves
df = FASTOutputFile(r'C:\Jonah\Uni\PhD\OESI\OpenFAST_cases\IEA-15_PlatformOnly_LINEAR\Waves\main.out').toDataFrame() # With waves
eta = df['Wave1Elev_[m]'].to_numpy()
time = df['Time_[s]'].to_numpy()

# # When only using 1DOF in OpenFAST (Heave)
# stiff_term = -0.0988420953056  # -M^-1 * K from OpenFAST
# damp_term  = -0.00406290937617  # -M^-1 * C from OpenFAST

# # 2. Extract the positive standalone terms (removing the negative sign)
# heave_inv_M_K = -stiff_term     # Reverts to +0.0988420953056
# heave_inv_M_C = -damp_term      # Reverts to +0.00406290937617
# # 3. Initialize the 6x6 standalone matrices with stable dummy values
# # (This keeps unneeded DOFs structurally locked and prevents division-by-zero errors)
# inv_M_K = np.eye(6) * 1e5       # High dummy stiffness for off-DOFs
# inv_M_C = np.eye(6) * 1e3       # High dummy damping for off-DOFs
# # 4. Inject the positive OpenFAST values into the Heave position 
# # (Assuming index 2 is your Heave/3rd DOF)
# inv_M_K[2, 2] = heave_inv_M_K
# inv_M_C[2, 2] = heave_inv_M_C

# Set up A and B matrix
doF = 6
# A = np.block([[np.zeros((doF,doF)), np.eye(doF), np.zeros((doF,n_r)), np.zeros((doF,n_exc))],
#               [-inv_M_K, -inv_M_C, -M_inv @ C_r, M_inv @ C_exc],
#               [np.zeros((n_r,doF)), B_r, A_r, np.zeros((n_r,n_exc))],
#               [np.zeros((n_exc,doF)), np.zeros((n_exc,doF)), np.zeros((n_exc,n_r)), A_exc]
#               ])

A = np.block([
    [np.zeros((n_dofs, n_dofs)), np.eye(n_dofs),             np.zeros((n_dofs, n_r)), np.zeros((n_dofs, n_exc))],
    [-inv_M_K,                   -inv_M_C,                   -M_inv @ C_r,            M_inv @ C_exc            ],
    [np.zeros((n_r, n_dofs)),    B_r,                        A_r,                     np.zeros((n_r, n_exc))   ],
    [np.zeros((n_exc, n_dofs)),  np.zeros((n_exc, n_dofs)),  np.zeros((n_exc, n_r)),  A_exc                    ]
])

B = np.block([[np.zeros((doF,1))],
              [np.zeros((doF,1))],
              [np.zeros((n_r,1))],
              [B_exc]
              ])

### ALTERNATIVE: Use A and B matrix from OpenFAST .lin file
lin = FASTLinearizationFile(r'C:\Jonah\Uni\PhD\OESI\OpenFAST_cases\IEA-15_PlatformOnly_LINEAR\No_waves\main.1.lin')
A_OF = lin['A']
A_OF = A_OF#[0:110,0:110]
# diff = A - A_OF
# A = A_OF
# B_all = lin['B']
# wave_input_idx = 6 
# B = B_all[:, [wave_input_idx]]#[0:110]

sysl = LinearStateSpace(A=A, B=B)
sysl.setInputTimeSeries(time,eta)
q0 = np.zeros(110) #np.zeros(100) when using OpenFAST matrices
q0[2] = -0.2613634003859304 #q0[0] when using OpenFAST matrices
# q0 = np.zeros(100)
# q0[0] = -0.2613634003859304
sysl.setStateInitialConditions(q0)
# print(sysl)
resln, _ = sysl.integrate(time, method='LSODA')
# resld, _ = sysl.integrate(time, method='impulse')

# =============================================================================
# --- Extract and Post-Process Results ---
# =============================================================================

# 1. Access the raw states array using SciPy's .y attribute
# Shape will be (n_states, n_timesteps) -> (116, n_timesteps)
x_states = resln.y  

# 2. Extract your 6 physical DOFs (the first 6 rows)
q_estimated = x_states[0:doF, :]  # Shape: (6, n_timesteps)

# Optional: Extract physical velocities (the next 6 rows)
q_dot_estimated = x_states[doF:2*doF, :]


# # --- Plots


plt.figure(figsize=(8, 2))
plt.plot(time,eta)
plt.grid(True, linestyle=':', alpha=0.6)
plt.xlabel('Time [s]')
plt.ylabel('Wave Elevation [m]')

# Remove mean in OpenFAST data for the plot
df['PtfmHeave_[m]'] = df['PtfmHeave_[m]'].values - np.mean(df['PtfmHeave_[m]'].values)

# Plot Heave only
plt.figure(figsize=(10, 5))
plt.plot(time, x_states[2,:],label=f'Estimated', color='royalblue', linewidth=1.5) #x_states[0,:] when using OpenFAST matrices
# plt.plot(time, x_states[0,:],label=f'Estimated', color='royalblue', linewidth=1.5)
plt.plot(time, df['PtfmHeave_[m]'].values,label='OpenFAST', color='crimson', linestyle='--', linewidth=1.2)
plt.legend(loc='upper right')
plt.xlabel('Time [s]')
plt.ylabel('Heave [m]')
plt.grid(True, linestyle=':', alpha=0.6)

# # Define labels and units for each degree of freedom
# dof_labels = ['Surge', 'Sway', 'Heave', 'Roll', 'Pitch', 'Yaw']
# dof_units  = ['[m]', '[m]', '[m]', '[rad]', '[rad]', '[rad]']

# # Create a figure with 6 stacked subplots sharing the same X-axis (Time)
# fig, axes = plt.subplots(6, 1, figsize=(10, 14), sharex=True)

# for i in range(doF):
#     # Plot the integrated state-space displacement trajectory
#     axes[i].plot(time, q_estimated[i, :], label=f'Estimated', color='royalblue', linewidth=1.5)
    
#     # -------------------------------------------------------------------------
#     # Optional Validation: Plot the OpenFAST ground truth if it exists in your df
#     # Map your 6 DOFs to the standard OpenFAST platform output channels
#     openfast_channels = ['PtfmSurge_[m]', 'PtfmSway_[m]', 'PtfmHeave_[m]', 'PtfmRoll_[deg]', 'PtfmPitch_[deg]', 'PtfmYaw_[deg]']
#     of_ch = openfast_channels[i]
    
#     if of_ch in df.columns:
#         # Note: OpenFAST outputs angular DOFs in DEGREES. 
#         # If your state-space matrices use RADIANS, convert OpenFAST here:
#         if i >= 3:  # Roll, Pitch, Yaw indices
#             of_data = df[of_ch].values * (np.pi / 180.0) # Convert deg to rad
#             lbl_unit = '[rad]'
#         else:
#             of_data = df[of_ch].values
#             lbl_unit = '[m]'
            
#         axes[i].plot(time, of_data, label='OpenFAST', color='crimson', linestyle='--', linewidth=1.2)
#         axes[i].set_ylabel(f'{dof_labels[i]} {lbl_unit}')
#     else:
#         # Fallback if channel names don't match your file exactly
#         axes[i].set_ylabel(f'{dof_labels[i]} {dof_units[i]}')
#     # -------------------------------------------------------------------------

#     axes[i].grid(True, linestyle=':', alpha=0.6)
#     axes[i].legend(loc='upper right')

# # Label the shared x-axis on the bottom plot
# axes[-1].set_xlabel('Time [s]')

# # Format layout
# fig.suptitle('Platform Structural Trajectory (q) Validation', fontsize=14, fontweight='bold', y=0.98)
# plt.tight_layout()
# plt.show()