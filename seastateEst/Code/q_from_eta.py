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
from seastateEst import FRad

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

# --- Structure mass
filepath_mass = r"C:\Jonah\Uni\PhD\OESI\OpenFAST_cases\IEA-15_Full_LINEAR_Falling\main.1.lin"
M_s, __ = linData.rigidBodyMass(filepath_mass)
# Added mass at infinite frequency
filepath_Ainf = r"C:\Jonah\Uni\PhD\OESI\OpenFAST_cases\IEA-15_Full\HydroData\IEA-15-240-RWT-UMaineSemi.1"
A_inf = FRad.Ainf(filepath_Ainf, WtrDens=1025)
# Total mass
M = M_s + A_inf
# Filter for active dofs
active_dof = list(dof.values())
M = M[np.ix_(active_dof, active_dof)]
# Compute inverse mass
M_inv = np.linalg.inv(M)

# --- Stiffness (K) and damping (C)
filepath_KC = r'C:\Jonah\Uni\PhD\OESI\OpenFAST_cases\IEA-15_PlatformOnly_LINEAR\No_waves\main.1.lin' #TODO chose from what OpenFAST simulation I really want K and C
inv_M_K = linData.invMK(filepath_KC, dof)   # M^(-1)*K
inv_M_C = linData.invMC(filepath_KC, dof)   # M^(-1)*C

# --- SS and SS_EXCTN matrices
# NAWEA 2026 cases:
# A_r, B_r, C_r, n_r = read_ssFiles.read_ss(r'C:\Jonah\Uni\PhD\OESI\OpenFAST_cases\IEA-15_PlatformOnly_LINEAR\No_waves\HydroData\IEA-15-240-RWT-UMaineSemi.ss')
# A_exc, B_exc, C_exc, n_exc = read_ssFiles.read_ssexctn(r'C:\Jonah\Uni\PhD\OESI\OpenFAST_cases\IEA-15_PlatformOnly_LINEAR\No_waves\HydroData\IEA-15-240-RWT-UMaineSemi.ssexctn')

# Comparison simple lin model with OF
A_r, B_r, C_r, n_r = read_ssFiles.read_ss(r'C:\Jonah\Uni\PhD\OESI\OpenFAST_cases\IEA-15_PlatformOnly_LINEAR\Comparison_Simple_lin_model\HydroData\marin_semi.ss')
A_exc, B_exc, C_exc, n_exc = read_ssFiles.read_ssexctn(r'C:\Jonah\Uni\PhD\OESI\OpenFAST_cases\IEA-15_PlatformOnly_LINEAR\Comparison_Simple_lin_model\HydroData\marin_semi.ssexctn')


# Filter for active dofs ONLY if this was not done when creating the ss and ssexctn files
B_r = B_r[:, active_dof]
C_r =  C_r[active_dof, :]
C_exc =  C_exc[active_dof, :]

# --- Load wave elevation time series from OpenFAST
# df = FASTOutputFile(r'C:\Jonah\Uni\PhD\OESI\OpenFAST_cases\IEA-15_PlatformOnly_LINEAR\No_waves\main.out').toDataFrame() # Heave only no waves
# df = FASTOutputFile(r'C:\Jonah\Uni\PhD\OESI\OpenFAST_cases\IEA-15_PlatformOnly_LINEAR\Waves\main.out').toDataFrame() # Heave only with waves
# df = FASTOutputFile(r'C:\Jonah\Uni\PhD\OESI\OpenFAST_cases\IEA-15_PlatformOnly_LINEAR\No_waves_pitch\main.out').toDataFrame() # Pitch only no waves
# df = FASTOutputFile(r'C:\Jonah\Uni\PhD\OESI\OpenFAST_cases\IEA-15_PlatformOnly_LINEAR\Waves_pitch\main.out').toDataFrame() # Pitch only with waves
df = FASTOutputFile(r'C:\Jonah\Uni\PhD\OESI\OpenFAST_cases\IEA-15_PlatformOnly_LINEAR\Comparison_Simple_lin_model\main.out').toDataFrame() # Comparison simple lin model with OF
n_r = 0

eta = df['Wave1Elev_[m]'].to_numpy()
time = df['Time_[s]'].to_numpy()

# A = np.block([
#     [np.zeros((n_dofs, n_dofs)), np.eye(n_dofs),             np.zeros((n_dofs, n_exc)), np.zeros((n_dofs, n_r))],
#     [-inv_M_K,                   -inv_M_C,                   M_inv @ C_exc,             -M_inv @ C_r           ],
#     [np.zeros((n_exc, n_dofs)),  np.zeros((n_exc, n_dofs)),  A_exc,                     np.zeros((n_exc, n_r)) ],
#     [np.zeros((n_r, n_dofs)),    B_r,                        np.zeros((n_r, n_exc)),    A_r                    ]
# ])

B = np.block([[np.zeros((n_dofs,1))],
              [np.zeros((n_dofs,1))],
              [B_exc],
              [np.zeros((n_r,1))]
              ])


# --- ALTERNATIVE: Use A and B matrix from OpenFAST .lin file
lin = FASTLinearizationFile(r'C:\Jonah\Uni\PhD\OESI\OpenFAST_cases\IEA-15_PlatformOnly_LINEAR\Comparison_Simple_lin_model\main.1.lin') # Comparison simple lin model with OF
# lin = FASTLinearizationFile(r'C:\Jonah\Uni\PhD\OESI\OpenFAST_cases\IEA-15_PlatformOnly_LINEAR\No_waves\main.1.lin') # Heave only
# lin = FASTLinearizationFile(r'C:\Jonah\Uni\PhD\OESI\OpenFAST_cases\IEA-15_PlatformOnly_LINEAR\No_waves_pitch\main.1.lin') # Pitch only
A_OF = lin['A']
A_OF = A_OF
# diff = A - A_OF
A = A_OF
# B_all = lin['B']
# wave_input_idx = 6 
# B = B_all[:, [wave_input_idx]]

# q0 = np.zeros(100)
q0 = np.zeros(8)

if dof["Heave"]:
    # q0[0] = -0.2613634003859304 # Heave NAWEA 2026 case
    q0[0] = -0.009784794055073265 # Comparison simple lin model with OF


if dof["Pitch"]:
    A[1,1] = A[1,1] * 0.3 # THIS IS A FUDGE FACTOR for pitch damping: TODO REMOVE!
    # q0[0] = 2.5 * (np.pi / 180.0) # Pitch

sysl = LinearStateSpace(A=A, B=B)
sysl.setInputTimeSeries(time,eta)

sysl.setStateInitialConditions(q0)
# # print(sysl)
resln, _ = sysl.integrate(time, method='LSODA')

# =============================================================================
# --- Extract and Post-Process Results ---
# =============================================================================

# 1. Access the raw states array using SciPy's .y attribute
# Shape will be (n_states, n_timesteps) -> (116, n_timesteps)
x_states = resln.y  

# 2. Extract activated dofs
q_estimated = x_states[0:n_dofs, :]

# # Optional: Extract physical velocities
# q_dot_estimated = x_states[n_dofs:2*n_dofs, :]


# --- Plots
plt.figure(figsize=(10, 3))
plt.plot(time,eta)
plt.grid(True, linestyle=':', alpha=0.6)
plt.xlabel('Time [s]')
plt.ylabel('Wave Elevation [m]')

if dof["Heave"]:
    # Remove mean in OpenFAST data for the plot
    df['PtfmHeave_[m]'] = df['PtfmHeave_[m]'].values - np.mean(df['PtfmHeave_[m]'].values)

    # Plot Heave only
    plt.figure(figsize=(10, 3))
    plt.plot(time, q_estimated.flatten(),label=f'Estimated', color='royalblue', linewidth=1.5) #TODO: If multiple DOFs are enabled this will plot the wrong thing
    plt.plot(time, df['PtfmHeave_[m]'].values,label='OpenFAST', color='crimson', linestyle='--', linewidth=1.2)
    plt.legend(loc='upper right')
    plt.xlabel('Time [s]')
    plt.ylabel('Heave [m]')
    plt.grid(True, linestyle=':', alpha=0.6)

if dof["Pitch"]:
    # Plot Pitch only
    plt.figure(figsize=(10, 3))
    plt.plot(time, q_estimated.flatten(),label=f'Estimated', color='royalblue', linewidth=1.5) #TODO: If multiple DOFs are enabled this will plot the wrong thing
    plt.plot(time, df['PtfmPitch_[deg]'].values * (np.pi / 180.0),label='OpenFAST', color='crimson', linestyle='--', linewidth=1.2)
    plt.legend(loc='upper right')
    plt.xlabel('Time [s]')
    plt.ylabel('Pitch [rad]')
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