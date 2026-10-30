# Enable Autoreload
import notebook_setup # type: ignore
notebook_setup.enable_autoreload()

# General packages
import numpy as np
import matplotlib.pyplot as plt

# OpenFAST toolbox
from openfast_toolbox.io import FASTOutputFile # type: ignore

# Seastate Estimator
from seastateEst import FRad

### --- Load OpenFAST .out and HD.sum file
# v5.0 files
dfOut = FASTOutputFile(r'C:\Jonah\Uni\PhD\OESI\seastateEst\OpenFAST_cases\IEA-15_PlatformOnly_v5.0\main.out').toDataFrame()
filepath_HDsum = r'C:\Jonah\Uni\PhD\OESI\seastateEst\OpenFAST_cases\IEA-15_PlatformOnly_v5.0\main.HD.sum'
# v4.1.2 files
# dfOut = FASTOutputFile(r'C:\Jonah\Uni\PhD\OESI\OpenFAST_cases\IEA-15_PlatformOnly\main.out').toDataFrame()
# filepath_HDsum = r'C:\Jonah\Uni\PhD\OESI\OpenFAST_cases\IEA-15_PlatformOnly\main.HD.sum'
dfHDsum = FRad.RdtnKrnl(filepath_HDsum)

# Velocity and acceleration channels in dfOut
velocity_channels = ['PtfmTVxi_[m/s]', 'PtfmTVyi_[m/s]', 'PtfmTVzi_[m/s]', 'PtfmRVxi_[deg/s]', 'PtfmRVyi_[deg/s]', 'PtfmRVzi_[deg/s]']
accel_channels = ['PtfmTAxi_[m/s^2]', 'PtfmTAyi_[m/s^2]', 'PtfmTAzi_[m/s^2]', 'PtfmRAxi_[deg/s^2]', 'PtfmRAyi_[deg/s^2]', 'PtfmRAzi_[deg/s^2]']

# Transform into a (6xN) numpy array
velocity_history = dfOut[velocity_channels].values.T
velocity_history[3:6, :] = np.radians(velocity_history[3:6, :]) # Convert Roll, Pitch and Yaw velocities from deg/s to rad/s
accel_history = dfOut[accel_channels].values.T
accel_history[3:6, :] = np.radians(accel_history[3:6, :])

time = dfOut['Time_[s]'].to_numpy()
# Fmem_surge = dfOut['B1RdtFxi_[N]'].to_numpy()
# Fmem_sway = dfOut['B1RdtFyi_[N]'].to_numpy()
# Fmem_heave = dfOut['B1RdtFzi_[N]'].to_numpy()
# Fmem_roll = dfOut['B1RdtMxi_[N-m]'].to_numpy()
# Fmem_pitch = dfOut['B1RdtMyi_[N-m]'].to_numpy()
# Fmem_yaw = dfOut['B1RdtMzi_[N-m]'].to_numpy()

# Added mass at infinite frequency
filepath_Ainf = r"C:\Jonah\Uni\PhD\OESI\seastateEst\OpenFAST_cases\IEA-15_Full\HydroData\IEA-15-240-RWT-UMaineSemi.1"
A_inf = FRad.Ainf(filepath_Ainf, WtrDens=1025)

# Load the nasty radiation_debug text file
data = np.loadtxt(r'C:\Jonah\Uni\PhD\OESI\seastateEst\OpenFAST_cases\IEA-15_PlatformOnly_v5.0\radiation_debug.txt', skiprows=1)
time_txt = data[:, 0] # Slice out the Time column (Column index 0)
F_rdtn_txt = data[:, 1:7] # Slice out the 6 elements of F_Rdtn (Column indices 1 through 6)
F_rdtn_txt = F_rdtn_txt.T # Match the shape

### --- Calculate the convolution
# Fmem = FRad.Conv_Rdtn(dfHDsum, velocity_history, tmax=60)
# print('Fmem at current time tn:\n', Fmem)


#### --- Extract all time steps and setup storage
N_steps = len(time)
dt = time[1] - time[0]  # Assuming a constant time step
Fmem_history = np.zeros((6, N_steps))
Frad_history = np.zeros((6, N_steps))
F_added_mass = np.zeros((6, N_steps))

print("Running convolution validation loop...")

# --- Loop through time steps to simulate the growing history matrix
for n in range(N_steps):
    # 1. Slice the velocity history up to the current time step (1 to n+1 steps)
    # This matrix shape starts at (6, 1) and grows to (6, n+1)
    current_velocity_slice = velocity_history[:, :n+1]
    
    # 2. Compute Fmem using only the historical data available up to step n
    Fmem_at_tn = FRad.Conv_Rdtn(dfHDsum, current_velocity_slice, RdtnTMax=60)
    
    # 3. Store the result (assuming Conv_Rdtn returns a 6-element vector or array)
    Fmem_history[:, n] = Fmem_at_tn.squeeze()

    # 4. Pull the instantaneous acceleration vector for this specific time step (Shape: (6,))
    current_acc = accel_history[:, n]
    
    # 5. Compute the inertial added mass term: -Ainf @ acceleration
    F_added_mass[:, n] = -A_inf @ current_acc

# 6. Total Radiation Force = Memory Force + Added Mass Force
Frad_history = Fmem_history + F_added_mass

print("Calculation complete. Plotting validation results...")

# %%
# --- Plots
# 1. Create a figure with 3 stacked subplots sharing the same X-axis
fig, axs = plt.subplots(6, 1, figsize=(10, 18), sharex=True)

# Data mapping configuration for looping sequentially from DOF 0 to 5
openfast_loads = [
    dfOut['B1RdtFxi_[N]'].to_numpy(),    # 0: Surge
    dfOut['B1RdtFyi_[N]'].to_numpy(),    # 1: Sway
    dfOut['B1RdtFzi_[N]'].to_numpy(),    # 2: Heave
    dfOut['B1RdtMxi_[N-m]'].to_numpy(),  # 3: Roll
    dfOut['B1RdtMyi_[N-m]'].to_numpy(),  # 4: Pitch
    dfOut['B1RdtMzi_[N-m]'].to_numpy()   # 5: Yaw
]

dof_names = ['Surge', 'Sway', 'Heave', 'Roll', 'Pitch', 'Yaw']
openfast_labels = ['B1RdtFxi', 'B1RdtFyi', 'B1RdtFzi', 'B1RdtMxi', 'B1RdtMyi', 'B1RdtMzi']

# 2. Loop through and plot each of the 6 Degrees of Freedom (DOF)
for i in range(6):
    # Determine the appropriate unit text based on the DOF type
    unit = '[N]' if i < 3 else '[N-m]'
    
    # axs[i].plot(time, openfast_loads[i], label=f'OpenFAST Frad ({openfast_labels[i]})', color='black', linewidth=2)
    # axs[i].plot(time, Frad_history[i, :], label=f'My Script (Frad {dof_names[i]})', color='crimson', linestyle='--')
    axs[i].plot(time, F_rdtn_txt[i,:], label=f'OpenFAST (Fmem)', color='black', linewidth=2)
    axs[i].plot(time, Fmem_history[i, :], label=f'My Script (Fmem {dof_names[i]})', color='crimson', linestyle='--')
    
    # Subplot styling
    axs[i].set_ylabel(f'Radiation Load {unit}')
    axs[i].set_title(f'{dof_names[i]} Radiation Verification')
    axs[i].legend(loc='upper right')
    axs[i].grid(True)
    # axs[i].set_xlim([0, 10])

# 3. Overall layout adjustments
axs[5].set_xlabel('Time [s]')  # X-label only on the absolute bottom plot
fig.suptitle('Platform 6-DOF Radiation Load Verification (IEA-15MW Platform Only)', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.show()


# %%
