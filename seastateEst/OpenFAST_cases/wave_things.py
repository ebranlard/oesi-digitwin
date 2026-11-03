# Enable Autoreload
import notebook_setup # type: ignore
notebook_setup.enable_autoreload()

# Load packages
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from openfast_toolbox.io import FASTOutputFile # type: ignore
from openfast_toolbox.io import FASTLinearizationFile # type: ignore

# Set LaTeX font for figures
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'

# def rigidBodyMass(filepath):
#     """
#     This function finds the rigid body 6x6 mass matrix and it's inverse from an OpenFAST .lin file.
#     This mass contains the infinite-frequency hydrodynamic added mass A_inf if Hydrodyne was activated.
#     """

#     lin = FASTLinearizationFile(filepath)

#     # Extract the state-space matrices as numpy arrays
#     A = lin['A']
#     B = lin['B']
#     C = lin['C']
#     D = lin['D']

#     # Extract the lists of descriptions
#     input_descriptions = lin['u_info']['Description']
#     state_descriptions = lin['x_info']['Description']

#     # Find input columns starting with "ED Platform"
#     # We match the first 6 occurrences corresponding to Surge, Sway, Heave, Roll, Pitch, Yaw forces
#     input_indices = [
#         i for i, description in enumerate(input_descriptions) 
#         if description.strip().startswith("ED Platform")
#     ][:6]

#     # Find state rows starting with "ED First time derivative" 
#     # AND containing "Platform" or "Ptfm" to isolate the velocity states
#     state_deriv_indices = [
#         i for i, description in enumerate(state_descriptions) 
#         if description.strip().startswith("ED First time derivative") and ("Platform" in description or "Ptfm" in description)
#     ]

#     # Extract the 6x6 B matrix subset
#     B_subset = B[np.ix_(state_deriv_indices, input_indices)]

#     # Find the rigid body 6x6 mass matrix by inverting B_subset
#     M = np.linalg.inv(B_subset)
#     # print(input_indices)
#     # print(state_deriv_indices)
#     return M, B_subset

########## Testing the code

# Load the .lin file

# filepath = r"C:\Jonah\Uni\PhD\OESI\OpenFAST_cases\IEA-15_PlatformOnly_LINEAR_Falling\main.1.lin"
# filepath = r"C:\Jonah\Uni\PhD\OESI\OpenFAST_cases\IEA-15_Full_LINEAR_Falling\main.1.lin"

# M, M_inv = rigidBodyMass(filepath)

#### For NAWEA 2026 abstract

# Read the .out or .outb file
#df = FASTOutputFile('IEA-15_PlatformOnly/main.out').toDataFrame()
#df = FASTOutputFile('IEA-15_PlatformOnly_LINEAR/main.out').toDataFrame()
df = FASTOutputFile('IEA-15_Full/main.out').toDataFrame()
print(df.columns.values)

# Convert wave elevation to a NumPy array
Wave1Elev = df['Wave1Elev_[m]'].to_numpy()

plt.figure()
df.plot(x='Time_[s]', y='Wave1Elev_[m]', grid=True)
plt.ylabel('Wave1Elev_[m]')
plt.show()

plt.figure()
df.plot(x='Time_[s]', y='PtfmTAzi_[m/s^2]', grid=True)
plt.ylabel('Heave Accel')
plt.show()

# 1. Compute PSD
dt = df['Time_[s]'].diff().mean()  # Average time step
fs = 1 / dt                    # Sampling frequency (Hz)

freqs, psd_elev = signal.welch(df['Wave1Elev_[m]'], fs=fs, nperseg=4096) # nperseg defines the frequency resolution; higher means more detail but more noise: 2048, 4096, 8192
freqs, psd_accel = signal.welch(df['PtfmTAzi_[m/s^2]'], fs=fs, nperseg=4096)

# 2. Plot the results
plt.figure(figsize=(10, 5))
plt.plot(freqs, psd_elev)
plt.title('Wave Elevation Spectrum')
plt.xlabel('Frequency (Hz)')
plt.ylabel('PSD ($m^2/Hz$)')
plt.xlim(0, 0.5)  # Focus on physical wave frequencies (usually < 0.5 Hz)
plt.grid(True)
plt.show()

plt.figure(figsize=(10, 5))
plt.plot(freqs, psd_accel)
plt.title('Heave accel Spectrum')
plt.xlabel('Frequency (Hz)')
plt.ylabel('PSD ($m^2/Hz$)')
plt.xlim(0, 0.5)  # Focus on physical wave frequencies (usually < 0.5 Hz)
plt.grid(True)
plt.show()

# 1. Normalize the PSDs
psd_elev_norm = psd_elev / np.max(psd_elev)
psd_accel_norm = psd_accel / np.max(psd_accel)

# 2. Plot both on the same axis
plt.figure(figsize=(10, 6))

plt.plot(freqs, psd_elev_norm, label='Wave Elevation', linewidth=1.5)
plt.plot(freqs, psd_accel_norm, label='Heave Acceleration', linewidth=1.5, linestyle='--')

plt.xlabel('Frequency (Hz)', fontsize=18)
plt.ylabel('Normalized PSD (-)', fontsize=18) # No units since it is normalized
plt.xlim(0, 0.5) 
plt.ylim(0, 1.1) # Buffer at the top for clarity
plt.grid(True, which='both', linestyle=':', alpha=0.5)
plt.legend(fontsize=18)
plt.tick_params(axis='both', labelsize=14)
# plt.savefig(r"C:\Jonah\Uni\PhD\Conferences\NAWEA_2026\OESI\figures\spectrum.pdf", dpi=300, bbox_inches='tight')
plt.show()