"""
This library includes the functions associated with computing the excitation Force in Cummins Equaion.
"""
import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
from openfast_toolbox.io import FASTOutputFile # type: ignore

def complex_WaveElev(eta_t, t):
    """
    This function computes the complex wave spectrum from a time series of wave elevation.

    INPUTS:
          eta_t: Wave elevation time series [m]
          t: Time [s]
    """
    # Time step
    dt = t[1] - t[0]

    # Compute the complex wave spectrum
    eta_omega = np.fft.fft(eta_t) * dt
    omega_fft = np.fft.fftfreq(len(eta_t), d=dt) * 2 * np.pi  # rad/s

    return eta_omega, omega_fft

def ExctnKrnl(filepath, WtrDens=1025, g=9.80665, beta=0, WAMITULEN=1):
    '''
    Loads a WAMIT .3 file and returns sorted frequencies (rad/s) and 
    the dimensional complex excitation FRF K_exc(omega) for all 6 DOFs.

    INPUTS:
          filepath: Filepath to WAMIT .3 file
          WtrDens: Water density (kg/m^3)
          g: Gravitational acceleration (m/s^2)
          beta: Wave heading angle (deg)
          WAMITULEN: Characteristic body length scale used to redimensionalize WAMIT output (meters)
    '''

    # Open the WAMIT .3 file
    # Columns: Period (s), BETA (Wave Heading) (deg), I (DOF 1-6), Mod(Xi) (Exciting force magnitude), PhaXi (Phase of force in deg), ReXi (Real part of force), ImXi (Imag part of force)
    wamit_data = pd.read_csv(filepath, delim_whitespace=True, names=["per", "beta", "i", "ModXi", "PhaXi", "ReXi", "ImXi"])


    # Convert Period (per > 0) to angular frequency omega [rad/s]
    wamit_data["omega"] = np.where(
        wamit_data["per"] > 0,
        2 * np.pi / wamit_data["per"],
        np.abs(wamit_data["per"])
    )

    # Filter for the chosen wave heading (beta)
    df = wamit_data[np.isclose(wamit_data["beta"], beta)].copy()

    # Build nondimensional complex FRF (=K_exc(omega))
    df["FRF"] = df["ReXi"] + 1j * df["ImXi"]

    # Pivot so rows = sorted frequencies (rad/s) and columns = DOFs 1..6
    # pivoted = df.pivot(index="per", columns="i", values="FRF").sort_index()
    pivoted = df.pivot(index="omega", columns="i", values="FRF").sort_index()

    omega_Krnl = pivoted.index.to_numpy()
    K_exc = pivoted.to_numpy()

    # Redimensionalize: rho * g * L^2 for forces (DOFs 1-3), rho * g * L^3 for moments (DOFs 4-6) (SEE WAMIT Manual Sec 3.3)
    K_exc[:, :3] *= WtrDens * g * (WAMITULEN ** 2)
    K_exc[:, 3:] *= WtrDens * g * (WAMITULEN ** 3)
    return K_exc, omega_Krnl


def IFFT_Exctn(K_exc, omega_wamit, eta_omega, omega_fft, dt):
    """
    Interpolates K_exc onto the two-sided FFT frequency grid, multiplies by
    eta(omega), and computes the time-domain wave excitation force via IFFT.

    INPUTS:
          K_exc: Complex FRF array from WAMIT, shape (N_wamit, 6)
          omega_wamit: Sorted WAMIT frequencies [rad/s], shape (N_wamit,)
          eta_omega: Two-sided complex wave spectrum, shape (N_time,)
          omega_fft: Two-sided FFT frequencies [rad/s], shape (N_time,)
          dt: Time step [s]

    RETURNS:
          F_exc_t: Time-domain wave excitation forces/moments, shape (N_time, 6)
    """
    N_time = len(omega_fft)
    K_interp = np.zeros((N_time, 6), dtype=complex)

    # Identify positive (including 0 Hz) and negative frequency indices
    pos_mask = omega_fft >= 0
    neg_mask = omega_fft < 0

    for dof in range(6):
        # 1. Interpolate Real and Imaginary parts onto positive FFT frequencies
        # left=0.0, right=0.0 sets forces to zero outside WAMIT's frequency range
        re_pos = np.interp(omega_fft[pos_mask], omega_wamit, K_exc[:, dof].real, left=0.0, right=0.0)
        im_pos = np.interp(omega_fft[pos_mask], omega_wamit, K_exc[:, dof].imag, left=0.0, right=0.0)
        K_interp[pos_mask, dof] = re_pos + 1j * im_pos

        # 2. Enforce Hermitian symmetry for negative frequencies: K(-w) = conj(K(+w))
        re_neg = np.interp(-omega_fft[neg_mask], omega_wamit, K_exc[:, dof].real, left=0.0, right=0.0)
        im_neg = np.interp(-omega_fft[neg_mask], omega_wamit, K_exc[:, dof].imag, left=0.0, right=0.0)
        K_interp[neg_mask, dof] = re_neg - 1j * im_neg  # Note the minus sign for conjugate

    # 3. Multiply in the frequency domain: F(w) = K(w) * eta(w) for all 6 DOFs
    F_exc_omega = K_interp * eta_omega[:, np.newaxis]

    # 4. Inverse FFT back to time domain (divide by dt to undo the dt scaling in complex_WaveElev)
    F_exc_t = np.fft.ifft(F_exc_omega / dt, axis=0).real

    return F_exc_t

### --- Testing the code

# filepath_WAMIT = r'C:\Jonah\Uni\PhD\OESI\seastateEst\OpenFAST_cases\IEA-15_PlatformOnly_v5.0\HydroData\marin_semi.3'
# K_exc, omega_Kernel = ExctnKrnl(filepath_WAMIT)

# # --- Load measurements from OpenFAST
# filepath_OF = r'C:\Jonah\Uni\PhD\OESI\seastateEst\OpenFAST_cases\IEA-15_PlatformOnly_v5.0'
# dfOut = FASTOutputFile(os.path.join(filepath_OF, 'main.out')).toDataFrame()

# # Load your time series from the OpenFAST .out file
# eta_t = dfOut['Wave1Elev_[m]'].to_numpy()  # Wave elevation time series
# time = dfOut['Time_[s]'].to_numpy()

# eta_omega, omega_fft = complex_WaveElev(eta_t, time)

# # Compute dt and run the frequency-domain multiplication + IFFT
# dt = time[1] - time[0]

# F_exc_t = IFFT_Exctn(K_exc, omega_Kernel, eta_omega, omega_fft, dt)

# F_heave_OF = dfOut['B1WvsF1zi_[N]'].to_numpy()
# diff_heave = F_exc_t[:, 2] - F_heave_OF
# plt.figure()
# plt.plot(time, diff_heave)
# plt.ylim([-0.3e6, 0.3e6])
# ### --- Plots

# # Define channel names, labels, and DOF indices for all 6 DOFs
# dof_configs = [
#     ("B1WvsF1xi_[N]",   0, "Surge Force [N]"),
#     ("B1WvsM1xi_[N-m]", 3, "Roll Moment [N·m]"),
#     ("B1WvsF1yi_[N]",   1, "Sway Force [N]"),
#     ("B1WvsM1yi_[N-m]", 4, "Pitch Moment [N·m]"),
#     ("B1WvsF1zi_[N]",   2, "Heave Force [N]"),
#     ("B1WvsM1zi_[N-m]", 5, "Yaw Moment [N·m]"),
# ]

# fig, axes = plt.subplots(3, 2, figsize=(16, 9), sharex=True)
# axes = axes.flatten()

# for ax, (ch_name, dof_idx, ylabel) in zip(axes, dof_configs):
#     if ch_name in dfOut.columns:
#         ax.plot(time, dfOut[ch_name], label=f"OpenFAST ({ch_name.split('_')[0]})", linewidth=1.5)
#     else:
#         print(f"Warning: {ch_name} not found in dfOut. Available wave channels: {[c for c in dfOut.columns if 'Wvs' in c]}")
        
#     ax.plot(time, F_exc_t[:, dof_idx], "--", label="Python F_exc", linewidth=1.2)
#     ax.set_ylabel(ylabel)
#     ax.legend(loc="upper right")
#     ax.grid(True)

# axes[4].set_xlabel("Time [s]")
# axes[5].set_xlabel("Time [s]")

# plt.tight_layout()
# plt.show()

# plt.figure()
# plt.plot(frequencies,eta_omega)
# plt.xlim([0, 10])