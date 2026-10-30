# Enable Autoreload
import notebook_setup # type: ignore
notebook_setup.enable_autoreload()

# Load packages
import numpy as np
import matplotlib.pyplot as plt
from seastateEst import linData
from seastateEst import read_ssFiles
from seastateEst import wave_filter

# Set LaTeX font for figures
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'


### --- Testing the code

### Load rigid body mass
# filepath = r"C:\Jonah\Uni\PhD\OESI\OpenFAST_cases\IEA-15_PlatformOnly_LINEAR_Falling\main.1.lin"
filepath = r"C:\Jonah\Uni\PhD\OESI\OpenFAST_cases\IEA-15_Full_LINEAR_Falling\main.1.lin"

M, M_inv = linData.rigidBodyMass(filepath)

### Read SSexctn and SS files
A_r, B_r, C_r, n_r = read_ssFiles.read_ss(r'C:\Jonah\Uni\PhD\OESI\Starting_documents\SS_Excitation_Fitting\alpha-Release\HydroData\OC4 Semi\marin_semi.ss')
A_exc, B_exc, C_exc, n_exc = read_ssFiles.read_ssexctn(r'C:\Jonah\Uni\PhD\OESI\Starting_documents\SS_Excitation_Fitting\alpha-Release\marin_semi.ssexctn')

### Testing wave filter
T_p = 10.0
H_s = 2.5
omega_0 = 2 * np.pi / T_p
f_0 = 1/T_p

zeta = 0.12

f = np.linspace(1e-5,0.5,1000)

S_eta = wave_filter.waveFilter2Ord(zeta, H_s, T_p, f)
S_PM = wave_filter.PM_waves(H_s,T_p,f)
S_JONSWAP = wave_filter.JONSWAP_waves(H_s,T_p,f)

plt.figure()
plt.plot(f, S_PM, label = 'Pierson Moskowitz')
plt.plot(f, S_JONSWAP, label = 'JONSWAP')
plt.plot(f, S_eta, label = '2nd order filter')
plt.legend()
plt.xlabel('$f$ in Hz')
plt.ylabel('PSD in m$^2$s')