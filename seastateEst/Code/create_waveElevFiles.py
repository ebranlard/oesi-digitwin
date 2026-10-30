import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.signal import chirp

### --- Parameters
dt = 0.1
tmax = 800
p1 = 5.5
p2 = 10
a1 = 0.1
a2 = 0.29

### --- Create time series
time = np.arange(0,tmax,dt)
omega1 = 2*np.pi / p1
omega2 = 2*np.pi / p2

## - Regular sin wave
Wave1Elev = a2 * np.sin(omega2 * time)

## - Two frequencies
# Wave1Elev = a1 * np.sin(omega1 * time) + a2 * np.cos(omega2 * time)

## - 4 highly distinct, non-matching periods
# w1 = 2 * np.pi / 4.1
# w2 = 2 * np.pi / 7.3
# w3 = 2 * np.pi / 11.9
# w4 = 2 * np.pi / 16.5
# Wave1Elev = (1.5 * np.sin(w1 * time) + 
#              1.0 * np.cos(w2 * time) + 
#              2.0 * np.sin(w3 * time) + 
#              0.8 * np.cos(w4 * time))

## - Stepfunction
# Wave1Elev = np.where(time < 50, 0.0, 2.0)

## - Chirp
# Wave1Elev = a1 * chirp(time, f0=0.02, t1=tmax, f1=0.5, method='linear')

## - Multi-Frequency Wave Packets
# center_time = 100
# sigma = 50.0  # Controls the width of the packet
# envelope = np.exp(-((time - center_time) ** 2) / (2 * sigma ** 2))
# Wave1Elev = a2 * np.sin(omega2 * time) * envelope

## - Sudden frequency change
# current_omega = np.where(time < 100, omega1, omega2)
# Wave1Elev = a1 * np.sin(current_omega * time)

### --- Ensure always starts at 0 and has zero mean
# 1. Find where consecutive elements have different signs
# np.signbit returns True for negative values, False for positive/zero
sign_changes = np.diff(np.signbit(Wave1Elev))

# 2. Get the index of the very first sign change
first_crossing_idx = np.where(sign_changes)[0][0] + 1

# 3. Slice the arrays from the first crossing onward
# Wave1Elev = Wave1Elev[first_crossing_idx:]
# time = time[first_crossing_idx:]

### --- Write to .Elev file
filename = 'MyWaves_Start0'
output_folder = r'C:\Jonah\Uni\PhD\OESI\OpenFAST_cases'

with open(os.path.join(output_folder, filename + '.Elev'), "w") as f:
    f.write("      Time\t Wave1Elev\n     (sec)\t       (m)\n")
    for t, elev in zip(time, Wave1Elev):
        f.write(f"    {t:.4f}\t {elev:.3e}\n")

### --- Plot
plt.figure()
plt.plot(time,Wave1Elev)
plt.xlabel('Time (s)')
plt.ylabel('Wave Elevation (m)')
plt.xlim([0, 200])