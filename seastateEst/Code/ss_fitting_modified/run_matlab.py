import os
import numpy as np
import shutil
from scipy.io import loadmat
import matlab.engine
import matplotlib.pyplot as plt

def startMatlabEngine():
    ### To Share MATLAB engine open MATLAB and type in the MATLAB command window: matlab.engine.shareEngine

    # Check for already open and shared MATLAB sessions
    print("Searching for an open MATLAB session...")
    shared_engines = matlab.engine.find_matlab()

    if shared_engines:
        print(f"Found open MATLAB session: '{shared_engines[0]}'. Connecting instantly...")
        # Connect to the first open MATLAB instance found
        eng = matlab.engine.connect_matlab(shared_engines[0])
    else:
        print("No open MATLAB found. Starting a new background engine... (takes a few seconds)")
        # Fallback to starting a new background one if none are open
        eng = matlab.engine.start_matlab()

    # Tell MATLAB to look inside your current directory for scripts
    current_directory = os.getcwd()
    eng.addpath(current_directory)

    print("MATLAB engine started successfully. Running script...")
    return eng

def stopMatlabEngine(eng):
    print("Disconnecting from MATLAB...")
    eng.quit()

def RedOrder_SSexc(n, options):
    '''
    This function filters the order vector based on active degrees of freedom and returns a matlab double.
    '''
    # Remove entries with index where DOF is set to 0
    n = [order for dof, order in zip(list(options["DoF"][0]), n) if dof == 1]

    # Convert to a matlab double
    n = matlab.double(n)
    return n

def recreate_matlab_plot(iteration_index=1):
    # Load the generated mat file
    mat_file = os.path.join("plot_data", f"plot_data_{iteration_index}.mat")
    try:
        data = loadmat(mat_file)
    except FileNotFoundError:
        print(
            f"Error: {mat_file} not found. Ensure the MATLAB function ran successfully."
        )
        return

    # Extract arrays and flatten them to 1D vectors
    t_orig = data["time_original"].flatten()
    y_orig = data["irf_original"].flatten()

    t_high = data["time_high_ss"].flatten()
    y_high = data["irf_high_ss"].flatten()

    t_low = data["time_low_ss"].flatten()
    y_low = data["irf_low_ss"].flatten()

    # Extract string metadata
    order = str(data["ss_order"][0])
    ylabel = str(data["ylabel_str"][0])
    # title_text = str(data["title_str"][0])
    title_text = str(data['title_str']).strip()

    # Recreate the precise plot layout
    plt.figure(figsize=(12, 5.5))

    # Match lines styles and colors from the original MATLAB image
    plt.plot(t_orig, y_orig, label="Original IRF", color="#1f77b4", linewidth=2)
    plt.plot(
        t_high,
        y_high,
        label="high order ss approx",
        color="#d62728",
        linewidth=2,
    )
    plt.plot(
        t_low,
        y_low,
        label=f"{order} th order ss approx",
        color="#ff7f0e",
        linewidth=2,
    )

    # Add styling details matching the figure
    plt.title(f"${title_text}$", fontsize=14)
    plt.xlabel("Time [s]")
    plt.ylabel(ylabel)
    plt.xlim(0, 40)  # Bound to 40 seconds
    plt.grid(True, linestyle=":", alpha=0.6)
    plt.legend(loc="upper right")

    # Clean up background style
    plt.tight_layout()
    plt.show()

### --- ### Run the script. 
## Inputs for SSRad
options_rad = {
    "FileName": r"C:\Jonah\Uni\PhD\OESI\seastateEst\Code\ss_fitting_modified\HydroData\marin_semi", # Wamit Files name and Location
    'DoF': matlab.double([1, 1, 1, 1, 1, 1]),   # Degree of Freedom enabled following the order Surge, Sway, Heave, Roll, Pitch and Yaw
    'twr': matlab.double([0.0, 2.5]),           # Typical local frequencies range (rad/s)
    'wwf': 1.0,                                 # Frequencies weighting factors
    'tfi': 2,                                   # Identification Method: 1-Freq Ident; 2-FDI; 3-Time-Domain LS method; 4-Time-Domain Realization theory
    'fit': 0.95,                                # Fit required for parametric model
    'ppmf': 1,                                  # Plot the parametric model fit (0 or 1)
    'fmt': 1                                    # Manual (1) or Automatic (0) order reduction (only relevant for method 4)
}

# eng = startMatlabEngine()
# eng.ss_fitting(options_rad,nargout=0) # nargout=0 is required for scripts or functions that do not return output variables.
# stopMatlabEngine(eng)


## Inputs for SSexc
tc = 12

options_exc = {
    "FileName": r"C:\Jonah\Uni\PhD\OESI\seastateEst\Code\ss_fitting_modified\HydroData\marin_semi", # Wamit Files name and Location
    "DoF": matlab.double([1, 1, 1, 1, 1, 1]),   # Degree of Freedom enabled
    "WaveTMax": 4600,                           # Analysis time for incident wave calculations (sec) [unused when WaveMod=0] [determines WaveDOmega=2Pi/WaveTMax in the IFFT]
    "dT": 0.2,                                  # Time step (should match FAST)
    "Fit": 0.99,                                # Fit required for the parametric model (max. recomended 0.99)
    "ppmf": 1,                                  # Plot the parametric model fit (0 or 1)
    "fmt": 1,                                   # Manual (1) or Automatic(0) order reduction (Only if method=4)
    "maxiter": 15,                              # maximum number of iterations for the Automatic method
    "WtrDen": 1025.0,                           # Water density (kg/m^3)
    "Grav": 9.80665,                            # Gravitational acceleration (m/s^2)
    "WaveDir": 0.0,                             # Incident wave propagation heading direction (degrees) [unused when WaveMod=0]
    "OutFileName": "marin_semi"                 # Output file name
}
base_n = [20, 20, 20, 20, 20, 20] # Order of the fit for each DoF
manual_order_SSexc = RedOrder_SSexc(base_n, options_exc) # Order of the fit filtered as needed for input of the MATLAB function

eng = startMatlabEngine()
eng.SS_Fitting_Wave_Excitation(options_exc,tc,manual_order_SSexc,nargout=0) # nargout=0 is required for scripts or functions that do not return output variables.
stopMatlabEngine(eng)

# Plot the results from the MATLAB function using the saved .mat file
recreate_matlab_plot(iteration_index=3)

### --- Running as a loop
# source_file_dir = r"C:\Jonah\Uni\PhD\OESI\Code\ss_fitting_modified\HydroData"
# dest_dir = r"C:\Jonah\Uni\PhD\OESI\OpenFAST_cases\Parametric_study\Lock_Rad\HydroData"

# ### --- While loop where the fit is varied and DoFs have different # states
# options_exc["Fit"] = 0.2    # Define starting value
# step = 0.1

# options_exc["fmt"] = 0      # Ensure automatic order reduction is enabled
# eng = startMatlabEngine()

# while options_exc["Fit"] <= 0.99:
#     # Run Matlab script
#     eng.SS_Fitting_Wave_Excitation(options_exc,tc,manual_order_SSexc,nargout=0) # nargout=0 is required for scripts or functions that do not return output variables.

#     # Read number of states from file
#     file_path = os.path.join(options_exc["FileName"] + ".ssexctn")
#     with open(file_path, 'r') as f:
#         lines = f.readlines()
#     m = [int(x) for x in lines[4].split('%')[0].split()]    # states per degree of freedom

#     suffix_str = "_".join(str(x) for x in m)
#     suffix_str = f"Fit{options_exc['Fit']:.1f}_{suffix_str}"    # Add the fit to file name

#     # Move and rename file to wanted destination
#     default_output = os.path.join(source_file_dir, "marin_semi.ssexctn")           
#     target_output = os.path.join(dest_dir, f"marin_semi_{suffix_str}.ssexctn")

#     if os.path.exists(default_output):
#         shutil.copy(default_output, target_output)

#     # Increase Fit by the step size
#     options_exc["Fit"] += step

# stopMatlabEngine(eng)


### --- For loop where all DoFs have same # states
# eng = startMatlabEngine()

# for i in range(1, 6):
#     base_n = [2*i, 2*i, 2*i, 2*i, 2*i, 2*i]                     # Order of the fit for each DoF
#     manual_order_SSexc = RedOrder_SSexc(base_n, options_exc)    # Order of the fit filtered as needed for input of the MATLAB function

#     # Run Matlab script
#     eng.SS_Fitting_Wave_Excitation(options_exc,tc,manual_order_SSexc,nargout=0) # nargout=0 is required for scripts or functions that do not return output variables.

#     # Move and rename file to wanted destination
#     suffix_str = "_".join(str(x) for x in base_n)
#     default_output = os.path.join(source_file_dir, "marin_semi.ssexctn")           
#     target_output = os.path.join(dest_dir, f"marin_semi_{suffix_str}.ssexctn")

#     if os.path.exists(default_output):
#         shutil.copy(default_output, target_output)

# stopMatlabEngine(eng)


