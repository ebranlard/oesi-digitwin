# Enable Autoreload
import notebook_setup # type: ignore
notebook_setup.enable_autoreload()

import os
import glob
import subprocess
import shutil
import matplotlib.pyplot as plt

import openfast_toolbox.io as io # type: ignore
# from ss_fitting_modified import run_matlab

def files_needed(base_path, case_suffix, lock):
    '''
    This function creates temporary PotFiles for current iteration.
    case_suffix:    Can be a string like '10_10_5_8_8_3' or list like [10, 10, 5, 8, 8, 3] indicating the number of states per Dof
    lock='ss':      Each case has a unique .ss file; .ssexctn is shared.
    lock='ssexctn': Each case has a unique .ssexctn file; .ss is shared.
    '''
    # Make sure case_suffix always is a string of this format 10_10_5_8_8_3, no matter how it was passed on 
    if isinstance(case_suffix, (list, tuple)):
        case_suffix = "_".join(str(x) for x in case_suffix)
    else:
        case_suffix = str(case_suffix)

    # Names of files needed current_files['']
    current_files = dict()
    current_files['ssexctn'] = os.path.join(base_path, "HydroData", f"marin_semi_{case_suffix}.ssexctn")
    current_files['ss'] = os.path.join(base_path, "HydroData", f"marin_semi_{case_suffix}.ss")
    current_files['hst']     = os.path.join(base_path, "HydroData", f"marin_semi_{case_suffix}.hst")
    current_files['1']       = os.path.join(base_path, "HydroData", f"marin_semi_{case_suffix}.1")
    current_files['3']       = os.path.join(base_path, "HydroData", f"marin_semi_{case_suffix}.3")

    # Define where the master files live
    shared_ssexctn_master = os.path.join(base_path, "HydroData", "marin_semi.ssexctn")
    shared_ss_master = os.path.join(base_path, "HydroData", "marin_semi.ss")
    shared_hst_master = os.path.join(base_path, "HydroData", "marin_semi.hst")
    shared_1_master = os.path.join(base_path, "HydroData", "marin_semi.1")
    shared_3_master = os.path.join(base_path, "HydroData", "marin_semi.3")

    # Copy the WAMIT files, these are always shared
    shutil.copy(shared_hst_master, current_files['hst'])
    shutil.copy(shared_1_master, current_files['1'])
    shutil.copy(shared_3_master, current_files['3'])

    # Copy either .ss or .ssexctn file depending on which file type is locked
    if lock == 'ssexctn':
        shutil.copy(shared_ssexctn_master, current_files['ssexctn'])
    elif lock == 'ss':
        shutil.copy(shared_ss_master, current_files['ss'])
    else:
        print(f"Error: Unknown lock parameter '{lock}'. Must be 'ss' or 'ssexctn'.")
        exit(1)

    return current_files

def delete_files(fileName, lock):
    '''
    This function deletes the temporary PotFiles with the correct name, if .ssexctn is to be reused.
    fileName: dict()
    '''
    # Remove temporary WAMIT files
    if os.path.exists(fileName['hst']):
        os.remove(fileName['hst'])
    if os.path.exists(fileName['1']):
        os.remove(fileName['1'])
    if os.path.exists(fileName['3']):
        os.remove(fileName['3'])

    # Remove temporary .ss or .ssexctn file depending on which file type is locked
    if lock == 'ssexctn':
        if os.path.exists(fileName['ssexctn']):
            os.remove(fileName['ssexctn'])
    elif lock == 'ss':
        if os.path.exists(fileName['ss']):
            os.remove(fileName['ss'])
    else:
        print(f"Error: Unknown lock parameter '{lock}'. Must be 'ss' or 'ssexctn'.")
        exit(1)
    return

def executeOpenFAST(openfast_exe, fst_file, base_path):
    '''
    This function is used to execute OpenFAST
    openfast_exe:   filepath of the .exe
    fst_file:       filepath of the main.fst
    '''
    print("\n--- Executing OpenFAST ---")
    print(f"Running command: {openfast_exe} {fst_file}")

    try:
        # Runs the executable and waits for it to complete
        subprocess.run([openfast_exe, fst_file], check=True, cwd=base_path)
        print("OpenFAST execution finished successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error: OpenFAST failed with return code {e.returncode}")
        exit(1)
    except FileNotFoundError:
        print(f"Error: Could not find the executable '{openfast_exe}'. Ensure it is in your PATH or specify the full path.")
        exit(1)
    return

def rename_outputs(base_path, case_suffix):
    '''
    This function finds default main.* outputs and renames them to match the case number.
    '''
    for ext in ['.out', '.outb', '.sum']:
        default_output = os.path.join(base_path, f"main{ext}")           
        target_output = os.path.join(base_path, f"main_{case_suffix}{ext}")
        
        if os.path.exists(default_output):
            shutil.move(default_output, target_output)
            print(f"Saved output file as: {target_output}")
    return

### --- Running the code

# --- File and Path Configurations ---
# base_path = r'C:\Jonah\Uni\PhD\OESI\OpenFAST_cases\Test_forrunningPython_DONOTDELETE'
base_path = r'C:\Jonah\Uni\PhD\OESI\seastateEst\OpenFAST_cases\Parametric_study\Lock_Rad'       # base_path is where the OpenFAST input files are

openfast_exe = os.path.join(base_path, "openfast_x64.exe")           # OpenFAST executable
fst_file = os.path.join(base_path, "main.fst")                       # Main OpenFAST configuration driver file
hydrodyn_file = os.path.join(base_path, "HydroDyn.dat")              # Primary HydroDyn input file

# Ensure correct settings in HydroDyn file
hd_data = io.FASTInputFile(hydrodyn_file)
hd_data['ExctnMod'] = 2
hd_data['RdtnMod'] = 2

# Reuse .ss or .ssexctn file (this file only needs to exist once and will be copied and deleted on each iteration)
# lock = 'ssexctn'
lock = 'ss'

### --- Define the list of case suffixes to run
### Option 1: Create list of cases manually
cases = ["10_15_7_12_6_15",] # 15_15_10_15_10_15 is the 99% fit

### Option 2: Create list of cases automatically based on what is in the folder
# search_path = os.path.join(base_path, "HydroData", "*.ssexctn")
# files = glob.glob(search_path)

# cases = []
# for file_path in files:
#     # Get just the filename (e.g., "marin_semi_2_2_2_2_2_2.ssexctn")
#     filename = os.path.basename(file_path)

#     # Remove the extension (leaves "marin_semi_2_2_2_2_2_2")
#     name_without_ext = os.path.splitext(filename)[0]

#     # Split by underscore to isolate the components
#     parts = name_without_ext.split("_")

#     # Check if the 7th element from the right starts with "Fit"
#     if len(parts) >= 7 and parts[-7].startswith("Fit"):
#         case_string = "_".join(parts[-7:])  # Extracts "Fit0.3_2_2_2_2_2_2"
#     else:
#         case_string = "_".join(parts[-6:])  # Backward compatibility: extracts "2_2_2_2_2_2"

#     cases.append(case_string)
# cases.sort()

## Run OpenFAST in a loop
for suffix in cases:
    # Check if the individual case item is a list or tuple
    if isinstance(suffix, (list, tuple)):
        suffix_str = "_".join(str(x) for x in suffix)
    else:
        suffix_str = str(suffix)
    # ----------------------------------

    print(f"\n=========================================")
    print(f"STARTING SIMULATION CASE: {suffix_str}")
    print(f"=========================================")

    # Update the PotFile root name for the current iteration
    hd_data['PotFile'] = f'"HydroData/marin_semi_{suffix_str}"'
    hd_data.write(hydrodyn_file)

    # Create temporary files using the standardized suffix string
    current_files = files_needed(base_path, suffix_str, lock)

    # Run OpenFAST
    executeOpenFAST(openfast_exe, fst_file, base_path)
    
    # Rename outputs using the base path and suffix string
    rename_outputs(base_path, suffix_str)

    # Delete temporary files
    delete_files(current_files, lock)

### Load OpenFAST .out files
fast_data = {}
for case in cases:
    # Construct the dynamic file path
    file_path = os.path.join(base_path, f"main_{case}.out")
    
    # Load the file and store it in the dictionary
    try:
        fast_data[case] = io.FASTOutputFile(file_path).toDataFrame()
        print(f"Successfully loaded: {case}.out")
    except Exception as e:
        print(f"Error loading {case}: {e}")

### --- Plots

### 6 DOF plot
# 1. Define the 6 DOFs with their OpenFAST channel names and display labels
dofs = {
    'Surge': ('PtfmSurge_[m]', 'Surge [m]'),
   #'Sway':  ('PtfmSway_[m]', 'Sway [m]'),
    'Heave': ('PtfmHeave_[m]', 'Heave [m]'),
    #'Roll':  ('PtfmRoll_[deg]', 'Roll [deg]'),
    'Pitch': ('PtfmPitch_[deg]', 'Pitch [deg]'),
   #'Yaw':   ('PtfmYaw_[deg]', 'Yaw [deg]')
}

# 2. Create a figure with #DOF subplots stacked vertically (up tp 6 rows, 1 column)
fig, axes = plt.subplots(len(dofs), 1, figsize=(12, 8), sharex=True)

# 3. Loop through each DOF and its assigned subplot axis
for i, (dof_name, (channel, y_label)) in enumerate(dofs.items()):
    ax = axes[i]
    
    # Loop through each case in your fast_data dictionary
    for case, data in fast_data.items():
        time = data['Time_[s]'].values
        dof_data = data[channel].values
        
        # Only add labels to the first subplot so the legend doesn't repeat 6 times
        label = f'OpenFAST - {case}' if i == 0 else ""
        ax.plot(time, dof_data, label=label, linewidth=1.2)
    
    # Format each individual subplot
    ax.set_ylabel(y_label)
    ax.grid(True, linestyle=':', alpha=0.6)

# 4. Global figure formatting
axes[0].legend(loc='upper right', bbox_to_anchor=(1.0, 1.35), ncol=3) # Put legend above the top plot
axes[-1].set_xlabel('Time [s]') # Label only the bottom X-axis
plt.suptitle('Platform 6-DOF Comparison Across Cases', y=0.98, fontsize=14, fontweight='bold')

plt.tight_layout()
plt.show()