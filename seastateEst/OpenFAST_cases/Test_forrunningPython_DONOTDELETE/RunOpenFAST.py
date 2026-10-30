import os
import subprocess
import shutil
import openfast_toolbox.io as io # type: ignore

def files_needed(case_num, lock):
    '''
    This function creates temporary PotFiles for current iteration.
    lock='ss':      Each case has a unique .ss file; .ssexctn is shared.
    lock='ssexctn': Each case has a unique .ssexctn file; .ss is shared.
    '''
    # Names of files needed current_files['']
    current_files = dict()
    current_files['ssexctn'] = f"HydroData/marin_semi_{case_num}.ssexctn"
    current_files['ss'] = f"HydroData/marin_semi_{case_num}.ss"
    current_files['hst']     = f"HydroData/marin_semi_{case_num}.hst"
    current_files['1']       = f"HydroData/marin_semi_{case_num}.1"
    current_files['3']       = f"HydroData/marin_semi_{case_num}.3"

    # Define where the master files live
    shared_ssexctn_master = "HydroData/marin_semi.ssexctn"
    shared_ss_master = "HydroData/marin_semi.ss"
    shared_hst_master = "HydroData/marin_semi.hst"
    shared_1_master = "HydroData/marin_semi.1"
    shared_3_master = "HydroData/marin_semi.3"

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

def executeOpenFAST(openfast_exe, fst_file):
    '''
    This function is used to execute OpenFAST
    openfast_exe:   filepath of the .exe
    fst_file:       filepath of the main.fst
    '''
    print("\n--- Executing OpenFAST ---")
    print(f"Running command: {openfast_exe} {fst_file}")

    try:
        # Runs the executable and waits for it to complete
        subprocess.run([openfast_exe, fst_file], check=True)
        print("OpenFAST execution finished successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error: OpenFAST failed with return code {e.returncode}")
        exit(1)
    except FileNotFoundError:
        print(f"Error: Could not find the executable '{openfast_exe}'. Ensure it is in your PATH or specify the full path.")
        exit(1)
    return

def rename_outputs(case_num):
    '''
    This function finds default main.* outputs and renames them to match the case number.
    '''
    for ext in ['.out', '.outb', '.sum']:
        default_output = f"main{ext}"           # OpenFAST will create a main.out file on every run because my driver file is main.fst 
        target_output = f"main_{case_num}{ext}"
        
        if os.path.exists(default_output):
            shutil.move(default_output, target_output)
            print(f"Saved output file as: {target_output}")
    return

### --- Running the code

# --- 1. File and Path Configurations ---
openfast_exe = "openfast_x64.exe"           # OpenFAST executable
fst_file = "main.fst"                       # Main OpenFAST configuration driver file
hydrodyn_file = "HydroDyn.dat"              # Primary HydroDyn input file

# Ensure correct settings in HydroDyn file
hd_data = io.FASTInputFile(hydrodyn_file)
hd_data['ExctnMod'] = 2
hd_data['RdtnMod'] = 2

# Reuse .ss or .ssexctn file (this file only needs to exist once and will be copied and deleted on each iteration)
lock = 'ssexctn'
# lock = 'ss'

for i in range(0, 2):
    print(f"\n=========================================")
    print(f"STARTING SIMULATION CASE: i = {i}")
    print(f"=========================================")

    # Update the PotFile root name for the current iteration
    hd_data['PotFile'] = f'"HydroData/marin_semi_{i}"'
    hd_data.write(hydrodyn_file)

    # Create temporary files 
    current_files = files_needed(i, lock)

    # Run OpenFAST
    executeOpenFAST(openfast_exe, fst_file)
    rename_outputs(i)

    # Delete temporary files
    delete_files(current_files, lock)