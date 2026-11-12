# Load packages
import numpy as np
from openfast_toolbox.io import FASTOutputFile # type: ignore
from openfast_toolbox.io import FASTLinearizationFile # type: ignore

# def rigidBodyMass(filepath):
#     """
#     This function finds the rigid body 6x6 mass matrix and it's inverse from an OpenFAST .lin file.
#     This mass contains the infinite-frequency hydrodynamic added mass A_inf if Hydrodyne was activated.
#     LinInputs = 2 in main.fst for this to work
#     """

#     lin = FASTLinearizationFile(filepath)

#     # Extract the state-space matrices as numpy arrays
#     # A = lin['A']
#     B = lin['B']
#     # C = lin['C']
#     # D = lin['D']

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

def rigidBodyMass(filepath):
    """
    This function finds the rigid body mass matrix and its inverse from an OpenFAST .lin file.
    It automatically detects which DOFs are enabled by looking at the states, and maps them
    to the X/Y/Z Force/Moment input descriptions to guarantee a square matrix.

    LinInputs = 2 in main.fst for this to work
    """
    lin = FASTLinearizationFile(filepath)
    B = lin['B']

    input_descriptions = lin['u_info']['Description']
    state_descriptions = lin['x_info']['Description']

    # 1. Find the velocity state rows that are actually present
    state_deriv_indices = [
        i for i, desc in enumerate(state_descriptions) 
        if desc.strip().startswith("ED First time derivative") and ("Platform" in desc or "Ptfm" in desc)
    ]

    # Standard keywords OpenFAST uses for platform DOFs in the states
    dof_keywords = ['surge', 'sway', 'heave', 'roll', 'pitch', 'yaw']
    
    # 2. Extract which DOFs are active from the states
    active_dofs = []
    for idx in state_deriv_indices:
        desc_lower = state_descriptions[idx].lower()
        for kw in dof_keywords:
            if kw in desc_lower:
                active_dofs.append(kw)
                break

    # 3. Explicitly map active DOFs to your specific input string signatures
    input_indices = []
    for kw in active_dofs:
        if kw == 'surge':   match_str = "x force"
        elif kw == 'sway':  match_str = "y force"
        elif kw == 'heave': match_str = "z force"
        elif kw == 'roll':  match_str = "x moment"
        elif kw == 'pitch': match_str = "y moment"
        elif kw == 'yaw':   match_str = "z moment"
        
        idx_u = [
            i for i, desc in enumerate(input_descriptions)
            if desc.strip().startswith("ED Platform") 
            and match_str in desc.lower()
        ]
        
        if idx_u:
            input_indices.append(idx_u[0])
        else:
            raise ValueError(f"Could not find matching ED Platform input for signature: '{match_str}'")

    # Extract the square B matrix subset
    B_subset = B[np.ix_(state_deriv_indices, input_indices)]

    # Compute the mass matrix
    M = np.linalg.inv(B_subset)

    return M, B_subset

def invMK(filepath, dof_dict):
    """
    This function finds M^(-1)*K from an OpenFAST .lin file.
    It accounts for the activated dofs. The dofs must match what was set in the OpenFAST ED file.
    Note: The sign is positive. If -M^(-1)*K is needed, the minus must be added when the function is called.
    """
    # Count how many DOFs are set to True
    n_dofs = sum(dof_dict.values())
    
    lin = FASTLinearizationFile(filepath)
    A = lin['A']

    # Slice dynamically based on the number of active DOFs
    inv_M_K = -A[n_dofs : 2 * n_dofs, 0 : n_dofs]
    return inv_M_K

def invMC(filepath, dof_dict):
    """
    This function finds M^(-1)*C from an OpenFAST .lin file.
    It accounts for the activated dofs. The dofs must match what was set in the OpenFAST ED file.
    Note: The sign is positive. If -M^(-1)*C is needed, the minus must be added when the function is called.
    """
    # Count how many DOFs are set to True
    n_dofs = sum(dof_dict.values())
    
    lin = FASTLinearizationFile(filepath)
    A = lin['A']

    # Slice dynamically based on the number of active DOFs
    inv_M_C = -A[n_dofs : 2 * n_dofs, n_dofs : 2 * n_dofs]
    return inv_M_C

########## Testing the code
# filepath = r"C:\Jonah\Uni\PhD\OESI\OpenFAST_cases\IEA-15_PlatformOnly_LINEAR\No_waves\main.1.lin"
# # filepath = r"C:\Jonah\Uni\PhD\OESI\OpenFAST_cases\IEA-15_Full_LINEAR\main.1.lin"

# dof = {
#     "Surge": False,
#     "Sway": False,
#     "Heave": True,
#     "Roll": False,
#     "Pitch": False,
#     "Yaw": False,
# }

# inv_M_K = invMK(filepath, dof)
# inv_M_C = invMC(filepath, dof)