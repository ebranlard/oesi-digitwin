# Load packages
import numpy as np
from openfast_toolbox.io import FASTOutputFile # type: ignore
from openfast_toolbox.io import FASTLinearizationFile # type: ignore

def rigidBodyMass(filepath):
    """
    This function finds the rigid body 6x6 mass matrix and it's inverse from an OpenFAST .lin file.
    This mass contains the infinite-frequency hydrodynamic added mass A_inf if Hydrodyne was activated.
    """

    lin = FASTLinearizationFile(filepath)

    # Extract the state-space matrices as numpy arrays
    # A = lin['A']
    B = lin['B']
    # C = lin['C']
    # D = lin['D']

    # Extract the lists of descriptions
    input_descriptions = lin['u_info']['Description']
    state_descriptions = lin['x_info']['Description']

    # Find input columns starting with "ED Platform"
    # We match the first 6 occurrences corresponding to Surge, Sway, Heave, Roll, Pitch, Yaw forces
    input_indices = [
        i for i, description in enumerate(input_descriptions) 
        if description.strip().startswith("ED Platform")
    ][:6]

    # Find state rows starting with "ED First time derivative" 
    # AND containing "Platform" or "Ptfm" to isolate the velocity states
    state_deriv_indices = [
        i for i, description in enumerate(state_descriptions) 
        if description.strip().startswith("ED First time derivative") and ("Platform" in description or "Ptfm" in description)
    ]

    # Extract the 6x6 B matrix subset
    B_subset = B[np.ix_(state_deriv_indices, input_indices)]

    # Find the rigid body 6x6 mass matrix by inverting B_subset
    M = np.linalg.inv(B_subset)
    # print(input_indices)
    # print(state_deriv_indices)
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

# Old, assumes 6DOFs
# def invMK(filepath):
#     """
#     This function finds M^(-1)*K from an OpenFAST .lin file.
#     Note: The sign is positive. If -M^(-1)*K is needed, the minus must be added when the function is called.
#     """

#     lin = FASTLinearizationFile(filepath)

#     # Extract the state-space A matrix as a numpy array
#     A = lin['A']

#     # Retrieve M^(-1)*K
#     inv_M_K = -A[6:12, 0:6]
#     return inv_M_K

# Old, assumes 6DOFs
# def invMC(filepath):
#     """
#     This function finds M^(-1)*C from an OpenFAST .lin file.
#     Note: The sign is positive. If -M^(-1)*C is needed, the minus must be added when the function is called.
#     """

#     lin = FASTLinearizationFile(filepath)

#     # Extract the state-space A matrix as a numpy array
#     A = lin['A']

#     # Retrieve M^(-1)*C
#     inv_M_C = -A[6:12, 6:12]
#     return inv_M_C

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