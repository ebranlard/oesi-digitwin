import numpy as np

def read_ssexctn(file_path):
    """
    This function reads .ssexctn files and returns the matrices for the state-space model.
    """
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    # Metadata
    heading = float(lines[1].split('%')[0])                 # Wave heading angle, Splits at % and takes everything to the left, then turns it into a float
    tc = float(lines[2].split('%')[0])                      # time offset (tc)
    n_states = int(lines[3].split('%')[0])                  # total # of states
    m = [int(x) for x in lines[4].split('%')[0].split()]    # states per degree of freedom

    # A-matrix
    # always starts at line 5 and has dimensions n_states x n_states
    A = np.array([line.split() for line in lines[5 : 5 + n_states]], dtype=float)

    # B-matrix
    # has dimensions n_states x 1
    start_b = 5 + n_states
    B = np.array([line.split() for line in lines[start_b : start_b + n_states]], dtype=float)
    
    # C-matrix
    # has dimensions 6 x n_states
    # start_c = start_b + n_states
    # C = np.array([line.split() for line in lines[start_c : start_c + 6]], dtype=float)
    start_c = start_b + n_states
    c_rows = []
    for line in lines[start_c : start_c + 6]:
        split_line = line.split()
        # If a line is empty, it means that DOF has 0 states; replace with zeros
        if not split_line:
            c_rows.append([0.0] * n_states)
        else:
            c_rows.append(split_line)
            
    C = np.array(c_rows, dtype=float)

    print('Reading .ssexctn file with header:')
    print('heading: ',heading)
    print('tc: ',tc)
    print('n: ',n_states)
    print('m: ',m)
    return A, B, C, n_states, m

def read_ss(file_path):
    """
    This function reads .ss files and returns the matrices for the state-space model.
    """
    with open(file_path, 'r') as f:
        lines = f.readlines()

    # Metadata
    dof = [int(x) for x in lines[1].split('%')[0].split()]  # Enabled DOFs following the order Surge, Sway, Heave, Roll, Pitch and Yaw
    n_states = int(lines[2].split('%')[0])                  # Total # of radiation states
    m = [int(x) for x in lines[3].split('%')[0].split()]    # radiation states per degree of freedom
    
    # A-matrix
    # always starts at line 4 and has dimensions n_states x n_states
    A = np.array([line.split() for line in lines[4 : 4 + n_states]], dtype=float)

    # B-matrix
    # has dimensions n_states x 6
    start_b = 4 + n_states
    B = np.array([line.split() for line in lines[start_b : start_b + n_states]], dtype=float)

    # C-matrix
    # has dimensions 6 x n_states
    start_c = start_b + n_states
    C = np.array([line.split() for line in lines[start_c : start_c + 6]], dtype=float)

    print('Reading .ss file with header:')
    print('DOF: ',dof)
    print('n: ',n_states)
    print('m: ',m)
    return A, B, C, n_states

def reduceExctn(dof_dict, A, B, C, n_states, m):
    """
    Reduces wave excitation state-space matrices by mapping active 
    DOFs to their sequential internal state indices.
    
    INPUTS:
      dof_dict: Dictionary mapping DOF names to boolean activation flags: dof = {"Surge": True,"Sway": False,"Heave": True,"Roll": False,"Pitch": True,"Yaw": False,}
      A: Full state matrix (n_states x n_states)
      B: Full input matrix (n_states x 1)
      C: Full output matrix (6 x n_states)
      n_states: Total number of full internal states (sum of m)
      m: List of 6 integers representing states per DOF
    """
    active_dof_list = list(dof_dict.values())
    
    # 1. Build an index mask for the internal states
    state_indices = []
    current_idx = 0
    
    for i, m_i in enumerate(m):
        # If this mechanical DOF is active, keep its internal states
        if active_dof_list[i]:
            state_indices.extend(range(current_idx, current_idx + m_i))
        current_idx += m_i
        
    # 2. Slice the matrices using the tracked state indices
    # A is sliced symmetrically in both dimensions
    A_red = A[np.ix_(state_indices, state_indices)]
    
    # B rows correspond to internal states
    B_red = B[state_indices, :]
    
    # C rows correspond to 6 DOFs; columns correspond to internal states
    C_red = C[active_dof_list, :][:, state_indices]
    
    # 3. Track the new reduced state count
    n_red = len(state_indices)
    
    return A_red, B_red, C_red, n_red

############################
# Testing the code
# A_r, B_r, C_r, n_r = read_ss(r'C:\Jonah\Uni\PhD\OESI\Starting_documents\SS_Excitation_Fitting\alpha-Release\HydroData\OC4 Semi\marin_semi.ss')

# A_exc, B_exc, C_exc, n_exc, m_exc = read_ssexctn(r'C:\Jonah\Uni\PhD\OESI\seastateEst\Code\ss_fitting_modified\HydroData\marin_semi.ssexctn')
# A_exc, B_exc, C_exc, n_exc, m_exc = read_ssexctn(r'C:\Jonah\Uni\PhD\OESI\seastateEst\OpenFAST_cases\Parametric_study\Lock_Rad\HydroData\marin_semi_Fit0.7_6_10_4_4_6_6.ssexctn')


# dof = {
#     "Surge": True,
#     "Sway": False,
#     "Heave": True,
#     "Roll": False,
#     "Pitch": True,
#     "Yaw": False,
# }
# A_red, B_red, C_red, n_red = reduceExctn(dof, A_exc, B_exc, C_exc, n_exc, m_exc)