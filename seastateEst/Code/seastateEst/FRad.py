"""
This library includes the functions associated with computing the radiation Force in Cummins Equaion.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def Ainf(filepath, WtrDens):
    '''
    This function returns the infinite-frequency added mass matrix Ainf, given a WAMIT .1 file.
    The values in .1 are nondimensional, assuming WAMITULEN=1, they must be multiplied by the water density.
    WtrDens:    Water density
    '''
    # Open the WAMIT .1 file
    # Columns: Period, DOF_i, DOF_j, Added_Mass, Damping
    wamit_data = pd.read_csv(filepath, delim_whitespace=True, names=["per", "i", "j", "A", "B"])

    # Isolate the rows where PER == 0.0: These values are the entries for Ainf
    ainf_rows = wamit_data[wamit_data["per"] == 0.0]

    if ainf_rows.empty:
        raise ValueError("Could not find PER = 0.0 rows in the WAMIT file!")

   # Reconstruct the 6x6 A_inf matrix depending on i and j
    A_inf = np.zeros((6, 6))
    for _, row in ainf_rows.iterrows():
        i, j = int(row["i"]) - 1, int(row["j"]) - 1
        A_inf[i, j] = row["A"]

    # Dimensionalize A_inf
    A_inf = A_inf * WtrDens
    return A_inf

def RdtnKrnl(filepath):
    '''
    This function extracts the radiation kernel K of the impulse response function from an OpenFAST HD.sum file
    '''
    # Find the header line and read its column names
    with open(filepath, 'r') as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        if line.split()[:2] == ['n', 't']:
            header_line = i
            columns = line.split()
            break

    # Read numerical data, skipping the header and units lines
    df = pd.read_csv(
        filepath,
        sep=r'\s+',
        skiprows=header_line + 2,
        names=columns
    )

    # Keep only rows with a numeric time-step index
    df = df[pd.to_numeric(df['n'], errors='coerce').notna()]

    # Extract NumPy arrays
    # n = df['n'].to_numpy(dtype=int)
    # t = df['t'].to_numpy(dtype=float)
    # K = df.iloc[:, 2:].to_numpy(dtype=float)

    return df

def K_from_df(K_df, k):
    '''
    This function reconstructs the radiation kernel K at a given index n
    K_df:   Data frame whith entries: n[-], t[s], K(1,1)[kg/s^2], K(1,2)[kg/s^2], ..., K(2,1)[kg/s^2], ..., K(6,6)[kgm/s^2]
    k:   Index
    '''
    # Initialize an empty matrix
    K = np.zeros((6, 6))

    # Fill the matrix
    for i in range(1, 7):      # Loops 1 to 6 (Output DOFs)
        for j in range(1, 7):  # Loops 1 to 6 (Input DOFs)
            
            # Format 1: OpenFAST v5.0+ style -> "K(1,1)"
            v5_name = f"K({i},{j})"
            
            # Format 2: Pre-v5.0 style -> "K11"
            legacy_name = f"K{i}{j}"
            
            # Select whichever column name actually exists in the DataFrame
            if v5_name in K_df.columns:
                column_name = v5_name
            elif legacy_name in K_df.columns:
                column_name = legacy_name
            else:
                raise KeyError(f"Could not find a column named '{v5_name}' or '{legacy_name}' in the DataFrame.")
            
            # Extract the value at the target index
            K[i-1, j-1] = K_df.loc[k, column_name]
    return K

def Conv_Rdtn(K_df, qdot_history, RdtnTMax):
    '''
    This function computes the instantaneous 6-DOF memory force Fmem at the CURRENT time by performing the convolution.
    The function does not account for the negative sign Fmem has in Frad.
    K_df:           Data frame whith entries: n[-], t[s], K(1,1)[kg/s^2], K(1,2)[kg/s^2], ..., K(2,1)[kg/s^2], ..., K(6,6)[kgm/s^2]
    qdot_history:   Velocity time history. Numpy array of shape (6, Total_Sim_Steps)
    tmax:           Analysis time for wave radiation kernel calculations (sec)
    '''
    #TODO there might be a factor 2 that has to multiplied to Fmem because K from HD.sum is missing it.
    ### --- Initialize
    # Determine dt and make sure dt in K and qdot are identical
    t = K_df['t'].to_numpy(dtype=float) # TODO make sure dt in K and qdot match
    dt = t[1] - t[0]

    # MaxInd is the minimum of: length of K, length of qdot, steps within tmax
    MaxInd = min(K_df['n'].shape[0], qdot_history.shape[1], int(RdtnTMax/dt)) #TODO maybe add -1 for zero indexing in python
    print('MaxInd for convolution:', MaxInd)

    # Setting up qdot for specified time window
    qdot = qdot_history[:, -MaxInd:]

    # Loading K as a 3D array for the specified time window with shape: (Timestep, 6, 6)
    # Example: K(2,:,:) is the full 6x6 K at timestep 2
    K = K_df.iloc[0:MaxInd, 2:38].to_numpy(dtype=float).reshape(-1, 6, 6)

    # Initialize empty Fmem vector that will be filled for every time step k
    Fmem_k = np.zeros([6, MaxInd])

    ### --- Perform numerical convolution to determine the load contribution from wave radiation damping: 
    # First time step
    Fmem_k[:, 0] = -0.5 * K[MaxInd - 1] @ qdot[:, 0]

    # Last time step
    Fmem_k[:, MaxInd-1] = -0.5 * K[0] @ qdot[:, MaxInd-1]

    # Intermediate time steps
    # for k in range(1, MaxInd - 1):
    #     Fmem_k[:, k] = K[MaxInd - 1 - k] @ qdot[:, k]
    Fmem_k[:, 1:MaxInd-1] = -1 * np.einsum('kij,jk->ik', K[MaxInd-2:0:-1], qdot[:, 1:MaxInd-1])      # Much faster than the for loop for intermediate step

    # Compute the total memory force:
    Fmem = dt * np.sum(Fmem_k, axis=1, keepdims=True)

    return Fmem


# --- Testing the code

### Testing the convolution
# filepath_HDsum = r'C:\Jonah\Uni\PhD\OESI\OpenFAST_cases\IEA-15_Full_OpenFASTv5.0\main.HD.sum'
# filepath_HDsum = r'C:\Jonah\Uni\PhD\OESI\OpenFAST_cases\IEA-15_PlatformOnly\main.HD.sum'
# df = RdtnKrnl(filepath_HDsum)
# K = K_from_df(df, 1)

# qdot = np.array([
#     [10, 11, 12],
#     [20, 21, 22],
#     [30, 31, 32],
#     [40, 41, 42],
#     [50, 51, 52],
#     [60, 61, 62]
# ])

# Fmem = Conv_Rdtn(df, qdot, tmax=60)
# print('Fmem at current time tn:\n', Fmem)

# plt.figure()
# plt.plot(df['t'].to_numpy(dtype=float), df['K(3,3)'].to_numpy(dtype=float))
# # plt.plot(df['t'].to_numpy(dtype=float), df['K33'].to_numpy(dtype=float))
# all_K = df.iloc[:, 2:38].to_numpy(dtype=float).reshape(-1, 6, 6)

### Testing A_inf
# filepath = r"C:\Jonah\Uni\PhD\OESI\OpenFAST_cases\IEA-15_Full\HydroData\IEA-15-240-RWT-UMaineSemi.1"
# A_inf = Ainf(filepath, WtrDens=1025)
# print(A_inf)



