import numpy as np
from scipy.linalg import expm

def build_mixed_measurement_matrices(A, B, dof):
    """
    Generically constructs Yx and Yu matrices for a mixed configuration:
    - Active translational DOFs (Surge, Sway, Heave) -> ACCELERATION rows
    - Active rotational DOFs (Roll, Pitch, Yaw)     -> INCLINATION (Position) rows
    
    Handles matrices A and B that are ALREADY pre-reduced to the active DOFs.
    """
    active_dof_flags = list(dof.values())  # e.g., [True, False, True, False, True, False]
    n_active_dofs = sum(active_dof_flags)
    n_total_states = A.shape[0]

    Yx_rows = []
    Yu_rows = []
    
    # Extract the pre-reduced acceleration rows from A and B
    A_accel_block = A[n_active_dofs : 2*n_active_dofs, :]
    B_accel_block = B[n_active_dofs : 2*n_active_dofs, :]
    
    reduced_idx = 0  # Tracks the state index
    
    # Loop over the original 6 DOFs to evaluate the physical configuration
    for i in range(6):
        if active_dof_flags[i]:
            
            # --- Condition A: Translational DOFs (Surge=0, Sway=1, Heave=2)
            if i < 3:
                # Appends the corresponding row from the pre-reduced acceleration blocks
                Yx_rows.append(A_accel_block[reduced_idx, :])
                Yu_rows.append(B_accel_block[reduced_idx, :])
                
            # --- Condition B: Rotational DOFs (Roll=3, Pitch=4, Yaw=5)
            else:
                # Appends an identity selector row pointing directly to the position state index
                yx_row_pos = np.zeros(n_total_states)
                yx_row_pos[reduced_idx] = 1.0  
                
                Yx_rows.append(yx_row_pos)
                Yu_rows.append(np.zeros(n_active_dofs)) # Direct force throughput on position is zero
            
            # Always increment the state index tracker when a DOF is active
            reduced_idx += 1
            
    # Stack everything cleanly into final 2D NumPy arrays (Shape: 3 x n_total_states)
    Yx_gen = np.vstack(Yx_rows)
    Yu_gen = np.vstack(Yu_rows)
    
    return Yx_gen, Yu_gen

def discretize_exact(A, B, dt):
    n = A.shape[0]
    m = B.shape[1]

    AB = np.zeros((n + m, n + m))
    AB[:n, :n] = A
    AB[:n, n:] = B

    E = expm(AB * dt)

    Ad = E[:n, :n]
    Bd = E[:n, n:]

    return Ad, Bd