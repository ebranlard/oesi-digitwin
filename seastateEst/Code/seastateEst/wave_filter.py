import numpy as np

def PM_waves(H_s,T_p,f):
    """
    This function returns the Pierson Moskowitz spectrum.
    """
    f_p = 1/T_p                                       # peak frequency
    S_PM = 0.3125*H_s**2*f_p**4*f**(-5)*np.exp(-1.25*(f_p/f)**4)
    return S_PM

def JONSWAP_waves(H_s,T_p,f):
    """
    This function returns the JONSWAP spectrum.
    """
    f_p = 1/T_p                                       # peak frequency

    # Determine sigma
    sigma = np.where(f <= f_p, 0.07, 0.09)
    
    # Determine alpha
    alpha = np.exp(- (f-f_p)**2 / (2*sigma**2*f_p**2))

    # Determine gamma
    k = T_p/np.sqrt(H_s)
    if k <= 3.6:
        gamma = 5
    elif 3.6 < k <= 5:
        gamma = np.exp(5.75 - 1.15*k)
    else:
        gamma = 1
    
    # Determine C
    C = 1 - 0.287*np.log(gamma)

    # Compute Pierson Moskowitz spectrum
    S_PM = PM_waves(H_s,T_p,f)

    # Compute JONSWAP spectrum
    S_JONSWAP = C * S_PM * gamma**alpha
    
    return S_JONSWAP

def waveFilter2Ord(zeta, H_s, T_p, f):
    """
    This function computes the second order wave filter.
    The input f is just a vector like: f = np.linspace(1e-5,0.5,1000)
    """
    f_0 = 1/T_p

    # Find k_w by matching the peak of S_eta and S_JONSWAP
    Smax = JONSWAP_waves(H_s,T_p,f_0)
    k_w = 8*np.pi**2*zeta*f_0**2*np.sqrt(Smax)

    # Compute S_eta #TODO: might have to be done in terms of omega, not in terms of f
    S_eta = k_w**2 / (
        16*np.pi**4*((f_0**2-f**2)**2 + (2*zeta*f_0*f)**2)
    )
    return S_eta

## Testing

# # S_eta in terms of omega
# omega = np.linspace(0, 2, 100)
# S_eta = k_w**2 / (
#     (omega_0**2 - omega**2)**2 + (2*zeta*omega_0*omega)**2
# )