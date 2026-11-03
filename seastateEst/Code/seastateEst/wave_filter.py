import numpy as np
import matplotlib.pyplot as plt

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

def waveFilter2Ord(zeta, H_s, T_p, f=None):
    """
    This function computes the second order wave shaping filter.

    INPUTS:
          zeta: parameter, typically 0.05 <= zeta <= 0.15
          H_s: Significant wave height [m]
          T_p: Peak spectral period [s]
          f: vector like: f = np.linspace(1e-5,0.5,1000);   If f is None, S_eta will be returned as 0.
        OUTPUTS:
          S_eta: PSD of wave shaping filter [m^2 s]
          k_w: parameter [-]
          omega_0: Peak spectral angular frequency [rad/s]
    """
    f_0 = 1/T_p
    omega_0 = 2 * np.pi * f_0

    # Find k_w by matching the peak of S_eta and S_JONSWAP
    Smax = JONSWAP_waves(H_s,T_p,f_0)
    k_w = 8*np.pi**2*zeta*f_0**2*np.sqrt(Smax)

    # Compute S_eta(f) and check if f is provided #TODO: might have to be done in terms of omega, not in terms of f
    if f is None:
        S_etaf = 0.0
    else:
        # Compute S_eta in terms of f
        S_etaf = k_w**2 / (
            16*np.pi**4*((f_0**2-f**2)**2 + (2*zeta*f_0*f)**2)
        )

    return S_etaf, k_w, omega_0

## Testing

# H_s = 1.1
# T_p = 8.52
# f = np.linspace(1e-5,0.5,1000)
# zeta = 0.1

# S_JS = JONSWAP_waves(H_s,T_p,f)
# S_eta, k_w, omega_0 = waveFilter2Ord(zeta, H_s, T_p, f)

# plt.figure()
# plt.plot(f,S_JS, label='JS')
# plt.plot(f, S_eta, label='2nd order filter')
# plt.xlabel('f in Hz')
# plt.legend()