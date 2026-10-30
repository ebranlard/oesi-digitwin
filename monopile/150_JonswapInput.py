"""
Generate wave time series based on the Jonswap spectrum for different water depth

Note: we could rely on OpenFAST but this gives us some control over it and allows the generation of a WaveKin file if necessary.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from welib.hydro.spectra import jonswap
from welib.hydro.wavekin import elevation2d, wavenumber
from welib.tools.spectral import fft_wrap
from numpy.random import uniform, seed
seed(10)


SSCond=[]
SSCond.append({'Hs':8.1, 'Tp':12.7, 'h':34}) # Old
SSCond.append({'Hs':2.5, 'Tp':10  , 'h':34}) # Old
SSCond.append({'Hs':6  , 'Tp':12  , 'h':34}) # Old

for cond in SSCond:
    """ 
    h: wave elevation

    """
    # --- Parameters
    t  = np.linspace(0,600,6001) # time vector  [s]
    h  = cond['h']
    Hs = cond['Hs'] # Significant wave height [m]
    Tp = cond['Tp'] # Peak period [s]
    g  = 9.81 # Gravity[m/s2]  

    # --- Jonswap spectrum
    dt       = t[1]-t[0]                   # timestep [s]
    df       = 1/np.max(t)                 # step size for frequency
    fHighCut = 1/(dt)/2.                   # Highest frequency in calculations
    freq     = np.arange(df, fHighCut, df)
    S = jonswap(freq, Hs, Tp=Tp, g=9.81)

    # --- Keep less components
    Ikeep=S>np.max(S)/100
    freq=freq[Ikeep]
    S   =S [Ikeep]

    Ikeep=np.linspace(0,len(S)-1,20).astype(int)
    freq=freq[Ikeep]
    S   =S [Ikeep]
    df =freq[1]-freq[0]



    # --- Solve dispersion relation
    k = wavenumber(freq, h, g)

    # --- Compute wave elevation based on amplitudes and random phases
    eps = uniform(0,2*np.pi,len(freq))  # random phases between 0 and 2pi
    epsShift = eps + np.pi/2 
    a   = np.sqrt(2*S*df)  # wave amplitudes based on spectrum
    x = 0                  # longitudinal distance where wave is evaluated [m]
    eta = elevation2d(a, freq, k, eps, t, x)    

    # --- Write WaveKin file
#     WaveKinFile='Monopile/_Flex_Monopile/Waves/UserDefSimpleJonswap{}_new.wki'.format(str(int(h)))
#     with open(WaveKinFile, 'w') as f:
#         f.write('WAVEKININPUT\n')
#         f.write('#001\n')
#         f.write(' 1\n')
#         f.write('** Example for UserDef (SpecType=9, table of frequencies amplitude phases)\n')
#         f.write('   {}    WaterDepth\n'.format(h))
#         f.write(' 8.1966401 Pile Diameter\n')
#         f.write('       1     Number of displacements\n')
#         f.write('   0.000\n')
#         f.write('      26     Number of relative depth\n')
#         f.write(' {} \n'.format(' '.join(np.linspace(0,1,26).astype(str)) ))
#         f.write('    0.10     Time step\n')
#         f.write('  600.00     Duration\n')
#         f.write('USERDEF\n')
#         f.write('       {}    Number of frequency components, next is table of [freq, amp, phase]\n'.format(len(freq)))
#         for fi,ai,phi in zip(freq,a,epsShift):
#             f.write(' {} {} {} \n'.format(fi,ai,phi))
#         f.write('SLAMOFF')

    # --- Write HydroDyn Elev file
    SSCase = f'UserDefJonswap_Hs={Hs}_Tp={Tp}_h={int(h)}'
    print('Case:', SSCase)
    if h==34:
        HDElevFile =f'simulations_wt/Waves/{SSCase}.Elev'
    else:
        raise
    with open(HDElevFile, 'w') as f:
       f.write(' Time            Wave1Elev\n')
       f.write('(sec)                  (m)\n')
       for tt,etatt in zip(t,eta):
           f.write('{} {}\n'.format(tt,etatt))
    print('>>> Elv', HDElevFile)

    # --- Write DataFrame
    DFComp =f'simulations_wt/Waves/{SSCase}.csv'
    print('>>> CSV', DFComp)
    df = pd.DataFrame(data=np.column_stack((freq,a,eps)), columns=['Frequency_[Hz]','Amplitude_[m]','Phase_[rad]'])
    df.to_csv(DFComp, index=False)


    # --- Compute FFT of wave elevation
    f_fft, S_fft, Info = fft_wrap(t, eta, output_type='PSD', averaging='none')
         
    # --- Plots
    fig,axes = plt.subplots(2, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.30, wspace=0.20)

    ax=axes[0]
    ax.plot(t, eta)
    ax.tick_params(direction='in')
    ax.autoscale(enable=True, axis='both', tight=True)
    ax.set_xlabel('Time [s]')
    ax.set_ylabel(r'Wave elevation [m]')
    ax.set_title('Hydro - wave generation')

    ax=axes[1]
    ax.plot(f_fft, S_fft, '-', label='Generated')
    ax.plot(freq, S     , 'k', label='Jonswap')
    ax.legend()
    ax.set_xlabel('Frequency [Hz]')
    ax.set_ylabel(r'Spectral density [m$^2$ s]')
    ax.tick_params(direction='in')
    ax.autoscale(enable=True, axis='both', tight=True)


plt.show()


