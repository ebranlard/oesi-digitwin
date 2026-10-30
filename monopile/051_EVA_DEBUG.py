"""

"""

def basicEVA(A):
    # Basic EVA
    D,Q    = np.linalg.eig(A)
    Lambda = np.diag(D)
    v      = np.diag(Lambda)

    # Selecting eigenvalues with positive imaginary part (frequency)
    if True:
        Ipos = np.imag(v)>0
        Q = Q[:,Ipos]
        v = v[Ipos]

    # Frequencies and damping based on compled eigenvalues
    omega_0 = np.abs(v)              # natural cylic frequency [rad/s]
    freq_d  = np.imag(v)/(2*np.pi)   # damped frequency [Hz]
    zeta    = - np.real(v)/omega_0   # damping ratio
    freq_0  = omega_0/(2*np.pi)      # natural frequency [Hz]
    # Sorting
    if True:
        I = np.argsort(freq_0)
        freq_d = freq_d[I]
        freq_0 = freq_0[I]
        zeta   = zeta[I]
        Q      = Q[:,I]
    return freq_d, zeta, Q, freq_0 

# ---
# fstFile = 'simulations_wt/00_EVA/OF_F3T1S1_H0A0.fst'


fstFiles=[]
# fstFiles += ['C:/Users/ebranlard/Documents/Work/2024-10-OESI-Digitwin/DigiTwinMonopile/code5_wt/simulations/00_EVA_DEBUG/OF_F2_NoRNA.fst']

# fstFiles += ['simulations/00_EVA_DEBUG/OF_F2_NoRNA_NoCM.fst']
# fstFiles += ['simulations/00_EVA_DEBUG/OF_F2_NoRNA.fst']
# fstFiles += ['simulations/00_EVA_DEBUG/OF_F2T1_NoRNA.fst']

# fstFiles += ['simulations/00_EVA_DEBUG/OF_F2T0_NoRNA.fst']


# fstFiles += ['simulations/00_EVA/OF_F2T1S0_H0A0_NoRot.fst']
# fstFiles += ['simulations/00_EVA/OF_F3T1S0_H0A0_NoRot.fst'] # <<<<

PWD = 'C:/Users/ebranlard/Documents/Work/2024-10-OESI-Digitwin/DigiTwinMonopile/code5_wt/'
import numpy as np
import matplotlib.pyplot as plt
from welib.essentials import *
from welib.yams.models.MTNSB import FASTmodel2MTNSB
from welib.yams.models.TNSB_FAST import FASTmodel2TNSB
from welib.weio.fast_linearization_file import FASTLinearizationFile
from welib.kalman.kalman import BuildSystem_Linear_MechOnly 
# fstFiles += ['C:/Users/ebranlard/Documents/Work/2024-10-OESI-Digitwin/DigiTwinMonopile/code5_wt/simulations/00_EVA_DEBUG/OF_F3T0_NoRNA_CB0.fst']
# fstFiles += ['C:/Users/ebranlard/Documents/Work/2024-10-OESI-Digitwin/DigiTwinMonopile/code5_wt/simulations/00_EVA_DEBUG/OF_F3T0_NoRNA_CB2.fst']

# fstFiles += ['simulations/00_EVA/OF_F3T0S0_H0A0_NoRot.fst'] 
# fstFiles += ['simulations/00_EVA/OF_F3T1S0_H0A0_NoRot.fst'] 
fstFiles += ['simulations/00_EVA/OF_F2T1S0_H0A0_NoRot.fst'] 

# fstFiles += ['simulations/00_EVA/OF_F3T1S1_H0A0_NoRot_CB0.fst'] 
# fstFiles += ['simulations/00_EVA/OF_F3T1B100_H0A0_NoRot.fst'] 

# fstFiles += ['simulations/00_EVA_BIG/OF_F3T2B100_H0A0_NoRot.fst'] 
# fstFiles += ['simulations/00_EVA_BIG/OF_F3T2B100_H0A0_NoRot_CB0.fst'] 

# fstFiles += ['simulations/00_EVA_BIG/OF_F3T2B111_H0A0_NoRot.fst'] 
# fstFiles += ['simulations/00_EVA_BIG/OF_F4T4B111_H0A0_NoRot.fst'] 

# --- Debug twr
# fstFiles += ['simulations/00_EVA_DEBUG_Twr/OF_F0T2S0_H0A0_NoRot.fst'] 


# fstFiles += ['simulations/00_EVA_DEBUG_Twr/OF_F0T2B100_H0A0_NoRot_NoTilt.fst'] 

# fstFiles += ['simulations/00_EVA_DEBUG_Twr/OF_F0T2B100_H0A0_NoRot_NoTilt_NoGrav.fst'] 


for fstFile in fstFiles:
    print('----------------------------------------------------------------')
    print(f'------------------ FST FILE {fstFile} ---------------------')
    print('----------------------------------------------------------------')


    linFile = fstFile.replace('.fst','.1.lin')
    yamlFile = fstFile.replace('.fst','.SD.sumpy.yaml')

    # --- Loading a "YAMS" model
    shapes_sub=[]
    shapes_bld=[]
    shapes_twr=[]
    if 'F3' in fstFile or 'F2' in fstFile:
        shapes_sub=[0,4]
    elif 'F4' in fstFile:
        shapes_sub=[0,1, 3,4]

    if 'T1' in fstFile:
        shapes_twr=[0]
    elif 'T2' in fstFile:
        shapes_twr=[0,1]
    elif 'T4' in fstFile:
        shapes_twr=[0,1,2,3]

    fixedShaft=True
    if 'S1' in fstFile:
        fixedShaft=False

    if 'B100' in fstFile:
        shapes_bld=[0]
        fixedShaft=False
    elif 'B110' in fstFile:
        shapes_bld=[0,1]
        fixedShaft=False
    elif 'B111' in fstFile:
        shapes_bld=[0,1,2]
        fixedShaft=False
    elif 'B1' in fstFile:
        shapes_bld=[0]
        fixedShaft=False

    if len(shapes_sub)==0:
        WT =  FASTmodel2TNSB(fstFile, shapes_twr=shapes_twr, shapes_bld=shapes_bld, DEBUG=False, 
                             bStiffening=True, main_axis='z', 
#                              algo='OpenFAST', 
                              algo='', 
                             fixedShaft=fixedShaft,
                             assembly='auto'
                             ).WT
    else:
        WT = FASTmodel2MTNSB(fstFile, shapes_sub=shapes_sub, shapes_twr=shapes_twr, shapes_bld=shapes_bld, DEBUG=False, 
                             bStiffening=True, main_axis='z', 
                             fixedShaft=fixedShaft, 
                             algo='', 
#                              algo='OpenFAST', 
                             SD_bOverride=True,
                             SD_bCI=True,
    #                          FEM_method='cbeam').WT
                             FEM_method='full').WT



    print('------------ beamFEM')
#     if WT.SD._FEM is None:
#         WT.SD.init(TP=(0,0,15))
        #WT.SD.toYAML(yamlFile)
    #FEM = WT.SD.beamFEM(method='full')
    #print('FEM keys',FEM.keys())

    # --- 
    print('---------------------- OPENFAST LIN FILE ------------------------')
    lin = FASTLinearizationFile(linFile)
    dfs=lin.toDataFrame()
    #print(dfs['A'].keys())
    #printMat('AA',dfs['A'].values)
    #printMat('AA',dfs['A'].values, digits=1, xmin=1e-8)
    freq_d, zeta, Q, freq_0  = lin.eva()
    for i, (fd, f0, ze) in enumerate(zip(freq_d, freq_0, zeta)):
        print('Mode {:4d} {:7.3f} {:7.3f} '.format(i, f0, ze))
    freq_0_OF = freq_0
    zeta_OF = zeta
    print('----------------------------------------------------------------')

    print('----------------------YAMS  -------')
    q = np.zeros((len(WT.q0),1))
    WT.grd.setDOF(q) # <<<<< NOTE: unless the q_init was a good initial value, for now we use all zero 
    freq_d, zeta, Q, freq_0 = WT.grd.eva()
    A,B,C,D = BuildSystem_Linear_MechOnly(WT.MM, WT.DD, WT.KK) #, nP=len(sAug), nU=len(sInp), nY=len(sMeas), Fp=None)
    #rintMat('AA', A, digits=1, xmin=1e-8)
    for i, (fd, f0, ze) in enumerate(zip(freq_d, freq_0, zeta)):
        #print('Mode {:4d} {:7.3f} {:7.3f} '.format(i, f0, ze))
        print('Mode {:4d} {:7.3f} {:7.3f}'.format(i, f0, ze))
    print('----------------------------------------------------------------')
    print(WT.q0.T)

    print('------------------------------------------------------------------')
    print(' Mode (#) |       Frequency (Hz)       |     Damping ratio (-)')
    print('------------------------------------------------------------------')
    print('          |     YAMS     OF      Err.  |    YAMS     OF      Err.')
    for i, (fd, f0, ze, f1, z1) in enumerate(zip(freq_d, freq_0, zeta, freq_0_OF, zeta_OF)):
        print('  {:4d}    |  {:7.3f} {:7.3f} {:7.1f}%  | {:7.3f} {:7.3f} {:7.1f}%'.format(i+1, f0, f1, (f1-f0)/f1*100, ze, z1, (z1-ze)/z1*100))
    print('------------------------------------------------------------------')



    # --- Subset
#     cols = [c for c in dfs['A'].keys() if 'SDCraig' not in c]
#     AA = dfs['A'].loc[cols, cols]
#     freq_d, zeta, Q, freq_0  = basicEVA(AA)
#     for i, (fd, f0, ze) in enumerate(zip(freq_d, freq_0, zeta)):
#         print('Mode {:4d} {:7.3f} {:7.3f} '.format(i, f0, ze))


    print('----------------------------------------------------------------')
#     print('----------------------FULL FEM  -------')
#     for i, f0 in enumerate(WT.SD._FEM.freq):
#         print('Mode {:4d} {:7.3f}'.format(i, WT.SD._FEM.freq[i]))
#         if i>5:
#             break
#     print('----------------------GUYAN FEM -------')
#     for i, f0 in enumerate(WT.SD._FEM.f_G):
#         print('Mode {:4d} {:7.3f}'.format(i, WT.SD._FEM.f_G[i]))
#         if i>5:
#             break

    print('----------------------------------------------------------------')
