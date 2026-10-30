""" 
README: 
    This file can be run twice.

    First, without reconHydro and hydroShape, to export the pickle
         then we need to run 101_*.py to generate the shape function.

    Afterwards, we can rerun with reconHydro and the proper hydroShapec




Conclusions from Analyses of MT100 (OLD)
 - SIM doesnt matter too much
 - CB0 is bad because SubDyn doesn't account for "distributed inertial loads"
   It's SubDyn's fault, YAMS does better
 - SubDyn has more noise as you go up, doesn't not give zero shear at the top..

"""
#
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
# Local
from welib.essentials import *
import welib.weio as weio

from helper_functions_digimono import sec_plotFM, sec_plot, time_plot, mean_rel_err
from monopileWT import monopileOutputs, monopileSetupFromOpenFAST
from welib.hydro.wavekin import elevation2d


def prescribedRun(fstFile, compFile=None, tRange=[0,600], shapes_sub=[0,4], nUnderSamp=1, reconHydro=False, tRangeStats=None):
    #reconHydro = False # <<<<<<<<<<<<<<<<<<<<<<<<<<< FOR FIRST RUN # TODO make a dedicated script
    # reconHydro = True 


    hydroShape = None # <<<<<<<<<<<<<<<<<<<<< FOR FIRST RUN
    # hydroShape = f'_data/IEAMonoPile_HydroShapeFunction{SSCase}.csv' # NOTE: optional

    # --------------------------------------------------------------------------------}
    # --- Main script 
    # --------------------------------------------------------------------------------{
    # --- Script derived parameters
    Case = os.path.basename(fstFile.replace('.fst',''))
    INFO(f'Case: {os.path.dirname(fstFile)}/{Case} - ReconHydro: {reconHydro}')

    # --- Setting up YAMS model
    pST, pSS, pHD, Sys, WT, ref = monopileSetupFromOpenFAST(fstFile, shapes_sub=shapes_sub, TMIN=tRange[0], TMAX=tRange[1], nSubSample=nUnderSamp, compFile=compFile, tuneM=False, hydroShape=hydroShape, reconHydro=reconHydro)
    zDepth      = pHD['zDepth']
    F_secRef    = ref['F_sec']
    state_ref   = ref['state']
    state_d_ref = ref['state_d']
    vTime       = ref['df']['Time_[s]'].values
    # --------------------------------------------------------------------------------}
    # --- Local Simulation
    # --------------------------------------------------------------------------------{
    with TimeIt('Calc outputs...'):
        F_sec  = np.zeros((6,len(zDepth),len(vTime)))
        F_hydro  = np.zeros(len(vTime))
        p_hydro  = np.zeros((len(zDepth), len(vTime)))
        p_hnoac  = np.zeros((len(zDepth), len(vTime)))
        p_hyref  = np.zeros((len(zDepth), len(vTime)))
        p_recon  = np.zeros((len(zDepth), len(vTime)))
        eta_sim = pSS['eta']
        eta_dot = pSS['eta_dot']
        if 'phi' in pHD:
            p_recon = pHD['phi'][:, np.newaxis] * eta_dot[np.newaxis, :] # Reconstructed p_hydro based on shape function
        else:
            print('>>>>>>>>>>>>>>>> PHI NOT PRESENT, EXPORT CASE?')

        for it,t in enumerate(vTime):
            Fs, outD   = monopileOutputs(t, state_ref[:,it], pST=pST, pSS=pSS, pHD=pHD, it=it, xd=state_d_ref[:,it], reconHydro=reconHydro)
            Fs_, outD_ = monopileOutputs(t, state_ref[:,it], pST=pST, pSS=pSS, pHD=pHD, it=it, xd=state_d_ref[:,it], reconHydro=False)
            F_sec[:,:,it] = Fs[:,:]
            F_hydro[it]   = outD['F_hydro']
            p_hydro[:, it] = outD['p_hydro']
            p_hnoac[:, it] = outD['p_hydro_no_struct_acc']
            p_hyref[:, it] = outD_['p_hydro']

    # --------------------------------------------------------------------------------}
    # --- Postpro 
    # --------------------------------------------------------------------------------{
    if not reconHydro:
        print('>>>>>>>>>>>>>>>> EXPORT')
        import dill
        results = {
            'pST': pST, 'pSS': pSS, 'pHD': pHD, 'Sys': Sys, 'WT': WT, 'ref': ref,
            'zDepth': zDepth, 'vTime': vTime, 'F_sec': F_sec,
            'F_hydro': F_hydro, 'p_hydro': p_hydro, 'p_hnoac': p_hydro, 'eta_sim': eta_sim
        }
        postproFile=fstFile.replace('.fst','_prescribed.dpkl')
        with open(postproFile, 'wb') as f:
            dill.dump(results, f)
        print('Written: ',postproFile)

    # --------------------------------------------------------------------------------}
    # --- Plots 
    # --------------------------------------------------------------------------------{
    # --- Section loads Plot
    fig = sec_plotFM(vTime, zDepth, F_sec, F_secRef, stat='meanabs', tRange=tRangeStats)

    fig = sec_plot(vTime, zDepth, p_hydro/1000, ref=p_hnoac/1000, other=p_recon/1000, stat='meanabs', tRange=tRangeStats, label='p_hydro', otherLab='Reconstructed')


    # IZ = [len(zDepth)-1, int(3*len(zDepth)/4), int(len(zDepth)/2), int(len(zDepth)/4), 0]
    IZ = [len(zDepth)-1, int(3*len(zDepth)/6), int(2*len(zDepth)/6), int(1*len(zDepth)/6), 0]

    fig, axes = plt.subplots(len(IZ), 1, sharey=True, sharex=True, figsize=(8.4,7.5))
    fig.subplots_adjust(left=0.12, right=0.96, top=0.95, bottom=0.05, hspace=0.07, wspace=0.20)
    for ii, iz in enumerate(IZ):
        time_plot (vTime, F_secRef[0, iz, :]/1e6,  F_sec[0, iz, :]/1e6, f'z={zDepth[iz]:.0f}m', tRange=tRangeStats, ax=axes[ii], fig=fig)
        print(f'z {zDepth[iz]:5.0f}: ', end=''); mean_rel_err(vTime, F_secRef[0, iz, :], vTime, F_sec[0, iz, :], tRange=tRangeStats)
    axes[-1].set_xlabel('Time [s]')
    fig.suptitle('FxSec [MN]')


    IZ = [np.abs(zDepth - val).argmin() for val in np.linspace(0, zDepth[0], 5)]

    fig, axes = plt.subplots(len(IZ), 1, sharey=True, sharex=True, figsize=(8.4,8.0))
    fig.subplots_adjust(left=0.12, right=0.95, top=0.96, bottom=0.05, hspace=0.20, wspace=0.01)
    for ii, iz in enumerate(IZ):
        time_plot (vTime, ref=p_hyref[iz, :]/1000,  sim=p_hydro[iz, :]/1000, other=p_recon[iz, :]/1000, label=f'z={zDepth[iz]:.0f}m', tRange=tRangeStats, ax=axes[ii], fig=fig, refLab='Not Rec.', otherLab='Reconstructed')
        if ii==0:
            axes[ii].legend()
    #     print('z', zDepth[iz])
    #     mean_rel_err(vTime, F_secRef[0, iz, :], vTime, F_sec[0, iz, :], tRange=tRangeStats)
    axes[-1].set_xlabel('Time [s]')
    fig.suptitle('pHydro [kN]')



    fig, ax = time_plot (vTime, ref['eta'], eta_sim, 'Wave elevation [m]', tRange=tRangeStats)

    fig, ax = time_plot (vTime, ref['HydroFx'], F_hydro, 'Hydro Fx [N]', tRange=tRangeStats)

    fig, ax = time_plot (vTime, F_secRef[4,0,:]/1e6, F_sec[4,0,:]/1e6, label='Sea bed moment [MNm]', tRange=tRangeStats)


    print('Eta     :', end=''); mean_rel_err(vTime, ref['eta'], vTime, eta_sim, tRange=tRangeStats)
    print('M_bot   :', end=''); mean_rel_err(vTime, F_secRef[4,0,:], vTime, F_sec[4,0,:], tRange=tRangeStats)

    plt.show()


if __name__ == '__main__':

    # --------------------------------------------------------------------------------}
    # --- Script parameters
    # --------------------------------------------------------------------------------{

    # --- Monopile Decay no Hydro - NOTE: Best is CB6
    # fstFile     = 'simulations/MT100/01_Decays/OF_NoHydro.fst'
    # fstFile     = 'simulations/MT100/01_Decays/OF_NoHydro_CB0_NoSim.fst'
    # fstFile     = 'simulations/MT100/01_Decays/OF_NoHydro_CB0.fst'
    # fstFile     = 'simulations/MT100/01_Decays/OF_NoHydro_CB2_NoSim.fst'
    # fstFile     = 'simulations/MT100/01_Decays/OF_NoHydro_CB2.fst'
    # fstFile     = 'simulations/MT100/01_Decays/OF_NoHydro_CB6_NoSim.fst'
    # fstFile     = 'simulations/MT100/01_Decays/OF_NoHydro_CB6.fst'
    # tRangeStats=[10.2,14.3];

    # --- Monopile Reg Wave
    #fstFile = 'simulations/05_RegWave/OF.fst'; tRangeStats=[40,100];
    #fstFile = 'simulations/05_RegWave/OF_F2_NoRNA.fst'; tRangeStats=[40,100];
    #fstFile = 'simulations/05_RegWave/OF_F3_NoRNA.fst'; tRangeStats=[40,100];
    #fstFile = 'simulations/05_RegWave/OF_F3_RNA.fst'; tRangeStats=[40,100];
    #fstFile = 'simulations/05_RegWave/OF_F3T0_NoRNA.fst'; tRangeStats=[40,100];

    # --- Monopile Reg Wave No Drag
    # fstFile = '../simulations/MT100/05_RegWave/OF_NoDrag.fst';
    # tRangeStats=[40,100];

    # # --- Monopile Reg Wave CB0
    # fstFile = '../simulations/MT100/05_RegWave/OF_CB0.fst';
    # flxFile = '../code3_fixedbottom/Monopile/_Flex_Monopile/res/Flex_MT100_RegularWave_2Shapes.res';
    # tRangeStats=[40,100];

    # --- Monopile Reg Wave CrazyCA
    # fstFile = '../simulations/MT100/05_RegWave/OF_CrazyCA.fst';
    # tRangeStats=[40,100];


    # --- MONOPILE JONSWAP LONG - _HS=2.5_TP=10.
    # SDCase='F3T1S1_H1A0'; SSCase='_Hs=2.5_Tp=10'; fstFile = f'simulations/06_Jonswap/OF_{SDCase}{SSCase}.fst'; tRange=[0,600]; tRangeStats=[35,600];
    # compFile =f'simulations/Waves/UserDefJonswap{SSCase}_h=34.csv'


    # --- MONOPILE JONSWAP LONG - _HS=8.1_TP=12.7
    # SDCase='F2T1S0_H1A0'; SSCase='_Hs=8.1_Tp=12.7'; fstFile = f'simulations/06_Jonswap/OF_{SDCase}{SSCase}.fst'; tRange=[0,600]; tRangeStats=[35,600];
    # compFile =f'simulations/Waves/UserDefJonswap{SSCase}_h=34.csv'

    # SDCase='F2T1S1_H1A0'; SSCase='_Hs=8.1_Tp=12.7'; fstFile = f'simulations/06_Jonswap/OF_{SDCase}{SSCase}.fst'; tRange=[0,600]; tRangeStats=[35,600];
    # compFile =f'simulations/Waves/UserDefJonswap{SSCase}_h=34.csv'

    # --- Works well
    #SDCase='F3T0_NoRNA_H1A0'; SSCase='_Hs=8.1_Tp=12.7'; fstFile = f'simulations/06_Jonswap/OF_{SDCase}{SSCase}.fst'; tRange=[0,600]; tRangeStats=[35,600];
    #compFile =f'simulations/Waves/UserDefJonswap{SSCase}_h=34.csv'
    fstFile = f'simulations/06_Jonswap/OF_F3T0_NoRNA_H1A0_Hs=8.1_Tp=12.7.fst'; tRange=[0,600]; tRangeStats=[35,600];
    compFile =f'simulations/Waves/UserDefJonswap_Hs=8.1_Tp=12.7_h=34.csv'
    compFile =None

    # SDCase='F3T1S1_H1A0'; SSCase='_Hs=8.1_Tp=12.7'; fstFile = f'simulations/06_Jonswap/OF_{SDCase}{SSCase}.fst'; tRange=[0,600]; tRangeStats=[35,600];
    # compFile =f'simulations/Waves/UserDefJonswap{SSCase}_h=34.csv'






    # --- Monopile Jonswap 
    # fstFile     = '../simulations/MT100/06_Jonswap/OF.fst';
    # flxFile = '../code3_fixedbottom/Monopile/_Flex_Monopile/res/Flex_MT100_Jonswap_2Shapes.res';
    # compFile   = '../code3_fixedbottom/Monopile/Waves/UserDefSimpleJonswap50.csv'
    # tRangeStats=[35,100];

    # --- Monopile Jonswap No Drag
    # fstFile     = '../simulations/MT100/06_Jonswap/OF_NoDrag.fst';
    # compFile   = '../code3_fixedbottom/Monopile/Waves/UserDefSimpleJonswap50.csv'
    # tRangeStats=[35,100]; 

    # # --- Monopile Jonswap CB0
    # fstFile     = '../simulations/MT100/06_Jonswap/OF_CB0.fst';
    # compFile   = '../code3_fixedbottom/Monopile/Waves/UserDefSimpleJonswap50.csv'
    # tRangeStats=[35,100]; 

    # # --- Monopile Jonswap CrazyCA
    # fstFile     = '../simulations/MT100/06_Jonswap/OF_CrazyCA.fst';
    # compFile   = '../code3_fixedbottom/Monopile/Waves/UserDefSimpleJonswap50.csv'
    # tRangeStats=[35,100]; 


    prescribedRun(fstFile, compFile=compFile, tRange=tRange, nUnderSamp=30, tRangeStats=tRangeStats)
    # 

