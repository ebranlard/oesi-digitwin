#
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
# Local
from welib.essentials import *
import welib.weio as weio

from welib.yams.windturbine import FASTWindTurbine
from welib.yams.section_loads_WT import YAMSSectionLoadCalculator

# --------------------------------------------------------------------------------}
# --- Script parameters
# --------------------------------------------------------------------------------{
tRange=[0,600]
compFile = None
shapes_sub = [0,4]
dt_resample=0.2
topLoadMethod='OFAero'

tRange = [0,10]
tRange = None
# topLoadMethod='OF'

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
# fstFile = 'simulations/05_RegWave/OF.fst'; tRangeStats=[40,100];
# fstFile = 'simulations/05_RegWave/OF_F2_NoRNA.fst'; tRangeStats=[40,100];
# fstFile = 'simulations/05_RegWave/OF_F3_NoRNA.fst'; tRangeStats=[40,100];
# fstFile = 'simulations/05_RegWave/OF_F3_RNA.fst'; tRangeStats=[40,100];
# fstFile = 'simulations/05_RegWave/OF_F3T0_NoRNA.fst'; tRangeStats=[40,100];

# --- Full turbine
fstFile = 'simulations/06_Jonswap/OF_F3T1S1_H1A1_Hs=8.1_Tp=12.7.fst'; tRangeStats=[40,100];
compFile = 'simulations/Waves/UserDefJonswap_Hs=8.1_Tp=12.7_h=34.csv'

# --------------------------------------------------------------------------------}
# --- Main script 
# --------------------------------------------------------------------------------{
# --- Script derived parameters
Case = os.path.basename(fstFile.replace('.fst',''))
outFigDir = os.path.join(os.path.dirname(fstFile),  'figs/')
outFile   = os.path.join(outFigDir, os.path.basename(fstFile).replace('.fst',''))
os.makedirs(outFigDir, exist_ok=True)
INFO(f'Case: {os.path.dirname(fstFile)}/{Case}')

df = weio.read(fstFile.replace('.fst','.outb')).toDataFrame()

WT = FASTWindTurbine(fstFile, algo='OpenFAST', SD_FEM_method='cbeam', subShapes=shapes_sub, verbose=True, fixedShaft=False).WT
if WT.pSS is not None:
    NOTE('Setting Components', compFile)
    WT.SS_setComponents(compFile)
#     NOTE('Setting Compute Eta')
#     WT.SS_computeEta(dfRef['Time_[s]'])
#     if hydroShapeFile is not None:
#         WT.HD_setShapeFunction(hydroShapeFile)
else:
    F_secRef = None


YSL = YAMSSectionLoadCalculator(WT=WT)

if topLoadMethod=='OFAero':
    df['Fadd_R_xh'] = df['RtAeroFxh_[N]']
    df['Fadd_R_yh'] = df['RtAeroFyh_[N]'] #*0
    df['Fadd_R_zh'] = df['RtAeroFzh_[N]'] #*0
    df['Madd_R_xh'] = df['RtAeroMxh_[N-m]']
    df['Madd_R_yh'] = df['RtAeroMyh_[N-m]']   # TODO TODO Quite Important
    df['Madd_R_zh'] = df['RtAeroMzh_[N-m]']   # TODO TODO Quite Important
    useTopLoadsFromDF=False
elif topLoadMethod=='OFAero_s':
    df['Fadd_R_xs'] = df['RtAeroFxh_[N]']
    df['Fadd_R_ys'] = df['RtAeroFyh_[N]']
    df['Fadd_R_zs'] = df['RtAeroFzh_[N]']
    df['Madd_R_xs'] = df['RtAeroMxh_[N-m]']
    df['Madd_R_ys'] = df['RtAeroMyh_[N-m]']
    df['Madd_R_zs'] = df['RtAeroMzh_[N-m]']
    useTopLoadsFromDF=False
elif topLoadMethod=='OF':
    useTopLoadsFromDF=True
else:
    raise NotImplementedError()

dfOut, sections = YSL.fromDF(df, useTopLoadsFromDF=useTopLoadsFromDF, useInterfaceLoadsFromDF=False, dt_resample=dt_resample, tRange=tRange)
dfOut.export(fstFile.replace('.fst', '.YAMS.outb'))

# --- Thinkgs in
# fig = YSL.plot_comp(sig='Wave1Elev_[m]', ylabel='Wave elevation [m]', figFilename=outFile+'_SL_MNP_Eta.png', tRange=tRange)
# fig.axes[0].plot(WT.pSS['eta_time'], WT.pSS['eta'], ':', label='From SS')
# fig.savefig(outFile + '_SL_MNP_Eta.png')

# fig = YSL.plot_comp(sig='HydroFxi_[N]', ylabel='Hydro Fx [N]', figFilename=outFile+'_SL_MNP_Eta.png', tRange=tRange)


fig, axes = plt.subplots(6, 1, sharey=False, figsize=(6.4,9.8))
fig.subplots_adjust(left=0.16, right=0.98, top=0.99, bottom=0.06, hspace=0.30, wspace=0.20)
YSL.plot_comp(sig = 'YawBrFxp_[kN]', ax   = axes[0])
YSL.plot_comp(sig = 'YawBrFyp_[kN]', ax   = axes[1])
YSL.plot_comp(sig = 'YawBrFzp_[kN]', ax   = axes[2])
YSL.plot_comp(sig = 'YawBrMxp_[kN-m]', ax = axes[3])
YSL.plot_comp(sig = 'YawBrMyp_[kN-m]', ax = axes[4])
YSL.plot_comp(sig = 'YawBrMzp_[kN-m]', ax = axes[5])

# YSL.plot_tower_accelerations(IsecTwr=[0,4,8], component='x', figFilename=outFile+'_SL_TWR_Acc.png')
# 
# zDepth = YSL.sec['monopile']['z']
# IZ = [int(3*len(zDepth)/6)-5, int(2*len(zDepth)/6), int(1*len(zDepth)/6), 0]
# YSL.plot_monopile_section_loads_stats(IZ=IZ, component=0, figFilename=outFile+'_SL_MNP_SL_STATS.png', tRange=tRange, stat='meanabs')
# YSL.plot_monopile_section_loads      (IZ=IZ, component=0, figFilename=outFile+'_SL_MNP_SL.png', tRange=tRange)
# 
# fig = YSL.plot_comp(sig='M1N1MKye_[N*m]', ylabel='Sea bed moment [MNm]', figFilename=outFile+'_SL_MNP_M_SeaBed.png', tRange=tRange, scale=1e6)
# 
# 
YSL.plot_tower_section_loads(IsecTwr=[0,4,8], component='MLyt', figFilename=outFile + '_SL_TWR_SL.png', printStats=True)
YSL.plot_section_loads_stats(IZ=None, component=0, figFilename=None, tRange=tRange, stat='meanabs')



plt.show()
# 

