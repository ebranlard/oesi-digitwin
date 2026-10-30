"""

"""

import numpy as np
import matplotlib.pyplot as plt
from welib.essentials import *

# from welib.yams.FNSB_FAST import FASTmodel2FNSB
# from welib.yams.windturbine import FASTWindTurbine
from welib.yams.models.MTNSB import FASTmodel2MTNSB
# from welib.system.mech_system import MechSystem
from welib.hydro.wavekin import wavenumber, elevation2d
from welib.weio.fast_linearization_file import FASTLinearizationFile

setFigurePath('../article4_omae2026/figs/')
setFigureFont(13)

# ---
# fstFile = 'simulations_wt/00_EVA/OF_F3T1S1_H0A0.fst'
fstFile = 'simulations/00_EVA/OF_F2T1S0_H0A0_NoRot.fst'
linFile = fstFile.replace('.fst','.1.lin')
export=False
# export=True

# --- Loading a "YAMS" model
WT = FASTmodel2MTNSB(fstFile, shapes_sub=[0,4], shapes_twr=[0], shapes_bld=[], DEBUG=False, 
                     bStiffening=True, main_axis='z', fixedShaft=True, algo='OpenFAST').WT


# --------------------------------------------------------------------------------}
# --- Just the Guyan modes 
# --------------------------------------------------------------------------------{
# SD = WT.SD
# 
# FEM = SD.beamFEM()
# print(SD)
# 
# nModesPlot=4
# x=FEM['xNodes'][0,:]
# Q=FEM['Q']
# 
# # fig,axes = plt.subplots(1, 2, sharex=False, sharey=True, figsize=(6.4,3.9)) # (6.4,4.8)
# 
# fig = plt.figure(figsize=(6.4,3.8))
# # gs = fig.add_gridspec(1, 3, width_ratios=[2, 1, 1])
# gs = fig.add_gridspec(1, 4, width_ratios=[5, 1, 5, 5])
# ax1 = fig.add_subplot(gs[0, 0])
# ax2 = fig.add_subplot(gs[0, 2])
# ax3 = fig.add_subplot(gs[0, 3])
# fig.subplots_adjust(left=0.12, right=0.94, top=0.94, bottom=0.13, hspace=0.27, wspace=0.23)
# 
# 
# fs=12
# # ax.plot(x, Q[1::6,0]  ,'--'  , label=r'$u_y$')
# ax1.plot(Q[2::6,0], (x-np.min(x))/(np.max(x)-np.min(x)), '-', c=fColrs(1))
# ax1.set_title('FEM 1', fontsize=fs)
# 
# z= WT.fnd.s_span
# z /=np.max(z)
# ax2.plot(WT.fnd.PhiU[0][0], z, '-'  , c=fColrs(1),  label='Guyan - Surge')
# ax2.set_title('Guyan 1 (Surge)', fontsize=fs)
# ax3.plot(WT.fnd.PhiU[1][0], z, '--' , c=fColrs(2),  label='Guyan - Pitch')
# ax3.set_title('Guyan 2 (Pitch)', fontsize=fs)
# 
# 
# ax1.set_xlabel(r'$\Phi_1(z)$ [-]')
# ax2.set_xlabel(r'$\Phi_1(z)$ [-]')
# ax3.set_xlabel(r'$\Phi_2(z)$ [-]')
# for ax in [ax1, ax2, ax3]:
#     #ax.set_xlabel('Shape function [-]')
#     ax.set_ylim([0,1])
#     ax.tick_params(direction='in', top=True, right=True)
# 
# # ax2.set_yticklabels([])
# ax3.set_yticklabels([])
# 
# ax1.set_ylabel('Vertical coordinate [-]')
# # plt.tight_layout()
# 
# fig._title = 'StructuralShapeFunctions' 
# plt.show()


# --------------------------------------------------------------------------------}
# --- TRY EVA and Shape function Concat
# --------------------------------------------------------------------------------{


# --- 
print('----- OPENFAST LIN FILE ----')
lin = FASTLinearizationFile(linFile)
dfs=lin.toDataFrame()
print('keys:',dfs.keys())
printMat('AA',dfs['A'].values)
freq_d, zeta, Q, freq_0  = lin.eva()
for i, (fd, f0, ze) in enumerate(zip(freq_d, freq_0, zeta)):
    print('Mode {:4d} {:7.3f} {:7.3f} '.format(i, f0, ze))
print('-------------')




# --------------------------------------------------------------------------------}
# --- TRY EVA and Shape function Concat
# --------------------------------------------------------------------------------{
# --- Eigenvalue analysis 
norm='tip'
# print(WT.q0)
nDOF = len(WT.q0)
q = WT.q0.values*0
# print(WT.q0.values*0)
q = np.zeros((nDOF,1))


FEM = WT.SD.beamFEM()

WT.grd.setDOF(q) # <<<<< NOTE: unless the q_init was a good initial value, for now we use all zero 


freq_d, zeta, Q, freq_0 = WT.grd.eva()

print('-----YAMS  -------')

printMat('MM',WT.MM)
printMat('KK',WT.KK)
printMat('DD',WT.DD)

from welib.kalman.kalman import BuildSystem_Linear_MechOnly 
# A,B,C,D = BuildSystem_Linear_MechOnly(WT.MM, WT.DD, WT.KK, nP=len(sAug), nU=len(sInp), nY=len(sMeas), Fp=None)
A,B,C,D = BuildSystem_Linear_MechOnly(WT.MM, WT.DD, WT.KK) #, nP=len(sAug), nU=len(sInp), nY=len(sMeas), Fp=None)
printMat('AA', A)




print('-----YAMS  -------')
for i, (fd, f0, ze) in enumerate(zip(freq_d, freq_0, zeta)):
    print('Mode {:4d} {:7.3f} {:7.3f}'.format(i, f0, ze))
print('-------------')






nMaxModes = min(4, len(freq_0))
# 
# --- Modes
nModesPlot=3
modes = WT.grd.modes(norm=norm, main_axis='z', nBodies=3, Q_in=Q)
nMaxModes = min(nModesPlot, len(modes))

fig,axes = plt.subplots(1,nMaxModes, sharey=False, figsize=(12.4,4.8)) # (6.4,4.8)
fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
print(len( modes[0][0,1:]))
for im, (m,ax) in enumerate(zip(modes,axes)):
#     ax.plot(U_th[im,:], x_th, 'k--',label='Theory')
    # NOTE: first point is "reference/ground" frame (0,0,0)
    #ax.plot(m[1,1:]/m[1,-1], m[0,1:], label='y')
    ax.plot(m[0,1:]/np.max(m[0,1:]), m[2,1:], label='z')
    #ax.plot(m[0,1:], m[2,1:], label='z')
#     ax.set_aspect('equal', 'box')
    ax.set_title('Mode {}'.format(im+1))
ax.set_xlabel('')
ax.set_ylabel('')
ax.legend()
# 
# 










if export:
    export2pdf()


plt.show()
