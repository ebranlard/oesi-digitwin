import numpy as np
import pandas as pd
import os
# Local
import welib.weio as weio
from welib.fast.FASTLin import *
from welib.essentials import *

# --- Script parameters
bLoadAll=True # <<<< Create the dump file
# bLoadAll=False 
bExport=True
Method='mean'; WS=None
Method='ws'; WS=4 # For a given WS


main_file   = 'OF_F0T4S1B3.fst'; work_dir   = 'simulations/_04_LinF0T4S1B3/'
# main_file   = 'OF_F2T4S1B3.fst'; work_dir   = 'simulations/_04_LinF2T4S1B3/'


main_file   = 'OF_F0T4S1B3.fst'; work_dir   = 'C:/Work/OuiWind/_04_LinF0T4S1B3_SvDSimp-EDMM2/'

sSim = work_dir.split('_04_Lin')[1]
nDOF=2
dumpFile  = f'_dump/FASTLin_dump_{sSim}.dat'
LinFolder = work_dir
OutFile   = f'data/FASTLin_ABCD_{sSim}_Method={Method}.dat'

# --- Read all the linearizations from FAST
if bLoadAll:
    FL=FASTLin(folder=LinFolder, prefix='ws_', verbose=True)
    FL.save(filename=dumpFile, verbose=True)
else:
    FL = FASTLin.from_pickle(dumpFile, verbose=True)

# import pdb; pdb.set_trace()
# FL= FASTLin(folder=LinFolder, prefix=f'ws_{WS}', nLin=1)

# --- Mean matrices
if Method=='mean':
    kwargs = {}
elif Method=='ws':
    kwargs = {'WS':[WS]}
else:
    raise NotImplementedError()

# with TimeIt('Computing stats...'):
#     A,B,C,D = FL.average(**args)
# try:
#     M = FL.stats('M')[0]
# except:
#     FAIL('M is missing')
# 
# print(FL.OP_Data[0].Data[0])
# print('A0\n',FL.OP_Data[0].Data[0].A)
# print('Amean\n',A)
# A=FL.OP_Data[0].Data[0]['A']
# B=FL.OP_Data[0].Data[0]['B']
# C=FL.OP_Data[0].Data[0]['C']
# D=FL.OP_Data[0].Data[0]['D']


sX_DOF  = ['qt1FA_[m]', 'psi_rot_[rad]', 'd_qt1FA_[m/s]', 'd_psi_rot_[rad/s]']
sED_DOF = ['7_TwFADOF1' ,'13_GeAz']
sU_DOF  = ['HubFxN1_[N]','Qgen_[Nm]','PitchColl_[rad]']
# sY_DOF  = ['NcIMUTAxs_[m/s^2]', 'RotSpeed_[rpm]','SvDGenTq_[kNm]', 'BPitch1_[deg]']
sY_DOF  = ['NcIMUTAxs_[m/s^2]', 'RotSpeed_[rpm]','SvDGenTq_[kNm]', 'BPitch1_[deg]']

Ar, Br, Cr, Dr, Mr = FL.average_subset(sX_sel=sX_DOF, sU_sel=sU_DOF, sY_sel=sY_DOF, sE_sel=sED_DOF, **kwargs)


# import pdb; pdb.set_trace()
#     def average_subset(self, sX_sel=None, sU_sel=None, sY_sel=None, sE_sel=None, WS=None, exportFile=None, baseDict=None):

# Ar, Br, Cr, Dr, Mr = extractABCD_forModel(A, B, C, D, model='T1S1', M=None)
# print(Ar)

