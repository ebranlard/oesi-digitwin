import numpy as np
import os
import glob
import matplotlib.pyplot as plt
import welib.fast.fastlib as fastlib
import welib.weio as weio
from welib.essentials import *


main_file   = 'OF_F0T4S1B3.fst'; work_dir   = 'simulations/_03_PC_F0T4S1B3/' ; sSim = 'F0T4S1B3'
main_file   = 'OF_F4T4S1B3.fst'; work_dir   = 'simulations/_03_PC_F4T4S1B3/' ; sSim = 'F4T4S1B3'
main_file   = 'OF_F0T4S1B3_SvDSimp.fst'; work_dir   = 'simulations/_03_PC_F0T4S1B3_SvDSimp/' ; sSim = 'F0T4S1B3_SvDSimp'

# --- Parameters for this script
fastfiles = glob.glob(os.path.join(work_dir, 'ws*.fst'))
print(f'FST file #{1:5d}: {fastfiles[0]}')
print(f'FST file #{len(fastfiles):5d}: {fastfiles[-1]}')


# --- Simple Postprocessing

# "RtAeroCp" "RtAeroCq" "RtAeroCt" "RtAeroPwr" "RtArea" "RtSkew" "RtSpeed" "RtTSR" "RtVAvgxh" "RtVAvgyh" "RtVAvgzh"
# ColMap = {'WS_[m/s]':'Wind1VelX_[m/s]','Pitch_[deg]':'BldPitch1_[deg]','RPM_[rpm]':'RotSpeed_[rpm]','P_[kW]':'RotPwr_[kW]'}
ColMap = {'Wind_Speed_[m/s]':'Wind1VelX_[m/s]','Pitch_[deg]':'BldPitch1_[deg]','Rotor_Speed_[rpm]':'RotSpeed_[rpm]','Mech_Power_[kW]':'RotPwr_[kW]'}
ColMap['Torque_[kNm]']           = 'RotTorq_[kN-m]'
ColMap['Thrust_[N]']            = '{RotThrust_[kN]}*1000'
ColMap['Generator_Torque_[kNm]'] = 'GenTq_[kN-m]'
ColMap['Elec_Power_[kW]']        = 'GenPwr_[kW]'


outFiles = [os.path.splitext(f)[0]+'.outb' for f in fastfiles]
avg_results = fastlib.averagePostPro(outFiles, avgMethod='periods', avgParam=2, ColMap=ColMap, ColSort='Wind_Speed_[m/s]')
avg_results.drop('Time_[s]', axis=1, inplace=True)
print('>> Average results:')
print(avg_results)

avg_results.to_csv('_data/PowerCurve_{}.csv'.format(sSim), index=False)
