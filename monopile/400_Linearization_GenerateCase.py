import numpy as np
import os
import welib.weio as weio
import welib.fast as fastlib
from welib.essentials import *


# --- Parameters for this script
ref_dir     = 'simulations/04_linTemplate/'
FAST_EXE    = 'simulations/of42.exe'
# FAST_EXE    = 'simulations/of421-EDMM.exe'
# trim_points = 'simulations/IEA-22-280-RWT/IEA-22-280-RWT_OperOpenFAST.csv'

main_file   = 'OF_F0T4S1B3.fst'; work_dir   = 'simulations/_04_LinF0T4S1B3/'; trim_points = '_data/PowerCurve_F0T4S1B3.csv'
# main_file   = 'OF_F2T4S1B3.fst'; work_dir   = 'simulations/_04_LinF2T4S1B3/'; trim_points = '_data/PowerCurve_F4T4S1B3.csv'
main_file   = 'OF_F0T4S1B3_SvDSimp.fst'; work_dir   = 'simulations/_04_LinF0T4S1B3_SvDSimp/' ; sSim = 'F0T4S1B3_SvDSimp'


# --- Defining the parametric study  (list of dictionnaries with keys as FAST parameters)
Cases=weio.read(trim_points).toDataFrame()
# Cases=Cases.iloc[:-6:1,:]
# Adding 0 wind
# Cases=Cases.append(Cases.iloc[0,:])
# Cases.iloc[-1,:] = 0
# Cases['Wind [m/s]'].loc[-1] = 0.0001
# Cases['Pitch [deg]'].loc[-1] = Cases['Pitch [deg]'].loc[0]

print(Cases)


Tmin       = 1500
nPerPeriod = 6

PARAMS     = []
for i,c in Cases.iterrows():
    ws    = c['Wind_Speed_[m/s]']
    rpm   = c['Rotor_Speed_[rpm]']
    pitch = c['Pitch_[deg]']
#     tt    = c['TTDspFA_[m]']
    Omega = rpm/60*2*np.pi
    if abs(Omega)<0.001:
        LinTimes = [Tmin]
        Tmax     = Tmin+1
    else:
        T = 2*np.pi/Omega
        LinTimes = np.linspace(Tmin,Tmin+T,nPerPeriod+1)[:-1];
        Tmax       = Tmin+1.01*T
    print('rpm',rpm,'Omega',Omega,'T',T)
        
    p=dict()
    p['__name__']       = 'ws_{:.0f}'.format(ws)
    p['FAST|TMax']         = Tmax
    p['FAST|TStart']       = 0
    p['FAST|Linearize']    = True
    p['FAST|NLinTimes']    = len(LinTimes)
    p['FAST|LinTimes']     = list(LinTimes)
    p['FAST|LinInputs']    = 2                   # <<<<<<<<<<<<<<<<<< To get full linearizations
    p['FAST|LinOutputs']   = 1                   # <<<<<<<<<<<<<<<<<< To get full linearizations
    p['FAST|OutFmt']       = '"ES18.11E2"'
#     p['FAST|CompServo']  = 0
    if abs(ws)<0.001:
        p['FAST|CompAero']    = 0
        p['FAST|CompInflow']  = 0
    # --- Trim
    p['FAST|CalcSteady'] = True
    if ws>=12:
        p['FAST|TrimCase']   = 3    # 3:Pitch
    else:
        p['FAST|TrimCase']   = 2    # 2:torque
    p['FAST|TrimTol']    = 0.001 # default 0.001
    p['FAST|TrimGain']   = 0.01 # default 0.01
    p['EDFile|GenDOF']   = True

    p['InflowFile|WindType'] = 1
    p['InflowFile|HWindSpeed'] = ws
    p['InflowFile|PLexp']  = 0.0                 # <<<<<<<<<<<<<<<<<
    p['EDFile|BlPitch(1)'] = pitch
    p['EDFile|BlPitch(2)'] = pitch
    p['EDFile|BlPitch(3)'] = pitch
    p['EDFile|RotSpeed']   = rpm
#     p['EDFile|TTDspFA']    = tt
    p['EDFile|FlapDOF1'] = False
    p['EDFile|FlapDOF2'] = False
    p['EDFile|EdgeDOF']  = False
    p['EDFile|TeetDOF']  = False
    p['EDFile|DrTrDOF']  = False
#     p['EDFile|GenDOF']   = False
    p['EDFile|YawDOF']   = False



    PARAMS.append(p)
# --- Generating all files in a workdir
fastfiles=fastlib.templateReplace(PARAMS,ref_dir, outputDir=work_dir, removeRefSubFiles=True, main_file=main_file)

# --- Creating a batch script just in case
fastlib.writeBatch(os.path.join(work_dir,'_Linearization_RUN_ALL.bat'),fastfiles,fastExe=FAST_EXE, nBatches=4)
# --- Running the simulations
# fastlib.run_fastfiles(fastfiles, fastExe=FAST_EXE, parallel=True, showOutputs=True, nCores=2)
# 
# 
