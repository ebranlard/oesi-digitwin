import numpy as np
import os
import welib.weio as weio
import welib.fast as fastlib
from welib.essentials import *


# --- Parameters for this script
ref_dir     = 'simulations/03_aeroTemplate/'
trim_points = 'simulations/IEA-22-280-RWT/IEA-22-280-RWT_OperOpenFAST.csv'
FAST_EXE    = 'simulations/of42.exe'
# FAST_EXE    = 'simulations/of421-EDMM.exe'

main_file   = 'OF_F0T4S1B3.fst'; work_dir   = 'simulations/_03_PC_F0T4S1B3/' ; sSim = 'F0T4S1B3'
# main_file   = 'OF_F4T4S1B3.fst'; work_dir   = 'simulations/_03_PC_F4T4S1B3/' ; sSim = 'F4T4S1B3'
main_file   = 'OF_F0T4S1B3_SvDSimp.fst'; work_dir   = 'simulations/_03_PC_F0T4S1B3_SvDSimp/' ; sSim = 'F0T4S1B3_SvDSimp'


Cases = weio.read(trim_points).toDataFrame()
print(Cases)
print('Keys:',Cases.keys())

# --- Defining the parametric study  (list of dictionnaries with keys as FAST parameters)
WS    = Cases['Wind_Speed_[m/s]'].values #[1:-12:2]
RPM   = Cases['Rotor_Speed_[rpm]'].values#[1:-12:2]
PITCH = Cases['Pitch_[deg]'].values#[1:-12:2]
BaseDict = {'TMax': 200}
PARAMS=[]
for i,(wsp,rpm,pitch) in enumerate(zip(WS, RPM, PITCH)): # NOTE: same length of WS and RPM otherwise do multiple for loops
    p=BaseDict.copy()
    p['OutFileFmt']  = 2
    p['FAST|Linearize']        = False
    p['FAST|CalcSteady']       = False
    p['EDFile|RotSpeed']       = rpm
    p['EDFile|BlPitch(1)']     = pitch
    p['EDFile|BlPitch(2)']     = pitch
    p['EDFile|BlPitch(3)']     = pitch
    p['EDFile|BlPitch(3)']     = pitch
    p['InflowFile|HWindSpeed'] = wsp
    p['InflowFile|WindType']   = 1 # Setting steady wind
    p['__name__']='ws{:02d}'.format(int(p['InflowFile|HWindSpeed']))
    PARAMS.append(p)
    i=i+1
# --- Generating all files in a workdir
fastfiles=fastlib.templateReplace(PARAMS, ref_dir, outputDir=work_dir, removeRefSubFiles=True, main_file=main_file, oneSimPerDir=False)
print(f'FST file {0:5d}: {fastfiles[0]}')
print(f'FST file {len(fastfiles):5d}: {fastfiles[-1]}')

# --- Creating a batch script just in case
fastlib.writeBatch(os.path.join(work_dir,'_RUN_ALL.bat'),fastfiles,fastExe=FAST_EXE, nBatches=4, echo=False)
# --- Running the simulations
fastlib.run_fastfiles(fastfiles, fastExe=FAST_EXE, parallel=True, showOutputs=False, nCores=4)

# --- Simple Postprocessing
outFiles = [os.path.splitext(f)[0]+'.outb' for f in fastfiles]
avg_results = fastlib.averagePostPro(outFiles,avgMethod='periods',avgParam=1, ColMap = {'Wind_Speed_[m/s]':'Wind1VelX_[m/s]','Pitch_[deg]':'BldPitch1_[deg]','Rotor_Speed_[rpm]':'RotSpeed_[rpm]','Mech_Power_[kW]':'RotPwr_[kW]'},ColSort='Wind_Speed_[m/s]')
print('>> Average results:')
print(avg_results)

avg_results.to_csv('PowerCurve_{}.csv'.format(sSim))
