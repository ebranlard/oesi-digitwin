import numpy as np
import os
import matplotlib.pyplot as plt

from welib.fast.case_gen import CPCT_LambdaPitch
import welib.fast.postpro as postpro
#import pyFAST.case_geneneration.case_gen as case_gen
#import pyFAST.input_output.postpro as postpro

# Get current directory so this script can be called from any location

def CPLambdaRun():
    FAST_EXE   = 'simulations/of42.exe'               # Location of a FAST exe (and dll)
    ref_dir    = 'simulations/IEA-22-280-RWT/'        # Folder where the fast input files are located (will be copied)
    main_file  = 'OF.fst'                              # Main file in ref_dir, used as a template
    workDir    = 'simulations/_CPLambdaPitch/'
    exportBase = 'simulations/IEA-22-280-RWT/IEA-22-280-RWT_CaseGen_'
    TMax=10

    # --- Computing CP and CT matrices for range of lambda and pitches
    Lambda = np.linspace(2 ,14,20)
    Pitch  = np.linspace(-5,45,20)

    Lambda = np.linspace(2 ,14,3)
    Pitch  = np.linspace(-5,45,3)

    baseDict = {'CompServo':0, 'CompSeaSt':0, 'CompHydro':0, 'CompSub':0}

    CPCT_LambdaPitch(ref_dir, main_file, Lambda, Pitch, 
                     fastExe=FAST_EXE, showOutputs=False, nCores=5, workDir=workDir,
                     TMax=TMax, baseDict=baseDict, bStiff=True,
                     exportFmt = 'ROSCO', exportBase=exportBase, 
                     reRun = False, # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                     plot=True, verbose=True)



if __name__=='__main__':
    CPLambdaRun()
    plt.show()

