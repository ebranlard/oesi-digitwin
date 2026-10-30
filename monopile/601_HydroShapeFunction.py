"""
NOTE: the shape functions can be put into two groups:

- Group1: the Peak ratio, Least Square and Analytical Sum, are similar
    These are Linear Estimates.
    They effectively represent the "Average" or "Expected" shape. 
    Because the sum is linear, the high-frequency components (which exist at the surface
    but are zero at the bed) "dilute" the average, pulling the seabed value down toward 0.15.

    >>> Use it for time-domain reconstruction (Expansion)


- Group2: The RMS Ratio, and Analytical RMS are almost identical. 
    These are Energy-Based Estimates. 
    They use the square of the amplitudes ($a_i^2$), they weight the "Energy" of the sea state. 
    Since energy is concentrated at the peak period ($T_p$), and long waves carry their energy 
    deep into the water column, these methods give more weight to the waves that actually 
    reach the seabed, resulting in the higher 0.23 value.

    >>> Use it for bulk loading


"""
import re
import numpy as np
import matplotlib.pyplot as plt
from welib.essentials import *
from welib.hydro.hydro_shape_function import hydro_shape_function, hydro_shape_function_from_TS, verify_hydro_reconstruction
from welib.hydro.hydro_shape_function import hydro_shape_from_pickle

if __name__ == '__main__':


    # --------------------------------------------------------------------------------}
    # --- Params 
    # --------------------------------------------------------------------------------{
    setFigurePath('../article6_journal/')
    setFigureFont(13)
    # export =True
    export =False

    IPlot=[]
#     IPlot=[0 ] # DEBUG
    IPlot+=[1 ]  # Article plot
    DEPTH = 34
    tRange=[0,600]; tRangeStats=[35,600]; 
    tRange=[35,600]; 

    # --- _HS=2.5_TP=10
    SSCase=f'_Hs=2.5_Tp=10';   fstFile = 'simulations/06_Jonswap/OF_F2T1S0_H1A0_Hs=2.5_Tp=10.fst'; 
    # --- _HS=8.1_TP=12.7
    SDCase='F3T0_NoRNA_H1A0'; SSCase='_Hs=8.1_Tp=12.7'; fstFile = f'simulations/06_Jonswap/OF_{SDCase}{SSCase}.fst';
    # SSCase= '_Hs=8.1_Tp=12.7'; fstFile = 'simulations/06_Jonswap/OF_F2T1S0_H1A0_Hs=8.1_Tp=12.7.fst'; tRange=[35,600]; 
    # SSCase=f'_Hs=2.5_Tp=10';   fstFile = 'simulations/06_Jonswap/OF_F3T1S1_H1A0_Hs=2.5_Tp=10.fst'; tRange=[35,600]; 


    matchSS = re.search(r"Hs=[\d\.]+_Tp=[\d\.]+", fstFile)
    SSCase = matchSS.group(0) if matchSS else None
    shapeFile = f'_data/IEAMonoPile_HydroShapeFunction_{SSCase}.csv'

    hydro_shape_from_pickle(fstFile, IPlot=IPlot, export=export, shapeFile=shapeFile)





    if export and 1 in IPlot:
        export2pdf()



    plt.show()
