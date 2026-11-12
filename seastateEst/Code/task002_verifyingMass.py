# Enable Autoreload
import notebook_setup # type: ignore
notebook_setup.enable_autoreload()

# General packages
import numpy as np
import matplotlib.pyplot as plt
import os

# OpenFAST toolbox
from openfast_toolbox.io import FASTOutputFile # type: ignore
from openfast_toolbox.io import FASTInputFile # type: ignore

# Welib
from welib.kalman import kalman

# Seastate Estimator
from seastateEst import FRad
from seastateEst import linData
from seastateEst import read_ssFiles
from seastateEst import wave_filter
from seastateEst import kalman_helper

# --- Structure mass
filepath_mass = r"C:\Jonah\Uni\PhD\OESI\seastateEst\OpenFAST_cases\IEA-15_Full_LINEAR_Falling\main.1.lin"
M_s, M_inv = linData.rigidBodyMass(filepath_mass)
# Added mass at infinite frequency
filepath_Ainf = r"C:\Jonah\Uni\PhD\OESI\seastateEst\OpenFAST_cases\IEA-15_Full\HydroData\IEA-15-240-RWT-UMaineSemi.1"
A_inf = FRad.Ainf(filepath_Ainf, WtrDens=1025)
# Total mass
M = M_s + A_inf

# Alternatively, read mass including added mass from a linearization file
M_alt, M_altInv = linData.rigidBodyMass(r'C:\Jonah\Uni\PhD\OESI\seastateEst\OpenFAST_cases\IEA-15_PlatformOnly_v5.0\main_lin.1.lin')