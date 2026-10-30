from welib.fast.subdyn import SubDyn

sdfile = 'simulations/00_EVA_DEBUG/SD_CB0_NoCM.dat'

SD = SubDyn(sdfile)
SD.init(TP=(0,0,15))

# FEM1 = SD.beamFEM(method='cbeam', force=True)
# FEM2 = SD.beamFEM(method='full', force=True)
# # print(FEM2)
# Q_G1,_Q_CB1, df_G1, df_CB1, Modes_G1, Modes_CB1, CB1 = SD.beamModes(nCB=0, method='cbeam')
# Q_G2,_Q_CB2, df_G2, df_CB2, Modes_G2, Modes_CB2, CB2 = SD.beamModes(nCB=0, method='full')

# p1 = SD.toYAMSData(method='cbeam', equispacing=False)
p1, damp_zeta1, RayleighCoeff1, DampMat1, df_G1 = SD.toYAMSData(method='cbeam' , equispacing=False)
p2, damp_zeta2, RayleighCoeff2, DampMat2, df_G2 = SD.toYAMSData(method='full' , equispacing=False)
df_G1.to_csv('_Guyan1.csv',index=False)
df_G2.to_csv('_Guyan2.csv',index=False)

# print(df_G1)
# print(df_G2)
# 

# for i, f0 in enumerate(FEM1['freq']):
#     print('Mode {:4d} {:7.3f} {:7.3f}'.format(i, FEM1['freq'][i], FEM2['freq'][i]))
#     if i>5:
#         break

# for i, f0 in enumerate(FEM['freq']):
#     print('Mode {:4d} {:7.3f} {:7.3f}'.format(i, FEM['freq'][i], SD._FEM.freq[i]))
#     if i>5:
#         break

