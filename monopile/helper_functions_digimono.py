import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
# Local 
import welib.weio as weio
from welib.essentials import *
import pickle
from welib.tools.colors import MW_Orange, fColrs

# 
COLRS   = [fColrs(1), MW_Orange ,fColrs(2), fColrs(3)]
try:
    cmap = matplotlib.colormaps['viridis']
except:
    cmap = matplotlib.cm.get_cmap('viridis')


CMAP = [(cmap(v)[0],cmap(v)[1],cmap(v)[2],0.9) for v in np.linspace(0,1,5+1)]
COLRS=[CMAP[1], CMAP[4]]
colRef=COLRS[0]
colSim=COLRS[1]
LWRef=2.4
LWSim=1.5
# 
np.set_printoptions(linewidth=300, precision=4)
pd.set_option("display.max_rows", 999)
pd.set_option("display.precision", 4)


def meanabs(x, **kwargs):
    return np.mean(np.abs(x), **kwargs)

def time_plot(t, ref=None, sim=None, label='', other=None, fig=None, ax=None, tRange=None, refLab='OpenFAST', otherLab='Other'):
    figNotProvided = fig is None
    if fig is None:
        fig=plt.figure()
        fig.subplots_adjust(left=0.18, right=0.94, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        ax = fig.add_subplot(111)
    if ref is not None:
        ax.plot(t, ref,'-'  , c=colRef, lw=LWRef, label=refLab)
    if sim is not None:
        ax.plot(t, sim,'--' , c=colSim, lw=LWSim, label='YAMS')
    if other is not None:
        ax.plot(t, other,':' , label=otherLab)
    if tRange is not None:
        ax.axvline(x = tRange[0] , color = 'k', ls='--')
        ax.axvline(x = tRange[1], color = 'k', ls='--')
    if figNotProvided:
        ax.legend()
        ax.set_xlabel('Time [s]')
    ax.set_ylabel(label)
    ax.tick_params(direction='in')
    return fig, ax


def sec_plotFM(vTime, zDepth, F_sec, F_secRef=None, stat='mean', tRange=None, label='Section Force', other=None, otherLab='Other'):
    if tRange is None:
        IT = np.arange(0,F_sec.shape[1])
    else:
        IT = np.logical_and(vTime>tRange[0], vTime<tRange[1])
        if len(IT)==0:
            IT = np.arange(0, sec.shape[1])
    fstat = {'mean':np.mean, 'max':np.max, 'std':np.std, 'meanabs':meanabs}[stat]
    if F_sec is not None:
        Fx_sec0 = fstat(F_sec   [0,:,IT],axis=0) # 3, nSpan nT
        My_sec0 = fstat(F_sec   [4,:,IT],axis=0)
    if F_secRef is not None:
        Fx_sec1 = fstat(F_secRef[0,:,IT],axis=0) # 3, nSpan nT
        My_sec1 = fstat(F_secRef[4,:,IT],axis=0)
    if other is not None:
        Fx_sec2 = fstat(other[0,:,IT],axis=0) # 3, nSpan nT
        My_sec2 = fstat(other[4,:,IT],axis=0)

    fig,axes = plt.subplots(1, 2, sharey=True, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.14, right=0.95, top=0.95, bottom=0.12, hspace=0.20, wspace=0.20)
    ax=axes[0]
    if F_secRef is not None:
        ax.plot(Fx_sec1/1e6, zDepth,  '-' , c=colRef, lw=LWRef, label='OpenFAST')
    if F_sec is not None:
        ax.plot(Fx_sec0/1e6, zDepth,  '--', c=colSim, lw=LWSim, label='YAMS')
    if other is not None:
        ax.plot(Fx_sec2/1e6, zDepth,  ':' ,label=otherLab)

    ax.tick_params(direction='in')

    ax.set_xlabel(label)
    ax.set_ylabel('Vertical position [m]')

    ax=axes[1]
    if F_secRef is not None:
        ax.plot(My_sec1/1e6, zDepth, '-',  c=colRef, lw=LWRef, label='OpenFAST')
    if F_sec is not None:
        ax.plot(My_sec0/1e6, zDepth, '--', c=colSim, lw=LWSim, label='YAMS')
    if other is not None:
        ax.plot(My_sec2/1e6, zDepth,  ':' ,label=otherLab)
    ax.set_xlabel('Section Moment My [MNm]')
    ax.tick_params(direction='in')
    if F_sec is not None:
        ax.legend()
    fig.suptitle(stat)
    return fig

def sec_plot(vTime, zDepth, sec, ref=None, other=None, stat='mean', tRange=None, label='Section Force', fig=None, refLab='OpenFAST', otherLab='Other'):
    if tRange is None:
        IT = np.arange(0, sec.shape[1])
    else:
        IT = np.logical_and(vTime>tRange[0], vTime<tRange[1])
        if len(IT)==0:
            IT = np.arange(0, sec.shape[1])
    fstat = {'mean':np.mean, 'max':np.max, 'std':np.std, 'meanabs':meanabs}[stat]
    sec0 = fstat(sec[:,IT], axis=1)
    if ref is not None:
        sec1 = fstat(ref[:,IT], axis=1)
    if other is not None:
        sec2 = fstat(other[:,IT], axis=1)

    if fig is None:
        fig, ax = plt.subplots(1, 1, sharey=True, figsize=(6.4,4.8)) # (6.4,4.8)
        fig.subplots_adjust(left=0.14, right=0.95, top=0.95, bottom=0.12, hspace=0.20, wspace=0.20)
    else:
        ax=fig.axes[0]

    if ref is not None:
        ax.plot(sec1, zDepth,  '-' , c=colRef, lw=LWRef, label=refLab)
    ax.plot(sec0,     zDepth,  '--', c=colSim, lw=LWSim, label='YAMS')
    if other is not None:
        ax.plot(sec2, zDepth,  '-' , c=colRef, lw=LWRef, label=refLab)
    ax.tick_params(direction='in')
    ax.legend()

    ax.set_xlabel(label)
    ax.set_ylabel('Vertical position [m]')
    fig.suptitle(stat)
    return fig




# def save(dat,filename):
#     with open(filename,'wb') as f:
#         pickle.dump(dat,f)
# 
# def load(filename):
#     with open(filename,'rb') as f:
#         dat=pickle.load(f)
#     return dat

# def Leq(t,y,m=5):
#     return eq_load(y, m=m, neq=t[-1]-t[0])[0][0]
# 
# def SNR(y):
#     return np.mean(y**2)/np.std(y)**2
# 
# def moving_average(a, n=3) :
#     """ 
#     perform moving average, return a vector of same length as input
#     """
#     a=a.ravel()
#     a = np.concatenate(([a[0]]*(n-1),a)) # repeating first values
#     ret = np.cumsum(a, dtype=float)
#     ret[n:] = ret[n:] - ret[:-n]
#     ret=ret[n - 1:] / n
#     return ret
# 
# 
def mean_rel_err(t1,y1,t2,y2, tRange=None):
    if tRange is not None: 
        b=np.logical_and(t1>tRange[0], t1<tRange[1])
        if sum(b)>0:
            t1 = t1[b]
            y1 = y1[b]
            b=np.logical_and(t2>tRange[0], t2<tRange[1])
            t2 = t2[b]
            y2 = y2[b]

    if len(y1)!=len(y2):
        y2=np.interp(t1,t2,y2)
    # Method 1 relative to mean
    ref_val = np.mean(y1)
    meanrelerr0=np.mean(np.abs(y1-y2)/ref_val)*100 
    #print('Mean rel error {:7.2f} %'.format( meanrelerr0))
    # Method 2 scaling signals
    Min=min(np.min(y1), np.min(y2))
    Max=max(np.max(y1), np.max(y2))
    y1=(y1-Min)/(Max-Min)+0.001
    y2=(y2-Min)/(Max-Min)+0.001
    meanrelerr=np.mean(np.abs(y1-y2)/np.abs(y1))*100 
    print('Mean rel error {:7.2f} %'.format( meanrelerr))
    return meanrelerr,meanrelerr0

# 
# def mean_rel_err(t1,y1,t2,y2):
#     if len(y1)!=len(y2):
#         y2=np.interp(t1,t2,y2)
#     # Method 1 relative to mean
#     ref_val = np.mean(y1)
#     meanrelerr0=np.mean(np.abs(y1-y2)/ref_val)*100 
#     print('Mean rel error {:7.2f} %'.format( meanrelerr0))
#     # Method 2 scaling signals
#     Min=min(np.min(y1), np.min(y2))
#     Max=max(np.max(y1), np.max(y2))
#     y1=(y1-Min)/(Max-Min)+0.001
#     y2=(y2-Min)/(Max-Min)+0.001
#     meanrelerr=np.mean(np.abs(y1-y2)/np.abs(y1))*100 
#     print('Mean rel error {:7.2f} %'.format( meanrelerr))
#     return meanrelerr,meanrelerr0
# 
# def findsignal(sL, pat):
#     I=[]
#     for i,s in enumerate(sL):
#         if s.find(pat)>=0:
#             print(i,s)
#             I.append(i)
#     return I
# # --------------------------------------------------------------------------------}
# # --- Measurements 
# # --------------------------------------------------------------------------------{
# def to_global(s,theta_s):
#     xg=  s*np.cos(theta_s*np.pi/180)
#     yg= -s*np.sin(theta_s*np.pi/180)
#     return xg,yg
# 
# def to_nac(s,theta_s,theta_n):
#     xg,yg=to_global(s,theta_s)
#     xn = xg*np.cos(theta_n*np.pi/180) + yg*np.sin(theta_n*np.pi/180)
#     yn = yg*np.cos(theta_n*np.pi/180) - xg*np.sin(theta_n*np.pi/180)
#     return xn,yn
# 
# def to_nac_combine(s1,theta_s1,s2,theta_s2,theta_n):
#     x1,y1 = to_nac(s1, theta_s1, theta_n)
#     x2,y2 = to_nac(s2, theta_s2, theta_n)
#     return x1+x2,y1+y2
# 
# def tdms_to_csv(filename, TimeFrame=None):
#     nGear=78.68
#     GenEff=0.95
# 
#     df=weio.TDMSFile(filename).toDataFrame()
#     if TimeFrame is not None:
#         df=df[np.logical_and(df['Time_[s]']>TimeFrame[0], df['Time_[s]']<TimeFrame[1])]
#         df['Time_[s]'] -=TimeFrame[0]
# 
# 
#     df['TorqueLSS_[kNm]']    = (df['Mainshaft_Downwind_Torque'] - 1.1511) * 1274 # [kNm]
#     df['TorqueHSS_[kNm]']    = df['TorqueLSS_[kNm]']/nGear
#     df['GenTq_[kNm]'] = df['Active Power']/(df['HSS RPM']*2*np.pi/60)/GenEff
#     df['TwrBsM1_[kNm]']     = (df['Tower Base Bend 1'] -(-2.5243E-05) ) * 4.9136E+07 # [kNm]
#     df['TwrBsM2_[kNm]']     = (df['Tower Base Bend 2'] -(-3.0322E-04) ) * 4.8912E+07 # [kNm]
# 
#     df.rename(columns={'WindSpeed_80m':'WS80_[m/s]'}, inplace=True)
#     df.rename(columns={'Wind_Direction_38m':'WD38_[deg]'}, inplace=True)
#     df.rename(columns={'WD1_87m':'WD87_[deg]'}, inplace=True)
#     df.rename(columns={'LSS RPM':'RotSpeed_[rpm]'}, inplace=True)
#     df.rename(columns={'HSS RPM':'GenSpeed_[rpm]'}, inplace=True)
#     df.rename(columns={'TowerTopDCAccel_NS':'TTAcc_NS_[g]'}, inplace=True)
#     df.rename(columns={'TowerTopDCAccel_EW':'TTAcc_EW_[g]'}, inplace=True)
#     df.rename(columns={'Pitch_Blade1':'BldPitch1_[deg]'}, inplace=True)
#     df.rename(columns={'Pitch_Blade2':'BldPitch2_[deg]'}, inplace=True)
#     df.rename(columns={'Yaw_Encoder':'Yaw_[deg]'}, inplace=True)
#     df.rename(columns={'AzimuthAO':'AzimuthAO_[deg]'}, inplace=True)
#     df.rename(columns={'Azimuth':'Azimuth_[deg]'}, inplace=True)
#     df.rename(columns={'Temp1':'AirTemp_[C]'}, inplace=True)
#     df.rename(columns={'Air_Press_1':'Patm_[kPa]'}, inplace=True)
# 
#     DropCols= [ 'Scan Errors','Late Scans', 'OPC_AI_In_GridMonRealPowerAct', 'OPC_AI_In_PitchAngleCurrent1', 'OPC_CMD_Close_Brake', 'OPC_DI_In_TBDiscBrakeClosed', 'OPC_In_RotorSpd', 'OPC_In_WindSpd', 'OPC_OpCtl_TurbineFullState', 'OPC_OpCtl_TurbineOperationalState', 'OPC_STATE_FAULT', 'OPC_TurbineStatusSCADA', 'OPC_Yaw_AsPublished', 'WD_Mod_Active', 'WD_Nacelle', 'WD_Nacelle_Mod',
#     'MS Excel Timestamp', 'LabVIEW Timestamp', 'Accel_Generator_Front', 'Accel_Generator_Back', 'Accel_HSS_rotor', 'Accel_HSS_generator', 'Accel_TorqueArm1_H', 'Accel_TorqueArm1_V', 'Accel_TorqueArm2_H', 'Accel_TorqueArm2_V', 'Accel_Planetary', 'Accel_TorqueArm1_Frame']
#     DropCols+= ['Tower Base Bend 1','Tower Base Bend 2','Mainshaft_Downwind_Torque']
#     DropCols+= ['AzimuthAO_[deg]']
#     DropMiss = [c for c in DropCols if c not in df.columns]
#     DropCols = [c for c in DropCols if c in df.columns]
#     print('Missing columns:',DropMiss)
#     df=df.drop(columns=DropCols)
# 
# 
#     # --- NOTE wind direction is counted negative along z
#     # yaw is -200 i.e. +160 deg,
#     # the axis pointing upwind is thus 160 deg (clockwise, wind convention)
#     # the axis pointing downwind is thus -20 deg (clockwise, wind convention)
#     theta_n=df['Yaw_[deg]']+180 # Offset of 180 since x axis is pointing downstream
#     
#     # NOTE: acceleration always in nacelle coordinate
#     accFA,accSS = to_nac_combine(df['TTAcc_NS_[g]'] , 186, df['TTAcc_EW_[g]'] , 96, theta_n=theta_n + 20) # NOTE offset
#     df['TTAcc_FA_[g]'] = accFA
#     df['TTAcc_SS_[g]'] = accSS
#     #M1,M2     = to_nac_combine(df['TwrBsM1_[kNm]'], 6, df['TwrBsM2_[kNm]'], 96 , theta_n=theta_n)
#     M1,M2     = to_nac_combine(df['TwrBsM1_[kNm]'], 276, df['TwrBsM2_[kNm]'], 186 , theta_n=theta_n-12.5) # NOTE offset
#     df['TwrBsM_FA_[kNm]'] = M1
#     df['TwrBsM_SS_[kNm]'] = M2
# 
# 
# #     for theta_n in np.arange(-90,90,10):
# #         sTheta='{}'.format(theta_n)
# #         accFA,accSS = to_nac_combine(df['TTAcc_NS_[g]'] , 186, df['TTAcc_EW_[g]'] , 96, theta_n=theta_n)
# #         df['TTAcc_FA_[g]'+sTheta] = accFA
# #         df['TTAcc_SS_[g]'+sTheta] = accSS
# # 
# #         #MFA,MSS     = to_nac_combine(df['TwrBsM1_[kNm]'], 276, df['TwrBsM2_[kNm]'], 186 , theta_n=theta_n)
# #         M1,M2     = to_nac_combine(df['TwrBsM1_[kNm]'], 276, df['TwrBsM2_[kNm]'], 186 , theta_n=theta_n)
# #         MFA=M1
# #         MSS=M2
# #         df['TwrBsM_FA_[kNm]'+sTheta] = MFA
# #         df['TwrBsM_SS_[kNm]'+sTheta] = MSS
# #         M1,M2     = to_nac_combine(df['TwrBsM1_[kNm]'], 276, df['TwrBsM2_[kNm]'], 6 , theta_n=theta_n)
# #         MFA=M2
# #         MSS=M1
# #         df['TwrBsM_FA2_[kNm]'+sTheta] = MFA
# #         df['TwrBsM_SS2_[kNm]'+sTheta] = MSS
# 
#     df.to_csv(filename.replace('.tdms','.csv'),index=False)
# 
# 
# 
