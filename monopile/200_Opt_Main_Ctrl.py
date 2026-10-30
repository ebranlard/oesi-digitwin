import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
import welib.weio as weio
from welib.fast.servodyn import ServoDyn
# from welib.tools.curve_fitting import gentorque
# from welib.tools.figure import *
# from welib.tools.signal_analysis import applyFilterDF


oper = weio.read('simulations/IEA-22-280-RWT/IEA-22-280-RWT_OperOpenFAST.csv').toDataFrame()
oper['Generator_Torque_[Nm]'] = oper['Generator_Torque_[kNm]']*1000
oper['Generator_Speed_[rpm]'] = oper['Rotor_Speed_[rpm]']


# svd = FASTInputFile('simulations/IEA-22-280-RWT/IEA-22-280-RWT_ServoDyn_Simple.dat')
# svd = ServoDyn('simulations/IEA-22-280-RWT/IEA-22-280-RWT_ServoDyn_Simple.dat')
svd = ServoDyn('C:/Work/_libs/welib/welib/kalman/examples/_simulations/IEA-22-280-RWT/IEA-22-280-RWT_ServoDyn_Simple.dat')


svd.VS_print()

rpm = np.linspace(2, 8, 30)
df = svd.VS_DataFrame(rpm=rpm)

yfit, fitter = svd.VS_fit(oper['Generator_Speed_[rpm]'], oper['Generator_Torque_[Nm]'], verbose=True)
# svd.VS_fit(df['Generator_Speed_[rpm]'], df['Generator_Torque_[Nm]'])


ax = svd.VS_plot(rpm=rpm, label='ServoDyn VS')
# print(oper.keys())
ax.plot(oper['Generator_Speed_[rpm]'], oper['Generator_Torque_[kNm]'], label='Oper')
ax.legend()


plt.show()






























raise

def powerCurvePostProWindStep(outfile, ColMap={}, tWindow=20, WS=None, KeepAll=False, sWS='WS_[m/s]', sTime='Time_[s]', nValsMin=20):
    def renameCol(x):
        for k,v in ColMap.items():
            if x==v:
                return k
        return x

    # --- Reading wind and output file
    print('Reading input...')
    df = weio.FASTOutputFile(outfile).toDataFrame()
    # --- Applying column map
    if len(ColMap)>0:
        for k,v in ColMap.items():
            if v=='':
                v=k
            if v not in df.columns.values:
                raise Exception('Column `{}` not found in output file, cannot apply column mapping'.format(v))
        df.rename(columns=renameCol,inplace=True)
    # --- Our main values
    for v in [sWS, sTime]:
        if v not in df.columns.values:
            raise Exception('Column `{}` need to be present in data. Please provide a column map to define it'.format(v))
    vWS = np.round(df[sWS].values,2)
    t   = df[sTime].values
    # --- Finding unique WS values if not provided
    if WS is None:
        WS=np.unique(np.round(vWS,2))
    # --- Looping on wind values, and averaging data over a window
    tol =0.01
    Iall = np.arange(len(t))
    result=None
    for i,ws in enumerate(WS):
        b = abs(vWS-ws)<tol
        I = Iall[b]
        if len(I)>nValsMin:
            iEnd = I[-1]
            tEnd = t[iEnd]
            tStart = t[iEnd]-tWindow
            IWindow=(t>tStart) & (t<tEnd)
            print('WS:{:5.2f} - t: {:6.1f}-{:6.1f} - nValues:{}'.format(ws,tStart,tEnd,len(IWindow)))
            MeanValues = pd.DataFrame(df[IWindow].mean()).transpose()
            if i==0:
                result = MeanValues.copy()
            else:
                result=result.append(MeanValues, ignore_index=True)

    # --- Sorting by WS
    result[sWS]=np.round(result[sWS],2)
    result.sort_values([sWS],inplace=True)
    result.reset_index(drop=True,inplace=True) 
    # --- Eliminating unnecessary columns
    if len(ColMap)>0 and (not KeepAll):
        ColNames = list(ColMap.keys())
    else:
        ColNames = result.columns.values
    result = result[ColNames]
    
    # --- Rounding
    result=result.round(4)
    return result 



setFigureFont(15)




outFiles = ['SWT-3p6-130-WindStep-CtrlSimp/SWT-3p6-130_WindStep.outb']
# outFiles +=['SWT-3p6-130-ROSCO_Model/SWT-3p6-130_WindStep_Offshore.outb']
# outFiles +=['SWT-3p6-130-ROSCO_Model/SWT-3p6-130_WindStep.OffshoreController.outb']
# outFiles +=['SWT-3p6-130-ROSCO_Model/SWT-3p6-130_WindStep_Offshore.Off_om010.out']
# outFiles +=['SWT-3p6-130-ROSCO_Model/SWT-3p6-130_WindStep_PowerTracking_Onshore.outb']
# outFiles +=['SWT-3p6-130-ROSCO_Model/SWT-3p6-130_WindStep_Offshore_Long.Ctrlv1.outb']
# outFiles +=['SWT-3p6-130-ROSCO_Model/SWT-3p6-130_WindStep_Offshore_Long.Ctrlv2.outb']
outFiles +=['SWT-3p6-130-ROSCO_Model/SWT-3p6-130_WindStep_Offshore_Long.outb']
measFile = '../../TetraSpar_data/Oper.csv'

P_rated_elec  = 3.6e6  #[W]
DT_Efficiency = 1.00  # drivetrain efficiency
nGear         = 1     # Gear ratio
rpm_rated     = 12.30 # TODO low speed shaft RPM rated
rpm_rated2    = 13.30 # TODO low speed shaft RPM rated
rho           = 1.225
R             = 65.0852
lambda_opt    = 9.4714 # Tip speed ratio at "optimum" 
Cp_opt        = 0.4214 # Power coefficient at "optimum"

SpdGenOn  = 6.6 # HSS [rpm]
SlPc      = 3                    # Rated generator slip percentage in Region 2 1/2 for simple variable-speed generator control (%) 

# --- Derived params
RPM2RADS=   2*np.pi/60;
RADS2RPM=1/(2*np.pi/60);

omega_hss_rated = rpm_rated * nGear * RPM2RADS
print('omega_hss_rated : ',omega_hss_rated, 'rad/s')

k     = DT_Efficiency * 0.5*rho*np.pi*R**5/lambda_opt**3 *Cp_opt   # At the rotor, in Nm/(rad/s)^2
k_hss = k/(nGear**3)                               # at the HSS in Nm/(rad/s)^2
k_rpm = k/(nGear**3*RADS2RPM**2)                    # at the HSS in Nm/rpm^2

# Rgn2K  = k_rpm
Rgn2K  = 12222.0 # % Generator torque constant in Region 2 for simple variable-speed generator control (HSS side) (N-m/rpm^2) 
print('Region 2 cst: ', np.around(k_rpm,6),'(from lambda_opt)')
print('Region 2 cst: ', np.around(Rgn2K,6),'(used)')
print('Region 2 cst: ', np.around(Rgn2K*RADS2RPM**2,6),'(in n-m/(rad/s^))')

# --- Derived parameters
P_aero = P_rated_elec / DT_Efficiency
Q_aero = P_aero / (rpm_rated2*RPM2RADS)


RtGnSp = rpm_rated * nGear # Rated generator speed for simple variable-speed generator control (HSS side) (rpm) 
RtTq= Q_aero*DT_Efficiency # Rated generator torque/constant generator torque in Region 3 for simple variable-speed generator control (HSS side) (N-m) [used only when VSContrl=1]



# --- "model" based on params above
GenSpeed_mdl = np.linspace(SpdGenOn*0.5,RtGnSp*1.3,200)                # RPM, HSS
GenTq_mdl    = gentorque(GenSpeed_mdl, (RtGnSp, RtTq  , Rgn2K , SlPc , SpdGenOn))

# --- Based on simulation
outFile=outFiles[0]

# df = weio.read(outFile).toDataFrame()
df = powerCurvePostProWindStep(outFile, ColMap={}, tWindow=40, WS=None, KeepAll=False, sWS='Wind1VelX_[m/s]', sTime='Time_[s]', nValsMin=20)
GenSpeed_sim = df['RotSpeed_[rpm]'].values
GenTq_sim    = df['GenTq_[kN-m]'].values*1000
WS_sim       = df['Wind1VelX_[m/s]'].values
Pitch_sim    = df['BldPitch1_[deg]'].values
Power_sim    = df['GenPwr_[kW]'].values


df2 = powerCurvePostProWindStep(outFiles[1], ColMap={}, tWindow=20, WS=None, KeepAll=False, sWS='Wind1VelX_[m/s]', sTime='Time_[s]', nValsMin=20)
GenSpeed_sim2= df2['RotSpeed_[rpm]'].values
GenTq_sim2   = df2['GenTq_[kN-m]'].values*1000
WS_sim2      = df2['Wind1VelX_[m/s]'].values
Pitch_sim2   = df2['BldPitch1_[deg]'].values
Power_sim2   = df2['GenPwr_[kW]'].values


# --- Based on measurements
df = weio.read(measFile).toDataFrame()
GenSpeed_mea = df['RPM_[rpm]'].values
GenTq_mea    = df['Qgen_[kNm]'].values*1000
WS_mea       = df['WS_[m/s]'].values
Pitch_mea    = df['Pitch_[deg]'].values
Power_mea    = df['P_[kW]'].values

# --- Fit based on data (TODO)
# bounds='RtGnSp=(10.4,11.8) , RtTq=(4.9e6,5.2e6), Rgn2K=(10000.0,30000) ,SlPc=(0,15) , SpdGenOn=(1,6)'
# p0 = [RtGnSp, RtTq, Rgn2K, SlPc, SpdGenOn]
#         y_fit, pfit, fitter = model_fit('fitter: gentorque', GenSpeed_mea, GenTq_mea)
# y_fit, pfit, fitter = model_fit('predef: gentorque', GenSpeed_mea, GenTq_mea) #, bounds=bounds, p0=p0)

# fitter = GeneratorTorqueFitter(GenSpeed_mea, GenTq_mea, p0=p0, bounds=bounds)
# fitter = GeneratorTorqueFitter(GenSpeed_sim, GenTq_sim, p0=p0, bounds=bounds)
# print(fitter)
# 
# fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.6)) # (6.4,4.8)
# fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
# ax.plot(GenSpeed_mdl, GenTq_mdl     ,'k-', label='Model')
# ax.plot(GenSpeed_mea, GenTq_mea     ,'.', label='Measurements')
# ax.plot(GenSpeed_sim, GenTq_sim     ,'--', label='Simple Torque Ctrl')
# ax.plot(GenSpeed_sim2, GenTq_sim2     ,'--', label='ROSCO')
# 
# 
# # ax.plot(GenSpeed_mdl, fitter.model['fitted_function'](GenSpeed_mdl),'--', label='Fit')
# # ax.plot(x, y_fit ,'-', label='')
# # ax.plot(x, fitter.model['fitted_function'](x) ,'.', label='')
# ax.set_xlabel('Generator speed [rpm]')
# ax.set_ylabel('Generator Torque [kNm]')
# ax.legend()
# ax.tick_params(direction='in')
# fig.savefig('ControllerQGen.png')

# 
# fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
# fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
# # ax.plot(GenSpeed_mea*RPM2RADS, GenTq_mea     ,'.', label='Measurements')
# ax.plot(GenSpeed_mdl*RPM2RADS, GenTq_mdl     ,'-', label='Model')
# # ax.plot(GenSpeed_sim*RPM2RADS, GenTq_sim     ,'--', label='Simulation')
# # ax.plot(GenSpeed_mdl, fitter.model['fitted_function'](GenSpeed_mdl),'--', label='Fit')
# # ax.plot(x, y_fit ,'-', label='')
# # ax.plot(x, fitter.model['fitted_function'](x) ,'.', label='')
# ax.set_xlabel('Generator speed [rad/s]')
# ax.set_ylabel('Generator Torque ')
# ax.legend()
# ax.tick_params(direction='in')
# 

# --- Wrt WS
XLIM=[0,19] 
fig,axes = plt.subplots(2, 3, sharey=False, figsize=(16.8,4.6)) # (6.4,4.8)
fig.subplots_adjust(left=0.07, right=0.98, top=0.95, bottom=0.11, hspace=0.34, wspace=0.32)
ax = axes[0,0]
ax.plot(WS_mea, Power_mea     ,'.', label='Measurements')
ax.plot(WS_sim, Power_sim     ,'--', label='Simple Torque Ctrl')
ax.plot(WS_sim2, Power_sim2   ,'--', label='ROSCO')
ax.set_xlabel('Wind speed [m/s]')
ax.set_ylabel('Generator Power [kW]')
ax.set_ylim([0,4000])
ax.legend(loc='lower right', fontsize=12)
ax.set_xlim(XLIM)
ax = axes[0,1]
ax.plot(WS_mea, GenTq_mea     ,'.', label='Measurements')
ax.plot(WS_sim, GenTq_sim     ,'--', label='Simple Torque Ctrl')
ax.plot(WS_sim2, GenTq_sim2    ,'--', label='ROSCO')
ax.set_xlabel('Wind speed [m/s]')
ax.set_ylabel('Generator Torque [kNm]')
ax.set_xlim(XLIM)

ax = axes[0,2]
ax.plot(GenSpeed_mdl, GenTq_mdl     ,'k-', label='Simple Torque Model')
ax.plot(GenSpeed_mea, GenTq_mea     ,'.', label='Measurements')
ax.plot(GenSpeed_sim, GenTq_sim     ,'--', label='Simple Torque Ctrl')
ax.plot(GenSpeed_sim2, GenTq_sim2   ,'--', label='ROSCO')
ax.legend(loc='lower center', bbox_to_anchor=(0.5,-1.35))
ax.set_xlabel('Generator speed [rpm]')
ax.set_ylabel('Generator Torque [kNm]')
ax.set_xlim([0,15])
ax.set_xlim(XLIM)

ax = axes[1,0]
ax.plot(WS_mea, GenSpeed_mea     ,'.', label='Measurements')
ax.plot(WS_sim, GenSpeed_sim     ,'--', label='Simple Torque Ctrl')
ax.plot(WS_sim2, GenSpeed_sim2     ,'--', label='ROSCO')
ax.set_xlabel('Wind speed [m/s]')
ax.set_ylabel('Generator Speed [rpm]')
ax.set_ylim([0,15])
ax.set_xlim(XLIM)

ax = axes[1,1]
ax.plot(WS_mea, Pitch_mea     ,'.', label='Measurements')
ax.plot(WS_sim, Pitch_sim     ,'--', label='Simple Torque Ctrl')
ax.plot(WS_sim2, Pitch_sim2    ,'--', label='ROSCO')
ax.set_xlabel('Wind speed [m/s]')
ax.set_ylabel('Pitch [deg]')
ax.set_xlim(XLIM)

ax = axes[1,2]
ax.axis('off')

for ax in axes.flatten():
    ax.tick_params(direction='in')

fig.savefig('ControllerOperAll.png')



print('SpdGenOn : {:20.3f} [rpm]'      .format(SpdGenOn))
print('VS_RtGnSp: {:20.3f} [rpm]'      .format(RtGnSp),', ',np.around(RtGnSp*RPM2RADS,4),'rad/s')
print('VS_RtTq  : {:20.3f} [Nm]'       .format(RtTq))
print('VS_Rgn2K : {:20.3f} [N-m/rpm^2]'.format(Rgn2K))
print('VS_SlPc  : {:20.3f} [%]'        .format(SlPc))
print('')
print('VS_Rgn2K : {:20.3f} [N-m/(rad/s)^2]'.format(Rgn2K*RADS2RPM**2))
print('')
RtTqCheck = Rgn2K*RtGnSp**2
print('TqCheck  : {:20.3f} [Nm]'.format(RtTqCheck))
if RtTqCheck>RtTq:
    print('[WARN] Rgn2K*RtGnSp**2 ({})  > RtTq ({})'.format(RtTqCheck, RtTq))




plt.show()

if __name__ == '__main__':
    pass
