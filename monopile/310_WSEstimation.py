import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
import welib.weio as weio
from welib.tools.colors import python_colors, fColrs, lighten_color
# from welib.tools.tictoc import Timer
# from welib.tools.stats import rsquare, mean_rel_err, comparison_stats
from welib.ws_estimator.tabulated import TabulatedWSEstimator
# from welib.ws_estimator.tabulated_floating import *


fstFilename = 'simulations/03_aero/OF_F0T4S1B3_TurbHat.fst'

aeroMapFile = 'simulations/IEA-22-280-RWT/IEA-22-280-RWT_Cp_Ct_Cq.rpf'
# aeroMapFile = 'simulations/IEA-22-280-RWT/IEA-22-280-RWT_CaseGen__CPCTCQ.rpf'
operFile    = 'simulations/IEA-22-280-RWT/IEA-22-280-RWT_OperOpenFAST.csv'

wse = TabulatedWSEstimator(fstFile=fstFilename, operFile = operFile, aeroMapFile=aeroMapFile)
print(wse)

# wse.operPlot()
# plt.show()
# 
# raise


# --- Time series estimation

relaxation=0.9 # To better see failures use 0
tRange = [0,600]
# tRange = [300,305]
method='crossing-oper'
method='crossing'
# method='oper'
debug=False


df = wse.estimateTimeSeriesFromOF(fstFilename.replace('.fst','.outb'), relaxation=relaxation, tRange=tRange, method=method, debug=debug)
wse.plotTimeSeriesEstimation()














plt.show()

print(wse)






if __name__ == '__main_____':
    tRange = None
    tRange = [0,600]
    # tRange = [300,305]
    method='crossing-oper'
    # method='oper'
    debug=False
    relaxation=0.0 # To better see failures use 0


    for fstFilename in fstFilenames:
        # ---
        caseName = os.path.splitext(os.path.basename(fstFilename))[0]


        # pkl = weio.read(pklFilename)
        # print(pkl)

        wse = TabulatedWSEstimatorFloating(fstFile=fstFilename, pickleFile=pklFilename)
        # fig = wse.plotValidSamples()
        # plt.show()

        # --- DEBUG 1
        if -1 in IStudy:
            # --- DEBUG 1 (need clipping)
            omega= 1.4076096156310802; pitch=15.007409964065339; phiy= 0.0
            #pitch=2.081; omega=0.623; phiy=0.000;
            vWS, vQ = wse.TorqueAt(omega=omega, pitch=pitch, phiy=phiy)
            print(vQ)

            wse.debugOP(omega, pitch ,phiy)




        if 0 in IStudy:
            # --- Load .outb file and estimate wind speed based on OpenFAST
            # --- 
        #     fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
        #     fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        #     for om in np.arange(0.1,2, 0.2):
        #         WS,Q = wse.TorqueAt(Pitch=0, Omega=om)
        #         ax.plot(WS, Q, label='Omega = {:3.1f}'.format(om))
        #     ax.set_xlabel('WS [m/s]')
        #     ax.set_ylabel('Q [N]')
        #     ax.set_xlim([0,30])
        #     ax.set_ylim([0,2e7])
        #     ax.legend()
        #     plt.show()

            wse = TabulatedWSEstimatorFloating(fstFile=fstFilename, pickleFile=pklFilename)
            print(wse)

            # --- Loading data
            if onshoreCtrl:
                df =wse.estimateTimeSeriesFromOF(fstFilename.replace('.fst','.CtrlOnshore.outb'), relaxation=relaxation, tRange=tRange, method=method, debug=debug)
            else:
                df =wse.estimateTimeSeriesFromOF(fstFilename, relaxation=relaxation, tRange=tRange, method=method, debug=debug)

            # --- Plot
            fig,axes = plt.subplots(3, 1, sharex=False, figsize=(13.4,7.0)) # (6.4,4.8)
            fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
            
            Ylim1 = [0,20] # WS
            Ylim2 = [0,4.5e6] # Q
            Ylim3 = [-11,25] # Oper


            # Where the data is invalid
            bInv = np.logical_or.reduce((df['Omega_[rad/s]'] <wse.omega[0],   df['Omega_[rad/s]']  >wse.omega[-1]))
            bInv = np.logical_or.reduce((df['Pitch_[deg]']   <wse.pitch[0],   df['Pitch_[deg]']    >wse.pitch[-1], bInv))
            bInv = np.logical_or.reduce((df['PtfmPitch_[deg]']<wse.phiy[0],   df['PtfmPitch_[deg]']>wse.phiy[-1] , bInv))
            bInv = np.logical_or.reduce((df['WS_ref_[m/s]']   <wse.WS[0],     df['WS_ref_[m/s]']   >wse.WS[-1]   , bInv))
            t = df['Time_[s]'].values 


            bInv2 = np.isnan(df['Qaero_est_[N]'])
            #with Timer('ValidValues'):
            #    bInv2 =  ~wse.validValues(df['WS_ref_[m/s]'], df['Omega_[rad/s]'], df['Pitch_[deg]']  , df['PtfmPitch_[deg]'])
            #import pdb; pdb.set_trace()
            bInv3 = np.logical_or(bInv, bInv2)


            bInv = bInv3
            b    = ~bInv3

            # WS plot
            ax=axes[0]
            #ax.fill_between(t, Ylim1[0], Ylim1[1], where=bInv2, alpha=0.1, color=python_colors(1))
            ax.fill_between(t, Ylim1[0], Ylim1[1], where=bInv, alpha=0.1, color=(0.5,0.5,0.5))

            stats, sStats = comparison_stats(t[b], df['WS_ref_[m/s]'].values[b], t[b], df['WS_est_[m/s]'].values[b])
            ax.text(2,Ylim1[0]+(Ylim1[1]-Ylim1[0])*0.89, sStats, fontsize=11 )

            ax.plot(df['Time_[s]'], df['WS_ref_[m/s]'],     color=fColrs(1), label='OpenFAST')
            ax.plot(df['Time_[s]'], df['WS_est_[m/s]'], ':',color=fColrs(4), label='Estimated')
            ax.set_ylabel('Wind speed [m/s]')
            ax.axhline(y = wse.WS[0 ], color=python_colors(0), linestyle = '--', lw=0.5)
            ax.axhline(y = wse.WS[-1], color=python_colors(0), linestyle = '--', lw=0.5)
            ax.set_ylim(Ylim1)
            ax.set_title(caseName.replace('_',' '))

            ax=axes[1]
            #ax.fill_between(t, Ylim2[0], Ylim2[1], where=bInv2, alpha=0.1, color=python_colors(1))
            ax.fill_between(t, Ylim2[0], Ylim2[1], where=bInv, alpha=0.1, color=(0.5,0.5,0.5))

            stats, sStats = comparison_stats(t[b], df['Qaero_ref_[N]'].values[b], t[b], df['Qaero_est_[N]'].values[b])
            ax.text(2,Ylim2[0]+(Ylim2[1]-Ylim2[0])*0.89, sStats, fontsize=11 )

            ax.plot(df['Time_[s]'], df['Qaero_ref_[N]' ]    , color=fColrs(1),   label='OpenFAST')
            ax.plot(df['Time_[s]'], df['Qaero_est_[N]' ],':', color=fColrs(4),  label='Estimated')
            #ax.plot(df['Time_[s]'], df['Qaero_eval_[N]'], '--', label='Evaluated')
            ax.set_ylabel('Qaero [N]')
            ax.set_ylim(Ylim2)
            ax.legend(loc='center left')

            ax=axes[2]
            #ax.fill_between(t, Ylim3[0], Ylim3[1], where=bInv2, alpha=0.1, color=python_colors(1))
            ax.fill_between(t, Ylim3[0], Ylim3[1], where=bInv, alpha=0.1, color=(0.5,0.5,0.5))

            colrs=[python_colors(0), python_colors(1), python_colors(2)]
            colrs=[fColrs(1), lighten_color(fColrs(1),0.3), lighten_color(fColrs(1),0.6)]

            ax.plot(df['Time_[s]'], df['Omega_[rad/s]']*30/np.pi, '-',c=colrs[0], label='Omega [rpm]')
            ax.plot(df['Time_[s]'], df['Pitch_[deg]']           , '--',c=colrs[1], label='Pitch [deg]')
            ax.plot(df['Time_[s]'], df['PtfmPitch_[deg]']       , '-.',c=colrs[2], label='PhiY [deg]')
# 
#             ax.axhline(y = wse.omega[0 ]*30/np.pi, color=colrs[0], linestyle = '--', lw=0.5)
#             ax.axhline(y = wse.omega[-1]*30/np.pi, color=colrs[0], linestyle = '--', lw=0.5)
#             ax.axhline(y = wse.pitch[0 ],          color=colrs[1], linestyle = '--', lw=0.5)
#             ax.axhline(y = wse.pitch[-1],          color=colrs[1], linestyle = '--', lw=0.5)
#             ax.axhline(y = wse.phiy[0 ],           color=colrs[2], linestyle = '--', lw=0.5)
#             ax.axhline(y = wse.phiy[-1],           color=colrs[2], linestyle = '--', lw=0.5)

            ax.set_ylim(Ylim3)
            ax.set_xlabel('Time [s]')
            ax.legend()




            for ax in axes.flatten():
                ax.set_xlim([0,600])
                ax.tick_params(direction='in')

# 
            if onshoreCtrl:
                fig.savefig('_figs/WSEPres_{}_onshoreCtrl.png'.format(caseName))
            else:
                fig.savefig('_figs/WSEPres_{}.png'.format(caseName))



        if 1 in IStudy:



            relaxations=[0,0.3,0.5,0.7]
            Cases=['DBEMT','NoDBEMT','DBEMTStiff']

            for relaxation in relaxations:
                for Case in Cases:

                    base='../data/turbines/NREL5MW'
                    if Case=='DBEMT':
                        Main='../TurbineModel/NREL5MW_SimpleSimulations/TurbWindStep_AllDOF.fst'; 
                    elif Case=='NoDBEMT':
                        Main='../TurbineModel/NREL5MW_SimpleSimulations/TurbWindStep_AllDOF_NoDBEMT.fst'; 
                    elif Case=='DBEMTStiff':
                        Main='../TurbineModel/NREL5MW_SimpleSimulations/TurbWindStep_Stiff.fst'; 

                    # --- Turbine data and estimator
                    wse = TabulatedWSEstimator(fst_file=Main)
                    wse.load_files(base=base,suffix='')
                    print(wse)


                    # --- Loading data
                    df = weio.read(Main.replace('.fst','.outb')).toDataFrame()
                #     df=df[df['Time_[s]']<100]
                    time      = df['Time_[s]']
                    #TTvx      = df['NcIMUTVxs_[m/s]']
                    ws_ref    = df['RtVAvgxh_[m/s]'] # Rotor avg
                    pitch     = df['BldPitch1_[deg]']
                    Qaero_ref = df['RtFldMxh_[N-m]']
                    omega     = df['RotSpeed_[rpm]']*2*np.pi/60 # rad/s
                    lambda_ref=omega*wse.R/ws_ref

                    Qaero = wse.Torque(ws_ref,pitch,omega)

                    # ----
                    print('Estimating...')
                    ws_est=np.zeros(omega.shape)
                    print(ws_est.shape)
                    WS0=ws_ref[0]*0.9 # estimate previous time step
                    for i,(Qa,p,o,ws0) in enumerate(zip(Qaero_ref,pitch,omega,ws_ref)):
                        ws_hat=wse.estimate(Qa, p,o, WS0,relaxation=relaxation)
                        ws_est[i]=ws_hat
                        WS0=ws_hat
                        if np.mod(i,1000)==0:
                            print(i,len(ws),'{:4.1f} {:4.1f}'.format(ws_ref[i],ws[i]))

                    Qaero2 = wse.Torque(ws_est,pitch,omega)

                    # --- Export
                    M=np.column_stack((time,ws_ref,ws,Qaero_ref,Qaero,Qaero2,omega,pitch))
                    header='time,ws_ref,ws,Qaero_ref,Qaero,Qaero2,omega,pitch'
                    np.savetxt('../data/wsest/'+Case+'_{:.1f}.csv'.format(relaxation),M,delimiter=',',header=header)

            # 
        if 2 in IStudy:
            # --- Relative error for each estimations
            files= glob.glob('../data/wsest/*.csv')
            print(files)
            for f in files:
                df=weio.read(f).toDataFrame()
                df=df[df['time']>50]
                df=df[df['time']<600]
                print(f)
                mean_rel_err(df['time'],df['ws_ref'],df['time'],df['ws'])

        if 3 in IStudy:
            base='./SWT-3p6-130'
            Main='../TurbineModel/NREL5MW_SimpleSimulations/TurbWindStep_AllDOF_NoDBEMT.fst'; 
            wse = TabulatedWSEstimator(fst_file=Main)
            wse.load_files(base=base,suffix='')
            print(wse)

        # 
plt.show()

