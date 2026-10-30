import os
import numpy as np

from numpy import trapezoid
from welib.hydro.morison import monopileHydroLoads1D
from welib.yams.section_loads import beamSectionLoadsFromShapeFunctions
from welib.essentials import *

import welib.weio as weio
from welib.yams.models.MTNSB import FASTmodel2MTNSB
from welib.yams.windturbine import FASTWindTurbine
from welib.system.mech_system import MechSystem
from welib.hydro.wavekin import wavenumber, elevation2d, wave_components



# --------------------------------------------------------------------------------}
# --- Monopile forces
# --------------------------------------------------------------------------------{
def monopileGF(t, q, qd, pST, pSS, pHD, qdd=None, reconHydro=False, nOut=[0], eta=None, eta_dot=None):
    """ 
    Compute generalized forces in degrees of freedom.
    INPUTS:
     - q : mechanical degrees of freedom (e.g. q=[q_s, q_d])
     - qd: velocities of mechanical degrees of freedom (e.g. qd=[qd_s, qd_d])
     - pSS: dictionary of sea state parameters
           keys: 'ap', 'fp', 'kp', 'epsp', 'WaterDepth', 'rho'
     - pHD: dictionary of hydrodynamic parameters
           keys: 'zDepth', 'D', 'CD', 'CM'
     - pST: dictionary of structural parameters
           keys: 'z', 'PhiU'

    NOTE: for degrees of freedom, not first order state vector
        q    = q[0:self.nDOF]
        qdot = q[self.nDOF:]

    Return:
     - generalized forces GF=shape(x)
     - outF : additional outputs
    """
    q  = q.flatten()
    qd = qd.flatten()
    GF = np.zeros_like(q) # Generalized force

    # Structural velocity
    PhiU     = pST['PhiU']
    ndim, nz = PhiU[0].shape # Number of dimensions (3) and spanwise positions
    v_struct = np.zeros((ndim, nz))
    for j, PhiUj in enumerate(PhiU):
        v_struct  += qd[j] * PhiUj
    # Structural Acceleratin
    a_struct = np.zeros((ndim, nz))
    if qdd is not None:
        for j, PhiUj in enumerate(PhiU):
            a_struct  += qdd[j] * PhiUj
    xdist=0

    # --- Hydro force
    if pSS is None:
        outH=None
    else:
        if eta_dot is None:
            eta     = np.interp(t, pSS['eta_time'], pSS['eta'])
            eta_dot = np.interp(t, pSS['eta_time'], pSS['eta_dot'])
        if reconHydro:
            outH=dict()
            nOut[0]+=1
            if nOut[0]<10:
                print('recon', t)
            p_hydro = pHD['phi'] * eta_dot

            # Limit calculations to below water and above sea bed
            z = pHD['zDepth']
            h = -np.min(z) # approx 
            bWet = np.logical_and(z<=0, z>=-h)  # TODO eta
            zWet = z[bWet].flatten()
            z_ref=np.min(zWet) 

            outH['p_hydro_no_struct_acc']  = p_hydro.copy()
            if qdd is not None:
                A = (np.pi * pHD['D']**2 / 4).ravel()
                p_AM_str  = -pSS['rho'] * pHD['Ca'].ravel() * A * a_struct[0,:].ravel()  # Added-mass (relative acceleration)
                p_AM_str[~bWet]  = 0 # Important
                p_hydro += p_AM_str

            outH['p_hydro']  = p_hydro
            outH['dM_hydro'] = p_hydro * (z-z_ref)
            outH['F_hydro']    = trapezoid(p_hydro[bWet]               , zWet) # [N]
            outH['M_sb_hydro'] = trapezoid(p_hydro[bWet] * (zWet-z_ref), zWet) # Sea bed moment [Nm]

        else:
            outH = monopileHydroLoads1D(t, pSS['ap'], pSS['fp'], pSS['kp'], pSS['epsp'], pSS['WaterDepth'], pHD['zDepth'], xdist, D=pHD['D'], Cd=pHD['Cd'], Cp=pHD['Cp'], Ca=pHD['Ca'] , rho = pSS['rho'], u_struct=v_struct[0,:], a_struct=a_struct[0,:])

        outH['eta']     = eta
        outH['eta_dot'] = eta_dot

        GF_s = np.trapezoid(outH['p_hydro']*PhiU[0][0,:], pHD['zDepth']) # Generalized force in surge
        GF_p = np.trapezoid(outH['p_hydro']*PhiU[1][0,:], pHD['zDepth']) # Generalized force in pitch
        GF[0]+=GF_s
        GF[1]+=GF_p

    outF=outH
#     outF['p_hydro'] = outH['p_hydro']
#     outF['F_hydro'] = outH['F_hydro']
#     outF['M_hydro'] = outH['M_sb_hydro']

    return GF, outF


# --------------------------------------------------------------------------------}
# --- Calc outputs 
# --------------------------------------------------------------------------------{
def monopileOutputs(t, x, pST, pSS, pHD, it=None, xd=None, Sys=None, reconHydro=False):
    """ 
    Compute Outputs:
      - Section Forces

    INPUTS:
     - x : full state vector (e.g. x=[q_s, q_p, qd_s, qd_p])

     - pSS: dictionary of sea state parameters
           keys: 'ap', 'fp', 'kp', 'epsp', 'WaterDepth', 'rho'
     - pHD: dictionary of hydrodynamic parameters
           keys: 'zDepth', 'D', 'CD', 'CM'
     - pST: dictionary of structural parameters
           keys: 's_span', 'PhiU', 'PhiV', 'm'


     - First order state vector: x:[q_s, q_p, qd_s, qd_p] (shape 2*nDOF)
    """
    # --- Local variables, allocs
    nDOF = len(x)//2
    x = x.reshape(-1,1) # column vector of states
    q = x[:nDOF]
    qd = x[nDOF:]


    # Velocity - Acceleration
    if xd is None:
        xd  = Sys.dqdt (t, q=x.flatten()) # Evaluate acceleration (solve_ivp doesn't store it)
    qdd = xd[nDOF:]


    # --- External forces (general and non general)
    # IMPORTANT:
    #  - If you provide qdd, then the added mass from the structural acceleration will be in p_hydro
    #    In that case, we should not provide m_hydro to the section loads (to avoid double counting)
    GF, outF = monopileGF(t, q, qd, pST=pST, pSS=pSS, pHD=pHD, qdd=qdd, reconHydro=reconHydro)


    # Main Inputs for section loads
    if pSS is not None:
        p_ext = outF['p_hydro']
    else:
        p_ext = np.zeros((3, len(pST['m'])))
    a_ext = np.array((0.,0.,-pST['gravity'])) # external acceleration (gravity/earthquake)

    # Loads at Foundation top, in foundation coordinates
    # NOTE: see other scripts from "code3_fixedbottom" for how to do this properly.
    F_top = np.array((0.,0.,0.))
    M_top = np.array((0.,0.,0.))
    M_lumped = None
#     M_lumped = np.zeros_like(pST['s_span'])
#     M_lumped[-1]= 1500000

    # Section loads
    # NOTE: only provide m_hydro if the acceleration was not provided and part of p_hydro
    m_hydro = np.zeros(len(pST['m']))
    #if pSS is not None:
    #    m_hydro = (pSS['rho'] * pHD['Ca'] * np.pi * pHD['D']**2 / 4).flatten()
    #    m_hydro[pST['z']>0] = 0 # Very important, TODO wave elevation
    F_sec, M_sec, outD =  beamSectionLoadsFromShapeFunctions(q, qd, qdd, p_ext, F_top, M_top, pST['s_span'], pST['PhiU'], pST['PhiV'], PhiK = pST['PhiK'], m=pST['m'],  a_ext=a_ext, debug=False,
                                                             m_hydro=m_hydro, M_lumped=M_lumped)
    F_sec = np.vstack((F_sec, M_sec))
    outD.update(outF)

    if not reconHydro:
        #outD['p_hydro_no_struct_acc'] = outF_no_Acc['p_hydro']
        p_no_acc = outF['p_FK']+ outF['p_AM_wav']+ outF['p_drag']
        outD['p_hydro_no_struct_acc'] = p_no_acc
    else:
        pass
        # p_hydro_no_struct_acc already returned in that case in 

    #if t>10:
    #   oo = outF_no_Acc
    #   PT = oo['p_hydro']
    #   import pdb; pdb.set_trace()
    return F_sec, outD



# --------------------------------------------------------------------------------}
# ---  
# --------------------------------------------------------------------------------{
def monopileSetupFromOpenFAST(fstFile, shapes_sub=[0,4], TMIN=0.0, TMAX=None, nSubSample=1, compFile=None, tuneM=False, hydroShape=None, reconHydro=False):

    # --- Loading a "YAMS" model
    WT = FASTmodel2MTNSB(fstFile, shapes_sub=shapes_sub, shapes_twr=[], shapes_bld=[], DEBUG=False, bStiffening=True, main_axis='z', fixedShaft=True).WT

    zBeam = np.array(sorted(WT.fnd.SD.pointsMN['z'].unique()))
    zDepth = WT.fnd.s_span - WT.WtrDpth
    zDepth = zBeam
    if not np.array_equal(zBeam, zDepth):
        raise Exception('Problem in model, both z are different')

    # --- Structure
    pST         = dict()
    pST['PhiU']   = WT.fnd.PhiU
    pST['PhiV']   = WT.fnd.PhiV
    pST['PhiK']   = WT.fnd.PhiK
    pST['m']      = WT.fnd.m
    pST['s_span'] = WT.fnd.s_span
    pST['gravity'] = WT.gravity

    pST['MM']   = WT.MM.copy()
    pST['KK']   = WT.KK.copy()
    pST['DD']   = WT.DD.copy()
    pST['z']    = WT.fnd.s_span - WT.WtrDpth
    pST['WT']   = WT # Big HACK
    # --- Sea state and hydrodynamics
    pSS = None
    pHD = None
    if WT.Hydro:
        # --- Sea state, wave kinematics
        pSS = dict()
        pSS['rho']        = WT.WtrDens      # [kg/m3]
        pSS['WaterDepth'] = WT.WtrDpth      # [m]
        if compFile is not None:
            dfComp = weio.read(compFile).toDataFrame()
            pSS['ap']   = dfComp['Amplitude_[m]']
            pSS['fp']   = dfComp['Frequency_[Hz]']
            pSS['epsp'] = dfComp['Phase_[rad]']
            # --- Wave Kinematics
            pSS['kp'] = wavenumber(pSS['fp'], pSS['WaterDepth'], pST['gravity']) # Wave numbers
            pSS['compFile'] = compFile
#         else:
#             print('>>> Using Wave of amplitude 3 and period 12 for now')
#             pSS['ap']   = np.array([3])     # Amplitudes
#             pSS['fp']   = np.array([1/12])  # frequencies [Hz]
#             pSS['epsp'] = np.array([np.pi]) # Deterministic phas


        # ---  Hydrodynamics 
        HD = WT.HD
        try:
            cprop = HD.getTab('SectionPropCyl')
        except:
            cprop = HD.getTab('SectionProp')
        D_ = cprop['PropD'].values[0]
        sprop = HD.getTab('SmplPropCyl')
        Cp_ = sprop['SimplCp'].values[0]
        Cd_ = sprop['SimplCd'].values[0]
        Ca_ = sprop['SimplCa'].values[0]
        print(f'HD props: Cp={Cp_} Cd={Cd_} Ca={Ca_} D={D_} ')
        PlaceHolderOnes = np.ones([len(zDepth),1])
        pHD = dict()
        pHD['zDepth'] = zDepth
        pHD['D']      = D_  * PlaceHolderOnes
        pHD['Cd']     = Cd_ * PlaceHolderOnes
        pHD['CM']     = (Ca_+Cp_) *PlaceHolderOnes
        pHD['Ca']     = Ca_ * PlaceHolderOnes
        pHD['Cp']     = Cp_ * PlaceHolderOnes
        pHD['m_hydro']= pSS['rho'] * np.pi / 4 * pHD['D'].ravel()**2 * pHD['Ca'].ravel()
            
        # Generalized hydro mass matrix 
        GM_hydro = np.zeros((len(WT.fnd.PhiU),len(WT.fnd.PhiU)))
        bWet = zDepth<=0
        for i, phi in enumerate(WT.fnd.PhiU):
            phi_x = phi[0,:]
            GM_hydro[i,i] = np.trapezoid(pHD['m_hydro'][bWet] * phi_x[bWet]**2, zDepth[bWet])

        pHD['GM_hydro']= GM_hydro
            
        print('GM_hydro:\n', GM_hydro)
        for i, phi in enumerate(WT.fnd.PhiU):
            # NOTE: diagonal only?
            WT.MM[i,i]+=GM_hydro[i,i]

    if tuneM:
        raise Exception('Removed ?')
        # TODO remove this, it's application specific
        # Tuning
        factM1=1.009
        factM2=1.003
        WT.MM[0,0]*=factM1
        WT.MM[1,1]*=factM2
        WT.KK[0,0]*=factM1
        WT.KK[1,1]*=factM2


    # --- Load hydroShape
    if hydroShape is not None:
        dfH = weio.read(hydroShape).toDataFrame()
        if not np.array_equal(dfH['z_[m]'],pHD['zDepth']):
            raise Exception('z depth different with shape function and hydrodyn')
        pHD['phi']  = dfH['phi_[Ns/m^2]'].values # Used to be called k_h_z
        pHD['phit'] = dfH['phit_[-]'].values
        bWet = zDepth<=0
        k_h = np.zeros(len(WT.fnd.PhiU))
        for i, phi in enumerate(WT.fnd.PhiU):
            phi_x = phi[0,:]
            k_h[i] = np.trapezoid(pHD['phi'][bWet] * phi_x[bWet], zDepth[bWet])
        print('k_h     : ', k_h)
        pHD['k_h'] = k_h


    # --- OPENFAST as a reference simulation
    outFile = fstFile.replace('.fst','.outb')
    if not os.path.exists(outFile):
        outFile = fstFile.replace('.fst','.out')
    if not os.path.exists(outFile):
        WARN('Cannot setup reference, out file does not exist: '+outFile)
        ref = None
    else:
        ref = dict()
        # --- Reading input loads
        df = weio.read(outFile).toDataFrame()

        import matplotlib.pyplot as plt
        if WT.Hydro and compFile is None:
            NOTE('Computing wave component from eta ref')
            pSS_comp = wave_components(df['Wave1Elev_[m]'], df['Time_[s]'], water_depth = pSS['WaterDepth'], g=pST['gravity'], eta_out=True, test=True, plot=True, aThreshold=1e-4)

        df = df.iloc[::nSubSample] # SubSampling for shorter comp time
        df = df[df['Time_[s]']<TMAX]   # Limiting time
        df = df[df['Time_[s]']>=TMIN]   # Limiting time
        vTime = df['Time_[s]'].values
        print('nSteps  : ',len(vTime))
        # --- State
        if shapes_sub==[0,4]:
            ref['state']   = df[ ['Q_Sg_[m]', 'Q_P_[rad]', 'QD_Sg_[m/s]', 'QD_P_[rad/s]']].values.T
            ref['state_d'] = df[ ['QD_Sg_[m/s]', 'QD_P_[rad/s]', 'QD2_Sg_[m/s^2]', 'QD2_P_[rad/s^2]']].values.T
        elif shapes_sub==[0]:
            ref['state']   = df[ ['Q_Sg_[m]',    'QD_Sg_[m/s]']].values.T
            ref['state_d'] = df[ ['QD_Sg_[m/s]', 'QD2_Sg_[m/s^2]']].values.T
        else:
            raise NotImplementedError()
        ## --- Nicknames
        ref['eta']         = df['Wave1Elev_[m]']
        ref['HydroFx']     = df['HydroFxi_[N]'].values
        # --- Section loads Ref
        # SubDyn secton outputs
        zBeam, F_sec, r_sec =  WT.fnd.SD.beamSecOutputs(df, verbose=False)
        ref['df']    = df
        ref['z']     = zBeam
        ref['F_sec'] = F_sec
        ref['r_sec'] = r_sec

    if ref is not None and WT.Hydro:
        #print('Wave freq {}  amplitude {} '.format(np.max(pSS['fp']),np.max(pSS['ap'])))
        if compFile is None:
            #pSS_ = wave_components(ref['eta'], vTime, water_depth = pSS['WaterDepth'], g=pST['gravity'], eta_out=True, test=False, plot=True, aThreshold=1e-5)
            #plt.show()
            #eta_sim = elevation2d(pSS['ap'], pSS['fp'], pSS['kp'], pSS['epsp'], vTime, x=0)
            pSS.update(pSS_comp)
            # We recompute wave elevation because time has changed.
            eta_sim = elevation2d(pSS['ap'], pSS['fp'], pSS['kp'], pSS['epsp'], vTime, x=0)

        else:
            eta_sim = elevation2d(pSS['ap'], pSS['fp'], pSS['kp'], pSS['epsp'], vTime, x=0)



        if len(vTime)>0:
            dt = vTime[1]-vTime[0]
            eta_dot = np.gradient(eta_sim, dt) # Velocity state
        else:
            print('[WARN] time is empty, eta dot is zero')
            eta_dot = eta_sim*0
        pSS['eta_time'] = vTime
        pSS['eta']      = eta_sim
        pSS['eta_dot']  = eta_dot

    #if reconHydro:
        # Safety
        #del pSS['ap']
        #del pSS['fp']
        #del pSS['kp']
        #del pSS['epsp']


    # --- Setup Sys
    Sys = MechSystem(WT.MM, WT.DD, WT.KK)
    Sys.setStateInitialConditions(WT.z0.values)
    # --- Set external loads time series
    ForceFunction = lambda t,q,qd : monopileGF(t, q, qd, pST=pST, pSS=pSS, pHD=pHD, reconHydro=reconHydro)[0]
    Sys.setForceFunction( ForceFunction )



    return pST, pSS, pHD, Sys, WT, ref


if __name__ == '__main__':
    SSCase='_Hs=2.5_Tp=10_h=50'  ; fstFile = f'../simulations/MT100/06_Jonswap/OF_Long{SSCase}.fst'; tRange=[0,600]; tRangeStats=[35,600]; 

    pST, pSS, pHD, Sys, WT, ref = monopileSetupFromOpenFAST(fstFile, shapes_sub=[0,4], TMIN=0.0, TMAX=None, nSubSample=1, compFile=None, tuneM=False, hydroShape=None, reconHydro=False)

    printMat(WT.MM    , var='MMa', digits=6)
    printMat(pST['MM'], var='MM' , digits=6)
    printMat(pST['KK'], var='KK' , digits=6)
    printMat(pST['DD'], var='DD' , digits=6)


