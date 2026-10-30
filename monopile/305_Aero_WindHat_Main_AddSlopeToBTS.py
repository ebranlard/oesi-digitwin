"""
Creates a turbulence time series with a step of wind

"""

import matplotlib.pyplot as plt
import numpy as np
from welib.weio.turbsim_file import TurbSimFile


if __name__ == '__main__':
    btsFile = 'simulations_wt/Wind/TurbSimTS_IEA22.bts'
    btsOut  = btsFile.replace('_HatGen.bts', '_.bts')

    bts = TurbSimFile(btsFile)
    print(bts)
    dfs= bts.toDataFrame()

    u_m = dfs['ZMidLine']['u_[m/s]'].mean()
    t_new  = bts.t
    Tmax = t_new[-1]+bts.dt


    #u0 = [0,4  ,14 ,12, 0  ,0]
    #t0 = [0,150,400,550,570, Tmax]
    u0 = np.array([0 , 0   ,  8 , 9   , 8   , 2   , 0])-1
    t0 =          [0 , 150 , 250, 400 , 550 , 570 , Tmax]
    u_new = np.interp(t_new, t0, u0)

    bts.superimposeTimeSeries(t_new, u_new)
    print(bts)
    dfs= bts.toDataFrame()
    bts.write(btsOut)

    fig, ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8))
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    ax.plot(dfs['ZMidLine']['t_[s]'], dfs['ZMidLine']['u_[m/s]'])
    ax.plot(t_new, u_new , label='')
    ax.plot(t_new, u_new+u_m , label='')
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.legend()

    plt.show()
