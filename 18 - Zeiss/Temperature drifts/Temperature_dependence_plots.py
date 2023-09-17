# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 16:38:48 2023

@author: ibaldoni
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def plot_frequency_fluctuation(data, title, filename, frequency_division=30, time_slice=slice(None),Savefig=True):
    f2 = data.f2.astype(float) * frequency_division   
    f2_mean = np.mean(f2)
    f2 = f2 - f2[0]
    f2 = f2*1e-3*.960
    
    u2 = data.u2.astype(float)
    time = data.time.astype(float)
    # time_slice = slice(None, time_slice.stop if time_slice.stop else onesecond)

    fig, ax = plt.subplots(2)
    ax[0].plot(time, u2,color='blue',label='Voltage')
    ax[0].set_ylabel('Voltage [V]')
    ax[0].grid()
    # ax[0].set_xlabel('Time [s]')
    # ax[0].legend()
    # ax[0].set_title(title)
    
    ax[1].plot(time, f2,color='green',label='Data')
    ax[1].set_xlabel('Time [s]')
    ax[1].set_ylabel('Length variation [pm]')
    ax[1].grid()
    # ax[1].set_title(title)
    # ax[1].legend()
    if Savefig: plt.savefig(filename, dpi=300)



Folder = 'C:\\Users\\ibaldoni\\Desktop\\Projekte und Verknüpfungen\\01 - Zeiss Projekte\\1 - ADAM-Carla-Smile\\1 - Measurements\\Results2\\'

Option = [  'raw-02-0002_Hysterese_Bias5A5f0_1-UI.dat',
            'raw-01-0001_PolarisationUeff0UPeak10-UI.dat',
            'raw-08-0008_StepsA_140s-UI.dat',
            'raw-047-0002_Hyst_Bias5A5f1-UI.dat',
            'raw-214-0040_Hyst_Bias10A5f40-UI.dat']


p=0

for option in Option:
    with open(Folder+option) as f:
        med = pd.read_table(f, sep='\t', skiprows=50, names=['time', 'f1', 'u1', 'f2', 'U2', 'U1', 'u2', 'uawg'])

    
    plt.style.use('seaborn-dark')
    # Plot Parameters
    with plt.style.context({
        'xtick.labelsize': 15,
        'ytick.labelsize': 15,
        'font.size': 15,
        'figure.figsize': (8,6),
        'grid.alpha': 0.75,
        'figure.constrained_layout.use':True
    }):
        plot_frequency_fluctuation(med, '', 'Temperature_dependence_'+str(p)+'.png',Savefig = False)
    p+=1
