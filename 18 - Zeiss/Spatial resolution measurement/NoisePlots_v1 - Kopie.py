# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 16:38:48 2023

@author: ibaldoni
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def plot_frequency_fluctuation(data, title, filename, frequency_division=30, time_slice=slice(None)):
    f2 = data.f2.astype(float) * frequency_division
    f2_mean = np.mean(f2)
    f2 = f2 - f2_mean
    time = data.time.astype(float)
    time_slice = slice(None, time_slice.stop if time_slice.stop else onesecond)
    plt.figure()
    plt.plot(time[time_slice], f2[time_slice],color='green')
    plt.xlabel('Time [s]')
    plt.ylabel('Frequency [Hz]')
    plt.title(title)
    plt.grid()
    plt.savefig(filename, dpi=300)

    

if __name__ == '__main__':
        
    folder = r'C:\Users\ibaldoni\Documents\In Github\DiverseSimulationsTools\18 - Zeiss\Spatial resolution measurement\\'
    Data = ['With_PDH_lock.dat','With_PDH_lock_v1.dat', 'Without_PDH_lock.dat']
    
    onesecond = 4200 #large number for getting all the data
    
    for fileName in Data:
        
        Label = ['Frequency fluctuation with PDH lock' if 'Without' not in fileName else 'Frequency fluctuation without PDH lock'][0]
    
        with open(folder+fileName) as f:
            med = pd.read_table(f, sep='\t', skiprows=28, names=['time', 'f1', 'u1', 'f2', 'U2', 'U1', 'u2', 'uawg'])
        
    
        # Plot Parameters
        with plt.style.context({
            'xtick.labelsize': 15,
            'ytick.labelsize': 15,
            'font.size': 15,
            'figure.figsize': (8,6),
            'grid.alpha': 0.75
        }):
            plot_frequency_fluctuation(med, Label, 
                                        f'{Label}_1s.png', time_slice=slice(onesecond))
        