# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 16:38:48 2023

@author: ibaldoni
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def calculate_averages(data):
    """
    Takes a list of data points and returns the average maximum and the average minimum of
    each window of 5000 points in the data.
    """
    window_size = 5000
    num_windows = len(data) // window_size
    max_averages = []
    min_averages = []
    
    for i in range(num_windows):
        window_start = i * window_size
        window_end = window_start + window_size
        window_data = data[window_start:window_end]
        max_averages.append(np.max(window_data))
        min_averages.append(np.min(window_data))
    
    return np.mean(max_averages), np.mean(min_averages)

def plot_frequency_fluctuation(data, title, filename, frequency_division=30, time_slice=slice(None)):
    f2 = data.f2.astype(float) * frequency_division
    f2_mean = np.mean(f2)
    f2 = f2 - f2_mean
    f2 = f2*1e-3
    time = data.time.astype(float)
    time_slice = slice(None, time_slice.stop if time_slice.stop else onesecond)
    
    max_avg, min_avg = calculate_averages(f2)
    print(f"Average maximum: {max_avg}")
    print(f"Average minimum: {min_avg}")
    
    max_med, min_med = calculate_averages(f2)
    print(f"Median maximum: {max_med}")
    print(f"Median minimum: {min_med}")
    

    plt.figure()
    plt.plot(time[time_slice], f2[time_slice],color='green',label='Data')
    plt.xlabel('Time [s]')
    plt.ylabel('Frequency [kHz]')
    plt.hlines(max_med,time[time_slice][0],time[time_slice].iloc[-1],color='red',linestyle='dashed',label=f'±{max_avg:.3}kHz')
    plt.hlines(min_med,time[time_slice][0],time[time_slice].iloc[-1],color='red',linestyle='dashed')
    plt.title(title)
    plt.grid()
    plt.legend()
    # plt.savefig(filename, dpi=300)
    
    

if __name__ == '__main__':
        
    folder = r'C:\Users\ibaldoni\Documents\In Github\DiverseSimulationsTools\18 - Zeiss\Spatial resolution measurement\\'
    Data = ['With_PDH_lock.dat','With_PDH_lock_v1.dat']
    
    onesecond = 500000 #large number for getting all the data
    
    for fileName in Data:
    
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
            plot_frequency_fluctuation(med, 'Laser frequency fluctuation with PDH lock', 
                                       'Laser_w_lock_1s.png', time_slice=slice(onesecond))
        