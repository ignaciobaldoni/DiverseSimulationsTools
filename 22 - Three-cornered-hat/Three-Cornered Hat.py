# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 15:13:22 2023

@author: ibaldoni
"""


import pandas as pd
import numpy as np
import sys
sys.path.append(r'//menloserver/MFS/99-Data_Warehouse/02-User_Folders-Public/i.baldoni/python Util functions')
from Plot_aux_functions import Plot,Plot_parameters # Import all functions from the script
import matplotlib.pyplot as plt
# Plot_parameters('seaborn-dark')


def getValuesfromCSV(phase_noise):
    phase_noise = phase_noise.drop(['phase_noise_raw','3','4','5','6'], 1)
    phase_noise = phase_noise.reset_index()
    try:
        phase_noise = phase_noise.drop(['level_0','level_1'], 1)
    except:
        pass
    psd   = phase_noise.iloc[::2][:-1]
    frequency = phase_noise.iloc[1::2]
    frequency = np.asarray(frequency).squeeze()
    psd = np.asarray(psd).squeeze()
    
    print(len(psd),len(frequency))
    
    return frequency, psd


def read_phase_noise_data(fileName):
    """Reads phase noise data from a file and returns a pandas dataframe."""
    headers = ['Frequency', 'PSD']
    if fileName[-4:] == '.csv': 
        headers = ['phase_noise_raw','phase_noise_avg','3','4','5','6']        
        phase_noise = pd.read_csv(fileName, sep=',', names=headers, skiprows=59)
        # print(phase_noise.head())
        frequency, psd =  getValuesfromCSV(phase_noise)
        
    if fileName[-4:] == '.txt': 
        phase_noise = pd.read_table(fileName, sep=';', names=headers, skiprows=1)
    
        frequency = phase_noise.Frequency[phase_noise.PSD.values<-1]
        psd = phase_noise.PSD[phase_noise.PSD.values<-1]
        
        
    return frequency, psd


# Data readout for three-cornered hat measurements
file = ["", "", ""]
file[0] = "20210419-142704_Data_PSD.txt"
file[1] = "20210419-142704_Data_PSD.txt"#"//menloserver/MFS/03-Operations/02-DCP/03-Entwicklungsprojekte/9556-COSMIC/52-Messergebnisse/20210510_PhaseNoiseRIOwithUSC/MeasR_0001.csv"
file[2] = "20210419-142704_Data_PSD.txt"#"//menloserver/MFS/03-Operations/02-DCP/03-Entwicklungsprojekte/9556-COSMIC/52-Messergebnisse/20210504_Rio_Laser_PhaseNoise/Trace_0001.csv"


for i in range(0,3):
    print(i)
    
    frequency,psd = read_phase_noise_data(file[i])
    psd = psd-i*1e2/frequency
    
    # print(frequency.head(),psd.head())
    plt.figure(1)
    Plot(frequency,psd, num=1,log_scale=True,label=f'file{i}')
    # plt.show()
    
    if i==0: psd1 = psd
    if i==1: psd2 = psd
    if i==2: psd3 = psd

# Calculate the 3CH formula
delta_psd12 = np.abs(psd1 - psd2)
delta_psd23 = np.abs(psd2 - psd3)
delta_psd31 = np.abs(psd3 - psd1)
abs_psd1 = np.sqrt((delta_psd12**2 + delta_psd31**2 - delta_psd23**2)/2)
abs_psd2 = np.sqrt((delta_psd12**2 + delta_psd23**2 - delta_psd31**2)/2)
abs_psd3 = np.sqrt((delta_psd23**2 + delta_psd31**2 - delta_psd12**2)/2)

Plot(frequency,delta_psd12, num=2,log_scale=True)
Plot(frequency,delta_psd23, num=2,log_scale=True)
Plot(frequency,delta_psd31, num=2,log_scale=True)
