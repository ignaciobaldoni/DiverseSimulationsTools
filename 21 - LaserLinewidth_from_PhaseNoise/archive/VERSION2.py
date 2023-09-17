# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 23:46:00 2023

@author: ibaldoni
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 22:52:34 2023

@author: ibaldoni
"""

import sys
sys.path.append(r'\\menloserver\MFS\99-Data_Warehouse\01-User_Folders-Private\i.baldoni\2023 - Ultrastable Microwaves\03 - Simulations')
from Plot_aux_functions import Plot,Plot_parameters # Import all functions from the script
# import matplotlib.pyplot as plt
Plot_parameters('seaborn-dark')

import pandas as pd
import numpy as np
import scipy.integrate as integrate
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
# from scipy.signal import gaussian

dB_3 = 3 #To ensure double sideband

def read_phase_noise_data(fileName, filetype = '.txt'):
    """Reads phase noise data from a file and returns a pandas dataframe."""
    headers = ['Frequency', 'PSD']
    if filetype == '.csv': 
        headers = ['phase_noise_raw','phase_noise_avg','3','4','5','6']        
        phase_noise = pd.read_csv(fileName, sep=',', names=headers, skiprows=59)
        
    if filetype == '.txt': 
        phase_noise = pd.read_table(fileName, sep=';', names=headers, skiprows=1)
        
        
    return phase_noise

def filter_psd(phase_noise):
    """Filters the PSD values in the phase_noise dataframe."""
    y_filtered = psd #phase_noise.PSD.values

    # y_filtered[-1500:] = -85
    return y_filtered

def create_interp_func(frequency, y_filtered):
    """Creates an interpolation function from the frequency and filtered PSD data."""

    interp_func = interp1d(frequency, y_filtered, kind='linear', fill_value='extrapolate')
    return interp_func

def calculate_S_df_f(L, frequencies):
    """Calculates S_df_f from the PSD values."""
    S_phi = 2 * 10**((L+dB_3)/10)
    S_df_f = S_phi * (frequencies)**2 # Does 2pi goes here?

    return S_df_f


def calculate_integral(frequencies, interp_func):
    """Calculates the integral using the interpolated PSD values."""
    
    L = interp_func(frequencies)
    # plt.plot(frequencies, L)
    S_df_f = calculate_S_df_f(L, frequencies)
    threshold = 8 * np.log(2) * frequencies / (np.pi**2)
    integrand = np.heaviside(S_df_f - threshold, 1) * S_df_f
    
    return integrate.simps(integrand, frequencies), threshold, S_df_f

def calculate_FWHM(integral):
    """Calculates the FWHM from the integral value."""
    return np.sqrt(8 * np.log(2) * integral)

def getValuesfromCSV(phase_noise):
    phase_noise = phase_noise.drop(['phase_noise_raw','3','4','5','6'], 1)
    phase_noise = phase_noise.reset_index()
    phase_noise = phase_noise.drop(['level_0','level_1'], 1)
    psd   = phase_noise.iloc[::2][:-1]
    frequency = phase_noise.iloc[1::2]
    frequency = np.asarray(frequency).squeeze()
    psd = np.asarray(psd).squeeze()
    
    return frequency, psd




if __name__ == '__main__':
    fileName = '20210419-142704_Data_PSD.txt'
    
    # fileName = "//menloserver/MFS/03-Operations/02-DCP/03-Entwicklungsprojekte/9556-COSMIC/52-Messergebnisse/20210510_PhaseNoiseRIOwithUSC/MeasR_0001.csv"
  

    phase_noise = read_phase_noise_data(fileName,filetype=fileName[-4:])
    

    if '.csv' in fileName:
        print('Estamos bien')
        frequency, psd =  getValuesfromCSV(phase_noise)
        
    if '.txt' in fileName:
        frequency = phase_noise.Frequency[phase_noise.PSD.values<-1]
        psd = phase_noise.PSD[phase_noise.PSD.values<-1]
           
    
    y_filtered = filter_psd(psd) 

    interp_func = create_interp_func(frequency, y_filtered)
    integral, threshold, S_df_f = calculate_integral(frequency, interp_func)
    FWHM = calculate_FWHM(integral)

    print(FWHM, 'Hz')
    
    plt.loglog(frequency, S_df_f)
    plt.plot(frequency, threshold)
    
    Plot(frequency,S_df_f,num=1,label='S$_{\delta f}(f)$',yLabel='PSD [Hz$^2$/Hz]')
    Plot(frequency,threshold,num=1,label='beta line',yLabel='PSD [Hz$^2$/Hz]')

    
    Plot(frequency,y_filtered,log_scale=True, yLabel='PSD [dBc/Hz]',num=2, 
          label = f'Phase Noise from ORS\nLinewidth = {FWHM:.4} Hz')
    
    
    
    
