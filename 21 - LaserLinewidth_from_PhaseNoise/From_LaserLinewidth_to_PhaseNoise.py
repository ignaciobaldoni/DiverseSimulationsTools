# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 22:52:34 2023

@author: ibaldoni
"""

#%% Personal functions 
import sys
sys.path.append(r'\\menloserver\MFS\99-Data_Warehouse\01-User_Folders-Private\i.baldoni\python Util functions')
from Plot_aux_functions import Plot,Plot_parameters, add_grids # Import all functions from the script
from util_Functions import units
# import matplotlib.pyplot as plt
Plot_parameters('seaborn-dark')

#%% import python modules
import pandas as pd
import numpy as np
import scipy.integrate as integrate
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt


dB_3 = 0 # = 3 To ensure double sideband



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
    S_df_f = calculate_S_df_f(L, frequencies)
    threshold = 8 * np.log(2) * frequencies / (np.pi**2)
    integrand = np.heaviside(S_df_f - threshold, 1) * S_df_f
    
    return integrate.simps(integrand, frequencies), threshold, S_df_f

def calculate_FWHM(integral):
    """Calculates the FWHM from the integral value."""
    return np.sqrt(8 * np.log(2) * integral)


def read_phase_noise_data(fileName):
    """Reads phase noise data from a file and returns a pandas dataframe."""
    file_extension = fileName.split('.')[-1]
    
    if file_extension == 'csv':
        print('Reading CSV file...')
        phase_noise = ''

        
    elif file_extension == 'txt':
        print('Reading TXT file...')
        headers = ['Frequency', 'PSD']
        phase_noise = pd.read_table(fileName, sep=';', names=headers, skiprows=1)
    else:
        print('Error: Invalid file extension.')
        return None
    
    return phase_noise

def getValuesfromCSV(fileName, sep=','):
    headers = ['Frequency','PSD']        
    phase_noise = pd.read_csv(fileName, sep=sep, names=headers, skiprows=61)
    psd = phase_noise['PSD']
    frequency = phase_noise['Frequency']
    
    return frequency, psd


    
    

if __name__ == '__main__':
    fileName = '20210419-142704_Data_PSD.txt'
    # fileName = 'pn1avg.csv'
    
    phase_noise = read_phase_noise_data(fileName)

    
    if phase_noise is not None:
        if 'csv' in fileName:
            frequency, psd =  getValuesfromCSV(fileName)
        elif 'txt' in fileName:
            frequency = phase_noise.loc[phase_noise['PSD'] < -1, 'Frequency']
            psd = phase_noise.loc[phase_noise['PSD'] < -1, 'PSD']
        else:
            print('Error: Invalid file extension.')
            frequency, psd = None, None

    
    psd = -80+np.linspace(-80,0,len(frequency))**-1
        
        
        
   
    interp_func = create_interp_func(frequency, psd)
    integral, threshold, S_df_f = calculate_integral(frequency, interp_func)
    FWHM = calculate_FWHM(integral)
    print(FWHM)
    
    FWHM_unit,factor = units(FWHM)
    print(FWHM/factor, str(FWHM_unit)+'Hz')
    

    Plot(frequency,S_df_f,num=1,label='S$_{\delta f}(f)$',yLabel='PSD [Hz$^2$/Hz]', loglog=True)
    add_grids()
    Plot(frequency,threshold,num=1,label='beta line',yLabel='PSD [Hz$^2$/Hz]')    
    add_grids()
    # Plot(frequency,psd,log_scale=True, yLabel='PSD [dBc/Hz]',num=2, 
    #       label = f'Phase Noise data\nLinewidth = {FWHM:.4} Hz')

    Plot(frequency,psd,log_scale=True, yLabel='PSD [dBc/Hz]',num=2, 
          label = f'Phase Noise data\nLinewidth = {FWHM/factor:.4} {FWHM_unit}Hz')
    
    add_grids()    
