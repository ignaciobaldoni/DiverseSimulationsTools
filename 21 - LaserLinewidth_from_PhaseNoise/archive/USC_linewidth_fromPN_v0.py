# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 22:52:34 2023

@author: ibaldoni
"""

import sys
sys.path.append(r'\\menloserver\MFS\99-Data_Warehouse\01-User_Folders-Private\i.baldoni\2023 - Ultrastable Microwaves\3 - Simulations')
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

def read_phase_noise_data(fileName):
    """Reads phase noise data from a file and returns a pandas dataframe."""
    headers = ['Frequency', 'PSD']
    phase_noise = pd.read_table(fileName, sep=';', names=headers, skiprows=1)
    return phase_noise

def filter_psd(phase_noise):
    """Filters the PSD values in the phase_noise dataframe."""
    y_filtered = phase_noise.PSD.values
    # y_filtered[-1500:] = -85
    return y_filtered

def create_interp_func(frequency, y_filtered):
    """Creates an interpolation function from the frequency and filtered PSD data."""
    interp_func = interp1d(frequency, y_filtered, kind='linear', fill_value='extrapolate')
    return interp_func

def calculate_S_df_f(L, frequencies):
    """Calculates S_df_f from the PSD values."""
    S_phi = 2 * 10**(L+dB_3/10)
    S_df_f = S_phi * (frequencies)**2 # Does 2pi goes here?
    return S_df_f


def calculate_integral(frequencies, interp_func):
    """Calculates the integral using the interpolated PSD values."""
    
    L = interp_func(frequencies)
    S_df_f = calculate_S_df_f(L, frequencies)
    threshold = 8 * np.log(2) * frequencies / (np.pi**2)
    integrand = np.heaviside(S_df_f - threshold, 1) * S_df_f
    return integrate.simps(integrand, frequencies)

def calculate_FWHM(integral):
    """Calculates the FWHM from the integral value."""
    return np.sqrt(8 * np.log(2) * integral)

if __name__ == '__main__':
    fileName = '20210419-142704_Data_PSD.txt'
    
    
    phase_noise = read_phase_noise_data(fileName)
    
    # phase_noise = phase_noise
    phase_noise = phase_noise[phase_noise.PSD.values<-1]
    print(phase_noise.head())
    frequency = phase_noise.Frequency.values 
    y_filtered = filter_psd(phase_noise)
    interp_func = create_interp_func(frequency, y_filtered)
    integral = calculate_integral(frequency, interp_func)
    FWHM = calculate_FWHM(integral)

    print(FWHM, 'Hz')
    
    Plot(frequency,y_filtered,log_scale=True, yLabel='PSD [dBc/Hz]', 
         label = f'Phase Noise from ORS\nLinewidth = {FWHM:.4} Hz')
    
    plt.plot()
    # Plot(frequency, interp_func(frequency),log_scale=True, yLabel = 'Interpolation')
