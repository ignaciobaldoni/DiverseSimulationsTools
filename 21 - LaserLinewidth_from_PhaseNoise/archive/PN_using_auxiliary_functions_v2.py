# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 12:49:46 2023

@author: ibaldoni
"""

# import sys
# # Add the directory containing the script to the module search path
# sys.path.append('/path/to/directory')
# # Import all functions from the script
# from script_name import *
# # Call the functions as usual
# result = function_name(argument)

import pandas as pd
from Plot_aux_functions import *
import matplotlib.pyplot as plt
Plot_parameters('default')

import numpy as np
import scipy.integrate as integrate
from scipy.interpolate import interp1d
from scipy.signal import gaussian

fileName = '20210419-142704_Data_PSD.txt'
headers = ['Frequency','PSD']
phase_noise = pd.read_table(fileName,sep=';',names=headers,skiprows=1)


PSD = phase_noise.PSD
Frequency = phase_noise.Frequency
y_filtered = PSD

y_filtered[-500:]=-91

def apply_filter():
    # Define the width of the Gaussian filter in Hz (you may need to adjust this value based on your data)
    sigma = 1
    # Define the length of the filter in number of data points
    filter_length = int(0.5 * sigma / (Frequency[1] - Frequency[0]))
    # Create the Gaussian filter
    gaussian_filter = gaussian(filter_length, sigma / (Frequency[1] - Frequency[0]))
    # Normalize the filter so that the sum of its values is 1
    gaussian_filter /= np.sum(gaussian_filter)
    # Apply the filter to the data
    y_filtered = np.convolve(PSD, gaussian_filter, mode='same')



# Create an interpolation function
interp_func = interp1d(Frequency, y_filtered, kind='linear',fill_value='extrapolate')

# Evaluate the interpolated function over the same frequency range as the original data
interp_psd = interp_func(Frequency)

# # Plot the results
# fig, ax = plt.subplots()
# ax.semilogx(Frequency, PSD, label='Original PSD')
# ax.semilogx(Frequency, y_filtered, label='Filtered PSD')
# ax.semilogx(Frequency, interp_psd, label='Interpolated function')
# ax.set_xlabel('Frequency [Hz]')
# ax.set_ylabel('PSD [dBc/Hz]')
# ax.legend()

def PN_to_FN(frequencies):
    L = interp_func(frequencies)
    
    S_phi = 2 * 10**(L/10)
    S_df_f = S_phi * (2 * np.pi * frequencies)**2
    
    plt.plot(frequencies,S_df_f)

plt.figure
PN_to_FN(Frequency)
    
    


def Integral(frequencies):    
    L = interp_func(frequencies)
    
    S_phi = 2 * 10**(L/10)
    S_df_f = S_phi * (2 * np.pi * frequencies)**2
    
    return np.heaviside(S_df_f - 8*np.log(2)*frequencies/(np.pi**2),1)*S_df_f


def solveIntegral():
    return integrate.quad(lambda f: Integral(f), -0.1, 1e6) 

# plt.figure()
plt.semilogx(Frequency,Integral(Frequency),Frequency,8*np.log(2)*Frequency/(np.pi**2))





A = solveIntegral()[0]
# print(A)
FWHM = np.sqrt(8*np.log(2)*A)
print(FWHM,'Hz')
