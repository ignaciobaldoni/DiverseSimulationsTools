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

# phase_noise = phase_noise[:len(phase_noise)-250]



PSD = phase_noise.PSD
Frequency = phase_noise.Frequency
y_filtered = PSD



# # Define the width of the Gaussian filter in Hz (you may need to adjust this value based on your data)
# sigma = 50
# # Define the length of the filter in number of data points
# filter_length = int(5 * sigma / (Frequency[1] - Frequency[0]))
# # Create the Gaussian filter
# gaussian_filter = gaussian(filter_length, sigma / (Frequency[1] - Frequency[0]))
# # Normalize the filter so that the sum of its values is 1
# gaussian_filter /= np.sum(gaussian_filter)
# # Apply the filter to the data
# y_filtered = np.convolve(PSD, gaussian_filter, mode='same')


# Create an interpolation function
interp_func = interp1d(Frequency, y_filtered, kind='cubic',fill_value='extrapolate')

# Evaluate the interpolated function over the same frequency range as the original data
interp_psd = interp_func(Frequency)

# Plot the results
fig, ax = plt.subplots()
ax.semilogx(Frequency, PSD, label='Original PSD')
ax.semilogx(Frequency, y_filtered, label='Filtered PSD')
ax.semilogx(Frequency, interp_psd,'.-', label='Interpolated function')
ax.set_xlabel('Frequency [Hz]')
ax.set_ylabel('PSD [dBc/Hz]')
ax.legend()


# What we measure is: L(f) = 10*np.log10(0.5*PhaseNoise(f)) = PSD
    
# F(f) = 10^(L(f)/10) / (2 * pi * f)^2

    



def from_Lf_to_PhaseNoise(f):
    psd = interp_func(f)   
    PhaseNoise = 2*10**(psd/10)

    FrequencyNoise = f**2*PhaseNoise
    
    plt.figure()
    plt.loglog(f,FrequencyNoise)

    return PhaseNoise

from_Lf_to_PhaseNoise(Frequency) 
def Integral(f):    
    PhaseNoise = 1#from_Lf_to_PhaseNoise(f)
    FrequencyNoise = f**2*PhaseNoise
    
    return np.heaviside(FrequencyNoise - 8*np.log(2)*f/(np.pi**2),1)*FrequencyNoise 



def solveIntegral():
    return integrate.quad(lambda f: Integral(f), 0, 1e4) #np.inf


A = solveIntegral()[0]

print(A)

FWHM = np.sqrt(8*np.log(2)*A)

print(FWHM*1e-6,'MHz')





# # What we measure is: L(f) = 10*np.log10(0.5*PhaseNoise(f)) = PSD

# def fNoise(PSD,f):
#     PhaseNoise = 2*10**(PSD/10)
#     FrequencyNoise = 1
#     # w = 2*np.pi*f
#     # FrequencyNoise = 1/(w**2) * PhaseNoise
    
#     return FrequencyNoise
 
# fNoise(interp_psd,Frequency)


# def Integral(f):
#     FrequencyNoise = fNoise(interp_func(f),f)
#     print(np.heaviside(FrequencyNoise - 8*np.log(2)*f/(np.pi**2),1)*FrequencyNoise )
#     return np.heaviside(FrequencyNoise -\
#                         8*np.log(2)*f/(np.pi**2),1)*FrequencyNoise 

# def solveIntegral():
#     return integrate.quad(lambda f: Integral(f), 0, 1e6) #np.inf

# t0 = 0
# A = solveIntegral()[0]

# print(A)

# FWHM = np.sqrt(8*np.log(2)*A)

# print(FWHM,'Hz')
