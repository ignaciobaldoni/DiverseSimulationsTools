# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 09:55:21 2024

@author: ibaldoni
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import matplotlib.mlab as mplm
import warnings
warnings.filterwarnings("ignore")

sys.path.append(r'\\menloserver\MFS\99-Data_Warehouse\02-User_Folders-Public\i.baldoni\python Util functions')
from Plot_aux_functions import Plot_parameters, add_grids
Plot_parameters(width=10)

folder = r'C:\Users\ibaldoni\Documents\In Github\DiverseSimulationsTools\18 - Zeiss\Spatial resolution measurement\\'
File = r'With_PDH_lock.dat'
# File = r'With_PDH_lock_v1.dat'
# File = r'Without_PDH_lock.dat'

Color = 'g'
savefigs = False


data = pd.read_table(folder+File, 
                      sep='\t', 
                      skiprows=28, 
                      names=['time', 'f1', 'u1', 'f2', 'U2', 'U1', 'u2', 'uawg'], 
                      encoding='unicode_escape')

frequency_division = 30
frequency_fluctuations = data.f2.astype(float) * frequency_division
f2_mean = np.mean(frequency_fluctuations)
frequency_fluctuations = frequency_fluctuations - f2_mean


time = data.time.astype(float)
num_samples = len(frequency_fluctuations)
time_interval = time.iloc[1]-time.iloc[0]#(time.iloc[-1] - time.iloc[0]) / (num_samples - 1)
# Sampling frequency is the reciprocal of the time interval
sampling_freq = 1 / time_interval

Max_freq_fluct = np.max(frequency_fluctuations)

beat_norm = frequency_fluctuations/Max_freq_fluct

print('Time interval:\t\t',time_interval*1e3,'ms')
print('Sampling frequency:\t', sampling_freq*1e-3,'kHz')

fig, ax = plt.subplots(2,1)
ax[0].plot(time,frequency_fluctuations*1e-3,color=Color)
ax[0].set_ylabel('Freq. fluct. [kHz]')
ax[0].set_xlabel('Time [s]')
add_grids(ax[0])

# Compute Power Spectral Density (PSD)
psd, frequencies = ax[1].psd(beat_norm,#frequency_fluctuations, 
                              NFFT=len(time),
                              Fs=sampling_freq,color=Color) 
                              #sides = 'twosided', 
                              #window=np.hanning(num_samples))
                             
if savefigs: plt.savefig('FrequencyFluctuation_and_PSD.png')

# # Convert PSD to dB/Hz
# psd_dB = 10 * np.log10(psd)
# # Compute phase noise in dBc/Hz
# phase_noise_dBc_Hz = psd_dB - 10 * np.log10(frequencies ** 2)

# # Plot PSD and phase noise
# plot_combined = False
# if plot_combined:
#     fig, ax = plt.subplots(2,1)
#     ax[0].loglog(frequencies, np.sqrt(psd.real*Max_freq_fluct**2),color=Color)
#     ax[0].set_ylabel('Frequency noise [Hz/\u221AHz]')
#     add_grids(ax[0])
#     ax[0].set_ylim([1e-1,10e5])

    
#     # Plot Phase Noise
#     ax[1].semilogx(frequencies, phase_noise_dBc_Hz,color=Color)
#     ax[1].set_xlabel('Frequency Offset [Hz]')
#     ax[1].set_ylabel(r'$\mathcal{L}$ [dBc/Hz]')
#     add_grids(ax[1])
#     ax[1].set_ylim([-150,40])

    # if savefigs: plt.savefig('PlotsCombined.png')
# else:
    
#     plt.figure()
#     plt.loglog(frequencies, np.sqrt(psd.real*Max_freq_fluct**2),color=Color)
#     plt.ylabel('Frequency noise [Hz/\u221AHz]')
#     plt.xlabel('Frequency Offset [Hz]')
#     add_grids()
#     plt.ylim([1e-1,10e5])
#     if savefigs: plt.savefig('FrequencyNoise.png')
    
#     # Plot Phase Noise
#     plt.figure()
#     plt.semilogx(frequencies, phase_noise_dBc_Hz,color=Color)
#     plt.xlabel('Frequency Offset [Hz]')
#     plt.ylabel(r'Phase Noise $\mathcal{L}$ [dBc/Hz]')
#     add_grids()
#     plt.ylim([-150,40])
#     if savefigs: plt.savefig('PhaseNoise.png')



