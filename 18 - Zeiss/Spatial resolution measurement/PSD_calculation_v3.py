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

# Color = 'g'
Color = '#053061'
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
time_interval = time.iloc[1]-time.iloc[0]
#(time.iloc[-1] - time.iloc[0]) / (num_samples - 1)
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


fffluctu=mplm.csd(beat_norm, beat_norm, Fs=sampling_freq,
                  detrend='mean',scale_by_freq='True',NFFT=len(time)) # result in 1/Hz


S_f = fffluctu[0].real

sqrtfflu=np.sqrt(S_f*f2_mean**2) # conversion to Hz²/Hz then Hz/sqrt(Hz)
ffreq = fffluctu[1] # Fourier frequencies

time_cst=int(1/sampling_freq*1000)
pn1, = plt.loglog(ffreq,sqrtfflu,label="FXE %i ms data" %time_cst, color = Color)
ax[1].set_ylabel('Frequency noise [Hz/\u221AHz]')
ax[1].legend()
ax[1].set_ylim([1e2,10e8])
ax[1].set_xlabel('Frequency Offset [Hz]')
add_grids(ax[1])
if savefigs: plt.savefig('PlotsCombined.png')

plt.figure()
phase_noise_dBc_Hz = 10*np.log10(S_f) - 10 * np.log10(ffreq ** 2)
plt.semilogx(ffreq,phase_noise_dBc_Hz, color = Color)
plt.ylim([-150,40])
plt.ylabel(r'$\mathcal{L}$ Phase Noise [dBc/Hz]')
plt.xlabel('Frequency Offset [Hz]')
add_grids()
if savefigs: plt.savefig('PhaseNoise.png')


##########       Modified allan deviation calculation       ###########

sys.path.append(r'\\menloserver\mfs\99-Data_Warehouse\01-User_Folders-Private\i.baldoni\05 - Data processing\00 - Allan deviation collection\14-Compute ADEV, Phase Noise, etc from FXE data')
from AllanDeviationFunctions import calcAllanDev
 
wavelength = 1542.14e-9   

# f0 = (speed_of_light/wavelength) 
f0 = 0*f2_mean

# Define the acquisition frequency in Hz (frequency rate)
freq_rate = 1000 

#### m = Modified || h = Hadamard || a = Allan deviation || oa = Overlapping
AllanType = 'm'

calcAllanDev(frequency_fluctuations, 
             channels=[1], 
             channel_names=['SMILE - ADAM'], 
             freq_rate = freq_rate, 
             TCH = False,
             f_0 = f0,
             Dedrift = False, 
             wavelength = wavelength, 
             Type = AllanType, 
             saveFigs = False,
             Taus = [0.01, 0.10, 1.0, 10.0], 
             normalized_data = False)

plt.ylim([10e-14,10e-5])



from AllanDeviationFunctions import PSD_Noise_from_Counter


PSD_Noise_from_Counter(beat_norm,
                       channels=[4],
                       channel_names = ['SMILE - ADAM'],
                           wavelength = 1542.14e-9, 
                           freq_rate = 1e3, f_0 = f2_mean,
                           Requirements = False, 
                           normalized_data= True,
                           Dedrift = False)