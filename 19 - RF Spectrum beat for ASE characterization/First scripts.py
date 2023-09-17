# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 13:07:19 2022

@author: ibaldoni
"""

import numpy as np
import matplotlib.pyplot as plt

E0 = 1
signal1Frequency = 7


samplingFrequency   = 10000;

# samplingInterval       = 1 / samplingFrequency; 
# print(samplingInterval)

beginTime           = 0; 
endTime             = 1; 

# t = np.arange(beginTime, endTime, samplingInterval);
t = np.linspace(beginTime, endTime,samplingFrequency)

noise = 0.
E1 = np.sin(2*np.pi*signal1Frequency*t)+noise*np.sin(2*np.pi*(signal1Frequency+1)*t)

# plt.plot(t,5+E1,label='E1')

fourierTransform = np.fft.fft(E1)/len(E1)      # Normalize amplitude
fourierTransform = fourierTransform[range(int(len(E1)/2))] # Exclude sampling frequency


tpCount     = len(E1)
print(tpCount)
values      = np.linspace(0,tpCount/2,tpCount/2)# np.arange(int(tpCount/2))
timePeriod  = tpCount/samplingFrequency
timePeriod  = 10
print(timePeriod)
# frequencies = np.fft.fftfreq(int(tpCount/2),d=samplingFrequency)
frequencies = values/timePeriod

Frequencies = values/timePeriod
# Frequency domain representation
plt.figure()
plt.title('Fourier transform depicting the frequency components')
plt.plot(frequencies, abs(fourierTransform),'-o',label = 'Frequency of each laser')
plt.xlabel('Frequency')
plt.legend()
plt.ylabel('Amplitude')
plt.xlim([0,10])

