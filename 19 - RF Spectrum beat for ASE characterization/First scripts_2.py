# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 10:03:19 2022

@author: ibaldoni
"""
# Noise bandwidths

import numpy as np
import matplotlib.pyplot as plt


eta = 1
h = 6.62607015e-34      # m2 kg / s
nu = 194e12             # Hz
Pin = 100e-6            # Watts
q = 1.602176634e-19     # Coulombs
P_out = 10e-3           # Watts

i_in = Pin * q/(h*nu)
print('i_in:\t',i_in,'A')

Be = 1e9 # Hz Bandwidth

n_sp = 0.9
N = 1
n2 = 0.59
n1 = N-n2
n_sp = n2/(n2-n1)
print('n_sp:\t',n_sp)

Gain = 10*np.log10(P_out/Pin)
# Gain = 10
G = 10**(Gain/10)

print('Gain:\t',G,'dB')
DeltaNu = 12e12

i_sp = n_sp*(G-1)*q*DeltaNu

print('i_sp:\t',i_sp,'A')

S = (G*i_in)**2
print('S:\t\t',S)

N_shot = 2*Be*q*(G*i_in+i_sp)
print('N_shot:\t',N_shot)

N_s_sp = 4*G*i_in*i_sp*Be/DeltaNu
print('N_s_sp:\t',N_s_sp)

N_sp_sp = i_sp**2*Be*(2*DeltaNu-Be)/(DeltaNu**2)
print('N_s_sp:\t',N_s_sp)

SNR_Out = S/(N_shot+N_s_sp+N_sp_sp)
SNR = 10*np.log10(SNR_Out)
print('SNR:\t',SNR,'dB')

# print(10*np.log10((G*i_in/i_sp)**2))

print('------------------------------------------------')
E0 = 1
signal1Frequency = 7
signal2Frequency = 4

Beat_frequency = np.abs(signal1Frequency-signal2Frequency)


# Fs/N
# Fs sample rate
# N is size of FFT

samplingFrequency   = 100;
# samplingRate = 44100
# At what intervals time points are sampled
samplingInterval       = 1 / samplingFrequency; 
print(samplingInterval)
# Begin time period of the signals
beginTime           = 0; 
endTime             = 10; 

t = np.arange(beginTime, endTime, samplingInterval);
print(len(t))

# t = np.linspace(0,10,1000)

# E1 = E0*np.exp(1j*omega1*t)
# E2 = E0*np.exp(1j*omega2*t)
noise = 0.
E1 = np.sin(2*np.pi*signal1Frequency*t)+noise*np.sin(2*np.pi*(signal1Frequency+1)*t)
E2 = np.sin(2*np.pi*signal2Frequency*t)+noise*np.sin(2*np.pi*(signal2Frequency+1)*t)
Beat = np.sin(2*np.pi*Beat_frequency*t) 

Beat_note = E1+E2



plt.plot(t,5+E1,label='E1')
plt.plot(t,E2+2.5,label='E2')
plt.plot(t,Beat-5,label='Beat')
plt.plot(t,Beat_note-1.5,label='Beat note')
plt.legend()

# Fourier transform

Beat_note = Beat_note/np.max(Beat_note)
Frequ_spec = np.fft.fft(Beat_note)

# Frequency domain representation
fourierTransform = np.fft.fft(Beat_note)/len(Beat_note)      # Normalize amplitude
fourierTransform = fourierTransform[range(int(len(Beat_note)/2))] # Exclude sampling frequency

print(len(fourierTransform))

beat_fourierTransform = np.fft.fft(Beat)/len(Beat)      # Normalize amplitude
Beat_fourierTransform = beat_fourierTransform[range(int(len(Beat)/2))] # Exclude sampling frequency


tpCount     = len(Beat_note)
print(tpCount)
values      = np.linspace(0,tpCount/2,tpCount/2)# np.arange(int(tpCount/2))
timePeriod  = tpCount/samplingFrequency
print(timePeriod)
# frequencies = np.fft.fftfreq(tpCount,d=samplingFrequency)
frequencies = values/timePeriod

# Frequency domain representation
plt.figure()
plt.title('Fourier transform depicting the frequency components')
plt.plot(frequencies, abs(fourierTransform),label = 'Frequency of each laser')
plt.plot(frequencies, abs(Beat_fourierTransform),label = 'Beat between lasers')
plt.xlabel('Frequency')
plt.legend()
plt.ylabel('Amplitude')
plt.xlim([0,10])

