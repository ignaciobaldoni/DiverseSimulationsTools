# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 13:56:58 2022

@author: ibaldoni
"""

# https://pysdr.org/content/sampling.html
# https://colab.research.google.com/github/varun19299/wireless-lab-iitm/blob/notebooks/11-visualising-rf-spectrum.ipynb

import numpy as np
import matplotlib.pyplot as plt


eta = 1
h = 6.62607015e-34      # m2 kg / s
nu = 194e12             # Hz
Pin = 100e-6            # Watts
q = 1.602176634e-19     # Coulombs
P_out = 1e-3            # Watts

i_in = Pin * q/(h*nu)

print('i_in:\t',i_in,'A')


Be = 1e9 # Hz Bandwidth

N = 1
n2 = 0.519
n1 = N-n2
n_sp = n2/(n2-n1)
print('n_sp:\t',n_sp)

Gain = 10*np.log10(P_out/Pin)
# Gain = 10
G = 10**(Gain/10)

i_out = G*Pin
print('i_out:\t',i_out,'A')

print('Gain:\t',G,'dB')
DeltaNu = 0.6e9

i_sp = n_sp*(G-1)*q*DeltaNu

print('i_sp:\t',i_sp,'A')

S = (G*i_in)**2
print('S:\t\t',S,'A²')

N_shot = 2*Be*q*(G*i_in+i_sp)
print('N_shot:\t',N_shot,'A²')

N_s_sp = 4*G*i_in*i_sp*Be/DeltaNu
print('N_s_sp:\t',N_s_sp,'A²')

N_sp_sp = i_sp**2*Be*(2*DeltaNu-Be)/(DeltaNu**2)
print('N_sp_sp:',N_sp_sp,'A²')

SNR_Out = S/(N_shot+N_s_sp+N_sp_sp)

Noise = (N_shot+N_s_sp+N_sp_sp)
SNR = 10*np.log10(SNR_Out)
print('SNR:\t',SNR,'dB')

Noise_sources = [N_s_sp, N_shot,N_sp_sp]
plt.figure()
plt.plot(10*np.log10(Noise_sources),'o')

NF = 3
sigma_ASE = np.sqrt(2*G**2*NF*h*nu*Pin*Be)
print('Sigma:\t',sigma_ASE)

# print(10*np.log10((G*i_in/i_sp)**2))

# freqss = [50e1]
# for i in freqss:
#     center_freq = i 
#     Frequency = center_freq
    
#     # Number of samples taken per second, is simply 1/T = Fs
#     # According DSP theory --> sample at twice the frequency of the signal in order to remove ambiguities
#     # Sample rate must be “at least twice the frequency of the maximum frequency component”
#     Fs = 3e6 # sample rate [Hz]
#     Ts = 1/Fs # sample period
#     # print(Ts)
#     N = 1024 # number of samples to simulate
    
#     t = Ts*np.arange(N)
#     x = np.exp(1j*2*np.pi*Frequency*t) # simulates sinusoid at x Hz
        
#     x = x[0:N] # we will only take the FFT of the first 1024 samples, see text below
#     PSD = (np.abs(np.fft.fft(x))/N)**2
#     PSD_log = 10.0*np.log10(PSD)
#     PSD_shifted = np.fft.fftshift(PSD_log)
    
#     f = np.arange(Fs/-2.0, Fs/2.0, Fs/N) # start, stop, step.  centered around 0 Hz
#     f += center_freq # now add center frequency
    
#     RBW = Fs/len(np.fft.fft(x))
#     print(RBW,'Hz')
    
#     plt.figure()
#     plt.plot(f, PSD_shifted,label = 'RBW = %s Hz' % RBW)
#     plt.grid()
#     plt.legend()
#     plt.show()


print('---------------------------- With noise ------------------------------')

Fs = 20*DeltaNu # sample rate
Ts = 1/Fs # sample period
N = 2*2048 # number of samples to simulate
Frequency = DeltaNu

t = Ts*np.arange(N)
x = 0.00015*np.exp(1j*2*np.pi*Frequency*t) # simulates sinusoid at 50 Hz

n = (np.random.randn(N) + 1j*np.random.randn(N))/np.sqrt(2) # complex noise with unity power
noise_power = 2.

# Noise = np.ones(N)*Noise**noise_power
# Noise = n * np.sqrt(noise_power)
# Noise = Noise* n * np.sqrt(noise_power)


r = x + 0* Noise #

PSD = (np.abs(np.fft.fft(r))/N)**2
PSD_log = 10.0*np.log10(PSD)
PSD_shifted = np.fft.fftshift(PSD_log)

f = np.arange(Fs/-2.0, Fs/2.0, Fs/N) # start, stop, step

# plt.figure()
# plt.plot(f, PSD_shifted)
# plt.xlabel("Frequency [Hz]")
# plt.ylabel("Magnitude [dB]")
# plt.grid(True)
# plt.show()

# Power on RF spectrum
Psignal = 10E-3 #[mW]
quantum_efficiency = 0.8
R = 50
Total_optical_power = 1
P = 10*np.log10(R/2*(Psignal*quantum_efficiency)**2/Total_optical_power)
print(P)