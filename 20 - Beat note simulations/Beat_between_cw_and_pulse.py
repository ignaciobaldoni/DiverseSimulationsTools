# -*- coding: utf-8 -*-
"""
Created on Mon May 22 12:16:04 2023

@author: ibaldoni
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import butter, lfilter#, fftshift, fftfreq

# CW laser parameters
cw1_frequency = 1e15  # Frequency of the first CW laser in Hz
cw2_frequency = cw1_frequency + 100e6  # Frequency difference of 100 MHz

# Fiber coupler parameters
coupling_ratio = 0.5  # Coupling ratio of the fiber coupler

# Low-pass filter parameters
cutoff_frequency = 250e6  # Cutoff frequency of the low-pass filter in Hz

# Time parameters
time_resolution = 1e-12  # Time resolution in seconds
duration = 1e-6  # Signal duration in seconds
time = np.arange(0, duration, time_resolution)

# CW laser signals
cw1_signal = np.sin(2 * np.pi * cw1_frequency * time)
gaussian_pulse = np.exp(-(time - duration / 2) ** 2 / (2 * (1e-12) ** 2))
cw2_signal = np.sin(2 * np.pi * cw2_frequency * time) * gaussian_pulse

# Fiber coupling
sum_signal = coupling_ratio * cw1_signal + coupling_ratio * cw2_signal
difference_signal = coupling_ratio * cw1_signal - coupling_ratio * cw2_signal

# Low-pass filter
sampling_rate = 1 / time_resolution
normalized_cutoff = cutoff_frequency / (0.5 * sampling_rate)
b, a = butter(4, normalized_cutoff, btype='low', analog=False, output='ba')
filtered_signal = lfilter(b, a, difference_signal)

# Frequency analysis
frequency = np.fft.fftshift(np.fft.fftfreq(len(time), d=time_resolution))
difference_spectrum = np.abs(np.fft.fftshift(np.fft.fft(difference_signal))) ** 2
filtered_spectrum = np.abs(np.fft.fftshift(np.fft.fft(filtered_signal))) ** 2

# Plotting
plt.figure(figsize=(12, 6))

# Time-domain plot
plt.subplot(121)
plt.plot(time, difference_signal, label='Difference Signal')
plt.plot(time, filtered_signal, label='Filtered Signal')
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.title('Difference Signal and Filtered Signal')
plt.legend()
plt.grid(True)

# Frequency-domain plot
plt.subplot(122)
plt.plot(frequency, difference_spectrum, label='Difference Spectrum')
plt.plot(frequency, filtered_spectrum, label='Filtered Spectrum')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Power')
plt.title('Difference Spectrum and Filtered Spectrum')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()
