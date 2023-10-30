# -*- coding: utf-8 -*-
"""
Created on Thu Oct 12 13:49:36 2023

@author: ibaldoni
"""
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt

# Define filter parameters
order = 2  # Filter order
cutoff_frequency_hz = 20000.0  # Cutoff frequency in Hz
final_frequency = 50e3


# Sample rate in Hz
sample_rate = final_frequency*2
n_points = 10000


# Normalize the cutoff frequency to the Nyquist frequency
cutoff_frequency_normalized = cutoff_frequency_hz / (sample_rate / 2)

# Create a low-pass Butterworth filter
b, a = signal.butter(order, cutoff_frequency_normalized, 'low')

# Compute the frequency response
w, h = signal.freqz(b, a, worN= n_points)

# Calculate the frequency axis in Hz
frequencies = w / (2 * np.pi) * sample_rate

# Plot the frequency response
plt.figure()
plt.plot(frequencies, 20 * np.log10(abs(h)))
plt.title("Photodiode frequency response (sim.)")
plt.xlabel("Frequency [MHz]")
plt.ylabel("Gain [dB]")
plt.ylim([-80,5])


