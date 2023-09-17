# -*- coding: utf-8 -*-
"""
Created on Wed Sep  6 11:19:15 2023

@author: ibaldoni
"""

import matplotlib.pyplot as plt
import numpy as np

Finesse = 1e4
cwLinewidth = 15e3

cwLocked = cwLinewidth/(2*np.pi*Finesse)

print(f'{cwLocked:.5f}Hz')


planck  = 6.62607015e-34
c       = 299792458.0
Length  = 0.1866
optical_power_carrier = 1e-3
wavelength = 1542.14-9




  
S_f = (np.sqrt(planck*c**3)/8)*1/(Finesse*Length*np.sqrt(wavelength*optical_power_carrier))
print(S_f,'Hz/sqrt(Hz)')
phase_noise = 20*np.log10(S_f)
print(phase_noise,'dBc/Hz')


S_L = (np.sqrt(planck*c*wavelength)/8)*1/(Finesse*np.sqrt(optical_power_carrier))
print(S_L,'Hz/sqrt(Hz)')

phase_noise = 20*np.log10(S_L)
print(phase_noise,'dBc/Hz')
