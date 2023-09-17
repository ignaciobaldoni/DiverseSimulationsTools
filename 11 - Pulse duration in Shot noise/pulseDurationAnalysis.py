# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 09:58:54 2023

@author: ibaldoni

Pulse duration simple analysis from RP Photonics
"""


import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append(r'\\menloserver\MFS\99-Data_Warehouse\02-User_Folders-Public\i.baldoni\python Util functions')
from Plot_aux_functions import Plot_parameters, add_grids
Plot_parameters(width=7)
from util_Functions import shot_noise_pulsed

t0_fs = 3174
t0 = t0_fs*1e-15

fiberLength = 1

GVD_fiber = -25509e-30
beta = GVD_fiber


D2 = GVD_fiber * fiberLength




Tau = t0*np.sqrt(1+(4*np.log(2)*D2/(t0**2))**2)

Tau = 4000E-15

print('Pulse duration:',Tau*1e15,'fs')


pulse_duration  = Tau

Optical_power = 9*(1e-3)
MW_frequency = 10e9
Responsitivity = 0.55
AM2PM_suppression = 30


SN_am_dBc_Hz, SN_pm_dBc_Hz, SN_AM_to_PM_suppression = shot_noise_pulsed(
                                    Optical_Power = Optical_power, 
                                    Responsitivity = Responsitivity,
                                    Harmonic_number = MW_frequency,
                                    tau=Tau,
                                    AM_to_PM_suppression=AM2PM_suppression)




noise_labels = [
    # ('Shot noise (cw)', shot_noise_cw_dBcHz),
    ('Shot noise (AM)', SN_am_dBc_Hz),
    ('Shot noise (PM)', SN_pm_dBc_Hz),
    ('Shot noise (AM to PM)', SN_AM_to_PM_suppression)]

# frequency_offset = np.arange(100,1e7,10)
# for label, value in noise_labels:
#     plt.semilogx(frequency_offset, np.ones_like(frequency_offset) * value,
#                  label=f'{label} = {value:.2f} dBc/Hz', linestyle='--')
# plt.legend()
    


pulses = [40,50,100,200,250,500,1000,2000,3000,4000,5000,10000]
Powers = [1,10,20,50,100]

noise = []
num = 0

for Optical_power in Powers:
    for pulse_duration in pulses:
    
        _, SN_pm_dBc_Hz, _ = shot_noise_pulsed(
                                            Optical_Power = Optical_power*1e-3, 
                                            Responsitivity = Responsitivity,
                                            Harmonic_number = MW_frequency,
                                            tau=pulse_duration*1e-15,
                                            AM_to_PM_suppression=AM2PM_suppression)
        noise.append(SN_pm_dBc_Hz)
    
        
    plt.semilogx(pulses,noise,'-o',label=f'{Optical_power} mW')
    plt.legend()
    noise = []
    num += 1

plt.xlabel('Pulse duration [fs]')
plt.ylabel('Phase Noise [dBc/Hz]')
add_grids()
plt.ylim([-250,-150])

    


# interleaver = 60 m??????


