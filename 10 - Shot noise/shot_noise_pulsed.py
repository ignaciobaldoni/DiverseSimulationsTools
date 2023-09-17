# -*- coding: utf-8 -*-
"""
Created on Thu May  4 14:13:55 2023

@author: ibaldoni
"""

# Shot noise
import numpy as np

def shot_noise_cw(nu  = 194.4e12, Optical_Power = 1e-3, Responsitivity = 0.9):
    
    # Natural constants
    electric_charge = 1.602e-19
    planck_constant = 6.62607015e-34
    Load_resistance = 50
    
    photocurrent = Responsitivity * Optical_Power 
    print(f'Photocurrent \t= {photocurrent*1E3} mA')
    
    shot_noise_current = 2*electric_charge*Responsitivity*Optical_Power
   
    shot_noise_power = shot_noise_current*Load_resistance
    # print(f'shot noise power at cw = {shot_noise_power} W')

    shot_noise_power_dBc_Hz = 10*np.log10(shot_noise_power*1e3)
    print(f'Shot noise (cw) = {shot_noise_power_dBc_Hz:.2f} dBc/Hz')
    
    return shot_noise_power_dBc_Hz, shot_noise_power
    

def shot_noise_pulsed(nu  = 194.4e12, Optical_Power = 1e-3, Responsitivity = 0.9, MW_signal=10e9, 
                      tau = 100e-15, Harmonic_number = 10e9):
    
    '''Following the paper from Quinlan et al. (2013)
    Analysis of shot noise in the detection
        of ultrashort optical pulse trains
    '''
    
    carrier_power_in_dBm = 10*np.log10(Optical_Power*1e3)
    # print(f'Power in dBm = {carrier_power_in_dBm}')
    
    tau_g = tau/(2*(np.sqrt(np.log(2))))
    # print(f'tau_g = {tau_g*1e15} fs')
    
    # Natural constants
    electric_charge = 1.602e-19
    planck_constant = 6.62607015e-34
    Load_resistance = 50
    
    photocurrent = Responsitivity * Optical_Power 
    print(f'Photocurrent \t= {photocurrent*1E3} mA')
    
    
    shot_noise_power = (electric_charge*(2*np.pi*Harmonic_number*tau_g)**2)/(2*photocurrent)
   
    # print(f'shot noise power = {shot_noise_power} W')

    shot_noise_power_dBc_Hz = 10*np.log10(shot_noise_power)#*1e3)-carrier_power_in_dBm
    print(f'Shot noise (PN) = {shot_noise_power_dBc_Hz:.2f} dBc/Hz')
    
    L_AM = electric_charge/photocurrent
    LAM = 10*np.log10(L_AM)#*1e3)
    print(f'Shot noise (AM) = {LAM:.2f} dBc/Hz')
    
    return shot_noise_power_dBc_Hz, shot_noise_power


shot_noise_level_dBcHz, _ = shot_noise_pulsed(Optical_Power = 8e-3, 
                                              Responsitivity=0.3,
                                              Harmonic_number= 10e9,
                                              tau=75e-15)






















# def shot_noise_cw(nu  = 194.4e12, P = 1e-3, R = 0.9):
    
#     # Natural constants
#     electric_charge = 1.602e-19
#     planck_constant = 6.62607015e-34
#     Load_resistance = 50
    
#     shot_noise_current = 2*electric_charge*R*P
   
#     shot_noise_power = shot_noise_current*Load_resistance
#     print(f'shot noise power = {shot_noise_power} W')

#     shot_noise_power_dBc_Hz = 10*np.log10(shot_noise_power*1e3)
#     print(f'shot noise power = {shot_noise_power_dBc_Hz} dBc/Hz')
    
#     return shot_noise_power_dBc_Hz, shot_noise_power
    

# def shot_noise_pulsed(nu  = 194.4e12, P = 1e-3, R = 0.9, MW_signal=10e9, 
#                       tau = 100e-15, repRate = 100e6):
    
#     '''Following the paper from Quinlan et al. (2013)
#     Analysis of shot noise in the detection
#         of ultrashort optical pulse trains
#     '''
    
#     tau_g = tau/(2*(np.sqrt(np.log(2))))
#     print(f'tau_g = {tau_g*1e15} fs')
    
#     # Natural constants
#     electric_charge = 1.602e-19
#     planck_constant = 6.62607015e-34
#     Load_resistance = 50
    
#     photocurrent = R*P  #0.01 
    
#     L_AM = electric_charge/photocurrent
#     LAM = 10*np.log10(L_AM)
#     print(f'LAM = {LAM} dBc/Hz')
    
    
#     print(f'Photocurrent = {photocurrent*1E3} mA')
#     Harmonic_number = 10E9
    
#     shot_noise_current = 2*(electric_charge*(np.pi*Harmonic_number*tau_g)**2)/(photocurrent)
   
#     shot_noise_power_in_dBm = shot_noise_current#*Load_resistance
#     # print(f'shot noise power = {shot_noise_power} W')

#     shot_noise_power_dBc_Hz = 10*np.log10(shot_noise_power_in_dBm)
#     print(f'shot noise power = {shot_noise_power_dBc_Hz} dBc/Hz')
    
    
    
#     return shot_noise_power_dBc_Hz, LAM
    

# shot_noise_pulsed(P = 16e-3, R = 0.94, tau = 12000e-15, MW_signal=12e9)

# # # shot_noise_cw(P = 8e-3, R = 0.3)
# # pulse_durations = np.arange(10,1500,100)
# # lista = []
# # for i in pulse_durations:
# #     shot_noise_power_dBc_Hz, Lam = shot_noise_pulsed(P = 8e-3, R = 0.3, tau = i*1e-15)
# #     lista.append(shot_noise_power_dBc_Hz)

# # Tr = 1/10e9
# # Tr1 = 1/1000e9
# # import matplotlib.pyplot as plt
# # plt.semilogx(pulse_durations*1e-15/Tr, -Lam+lista,'.-')
# # plt.semilogx(pulse_durations*1e-15/Tr1, -Lam+lista,'.-')
    
    