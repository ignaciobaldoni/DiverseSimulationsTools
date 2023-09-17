# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 17:48:36 2021

edited on 05.05.2023

@author: ibaldoni
"""



import numpy as np

# Define input parameters
wavelength = 1050.25e-9  # meters
dc_bias_current = 4e-3  # amps OR VOLTS?? Originally was 4
characteristic_impedance = 1100  # ohms

# Define physical constants
speed_of_light = 299792458  # meters per second
planck_constant = 6.62607015e-34  # joule seconds
electron_charge = 1.602176634e-19  # coulombs

# Calculate derived parameters
angular_frequency = speed_of_light / wavelength
small_signal_ext_quant_efficiency = electron_charge / (planck_constant * angular_frequency)
carrier_lifetime = 0.8 / small_signal_ext_quant_efficiency
optical_input_power = dc_bias_current / (characteristic_impedance * 0.8)

# Calculate RIN
current_noise_spectral_density = small_signal_ext_quant_efficiency ** 2 *\
                                carrier_lifetime * optical_input_power * \
                                (1 - carrier_lifetime) * 2 * planck_constant *\
                                angular_frequency + small_signal_ext_quant_efficiency ** 2 * carrier_lifetime ** 2 * optical_input_power * 2 * planck_constant * angular_frequency
voltage_noise_spectral_density = characteristic_impedance ** 2 * current_noise_spectral_density
rin = 10 * np.log10(voltage_noise_spectral_density / dc_bias_current ** 2)

# Print result
print(f"RIN: {rin:.2f} dB")



###############################################################################

c=299789452.0
wo = c/1050.25e-9
hbar = 1.054571817e-34

electron = 1.602176634e-19

R0 = electron/(hbar*wo)

tau = 0.8/R0
Z = 1100
DC = 4
Pin = DC / (Z*0.8)


S_i = R0**2 * tau*Pin * (1-tau) * 2*hbar*wo+\
    R0**2*tau**2*Pin * 2*hbar*wo
    
Sv = Z**2*S_i

Rin = 10*np.log10(Sv/(DC**2))

print(Rin)
