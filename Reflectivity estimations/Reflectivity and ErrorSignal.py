# -*- coding: utf-8 -*-
"""
Created on Sat Feb 11 19:39:14 2023

@author: ignacio
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sp
import cmath
import itertools
from typing import List

# Plot Parameters
plt.rcParams.update({
    'xtick.labelsize': 15,
    'ytick.labelsize': 15,
    'font.size': 15,
    'figure.figsize': (8,6),
    'savefig.dpi': 300,
    'grid.alpha': 0.75,
    'figure.constrained_layout.use':True
})

plt.style.use('seaborn-dark')

# Define Constants
SPEED_OF_LIGHT = 299792458.0     # Speed of light [m/s]
VACUUM_REFRACTIVE_INDEX = 1.00             # Refractive index of vacuum
AMPLITUDE_OF_ELECTRIC_FIELD = 1
MODULATION_DEPTH = 0.09
    

def def_plots(r1,r2, omega, error_signal, transmission, resonance, finesse, linewidth,normalization=False):
    error_signal_label = f'Error signal\n{SIDEBAND_MODULATION:.2f} MHz Mod'
    finesse_str = f'{finesse:.1f}'
    linewidth_str = f'{linewidth:.2f}'
    resonance_label = f'Resonance (norm.)'
    transmission_label = f'Visibility\nFinesse = {finesse_str}\nLinewidth = {linewidth_str} MHz'
    title_label = f'R1: {r1*100:.1f}% || R2: {r2*100:.1f}%'

    fig, axes = plt.subplots(nrows=2, ncols=1, gridspec_kw={'hspace': 0.10})
    
    axes[0].plot(omega, error_signal, label=error_signal_label)
    if normalization:
        axes[0].text(omega[0],0.5,'Normalized')
    axes[0].set_ylim([-1.1,1.1])
    axes[0].legend()
    axes[0].set_title(title_label)
           
    # axes[1].plot(omega, resonance, label=resonance_label)
    axes[1].plot(omega, transmission, label=transmission_label)  
    axes[1].set_ylim([-0.1,1.1])
    axes[1].legend()

    for ax in axes:
        ax.grid(True)
        ax.set_xlabel('Frequency sweeping [MHz]')
        ax.set_ylabel('Intensity')
    
    plt.savefig(f'Error_signal_visibility_{r1*100}_{r2*100}.png')

def after_plots(num, r1,r2, variable, variable_str):
    label = f' $R_2$ = {r1*100:.2f}%'
    unit = '[pm]' if variable_str == 'Error' else ''
    plt.figure(num)
    plt.plot(r2*100, variable, label=label)
    plt.xlabel('Sample reflectivity [%]')
    plt.ylabel(f'{variable_str} {unit}')
    plt.title(f'{variable_str} vs $R_2$')
    plt.grid(True, alpha=0.75)
    plt.legend(loc='best', bbox_to_anchor=(0.35, 0.42))
    plt.savefig(f'{variable_str}_{r1:.0f}.png', dpi=300)
    
    
# Calculate cavity parameters
def calc_cavity_params(r1, r2, show_cav_params=True): 
    finesse = 2 * np.pi / -np.log(r1 * r2)
    fsr = SPEED_OF_LIGHT / (2 * VACUUM_REFRACTIVE_INDEX * CAVITY_LENGTH) * 1E-6
    linewidth = fsr / finesse
    roundtrip_time = 1 / (fsr * 1e6)
    decay_time = roundtrip_time / -np.log(r1 * r2)
    Minimum_laser_linewidth = LASER_WAVELENGTH**2 / (2*np.pi * decay_time * finesse) 
    
    if show_cav_params:
        print('FSR:\t\t',fsr, 'MHz')
        print('Finesse:\t',finesse)
        print('Linewidth:\t',linewidth,'MHz')
        print('RT time:\t',roundtrip_time*1e9,'ns')
        print('Decay time:\t',decay_time*1e9,'ns')
        print('Min. laser linewidth:\t',Minimum_laser_linewidth*1e12,'pm')
        
    return fsr, finesse, linewidth, roundtrip_time, decay_time, Minimum_laser_linewidth





def FPC(r1,r2,omega,FSR):
    return (-np.sqrt(r1) + np.sqrt(r2) * np.exp(1j * 2 * np.pi * omega / FSR)) \
           / (1 - np.sqrt(r1 * r2) * np.exp(1j * 2 * np.pi * omega / FSR))




def run_simulation(R1: List[float], R2: np.ndarray, plot_visibility: bool = False,
                   plot_finesse: bool = False, plot_error: bool = False, plotPDH_errorSignal: bool = False,
                   normalization: bool = False, plot_phase: bool = False, show_cav_params: bool = False,
                   sign: int = -1):

    visibility_arr = np.zeros_like(R2)
    finesse_arr = np.zeros_like(R2)
    error_arr = np.zeros_like(R2)
    SNR_arr = np.zeros_like(R2)
    i = 0
    for r1, r2 in itertools.product(R1, R2):
        FSR, finesse, linewidth, _, _, _ = calc_cavity_params(r1, r2, show_cav_params=False)
        omega = np.linspace(-FSR * .05 - SIDEBAND_MODULATION, FSR * .05 + SIDEBAND_MODULATION, 2 ** 14)
        
        # omega = np.linspace(5, FSR * .05 + SIDEBAND_MODULATION, 2 ** 14)
        
        Pc = sp.jv(0, MODULATION_DEPTH) ** 2 * np.abs(AMPLITUDE_OF_ELECTRIC_FIELD)
        Ps = sp.jv(1, MODULATION_DEPTH) ** 2 * np.abs(AMPLITUDE_OF_ELECTRIC_FIELD)
        
        
        E_ref = AMPLITUDE_OF_ELECTRIC_FIELD *(FPC(r1,r2,omega,FSR) * sp.jv(0,MODULATION_DEPTH)*np.exp(1j*omega) +\
                    FPC(r1,r2,omega+SIDEBAND_MODULATION,FSR) * sp.jv(1,MODULATION_DEPTH)*np.exp(1j*(omega+SIDEBAND_MODULATION)) + \
                    - FPC(r1,r2,omega-SIDEBAND_MODULATION,FSR) * sp.jv(1,MODULATION_DEPTH)*np.exp(1j*(omega-SIDEBAND_MODULATION)) )
            
        P_ref = np.abs(E_ref)**2
        
        error_signal = sign*(-2 * np.sqrt(Pc*Ps) * np.imag(FPC(r1,r2,omega,FSR)*np.conj(FPC(r1,r2,omega+SIDEBAND_MODULATION,FSR))-\
                        np.conj(FPC(r1,r2,omega,FSR))*FPC(r1,r2,omega-SIDEBAND_MODULATION,FSR))  )  
        transmission = np.abs(FPC(r1,r2,omega,FSR))**2
        Acirc = 1 / np.abs(1 - r1 * r2 * np.exp(-1j * 2 * np.pi * omega / FSR)) ** 2
        resonance = Acirc / np.max(Acirc)
        visibility = np.min(transmission)
        
        # This values are obtained experimentally from the used system. 
        error_kHz = linewidth * 1e6 * 300e3 /34.63e6
        error_pm  = error_kHz* 0.960/1e3 # To pm/kHz
        error_arr[i] = error_pm        
        visibility_arr[i] = 1-visibility
        
        
        depth = 1-visibility
        T2 = 1-np.sqrt(r1*r2)
        G = 4*np.sqrt(r1*r2)*1/((1-np.sqrt(r1*r2))**2)
        I_signal = T2*G/(1-np.sqrt(r1*r2))*10e-3
        SNR_arr[i] = np.log10(I_signal*np.sqrt(depth/(1-depth)))
        
        
        finesse_arr[i] = finesse

        if plotPDH_errorSignal:
            if len(R2)<5:
                def_plots(r1, r2, omega, error_signal, transmission, resonance,finesse,linewidth)

            else:
                print('Too many PDH error signals to plot!')
        i+=1
        
        if r2 == R2[-1]:
            if plot_finesse: after_plots(2,r1,R2, finesse_arr,'Finesse')
            if plot_visibility: after_plots(3,r1,R2, visibility_arr,'Visibility')
            if plot_visibility: after_plots(5,r1,R2, SNR_arr,'SNR')
            if plot_error: after_plots(4,r1,R2, error_arr,'Error')
            i=0
            
            
        
if __name__ == "__main__":        
    
    plot_Finesse            = False
    plot_Error              = False
    plot_Visibility         = True
    plot_PDH_errorSignal    = False

    CAVITY_LENGTH           = 0.1866          # Cavity Length in Meters
    SIDEBAND_MODULATION     = 20              # Sideband modulation [MHz]      
    LASER_WAVELENGTH        = 1542.14e-9
    
    R1 = [0.85,
           0.94,
           0.95,
           0.990,
           0.999]

    
    R2 = np.arange(0.80,1,0.001) 

    
    run_simulation(R1,R2,
                   plot_finesse=plot_Finesse,
                   plot_error=plot_Error,
                   plot_visibility=plot_Visibility,
                   plotPDH_errorSignal = plot_PDH_errorSignal)