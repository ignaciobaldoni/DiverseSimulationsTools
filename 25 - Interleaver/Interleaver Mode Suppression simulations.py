# -*- coding: utf-8 -*-
"""
Created on Wed Oct 11 16:08:00 2023

@author: ibaldoni

Based on eq. [3], [4] and [5] from 
2022 - Franklyn Quinlan - The photodetection of ultrashort optical pulse trains
for low noise microwave signal generation
"""

import numpy as np
import sys
sys.path.append(r'\\menloserver\MFS\99-Data_Warehouse\02-User_Folders-Public\i.baldoni\python Util functions')
from Plot_aux_functions import Plot_parameters 

import matplotlib.pyplot as plt
import math
from scipy import signal
from scipy.signal import find_peaks
from FSWP_csvDataReader import FSWP_Spectrum_csvReader
import scipy.constants as c
from math import copysign
import pandas as pd

### Constants in the script 
Plot_parameters(width=8)
todBm = 1e3
noise_floor = -78

frequency = np.linspace(1e6, 50e9, 200000)
step = 25e6#*1.25
frequency = np.arange(0,50e9+step,step)



###############################################################################
###############################################################################
def Gf(tau_d,f,loss_factor = 1/2): 
    return loss_factor*(1+np.exp(1j*tau_d*2*np.pi*f))

def MW_Power_spectrum(G_total):
    return np.abs(G_total)**2

def tau_interleaver(RepRate,stage):
    return (1/RepRate)/(2**stage)

def Gf_realCase(tau_d,f,tau1=0,tau2=0,tau_in=0,
                Alpha=1/2,Beta=1/2, 
                coupler_loss_factor=0.9):
    
    # 0.9 because the insertion loss of the coupler usually is 0.5 dB
    
    GF = Alpha*np.exp(np.pi**2*f**2*(tau1**2-tau_in**2))+\
        Beta*np.exp(np.pi**2*f**2*(tau2**2-tau_in**2))*\
        np.exp(1j*tau_d*2*np.pi*f)
    
    return coupler_loss_factor*GF

###############################################################################

def find_local_extrema(amplitudes, frequencies, repRate = 125, stage_number=4,
                       plots = False):

    # Find local maxima and minima using scipy's find_peaks
    maxima, _ = find_peaks(amplitudes,distance=50)#,distance=repRate*2**stage_number*2)#,height=-11)
    minima, _ = find_peaks(-amplitudes,distance=4)

    # Plot the function and mark the local maxima and minima
    if plots:
        plt.plot(frequencies[maxima], amplitudes[maxima], 'ro', label='Local Maxima')
        plt.plot(frequencies[minima], amplitudes[minima], 'go', label='Local Minima')
        plt.legend(fontsize = 12,#loc = 'upper center'
                   )

    return minima, maxima, amplitudes[minima],amplitudes[maxima]

###############################################################################

def P_u(frequency, R = 50, responsitivity_pd = .55, 
        optical_power = 6.5e-3,
        Cutoff_frequency_hz=10000,
        plots = False):
    
    
    def PD_TF():        
        # Define filter parameters
        order = 2  # Filter order
        final_frequency = 70000 #cutoff_frequency_hz*2
    
        # Sample rate in Hz
        sample_rate = final_frequency*2
        n_points = len(frequency)
    
        # Normalize the cutoff frequency to the Nyquist frequency
        cutoff_frequency_normalized = Cutoff_frequency_hz / (sample_rate / 2)
        # print('Cutoff f:', cutoff_frequency_normalized*final_frequency,'GHz')
    
        # Create a low-pass Butterworth filter
        b, a = signal.butter(order, cutoff_frequency_normalized, 'low')
    
        # Compute the frequency response
        w, H_f = signal.freqz(b, a, worN= n_points)
    
        # Calculate the frequency axis in Hz
        frequencies = w / (2 * np.pi) * sample_rate
    
        if plots:
            # Plot the frequency response
            plt.plot(frequencies/1e3, 20 * np.log10(abs(H_f)))
            plt.title("Photodiode frequency response (sim.)")
            plt.xlabel("Frequency [GHz]")
            plt.ylabel("Gain [dB]")
            plt.ylim([-80,5])
        
        # m = 0.01 # Modulation index
        # cw_modulated =  (0.5*m**2)*(Iavg**2)*R*(np.abs(H_f)**2)    
        
        return H_f
    
    H_f = PD_TF()
    
    Iavg = optical_power*responsitivity_pd
    P_mu_pulse = 2*(Iavg**2)*(np.abs(H_f)**2)*R   

    return P_mu_pulse # This is in Watts

###############################################################################
        
def dispersion(stage, tau0_fs):
    
    def TAU(t0, fiberLength):
        GVD_fiber = -25509e-30
        D2 = GVD_fiber * fiberLength
        return t0 * (1 + (4 * np.log(2) * D2 / t0**2)**2)**0.5
    
    fiber_lengths = {
        1: (1.25, 0.43),
        2: (0.8, 0.392),
        3: (0.7, 0.496),
        4: (0.7, 0.598),
        5: (0.7, 0.649),
        6: (0.7, 0.675)
    }
    
    if stage in fiber_lengths:
        fiberLength1, fiberLength2 = fiber_lengths[stage]
        Tau1 = TAU(tau0_fs, fiberLength1)
        Tau2 = TAU(tau0_fs, fiberLength2)
        
        # print(f'{(fiberLength1 - fiberLength2):.4f} cm')
        # print(f'{((1-Tau2/Tau1)*100):.5f} % between tau1 and tau2')
        return Tau1, Tau2
    else:
        print("Invalid stage")
        return None, None
    
      
        
###############################################################################
###############################################################################

if __name__ == "__main__":    
    
    number_of_stages    = 4
    f0_reprate          = 125e6
    coupler_factor      = 0.5
    Cutoff_Freq_MHz_PD  = 20000
    Tau_In              = 1.522e-12
    plots               = True
    saveplots           = False
    
    
    
    first_harmonic = (f0_reprate*1e-6)*(2**number_of_stages)
    print(f'Number of stages: {number_of_stages} \n'+\
          f'Expected first frequency = {first_harmonic:.2f} MHz')
    print('----------------------------------------------------------------')
    
    
    
    
    
    gf_product = 1.0
    for i in range(1, number_of_stages+1):
        
        
        Tau_d = tau_interleaver(f0_reprate, i)
#        error = 0.0*Tau_d
#        Tau_d = Tau_d-error
        ds = Tau_d/1.4682*c.c
 ##############################################       
        d_ds = 1*150e-6#500e-6
 ##############################################       
        t_d_ds = d_ds*1.4682/c.c
        Tau_d = Tau_d+t_d_ds
        
        # Let's consider dispersion in every stage.
        Tau1, Tau2 = dispersion(i,Tau_In)              

        # ### Pulse duration does not seem to affect greatly the suppression of 
        # ### the unwanted modes
        # Tau1 = Tau_In+.27e-12*(number_of_stages-i)
        # Tau2 = Tau1*0.9
        

# The power distribution in the interleaver stage does have a clear influence. 
# A random value between 0.45 and 0.55 is stablished to account for unpredicted
# losses because of the coupler manufacturing issues or due to splicing.  

        ALPHA = np.random.uniform(0.45, 0.55)     
        
        if i == 1: ALPHA = 0.5
        if i == 2: ALPHA = 0.5
        if i == 3: ALPHA = 0.5
        if i == 4: ALPHA = 0.5
        if i == 5: ALPHA = 0.5
        if i == 6: ALPHA = 0.5

        
        BETA = (1-ALPHA)
        
        # Here we implement the realistic case
        gf = Gf_realCase(Tau_d, 
                         frequency,
                         tau1=Tau1,
                         tau2=Tau2,
                         tau_in=Tau_In,
                         Alpha = ALPHA,
                         Beta = BETA)  
        
        
    # In case we want to see the "ideal" case, activate the following function.
        #Gf(Tau_d, frequency)
        
        gf_product *= gf
        
        
        print(f'----------- Stage {i} -----------')
        print(f'Coupling: arm1 {ALPHA:.2f} - arm2 {BETA:.2f}')
        print(F'Initial tau = {Tau_In*1e12:.4f} ps')
        
        Tau_In = Tau1
        
        
    
    
    print('----------------------------------------------------------------')   
    GTotal = gf_product*coupler_factor
    MW_spectrum = MW_Power_spectrum(GTotal)*\
        (P_u(frequency,Cutoff_frequency_hz=Cutoff_Freq_MHz_PD))
    GTotal_dBm = 10*np.log10(MW_spectrum*todBm)    
    #GTotal_dBm_fr = GTotal_dBm*(F_r(frequency))
    lf = len(frequency)
    GTDBMMWFR=[]
    for i in np.arange(0, lf,1):
        if frequency[i] % 125e6 == 0:
            GTDBMMWFR.append(GTotal_dBm[i]+6.0)
        else:
            GTDBMMWFR.append(-80.)
    
    
    if plots:
        plt.figure()
        plt.plot(frequency*1e-9, GTDBMMWFR,label='Simulated Microwave Signal',color = 'b')
        plt.ylabel('Power [dBm]')
        plt.xlabel('Frequency [GHz]')
        plt.xlim([0,50])
        plt.ylim([-60,5])
        
    
    
    ###########################################################################
    ######################## We call measured data ########################
    ###########################################################################
    
    folder = r'\\menloserver\MFS\03-Operations\02-DCP\01-DCP_Management\01-Gruppenmanagement\UMS\06_Measurements TBD\20230929 - Interleaver five stages\1 - Raw Data\\'
    
    files = [r"4 stages_50GHz.CSV"]
    
    for i in files:
        
        frequency, phaseNoise = FSWP_Spectrum_csvReader(folder+i,number_of_traces=1)
            
        plt.plot(frequency*1e-9, phaseNoise+6.0,alpha=0.5, color = 'r', label='Measured Microwave Signal')
        #plt.plot(frequency*1e-9, phaseNoise+6.0, label='Measured Microwave Signal')       
        plt.xlabel('Frequency [GHz]')
        plt.ylabel('Amplitude [dBm]')
        #plt.title('11b)                          five stages adjusted',loc='left',fontsize=15)
        plt.title('12b)',fontsize=20,loc = 'left') 
        #plt.title('$\Delta$L=500 $\mu$m and 5% coupler imbalance',fontsize=15)
        plt.title('five stages adjusted',fontsize=15)
        
        
    plt.legend(loc='lower left',framealpha=1.0)
    plt.grid(visible=True, which='both', color='lightgrey', linestyle='-')
    plt.ylim(-60,5)
    plt.xlim(0,50)
    plt.show()
    if saveplots: plt.savefig('SimInterleaver_5050_wholeSpan.png')

    

# # The difference between the Frequency Response of the PD and the plotted value
# # is due to the losses considered in the couplers, splices and coupling ratios 
