# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 14:46:30 2023

@author: ignacio

[1] Nagourney-W. Quantum Electronics for Atomic-Physics (2010)

This script calculates the theoretical resolution value of Optical Reference 
according to the finesse of the cavity using control theory.
However, there is an issue with the calculation of the proportional gain (K)
which causes a mismatch with the experimental values. 
Further investigation and refinement of the algorithm is necessary to improve 
the accuracy of the calculations and correct the theoretical approximations. 
The cavity is locked using a Pound-Drever-Hall loop (PDH Lock)

"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import control
import scipy.special as sp

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


def PDH_Error(R1,R2,P,I, B=0, ErrorSignal_Voltage = 0.1, plot_number=2,
              f_dev = 10, Amp = 210e6,Laser_drift = False,Cavity_drift = False, 
              Bode_Plot = False, plots_response = False):
    
    # definition of fixed quantities
    SPEED_OF_LIGHT  = 299792458.0     # Speed of light [m/s]
    CAVITY_LENGTH   = 0.186           # Cavity Length in Meters
    VACUUM_REFRACTIVE_INDEX = 1.00      
   
    def calc_cavity_params(r1, r2):
        
        ''' 
        This function returns cavity parameters according to the reflectivity
        of the mirrors
        '''        
        
        finesse = 2 * np.pi / -np.log(r1 * r2)
        fsr = SPEED_OF_LIGHT / (2 * VACUUM_REFRACTIVE_INDEX * CAVITY_LENGTH) 
        linewidth = fsr / finesse
        roundtrip_time = 1 / fsr
        t_storage = roundtrip_time / -np.log(r1 * r2)
           
        return fsr, finesse, linewidth, roundtrip_time, t_storage
    
    # And here we run it
    fsr,finesse,linewidth,roundtrip_time, τ_C = calc_cavity_params(R1, R2)

    # Simulation time (1s)
    T = .1
    t = np.arange(0,T,1e-4)
    
    # Definition of the transfer functions that compose the entire system   
    
    kk = np.sqrt(R1*R2)/(1-R1*R2)

    K = kk*ErrorSignal_Voltage/linewidth

    omega_n = (SPEED_OF_LIGHT*np.pi/CAVITY_LENGTH)
    TF_cavity = control.tf([K*omega_n**2], [1, 2 * np.pi * linewidth, omega_n**2])
     
    # Create the PI Control transfer function
    kp, ki = P, I #P/TI
    TF_PI = control.tf([ki, kp], [1, 0])
    
    # Create the laser transfer function
    K_L, τ_L = 150E6,  9.5e-8 
    TF_laser = control.tf([K_L], [0,1])    
    
    # Calculate the open loop transfer function
    G =  TF_PI * TF_laser * TF_cavity 
    
    DeltaCavityLaser = Amp*np.sin(2*np.pi*f_dev*t)+B*t
    F0 = control.tf([B],[1,0,0])
    F1 = control.tf([Amp*2*np.pi*f_dev],[1,0,4*np.pi**2*f_dev**2])
        
    # EQ. 11.118 from Chapter "Laser frequency stabilization and control systems" [1]
    F0_s =  F1/(1+G)    
    # F0_s = (F0*G + F1)/(1+G)
    
    # Simulation of the control loop being active
    t, Yi   = control.impulse_response(F0_s,t)
    
    FactorPDH = 1e3
    FactorCavity = 0.960
    Error = np.mean([np.abs(np.min(Yi)/FactorPDH*FactorCavity),np.max(Yi)/FactorPDH*FactorCavity])
    
    if plots_response: 
        # Plot the response
        TT = np.linspace(t[0],t[-1],100)        
        plt.plot(t, Yi/FactorPDH,label = 'Response - PDH Lock active')
        plt.xlabel('Time [s]')
        plt.ylabel('Frequency [kHz]')        
        plt.legend()
        
    return Error
    

if __name__ == "__main__":
    
    Proportional_gain   = 700
    Integration_gain    = 300  #Proportional_gain*1.2*10
    
    plot_number         = 2
    Bode_Plot           = False
    Laser_drift         = False
    Cavity_drift        = False 
    plots_response      = False
    
    # System configuration, obtained from the experimental values
    ErrorSignalPeakToPeakVoltage = 0.150 #[V]
    f_dev   = 10 #[Hz]
    Amp     = 210e6 #[Hz] 
    
    r2 = np.linspace(.80,.9999,50)
    
    r1 = [.85,.90,.95,.99,.993]    
    
    for R1 in r1:
        error = []
        for R2 in r2:        
        
            Error =   PDH_Error(R1,R2,Proportional_gain,Integration_gain,
                      Bode_Plot=Bode_Plot,
                      Laser_drift = Laser_drift,
                      plot_number=plot_number,
                      Cavity_drift=Cavity_drift,
                      ErrorSignal_Voltage = ErrorSignalPeakToPeakVoltage,
                      f_dev = f_dev, Amp = Amp,
                      plots_response = False)
            
            error = np.append(error,Error)
        plt.plot(r2*100,error,label='$R{_1}$='+str(R1*100)+'%',linestyle='dashed')
        plt.legend()
    plt.ylabel('Error [pm]')
    plt.xlabel('Sample reflectivity [%]')
    plt.ylim([-10,400])
    plt.grid()
    plt.title('Theoretical error vs. $R_2$')
    plt.savefig('Theoretical_error.png')
    
'''
To include the information regarding the slope of the error signal "m" in the 
transfer function G(s), you can introduce an additional factor "m" in the 
numerator of the transfer function. This additional factor can be interpreted 
as the gain of the system when the laser is scanned across the cavity resonance 
without any feedback loop implemented on the laser. Therefore, the modified 
transfer function G(s) can be expressed as:
G(s) = K*m / (s + pi * Δf)
Here, "m" is the slope of the error signal obtained when the laser is scanned 
across the cavity resonance without any feedback loop implemented on the laser, 
and "K" is a constant that depends on the physical parameters of the Fabry 
Perot Cavity.
It is important to note that "m" and "K" are not the same thing. "m" represents 
the gain of the system in the out-of-loop regime, while "K" represents the 
overall transfer function of the cavity and the control system when the laser 
is locked to the cavity resonance using a feedback loop. Therefore, even though 
both "m" and "K" have units of volts per hertz, they are not equal to each 
other.
'''