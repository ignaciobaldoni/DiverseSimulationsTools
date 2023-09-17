# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 12:00:00 2020

@author: ibaldoni

TRN calculation 
"""

import numpy as np
import matplotlib.pyplot as plt





def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = printEnd)
    # Print New Line on Complete
    if iteration == total: 
        print()    
    






w = np.linspace(1,1e4,10000) # Frequency of interest range

rho     = 3.29 * 10**3      # Density (kg m3)
n_o     = 1.996             # Refractive index
dn_dT   = 2.45 * 10**(-5)   # Thermo optic coefficient (K^-1)
kappa   = 30                # Thermal conductivity (W^m-1 K^-1)  
C_p     = 800               # Specific heat capacity (J kg^-1 K^-1)
kb      = 1.38 * 10**(-23)  # Boltzmann constant (J/K) 
c       = 299792458

Radius = [0.00025, 0.0005, 0.0025]#, 0.004]
for R in Radius:
#    R = 0.00025                  # Ring radius (m)
    T = 23 + 273           # Bath temperature (K) 
    l = 66                      # Orbital number (fundamental mode)
    wl = 1550 * 1e-9            # Wavelength lambda (m)
    V_eff = 2.8 * (wl/n_o)**3     # Effective mode volume 
    
    dr = 0.1
    dz = 0.1
    td = ((np.pi/4)**(1/3))*(rho*C_p*dr**2)/kappa
    
    m = l - 10
    p = l - m
    
    
    part2 = 1 + ((R**2*rho*C_p*w)/(kappa * 3**(5/3)))**(3/2)+1/6*((R**2*rho*C_p*w)/(kappa*8*l**(1/3)))**2
    
    # Homogeneous microresonator in an inﬁnite heat bath
    S_dt_1 = (kb * T**2 * R**2)/(12*kappa*V_eff)*part2**(-1)
    
    
    # Thermal decomposition method which does not take into account the interaction with the environment
    
#    part3 = np.sqrt(1/(2*p+1))*1/(R*np.sqrt(dz**2-dr**2))*1/(1+(w*td)**(3/4))**2
    
#    S_dt_2 = (kb*T**2)/np.sqrt(np.pi**3*kappa*rho*C_p*w)*part3
    
    f0 = c/wl 
    
    S_df = (f0 * 1/n_o * dn_dT)**2 * S_dt_1
    
    plt.figure(1)
    plt.plot(np.log10(w), np.log10(np.sqrt(S_df)), label = 'R = %s um' % str(R*1e6))
    plt.legend(loc = 'upper right')
    plt.hlines(0,0,np.log10(w[-1])+0.25,linewidth=0.3, linestyles='dashed')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Noise spectral density ( Hz / Hz ^ (1/2) )')
    
    # Conversion to phase noise S_theta
    
    S_theta = np.log10(S_df*np.pi**2/(w**2))
    
    plt.figure(2)
    plt.plot(np.log10(w), S_theta, label = 'R = %s um' % str(R*1e6))
    plt.legend(loc = 'upper right')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Phase noise ( dBc / Hz )')    
    
    
    
#    printProgressBar(0, n, prefix = 'Progress:', suffix = 'Complete', length = 50)
#    for it in range(1,len(detuning)):
#        printProgressBar(it + 1, n, prefix = 'Progress:', suffix = 'Complete', length = 50)
#        dOm_curr = detuning[it] # detuning value
#        sol[it] = r.integrate(r.t+t_st)