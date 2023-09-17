# -*- coding: utf-8 -*-
"""
Created on Tue Oct 19 14:53:43 2021

@author: ibaldoni
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy import interpolate

import scipy.optimize as scopt
import scipy.integrate as scint

import h5py

import time

study_thermal_triangle = False
simulation = True


def freq_wavelength(question,value,freq = None,wavelengh = None):
    if question == 'f':
#        print(c/lambda0**2 * value*1e-9,'GHz')
        return c/lambda0**2 * value
    elif question == 'w':
#        print((value*c/(freq0**2))*1e9,'nm')
        return (value*c/(freq0**2))

#read out of the file
fileName = 'C:\\Users\\ibaldoni\\Desktop\\D63_1_12GHz_F13_soliton_stepsGenerated light_38.h5' ##\\menloserver\MFS\03-Operations\02-DCP\03-Entwicklungsprojekte\9556-COSMIC\52-Messergebnisse\20210217_20GHz_multiple_drifts_tests\20210217_20GHz_temp_stab_27100deg_1500mA_TEC11000(lower_branch)_800Vpp_26Hz\3-Results\Frequency trace
dict_Result = h5py.File(fileName, 'r')

timeTrace = dict_Result['measurementData'] \
                            ['fiberResonatorTransmission_traceData'] \
                            ['time_axes']

FiberLoopCavity = dict_Result['measurementData'] \
                            ['fiberResonatorTransmission_traceData'] \
                            ['voltage_axes']
                            
ringResonatorTransmissionTrace \
    = dict_Result['measurementData'] \
                ['ringResonatorTransmission_traceData'] \
                ['voltage_axes']

ringResonator = ringResonatorTransmissionTrace[()]
timeTrace = timeTrace[()]
flc = FiberLoopCavity[()]
flc_2 = 1- flc




if simulation == True:    
    matriz = 550000# len(ringResonatorTransmissionTrace)-1 + 900000
    
    ##### ------------ Fundamental settings ------------ #####  
    
    c = 299792458.0
    center_pump = 1542.142*1e-9#1542.142*1e-9
    lambda0 = center_pump
    freq0 = c/lambda0
    
    pump2 = 1333e-9
    lambda2 = pump2
    
    forward_backward = False
    
    ##### ------------ Silicon nitride properties ------------ #####  
    dn = 2.45e-5
#    alpha = 3.3e-6
    exp_coeff = 2.5e-6
    n = 1.996
    density = 3.29e3 
    alpha = (1/n * dn + exp_coeff)
    Thermal_conductivity = 30
    specificHeat = 800
    print(Thermal_conductivity)
    ##### ------------ Chip dimensions and properties ------------ #####  
            
    diameter = 3.776e-3   #3.776e-3*.5   
    width = 1.5e-6          #1.3e-6
    height = 0.8e-6
    D = np.pi*diameter
    Volume = width * height * D
    mass = density * Volume
    
    Cp = specificHeat * mass
    K = Thermal_conductivity * np.pi * diameter#*  1e-6# V/L**2??
    print(K)
    ##### ------------  Tuning settings ------------ #####  
    
    deltaf = 900e6  # 700 MHz
    deltalambda = deltaf*lambda0**2/c
     
    begin = center_pump - deltalambda/6#1542.14199e-9
    end = center_pump + deltalambda/2 #+ deltalambda/2   #1542.24e-9
    
    rise_time = timeTrace[-1]-timeTrace[0]#0.1e-3  # in seconds
    
    ##### ------------  Other settings of the setup ------------ #####  
    
    linewidthHz = 40e6
    linewidth = (c/freq_wavelength('f',center_pump)**2)*linewidthHz
    
    PumpPower = 0.3725
    coupling = .55
    Q = freq_wavelength('f',center_pump)/linewidthHz
    Qabs = Q*2 #8.0e6
    OptPower = PumpPower * coupling * Q/Qabs
    
    PumpPower2 = 0.2e-24
    coupling = .55
    Q2 = freq_wavelength('f',pump2)/linewidthHz
    Qabs2 = Q*2 
    OptPower_sec = PumpPower2 * coupling * Q2/Qabs2
    
    
    ##### ------------  Magic numbers for simulation ------------ #####  
    
    const1 = OptPower/Cp 
    const2 = K/OptPower
    
    magic_number  = 0.061927939922785974  # For stability of finite differences method
    magic_number2 = 9.5859829171354e-08
    #matriz = np.int(rise_time*const1/magic_number) - 1
    
    print('Number of points:',matriz)
    dt = rise_time/(matriz+1)
    
#    if dt*const2<magic_number2:
#        matriz = np.int(rise_time*const2/magic_number2) - 1
#        print('New matrix calculation')
    
    if forward_backward == True:   
        forward_pump = np.linspace(begin,end,matriz+1)
        backward_pump = np.linspace(end,begin,matriz)
        pump_1 = np.concatenate((forward_pump, backward_pump), axis=0)
        kk = 2
    else:
        kk = 1
        pump_1 = np.linspace(begin,end,matriz+1)
    
    
    ##### ------------  Thermal triangle calculation ------------ #####  
        
    X = np.zeros(kk*matriz+1)
    X2 = np.zeros(kk*matriz+1)
    a = np.zeros((kk*matriz+1),dtype = 'complex')
    
    tt = np.linspace(timeTrace[0],timeTrace[-1],kk*matriz+1)
    
    y = pump_1
    pump = interpolate.interp1d(tt, y)
    
    const1 = OptPower/Cp #dfdt 
    const2 = K/OptPower
    
    const3 = OptPower_sec/Cp
    const4 = K/OptPower_sec
    
    X[0] = 0#dt*const1 * ((2*(pump(tt[0]) - lambda0)/linewidth)**2 + 1)**(-1) + \
           #dt*const3 * ((2*(pump2 - lambda2)/linewidth)**2 + 1)**(-1)
    X2[0] = 0
      
#    scale1 = int(len(tt)/4)
#    scale = [tt[0],tt[scale1],tt[scale1*2],tt[scale1*3],tt[-1]]
    
#    dt2 = scale1
    k1 = K
    k2 = 740
    Cp2 = 1.38
    dt2 = dt
    
    gammaT  = freq_wavelength('f',center_pump)/Q
    GammaT  = K
    gT      = -freq_wavelength('f',center_pump)/n*dn 
    gammaE  = freq_wavelength('f',center_pump)/Qabs * 1e-12
    A       = 0.5
    A2 = OptPower
    etaT = 0 
    
    
    
    detuning1 = np.linspace(-1000000,10000000,matriz)
    print(dt)
    
    for i in range(0,kk*matriz):
#        detuning =freq_wavelength('f',pump(tt[i])) - freq_wavelength('f',lambda0)
        detuning=detuning1[i]
#        time.sleep(0.1)

        a[i+1] = a[i]+dt*((1j*detuning-gammaT/2)*a[i]-1j*gT*X[i]*a[i]+1j*np.sqrt(gammaE)*A)
        
        
        
#        etaT    = 1e-12
        X[i+1] = X[i] + dt*(-GammaT * X[i] + etaT*np.real(np.conj(a[i])*a[i]))
        
        
#
#        X[i+1] =  X[i] + \
#                  dt*1/Cp*np.real(np.conj(a[i])*a[i]) * (
#                  ((2*(pump(tt[i]) - lambda0 * (1 + alpha*(X[i]+X2[i])))/linewidth)**2 + 1)**(-1) \
#                            - K/A2*X[i])
#                  
#        A2 = np.real(np.conj(a[i])*a[i])
#         

        

#        var_power = 1#(1+i*1e-5)
#        X[i+1] =  X[i] + \
#                  dt*const1*var_power * (
#                  ((2*(pump(tt[i]) - lambda0 * (1 + alpha*(X[i]+X2[i])))/linewidth)**2 + 1)**(-1) \
#                            - const2*X[i]) + \
#                  dt*const3*(
#                  ((2*(pump2 - lambda2 * (1 + alpha*X[i]))/linewidth)**2 + 1)**(-1) 
#                            - const4*X[i]\
#                            )
#         
##        if any(tt[i]==scale):
##            X2[i+1] = X2[i]+k1/Cp2*dt2*X[i]-k2*dt2/Cp2*X2[i]
##        else:
##            X2[i+1]=X2[i]
#        X2[i+1] = X2[i]+k1/Cp2*dt2*X[i]-k2*dt2/Cp2*X2[i]
                

    
    temp_v = np.round(np.max(X)-X[-1],4)
    print('Temperature variation = %s K'%temp_v)
    
    deltaLambda = lambda0 * (1/n * dn + 1/D * exp_coeff) * X
    freq_v = np.round(np.max(1e-9*freq_wavelength('f',deltaLambda))-1e-9*freq_wavelength('f',deltaLambda)[-1],4)
    print('Frequency variation = %s GHz'%freq_v)
    
    
    
    
    
    #### ------------  Plot results ------------ #####  
    
    x_plot = tt # (c/center_pump-c/pump(tt))*1e-9 
    
    fig, axs = plt.subplots(4,figsize=(12,8))
    fig.suptitle('Preliminary results plots')
    
    axs[0].plot(x_plot*1e3,X,label='Cp = %s' % specificHeat)
    axs[0].set_ylabel('$\Delta$T [K]')
#    axs[0].legend()
    
    axs[1].plot(x_plot*1e3,deltaLambda*1e9,label='const1 = %s \n const2 = %s' % (const1,const2))
    axs[1].set_ylabel('$\Delta \lambda$ [nm]')
#    axs[1].legend()
    
    axs[2].plot(x_plot*1e3,1e-9*freq_wavelength('f',deltaLambda),label='const1 = %s \n const2 = %s' % (const1,const2))
#    axs[2].set_xlabel('Time [ms]')
    axs[2].set_ylabel('$\Delta $f [GHz]')
#    axs[2].legend()
    
    axs[3].plot(x_plot*1e3,c/pump(tt)*1e-6,label='const1 = %s \n const2 = %s' % (const1,const2))
    axs[3].set_xlabel('Time [ms]')
    axs[3].set_ylabel('Laser detuning [MHz]')
    


if study_thermal_triangle == True:
    
    thermal_triangle_begin = np.where(ringResonator>0.09)[0][0]
    print(thermal_triangle_begin)
    thermal_triangle_finish = np.where(ringResonator==np.max(ringResonator))[0][-1]
    print(thermal_triangle_finish)
    
    
    plt.plot(timeTrace,ringResonatorTransmissionTrace)
    plt.plot(timeTrace[thermal_triangle_begin],ringResonator[thermal_triangle_begin],'o')
    plt.plot(timeTrace[thermal_triangle_finish],ringResonator[thermal_triangle_finish],'o')
    plt.plot(timeTrace, 1-flc)
    
    from scipy.signal import find_peaks
    peaks, _ = find_peaks(flc_2,distance = 4000, height=.535)
    
    freq = [0,180, 360,540,720]
    unit = timeTrace[peaks]
    #unit = peaks
    freq_calib = interpolate.interp1d(unit, freq,fill_value='extrapolate')
    
    plt.plot(timeTrace[peaks],flc_2[peaks],'o')
    
    plt.figure()
    plt.plot(freq_calib(timeTrace),ringResonator)
    
    width_thermal_triangle = freq_calib(timeTrace[thermal_triangle_finish])-freq_calib(timeTrace[thermal_triangle_begin])
    
    print('Thermal triangle = ',width_thermal_triangle,'MHz')
    print('Linewidth = ',linewidthHz*1e-6,'MHz')


    tauCharge = Q/(c/center_pump)
