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

#read out of the file
fileName = 'C:\\Users\\ibaldoni\\Desktop\\D63_1_12GHz_F13_soliton_stepsGenerated light_38.h5'
dict_Result = h5py.File(fileName, 'r')

timeTrace = dict_Result['measurementData'] \
                            ['fiberResonatorTransmission_traceData'] \
                            ['time_axes']

# get fiber loop cavity transmission trace
voltageTrace = dict_Result['measurementData'] \
                            ['fiberResonatorTransmission_traceData'] \
                            ['voltage_axes']
                            

ringResonatorTransmissionTrace \
    = dict_Result['measurementData'] \
                ['ringResonatorTransmission_traceData'] \
                ['voltage_axes']

matriz = len(ringResonatorTransmissionTrace)-1


c = 299792458.0
center_pump = 1542.142*1e-9
lambda0 = center_pump
freq0 = c/lambda0

forward_backward = False

##### ------------ Silicon nitride properties ------------ #####  
dn = 2.45e-5
alpha = 3.3e-6
n = 1.996
density = 3.29e3 
Thermal_conductivity = 30
specificHeat = 800

##### ------------ Chip dimensions and properties ------------ #####  
        
diameter = 3.776e-3   
width = 1.5e-6
height = 0.8e-6
D = np.pi*diameter
Volume = width * height * D
mass = density * Volume


def freq_wavelength(question,value,freq = None,wavelengh = None):
    if question == 'f':
        print(c/lambda0**2 * value*1e-9,'GHz')
        return c/lambda0**2 * value
    elif question == 'w':
        print((value*c/(freq0**2))*1e9,'nm')
        return (value*c/(freq0**2))



##### ------------  Tuning settings ------------ #####  

deltaf = 700e6  # 700 MHz
deltalambda = deltaf*lambda0**2/c
 
begin = center_pump - deltalambda/2#1542.14199e-9
end = center_pump + deltalambda/2   #1542.24e-9

rise_time = timeTrace[-1]#0.1e-3  # in seconds

##### ------------  Other settings of the setup ------------ #####  

linewidthHz = 40e6
linewidth = (c/freq_wavelength('f',center_pump)**2)*linewidthHz

PumpPower = 0.3
coupling = .55
Q = freq_wavelength('f',center_pump)/linewidthHz
Qabs = Q*2 #8.0e6

OptPower = PumpPower * coupling * Q/Qabs
Cp = specificHeat * mass
K = Thermal_conductivity * np.pi * diameter # V/L**2??

const1 = OptPower/Cp 
const2 = K/OptPower

magic_number  = 0.061927939922785974  # For stability of finite differences method
magic_number2 = 9.5859829171354e-08
#matriz = np.int(rise_time*const1/magic_number) - 1

print('Number of points:',matriz)
dt = rise_time/(matriz+1)

if dt*const2<magic_number2:
    matriz = np.int(rise_time*const2/magic_number2) - 1
    print('New matrix calculation')


sweep_speed = deltaf/rise_time


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
DeltaLambda = np.zeros(kk*matriz+1)
X[0] = 0
DeltaLambda[0] = 0
    
tt = np.linspace(timeTrace[0],timeTrace[-1],kk*matriz+1)
#tt = np.linspace(0, kk*rise_time, kk*matriz+1)



y = pump_1
pump = interpolate.interp1d(tt, y)

const1 = OptPower/Cp #dfdt 
const2 = K/OptPower
    
for i in range(0,kk*matriz):
    
     X[i+1] =  X[i] + \
              dt*const1 * (
              ((2*(pump(tt[i]) - lambda0 * (1 + alpha*X[i]))/linewidth)**2 + 1)**(-1) \
                        - const2*X[i]\
                        )

print(dt*const1)
print(dt*const2)

##### ------------  Plot results ------------ #####  
    
#Detuning = c/lambda0**2 * (pump-lambda0)*1e-6
x_plot = tt

plt.figure(1)    
plt.plot(x_plot*1e3,X,label='Cp = %s' % specificHeat)
plt.xlabel('Time [ms]')
plt.ylabel('$\Delta$T [K]')
plt.legend()

plt.figure(2)
plt.plot(x_plot*1e3,1e9*lambda0 * (1/n * dn + 1/D * alpha) * X,label='const1 = %s \n const2 = %s' % (const1,const2))
plt.xlabel('Time [ms]')
plt.ylabel('$\Delta \lambda$ [nm]')
plt.legend()







# Fit differential equation
y_data = ringResonatorTransmissionTrace  # Example data




def f(tt, const10,const20): 
    

    
    X = np.zeros(kk*matriz+1)    
    X[0] = 0
    
    for i in range(0,kk*matriz):
    
         X[i+1] =  X[i] + \
                  dt*const10 * (
                  ((2*(pump(tt[i]) - lambda0 * (1 + alpha*X[i]))/linewidth)**2 + 1)**(-1) \
                            - const20*X[i]\
                            )
                       
    return X


guess = [const1,const2]
popt, pcov = scopt.curve_fit(f, tt, y_data, p0=guess, maxfev=10000)


fit = f(tt, *popt)

plt.figure()
plt.plot(tt*1e3, y_data, 'ro', label='data')
plt.plot(tt*1e3, fit, 'b-', label='const1 =%s\n const2 = %s'%(popt[0],popt[1]))
plt.legend(loc='best')
