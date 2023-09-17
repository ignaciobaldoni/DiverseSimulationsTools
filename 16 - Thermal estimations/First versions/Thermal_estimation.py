# -*- coding: utf-8 -*-
"""
Created on Mon Oct 18 10:11:56 2021

@author: ibaldoni
"""


import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

c = 299792458.0
dn = 2.45e-5
alpha = 3.3e-6


density = 3.29e3 
diameter = 3.776e-3 
center_pump = 1542.142*1e-9
lambda0 = center_pump


width = 1.5e-6
height = 0.8e-6
n = 1.996

deltaT = np.linspace(-10,10,100)
x = deltaT

D = np.pi*diameter

# Wavelength variation as a function of temperature
DeltaLambda = lambda0 * (1/n * dn + 1/D * alpha) * x

#plt.plot(deltaT,DeltaLambda*1e9)
#plt.xlabel('Temperature variation (K)')
#plt.ylabel('Wavelength variation (nm)')



Volume = width * height * np.pi * diameter


mass = density * Volume

begin = center_pump - 0.000001e-9#1542.14199e-9
end = center_pump + 0.00001e-9   #1542.24e-9


linewidthHz = 100e6

linewidth = (1/194e12**2)*linewidthHz

PumpPower = 0.2
coupling = .55
Q = .5e6
Qabs = 1.0e6

K = 30
Cp = 800

OptPower = PumpPower * coupling * Q/Qabs
#specificHeat = 800
#Cp = specificHeat * mass
#K = 30 * np.pi * 3.776e-3

matriz = 5000
#forward_pump = np.linspace(begin,end,matriz+1)
#backward_pump = np.linspace(end,begin,matriz)
#pump = np.concatenate((forward_pump, backward_pump), axis=0)

pump_1 = np.linspace(begin,end,matriz+1)
####  --------- Hausgemacht --------- ####

rise_time = 10

from scipy import interpolate
tt = np.linspace(0, rise_time, matriz+1)
y = pump_1
pump = interpolate.interp1d(tt, y)


X = np.zeros(matriz+1)
DeltaLambda = np.zeros(matriz+1)
X[0] = 0
DeltaLambda[0] = 0

dt = rise_time/(matriz+1)

#dfdt = (700e6/1e0)

yo_lo_digo ='Si' 
if yo_lo_digo == 'Si':
    
    const1 = OptPower/Cp #dfdt 
    const2 = K/OptPower
    
#    real = (1 + (1/n * dn + 1/D * alpha)
    for i in range(0,matriz):
        
        X[i+1] =  X[i] + \
                   dt/Cp * (OptPower*((2*(
                            pump(tt[i]) - lambda0 * (1+alpha*X[i])
                            )/linewidth)**2 + 1)**(-1) \
                     - K*X[i]\
                     )
    
    
    #Detuning = c/lambda0**2 * (pump-lambda0)*1e-6
    plt.plot(tt,X)
    plt.xlabel('Time [s]')
    plt.ylabel('$\Delta$T [K]')
    
    plt.figure()
    plt.plot(tt,1e9*lambda0 * (1/n * dn + 1/D * alpha) * X)
    plt.xlabel('Time [s]')
    plt.ylabel('$\Delta \lambda$ [nm]')




yo_lo_digo ='Si' 
if yo_lo_digo == 'Si':
    # function that returns dy/dt
    def model(y,pump):
        dydt = 2/Cp * (OptPower / ((2*(
                            pump - lambda0 * (1 + (1/n * dn + 1/D * alpha)*y))/linewidth)**2 + 1) - K*y)
        
        return dydt
    
    # initial condition
    y0 = 0.01
    
    # time points
    pump = np.linspace(begin,end,matriz+1)
    
    # solve ODE
    y = odeint(model,y0,pump)
    
    # plot results
    plt.plot(pump,y)
    plt.xlabel('time')
    plt.ylabel('y(t)')
    plt.show()