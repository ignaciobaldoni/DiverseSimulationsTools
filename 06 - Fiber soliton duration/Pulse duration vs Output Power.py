# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 11:40:14 2022

@author: ibaldoni
"""
import numpy as np
import matplotlib.pyplot as plt

AvgPower = np.linspace(6,11,50)

p0 = np.linspace(500.1000,5000,50)

beta = 22958
gamma = .7904*1e-3
T0 = np.sqrt(beta/(p0*gamma))
repRate = 13e9

Avg = p0*T0*1e-15*repRate/0.88
#plt.figure(num=1)
# plt.plot(p0,Avg, '-o')
# plt.plot(Avg,T0, '-o')

Current_measured = [5000,6000,7000,8000,9000,10000]
Power_measured = np.array([6,7,8,9,10,11])
Time_meas = np.array([200,145,110,90,75,69])

Peak_power_meas = Power_measured*0.88/(Time_meas*1e-15*repRate)

# plt.figure()
#plt.plot(Power_measured,Time_meas,'-o')
#plt.plot(Power_measured,Peak_power_meas,'-o')
plt.plot(Peak_power_meas,Time_meas,'-o')

import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def Pulse_duration(x, a):
    return a * 1/np.sqrt(x)
plt.figure()
xdata = Peak_power_meas
ydata = Time_meas
plt.plot(xdata, ydata, 'bo-', label='data')

popt, pcov = curve_fit(Pulse_duration, xdata, ydata)#,sigma=[1,1,0.01,0.01,.501,.501])


plt.plot(xdata, Pulse_duration(xdata, *popt), 'r-',
         label='fit: a=%5.3f' % tuple(popt))

plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.show()