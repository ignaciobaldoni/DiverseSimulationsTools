# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 12:02:42 2020

@author: ibaldoni
"""

import numpy as np
import matplotlib.pyplot as plt

k = 1.3806504e-16 
T = 300
rho = 2.2
C = 6.7e6
R = 50e-6
Veff = 10e-9

dT2 = (k * T**2)/(rho * C * Veff)  

result = np.sqrt(dT2)
#print('Variance of temperature fluctuations in volume V = %s uK' % (result*1e6))


#print(1.45e-5 * 30e-6 / 1.405)

c = 299792458
f = c/1500e-9

#print(2.45e-5 * 60e-6 / 1.996 * f * 1e-3)


FSR = 19.6
power = 193e12/(FSR*1e9)

Power_to_dB = 10*np.log10(power**2)
#print(Power_to_dB)


Diameter    = 3.770e-3  # Diameter of the ring in m
alpha_SiN   = 3.6e-6    # Thermal expansion coefficient of Silicon Nitride

alpha_Si    = 2.6e-6    # Thermal expansion coefficient of Silicon (Waver) 
alpha_SiO   = 5.5e-7    # Thermal expansion coefficient of Silicon (Fused Silica)

lchip = 5e-3 # 5 mm
   
Temperatures = np.linspace(283,323,50)

T0 = 300


deltaT = [(i-T0) for i in Temperatures]
deltaT = np.asarray(deltaT)


def dn_Si3N4_dT(T):
    return 3.211e-7 - 1.990e-8 * T + 8.856e-10 * T**2 \
                    -1.375e-12 * T **3 - 1.105e-15 * T **4

print("Strong claim: n_g doesn't depend on the temperature")
dD_dT = alpha_Si*deltaT*Diameter

rep_rate = 12.13E9 # Hz
n_g = c/(np.pi*Diameter*rep_rate)

n_SiN = 1.9920688512833209
n = n_SiN

def frep_dL(deltaT,n):
    dL = alpha_Si*deltaT*Diameter
    new_D = Diameter + dL
    return c/(new_D*np.pi*n)

def dfrep_dT_v1(n_g,L,dn_dT,deltaT):
    return -((c/L)*(1/n_g**2)*dn_dT * deltaT)

plt.plot(deltaT,(rep_rate+dfrep_dT_v1(n,np.pi*Diameter,2.45e-5,deltaT))*1E-9,label="Thermorefractive contribution")
plt.plot(deltaT,(frep_dL(deltaT,n_g))*1E-9,label="Si Thermal expansion contribution")
plt.xlabel("Temperature variation $\Delta$T (°C)")
plt.ylabel("Frequency variation (MHz)")
plt.hlines(rep_rate*1E-9,deltaT[0],deltaT[-1],linewidth = 0.7, linestyle = "--")
plt.legend()

def dfrep_dT(Diameter):
    
    L = np.pi*Diameter
    dn_dT = 2.45e-5 
    n_SiN = 1.9920688512833209
    n = n_SiN
    dL_dT = alpha_Si*L
    k = 1   #16027 
    
    dfrep_dT = -(k*c)/((L*n)**2)*(n*dL_dT + L *dn_dT)
    
    return dfrep_dT
    
#plt.figure(3)
plt.plot(deltaT,(rep_rate+dfrep_dT(Diameter)*deltaT)*1E-9,
         label = "Both contributions df/dT")
plt.xlabel("Temperature variation $\Delta$T (°C)")
plt.ylabel("Resonance variation (GHz)")
plt.legend()

#plt.figure()
#plt.plot(Temperatures-273,dfrep_dT_v1(n,np.pi*Diameter,2.45e-5,deltaT)*1E-6,
#         label = "Thermo Refractive contribution")
#
#plt.plot(Temperatures-273,(rep_rate-frep_dL(deltaT,n))*1e-9,
#         label = "Thermal expansion contribution")
#plt.ylabel("Frequency variation (MHz)")
#
#plt.plot(Temperatures-273,(dfrep_dT(Diameter)*deltaT)*1E-6,
#         label = "Both effects combined")
#plt.legend()

