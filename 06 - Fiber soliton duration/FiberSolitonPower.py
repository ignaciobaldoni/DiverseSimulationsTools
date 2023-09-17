# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 13:25:47 2022
@author: ibaldoni
"""

import numpy as np
import matplotlib.pyplot as plt

plt.rc('xtick', labelsize=17) 
plt.rc('ytick', labelsize=17) 
plt.rcParams.update({'font.size': 17})
plt.rcParams['figure.figsize'] = (8, 6)

# Fiber parameters
n2      = 2.2E-20       # [m^2/W]
A       = 85E-12        # [m^2]        
disp    = 18e-6         #[ps/nm/km]

# System properties
wl      = 1542e-9       # [m]
T0      = 500e-15       # [s]
f_rep   = 12.1269e9      # [Hz]

# Physical constants
c       = 299789452.0   # [m/s]


#%% Runs the code

def beta2(disp,wl):
    return -disp*wl**2/(2*np.pi*c)

# print(beta2(disp,wl))

def gamma(wl):
    return 2*np.pi*n2/(A*wl)

# print(gamma(wl)*1e3)

def E_soliton(gamma,beta2,t0):
    return 2*abs(beta2(disp,wl))/(abs(gamma(wl)*t0)/1.665)

print(E_soliton(gamma,beta2,T0)*1e9,'nJ')

def Average_power(t0,f_rep,gamma,beta2):
    return np.round(E_soliton(gamma,beta2,t0)*f_rep,4)

print(Average_power(T0,f_rep,gamma,beta2),'W')

Peak_power_for_soliton = abs(beta2(disp,wl))/((T0/1.665)**2*gamma(wl))
print('Peak power for soliton', Peak_power_for_soliton,'W')       

SolitonN = (T0/1.665)*np.sqrt(Peak_power_for_soliton*gamma(wl)/abs(beta2(disp,wl)))
print('Soliton number: ',SolitonN)           


#%% Plotting
durations = np.linspace(50e-15,300e-15,100)
plt.plot(durations*1e15,Average_power(durations,f_rep,gamma,beta2),'-o')
plt.plot(79,Average_power(79e-15,f_rep,gamma,beta2),'o',
          label = '79 fs = %s W' % Average_power(79e-15,f_rep,gamma,beta2))
plt.plot(100,Average_power(100e-15,f_rep,gamma,beta2),'o',
          label = '100 fs = %s W' % Average_power(100e-15,f_rep,gamma,beta2))
plt.plot(150,Average_power(150e-15,f_rep,gamma,beta2),'o',
          label = '150 fs = %s W' % Average_power(150e-15,f_rep,gamma,beta2))
plt.plot(200,Average_power(200e-15,f_rep,gamma,beta2),'o',
          label = '200 fs = %s W' % Average_power(200e-15,f_rep,gamma,beta2))
plt.plot(250,Average_power(250e-15,f_rep,gamma,beta2),'o',
          label = '250 fs = %s W' % Average_power(250e-15,f_rep,gamma,beta2))
plt.plot(300,Average_power(300e-15,f_rep,gamma,beta2),'o',
          label = '300 fs = %s W' % Average_power(300e-15,f_rep,gamma,beta2))

plt.legend()
plt.xlabel('Aimed duration [fs]')
plt.ylabel('Needed average power [W]')
plt.grid()
