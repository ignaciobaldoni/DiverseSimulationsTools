# -*- coding: utf-8 -*-
"""
Created on Thu Nov 11 09:55:56 2021

@author: ibaldoni
"""

import numpy as np
import matplotlib.pyplot as plt
FigSize = (12,8)

c = 299792458.0
center_pump = 1050.25*1e-9

freq0 = c/center_pump
print('Center frequency:',freq0*1e-12,'THz')

repRate = 50e9   #50.3e9

q = freq0/repRate
print('Modes:',int(q))

dL  = 2.6e-6
dn  = 2.45e-5
n   = 1.9921
#L   = 3.776e-3*np.pi

L = c/(repRate*n)
print('L = ',L*1e3,'mm')

dfdT = -c*q*(1/(n*L)**2) * (n*L*dL+L*dn)

print('dfdT =',dfdT*1e-9,'GHz/K')


Current = [0.06,0.13,0.19]
Resonance = [freq0+repRate*3,
             freq0+repRate*2,
             freq0+repRate*1]

Resonance = [i*1e-12 for i in Resonance]


fig, ax = plt.subplots(figsize=FigSize)
ax.plot(Current, Resonance,'o-',label = '1050 nm resonances')
ax.set_xlabel('Current (A)',fontsize = 17)
ax.set_ylabel('Resonance position (THz)',color="b",fontsize = 17)
ax.tick_params(labelsize=17)
ax.grid()


center_pump = 1541.94*1e-9
freq2 = c/center_pump
Current = [0.01,0.13,0.2]
Resonance = [freq2+repRate*3,
             freq2+repRate*2,
             freq2+repRate*1]

ax2=ax.twinx()
# plt.figure(figsize=FigSize)
Resonance = [i*1e-12 for i in Resonance]
ax2.plot(Current, Resonance,'ro-',label = '1542 nm resonances')
ax2.set_xlabel('Current (A)',fontsize = 17)
ax2.set_ylabel('Resonance position (THz)',color="red",fontsize = 17)
ax2.tick_params(labelsize=17)
# ax.grid()
plt.legend()