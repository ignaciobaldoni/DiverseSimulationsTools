# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 09:22:41 2022

@author: ibaldoni
"""


import matplotlib.pyplot as plt


plt.rc('xtick', labelsize=17) 
plt.rc('ytick', labelsize=17) 
plt.rcParams.update({'font.size': 17})
plt.rcParams['figure.figsize'] = (10, 6)


import numpy as np

Power = [1.88,3,4.84,5.89,7.15,8.27,9.95,11.7,12.3,13,14.6]

duration = [263,217,181,177,190,203,167,106,90.4,71.6,69]

plt.plot(Power,duration, '-o',label='Input pulse = 20 mW, \n$T$ = 200 fs')
plt.xlabel('Power [W]')
plt.ylabel('Pulse duration [fs]')
#plt.title()
plt.legend()
plt.grid()


# Power = 13.4 W
Input_duration = [100,200,300,400,500,600,700,800,900,1000]

duration = [70,87.4,286,174,74.2,84.6,155,285,472,664]

plt.figure()
plt.plot(Input_duration,duration, '-o',label='no chirp')
plt.xlabel('Input pulse duration [fs]')
plt.ylabel('Pulse duration [fs]')
#plt.title()
plt.legend()
plt.grid()


# Power = 10 W
Input_duration = [100,200,300,400,500,600,700,1000]

duration = [121,175,120,81,101,169,284,775]
chirp_300 = [122,178,122,81.1,97.2,159,265,738]
chirp_1000 = [122,187,129,81.8,88.7,136,222,723]

plt.figure()
plt.plot(Input_duration,duration, '-o',label='no chirp')
# plt.plot(Input_duration,chirp_300, '-o',label='chirp 300 GHz/ps')
# plt.plot(Input_duration,chirp_1000, '-o',label='chirp 1000 GHz/ps')
plt.xlabel('Input pulse duration [fs]')
plt.ylabel('Pulse duration [fs]')
#plt.title()
plt.legend()
plt.grid()