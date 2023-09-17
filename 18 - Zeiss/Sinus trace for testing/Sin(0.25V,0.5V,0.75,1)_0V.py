# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 16:40:36 2023

@author: ibaldoni
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

## File from Zeiss as example
# File = '0002_Hysterese_Bias5A5f0_1.csv'

# Original = pd.read_csv(File,sep=',',names=['time','V'])
# Original=Original.drop(Original.index[0])  

# t = Original.time
# t = t.astype(float)
# V = Original.V
# V = V.astype(float)

# plt.plot(t,V,label='Zeiss')

### Trajcetory definition

t = np.arange(0,6500,.1)

for i in range(0,len(t)):
    t[i] = np.round(t[i],2)

t1 = 500
t2 = 2000
t3 = 3500
t4 = 5000
t11 = 754

freq = .1
bias1 = 0
bias = 5
t_Add = 250

A1 = .25/2
A2 = .5/2
A3 = .75/2
A4 = 1.00/2

funt = (A1+bias1+A1*np.cos(freq*t+np.pi-freq*t1))*np.heaviside(t-t1,1)*np.heaviside(-t+t1+t11,1)
funt2 = (A2+bias1+A2*np.cos(freq*t+np.pi-freq*t2))*np.heaviside(t-t2,1)*np.heaviside(-t+t2+t11,1)
funt3 = (A3+bias1+A3*np.cos(freq*t+np.pi-freq*t3))*np.heaviside(t-t3,1)*np.heaviside(-t+t3+t11,1)
funt4 = (A4+bias1+A4*np.cos(freq*t+np.pi-freq*t4))*np.heaviside(t-t4,1)*np.heaviside(-t+t4+t11,1)

FUN_5V = funt + funt2 + funt3 + funt4 + bias
plt.plot(t,FUN_5V,'-'); plt.ylabel('Voltage [V]'); plt.xlabel('Time [s]')

d = {'time': t, 'voltage': FUN_5V}
df = pd.DataFrame(data=d)
df.to_csv("SMILE_Trajectory_bias5V.csv",index=False)

FUN_5V_0V = (funt + funt2 + funt3 + funt4 + bias)*np.heaviside(t-t1+t_Add,1)*np.heaviside(-t+t4+t_Add+t11,1)
plt.figure()
plt.plot(t,FUN_5V_0V,'-'); plt.ylabel('Voltage [V]'); plt.xlabel('Time [s]')

d = {'time': t, 'voltage': FUN_5V_0V}
df = pd.DataFrame(data=d)
df.to_csv("SMILE_Trajectory_bias5V_0V.csv",index=False)