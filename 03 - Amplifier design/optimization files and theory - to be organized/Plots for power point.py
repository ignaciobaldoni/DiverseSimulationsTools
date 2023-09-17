# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 13:33:17 2022

@author: ibaldoni
"""

import matplotlib.pyplot as plt


plt.rc('xtick', labelsize=17) 
plt.rc('ytick', labelsize=17) 
plt.rcParams.update({'font.size': 17})
plt.rcParams['figure.figsize'] = (8,6)


import numpy as np

# Input = np.arange(0.01,0.260,0.030)
# SNR = [9.34,14.7,16.7,17.8,18.7,19.3,19.8,20.2,20.6]

# plt.plot(Input, SNR,'o-',label = 'SNR\nfixed pump\nfixed length')
# plt.grid()
# plt.legend(loc='lower right')
# plt.ylabel('SNR [dB]')
# plt.xlabel('Seed power [mW]')
# plt.title('-')


# Input = np.arange(50,440,50)
# Output = [7.25,12.8,16.8,19.7,22,23.9,25.4,26.7]
# Gain = [13.7,16.2,17.4,18.1,18.6,18.9,19.2,19.4]

# plt.figure()
# plt.plot(Input, Gain,'o-',label = 'Gain\nfixed seed\nfixed length')
# plt.grid()
# plt.legend(loc='lower right')
# plt.ylabel('Gain [dB]')
# plt.xlabel('Pump power [mW]')
# plt.title('Second Stage')



# Input = np.arange(10,260,30)
# Output = [217,289,302,307,310,311,313,313,314]
# Gain = [3.38,4.56,4.73,4.79,4.83,4.85,4.87,4.88,4.89]

# plt.figure()
# plt.plot(Input, Gain,'o-',label = 'Gain\nfixed seed\nfixed length')
# plt.grid()
# plt.legend(loc='lower right')
# plt.ylabel('Gain [dB]')
# plt.xlabel('Pump power [mW]')
# plt.title('First Stage')


Input = np.arange(0.01,0.260,0.030)
SNR = [10.4,16.4,18.7,20.4,21.5,22.4,23.1,23.8,24.3]

plt.figure()
plt.plot(Input, SNR,'o-',label = 'L = 0.5 m\nPump = 120 mW')
plt.grid()
plt.legend(loc='lower right')
plt.ylabel('SNR [dB]')
plt.xlabel('Seed power [mW]')


Input = np.arange(0.01,0.260,0.030)
SNR = [10.2,16.2,18.6,20.2,21.3,22.2,22.9,23.6,24.1]

plt.plot(Input, SNR,'o-',label = 'L = 0.5 m\nPump = 50 mW')
plt.grid()
plt.legend(loc='lower right')
plt.ylabel('SNR [dB]')
plt.xlabel('Seed power [mW]')
plt.grid()
plt.title('SNR vs. seed power')

Input = np.arange(0.01,0.260,0.030)
SNR = [8.67,14.7,17.1,18.7,19.8,20.7,21.5,22.1,22.7]

plt.plot(Input, SNR,'o-',label = 'L = 1.5 m\nPump = 250 mW')
plt.legend(loc='lower right')
plt.ylabel('SNR [dB]')
plt.xlabel('Seed power [mW]')

