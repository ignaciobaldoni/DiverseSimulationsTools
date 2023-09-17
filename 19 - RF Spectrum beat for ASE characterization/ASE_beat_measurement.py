# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 12:38:55 2022

@author: ibaldoni
"""

import numpy as np
import matplotlib.pyplot as plt

plt.close('all')
import pandas as pd

plt.rc('xtick', labelsize=17) 
plt.rc('ytick', labelsize=17) 
plt.rcParams.update({'font.size': 17})
plt.rcParams['figure.figsize'] = (14, 10)

import warnings
warnings.filterwarnings("ignore")

import scipy.integrate as integrate
import scipy.special as special
import time as dt

header_list = ["Frequency", "Trace"]
skipRow = 61#29

import os
directory_data = os.getcwd()+'\\raw_data'  
path, dirs, files = next(os.walk(directory_data))


exceptions_in = ['dBµV_Amplified_NOT_DWDM','dBµV_HigherASE','Beat Note']
p=0
for i in files[:]:
    if any((ele_in in i) for ele_in in exceptions_in):
        fileName =  directory_data+'\\'+str(i)
        print(fileName)
        # seed_input = fileName[14:23]


        med=pd.read_csv(fileName,skiprows = skipRow, names=header_list)
        med.Trace = med.Trace/np.max((np.abs(med.Trace)))
        # if 'dBµV' in i: med.Trace = med.Trace-10*np.log(50)-65
        
        center = np.argmax(med.Trace)
        med.Frequency = med.Frequency - med.Frequency[center]
        plt.plot(med.Frequency*1e-6,med.Trace,label=str(i),alpha=1*(1-p*0.3))
        # if 'PhaseNoise' in fileName:
            
        #     plt.plot(np.log10(med.Frequency),med.Trace,label=str(i))
        # else:
        #     plt.plot(med.Frequency,med.Trace,label=str(i))
        plt.legend()
        p+=1
        
plt.grid()
plt.ylabel('Power [dBm]')        
plt.xlabel('Frequency [MHz]')        


header_list = ["wavelength", "PSD"]
skipRow = 3

OSA = 'MPQ'
OSA_Calibration = 3


exceptions_in = ['ASE3']

p=0
plt.figure()
for i in files[:]:
    if any((ele_in in i) for ele_in in exceptions_in):
        fileName = 'raw_data\\'+ str(i)
        # print(fileName)
        seed_input = fileName[14:23]


        med=pd.read_csv(fileName,skiprows = skipRow, names=header_list)
        
        med = med[:1001]
    
        wavelength      = med['wavelength']
        
        wavelength = wavelength.astype(float)
        PSD             = med['PSD']
        PSD = PSD.astype(float)
        
        PSD = PSD + OSA_Calibration # For Severus: 13.4 
        fit_spectrum = True
        spectrum_plot = True            
                    
        if fit_spectrum == True:      
    #        cw laser contribution calculation 
            lin = 10**(PSD/10) # Power originally is on dBm/nm    
            
            # PSD = 10*np.log10(lin/np.max(lin))
            
            
            total=integrate.simps(lin,wavelength)

            total = np.round(total,3)

            print('Measured area:',total,'mW')
            
            

        if spectrum_plot == True:
            # plt.figure()
            
            plt.plot(wavelength,PSD,alpha=0.63)#,label='Power = %s mW'%total)
            print('Power measured comb = %s mW'%(total))
            
            plt.xlabel('Wavelength [nm]')
            plt.ylabel('Power [dBm/nm]')
           
            
            
            
            
            #if p==1: plt.title('ASE SPIROU Amplifier dif=%s dB' % Value)
            
            # plt.legend()
            Ylim = [-50,15]
            plt.ylim(Ylim)
            plt.xlim([1500,1600])
            p+=1
plt.grid()
            # plt.savefig('Spirou ASE 430 µW seed.png')
            # plt.savefig('Spirou ASE 10 mW seed.png')
# plt.savefig('Spirou ASE'+str(seed_input)+'.png')
# plt.savefig('Spirou ASE'+str(Opt_power)+'.png')
# plt.savefig('Comparison'+str(exceptions_in[0])+'and_low_seed.png')
    
          