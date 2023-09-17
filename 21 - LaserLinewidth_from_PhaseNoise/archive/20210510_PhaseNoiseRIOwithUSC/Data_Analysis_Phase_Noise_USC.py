# -*- coding: utf-8 -*-
"""
Created on Wed Apr 28 18:16:48 2021

@author: ibaldoni


PXA Phase noise measurement read out from USB Stick

"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import datetime as dt
import os.path

#%% 
plt.close('all')

skipRow = 59
figSize = (14,10)
Headers = ['1','2','3','4','5','6']

Phase_noise_files       = True
Phase_noise_in_radians  = False

phase_noise_at_ten_k = []

import os
directory_data = os.getcwd()  
path, dirs, files = next(os.walk(directory_data))




#%% Phase noise files
if Phase_noise_files == True:
    exceptions_not = ['.py','.png']
    exceptions_in = ['.txt']
    
  
    for i in files:

        if (all((ele not in i) for ele in exceptions_not) and
            all((ele_in in i) for ele_in in exceptions_in)):
            print(i)
            
            fileName =  i
            
            phase_noise = pd.read_csv(fileName,skiprows = skipRow, names = Headers)

            phase_noise = phase_noise.drop(['3','4','5','6'], 1)

            
            phase_noise = phase_noise.reset_index()
            print(phase_noise)
            phase_noise = phase_noise.drop(['level_0','level_1'], 1)

            
            Frequency   = phase_noise.iloc[::2]
            Frequency   = Frequency[1:]
            PSD         = phase_noise.iloc[1::2]

#            Frequency = phase_noise['2']
#            PSD = phase_noise['1']

            

            Name = ''
            Xlabel = 'PSD (dBc/Hz)'
            
            if Phase_noise_in_radians == True: 
                PSD_to_rad = np.sqrt(10**(2*PSD/10))
                #np.sqrt(np.sqrt(2*(10**(PSD/10))))
                PSD = PSD_to_rad*1e6 # urad
                Name = 'uradians' 
                Xlabel = 'PSD (urad/sqrt(Hz))'
     
#            tenK        = np.where(Frequency == 10000)[0][0]
#            PN_tenK     = np.round(PSD.iloc[tenK],2)
#            PN_tenK_avg = PSD.iloc[tenK]            
#            tenK_plot   = Frequency.iloc[tenK]
            
            PN_tenK = 'Unknown'
            Ylim = [-100,0]
            
            phase_noise_at_tenK = str(PN_tenK)+' dBc/Hz @ 10 kHz'
            
            plt.figure(figsize=figSize)
            plt.plot(PSD,Frequency,label = 'RIO Laser')
#            plt.plot(PSDDD,Freqeqe, label ='USC')
            
#            plt.plot(10000,-110,'ro')
            #plt.plot(tenK_plot,PN_tenK,'o',tenK_plot,PN_tenK_avg,'o')
            
            plt.ylabel(Xlabel,fontsize = 20)
            plt.xlabel('Frequency (Hz)', fontsize = 20)
            plt.ylim(Ylim)
            plt.tick_params(labelsize=17)
#            plt.text(10000,-40,phase_noise_at_tenK,fontsize = 17,
#                     color='black', bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=1'))
            plt.title(str(i), fontsize = 17)
            plt.xscale('log')
            plt.grid(True, which="both")
            plt.legend()
            
#            plt.savefig(fileName+'_USC_'+str(Name)+'.png')
#            plt.close()
            
            
#Freqeqe = Frequency
#PSDDD = PSD 

