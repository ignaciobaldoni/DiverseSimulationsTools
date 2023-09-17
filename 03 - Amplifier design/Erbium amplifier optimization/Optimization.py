# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 13:11:24 2022

@author: ibaldoni
"""

import matplotlib.pyplot as plt
import pandas as pd

plt.rc('xtick', labelsize=17) 
plt.rc('ytick', labelsize=17) 
plt.rcParams.update({'font.size': 15})
plt.rcParams['figure.figsize'] = (27, 16)

import numpy as np

sim_folder = 'sim_results\\'
'''
 This should be further automated. 
 x = Stage 1
 y = Stage 2
 
'''

# x = np.arange(1,2.6,0.1)
# y = np.arange(0.5,2.1,0.1)
# fileName = sim_folder+'5dB_22dB.dat'

# x = np.arange(0.25,1.65,0.1)
# y = np.arange(1,8.26,0.25)
# fileName = sim_folder+'22dB_5dB.dat'

# x = np.arange(0.25,1.61,0.1)
# y = np.arange(0.25,2.61,0.1)
# fileName = sim_folder+'22dB_22dB.dat'

# x = np.arange(0.25,1.65,0.1)
# y = np.arange(1,8.6,0.5)
# fileName = sim_folder+'5dB_5dB.dat'

x = np.arange(0.1,2.2,0.1) # Stage 1
y = np.arange(0.1,2.2,0.1) # Stage 2
fileName = sim_folder+'22dB_22dB_v2.dat'


Sim = pd.read_table(fileName,sep=',',names=['Stage1','Stage2','SNR',
                                            'Power','Pulse_duration'])

sim = pd.DataFrame(Sim)


Optim_power = 20
rslt_df = sim[(sim['SNR'] > 18.2135) &
          (sim['Power']>Optim_power*1e-3)]

print('\nResult dataframe :\n', rslt_df)


z = np.array(sim.SNR)
Power_ = np.round(np.array(sim.Power)*1e3,1)

Power = Power_.reshape(len(x), len(y))

Z = z.reshape(len(x), len(y))
plt.imshow(Z, interpolation='bilinear')


for (j,i),label in np.ndenumerate(Power):
    plt.text(i,j,'('+str(label)+')',fontsize = 14,ha='center',va='top')
# for (j,i),label in np.ndenumerate(Pulse):
#     plt.text(i,j,str(label),ha='center',va='center')
for (j,i),label in np.ndenumerate(Z):
    plt.text(i,j,str(label),fontsize = 14,ha='center',va='bottom')
    


positions = np.arange(0,len(y))
labels = np.round(y,1)
plt.xticks(positions, labels)

positions = np.arange(0,len(x))
labels = np.round(x,1)
plt.yticks(positions, labels)
plt.xlabel('Fiber length (2nd Stage) [m]')
plt.ylabel('Fiber length (1st Stage) [m]')
plt.title('SNR '+str(fileName[:-4]))
plt.tight_layout()
    
cbar = plt.colorbar()

for index, row in rslt_df.iterrows():
    b = np.abs(x - row['Stage1'])<0.05
    r = np.array(range(len(b)))
    c = np.abs(x - row['Stage2'])<0.05
    rr = np.array(range(len(b)))
    pts = np.array([r[b][0],rr[c][0]])
    plt.scatter(pts[1], pts[0], marker=".", color="red", s=200)
    
#plt.savefig('SNR '+str(fileName[:-4])+'.png')
plt.show()

