# -*- coding: utf-8 -*-
"""
Created on Tue Sep  5 17:33:44 2023

@author: ibaldoni
"""

import random
import matplotlib.pyplot as plt
import numpy as np
import allantools as at

def random_spike_generator(duration=10000, spike_interval=50, spike_height = 5):
    data = []
    for t in range(duration):
        if t % spike_interval == 0:
            # Create a spike. The magnitude and direction (up or down) is random.
            spike = random.choice([spike_height, -spike_height])*1e-12
            data.extend([spike] * spike_interval)
        else:
            # Append a random value between -1 and 1 for noise.
            noise = random.uniform(-1, 1)*1e-12
            data.append(noise)

    # Plot the data
    plt.figure(figsize=(15, 6))
    plt.plot(range(len(data)), data, color='b')
    plt.xlabel('Data Points')
    plt.ylabel('Value')
    plt.title('Random Spike Generator')
    plt.grid(True)
    plt.show()
    
    return data

if __name__ == "__main__":
    Data1 = random_spike_generator(spike_interval=1,spike_height = 1e-12)
    Data2 = random_spike_generator(spike_interval=10)
    Data3 = Data1*np.random.uniform(0,10,len(Data1))



    
    # Assuming you have your data in the 'data' list
    data = np.array(Data1)
    data2 = np.array(Data2)
    data3 = np.array(Data3)
    
    # Calculate the Allan Deviation
    (t2, ad, ade, adn) = at.mdev(data,data_type="freq", taus="decade")
    plt.loglog(t2, ad,'o-',label = 'Noiseless')
    (t2, ad, ade, adn) = at.mdev(data2,data_type="freq", taus="decade")
    plt.loglog(t2, ad,'o-',label = 'Spike every 10 s')
    (t2, ad, ade, adn) = at.mdev(data3,data_type="freq", taus="decade")
    plt.loglog(t2, ad,'o-',label = 'High random noise')
    

    # plt.loglog(t2, ad,'o-')
    plt.xlabel('Averaging Time (tau)')
    plt.ylabel('Allan Deviation')
    plt.title('Allan Deviation Plot')
    plt.grid(True)
    plt.legend()
    # plt.ylim([10e-14,10e-11])
    plt.show()
