# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 22:02:12 2023

@author: ibaldoni
"""


import numpy as np

S_df = np.array([1e7, 1e6, 1e5, 1e4, 5e3, 5e3, 5e3])
frequencies = np.array([10, 100, 1000, 1e4, 1e5, 1e6, 1e7])
P = 1 # Optical power in watts

S_phi =  S_df/(2 * np.pi * frequencies)**2
L = 10 * np.log10(0.5 * S_phi)

print(L)


import numpy as np

def S_df(frequencies):
    L = np.array([31.02610268, 1.02610268, -28.97389732, -58.97389732, -81.98419728, -101.98419728, -121.98419728])
    frequencies = np.array([10, 100, 1000, 1e4, 1e5, 1e6, 1e7])
    # P = 1e-3 # Optical power in watts
    
    S_phi = 2 * P**2 * 10**(L/10)
    S_df_f = S_phi * (2 * np.pi * frequencies)**2
    return S_df_f

print(S_df)
