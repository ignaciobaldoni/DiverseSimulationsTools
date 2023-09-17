# -*- coding: utf-8 -*-
"""
Created on Mon Oct 25 09:55:52 2021

@author: ibaldoni
"""

import numpy as np
import matplotlib.pyplot as plt

c = 299792458.0

w0 = 286e12     #1050e-9
wmu = 282e12    #542e-9

w00 = w0-wmu

D1 = 25.718e9*(2*np.pi)
D2 = 30.292e3*(2*np.pi)/2
D3 = -190.88*np.pi*2/6 

def w_mu(mu):
    
    return w0+D1*mu+D2*mu**2+D3*mu**3

# Enter the coefficients of the poly 
# in the array
coeff = [D3, D2, D1, w00]
print(np.roots(coeff))

mu = np.linspace(-100000,100000,1000000)
plt.plot(mu,w_mu(mu))