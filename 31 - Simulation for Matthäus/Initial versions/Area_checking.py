# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 17:02:16 2022

@author: ibaldoni
"""

import numpy as np
# from scipy.integrate import dblquad
import matplotlib.pyplot as plt

plt.rc('xtick', labelsize=17) 
plt.rc('ytick', labelsize=17) 
plt.rcParams.update({'font.size': 17})
plt.rcParams['figure.figsize'] = (14, 10)

# from matplotlib import cm
# from mpl_toolkits.mplot3d import Axes3D
# from math import pi
from scipy import integrate
# import time

plot3d = True
plot2d = True

def volume_Gaussian(amplitude, xy0, sigma_xy):
    return 2*np.pi*amplitude*sigma_xy[0]*sigma_xy[1]

def gaussian2D(x, y, amplitude, xy0, sigma_xy):
    x = x - xy0[0]
    y = y - xy0[1]
    return (amplitude*np.exp( -x**2/2/sigma_xy[0]**2 - y**2/2/sigma_xy[1]**2 ))**2                                                          

N = 600
X = np.linspace(-2, 2, N)
Y = np.linspace(-2, 2, N)
X, Y = np.meshgrid(X, Y)

# coord4 = 1.6
 # Gaussians
Z = 1+np.exp(-(X**2+Y**2))
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_surface(X, Y, Z)

coord2 = .0
coord1 = coord2
coord3 = 1.5

 # Gaussians
Z2 = coord1*X+coord2*Y+coord3
# fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_surface(X, Y, Z2)    
        
f = lambda z, y, x: coord1*x+coord2*y+coord3
Area1 = integrate.tplquad(f, -0.25,.25, 
                            lambda x: -.25,
                          lambda x: .25,
                          lambda x, y: 0, lambda x, y: 1+np.exp(-(x**2+y**2)))

AREA = Area1[0]
area = np.round(Area1[0],3)

print('Overlap:',area)