# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 17:02:16 2022

@author: ibaldoni
"""

import numpy as np
from scipy.integrate import dblquad
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D


def volume_Gaussian(amplitude, xy0, sigma_xy):
    return 2*np.pi*amplitude*sigma_xy[0]*sigma_xy[1]

def gaussian2D(x, y, amplitude, xy0, sigma_xy):
    x = x - xy0[0]
    y = y - xy0[1]
    return amplitude*np.exp( -x**2/2/sigma_xy[0]**2 - y**2/2/sigma_xy[1]**2 )                                                          

N = 600
X = np.linspace(-3, 3, N)
Y = np.linspace(-3, 4, N)
X, Y = np.meshgrid(X, Y)

mu = np.array([1, 10])
Sigma = np.array([ 1. , 1.])

mu2 = np.array([0, 0])
Sigma2 = np.array([ 1. , 1 ])

amp1 = 1
amp2 = 1


from scipy import integrate
f = lambda y, x: gaussian2D(x, y, amp1, mu, Sigma)
result1 = integrate.dblquad(f, -20, 20, 
                            lambda x: -20,
                          lambda x: 20)

from scipy import integrate
f = lambda y, x: gaussian2D(x, y, amp2, mu2, Sigma2)
result2 = integrate.dblquad(f, -20, 20, 
                            lambda x: -20,
                          lambda x: 20)
result1 = result1[0]
result2 = result2[0]



# Gaussians
Z = gaussian2D(X, Y, amp1, mu, Sigma)
Z2 = gaussian2D(X, Y, amp2, mu2, Sigma2)


# Two Gaussian
def A_minus_B(x, y, argsA, argsB):
    return gaussian2D(x, y, *argsA) - gaussian2D(x, y, *argsB)

argsA = (amp1, mu, Sigma)
argsB = (amp2, mu2, Sigma2)



plt.figure()
plt.contour(X, Y, Z)#, zdir='z')#, offset=-0.15, cmap=cm.viridis)
plt.contour(X, Y, Z2)#, zdir='z')#, offset=-0.15, cmap=cm.viridis)
import numpy as np
from matplotlib import pyplot as plt
from math import pi

u=mu2[0]     #x-position of the center
v=mu2[1]      #y-position of the center
a=2.     #radius on the x-axis
b=2   #radius on the y-axis

t = np.linspace(0, 2*pi, 100)
plt.plot( u+a*np.cos(t) , v+b*np.sin(t),'b')
plt.grid(color='lightgray',linestyle='--')


from scipy import integrate
f = lambda y, x: 1/result1*gaussian2D(x, y, amp1, mu, Sigma)
result = integrate.dblquad(f, -a+u, a+u, 
                            lambda x: v - np.sqrt(b**2*(1 - ((x-u)/a)**2)),
                          lambda x: v + np.sqrt(b**2*(1 - ((x-u)/a)**2)))

print('gaussian1',result)

u=mu[0]     #x-position of the center
v=mu[1]      #y-position of the center
a=2.     #radius on the x-axis
b=2        #radius on the y-axis



t = np.linspace(0, 2*pi, 100)
plt.plot( u+a*np.cos(t) , v+b*np.sin(t),'r' )
plt.grid(color='lightgray',linestyle='--')


from scipy import integrate
f = lambda y, x: 1/result2*gaussian2D(x, y, amp2, mu2, Sigma2)
result = integrate.dblquad(f, -a+u, a+u, 
                            lambda x: v - np.sqrt(b**2*(1 - ((x-u)/a)**2)),
                          lambda x: v + np.sqrt(b**2*(1 - ((x-u)/a)**2)))

print('gaussian2',result)