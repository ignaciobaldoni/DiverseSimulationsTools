# -*- coding: utf-8 -*-
"""
Created on Thu Jul 11 18:05:08 2019
@author: La Silla MPQ
"""

import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit

def n_SiO2(wl):
    # Sellmeier equation and coefficients for fused silica
    c = [0.6961663, 0.0684043**2, 
         0.4079426, 0.1162414**2, 
         0.8974794, 9.896161**2]
    wl = wl*1.0e6 # conversion into microns
    l2 = wl * wl
    p1 = (c[0]*l2)/(l2-c[1])
    p2 = (c[2]*l2)/(l2-c[3])
    p3 = (c[4]*l2)/(l2-c[5])
    return np.sqrt(1.+p1+p2+p3)


def n_Si3N4(wl):
    # Sellmeier equation and coefficients (fit is reasonably good down to 500 nm)
    c = [3.40789621e+00 ,1.72579862e-02, -4.62288837e-01, 
         1.72570350e-02, 9.26379962e+00, -1.84302988e+04]
    wl = wl*1.0e6 # conversion into microns
    l2 = wl * wl
    p1 = (c[0]*l2)/(l2-c[1])
    p2 = (c[2]*l2)/(l2-c[3])
    p3 = (c[4]*l2)/(l2-c[5])
    return np.sqrt(1.+p1+p2+p3)


def n_Si3N4_fit(wl, c0, c1, c2, c3, c4, c5):
    wl = wl*1.0e6 # conversion into microns
    l2 = wl * wl
    p1 = (c0*l2)/(l2-c1)
    p2 = (c2*l2)/(l2-c3)
    p3 = (c4*l2)/(l2-c5)
    return np.sqrt(1.+p1+p2+p3)


if __name__ == "__main__":
    #print W[:850], N[:850]
    W = np.loadtxt("wavelengths.txt")
    N = np.loadtxt("refractive_index.txt")
    p0 = [0.6961663, 0.0684043**2, 0.4079426, 0.1162414**2, 0.8974794, 9.896161**2]
    p1, pcov = curve_fit(n_Si3N4_fit, W[:600]*1e-9, N[:600], p0=p0)
    #p, pcov = curve_fit(n_Si3N4_fit, W, N, p0=p)
    N2= [n_SiO2(w*1e-9) for w in W]
    #N3= [n_Si3N4_fit(w*1e-9, p1[0],p1[1],p1[2],p1[3],p1[4],p1[5]) for w in W]
    N4= [n_Si3N4(w*1e-9) for w in W[:600]]
    plt.plot(W,N,label="Si3N4")
    plt.plot(W,N2,label="SiO2")
    #plt.plot(W,N3,label="Si3N4 fit")
    plt.plot(W[:600],N4,label="Si3N4 fit")
    plt.legend(loc=0)
    plt.show()
