# -*- coding: utf-8 -*-
"""
Created on Mon Oct 29 17:00:40 2018
@author: Rafael Probst
"""

#####################################
###         USER INPUT            ###
#####################################

input_file = "scale_1.50_5e-5.txt"
# plot_param = "n_eff (real)"
plot_param = "dispersion parameter"
# plot_param = "effective mode area"
# plot_param = "Volume"
#plot_param = "transmission [dB/m]"
#plot_param = "transmission (55 mm) [dB]"
label = input_file[11:-4]
print(label)
title = "900 nm height, 1.5 um width (slow)"

#####################################
###        PROGRAM CODE           ###
#####################################

import numpy as np
from matplotlib import pyplot as plt
from scipy import interpolate

### function definitions ###

def compute_dispersion_parameter(wl, n):
    """
    wl, n must be lists with vacuum wavelength (in m!)
    and refractive index, respectively.
    Output: dispersion parameter D in s/m^2
    multiply by 1e6 to obtain ps/(km*nm)
    """
    c0 = 299792458.
    # compute first derivative
    wl1 = [(wl[i+1]+wl[i])/2. for i in range(len(wl)-1)]
    n1 = [(n[i+1]-n[i])/(wl[i+1]-wl[i]) for i in range(len(wl)-1)]
    # compute second derivative
    wl2 = [(wl1[i+1]+wl1[i])/2. for i in range(len(wl1)-1)]
    n2 = [(n1[i+1]-n1[i])/(wl1[i+1]-wl1[i]) for i in range(len(wl1)-1)]
    # compute dispersion parameter -> return
    return [(-wl2[i]/c0)*n2[i] for i in range(len(wl2))]

def Volume28(wl,n):
    Vol = [2.8*(-wl[i]/n[i])**3 for i in range(len(wl))]
    
    return Vol
    

def find_zero_disp(wl, D):
    """
    This function finds the zero-dispersion points of the 
    dispersion parameter D as a function of wavelength wl.
    input: list of wavelengths wl, list of dispersion values D
    returns a list of all zero-dispersion values (ZDWs)
    wl can be in m, nm, or other, unit of ZDWs will be accordingly.
    """
    return interpolate.UnivariateSpline(wl,D,s=0).roots()

R = 3.776e-3
# main part of the program
data = np.loadtxt(input_file)
shift = 3 if "slow" in label else 0
wl   = [d[0] for d in data]
n    = [d[1+shift] for d in data]
t    = [d[2+shift] for d in data]
a    = [d[3+shift] for d in data]
D    = compute_dispersion_parameter(wl, n)
D    = [d*1.0e6/(2*np.pi)  for d in D] # convert to ps/nm/km
a    = [x*1.0e12 for x in a] # convert o µm^2

Vol = [(2.8*(wl[i]/n[i])**3)*1e12 for i in range(len(wl))]

V_2  = [2*np.pi*R*aa for aa in a]

# plotting
plt.figure()
wl_nm  = [w*1.0e9 for w in wl]
if plot_param == "dispersion parameter":
    wl_nm2 = wl_nm[1:-1]
    THz = [1e-3*299792458/w for w in wl_nm2]
    plt.plot(wl_nm2, D, label=label,linewidth = 2.0)
    plt.ylabel("$D_{2}$ / $2\pi$ [ps/(nm*km)]")
    plt.xlabel("Frequency (THz)")
    plt.grid()
    plt.legend()
    plt.title(title)

if plot_param == "Volume":
    plt.plot(wl_nm, Vol, label=label, color ='r')
#    plt.plot(wl_nm, V_2, label=label, color ='b')    
    plt.ylabel("effective Volume")
    plt.xlabel("Wavelength (nm)")

    # plt.ylim(-1500,175)
    # header = "zero-dispersion: "
    # for z in find_zero_disp(wl_nm2, D):
    #     header += str(round(z,1))+" nm "

if "transmission" in plot_param:
    plt.ylabel(plot_param)
    plt.ylim(-30,5)
    if "(" in plot_param:
        mm = float(plot_param.split("(")[1].split(" mm")[0])
        t = [mm*x/1000. for x in t]
    plt.plot(wl_nm, t, label=label)
    plt.xlabel("Wavelength (nm)")
        
if plot_param == "n_eff (real)":
    plt.plot(wl_nm, n, label=label)
    plt.ylabel("effective refractive index (real part)")
    plt.xlabel("Wavelength (nm)")

if plot_param == "effective mode area":
    plt.plot(wl_nm, a, label=label)
    plt.text.usetex=True
    plt.ylabel("effective mode area [${\mu}m^2$]")
    plt.xlabel("Wavelength (nm)")

plt.savefig(str(input_file[:-4])+str(plot_param)+'.png')