# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 16:04:00 2020

@author: ibaldoni

Results from JCM Wave calculations
"""

#####################################
###         USER INPUT            ###
#####################################


import numpy as np
from matplotlib import pyplot as plt
from scipy import interpolate
from scipy import integrate as integ

s = (146,3)
n_eff_l_T = np.zeros(s)



input_file = "scale_1.50_800nm.txt"

#input_file = "C:\\Users\\ibaldoni\\Desktop\\JCM_Results\\20200720_scale_1.50.txt"

#plot_param = 'Nada'
#plot_param = "n_eff (real)"
plot_param = "dispersion parameter"
#plot_param = "effective mode area"
#plot_param = "transmission [dB/m]"
#plot_param = "transmission (55 mm) [dB]"
label = "800 nm height, 1.5 um width (slow)"

#####################################
###        PROGRAM CODE           ###
#####################################

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

def find_zero_disp(wl, D):
    """
    This function finds the zero-dispersion points of the 
    dispersion parameter D as a function of wavelength wl.
    input: list of wavelengths wl, list of dispersion values D
    returns a list of all zero-dispersion values (ZDWs)
    wl can be in m, nm, or other, unit of ZDWs will be accordingly.
    """
    return interpolate.UnivariateSpline(wl,D,s=0).roots()


# main part of the program
data = np.loadtxt(input_file)
shift = 3 if "slow" in label else 0
wl   = [d[0] for d in data]
n    = [d[1+shift] for d in data]
t    = [d[2+shift] for d in data]
a    = [d[3+shift] for d in data]
D    = compute_dispersion_parameter(wl, n)
#    dnl_dT = compute_dnl_dT(wl, n, Temperatures)
D    = [d*1.0e6/(2*np.pi)  for d in D] # convert to ps/nm/km
a    = [x*1.0e12 for x in a] # convert o µm^2


#mu = np.linspace(-10,10,144)
#D = [d/2 for d in D]
#Dint = mu**2 
#plt.plot(Dint)
#Dint = integ.cumtrapz(D,wl[1:-1])
#plt.plot(wl[1:-2], Dint)


c = 299792458
    
l_carrier = 1542   
f_carrier =  c/(l_carrier*1e-9)*1e-12

# plotting
wl_nm  = [w*1.0e9 for w in wl]
if plot_param == "dispersion parameter":
    wl_nm2 = wl_nm[1:-1]
    THz = [1e-3*299792458/w for w in wl_nm2]
    plt.plot(THz, D, label=label)
    plt.vlines(f_carrier, min(D),max(D)+2,linewidth = 0.7, linestyle = '--', color = 'r')
    plt.ylabel("$D_{2}$ / $2\pi$ [ps/(nm*km)]")
    plt.xlabel("Frequency (THz)")


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

#%%
def group_refractive_index(n,wl_nm,dnl_dT, Temperatures):
    
    def compute_dn_dT(n, T):        
        # Para esto, necesito todos los resultados
        dn_dT = [(n[:][i+1]-n[:][i])/(T[i+1]-T[i]) for i in range(len(T)-1)]
     
        return dn_dT
    
    def compute_dnl_dT(wl, n, T):
        
        # Para esto, necesito todos los resultados        
        n1 = [(n[i+1]-n[i])/(wl[i+1]-wl[i]) for i in range(len(wl)-1)]        
        dnl_dT = [(n1[i+1]-n1[i])/(T[i+1]-T[i]) for i in range(len(T)-1)]
        return dnl_dT    
    
#    Lr =  V_eff/A_eff
    
    
    dn      = np.diff(n)
    dlambda = np.diff(wl_nm)
    dn_dl = dn/dlambda
    
    wl_nm1 = wl_nm[:-1]
    n_l = n[:-1]
    
    n_g = n_l - wl_nm1 * dn_dl
    
    
    import matplotlib.pyplot as plt
    from scipy import interpolate

    n_gc = interpolate.interp1d(wl_nm1, n_g)
    
#        dn_dl_int = interpolate.interp1d(wl_nm1, dn_dl)
    c = 299792458
    
    l_carrier = 1542   
    
    frep = 12*1e9 #20 GHz
    
    L = c/(n_gc(l_carrier)*frep)
    print('Length = %s mm' % (L*1000))
    dn_dT = 2.45 * 10**(-5)     # 1/K
    
#        kb = 1.380649e-23           # J/K    
#    k = 1.3806504e-16          # erg/K 

    dn2_dTdl = dnl_dT[0]
    
    dfrep_dT = c/L*(1/n_gc(l_carrier)**2)*(dn_dT - (l_carrier*1e-9)*dn2_dTdl) #* fruta    




    Plots = False
    if Plots==True:
        plt.figure()
        plt.plot(wl_nm, n)
        plt.xlabel('Wavelength (nm)')
        plt.ylabel("$n_{eff}$ $(\lambda)$", fontsize = 17)
        plt.title("Effective refractive index vs. Wavelength")
        plt.legend()
        
        plt.figure()
        plt.plot(wl_nm1, n_g)
        plt.plot(l_carrier,n_gc(l_carrier),'o', label = '%s nm' % l_carrier)
        plt.xlabel('Wavelength (nm)')
        plt.ylabel("$n_{g}$ $(\lambda)$", fontsize = 17)
        plt.title("Group refractive index vs. Wavelength")
        plt.legend()
        
        plt.figure()
        plt.plot(wl_nm1, dfrep_dT)
        plt.xlabel('Wavelength (nm)')
        plt.ylabel("$df_{rep} / dT $ [ Hz / K ]", fontsize = 17)
        plt.title("Repetition rate variation  vs. Wavelength")
        


   
    return wl_nm1, n_g
