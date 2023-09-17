# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 10:38:59 2020

@author: ibaldoni
NO TEMPERATURE

Results from JCM Wave calculations
"""

#####################################
###         USER INPUT            ###
#####################################


Temperatures = [300]

import numpy as np
from matplotlib import pyplot as plt
from scipy import interpolate

s = (146,3)
n_eff_l_T = np.zeros(s)
c = 299792458

for jj in range(0,len(Temperatures)):
    
    input_file = '//menloserver//MFS//99-Data_Warehouse//01-User_Folders-Private//i.baldoni//Arbeit//Python codes//Ignacio//JCM_Results//scale_1.50_300.txt'
    #"C:\\Users\\ibaldoni\\Desktop\\Python Codes\\Ignacio\\JCM_Results\\scale_1.50_"+str(Temperatures[jj])+".txt"

    #input_file = "C:\\Users\\ibaldoni\\Desktop\\JCM_Results\\20200720_scale_1.50.txt"
    
    plot_param = "n_eff (real)"
    #plot_param = "dispersion parameter"
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
    
    
    n_eff_l_T[:,jj] = n
    
    
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
    #
        # plt.ylim(-1500,175)
    #    header = "zero-dispersion: "
    #    for z in find_zero_disp(wl_nm2, D):
    #        header += str(round(z,1))+" nm "
    
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
    
    
    def group_refractive_index(n,wl_nm, Temperatures):
        
        def compute_dn_dT(n, T):        
            # Para esto, necesito todos los resultados
            dn_dT = [(n[:][i+1]-n[:][i])/(T[i+1]-T[i]) for i in range(len(T)-1)]
         
            return dn_dT
        
        def compute_dnl_dT(wl, n, T):
            
            # Para esto, necesito todos los resultados        
            n1 = [(n[i+1]-n[i])/(wl[i+1]-wl[i]) for i in range(len(wl)-1)]        
            dnl_dT = [(n1[i+1]-n1[i])/(T[i+1]-T[i]) for i in range(len(T)-1)]
            return dnl_dT            
        
        dn      = np.diff(n)
        dlambda = np.diff(wl_nm)
        dn_dl = dn/dlambda
        
        wl_nm1 = wl_nm[:-1]
        n_l = n[:-1]
        
        n_g = n_l - wl_nm1 * dn_dl
        
        
        import matplotlib.pyplot as plt
        from scipy import interpolate
    
        n_gc = interpolate.interp1d(wl_nm1, n_g)
        
        
        
##        dn_dl_int = interpolate.interp1d(wl_nm1, dn_dl)
#        
#        
#        l_carrier = 1542   
#        
#        frep = 12*1e9 #20 GHz
#        
#        L = c/(n_gc(l_carrier)*frep)
#        print('Length = %s mm' % (L*1000))
#        dn_dT = 2.45 * 10**(-5)     # 1/K
#        
##        kb = 1.380649e-23           # J/K    
#    #    k = 1.3806504e-16          # erg/K 
#    
#        dn2_dTdl = dnl_dT[0]
        l_carrier = 1542   
#        
        frep = 12.13*1e9 #20 GHz
#        
        L = c/(n_gc(l_carrier)*frep)
        print(L)
#        
#
#
#
#    
    
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

        return n_gc #return wl_nm1, n_g  
        
    
    
    
n_gc = group_refractive_index(n,wl_nm, Temperatures)

l_carrier = 1542   
frep = 12.13*1e9 #20 GHz
L = 0.011831550848885661

l_var = np.linspace(1541,1543,105)
f_rep_var = c/(n_gc(l_var)*L)

l_var_plot = c/(l_var*1e-9)*1e-12

fig, ax1 = plt.subplots()

ax1.plot(l_var_plot, (frep - f_rep_var)*1e-6,'o')
ax1.set_ylabel('Repetition rate variation (MHz)')
ax1.set_xlabel('Wavelength (THz)')
ax2 = ax1.twinx()
ax2.set_ylabel('Group index')  # we already handled the x-label with ax1
ax2.plot(l_var_plot, n_gc(l_var))
ax2.tick_params(axis='y')

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()

plt.suptitle('Repetition rate variation due to changes in the group index')    