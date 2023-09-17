# -*- coding: utf-8 -*-
"""
Created on Mon Oct 29 17:00:40 2018
@author: Rafael Probst
edited for Python 3.7: IB 2022
"""

#####################################
###         USER INPUT            ###
#####################################


import numpy as np

widths  = np.arange(0.1,1.05,0.05)

savePlot = True


for i in widths:
    folder = "./NL_20110305_12 (Spirou)\dispersion"
    input_file = folder + "\scale_{:.2f}.txt".format(i)
    #plot_param = "n_eff (real)"
    plot_param = "dispersion parameter"
    #plot_param = "effective mode area"
    #plot_param = "transmission [dB/m]"
    #plot_param = "transmission (40 mm) [dB]"
    label = "LP$_{01}$ slow"
    
    
    
    #####################################
    ###        PROGRAM CODE           ###
    #####################################
    
    from matplotlib import pyplot as plt
    from scipy import interpolate
    
    plt.figure()
    
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
    D    = [d*1.0e6  for d in D] # convert to ps/nm/km
    a    = [x*1.0e12 for x in a] # convert o µm^2
    
    # plotting
    wl_nm  = [w*1.0e9 for w in wl]
    if plot_param == "dispersion parameter":
        wl_nm2 = wl_nm[1:-1]
        plt.plot(wl_nm2, D, label=label)
        plt.ylabel("group-velocity dispersion [ps/(nm*km)]")
        plt.ylim(np.min(D)-100,np.max(D)+100)
        header = "zero-dispersion: "
        for z in find_zero_disp(wl_nm2, D):
            header += str(round(z,1))+" nm "
    
    if "transmission" in plot_param:
        plt.ylabel(plot_param)
        plt.ylim(-30,5)
        if "(" in plot_param:
            mm = float(plot_param.split("(")[1].split(" mm")[0])
            t = [mm*x/1000. for x in t]
        #plt.plot(wl_nm, t, label=label)
        plt.plot(wl_nm, t, label="LP$_{01}$ slow axis")
            
    if plot_param == "n_eff (real)":
        plt.plot(wl_nm, n, label=label)
        plt.ylabel("effective refractive index (real part)")
    
    if plot_param == "effective mode area":
        plt.plot(wl_nm, a, label=label)
        plt.text.usetex=True
        plt.ylabel("effective mode area [${\mu}m^2$]")
    
    fig = plt.gcf()
    fig.set_size_inches(7, 4.2, forward=True)
    #plt.xlim(350., 1750.)
    plt.xlim(600., 2600.)
    plt.xlabel("wavelength [nm]")
    plt.legend(loc="best", fontsize = 'small')
    plt.title('Width = {:.2f} nm'.format(i))
    plt.grid(True)
    plt.tight_layout()
    if savePlot: plt.savefig(folder +"/" +str(plot_param)+'_{:.2f} nm.png'.format(i))
