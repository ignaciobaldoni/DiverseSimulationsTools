# -*- coding: utf-8 -*-
"""
Created on Mon Oct 29 17:00:40 2018
@author: Rafael Probst
"""

#####################################
###         USER INPUT            ###
#####################################






heights = [880]
for i in heights:
    
    import os
    directory_data = os.getcwd()+'/dispersion_h'+str(i)+'nm'  
    path, dirs, files = next(os.walk(directory_data))
    
    Mx_before = 0
    for j in files:
        
        width=j[6:-4]
        input_file = directory_data+"/"+j
        
        #plot_param = "n_eff (real)"
        plot_param = "dispersion parameter"
        #plot_param = "effective mode area"
        #plot_param = "transmission [dB/m]"
        #plot_param = "transmission (55 mm) [dB]"
        label = str(i)+" nm height, "+str(width)+" um width (slow)"
        xLabel = 'wavelength' #or 'frequency'
        center_frequency = 1542.14 # [nm]
        
        #####################################
        ###        PROGRAM CODE           ###
        #####################################
        
        import numpy as np
        from matplotlib import pyplot as plt
        plt.rc('xtick', labelsize=17) 
        plt.rc('ytick', labelsize=17) 
        plt.rcParams.update({'font.size': 17})
        plt.rcParams['figure.figsize'] = (16, 8)
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
        D    = [d*1.0e6/(2*np.pi)  for d in D] # convert to ps/nm/km
        a    = [x*1.0e12 for x in a] # convert o µm^2
        
        # plotting
        wl_nm  = [w*1.0e9 for w in wl]
        if plot_param == "dispersion parameter":
            wl_nm2 = wl_nm[1:-1]
            THz = [1e-3*299792458/w for w in wl_nm2]
            if xLabel == 'frequency': plt.plot(THz, D, label=label) 
            if xLabel != 'frequency': plt.plot(wl_nm2, D, label=label)
            plt.ylabel("$D_{2}$ / $2\pi$ [ps/(nm*km)]")
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
                
        if plot_param == "n_eff (real)":
            plt.plot(wl_nm, n, label=label)
            plt.ylabel("effective refractive index (real part)")
        
        if plot_param == "effective mode area":
            plt.plot(wl_nm, a, label=label)
            plt.text.usetex=True
            plt.ylabel("effective mode area [${\mu}m^2$]")
        
        
        Max = max(Mx_before,int(np.max(D)))
        Min = 0
        fig = plt.gcf()
        #fig.set_size_inches(7, 4.2, forward=True)
        if xLabel == 'frequency': plt.xlim(170., 250.)
        if xLabel != 'frequency': plt.xlim(1000., 2000.)
        plt.ylim(Min, Max+5)
        plt.yticks(np.arange(Min, Max+5, 2))
        if xLabel == 'frequency': plt.xlabel("frequency [THz]") 
        if xLabel != 'frequency': plt.xlabel("wavelength [nm]")
        plt.legend(loc="lower right", fontsize = 'small')
        plt.grid(True)
        #plt.tight_layout()
        plt.show()
        Mx_before = Max
        
        print(wl_nm2[67],D[67])
        print(wl_nm2[68],D[68])
        
    if center_frequency: plt.vlines(center_frequency, Min, Max+5, 'r',
                                    alpha=0.25, label='Center frequency = %s'%center_frequency)
    plt.title('Dispersion parameter for height = %s nm' % i)
    #plt.savefig('Height = %s nm.png' % i)
    #plt.close()
    
    width = np.arange(2.25,2.5,.05)
    d2_1540 = [45.26024265090844,45.12863346157875,44.39649220986762,44.30782695239805,43.46550468037562]
    d2_1550 = [43.413373521870334,42.48716115684584,42.466113246687286,41.47622100680247,41.48164233206151]
    plt.figure()
    plt.plot(width,d2_1540,'o-',label = '1540 nm')
    plt.plot(width,d2_1550,'o-',label = '1550 nm')        
    plt.grid(True)
    plt.ylabel("$D_{2}$ / $2\pi$ [ps/(nm*km)]")
    plt.xlabel("Width [${\mu}$m]")
    plt.title('Dispersion parameter at specific wavelenths for different widths')
    plt.legend()