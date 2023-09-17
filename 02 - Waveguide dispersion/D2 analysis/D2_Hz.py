# -*- coding: utf-8 -*-
"""
Created on Mon Oct 29 17:00:40 2018
@author: Rafael Probst
"""

#####################################
###         USER INPUT            ###
#####################################


heights = [870,880,890]

D2_1550_880 = []    
D2_1540_880 = []

D2_1550_870 = []    
D2_1540_870 = []

D2_1550_890 = []    
D2_1540_890 = []


for i in heights:
    
    
    
    import os
    directory_data = os.getcwd()+'/dispersion_h'+str(i)+'nm'  
    path, dirs, files = next(os.walk(directory_data))
    
    Mx_before = 0
    for j in files:
        
        width=j[6:-4]
        input_file = directory_data+"/"+j
        
        #plot_param = "n_eff (real)"

        plot_param = 'Microresonator D2 parameter'

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
        
        def n_g1(n,wl_nm,w0):
            
        #    Lr =  V_eff/A_eff
            
            
            dn      = np.diff(n)
            dlambda = np.diff(wl_nm)
            dn_dl = dn/dlambda
            
            wl_nm1 = wl_nm[:-1]
            n_l = n[:-1]
            
            n_g = n_l - wl_nm1 * dn_dl
            
            n_gc = interpolate.interp1d(wl_nm1, n_g)
            
            print('ng = %s'%(n_gc(w0*1e-9)))

            return n_gc(w0*1e-9)
            
                 
        D1 = 2*np.pi*12.126e9 #[Hz]
             
        def beta1(n,wl,w0):
             c0 = 299792458.                     #[m/s]
             #omega0 = center_frequency*1e-9 #[Hz]
             return (c0/n_g1(n,wl,w0))**(-1)
        
        def beta2(wl,n):
            """
            wl, n must be lists with vacuum wavelength (in m!)
            and refractive index, respectively.
            Output: beta2 in s^2/m

            """
            c0 = 299792458.                     #[m/s]    
            Disp = np.array(compute_dispersion_parameter(wl, n))
            wl = np.array(wl)
            #print(Disp*(wl[1:len(wl)-1])**2)
            return -(wl[1:len(wl)-1])**2/(2*np.pi*c0)*Disp
        # ((wl**2)/(2*np.pi*c0))*compute_dispersion_parameter(wl, n)
        
        def compute_D2(wl,n):
            return -beta2(wl,n)*D1**2/beta1(n,wl,center_frequency)
        
        MR_d2 = compute_D2(wl, n)*1e-3
        
        if i == 880: D2_1550_880.append(MR_d2[68])
        if i == 880: D2_1540_880.append(MR_d2[67])

        if i == 870: D2_1550_870.append(MR_d2[68])
        if i == 870: D2_1540_870.append(MR_d2[67])
        
        if i == 890: D2_1550_890.append(MR_d2[68])
        if i == 890: D2_1540_890.append(MR_d2[67])
        
        
        wl_nm  = [w*1.0e9 for w in wl]
        if plot_param == "Microresonator D2 parameter":
            
            wl_nm2 = wl_nm[1:-1]
            THz = [1e-3*299792458/w for w in wl_nm2]
            if xLabel == 'frequency': plt.plot(THz, MR_d2, label=label) 
            if xLabel != 'frequency': plt.plot(wl_nm2, MR_d2, label=label)
            plt.ylabel("$D_{2}$ [kHz]")
            if xLabel == 'frequency': plt.xlabel("frequency [THz]") 
            if xLabel != 'frequency': plt.xlabel("wavelength [nm]")
            
            

        if plot_param == "n_eff (real)":
            plt.plot(wl_nm, n, label=label)
            plt.ylabel("effective refractive index (real part)")
        
        if plot_param == "effective mode area":
            plt.plot(wl_nm, a, label=label)
            plt.text.usetex=True
            plt.ylabel("effective mode area [${\mu}m^2$]")
        
    if center_frequency and i==800: plt.vlines(center_frequency, -300, 300, 'r',
                                    alpha=0.25, label='Center frequency = %s'%center_frequency)
    
    if center_frequency and i>=850: plt.vlines(center_frequency, 0, 350, 'r',
                                    alpha=0.25, label='Center frequency = %s'%center_frequency)
    
    
    
    
    plt.legend()
plt.grid()

plt.figure()
plt.plot(D2_1550_890,'o-',label='h = 890nm, 1550nm')
plt.plot(D2_1540_890,'o-',label='h = 890nm, 1540nm')
plt.plot(D2_1550_880,'o-',label='h = 880nm, 1550nm')
plt.plot(D2_1540_880,'o-',label='h = 880nm, 1540nm')
plt.plot(D2_1550_870,'o-',label='h = 870nm, 1550nm')
plt.plot(D2_1540_870,'o-',label='h = 870nm, 1540nm')
plt.legend()
plt.grid()