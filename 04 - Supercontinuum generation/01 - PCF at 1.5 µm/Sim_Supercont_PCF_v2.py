# -*- coding: utf-8 -*-
"""
Created on Wed Sep 26 17:14:24 2018
@author: La Silla MPQ
edited for Python 3.7: IB 2022
"""

import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
import pynlo, bisect, os
from scipy.interpolate import UnivariateSpline, interp1d
import pandas as pd  

"""
Usage:
Enter simulation parameters in the user area, run the program.
Wait for simulation to complete. Plot with -simulation results will open.
If saveTXT and/or savePlot is True, find results saved in "results" folder.

"""

opt = []
Taper_Length = np.linspace(0.038, 0.10,10) # [0.030,0.140]
Waist_Scale = [0.25]#np.linspace(0.10, 0.40,10) #[0.3,0.20]
Waist_Length = np.linspace(0.010, 0.0433,10) #[0.3,0.20]
Waist_scale = Waist_Scale
for taper_length in Taper_Length: 
    #for Waist_scale in Waist_Scale:
    for Waist_len in Waist_Length:
        #################################################
        ###                USER AREA                  ###
        #################################################
        
        # Normally, a code should be kept strictly in SI (mks) units. 
        # However, pynlo forces the code to be at least partly in 
        # other units (nm, THz, ps), so we go along with it.
        
        # pulse parameters
        shape       = 'sech2'       # pulse shape ('Gauss' or 'sech2')
        FWHM_ps     = 0.110         # pulse duration [ps]
        carrier_nm  = 1542          # pulse central wavelength [nm]
        frep_MHz    = 12126         # pulse repetition rate [MHz]
        pow_mW      = 8000          # average power [mW]
        GDD         = 0e-3          # linear pre-chirp [ps^2]
        Noise       = False         # Does not seem to change the spectrum at all
        noise_type  = 'sqrt_N_freq' #{one_photon_freq,sqrt_N_freq}
        
        # fiber parameters
        
        untap_len   = 0.010     # initial untapered length [m]
        taper_len   = taper_length     # 0.130 (down-)taper transition length [m]
        waist_len   = Waist_len     # waist length [m]
        uptap_len   = 0.020     # uptaper length [m]
        untap_scale = 1.00      # scale of fiber diam. rel. to original value
        waist_scale = Waist_Scale[0]       # .3 scale of fiber diameter at waist
        
        propa_len   = taper_length + 2*untap_len + waist_len + uptap_len   # 0.20 propagation length [m]
        
        workingDir  = "NL_20110305_12 (Spirou)" # path to working directory
        axis        = "slow"    # polarization ("slow" or "fast")
        
        # simulation parameters
        window_ps   = 10.0          # simulation window [ps]
        steps       = 200           # simulation steps
        points      = 20500        # 13000 simulation points
        low_nm      = 800          # start wavelength [nm]
        high_nm     = 2700          # stop wavelength [nm]
        saveTXT     = False         # save output spectrum as txt?
        savePlot    = saveTXT       # save the plot (will go full screen)
        
        cmap = plt.cm.get_cmap("jet")
        
        #################################################
        ###           FUNCTION DEFINITIONS            ###
        #################################################
        
        # In this section we define functions and the related global variables.
        # Several functions need to be provided to pynlo when doing taper simulations,
        # in order to specify the evolution of dispersion and refractive index
        # (indrectly the evolution of the taper diameter) along the fiber.
        # However, the requirements of pynlo to these functions are such that
        # we cannot pass all parameters that we would like to. This forces us to
        # set them externally with global variables. Actually, this is very bad
        # programming practice, but in this case it is the only way to make certain
        # things possible and to avoid repeting expensive operations with every function call.
        
        # transcribe z-positions (z_stops) of taper sections into list
        lengths = [untap_len, taper_len, waist_len, uptap_len]
        z_stops = [sum(lengths[:i+1]) for i in range(len(lengths))]
        
        # scaling of taper diameter as a function of position z along fiber
        def scale(z):
            if z <= z_stops[0]:              # initial untapered region
                return untap_scale
            if z_stops[0] < z <= z_stops[1]: # down-taper
                return interp1d(z_stops[:2],[untap_scale,waist_scale])(z)
            if z_stops[1] < z <= z_stops[2]: # taper waist
                return waist_scale
            if z_stops[2] < z <= z_stops[3]: # up-taper
                return interp1d(z_stops[2:],[waist_scale,untap_scale])(z)
            if z_stops[3] < z:               # final untapered region
                return untap_scale
        
        # get find the dispersion files
        dispFolder = workingDir + "/dispersion"
        dispFiles = [d for d in os.listdir(dispFolder) if d.endswith(".txt")]    
        # read dispersion files and store as list of arrays
        dispData  = [np.loadtxt(dispFolder+"/"+d) for d in dispFiles]
        # store taper scale for each dispersion file (from dispersion file name)
        dispScale = [float(d.split("_")[1].replace(".txt", "")) for d in dispFiles]
        # get wavelength scale of dispersion data, converted to nm
        dispWL_nm = [w*1e9 for w in dispData[0][:,0]]
        # get the (approx.) list index of the carrier wavelength
        cwl_index = bisect.bisect(dispWL_nm, carrier_nm)
        
        # check: all dispersion files must use the same wavelength scale
        for i in range(1,len(dispData)):
            # all wavelength array elements must have the same value as in the first file
            if not np.all(dispData[i][:,0]==dispData[0][:,0]): # otherwise abort with error
                ValueError("Dispersion files use different wavelength scales")
        
        # switch between slow axis and fast axis
        # if slow axis chosen: delete columns of fast axis from dispersion matrix
        if axis == "slow":
            dispData  = [np.delete(d, np.s_[1:4], axis=1) for d in dispData]
        # in case of fast axis: there is nothing to be done: the first few columns
        # contain the data for the fast axis, which will be used automatically
        elif axis != "fast": # if axis neither "fast" nor "slow": abort with error 
            ValueError("Select either fast or slow polarization axis")
        
        # function for pynlo to specify the evolution of dispersion along fiber
        x = -1  # work-around for next function to fix problem with pynlo
        def dispersion(z):  # z: position along fiber [m]
            if x>=0: z=x    # global parameter x sets z-value externally
            n_eff = []      # dispersion is specified in terms of effective index
            for i in range(len(dispWL_nm)):
                N = [dispData[j][i,1] for j in range(len(dispData))]
                k = min([3,len(dispData)-1]) # order of spline
                n = interp1d(dispScale,N,kind=k)(scale(z))
                n_eff.append(float(n))
            return np.array(dispWL_nm), np.array(n_eff)
        
        # function returns zero-dispersion wavelengths (ZDWs)
        def get_ZDW(wl_nm, D): # input D: dispersion parameter
            return UnivariateSpline(wl_nm,D).roots()
        
        # gamma: function that returns nonlinear coefficient [1/(W*m)]
        def gamma(z):       # z [m]: position along the fiber
            n2 = 2.74e-20   # n2 for fused silica [m^2/W]
            # A (list): effective mode area (at carrier wl) for different taper scales
            A = [dispData[j][cwl_index,3] for j in range(len(dispData))]
            # a: effective mode for scale at given position z along the fiber
            a = interp1d(dispScale,A)(scale(z))
            # compute nonliner coefficient from effective mode area a
            return 2*np.pi*n2/(carrier_nm*1e-9*a)
        
        def dB(num):
            # convert to dB
            return 10 * np.log10( num )
        
        def power(num):
            # convert amplitude (electric field) to power value
            return np.absolute( num )**2
        
        
        #################################################
        ###      PHYSICS SIMULATION WITH PYNLO        ###
        #################################################
        
        # Much of the code in this section has been developed on the basis of the
        # pnlo example "simple_SSFM.py", see pynlo folder or website:
        # https://pynlo.readthedocs.io/en/latest/example_simple.html
        # Additional functionality was added using the pynlo documentation:
        # https://pynlo.readthedocs.io/en/latest/pynlo.html
        
        # create the pulse (power doesn't matter because pulse energy will be set later)
        # if desired the input pulse can be chirped linearly/quadradically: GDD [ps^2], TOD [ps^3]
        if shape == "sech2":
            pulse = pynlo.light.DerivedPulses.SechPulse(power = 1, T0_ps = FWHM_ps/1.76,
                    center_wavelength_nm = carrier_nm, time_window_ps = window_ps, GDD = GDD,
                    TOD = 0, NPTS = points)
        elif shape == "Gauss":
            pulse = pynlo.light.DerivedPulses.GaussianPulse(power = 1, T0_ps = FWHM_ps,
                    center_wavelength_nm = carrier_nm, time_window_ps = window_ps, GDD = GDD,
                    TOD = 0, NPTS = points)
        
        # set the pulse energy
        pulse.set_epp((pow_mW*1.0e-3)/(frep_MHz*1.0e6))
        
        if Noise: 
            print('Noise implemented as %s'%noise_type) 
            pulse.add_noise(noise_type=noise_type)
        
        # create the fiber (betas and gamma are given but will be overwritten thereafter)
        fiber1 = pynlo.media.fibers.fiber.FiberInstance()
        fiber1.generate_fiber(propa_len, center_wl_nm=carrier_nm, betas=(-1.1e+01,-3.2e-02,2.3e-04),
                              gamma_W_m=0.5, gvd_units='ps^n/km', gain=-0.0) # gain = atten. coeff. [1/m]
        # supply external dispersion function dispersion(z), overwrites betas from previous step
        fiber1.set_dispersion_function(dispersion, dispersion_format='n')
        # provide external function gamma(z): specifies nonlinear coefficient along fiber
        fiber1.set_gamma_function(gamma)
        
        # Set up nonlinear interactions
        # local_error makes different result, smaller number slows down the simulation 
        # when the broadening happens.i.e. goes to higher resolution.
        evol = pynlo.interactions.FourWaveMixing.SSFM.SSFM(local_error=0.0001, USE_SIMPLE_RAMAN=False,
               disable_Raman = False, disable_self_steepening  = False)
        # compute pulse propagation
        Y, AW, AT, pulse_out = evol.propagate(pulse_in=pulse, fiber=fiber1, n_steps=steps, 
                               reload_fiber_each_step=True) # reloading fiber in case of taper
        
        
        
        
        ################################################
        ###    DATA PROCESSING
        ################################################
        
        # convert frequency [THz] to wavlength [nm]
        wl_nm = np.array([1e-3*299792458./f for f in pulse.F_THz])
        
        # logical values in this list mark the wavelengths within the region of interest
        iis   = np.logical_and(wl_nm > low_nm, wl_nm < high_nm)
        
        # get zero-dispersion wavelengths along the fiber
        ZDW = []
        for i,y in enumerate(Y):
            x = y # x is global parameter; dirty but necessary
            D = fiber1.Beta2_to_D(pulse) # this is pretty slow :(
            ZDW.append(get_ZDW(wl_nm[iis][::-1],D[iis][::-1]))
        
        # transform the data a bit
        xW = wl_nm[iis]
        inSpec  = pulse.AW[iis]
        outSpec = pulse_out.AW[iis]
        zW = dB( power(np.transpose(AW)[:, iis]) )
        zT = dB( power(np.transpose(AT)) )
        
        # convert to dB/nm
        # check conversion factor
        #print( xW)
        delta_nm  = [(xW[i-1]-xW[i+1])/2. for i in range(1,len(xW)-1,1)]
        delta_nm  = np.array( [delta_nm[0]] + delta_nm + [delta_nm[-1]] )
        dB_nm_in  = dB( power(inSpec) / delta_nm )
        dB_nm_out = dB( power(outSpec) / delta_nm )
        
        
        #################################################
        ###                PLOTTING                   ###
        #################################################
        
        # Forgot to close the plot window from the previous run?
        # --> Lets' close it now.
        plt.close()
        
        # set up plots
        fig = plt.figure(figsize=(12,12))
        ax0 = plt.subplot2grid((3,2), (0, 0), rowspan=1)
        ax1 = plt.subplot2grid((3,2), (0, 1), rowspan=1)
        ax2 = plt.subplot2grid((3,2), (1, 0), rowspan=2, sharex=ax0)
        ax3 = plt.subplot2grid((3,2), (1, 1), rowspan=2, sharex=ax1)
        
        # plot the data
        Y_mm = Y * 1e3  # convert distance to mm
        ax0.plot(xW,  dB_nm_out, color = 'blue')
        ax0.plot(xW,  dB_nm_in,     color = 'grey')
        ax1.plot(pulse_out.T_ps,  dB(power(pulse_out.AT)), color = 'blue')
        ax1.plot(pulse.T_ps,      dB(power(pulse.AT)),     color = 'grey')
        ax2.pcolormesh(xW, Y_mm, zW, vmin=np.max(zW)-40.0, vmax=np.max(zW),cmap=cmap)
        ax3.pcolormesh(pulse.T_ps, Y_mm, zT, vmin=np.max(zT)-40.0, vmax=np.max(zT),cmap=cmap)
        
        # draw lines indicating taper region
        for z in z_stops:
            ax2.axhline(z*1e3, linestyle="--", color="white")
            ax3.axhline(z*1e3, linestyle="--", color="white")
        
        # plot the zero-dispersion wavelengths
        for i,y in enumerate(Y_mm): # for each propagation step...
            if len(ZDW) > 0:        # (only proceed if there are any ZDWs)
                for z in ZDW[i]:    # ...plot each zero-dispersion point
                    ax2.plot(z, y, color="white", marker=".")
        
        # format axes
        ax0.set_ylim(-40, 10)
        ax1.set_ylim(  0, 40)
        ax1.set_xlim(np.min(pulse.T_ps), np.max(pulse.T_ps))
        ax2.set_xlim(low_nm, high_nm)
        ax2.set_ylim(0,Y_mm[-1])
        ax3.set_ylim(0,Y_mm[-1])
        ax0.set_ylabel('power [dB]')
        ax2.set_ylabel('propagation distance [mm]')
        ax2.set_xlabel('wavelength [nm]')
        ax3.set_xlabel('time [ps]')
        plt.tight_layout()
        plt.suptitle('waist = %s\ntaper length = %s m\nwaist len %s'
                     %(Waist_scale,taper_length,Waist_len),x=0.25)
        #plt.show()
        
        
        #################################################
        ###           WRITE RESULT TO FILE            ###
        #################################################
        
        if saveTXT or savePlot:
            
            # create a folder for the results, if it doesn't exist already
            folder = workingDir + "/results"    # path to results folder
            if not os.path.exists(folder): os.makedirs(folder)
            
            # create file name for the file to be added (contains date and time)
            now = str(dt.datetime.now()).split(".")[0]  # split to remove ms value
            now = now.replace(":",".").replace(" ", "_")# reformat date/time string
        
        if saveTXT:
            # header for the file contains parameters of the simulation
            header = now + "\n"  # write time/date also into header
            for e in ["shape", "FWHM_ps", "carrier_nm", "frep_MHz", "pow_mW", "GDD", "propa_len", 
                      "untap_len", "taper_len", "waist_len", "uptap_len", "untap_scale", 
                      "waist_scale", "workingDir", "axis", "window_ps", "steps", "points"]:
                header += e + " = " + str(eval(e)) + "\n"       # value of each parameter
            header += "\nwavelength [nm]\t power spectral density [dB]" # column headings
            
            # write the spectrum at the fiber output into the file
            result = np.transpose( np.array([xW[::-1],  dB_nm_out[::-1]]) )
            np.savetxt(folder+"/"+now+".txt", result, fmt="%f", delimiter="\t", header=header)
            print("Spectrum saved under " + folder + "/" + now + ".txt")
        
        
        if savePlot:
            # expand plot to full screen and save the figure as png
            #plt.get_current_fig_manager().window.showMaximized()
            plt.savefig(folder + "/" +str(pow_mW)+'mW_'+ now + ".png", dpi = 570)
            print("Plot saved under " + folder + "/" + now + ".png")
        
        plt.show()
        
        
        
        
        print('----------------------------------------------------')
        print('Pulse duration: \t %s fs'%(FWHM_ps*1e3))
        print('Repetition rate: \t %s GHz'%(frep_MHz*1e-3))
        print('Coupled power: \t\t %s W'%(pow_mW*1.25))
        print('GDD: \t\t\t\t %s ps^2'%(GDD))
        
        print('----------------------------------------------------')
        print('PCF Length:\t\t\t\t\t %s cm'%(propa_len*100))
        print('Initial untapered length:\t %s cm'%(untap_len*100))
        print('Tapered length:\t\t\t\t %s cm'%(taper_len*100))
        print('Waist length:\t\t\t\t %s cm'%(waist_len*100))
        print('Uptaper length:\t\t\t\t %s cm'%(uptap_len*100))
        print('Diameter prior tapering:\t %s um'%(untap_scale*125))
        print('Waist diameter in tapering:\t %s um'%(waist_scale*125))
        
        
        end_iteration = np.where(np.abs(xW - 1610)<0.45)[0][0]
        begin_iteration = np.where(np.abs(xW - 1500)<0.45)[0][0]
        
        
        hh = end_iteration 
        while xW[hh]<=xW[0]:
            if np.abs(dB_nm_out[hh]-(-30))<0.5:
                minimo = dB_nm_out[hh]
                #print(xW[hh])
                break
            else:
                hh-=1
        
        ll = begin_iteration
        while xW[ll]>=xW[-1]:
            #print(xW[ll])
            if np.abs(dB_nm_out[ll]-(-30))<0.5:
                minimo2 = dB_nm_out[ll]
                break
            else:
                ll+=1
        
        
            
        low_lambda = np.where(dB_nm_out == minimo2)[0][0]
        high_lambda = np.where(dB_nm_out == minimo)[0][0]
        #print(xW[low_lambda],dB_nm_out[low_lambda])
        #print(xW[high_lambda],dB_nm_out[high_lambda])
        
        SB_Range = xW[high_lambda]-xW[low_lambda]
        #print(xW[low_lambda],xW[high_lambda],SB_Range)
        
        result_opt = [taper_len,Waist_scale,Waist_len,xW[low_lambda],xW[high_lambda],SB_Range]
        
        opt.append(result_opt)


Opt = pd.DataFrame(opt)
new_column_names = ['Taper_length', 'Waist_scale','Waist_len','High_wavelength','Low_wavelength', 'SC_range']
Opt.columns = new_column_names


plt.rc('xtick', labelsize=17) 
plt.rc('ytick', labelsize=17) 
plt.rcParams.update({'font.size': 17})
plt.rcParams['figure.figsize'] = (14, 10)


now = str(dt.datetime.now()).split(".")[0]  # split to remove ms value
now = now.replace(":",".").replace(" ", "_")# reformat date/time string



if len(Waist_Scale) >= 3:
    
    k=0
    p=0
    plt.figure()
    for h in Taper_Length:
        
        Label = 'Taper length = %s m'%(np.round(h,3))
        plt.plot(Opt.Waist_scale[k+0:k+10],Opt.SC_range[k+0:k+10],'-o',label=Label,alpha=(1-p*0.1))
        k = k + 10
        p = 0*p+1
    
    plt.legend()
    plt.grid()
    plt.title('Uptaper = %s , waist length = %s, untap_length = %s'%(uptap_len,waist_len, untap_len))
    plt.xlabel('Waist scale')
    plt.ylabel('Broadening range [nm]')
    plt.savefig('Taper length, waist scale optimization'+ now + ".png", dpi = 570)
    
    
    
    k=0
    p=0
    plt.figure()
    for h in Taper_Length:
        
        Label = 'Taper length = %s m'%(np.round(h,3))
        plt.plot(Opt.Waist_scale[k+0:k+10],Opt.High_wavelength[k+0:k+10],'-o',label=Label,alpha=(1-p*0.1))
        k = k + 10
        p = 0*p+1
    
    plt.legend()
    plt.grid()
    plt.title('Uptaper = %s , waist length = %s, untap_length = %s'%(uptap_len,waist_len, untap_len))
    plt.xlabel('Waist scale')
    plt.ylabel('Higher wavelength accessed [nm]')
    plt.savefig('Taper length, waist scale opt_Higher wavelength accessed'+ now + ".png", dpi = 570)
    
    k=0
    p=0
    plt.figure()
    for h in Taper_Length:
        
        Label = 'Taper length = %s m'%(np.round(h,3))
        plt.plot(Opt.Waist_scale[k+0:k+10],Opt.Low_wavelength[k+0:k+10],'-o',label=Label,alpha=(1-p*0.1))
        k = k + 10
        p = 0*p+1
    
    plt.legend()
    plt.grid()
    plt.title('Uptaper = %s , waist length = %s, untap_length = %s'%(uptap_len,waist_len, untap_len))
    plt.xlabel('Waist scale')
    plt.ylabel('Lower wavelength accessed [nm]')
    plt.savefig('Taper length, waist scale opt_Lower wavelength accessed'+ now + ".png", dpi = 570)


else:
    g=len(Waist_Length)
    k=0
    p=0
    plt.figure()
    for h in Taper_Length:
        
        Label = 'Taper length = %s m'%(np.round(h,3))
        plt.plot(Opt.Waist_len[k+0:k+g],Opt.SC_range[k+0:k+g],'-o',label=Label,alpha=(1-p*0.1))
        k = k + g
        p = 0*p+1
    
    plt.legend()
    plt.grid()
    plt.title('Uptaper = %s , waist scale = %s, untap_length = %s'%(uptap_len,waist_scale, untap_len))
    plt.xlabel('Waist length [m]')
    plt.ylabel('Broadening range [nm]')
    plt.savefig('Taper length, waist scale optimization'+ now + ".png", dpi = 570)
    
    
    
    k=0
    p=0
    plt.figure()
    for h in Taper_Length:
        
        Label = 'Taper length = %s m'%(np.round(h,3))
        plt.plot(Opt.Waist_len[k+0:k+g],Opt.High_wavelength[k+0:k+g],'-o',label=Label,alpha=(1-p*0.1))
        k = k + g
        p = 0*p+1
    
    plt.legend()
    plt.grid()
    plt.title('Uptaper = %s , waist scale = %s, untap_length = %s'%(uptap_len,waist_scale, untap_len))
    plt.xlabel('Waist length [m]')
    plt.ylabel('Higher wavelength accessed [nm]')
    plt.savefig('Taper length, waist scale opt_Higher wavelength accessed'+ now + ".png", dpi = 570)
    
    k=0
    p=0
    plt.figure()
    
    for h in Taper_Length:
        
        Label = 'Taper length = %s m'%(np.round(h,3))
        plt.plot(Opt.Waist_len[k+0:k+g],Opt.Low_wavelength[k+0:k+g],'-o',label=Label,alpha=(1-p*0.1))
        k = k + g
        p = 0*p+1
    
    plt.legend()
    plt.grid()
    plt.title('Uptaper = %s , waist scale = %s, untap_length = %s'%(uptap_len,waist_scale, untap_len))
    plt.xlabel('Waist length [m]')
    plt.ylabel('Lower wavelength accessed [nm]')
    plt.savefig('Taper length, waist scale opt_Lower wavelength accessed'+ now + ".png", dpi = 570)
