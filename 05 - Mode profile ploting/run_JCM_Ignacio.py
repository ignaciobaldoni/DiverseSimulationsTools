
"""
If "requests" package is missing, open a command window and type:
easy_install pip
pip install requests
"""

##############################################
###              USER AREA                 ###
##############################################

# wavelength sweep, modes
center  = 1550.0e-9     # center of the wavelength scanning range. Scan will start from here
stop1   = 2000.0e-9     # Long-wavelength stop of the wavelength scan
stop2   =  550.0e-9     # Short-wavelengh stop of the wavelength scan #550.0e-9
delta   =   10.0e-9     # Step size of the wavelength scan
n_guess =     1.800     # average index of the desired modes at "center"
height  =     0.800     # waveguide height [um]
start   =     1.500     # Start value of width scan [um]
stop    =     1.490     # Stop value of width scan [um]
step    =     0.050     # Step size of width scan [um]



#n_Si3N4 = 1.4496309898590634
#
#n_SiO2 = 1.732509249636967

##############################################
###             PROGRAM CODE               ###
##############################################

import sys
import os
sys.path.append(os.path.join(os.getenv('JCMROOT'), 'ThirdPartySupport', 'Python'))
import jcmwave
jcmwave.info()
import numpy as np
import datetime
from scipy.interpolate import InterpolatedUnivariateSpline as spline

import matplotlib.pyplot as plt
from Mode_profile_reading import Mode_profile

###########################################
###          aux functions              ###
###########################################

def n_SiO2(wl,T,T0):
    # Sellmeier equation and coefficients for fused silica
    c = [0.6961663, 0.0684043**2, 
         0.4079426, 0.1162414**2, 
         0.8974794, 9.896161**2]
    wl = wl*1.0e6 # conversion into microns
    l2 = wl * wl
    p1 = (c[0]*l2)/(l2-c[1])
    p2 = (c[2]*l2)/(l2-c[3])
    p3 = (c[4]*l2)/(l2-c[5])

    f_T = 0.95e-5 # RIU/°C   
#    T0 = 293
    
    dn_SiO2_dT = 1.09e-5 #-1.167e-7 + 1.727e-8 * T + 1.861e-10 * T**2 \
                    #- 5.781e-12 * T**3 + 4.221e-16 * T**4
                    
    f_T = 0*dn_SiO2_dT 
    
    return np.sqrt(1.+p1+p2+p3) + f_T*(T-T0)


def n_Si3N4(wl,T,T0):
    # Sellmeier equation and coefficients (fit is reasonably good down to 500 nm)
    c = [3.40789621e+00 ,1.72579862e-02, -4.62288837e-01, 
         1.72570350e-02, 9.26379962e+00, -1.84302988e+04]
    wl = wl*1.0e6 # conversion into microns
    l2 = wl * wl
    p1 = (c[0]*l2)/(l2-c[1])
    p2 = (c[2]*l2)/(l2-c[3])
    p3 = (c[4]*l2)/(l2-c[5])
    
#    f_T = 2.45e-5 # RIU/°C    
#    T0 = 293
    
    dn_Si3N4_dT = 2.45e-5 # 3.211e-7 - 1.990e-8 * T + 8.856e-10 * T**2 \
                    #-1.375e-12 * T **3 - 1.105e-15 * T **4   
                    
    f_T = 0*dn_Si3N4_dT
    
#    np.sqrt(1.+p1+p2+p3) + f_T*(T-T0)

    return np.sqrt(1.+p1+p2+p3) +  f_T*(T-T0) 


def run_jcm(T_set, T_anterior, lambda_0, n_guess, width, height , modes=2):
    # set parameters
    keys = {'lambda0': lambda_0, 'permittivity_SiO2': n_SiO2(lambda_0, T_set,T_anterior)**2, 
            'permittivity_Si3N4': n_Si3N4(lambda_0,T_set,T_anterior)**2,
            'n_guess': n_guess, 'height': height, 'width': width, 'modes': modes}
    # run JCM mode solver
    results = jcmwave.solve('./project.jcmpt', keys=keys)
    print "guess was:", n_guess
    print "width", width, ": ", lambda_0*1e9, "nm\t", "completed"
    
    plt.figure()
    _,_, ax = Mode_profile(plot_mode = True,
                                 fit_mode = False,
                                 plot_eqA = False,
                                 other_plots = False)
                                 
    plt.show()
    plt.pause(0.5)
    return process_results(results, lambda_0)


def process_results(results, lambda_0):

    
    """
    Function to restructure the results of the mode solver.
    Computes the effective mode area and propagation loss.
    "results": results array as returned by the JCM mode solver
    "lambda_0": Vacuum wavelength for which the mode solver was operated.
    """
    R = [lambda_0]
    for i in range(len(results[1]['Integrate(E^2)']))[::-1]:
        n  = results[0]['eigenvalues']['effective_refractive_index'][i]
        I1 = results[1]['Integrate(E^2)'][i][1].real
        I2 = results[1]['Integrate(E^2)'][i][0].real
        I3 = results[2]['Integrate(E^4)'][i][1].real
        I4 = results[2]['Integrate(E^4)'][i][0].real
        a = ((I1+I2)**2)/(I3+I4)  # compute effective mode area
        g_ = np.exp(-4*np.pi*n.imag/lambda_0) # factor 4 instead of 2: intensity, not amplitude!
        g  = 10.*np.log10(g_) # convert to damping coefficient in dB/m
                  
        


        R += [n.real, g, a]
        
        
        
       
    return R


def save_results(width, results, T_set, folder="."):
    """
    Writes the simulation results for a given fiber scaling into a file.
    File name is automatic, contains scale.
    "folder": Working directory. Files are saved in subfolder "dispersion".
    "scale": Scale factor of the fiber structure (1.0 is original size).
    "results": Matrix that contains (processed) simulation results.
    """
    # Check if subfolder "dispersion" exists (create if necessary).
    if not os.path.exists(folder+"/dispersion"+str(T_set)): 
        os.makedirs(folder+"/dispersion"+str(T_set))
    # create column headers for output table. First column: wavenength
    header = "wavelength [m]\t"
    for i in range(1,len(results[0])-1,3):
        header += "effective refractive index\ttransmission [dB/m]\teffective mode area [m]\t"
    # automatically generate file name (+ rest of the path), save file
    path = folder + "/dispersion"+str(T_set)+ "/scale_" + "{:.2f}".format(width) + ".txt"
    np.savetxt(path, results, delimiter="\t", header=header)
    print "saved to:", path


def N_center(R):
    """
    Returns the center effective index for a range containing several modes.
    Useable as "n_guess" to find these modes with the mode solver.
    "R": results matrix.
    """
    Nc = []
    for r in R:
        N = [r[i] for i in range(1,len(r),3)]
        Nc.append( (min(N)+max(N))/2. )
    return Nc


def wavelengths(R):
    """ extracts the list of wavelengths from results array R"""
    return [r[0] for r in R]


###################################################
####            PROGRAM EXECUTION               ###
###################################################


#run_jcm(lambda_0, n_guess, width, height)


if __name__ == "__main__":
    
    # Record program start time. Will be used to monitor runtime.
    t1 = datetime.datetime.now()
    
    
    Temperatures = [0]#np.linspace(280,320,10)
    
    
    T0 = Temperatures[0]
    for T in Temperatures:
        
       
        
        print('Studied temperature %s °C' % (T-273))
        print("Temperatura anterior %s °C" % (T0-273))
        
        # Make a few small initial steps along the taper, Simulate using center wavelength.
        # These small taper steps will be used to start the spline extrapolation of the 
        # effective refractive index to get n_guess for the mode solver (at center wavelength)
        steps = [start+0.01, start+0.02]  # list of steps along the taper
        # "R2": traces evolution of the madal parameters along the taper (at center wavelength)
        R2  = [ run_jcm(T, T0, center, n_guess, steps[0], height) ]
        R2 += [ run_jcm(T, T0, center, n_guess, steps[1], height) ]
        
        # List index of center wavelength in results list.
        center_ind = int(round((center-stop2)/delta,0))
    
    
        # Start scanning along the taper (scale parameter to describe fiber down-scaling)
        for width in np.arange(start, stop-1e-5, -step):
            
            # Compute n_guess (at center wl) for new taper scale by spline extrapolation
            n_guess = spline(steps, N_center(R2), k=min(2,len(R2)-1))(width)
            
            # Run the mode solver twice at center wl and slightly detuned from center wl
            # "R": results matrix for the current fiber scale
            R  = [ run_jcm(T, T0, center, n_guess, width, height) ]
            R += [ run_jcm(T, T0, center+0.5e-9, n_guess, width, height) ]
            
            # Wavelength scan twoards long-wavelength stop
            for wl in np.arange(center+delta, stop1+1e-10, delta):
                # n_guess for each wavlength computed through spline extrapolation from old values
                n_guess = spline(wavelengths(R), N_center(R), k=min(2,len(R)-1))(wl)
                # run the mode solver and append result to result matrix "R"
                R.append( run_jcm(T,T0,wl, n_guess, width, height) )
            del R[1] # remove the extra data point that we needed for initial extrapolation
            
            # Wavelength scan towards short-wavelength stop
            for wl in np.arange(center-delta, stop2-1e-10, -delta):
                # n_guess for each wavlength computed through spline extrapolation from old values
                n_guess = spline(wavelengths(R), N_center(R), k=min(2,len(R)-1))(wl)
                # run the mode solver and append result to result matrix "R"
                R.insert(0, run_jcm(T,T0,wl, n_guess, width, height))
            
            # Save the result of the wavelength scan to a file
            save_results(width, R, T)
            # Update lists that track evolution along taper (for next extrapolation iteration)
            steps, R2 = [width] + steps, [R[center_ind]] + R2
            
#            T0 = T
        
        # record stop time of program, print start, stop and runtime to screen
    t2 = datetime.datetime.now()
    print t1, t2, t2-t1