# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 16:32:41 2022

@author: ibaldoni
"""

import numpy as np
import matplotlib.pyplot as plt

plt.rc('xtick', labelsize=17) 
plt.rc('ytick', labelsize=17) 
plt.rcParams.update({'font.size': 17})
plt.rcParams['figure.figsize'] = (14, 10)

import scipy.integrate as integrate
import time as dt


begin = dt.time()

Interf_AK = True
New_terms = False
AK_terms = True
Simulation_analysis = True

E0 = 1
omega = 194*2*np.pi # THz
TOD = 1

Duration = [0.25] #np.arange(1,3,0.50)
WIDTH = []

for duration in Duration:
    chirp = 0 #194*2
    
    n_points = 5000
    
    lim1 = -duration*6
    lim2 = -lim1
    
    delay = np.linspace(lim1,lim2,n_points)
    
    # For normalization
    factor8 = 1./(np.sqrt(2.*np.pi)*(duration/2.355))
    
    
    def gaussian(x, mu, sig):
        return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))
            #1./(np.sqrt(2.*np.pi)*sig)*np.exp(-np.power((x - mu)/sig, 2.)/2)
            
    def E_pulse(x): #2.355 x sigma = FWHM
        return gaussian(x, 0, duration/2.355)*np.exp(1j*(omega+chirp*x)*x)
    
    t = np.linspace(-duration*10,duration*10,n_points)
    plt.figure()
    plt.plot(t,E_pulse(t))
    plt.plot(t,gaussian(t, 0, duration/2.355))
    plt.grid()
    
    if Interf_AK == True:
        def IAC2(x,t0):
            return np.abs((E_pulse(x)+E_pulse(x-t0))**2)**2
        
        def IAC(t0):
            return integrate.quad(lambda x: IAC2(x,t0), -duration*100,duration*100) #-np.inf, np.inf)
           
        IAC_delay = []
        for i in delay:
            IAC_delay.append(IAC(i)[0])
            
        IAC_delay = np.array(IAC_delay)
    
        plt.figure()
        plt.plot(delay,IAC_delay*factor8#/np.max(IAC_delay)#*factor8
                  ,'g-',label = 'Interferometric Autocorrelation')
        plt.legend()
        plt.xlabel('Delay [ps]')
        plt.ylabel('Autocorrelation function')
        plt.grid()
    
    if AK_terms == True:
        def AC(x,t0):
            return np.abs(2*E_pulse(x)*E_pulse(x-t0))**2
        
        def IntensityAC(t0):
            return integrate.quad(lambda x: AC(x,t0), -duration*100,duration*100) #-np.inf, np.inf)
        
        # p = 0
        IntensityAC_delay = []#np.zeros(len(delay))
        for i in delay:
            IntensityAC_delay.append(IntensityAC(i)[0])
            
        IntensityAC_delay = np.array(IntensityAC_delay)
        

    
    if New_terms == True:
        def I2(x,t0):
            return np.abs(E_pulse(x)**2+E_pulse(x-t0)**2)**2
        
        def Interf(t0):
            return integrate.quad(lambda x: I2(x,t0), -duration*100,duration*100) #-np.inf, np.inf)
        
        
        Interf_delay = []
        for i in delay:
            Interf_delay.append(Interf(i)[0])
            
            
        Interf_delay = np.array(Interf_delay)   

    
    
    from scipy import interpolate
    def fwhm(x, y, k=10):
        
        """
        Determine full-with-half-maximum of a peaked set of points, x and y.
        Assumes that there is only one peak present in the datasset.  The function
        uses a spline interpolation of order k.
        """
        class MultiplePeaks(Exception): pass
        class NoPeaksFound(Exception): pass
    
        half_max = (np.amax(y)+np.amin(y))*0.5
        s = interpolate.splrep(x, y - half_max)
        roots = interpolate.sproot(s)
        
        if len(roots) > 2:
            return abs(roots[1] - roots[0])
            raise MultiplePeaks("The dataset appears to have multiple peaks, and "
                    "thus the FWHM can't be determined.")
        elif len(roots) < 2:
            raise NoPeaksFound("No proper peaks were found in the data set; likely "
                    "the dataset is flat (e.g. all zeros).")
        else:
            return abs(roots[1] - roots[0])
        
    
    if Simulation_analysis == True:
        
        magic_number = 20 
        magic_number2 = int(n_points/2)
        
        T = 2*np.pi/omega
        
        import copy
        Norm_AK = IAC_delay/np.max(IAC_delay)
        Frequ_spec = (np.fft.fft(Norm_AK))
        
        Frequ_spec2=copy.copy(Frequ_spec);
        
        ## Filters high frequencies 
        Frequ_spec2[8:len(Frequ_spec)]=0
    
        IAC=np.fft.ifft(Frequ_spec2) #intensity autocorrelation
        
        width=fwhm(delay, np.abs(IAC))    
        theory_width = fwhm(delay, IntensityAC_delay)
              
        Frequ_spec[0]=0
        F=np.argmax(Frequ_spec[magic_number:magic_number2])
        
        F += magic_number
    
        period=len(Norm_AK)/F
        print('Number of points that are in one period is', period)
        taxis=T*np.linspace(0, len(Norm_AK), len(Norm_AK))/period
        
        ind=np.argmax(Norm_AK)     
        taxis=taxis-taxis[ind] #Centers plot at maximum
        
        width=fwhm(taxis, np.abs(IAC))    
        theory_width = fwhm(taxis, IntensityAC_delay)
        
        # Autocorrelation FWHM = 1.41 x FWHM of gaussian pulse  
        print('Expected duration', duration,'ps')
        print('Algorithm width', width/1.41,'ps') 
        print('Theoretical width', theory_width,'ps')
        
        WIDTH.append(width/1.41)
    
    end = dt.time()
    print('----------- Code running for',end-begin,'s -----------------------')

    
WIDTH = np.array(WIDTH)
