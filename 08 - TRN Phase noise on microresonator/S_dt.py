# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 08:43:31 2020

@author: ibaldoni
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from Lectura_de_datos import Lectura

plots = True

# We read the JCM Wave results that we are interested

w = np.linspace(1e1,1e8,int(1e5))

D1 = 3.776e-3
D2 = 2.30e-3 #2.30mm
#D3 = 4e-3
Rs = [D1/2,D2/2]#, D3/2]        # m^2 

TRN_effect_on_detuning = -55

Width = 800
l_carrier = 1542.14

input_file = "scale_1.50_"+str(Width)+"nm.txt"

A_eff_micron2 = Lectura(input_file,Width,l_carrier,"effective mode area")
A_eff = A_eff_micron2 * 1e-12 

plt.close()

##### Del fit salen estos resultados para 1542.14nm
#Equivalent_width = 6.01536059634e-13
#A_eff   = Equivalent_width
#dr      = 5.20752721748e-07



for R in Rs:
    
    if R > 1.5e-3:
        Label = '12 GHz'      
        Color = 'r'
    else:
        Label = '20 GHz'
        Color = 'b'
        
        
    Fundamental_TRN_effect = [-84 if R == D1/2 else -80] [0]

#    dr = [np.abs(-4.1450e-11-4.0889e-7) if Width==800 else np.abs(-7.0853e-10-3.8896e-7)][0]  # Obtained from the mode profile of JCM Wave
    
    #### Del fit salen estos resultados para 1542.14nm
#    Equivalent_width = 6.01536059634e-13
#    A_eff   = Equivalent_width
    dr      = 5.20752721748e-07
    
    V_mode = 2*np.pi*R*A_eff
    
    rho = 3.29e3 # kg m^(-3)
    k = 30  # W m^(-1) K^(-1)
    C = 800 # J (kg K)^(-1)
    kb = 1.3806488e-23 #J K^(−1)
    T_set = 23
    T = 273 + T_set
    
    f0 = 194.1e12  # Resonance frequency
    n0 = 1.996  # Refractive index 
    dn_dT = 2.45e-5 # Thermo-optic coefficient
    
    factor = (f0/n0)*dn_dT
    
    td = np.pi**(1/3)*rho*C*dr**2/(4**(1/3)*k)
    
#    coeff1 = 9*np.sqrt(3)*dr**2/(8*(2*np.pi)**(2/3)*k*np.sqrt(td))
    
    coeff = 9*np.sqrt(3*td)/(8*np.pi)
    
    def S_dt(w):
        return 1/V_mode * (kb*T**2)/(np.sqrt(w)*rho*C) * coeff * ((1+(w*td)**(3/4))**2)**(-1)
#        return kb*T**2/V_mode * coeff1 * (np.sqrt(w)*(1+(w*td)**(3/4))**2)**(-1)
    
    S_df_w = factor**2 * S_dt(w)    
    S_df_f = 2*np.pi*S_df_w
    
    S_dphi_w = S_df_w*(2*np.pi/w)**2    
    S_dphi_f = 2*np.pi*S_dphi_w
    
    
    SSB_fund = 10*np.log10(0.5*S_dphi_f) + Fundamental_TRN_effect   
#    SSB_fund2 = 10*np.log10(S_dphi)-3 + Fundamental_TRN_effect   # 10.log10(Sf(f)) – 3dB [dBc/Hz]
    
    SSB_detun = 10*np.log10(0.5*S_dphi_f) + TRN_effect_on_detuning    
#    SSB_detun2 = 10*np.log10(S_dphi)-3 + TRN_effect_on_detuning    

    def SSB_fund_int(y,interes):
        x = np.log10(y/(2*np.pi))
        interpolacion = interpolate.interp1d(x, SSB_fund)
        interes_ = np.log10(interes)
        result = interpolacion(interes_)
        return result
    
    print(SSB_fund_int(w,1e4))


#%%    PLOTS     
    Linestyle='--'
    if R > 1.5e-3:
        Label = '12 GHz'      
        Color = 'r'
    else:
        Label = '20 GHz'
        Color = 'b'

    if plots == True:             
        plt.figure(num=1,figsize=(14,10))    
        plt.grid()
        plt.plot(w/(2*np.pi),np.sqrt(S_df_w)**2,linewidth = 2,label = Label, color=Color)
        plt.yscale('log')
        plt.xscale('log')
        plt.yticks(fontsize=25, rotation=0)
        plt.xticks(fontsize=25, rotation=0)        
        plt.xlabel('Frequency (Hz)',fontsize=23)
        plt.ylabel("${S_{\delta f}}^{1/2}$ (Hz / Hz$^{1/2}$)",fontsize=23)
        plt.legend(fontsize=17)
        plt.xlim(1e1,1e8)
    plt.grid()

    plt.figure(num=2, figsize=(14,10))   
    plt.grid()
    
    
    plt.plot(w/(2*np.pi),SSB_fund,linewidth = 2,label = Label+str(' Fundamental TRN'),color = Color)
    plt.plot(w/(2*np.pi),SSB_detun,linewidth = 1.2,label = Label+str(' TRN on detuning'),color = Color,linestyle=Linestyle)
#    plt.plot(w,SSB_fund2,linewidth = 2,label = Label+str(' Fundamental TRN'),color = Color)
#    plt.plot(1e4,SSB_fund_int(w,1e4),'o')
#    plt.plot(w,SSB_detun2,linewidth = 2,label = Label+str(' TRN on detuning'),color = Color,linestyle=Linestyle)
#    plt.yscale('log')
#    plt.hlines(-150,0,w[-1],linewidth = 0.5,linestyle=Linestyle)
    plt.xscale('log')
    plt.yticks(fontsize=25, rotation=0)
    plt.xticks(fontsize=25, rotation=0)
    plt.ylim(-160,0)

    plt.xlabel('Frequency (Hz)',fontsize=23)
    plt.ylabel("${SSB_{\phi}}$ (dBc / Hz)",fontsize=23)
    plt.legend(fontsize=17)
  
    def S_dt2(w):
        print('The another method')
        p = 2
        dz = 5.4*dr
        return kb*T**2/np.sqrt(np.pi**3 * k * rho * C * w)*np.sqrt(1/(2*p+1))*1/(R*np.sqrt(dz**2-dr**2))*1/(1+(w*td)**(3/4))**2
plt.grid()