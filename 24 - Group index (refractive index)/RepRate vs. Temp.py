# -*- coding: utf-8 -*-
"""
Created on Wed Jul 29 15:03:46 2020

@author: ibaldoni

More simple calculations
"""
import numpy as np
import matplotlib.pyplot as plt

c = 299792458

Diameter    = 3.770e-3  # Diameter of the ring in m
alpha_SiN   = 3.6e-6    # Thermal expansion coefficient of Silicon Nitride

alpha_Si    = 2.6e-6    # Thermal expansion coefficient of Silicon (Waver) 
alpha_SiO   = 5.5e-7    # Thermal expansion coefficient of Silicon (Fused Silica)

n_SiN = 1.9920688512833209
dn_dT = 2.45e-5 

rep_rate = c/(np.pi*Diameter*n_SiN)

l_pump = 1542   
f_pump =  c/(l_pump*1e-9)
#k = 16027
k = 1#int(f_pump/rep_rate)
print("k mode = %s" %k)


def dfrep_dT_v1(n,L,dn_dT,k):
    L = np.pi*L 
    return -((k*c/L)*(1/n**2)*dn_dT)

def frep_dL(rep_rate, deltaT,n,k):
    dL = alpha_Si*deltaT*Diameter
    new_D = Diameter + dL
    new_freq = k*c/(new_D*np.pi*n)
    dfl = new_freq - k*rep_rate
    dT = deltaT
    dfl_dT = dfl/dT
    return dfl_dT

def dfrepdT(Diameter,n,k, ng_int=1, Elastic = True, ThermoOptic = True):

        
    L = np.pi*Diameter
    dL_dT = alpha_Si*L
    
    
    dfrep_dT = -(k*c)/((L*n)**2)*(n*dL_dT + L *dn_dT)

    df_el = -(k*c)/((L*n)**2)*(n*dL_dT)

    df_to = -(k*c)/((L*n)**2)*(L *dn_dT)
    
    df_carrier = -(k*c)/((L*ng_int)**2)*(ng_int*dL_dT + L *dn_dT)
    df_el_carrier = -(k*c)/((L*ng_int)**2)*(ng_int*dL_dT)

    df_to_carrier = -(k*c)/((L*ng_int)**2)*(L *dn_dT)
    
    
    return dfrep_dT, df_el, df_to, df_el_carrier, df_to_carrier, df_carrier

T0 = 300
Temperatures = [301]

deltaT = [(i-T0) for i in Temperatures]
deltaT = np.asarray(deltaT)

df_dT = dfrep_dT_v1(n_SiN,Diameter,dn_dT,k)
print("Thermorefractive contribution\t = %s MHz/K"
      %round(df_dT*1e-6,4))

dfl_dT = frep_dL(rep_rate, deltaT,n_SiN,k)
print("Thermal expansion contribution\t = %s MHz/K"
      %round(dfl_dT[0]*1e-6,4))

dfrep_dT, _, _, _,_,_ = dfrepdT(Diameter,n_SiN,k,1)
print('Both contributions dfrep_dT\t\t = %s MHz/K'%round(dfrep_dT*1e-6,4))

print('Refractive index used: %s' % n_SiN)

def group_refractive_index(n,wl_nm, interes):
    
    dn      = np.diff(n)
    dlambda = np.diff(wl_nm)
    dn_dl = dn/dlambda
    
    wl_nm1 = wl_nm[:-1]
    n_l = n[:-1]
    
    n_g = n_l - wl_nm1 * dn_dl
    
    from scipy import interpolate

    n_gc = interpolate.interp1d(wl_nm1, n_g)
    n_interes = n_gc(interes)
    
    return n_g, n_interes




    
#%%   
    
    
Width = 800
label = "800 nm height, 1.5 um width (slow)"
input_file = "scale_1.50_"+str(Width)+"nm.txt"
data = np.loadtxt(input_file)
shift = 3 if "slow" in label else 0
wl   = [d[0] for d in data]
n    = [d[1+shift] for d in data]
t    = [d[2+shift] for d in data]
a    = [d[3+shift] for d in data]
a    = [x*1.0e12 for x in a] # convert o µm^2
# plotting
wl_nm  = [w*1.0e9 for w in wl]

n_g, n_gint = group_refractive_index(n,wl_nm, 1550)

dfrep2_dT, df_el, df_to, df_el_carrier, df_to_carrier, dfrep_dT_carrier = dfrepdT(Diameter,n_g,k, n_gint)

dfrep2_dT = dfrep2_dT*1e-6 #MHz/K
df_el = df_el*1e-6 
df_to = df_to*1e-6 

df_el_carrier = df_el_carrier*1e-6 
df_to_carrier = df_to_carrier*1e-6 
dfrep_dT_carrier = dfrep_dT_carrier*1e-6 

print('---------------------------------------------------------')

print("Thermorefractive contribution\t = %s MHz/K"
      %round(df_to_carrier,4))

print("Thermal expansion contribution\t = %s MHz/K"
      %round(df_el_carrier,4))

print('Both contributions dfrep_dT\t\t = %s MHz/K'%round(dfrep_dT_carrier,4))

print('And in the carrier = %s MHz/K'%round(dfrep_dT_carrier,4))



plt.figure()
plt.plot(wl_nm[0:-1],n_g, linewidth = 2)
plt.xlabel('Wavelength (nm)', fontsize = 17)
plt.ylabel('Group Index', fontsize = 17)
plt.yticks(fontsize=17, rotation=0)
plt.xticks(fontsize=17, rotation=0)
plt.grid()

fig, ax1 = plt.subplots()

ax2 = ax1.twinx()
lns1 = ax2.plot(wl_nm[0:-1],df_el, linewidth = 2, label = 'Thermo elastic contribution', linestyle='--')
lns2 = ax2.plot(wl_nm[0:-1],df_to, linewidth = 2, label = 'Thermo optical contribution', linestyle='--')
ax1.tick_params(direction='out', length=6, width=1, labelsize=17)
ax2.tick_params(direction='out', length=6, width=1, labelsize=17)
ax1.set_ylabel(r"Red line", fontsize=17)
ax2.set_ylabel(r"Dashed lines", fontsize=17)

#ax2.set_xlabel('Wavelength (nm)', fontsize = 17)
#ax2.ylabel("$df_{rep} / dT $", fontsize = 17)

lns3 = ax1.plot(wl_nm[0:-1],dfrep2_dT, linewidth = 2, label = 'Both contributions', color='r')
#plt.yticks(fontsize=17, rotation=0)
#plt.xticks(fontsize=17, rotation=0)
#ax1.legend(fontsize=17,loc='right')
#ax2.legend(fontsize=17,loc='center right')

# added these three lines
lns = lns1+lns2+lns3
labs = [l.get_label() for l in lns]
ax1.legend(lns, labs, fontsize=17, loc='center right')
plt.title("$df_{rep} / dT $ (MHz/K)", fontsize = 27)
ax1.set_xlabel('Wavelength (nm)', fontsize = 17)

plt.grid()

print(n_gint)




















