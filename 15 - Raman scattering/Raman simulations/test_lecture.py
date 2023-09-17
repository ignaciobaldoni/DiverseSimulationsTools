import matplotlib.pyplot as plt
import numpy as np
import PyCORe_main_cp as pcm

#plt.close('all')
Num_of_modes = 2**8
#
D2 = 17e6 # Hz #D2 = 19e6

mu = np.arange(-Num_of_modes/2,Num_of_modes/2)
Dint = 2*np.pi*(mu**2*D2/2)
#plt.plot(Dint)

dNu_ini = 0.1e8 #-.2e8
dNu_end = 3.5e8
nn = 2700#2500
ramp_stop = 0.99
dOm = 2*np.pi*np.concatenate([np.linspace(dNu_ini,dNu_end, int(nn*ramp_stop)),dNu_end*np.ones(int(np.round((1-ramp_stop)*nn)))])

dz = dOm[1]-dOm[0]

rel_T_detuning = 10e-10
slow_time = len(dOm)*rel_T_detuning


P0 = 0.049 ### 0.7W ## 1W
Pump = np.zeros(len(mu),dtype='complex')
Pump[0] = np.sqrt(P0)

c = 299792458.0
center_wavelength = 1550 # In nm
frequency = c / (center_wavelength*1e-9) #285e12#194e12

FSR = 50E9 #in Hz
FM2FM = -84 # Because frequency division


PhysicalParameters = {'n0' : 1.9,
                      'n2' : 2.4e-19,### m^2/W
                      'FSR' : FSR,#1000e9 ,
                      'w0' : 2*np.pi*frequency,
                      'width' : 1.5e-6,
                      'height' : 0.9e-6,#1.35e-6,
                      'kappa_0' : 14e6*2*np.pi,
                      'kappa_ex' :  0.998**mu*8e6,#0.998**mu*25e6*2*np.pi,
                      'Dint' : Dint}

simulation_parameters = {'slow_time' : slow_time,#2*1e-6,
                         'detuning_array' : dOm,
                         'Raman_term' : 20e-15, #20e-11
                         'fraction_Raman': 0.2,
                         'noise_level' : 1e-7,
                         'output' : 'map',
                         'absolute_tolerance' : 1e-8,
                         'relative_tolerance' : 1e-8,
                         'max_internal_steps' : 2000}


single_ring = pcm.Resonator(PhysicalParameters)

map2d=single_ring.Propagate_SplitStep(simulation_parameters, Pump,dt=1e-3/2)
#%%
fig = plt.figure()
ax = fig.add_subplot()
ax.set_xlabel("Normalized detuning")
ax.set_ylabel("Intracavity power, a.u.")
ax.plot(dOm*2/single_ring.kappa,np.mean(np.abs(map2d)**2,axis=1))
#plt.show()

pcm.Plot_MapNew(np.fft.fftshift(np.fft.fft(map2d,axis=1),axes=1),dOm*2/single_ring.kappa,frequency,FSR, D2, FM2FM)

#%%
#fig = plt.figure()
#ax = fig.add_subplot()
#ax.set_xlabel("Detuning number")
#ax.set_ylabel("Normalized detuning")
#ax.plot(dOm*2/single_ring.kappa)
#plt.show()
#
##%%
#fig = plt.figure()
#ax = fig.add_subplot()
#ax.set_xlabel("Mode number")
#ax.set_ylabel("Dint")
#ax.plot(mu,Dint,'o')
#plt.show()
#%%
#pcm.Plot_MapNew(np.fft.fftshift(np.fft.fft(map2d,axis=1),axes=1),dOm*2/single_ring.kappa,frequency)
#
#df_optical = 0.0033e12
#FM2FM = -84 # Because frequency division
#
#df_microwave = np.sqrt(10**(FM2FM/10))*df_optical
#print(str(df_microwave*1e-6)+' MHz')