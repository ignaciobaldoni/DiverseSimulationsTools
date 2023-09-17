import matplotlib.pyplot as plt
import numpy as np
import PyCORe_main_cp as pcm

plt.close('all')
print('KECOMO Chip: D63_1_12GHz_F13_C12')
Num_of_modes = 2**7
D2 = 19e6 # Hz #D2 = 17e6 (According to Junqiu's report)
D1 = 12.13e9            # (According to Junqiu's report)
mu = np.arange(-Num_of_modes/2,Num_of_modes/2)
Dint = 2*np.pi*(mu**2*D2/2)

dNu_ini = 0.5e8
dNu_end = 2e8
coupling_in     = 0.47
coupling_out    = 0.40
Through_coupling_efficiency = coupling_in * coupling_out
print(Through_coupling_efficiency)

Output_EDFA = 0.200 # Watts

nn = 2500
ramp_stop = 0.99
dOm = 2*np.pi*np.concatenate([np.linspace(dNu_ini,dNu_end, int(nn*ramp_stop)),dNu_end*np.ones(int(np.round((1-ramp_stop)*nn)))])
dOm_total = (dOm[-1]-dOm[0])*1e-6
print('Complete detuning: %s MHz' % (np.round(dOm_total,2)))

print(dOm_total*1e6/D1) if (dOm_total*1e6)/D1<(0.05*2*np.pi) else print('Approximation fails')


k_ex = 0.998**mu*1e6*2*np.pi * 8
kappa_0 = 14e6*(2*np.pi) # (According to Junqiu's report)
coupling = np.max(k_ex/(kappa_0+k_ex))
linewidth = (k_ex + kappa_0)/(2*np.pi)
print('Critical couple') if coupling == 0.5 else print('')
print('Undercoupled') if coupling < 0.5 else print('Overcoupled')
print('Linewidth = %s MHz'%np.round(np.mean(linewidth)*1e-6,2))
print('Finesse = %s' % np.round(D1/np.mean(linewidth)),2)

dz = dOm[1]-dOm[0]

rel_T_detuning = 10e-10
slow_time = len(dOm)*rel_T_detuning

PhysicalParameters = {'n0' : 1.9,
                      'n2' : 2.4e-19,### m^2/W
                      'FSR' : D1,#1000e9 ,
                      'w0' : 2*np.pi*194e12,
                      'width' : 1.5e-6,
                      'height' : 0.9e-6,#1.35e-6,
                      'kappa_0' : kappa_0,
                      'kappa_ex' : k_ex ,#0.998**mu*25e6*2*np.pi,
                      'Dint' : Dint}

simulation_parameters = {'slow_time' : slow_time,#2*1e-6,
                         'detuning_array' : dOm,
                         'Raman_term' : 0.,
                         'fraction_Raman': 0.,
                         'noise_level' : 1e-7,
                         'output' : 'map',
                         'absolute_tolerance' : 1e-8,
                         'relative_tolerance' : 1e-8,
                         'max_internal_steps' : 2000}

P0 = Output_EDFA*coupling_in    #0.059
Pump = np.zeros(len(mu),dtype='complex')
Pump[0] = np.sqrt(P0)

print(P0*1000,'mW')

single_ring = pcm.Resonator(PhysicalParameters)

map2d=single_ring.Propagate_SplitStep(simulation_parameters, Pump,dt=1e-3/2)

#%%
fig = plt.figure()
ax = fig.add_subplot()
ax.set_xlabel("Normalized detuning")
ax.set_ylabel("Intracavity power, a.u.")
ax.plot(dOm*2/single_ring.kappa,np.mean(np.abs(map2d)**2,axis=1))

#%%
pcm.Plot_MapNew(np.fft.fftshift(np.fft.fft(map2d,axis=1),axes=1),dOm*2/single_ring.kappa)