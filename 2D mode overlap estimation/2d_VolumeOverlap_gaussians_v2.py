# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 17:02:16 2022

@author: ignacio
"""

import numpy as np
import matplotlib.pyplot as plt

plt.rc('xtick', labelsize=17) 
plt.rc('ytick', labelsize=17) 
plt.rcParams.update({'font.size': 17})
plt.rcParams['figure.figsize'] = (14, 10)

from scipy import integrate
import time


plot3d = True
plot2d = True


def volume_Gaussian(amplitude, xy0, sigma_xy):
    return 2*np.pi*amplitude*sigma_xy[0]*sigma_xy[1]

def gaussian2D(x, y, amplitude, xy0, sigma_xy):
    x = x - xy0[0]
    y = y - xy0[1]
    return (amplitude*np.exp( -x**2/2/sigma_xy[0]**2 - y**2/2/sigma_xy[1]**2 ))**2                                                          


Theta = [0,0.00029,0.00058,0.00087,0.00116,0.00145,0.00174,0.00203,0.00232,0.00261,0.0029,0.00319,0.00348,0.00377,0.00406,0.00435]
Theta = [i*1e3 for i in Theta]

for run in range(0,1):
# First run
    print('First run')
    if run==0:
        w1 = 0.00024 #0.001
        w2 = 0.00024 #0.001
        Mu = [0,2.9E-05,5.8E-05,8.7E-05,0.000116,0.000145,0.000174,0.000203,0.000232,0.000261,0.00029,0.000319,0.000348,0.000377,0.000406,0.000435]
        Overlap=[]
    
    # Second run
    if run==1:
        w1 = 0.00019
        w2 = 0.00019
        Mu = [0,2.32E-05,4.64E-05,6.96E-05,9.28E-05,0.000116,0.0001392,0.0001624,0.0001856,0.0002088,0.000232,0.0002552,0.0002784,0.0003016,0.0003248,0.000348]
        Overlap=[]
    
    
    for i in Mu:
        print('position',i)
        mu = np.array([0,i])
        Sigma = np.array([ w1 , w1])
        
        mu2 = np.array([0, 0])
        Sigma2 = np.array([ w2 , w2 ])
        
        amp1 = 1/(2*np.pi*Sigma[0]*Sigma[1])#0.58134
        amp2 = 1/(2*np.pi*Sigma2[0]*Sigma2[1])#0.5814
        
        N = 1024
        X = np.linspace(np.min([mu,mu2])-3*np.max([Sigma,Sigma2]), np.max([mu,mu2])+3*np.max([Sigma,Sigma2]), N)
        Y = np.linspace(np.min([mu,mu2])-3*np.max([Sigma,Sigma2]), np.max([mu,mu2])+3*np.max([Sigma,Sigma2]), N)
        X, Y = np.meshgrid(X, Y)
        
        
        # Gaussians
        Z = gaussian2D(X, Y, amp1, mu, Sigma)
        Z2 = gaussian2D(X, Y, amp2, mu2, Sigma2)
        
        if plot3d == True:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax.plot_surface(X, Y, Z)
            ax.plot_surface(X, Y, Z2)
            # plt.savefig('3D_mu_'+str(i)+'.png')
               
        if plot2d == True:
            plt.figure()
            plt.contour(X, Y, Z)#, zdir='z')#, offset=-0.15, cmap=cm.viridis)
            plt.contour(X, Y, Z2)#, zdir='z')#, offset=-0.15, cmap=cm.viridis)
            # plt.savefig('mu_'+str(i)+'.png')
            
        
        
        u=mu2[0]     #x-position of the center
        v=mu2[1]      #y-position of the center
        a=Sigma2[0]*1.37     #radius on the x-axis
        b=Sigma2[1]*1.37   #radius on the y-axis
        
    
        f = lambda z, y, x: gaussian2D(x, y, amp1, mu, Sigma)
        Area1 = integrate.tplquad(f, -a+u, a+u, 
                                    lambda x: v - np.sqrt(b**2*(1 - ((x-u)/a)**2)),
                                  lambda x: v + np.sqrt(b**2*(1 - ((x-u)/a)**2)),
                                  lambda x, y: 0, lambda x, y: gaussian2D(x, y, amp2, mu2, Sigma2))
        
        
        
        AREA = Area1[0] 
        
        area = np.round(Area1[0],3)
        
        Overlap.append(AREA)
        time.sleep(0.01)

    
    

    
plt.figure(42)
plt.plot(Theta, Overlap/np.max(Overlap),'-bo')#,label='Spot size = %s µm' % (np.round(w1*1e6,2)))
plt.ylabel('Coupling')
plt.xlabel('Angular offset [mrad]')
plt.grid()
# plt.savefig('Coupling vs Angular Offset.png')