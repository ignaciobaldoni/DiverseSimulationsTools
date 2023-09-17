# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 17:02:16 2022

@author: ibaldoni
"""

import numpy as np
from scipy.integrate import dblquad
import matplotlib.pyplot as plt

plt.rc('xtick', labelsize=17) 
plt.rc('ytick', labelsize=17) 
plt.rcParams.update({'font.size': 17})
plt.rcParams['figure.figsize'] = (14, 10)

from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from math import pi
from scipy import integrate
import time


plot3d = False
plot2d = False


def volume_Gaussian(amplitude, xy0, sigma_xy):
    return 2*np.pi*amplitude*sigma_xy[0]*sigma_xy[1]

def gaussian2D(x, y, amplitude, xy0, sigma_xy):
    x = x - xy0[0]
    y = y - xy0[1]
    return (amplitude*np.exp( -x**2/2/sigma_xy[0]**2 - y**2/2/sigma_xy[1]**2 ))**2                                                          

# Mu=np.linspace(0,0.001,1)
# Mu=np.linspace(0,0.01,50)
Theta = [0,0.00029,0.00058,0.00087,0.00116,0.00145,0.00174,0.00203,0.00232,0.00261,0.0029,0.00319,0.00348,0.00377,0.00406,0.00435]
Theta = [i*1e3 for i in Theta]

for run in range(0,2):
# First run
    if run==0:
        w1 = 0.00024776067021345/2 #0.001
        w2 = 0.00024776067021345/2##0.001
        Mu = [0,2.9E-05,5.8E-05,8.7E-05,0.000116,0.000145,0.000174,0.000203,0.000232,0.000261,0.00029,0.000319,0.000348,0.000377,0.000406,0.000435]
        Overlap=[]
    
    # Second run
    if run==1:
        w1 = 0.000198208536567177/2
        w2 = 0.000198208536567177/2
        Mu = [0,2.32E-05,4.64E-05,6.96E-05,9.28E-05,0.000116,0.0001392,0.0001624,0.0001856,0.0002088,0.000232,0.0002552,0.0002784,0.0003016,0.0003248,0.000348]
        Overlap=[]
    
    
    for i in Mu:
        print('position',i)
        mu = np.array([0,i])
        Sigma = np.array([ w1 , w1])
        
        mu2 = np.array([0, 0])
        Sigma2 = np.array([ w2 , w2 ])
        
        amp1 = 1/(2*np.pi*Sigma[0]*Sigma[1])*0.57801#0.58134
        amp2 = 1/(2*np.pi*Sigma2[0]*Sigma2[1])*0.57801#0.5814
        
        N = 1024
        X = np.linspace(np.min([mu,mu2])-3*np.max([Sigma,Sigma2]), np.max([mu,mu2])+3*np.max([Sigma,Sigma2]), N)
        Y = np.linspace(np.min([mu,mu2])-3*np.max([Sigma,Sigma2]), np.max([mu,mu2])+3*np.max([Sigma,Sigma2]), N)
        X, Y = np.meshgrid(X, Y)
        
        
        # Gaussians
        Z = gaussian2D(X, Y, amp1, mu, Sigma)
        Z2 = gaussian2D(X, Y, amp2, mu2, Sigma2)
        
        if plot3d == True:
            fig = plt.figure()
            ax = fig.gca(projection='3d')
            ax.plot_surface(X, Y, Z)
            ax.plot_surface(X, Y, Z2)
        
        # Two Gaussian
        # def A_minus_B(x, y, argsA, argsB):
        #     return gaussian2D(x, y, *argsA) - gaussian2D(x, y, *argsB)
        
        argsA = (amp1, mu, Sigma)
        argsB = (amp2, mu2, Sigma2)
        
        
        if plot2d == True:
            plt.figure()
            plt.contour(X, Y, Z)#, zdir='z')#, offset=-0.15, cmap=cm.viridis)
            plt.contour(X, Y, Z2)#, zdir='z')#, offset=-0.15, cmap=cm.viridis)
            
        
        
        u=mu2[0]     #x-position of the center
        v=mu2[1]      #y-position of the center
        a=2*Sigma2[0]     #radius on the x-axis
        b=2*Sigma2[1]   #radius on the y-axis
        
        
        if plot2d == True:
            t = np.linspace(0, 2*pi, 100)
            plt.plot( u+a*np.cos(t) , v+b*np.sin(t),'b')
            plt.grid(color='lightgray',linestyle='--')
        
        
    
        f = lambda y, x: gaussian2D(x, y, amp1, mu, Sigma)
        Area1 = integrate.dblquad(f, -a+u, a+u, 
                                    lambda x: v - np.sqrt(b**2*(1 - ((x-u)/a)**2)),
                                  lambda x: v + np.sqrt(b**2*(1 - ((x-u)/a)**2)))
        
        print('gaussian1',Area1)
        
        u=mu[0]     #x-position of the center
        v=mu[1]      #y-position of the center
        a=2*Sigma[0]     #radius on the x-axis
        b=2*Sigma[1]      #radius on the y-axis
        
        
        
        
        if plot2d == True:
            t = np.linspace(0, 2*pi, 100)
            plt.plot( u+a*np.cos(t) , v+b*np.sin(t),'r' )
            plt.grid(color='lightgray',linestyle='--')
            plt.ylabel('Y Axis')
            plt.xlabel('X Axis')
        
        
        
        f = lambda y, x: gaussian2D(x, y, amp2, mu2, Sigma2)
        Area2 = integrate.dblquad(f, -a+u, a+u, 
                                    lambda x: v - np.sqrt(b**2*(1 - ((x-u)/a)**2)),
                                  lambda x: v + np.sqrt(b**2*(1 - ((x-u)/a)**2)))
        
        print('gaussian2',Area2)
        
        AREA = Area1[0] + Area2[0]
        
        area = np.round(Area1[0] + Area2[0],3)
        
        print('Overlap:',area)
        
        Overlap.append(AREA)
        #time.sleep(0.01)
        # plt.title('Overlap = %s'%area)
    
    
    # error = 0.01
    # plt.title(r'Overlap = %s $\pm$ %s'%(area,error))
    
plt.figure(1)
plt.plot(Theta, Overlap/np.max(Overlap),'-bo')#,label='Spot size = %s µm' % (np.round(w1*1e6,2)))
plt.ylabel('Coupling')
plt.xlabel('Angular offset [mrad]')

# plt.legend()
# plt.title('Comparison between two radius')
plt.grid()
# plt.savefig('Same spot size for both modes.png')