# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 15:44:52 2022

@author: ibaldoni
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

def solve(m1,m2,std1,std2):
  a = 1/(2*std1**2) - 1/(2*std2**2)
  b = m2/(std2**2) - m1/(std1**2)
  c = m1**2 /(2*std1**2) - m2**2 / (2*std2**2) - np.log(std2/std1)
  return np.roots([a,b,c])

Overlap = []
def gaussian(x, mu, sig):
    return np.abs(1./(np.sqrt(2.*np.pi)*(sig))*np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.))))**2

w2 = 0.00024776067021345 
Mu = [0,2.9E-05,5.8E-05,8.7E-05,0.000116,0.000145,0.000174,0.000203,0.000232,0.000261,0.00029,0.000319,0.000348,0.000377,0.000406,0.000435]

for i in Mu:
    
    m1 = .0
    std1 = w2
    m2 = i
    std2 = w2
    
    if m1==m2 and std2==std1:
        m1=m1+1e-6
        # Otherwise the code crashes for not finding the intersection point
    
    #Get point of intersect
    result = solve(m1,m2,std1,std2)
    
    #Get point on surface
    x = np.linspace(-np.max(Mu)*2,np.max(Mu)*2,10000)
    
    
    plt.figure()
    plot1=plt.plot(x,gaussian(x,m1,std1)/np.max(gaussian(x,m1,std1)))
    plot2=plt.plot(x,gaussian(x,m2,std2)/np.max(gaussian(x,m2,std2)))
    
    plot3=plt.plot(result,gaussian(result,m1,std1)/np.max(gaussian(x,m1,std1)),'o')
    plt.plot(w2,gaussian(w2,m1,std1)/np.max(gaussian(x,m1,std1)),'o')
    # plot3=plt.plot(0,0.2,'wo',label = 'Overlap:%s'%area)
    
    #Plots integrated area
    r = result[0]
    olap = plt.fill_between(x[x>r], 0, gaussian(x[x>r],m1,std1)/np.max(gaussian(x,m1,std1)),alpha=0.3)
    olap = plt.fill_between(x[x<r], 0, gaussian(x[x<r],m2,std2)/np.max(gaussian(x,m2,std2)),alpha=0.3)
    
    plt.show()
    
    from scipy import integrate
    x1 = lambda x: gaussian(x,m1,std1)
    x2 = lambda x: gaussian(x,m2,std2)
    
    part1 = integrate.quad(x2, 0, r)[0]
    part2 = integrate.quad(x1, r, 2*np.max(Mu))[0]
    
    area = part1 + part2
    # plt.title('Overlap (normalized): %s '%area)
    
    print('Overlap:',area)

    Overlap.append(area)    
    
    
# plt.figure(1)
Theta = [0,0.00029,0.00058,0.00087,0.00116,0.00145,0.00174,0.00203,0.00232,0.00261,0.0029,0.00319,0.00348,0.00377,0.00406,0.00435]
Theta = [i*1e3 for i in Theta]
plt.figure()
plt.plot(Theta, Overlap/np.max(Overlap),'-o',label='Spot size = %s µm' % (np.round(w2*1e6,2)))
plt.ylabel('Coupling')
plt.xlabel('Angular offset [mrad]')
plt.grid()
# plt.savefig('Same spot size for both modes.png')