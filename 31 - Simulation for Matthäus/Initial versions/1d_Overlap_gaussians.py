# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 15:44:52 2022

@author: ibaldoni
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
# norm.cdf(1.96)

def solve(m1,m2,std1,std2):
  a = 1/(2*std1**2) - 1/(2*std2**2)
  b = m2/(std2**2) - m1/(std1**2)
  c = m1**2 /(2*std1**2) - m2**2 / (2*std2**2) - np.log(std2/std1)
  return np.roots([a,b,c])

m1 = 3.0
std1 = 1.0
m2 = 5.0
std2 = 1.0

#Get point of intersect
result = solve(m1,m2,std1,std2)
r = result[0]

if m1==m2 and std2==std1:
    m1=m1+0.01
    # Otherwise the code crashes for not finding the intersection point

from scipy import integrate
x1 = lambda x: norm.pdf(x,m1,std1)
x2 = lambda x: norm.pdf(x,m2,std2)

part1 = integrate.quad(x1, r, 10)[0]
part2 = integrate.quad(x2, 0, r)[0]

area = part1+part2
print('Overlap:',area)


#Get point on surface
x = np.linspace(-2.5,9,10000)
plot1=plt.plot(x,norm.pdf(x,m1,std1))
plot2=plt.plot(x,norm.pdf(x,m2,std2))

plot3=plt.plot(result,norm.pdf(result,m1,std1),'o')
# plot3=plt.plot(0,0.2,'wo',label = 'Overlap:%s'%area)
plt.title('Overlap: %s '%area)
# plt.legend()

#Plots integrated area
olap = plt.fill_between(x[x>r], 0, norm.pdf(x[x>r],m1,std1),alpha=0.3)
olap = plt.fill_between(x[x<r], 0, norm.pdf(x[x<r],m2,std2),alpha=0.3)
plt.show()


