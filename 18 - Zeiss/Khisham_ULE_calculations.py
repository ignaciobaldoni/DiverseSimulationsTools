# -*- coding: utf-8 -*-
"""
Created on Thu Oct 12 15:34:07 2023
https://opg.optica.org/ao/abstract.cfm?uri=ao-62-13-3347
@author: ibaldoni
"""

##############################################################################


cte_glass_ppb = 0.9

cte_glass_ppm = cte_glass_ppb*1e-3
cte_glass = cte_glass_ppm*1e-6

dT = 10e-3
L = 10e-3

dL = cte_glass*dT*L
print(dL*1e12,'pm') 

##############################################################################

Bulk_Modulus = 34.1e9 #Pa
dP = 1e-5 
#Pa variation (10e-7 mBar)
dL_P = (dP/Bulk_Modulus)*L

print(dL_P*1e12,'pm')