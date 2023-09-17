# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 17:04:54 2020

@author: ibaldoni
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


wavelength = pd.read_fwf('wavelengths.txt')
n0 = pd.read_fwf("refractive_index.txt")

dn = n0.diff()
u = wavelength.values


df = pd.DataFrame(dn * u)

ng = n0 - df

plt.plot(wavelength, ng)
plt.plot(wavelength, n0)