# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 16:47:42 2020

@author: Ignacio
"""

import nazca as nd
import nazca.geometries as geom
import numpy as np
from scipy import spatial


c = 299792458
freq = 10 #GHz

wl = 1550*1e-9 #1550 nm

n_g = 2.5

chip_size = 5000 #in microns
taper_size = 200

diameter = c/(n_g*freq*1e9*np.pi)*1E6 # In um
print(diameter)

diameter = 3776


offset_waveguide    = 921
valor_cero          = chip_size*0.5 - offset_waveguide
left_border         = -chip_size*0.5
right_border        = -left_border
dist_num_text       = 10
dist_waveguide      = 82
gap                 = 0.6
width_waveguide     = 1.5
length_tec          = 300







chip = nd.Polygon(layer=36, points=geom.square(length
                                              =chip_size))
resonator = nd.Polygon(layer=17, points=geom.ring(radius=diameter*0.5, N=4000, width=1.5))



waveguide = nd.Polygon(layer=10, points=geom.rectangle(length=chip_size-2*taper_size, height = 1.5))

tapers1 = nd.Polygon(layer=16, points=geom.taper(length=taper_size, width1=0.25, width2=1.5, grow=0))
tapers2 = nd.Polygon(layer=16, points=geom.taper(length=taper_size, width1=1.5, width2=0.25, grow=0))

tec = nd.Polygon(layer=4, points=geom.square(length=length_tec))

chip.put(-2500,-2500)
tec.put(2500-length_tec,2500-length_tec)






coord_1waveguide= -valor_cero+2*width_waveguide
coord = [(0.,coord_1waveguide)]

# I create a new variable to make easier to find the minimal distance between the waveguide and the resonator
new_variable = np.round(resonator.points,10) 
tree = spatial.KDTree(new_variable)
min_point = tree.query(coord)

initial_pos = min_point[0][0]
print(initial_pos)

resonator.put(0,initial_pos+gap)

  

for i in range(0,9):
    
    position_y = -(valor_cero + i*dist_waveguide)
    
    waveguide.put(left_border+taper_size,position_y)
    tapers1.put(left_border, position_y)
    tapers2.put(right_border-taper_size, position_y)
    if i != 0:
        nd.text(text = str(float(i)), height = 27, layer = 4).put(left_border+taper_size*0.25, position_y+dist_num_text)
        nd.text(text = str(float(i)), height = 27, layer = 4).put(right_border-taper_size*0.75, position_y+dist_num_text)
    
f = nd.Font('cousine')   
message = "This chip layout was made by Ignacio." 
nd.text(text=message, height=70, layer=4).put(-2300, 2300)  

nd.export_gds(filename='Igna_1st_chip.gds')