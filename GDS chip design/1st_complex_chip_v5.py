# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 09:49:54 2020

@author: Ignacio

My first GDS desing of a more complex chip 

All values are in um
"""

import gdspy
from scipy import interpolate
       
if __name__ == "__main__":

    lib = gdspy.GdsLibrary()
    
    c = 299792458.
    chip_size = 5000.0

# Negative resist example
    width = 1.5     
    
    ng = 4 # Simulations results
    
    Rep_rate = 100
    Rep_rate_hz = Rep_rate*1e9
    
    ring_radius = c/(2*Rep_rate_hz*ng)*1e6 # In um
    ring_radius = 500.0
    print('Radius of the microresonator = %s um'% (ring_radius))
    print('Repetition rate = %s GHz' % Rep_rate)
    
    bend_radius = ring_radius
    taper_len = 200.0

    width_taper = 0.25    
    dist_waveguides = 50    
    offset_from_zero = 0.5*ring_radius
    
    gap = 0.165
    N_resonators = 2
            
    # Resonators in each row
    n_resonators = [i for i in range(N_resonators)]    
    n_filas = 2 # Number of rows in the chip  
    
    dist_entre_filas = ring_radius*4*0.9
    
    
#    Mis variables falopa
    corrimiento_en_x = 1.5*ring_radius
    posy_taper = 2.5*ring_radius    
    distancia_curv = dist_entre_filas #4*ring_radius
    Offset_ring = 0.5*(distancia_curv + corrimiento_en_x)

    centro1 = 0   # x offset of the resonators
    centro2 = ring_radius + gap + 0.5*width   # y offset of the resonators     
    
##  Total surface of the chip        
    chip = gdspy.Cell('Rectangle')
    chip.add(gdspy.Rectangle((-chip_size*0.5,-chip_size*0.5),
                             (chip_size*0.5,chip_size*0.5),
            layer = 2))
    
##  Microresonator    
    ring = lib.new_cell("NRing")
    ring.add(
        gdspy.Round((centro1, centro2), 
                    ring_radius, 
                    ring_radius - width, 
                    tolerance=0.001, 
                    layer = 0)
    )

##  Taper of the left side of the chip
    taper_left = lib.new_cell("NTaper_L")
    taper_left.add(gdspy.Path(width_taper, (0, offset_from_zero)).segment(taper_len, "+x", final_width=width, layer = 5))

##  Taper of the right side of the chip    
    taper_right = lib.new_cell("NTaper_R")
    taper_right.add(gdspy.Path(width_taper, (0, posy_taper + offset_from_zero
                                             )).segment(taper_len, "-x", final_width=width, layer = 5))
    
    c = lib.new_cell("Negative")
    
    c.add(chip)
    
    # Horizontal text with height 2.25
    htext = gdspy.Text("Company, IB, %s GHz, Radius = %s, Gap = %s,  Width = %s " % (Rep_rate, ring_radius, gap, width), 
                       70, (-chip_size*.5 + taper_len, chip_size*.5 - taper_len), layer = 4)
    
    c.add(htext)
     
    series = [0, -1, 1, -2, 2, -3, 3, -4, 4, -5, 5] # https://oeis.org/A001057
    
    for jj in range(n_filas):
        for i, gap in enumerate(n_resonators):                        

            posicion_y_en_el_chip = dist_entre_filas*series[jj]
                        
            dist_wg = -dist_waveguides*i
                        
            path = gdspy.FlexPath(#Position inicial del path!
                [(-chip_size*0.5+taper_len, 
                  dist_wg+posicion_y_en_el_chip+offset_from_zero)],
                width=width,
                corners="circular bend",
                bend_radius=bend_radius,
                gdsii_path=True,
                offset = 0,
                layer = 3)
                        
            # First coordinate defines the length of the channel
            path.segment(((i+1)* distancia_curv, 0), relative=True)
            
            # Second coordinate defines the length of the "going up" channel
            path.segment((0, posy_taper), relative=True)  
            
            ## Parte de arriba de la waveguide          
            path.segment((chip_size -2*taper_len - (i+1)* distancia_curv, 0), relative=True)
            
            c.add(path)
            c.add(gdspy.CellReference(ring, (i* distancia_curv + Offset_ring-
                                             chip_size*0.5, # X of the resonators
                                              dist_wg+posicion_y_en_el_chip+offset_from_zero)))
            
            text_l = gdspy.Text("N %s, Gap = %s" % (i, gap), 
                       10, (-chip_size*.5 + 0.5*taper_len, dist_wg+posicion_y_en_el_chip+offset_from_zero+10), layer = 4)
    
            c.add(text_l)
            
            text_r= gdspy.Text("N %s, Gap = %s" % (i, gap), 
                       10, (chip_size*.5 - taper_len, posicion_y_en_el_chip+\
                            posy_taper + offset_from_zero + 10 + dist_wg), layer = 4)
        
            c.add(text_r)

        
    ### Taper 1    
        c.add(gdspy.CellArray(taper_left, 1, len(n_resonators ),\
                              (0, -dist_waveguides),\
                              (-chip_size*0.5, posicion_y_en_el_chip)))
        
        
        
    ### Taper 2
        c.add(gdspy.CellArray(taper_right, 1, len(n_resonators ), \
                              (0, -dist_waveguides), \
                              (chip_size*0.5, posicion_y_en_el_chip)))
        
        
        
    
    # Save to a gds file and check out the output
    lib.write_gds("Double_ring_normal.gds")