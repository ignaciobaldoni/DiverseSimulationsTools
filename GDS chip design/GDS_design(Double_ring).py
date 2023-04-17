# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 09:49:54 2020

@author: Ignacio

GDS Design, rings with different radii and gaps

All values are in um
"""

import gdspy

def Ring(centro1,centro2,ring_radius,width,name_ring,i):   
    
    name_ring.add(gdspy.Round((centro1, centro2), 
                    ring_radius, 
                    ring_radius - width, 
                    tolerance=0.001, 
                    layer = 1))
    
def Offset(ring_radius):
#    offset=0.5*np.max([ring_radius1,ring_radius2])
    offset=0.5*ring_radius
    return offset

def Taper(taper_name, position):        
    taper_name.add(gdspy.Path(width_taper, (0, position + offset_from_zero
                                             )).segment(taper_len, "-x", final_width=width, layer = 5))
       
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
    ring_radius1 = 450.0
    ring_radius2 = 250 
    print('Radius of the microresonator = %s um'% (ring_radius))
    print('Repetition rate = %s GHz' % Rep_rate)
    
    taper_len = 200.0

    width_taper = 0.25    
    dist_waveguides = 50  
    
    offset_from_zero = Offset(ring_radius1)
    
    min_gap = 0.6
    
    n_filas = 3 # Number of rows in the chip  
    
    centro1 = 0   # x offset of the resonators
    
##  Total surface of the chip        
    chip = gdspy.Cell('Rectangle')
    chip.add(gdspy.Rectangle((-chip_size*0.5,-chip_size*0.5),
                             (chip_size*0.5,chip_size*0.5),
            layer = 2))
   

##  Taper of the left side of the chip
    taper_left = lib.new_cell("NTaper_L")
    taper_left.add(gdspy.Path(width_taper, (0, offset_from_zero)).segment(taper_len, "+x", final_width=width, layer = 5))

##  Taper of the right side of the chip    
    taper_right1 = lib.new_cell("NTaper_R1")
    taper_right2 = lib.new_cell("NTaper_R2")    
    
    c = lib.new_cell("Negative")
    
    c.add(chip)
    
    # Horizontal text with height 2.25
    htext = gdspy.Text("Company, IB, %s GHz, Radius = %s, Gap = %s,  Width = %s " % (Rep_rate, ring_radius1, min_gap, width), 
                       70, (-chip_size*.5 + taper_len, chip_size*.5 - taper_len), layer = 4)
    
    c.add(htext)
     
    series = [0, -1, 1, -2, 2, -3, 3, -4, 4, -5, 5]
    
    for jj in range(n_filas):
                
        if jj == 1:
            N_resonators = 2           
            # Resonators in the row
            n_resonators = [i for i in range(N_resonators)]  
        else:
            N_resonators = 4      
            # Resonators in the row
            n_resonators = [i for i in range(N_resonators)]  
            
        for i in range(N_resonators):  
            
            ##  Microresonator    
            ring1 = lib.new_cell("NRing_type1"+str(jj)+str(i))
            ring2 = lib.new_cell("NRing_type2"+str(jj)+str(i))
            
            gap = (1+i)*min_gap

            ring_radius = [ring_radius1 if jj==1 else ring_radius2][0]            
            ring = [ring1 if jj ==1 else ring2][0]

            centro2 = ring_radius + gap + 0.5*width
            
            bend_radius = ring_radius
            posy_taper = posy_taper = 2.5*ring_radius
            dist_entre_filas = ring_radius*3.4                                              
            corrimiento_en_x = 2.5*ring_radius       
            distancia_curv = 1.17*dist_entre_filas
            Offset_ring = 0.5*(distancia_curv + corrimiento_en_x)
            posicion_y_en_el_chip = 1.1*dist_entre_filas*series[jj]                        
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
            
            
            Ring(centro1,centro2,ring_radius,width,ring,i)           
            c.add(gdspy.CellReference(ring, (i* distancia_curv + Offset_ring-
                                             chip_size*0.5, # X of the resonators
                                              dist_wg+posicion_y_en_el_chip+offset_from_zero)))
            
            if jj == 1: 
                Taper(taper_right1,2.5*ring_radius1)
                c.add(gdspy.CellArray(taper_right1, 1, len(n_resonators), \
                              (0, -dist_waveguides), \
                              (chip_size*0.5, posicion_y_en_el_chip)))
            else:
                Taper(taper_right2,2.5*ring_radius2)
                c.add(gdspy.CellArray(taper_right2, 1, len(n_resonators), \
                              (0, -dist_waveguides), \
                              (chip_size*0.5, posicion_y_en_el_chip)))
            
            text_l = gdspy.Text("N %s, Gap = %s" % (i, gap), 
                       10, (-chip_size*.5 + 0.5*taper_len + i*50, dist_wg+posicion_y_en_el_chip+offset_from_zero+10), layer = 4)
    
            c.add(text_l)
            
            text_r= gdspy.Text("N %s, Gap = %s" % (i, gap), 
                       10, (chip_size*.5 - taper_len -  i*50, posicion_y_en_el_chip+\
                            posy_taper + offset_from_zero + 10 + dist_wg), layer = 4)     
            c.add(text_r)
                         
    ### Taper 1    
        c.add(gdspy.CellArray(taper_left, 1, len(n_resonators ),\
                              (0, -dist_waveguides),\
                              (-chip_size*0.5, posicion_y_en_el_chip)))

    # Save to a gds file and check out the output
    lib.write_gds("Double_ring.gds")