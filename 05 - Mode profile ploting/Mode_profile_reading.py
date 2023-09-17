# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 09:41:11 2020

@author: Administrator
"""
# https://docs.jcmwave.com/JCMsuite/html/PythonInterface/bf647454e36069fd16f1a7a35cf6a865.html
# https://docs.jcmwave.com/JCMsuite/html/PythonInterface/57af15bb25dce719980d910940c2cd04.html?version=4.0.3#af15bb25dce719980d910940c2cd04
# https://docs.jcmwave.com/JCMsuite/html/PythonInterface/849f4e5b5a742e774b22bb4811574000.html
# https://scipython.com/blog/non-linear-least-squares-fitting-of-a-two-dimensional-data/

#   For exporting files on jcmWave
#   https://docs.jcmwave.com/JCMsuite/html/ParameterReference/b6dad0616b1396ef63dd0fe7dbc17d36.html?version=4.0.3

import sys
import os
sys.path.append(os.path.join(os.getenv('JCMROOT'), 'ThirdPartySupport', 'Python'))
import jcmwave
jcmwave.info()

def Mode_profile(plot_mode = True,
                 fit_mode = False,
                 plot_eqA = False,
                 other_plots = False):
    
    import numpy as np
    import matplotlib.pyplot as plt
    
    
    cfb = jcmwave.loadcartesianfields('project_results/fieldbag_cartesian.jcm')
    amplitude = cfb['field'][0]    
    intensity = (amplitude.conj()*amplitude).sum(2).real 
    
    np.savetxt('intensity.csv', intensity, delimiter=',')
    np.savetxt('cfb_x.csv', cfb['X'], delimiter=',')
    np.savetxt('cfb_y.csv', cfb['Y'], delimiter=',')
    
    if plot_mode == True:
         # Plot output of JCM Wave
        fig, ax = plt.subplots(1)
        ax.pcolormesh(cfb['X'], cfb['Y'], intensity, shading='gouraud') 
        ax.axis('tight')
        ax.set_aspect('equal')
        ax.xaxis.major.formatter.set_powerlimits((-1, 0))
        ax.yaxis.major.formatter.set_powerlimits((-1, 0))        
    else:
        ax = 'No mode plotted'
         
    
        
        
    if fit_mode == True:
        
    
        # We perform the fit    //  Normalization for simpler fitting?
        intensity_norm = intensity/np.max(intensity)
        
        # For a 2D Fit with Python the array must be reshaped
        intensity_ = intensity_norm.ravel()
        
        # define Gaussian
        def twoD_Gauss((x,y),amplitude,x0,y0,sigma_x,sigma_y):
            x0=float(x0)
            y0=float(y0)
            gauss2d = amplitude*np.exp(-(((x-x0)**(2)/(2*sigma_x**(2))) + ((y-y0)**(2)/(2*sigma_y**(2)))))
            return gauss2d.ravel()
        
        # Fit with gaussian
        import scipy.optimize as opt
        x = cfb['X']
        y = cfb['Y']
        
        initial_guess = (1 , 0, 0, .3e-6, 0.3e-6)
        popt, pcov = opt.curve_fit(twoD_Gauss, (x,y), intensity_, p0=initial_guess)

       
        ## We calculate the volume under the gaussian to get the equivalent area
        from scipy.integrate import dblquad
        import numpy as np
        def Gauss_vol(y, x):
            'y must be the first argument, and x the second.'
            return twoD_Gauss((x, y), *popt)
        
        volume, err = dblquad(Gauss_vol, x[0][0], x[-1][0],
                           lambda x: y[0][0],
                           lambda x: y[0][-1])
        print('Normalized volume is %s: ' % volume)
        
        # When normalized it will be equal to 1, obviously...
        height = np.max(intensity_norm)  
        
        Equivalent_area = volume/height
        print('Equivalent area is %s: ' % volume)
        x_eqA = np.sqrt(Equivalent_area) # Side of the equivalent square, for plotting 
        
    #    We are interested in the mode radius. We use the FWHM but this 
    #    criteria should be thought better. Maybe define 1/e radius
        
        FWHM_x = 2*np.sqrt(2*np.log(2))*popt[3]
        FWHM_y = 2*np.sqrt(2*np.log(2))*popt[4]
    
    #   dr1 = np.sqrt((FWHM_x*0.5)**2 + (FWHM_y*0.5)**2)  # The radius of the circle 
    
        dr = 0.5*FWHM_x
    
  
        # Plotting     
        if plot_eqA == True:
             # Plot output of JCM Wave
            plt.pcolormesh(cfb['X'], cfb['Y'], intensity, shading='gouraud') 
            plt.axis('tight')
            plt.gca().set_aspect('equal')
            plt.gca().xaxis.major.formatter.set_powerlimits((-1, 0))
            plt.gca().yaxis.major.formatter.set_powerlimits((-1, 0))
            
            # Plot the equivalent area as a square
            plt.axes()
            rectangle = plt.Rectangle((-x_eqA/2,-x_eqA/2), x_eqA, x_eqA, 
                                      fc='magenta',ec="red", alpha=0.5)
            plt.gca().add_patch(rectangle)
            plt.axis('scaled')
    
            # Plots a circle with radius dr1        
    #        plt.axes()
    #        Circle = plt.Circle((0,0),dr1, fc='cyan',ec="blue", alpha=0.5)
    #        plt.gca().add_patch(Circle)
    #        plt.axis('scaled')
            
            # Plots a Ellipse with x = FWHM_x and y = FWHM_y
            from matplotlib.patches import Ellipse
            plt.axes()
            Elipse = Ellipse((0,0),FWHM_x,FWHM_y, fc='cyan',ec="blue", alpha=0.5)
            plt.gca().add_patch(Elipse)
            plt.axis('scaled')
                
        #% Plotting other nice results    
        if other_plots == True:
            # This import registers the 3D projection, but is otherwise unused.
            from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
            from matplotlib import cm
            #from matplotlib.ticker import LinearLocator, FormatStrFormatter
        
            fig = plt.figure()
            ax = fig.gca(projection='3d')
            
            # Plot the surface.
            surf = ax.plot_surface(cfb['X'], cfb['Y'], intensity_norm, cmap=cm.coolwarm,
                                   linewidth=0, antialiased=False)
            # Add a color bar which maps values to colors.
            fig.colorbar(surf, shrink=0.5, aspect=5)
            plt.show()
            
            data_fitted = twoD_Gauss((cfb['X'], cfb['Y']), *popt)
            data_fitted = data_fitted.reshape(len(x), len(y))
            
            fig = plt.figure()
            ax = fig.gca(projection='3d')
            surf2 = ax.plot_surface(x, y, data_fitted, cmap=cm.coolwarm,
                                       linewidth=0, antialiased=False)
            fig.colorbar(surf2, shrink=0.5, aspect=5)
            plt.show()
    else:
        Equivalent_area =   'Not calculated'
        dr              =   'Not calculated'
    
    return Equivalent_area, dr, ax

