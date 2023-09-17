import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import complex_ode,solve_ivp
import matplotlib.ticker as ticker
import matplotlib.colors as mcolors
from scipy.constants import pi, c, hbar
from matplotlib.widgets import Slider, Button, TextBox
from matplotlib.animation import FuncAnimation
import matplotlib.image as mpimg
from scipy.optimize import curve_fit
import time
from scipy.sparse import block_diag,identity,diags
import scipy.integrate as integrate

class Resonator:
    def __init__(self, resonator_parameters):
        #Physical parameters initialization
        self.n0 = resonator_parameters['n0']
        self.n2 = resonator_parameters['n2']
        self.FSR = resonator_parameters['FSR']
        self.w0 = resonator_parameters['w0']
        self.width = resonator_parameters['width']
        self.height = resonator_parameters['height']
        self.kappa_0 = resonator_parameters['kappa_0']
        self.kappa_ex_arr = np.fft.ifftshift(resonator_parameters['kappa_ex'])
        self.kappa_ex = self.kappa_ex_arr[0]
        self.Dint = np.fft.ifftshift(resonator_parameters['Dint'])
        #Auxiliary physical parameters
        self.Tr = 1/self.FSR #round trip time
        self.Aeff = self.width*self.height 
        self.Leff = c/self.n0*self.Tr 
        self.Veff = self.Aeff*self.Leff 
        self.g0 = hbar*self.w0**2*c*self.n2/self.n0**2/self.Veff
        self.gamma = self.n2*self.w0/c/self.Aeff
        self.kappa = self.kappa_0 + self.kappa_ex
        self.N_points = len(self.Dint)
        mu = np.fft.fftshift(np.arange(-self.N_points/2, self.N_points/2))
        def func(x, a, b, c, d):
            return a + x*b + c*x**2/2 + d*x**3/6
        popt, pcov = curve_fit(func, mu, self.Dint)
        self.D2 = popt[2]
        self.D3 = popt[3]

#        print(self.kappa/2/self.D2)

        
    def noise(self, a):
#        return a*np.exp(1j*np.random.uniform(-1,1,self.N_points)*np.pi)
        return a*(np.random.uniform(-1,1,self.N_points) + 1j*np.random.uniform(-1,1,self.N_points))

    #   Propagate Using the Step Adaptive  Method
    def Propagate_SAM(self, simulation_parameters, Pump, Seed=[0]):
        start_time = time.time()

        T = simulation_parameters['slow_time']
        abtol = simulation_parameters['absolute_tolerance']
        reltol = simulation_parameters['relative_tolerance']
        out_param = simulation_parameters['output']
        nmax = simulation_parameters['max_internal_steps']
        detuning = simulation_parameters['detuning_array']
        eps = simulation_parameters['noise_level']
        
        pump = Pump*np.sqrt(1./(hbar*self.w0))
        if Seed[0] == 0:
            seed = self.seed_level(Pump, detuning[0])*np.sqrt(2*self.g0/self.kappa)
        else:
            seed = Seed*np.sqrt(2*self.g0/self.kappa)
        ### renarmalization
        T_rn = (self.kappa/2)*T
        f0 = pump*np.sqrt(8*self.g0*self.kappa_ex/self.kappa**3)
        print('f0^2 = ' + str(np.round(max(abs(f0)**2), 2)))
        print('xi [' + str(detuning[0]*2/self.kappa) + ',' +str(detuning[-1]*2/self.kappa)+ ']')
        noise_const = self.noise(eps) # set the noise level
        nn = len(detuning)
        ### define the rhs function
        def LLE_1d(Time, A):
            A = A - noise_const#self.noise(eps)
            A_dir = np.fft.ifft(A)*len(A)## in the direct space
            dAdT =  -1*((self.kappa_ex_arr+self.kappa_0)/self.kappa + 1j*(self.Dint + dOm_curr)*2/self.kappa)*A + 1j*np.fft.fft(A_dir*np.abs(A_dir)**2)/len(A) + f0#*len(A)
            return dAdT
        
        t_st = float(T_rn)/len(detuning)
        r = complex_ode(LLE_1d).set_integrator('dop853', atol=abtol, rtol=reltol,nsteps=nmax)# set the solver
        r.set_initial_value(seed, 0)# seed the cavity
        sol = np.ndarray(shape=(len(detuning), self.N_points), dtype='complex') # define an array to store the data
        sol[0,:] = seed
        #printProgressBar(0, nn, prefix = 'Progress:', suffix = 'Complete', length = 50, fill='elapsed time = ' + str((time.time() - start_time)) + ' s')
        for it in range(1,len(detuning)):
            self.printProgressBar(it + 1, nn, prefix = 'Progress:', suffix = 'Complete,', time='elapsed time = ' + '{:04.1f}'.format(time.time() - start_time) + ' s', length = 50)
            #self.print('elapsed time = ', (time.time() - start_time))
            dOm_curr = detuning[it] # detuning value
            sol[it] = r.integrate(r.t+t_st)
            
        if out_param == 'map':
            return sol
        elif out_param == 'fin_res':
            return sol[-1, :] 
        else:
            print ('wrong parameter')
       
    def Propagate_SplitStep(self, simulation_parameters, Pump, Seed=[0], dt=1e-3):
        start_time = time.time()
        T = simulation_parameters['slow_time']
        out_param = simulation_parameters['output']
        detuning = simulation_parameters['detuning_array']
        eps = simulation_parameters['noise_level']
        TRaman = simulation_parameters['Raman_term']
        fR = simulation_parameters['fraction_Raman']
        #dt = simulation_parameters['time_step']#in photon lifetimes
        
        pump = Pump*np.sqrt(1./(hbar*self.w0))
        if Seed[0] == 0:
            seed = self.seed_level(Pump, detuning[0])*np.sqrt(2*self.g0/self.kappa)
        else:
            seed = Seed*np.sqrt(2*self.g0/self.kappa)
        ### renarmalization
        T_rn = (self.kappa/2)*T # T is the roundtrip time. T_rn is the slow time

        f0 = pump*np.sqrt(8*self.g0*self.kappa_ex/self.kappa**3)
        print('f0^2 = ' + str(np.round(max(abs(f0)**2), 2)))
        print('xi [' + str(detuning[0]*2/self.kappa) + ',' +str(detuning[-1]*2/self.kappa)+ ']')
        noise_const = self.noise(eps) # set the noise level
        nn = len(detuning)
        
        t_st = float(T_rn)/len(detuning)
        #dt=1e-4 #t_ph
        
        sol = np.ndarray(shape=(len(detuning), self.N_points), dtype='complex') # define an array to store the data
#        Raman_contribution = sol
        sol[0,:] = (seed)
        self.printProgressBar(0, nn, prefix = 'Progress:', suffix = 'Complete', length = 50)
        for it in range(1,len(detuning)):
            
            self.printProgressBar(it + 1, nn, prefix = 'Progress:', suffix = 'Complete,', time='elapsed time = ' + '{:04.1f}'.format(time.time() - start_time) + ' s', length = 50)
            dOm_curr = detuning[it] # detuning value
            t=0
            buf = sol[it-1]
            buf-=noise_const
#            Tint = np.arange(0,t_st,dt)
#            Tint = np.linspace(0,t_st,len(self.kappa_ex_arr))
            while t<t_st:
                buf_dir = np.fft.ifft(buf)## in the direct space
                f = np.fft.ifft(f0)*len(buf)
                ## First step
                #buf =buf + dt*(1j/len(buf)*np.fft.fft(buf_dir*np.abs(buf_dir)**2))
                #buf =np.fft.fft(np.exp(dt*(1j/len(buf)*np.fft.fft(buf_dir*np.abs(buf_dir)**2) ))*buf_dir)
                ## second step
                #buf = np.exp(-dt *(1+1j*(self.Dint + dOm_curr)*2/self.kappa + f0/buf)) * buf
                '''
                Nonlinear part of the Lugiato Lefever equation.
                The nonlinear part is the one that is solved in the TIME domain 
                It is necessary to Fourier transform back and forth because the
                linear step is made in the frequency domain while the nonlinear
                step is made in the time domain
                buf_dir = A 
                '''
                
#                D1      = 12.13e9 * 2 * np.pi / self.kappa  ## ??????
#                fR      = 0.4                         # 20%
#                TRaman  = Raman
#                D1 = self.FSR
#                tr = 2*np.pi/D1                
                
                field_fft = np.fft.ifft(np.abs(buf_dir) ** 2)
                omega = 2.*np.pi*np.fft.fftfreq(len(buf_dir),dt) # Siempre la misma
                out_field = np.fft.fft(-1j*omega*field_fft)
#                TRaman = 0
                ad_coeff = np.sqrt(self.kappa * 0.5) * (self.g0 * np.sqrt(self.D2))**(-1)
                # print(ad_coeff)
                ### Raman Term  2*np.pi*
                if TRaman != 0.0:                    
                    # Raman =  ad_coeff * 1j * self.g0 * fR * (- 2*np.pi* out_field * TRaman/T) #tr
                    Raman =   ad_coeff* 1j * self.g0 * fR * (- out_field * TRaman/T) #tr
                else:
                    Raman = 0.0                     
    
                buf_dir = np.fft.fft(np.exp(dt *(1j * np.abs(buf_dir)**2 + Raman + f/buf_dir )) * buf_dir)
                
#                Original = buf_dir = np.fft.fft(np.exp(dt *(1j *  np.abs(buf_dir) ** 2 + f/buf_dir )) * buf_dir)
            
                # Second step: Everything is in the Fourier space
                buf = np.exp(-dt *((self.kappa_ex_arr+self.kappa_0)/self.kappa+1j*(self.Dint + dOm_curr)*2/self.kappa )) * buf_dir
                #buf = np.fft.fft(buf_dir)
                t+=dt #Fast time
            sol[it] = buf

#            Raman_contribution[it] = Raman
            
        if out_param == 'map':
            return sol
        elif out_param == 'fin_res':
            return sol[-1, :] 
        else:
            print ('wrong parameter')
        
    def seed_level (self, pump, detuning):
        f_norm = pump*np.sqrt(1./(hbar*self.w0))*np.sqrt(8*self.g0*self.kappa_ex/self.kappa**3)
        detuning_norm  = detuning*2/self.kappa
        stat_roots = np.roots([1, -2*detuning_norm, (detuning_norm**2+1), -abs(f_norm[0])**2])
        ind_roots = [np.imag(ii)==0 for ii in stat_roots]
        res_seed = np.zeros_like(f_norm)
        res_seed[0] = abs(stat_roots[ind_roots])**.5/np.sqrt(2*self.g0/self.kappa)
        return res_seed
    
    def seed_soliton(self, pump, detuning):
        fast_t = np.linspace(-pi,pi,len(pump))*np.sqrt(self.kappa/2/self.D2)
        f_norm = np.sqrt(pump/(hbar*self.w0))*np.sqrt(8*self.g0*self.kappa_ex/self.kappa**3)
        detuning_norm  = detuning*2/self.kappa
        stat_roots = np.roots([1, -2*detuning_norm, (detuning_norm**2+1), -abs(f_norm[0])**2])
        
        ind_roots = [np.imag(ii)==0 for ii in stat_roots]
        B = np.sqrt(2*detuning_norm)
        return np.fft.fft(np.min(np.abs(stat_roots[ind_roots]))**.5 + B*np.exp(1j*np.arccos(2*B/np.pi/f_norm[0])*2)*np.cosh(B*fast_t)**-1)/np.sqrt(2*self.g0/self.kappa)/len(pump)
                
    def NeverStopSAM (self, T_step, detuning_0=-1, Pump_P=2., nmax=1000, abtol=1e-10, reltol=1e-9, out_param='fin_res'):
        self.Pump = self.Pump/abs(self.Pump)
        def deriv_1(dt, field_in):
        # computes the first-order derivative of field_in
            field_fft = np.fft.fft(field_in)
            omega = 2.*np.pi*np.fft.fftfreq(len(field_in),dt)
            out_field = np.fft.ifft(-1j*omega*field_fft)
            return out_field
        
        def deriv_2(dt, field_in):
        # computes the second-order derivative of field_in
            field_fft = np.fft.fft(field_in)
            omega = 2.*np.pi*np.fft.fftfreq(len(field_in),dt)
            field_fft *= -omega**2
            out_field = np.fft.ifft(field_fft)
            return out_field 
        
        def disp(field_in,Dint_in):
        # computes the dispersion term in Fourier space
            field_fft = np.fft.fft(field_in)
            out_field = np.fft.ifft(Dint_in*field_fft)     
            return out_field

        ### define the rhs function
        def LLE_1d(Z, A):
            # for nomalized
            if np.size(self.Dint)==1 and self.Dint == 1:
                 dAdt2 = deriv_2(self.TimeStep, A)
                 dAdT =  1j*dAdt2/2 + 1j*self.gamma*self.L/self.Tr*np.abs(A)**2*A - (self.kappa/2+1j*dOm_curr)*A + np.sqrt(self.kappa/2/self.Tr)*self.Pump*Pump_P**.5
            elif np.size(self.Dint)==1 and self.Dint == -1:
                 dAdt2 = deriv_2(self.TimeStep, A)
                 dAdT =  -1j*dAdt2/2 + 1j*self.gamma*self.L/self.Tr*np.abs(A)**2*A - (self.kappa/2+1j*dOm_curr)*A + np.sqrt(self.kappa/2/self.Tr)*self.Pump*Pump_P**.5
            else:  
                # with out raman
                Disp_int = disp(A,self.Dint)
                if self.Traman==0:
                    dAdT =  -1j*Disp_int + 1j*self.gamma*self.L/self.Tr*np.abs(A)**2*A - (self.kappa/2+1j*dOm_curr)*A + np.sqrt(self.kappa/2/self.Tr)*self.Pump*Pump_P**.5
                else:
                    # with raman
                    dAAdt = deriv_1(self.TimeStep,abs(A)**2)
                    dAdT =  -1j*Disp_int + 1j*self.gamma*self.L/self.Tr*np.abs(A)**2*A - (self.kappa/2+1j*dOm_curr)*A -1j*self.gamma*self.Traman*dAAdt*A + np.sqrt(self.kappa/2/self.Tr)*self.Pump*Pump_P**.5
            return dAdT
        
        r = complex_ode(LLE_1d).set_integrator('dopri5', atol=abtol, rtol=reltol,nsteps=nmax)# set the solver
        r.set_initial_value(self.seed, 0)# seed the cavity
                
        img = mpimg.imread('phase_space.png')
        xx = np.linspace(-1,5,np.size(img,axis=1))
        yy = np.linspace(11,0,np.size(img,axis=0))
        XX,YY = np.meshgrid(xx,yy)
                
        fig = plt.figure(figsize=(11,7))        
        plt.subplots_adjust(top=0.95,bottom=0.1,left=0.06,right=0.986,hspace=0.2,wspace=0.16)

        ax1 = plt.subplot(221)
        ax1.pcolormesh(XX,YY,img[:,:,1])
        plt.xlabel('Detuning')
        plt.ylabel('f^2')
        plt.title('Choose the region')
        plt.xlim(min(xx),max(xx))
        dot = plt.plot(detuning_0, Pump_P,'rx')
        
        
        ax2 = plt.subplot(222)
        line, = plt.plot(abs(self.seed)**2)
        plt.ylim(0,1.1)
        plt.ylabel('$|\Psi|^2$')
        
        ax3 = plt.subplot(224)
        line2, = plt.semilogy(self.mu, np.abs(np.fft.fft(self.seed))**2)
        plt.ylabel('PSD')
        plt.xlabel('mode number')
        ### widjets
        axcolor = 'lightgoldenrodyellow'

        resetax = plt.axes([0.4, 0.025, 0.1, 0.04])
        button = Button(resetax, 'Stop', color=axcolor, hovercolor='0.975')
        
        axboxf = plt.axes([0.1, 0.35, 0.1, 0.075])
        text_box_f = TextBox(axboxf, 'f^2', initial=str(Pump_P))
        
        axboxd = plt.axes([0.1, 0.25, 0.1, 0.075])
        text_box_d = TextBox(axboxd, 'Detuning', initial=str(detuning_0))
        
        Run = True
        def setup(event): 
            global Run
            Run = False   
        button.on_clicked(setup)
        
        def onclick(event): 
            if event.inaxes == ax1:
                ix, iy = event.xdata, event.ydata     
                text_box_d.set_val(np.round(ix,4))
                text_box_f.set_val(np.round(iy,4))
                ax1.plot([ix],[iy],'rx')
   

        fig.canvas.mpl_connect('button_press_event', onclick)
        
        while Run:
            dOm_curr = float(text_box_d.text) # get the detuning value
            Pump_P = float(text_box_f.text)
            Field = r.integrate(r.t+T_step)
            F_mod_sq = np.abs(Field)**2
            F_sp = np.abs(np.fft.fft(Field))**2
            line.set_ydata(F_mod_sq)
            line2.set_ydata(F_sp)
            ax2.set_ylim(0, max(F_mod_sq))
            ax3.set_ylim(min(F_sp),max(F_sp))
            plt.pause(1e-10)
        
    def printProgressBar (self, iteration, total, prefix = '', suffix = '', time = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
        """
        Call in a loop to create terminal progress bar
        @params:
            iteration   - Required  : current iteration (Int)
            total       - Required  : total iterations (Int)
            prefix      - Optional  : prefix string (Str)
            suffix      - Optional  : suffix string (Str)
            decimals    - Optional  : positive number of decimals in percent complete (Int)
            length      - Optional  : character length of bar (Int)
            fill        - Optional  : bar fill character (Str)
            printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
        """
        percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
        filledLength = int(length * iteration // total)
        bar = fill * filledLength + '-' * (length - filledLength)
        print('\r%s |%s| %s%% %s %s' % (prefix, bar, percent, suffix, time), end = printEnd)
        # Print New Line on Complete
        if iteration == total: 
                print()


              
            
def Plot_Map(map_data, zmin, zmax,dt=1,dz=1,colormap = 'cubehelix',z0=0):
    def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
        '''
        Function to offset the "center" of a colormap. Useful for
        data with a negative min and positive max and you want the
        middle of the colormap's dynamic range to be at zero
    
        Input
        -----
          cmap : The matplotlib colormap to be altered
          start : Offset from lowest point in the colormap's range.
              Defaults to 0.0 (no lower ofset). Should be between
              0.0 and `midpoint`.
          midpoint : The new center of the colormap. Defaults to 
              0.5 (no shift). Should be between 0.0 and 1.0. In
              general, this should be  1 - vmax/(vmax + abs(vmin))
              For example if your data range from -15.0 to +5.0 and
              you want the center of the colormap at 0.0, `midpoint`
              should be set to  1 - 5/(5 + 15)) or 0.75
          stop : Offset from highets point in the colormap's range.
              Defaults to 1.0 (no upper ofset). Should be between
              `midpoint` and 1.0.
        '''
        cdict = {
            'red': [],
            'green': [],
            'blue': [],
            'alpha': []
        }
    
        # regular index to compute the colors
        reg_index = np.linspace(start, stop, 257)
    
        # shifted index to match the data
        shift_index = np.hstack([
            np.linspace(0.0, midpoint, 128, endpoint=False), 
            np.linspace(midpoint, 1.0, 129, endpoint=True)
        ])
    
        for ri, si in zip(reg_index, shift_index):
            r, g, b, a = cmap(ri)
    
            cdict['red'].append((si, r, r))
            cdict['green'].append((si, g, g))
            cdict['blue'].append((si, b, b))
            cdict['alpha'].append((si, a, a))
    
        newcmap = mcolors.LinearSegmentedColormap(name, cdict)
        plt.register_cmap(cmap=newcmap)
    
        return newcmap


    def onclick(event):
        ix, iy = event.xdata, event.ydata
        x = int(np.floor(ix/dz))
        plt.suptitle('Chosen distance z = %f km'%ix, fontsize=20)
        ax.lines.pop(0)
        ax.plot([ix,ix], [0, dt*np.size(map_data,1)],'r')

        ax2 = plt.subplot2grid((4, 1), (2, 0))            
        ax2.plot(np.arange(0,dt*np.size(map_data,1),dt), abs(map_data[x,:])**2, 'r')
        ax2.set_ylabel('Power (W)')
        ax2.set_xlim(0, dt*np.size(map_data,1))        
        ax3 = plt.subplot2grid((4, 1), (3, 0))
        ax3.plot(np.arange(0,dt*np.size(map_data,1),dt), np.angle(map_data[x,:])/(np.pi),'b')
        if max( np.unwrap(np.angle(map_data[x,:]))/(np.pi)) - min( np.unwrap(np.angle(map_data[x,:]))/(np.pi))<10:
            ax3.plot(np.arange(0,dt*np.size(map_data,1),dt), np.unwrap(np.angle(map_data[x,:]))/(np.pi),'g')
        ax3.set_xlabel('Time (ps)')
        ax3.set_ylabel('Phase (rad)')
        ax3.set_xlim(0, dt*np.size(map_data,1))
        ax3.yaxis.set_major_locator(ticker.MultipleLocator(base=1.0))
        ax3.yaxis.set_major_formatter(ticker.FormatStrFormatter('%g $\pi$'))
        ax3.grid(True)
        plt.show()
        f.canvas.draw()
        
    f = plt.figure()
    ax = plt.subplot2grid((4, 1), (0, 0), rowspan=2)
    plt.suptitle('Choose the coordinate', fontsize=20)
    f.set_size_inches(10,8)
    Z,T = np.meshgrid( np.arange(zmin,zmax+dz*10,dz), np.arange(0, dt*np.size(map_data,1),dt))
#    orig_cmap = plt.get_cmap('viridis')
#    colormap = shiftedColorMap(orig_cmap, start=0., midpoint=.5, stop=1., name='shrunk')
    pc = ax.pcolormesh(Z, T, abs(np.transpose(map_data))**2, cmap=colormap)
    ax.plot([0, 0], [0, dt*np.size(map_data,1)-dt], 'r')
    ax.set_xlabel('Detuning')
    ax.set_ylabel('Time (ps)')
    ax.set_ylim(0, dt*np.size(map_data,1))
    ax.set_xlim(zmin, zmax)
    ix=z0
    x = int(np.floor(ix/dz))
    plt.suptitle('Chosen distance z = %f km'%ix, fontsize=20)
    ax.lines.pop(0)
    
    ax.plot([ix,ix], [0, dt*np.size(map_data,1)],'r')

    ax2 = plt.subplot2grid((4, 1), (2, 0))            
    ax2.plot(np.arange(0,dt*np.size(map_data,1),dt), abs(map_data[x,:])**2, 'r')
    ax2.set_ylabel('Power (W)')
    ax2.set_xlim(0, dt*np.size(map_data,1))        
    ax3 = plt.subplot2grid((4, 1), (3, 0))
    ax3.plot(np.arange(0,dt*np.size(map_data,1),dt), np.angle(map_data[x,:])/(np.pi),'b')
    if max( np.unwrap(np.angle(map_data[x,:]))/(np.pi)) - min( np.unwrap(np.angle(map_data[x,:]))/(np.pi))<10:
        ax3.plot(np.arange(0,dt*np.size(map_data,1),dt), np.unwrap(np.angle(map_data[x,:]))/(np.pi),'g')
    ax3.set_xlabel('Time (ps)')
    ax3.set_ylabel('Phase (rad)')
    ax3.set_xlim(0, dt*np.size(map_data,1))
    ax3.yaxis.set_major_locator(ticker.MultipleLocator(base=1.0))
    ax3.yaxis.set_major_formatter(ticker.FormatStrFormatter('%g $\pi$'))
    ax3.grid(True)
#    f.colorbar(pc)
    plt.subplots_adjust(left=0.07, bottom=0.07, right=0.95, top=0.93, wspace=None, hspace=0.4)
    f.canvas.mpl_connect('button_press_event', onclick)
    
def Plot_MapNew(map_data, detuning,frequency, FSR, D2, FM2FM,colormap = 'cubehelix'):
    
    dOm = detuning[1]-detuning[0]
    dt=1
    FSR = FSR*1e-12
    D2 = D2*1e-12
   
   
    Num_of_modes = map_data[0,:].size
    mu = np.arange(-Num_of_modes/2,Num_of_modes/2)
    def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
        '''
        Function to offset the "center" of a colormap. Useful for
        data with a negative min and positive max and you want the
        middle of the colormap's dynamic range to be at zero
    
        Input
        -----
          cmap : The matplotlib colormap to be altered
          start : Offset from lowest point in the colormap's range.
              Defaults to 0.0 (no lower ofset). Should be between
              0.0 and `midpoint`.
          midpoint : The new center of the colormap. Defaults to 
              0.5 (no shift). Should be between 0.0 and 1.0. In
              general, this should be  1 - vmax/(vmax + abs(vmin))
              For example if your data range from -15.0 to +5.0 and
              you want the center of the colormap at 0.0, `midpoint`
              should be set to  1 - 5/(5 + 15)) or 0.75
          stop : Offset from highets point in the colormap's range.
              Defaults to 1.0 (no upper ofset). Should be between
              `midpoint` and 1.0.
        '''
        cdict = {
            'red': [],
            'green': [],
            'blue': [],
            'alpha': []
        }
    
        # regular index to compute the colors
        reg_index = np.linspace(start, stop, 257)
    
        # shifted index to match the data
        shift_index = np.hstack([
            np.linspace(0.0, midpoint, 128, endpoint=False), 
            np.linspace(midpoint, 1.0, 129, endpoint=True)
        ])
    
        for ri, si in zip(reg_index, shift_index):
            r, g, b, a = cmap(ri)
    
            cdict['red'].append((si, r, r))
            cdict['green'].append((si, g, g))
            cdict['blue'].append((si, b, b))
            cdict['alpha'].append((si, a, a))
    
        newcmap = mcolors.LinearSegmentedColormap(name, cdict)
        plt.register_cmap(cmap=newcmap)
    
        return newcmap


    def onclick(event):
        
        ix, iy = event.xdata, event.ydata
        x = int(np.floor((ix-detuning.min())/dOm))
        max_val = (abs(map_data[x,:])**2).max()
        plt.suptitle('Chosen detuning '+r'$\zeta_0$'+ '= %f'%ix, fontsize=20)
        ax.lines.pop(0)
        ax.plot([ix,ix], [-np.pi, np.pi ],'r')

        ax2 = plt.subplot2grid((5, 1), (2, 0))            
        ax2.plot(phi, abs(map_data[x,:])**2/max_val, 'r')
        ax2.set_ylabel('Intracavity power [a.u.]')
        ax2.set_xlim(-np.pi,np.pi)
        ax2.set_ylim(0,1)        
        ax3 = plt.subplot2grid((5, 1), (3, 0))
        ax3.plot(phi, np.angle(map_data[x,:])/(np.pi),'b')
#        if max( np.unwrap(np.angle(map_data[x,:]))/(np.pi)) - min( np.unwrap(np.angle(map_data[x,:]))/(np.pi))<10:
#            ax3.plot(np.arange(0,dt*np.size(map_data,1),dt), np.unwrap(np.angle(map_data[x,:]))/(np.pi),'g')
        ax3.set_xlabel(r'$\varphi$')
        ax3.set_ylabel('Phase (rad)')
        ax3.set_xlim(-np.pi,np.pi)
        ax3.yaxis.set_major_locator(ticker.MultipleLocator(base=1.0))
        ax3.yaxis.set_major_formatter(ticker.FormatStrFormatter('%g $\pi$'))
        ax3.grid(True)
        
#        mu_freq = mu
#        mu_freq0 = 0
        
        
        
        mu_freq     = frequency*1e-12 + FSR*mu * D2* mu**2
        mu_freq0    = frequency*1e-12
#        mu_freq     = c/(mu_freq*1e3) # Wavelength
#        mu_freq0     = c/194e3
        
        mu_freq = mu
        mu_freq0 = 0
        
        

        spectra = 10*np.log10(abs(np.fft.fftshift(np.fft.fft(map_data[x,:])))**2/(abs(np.fft.fft(map_data[x,:]))**2).max())
        
        Spectra = spectra - np.min(spectra) 
 
        center_of_mass_x = np.sum(mu_freq*Spectra)/np.sum(Spectra)
        
        
        ax4 = plt.subplot2grid((5, 1), (4, 0))            
        ax4.plot(mu_freq,spectra)
        min_spectrum = np.min(spectra)
        max_spectrum = np.max(spectra)
        ax4.text(mu_freq0 + 3,np.mean(spectra),str(np.round(mu_freq0-center_of_mass_x,4))+' THz')
        
         
        
        ax4.vlines(mu_freq0,min_spectrum-3,max_spectrum+3,'g', linestyle = '--')
        ax4.vlines(center_of_mass_x,min_spectrum-3,max_spectrum+3,'r', linestyle = '--')
        ax4.set_ylabel('Spectrum, dB')
#        ax4.set_xlim(mu.min(),mu.max())
        #ax4.set_ylim(-100,3)   
        plt.show()
        f.canvas.draw()
        ### ------------------------------------------------------- ###
        df_optical = mu_freq0-center_of_mass_x
        
        df_microwave = np.sqrt(10**(FM2FM/10))*(df_optical*1e12)
        print('%s\t%s'%('delta optical','delta f microwave'))
        print('%s\t%s'%(str(np.round(df_optical,6))+' THz',str(np.round(df_microwave*1e-6,6))+' MHz'))
        ### ------------------------------------------------------- ###
        
    
    f = plt.figure()
    ax = plt.subplot2grid((5, 1), (0, 0), rowspan=2)
    plt.suptitle('Choose the detuning', fontsize=20)
    f.set_size_inches(10,8)
    phi = np.linspace(-np.pi,np.pi,map_data[0,:].size)
#    orig_cmap = plt.get_cmap('viridis')
#    colormap = shiftedColorMap(orig_cmap, start=0., midpoint=.5, stop=1., name='shrunk')
    pc = ax.pcolormesh(detuning, phi, abs(np.transpose(map_data))**2, cmap=colormap)
    ax.plot([0, 0], [-np.pi, np.pi], 'r')
    ax.set_xlabel('Detuning')
    ax.set_ylabel(r'$\varphi$')
    ax.set_ylim(-np.pi, np.pi)
    ax.set_xlim(detuning.min(),detuning.max())
    ix=0
    
    x = int(np.floor((ix-detuning.min())/dOm))
    max_val = (abs(map_data[x,:])**2).max()
    plt.suptitle('Chosen detuning '+r'$\zeta_0$'+ '= %f km'%ix, fontsize=20)
    ax.lines.pop(0)
    
    ax.plot([ix,ix], [-np.pi, np.pi ],'r')
    
    ax2 = plt.subplot2grid((5, 1), (2, 0))            
    ax2.plot(phi,abs(map_data[x,:])**2/max_val, 'r')
    ax2.set_ylabel('Intracavity power [a.u.]')
    ax2.set_xlim(-np.pi,np.pi)
    ax2.set_ylim(0,1)        
    ax3 = plt.subplot2grid((5, 1), (3, 0))
    ax3.plot(phi, np.angle(map_data[x,:])/(np.pi),'b')
#    if max( np.unwrap(np.angle(map_data[x,:]))/(np.pi)) - min( np.unwrap(np.angle(map_data[x,:]))/(np.pi))<10:
#        ax3.plot(np.arange(0,dt*np.size(map_data,1),dt), np.unwrap(np.angle(map_data[x,:]))/(np.pi),'g')
    ax3.set_xlabel(r'$\varphi$')
    ax3.set_ylabel('Phase (rad)')
    ax3.set_xlim(-np.pi,np.pi)
    ax3.yaxis.set_major_locator(ticker.MultipleLocator(base=1.0))
    ax3.yaxis.set_major_formatter(ticker.FormatStrFormatter('%g $\pi$'))
    ax3.grid(True)
    ax4 = plt.subplot2grid((5, 1), (4, 0))            
    ax4.plot(mu,10*np.log10(abs(np.fft.fftshift(np.fft.fft(map_data[x,:])))**2/(abs(np.fft.fft(map_data[x,:]))**2).max()), 'b.')
    ax4.set_ylabel('Spectrum, dB')
#    ax4.set_xlim(mu.min(),mu.max())
    #ax4.set_ylim(-50,3)        
#    f.colorbar(pc)
    plt.subplots_adjust(left=0.07, bottom=0.07, right=0.95, top=0.93, wspace=None, hspace=0.4)
    f.canvas.mpl_connect('button_press_event', onclick)
"""
here is a set of useful standard functions
"""
if __name__ == '__main__':
    print('PyCORe')
    