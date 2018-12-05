from imports import *
from functions import *
from parameters import *
import PAA_LISA
import NOISE_LISA

class WFE():
    def __init__(self,**kwargs):
        global labda, w0, E0, k, D, LA
        
        
        self.Ndata = kwargs.pop('Ndata',False)
        self.tele_control = kwargs.pop('tele_control','full control')
        self.side = kwargs.pop('side','l')
        self.speed_on = kwargs.pop('speed_on',True)
        self.simple=True
        self.jitter=[False,False]
        if self.Ndata==False:
            print('Please input NOISE_LISA object')
        else:
            labda = self.Ndata.data.labda
            w0 = self.Ndata.data.w0
            E0 = 1 #...adjust
            k = (2*np.pi)/labda
            D = self.Ndata.data.D
            LA = PAA_LISA.utils.la()
            self.pupil()
            self.scale=1
 

    def pupil(self,D_calc = D,Nbins=10):
        self.xlist = np.linspace(-D_calc,D_calc,Nbins)
        self.ylist = self.xlist
        self.Deltax = self.xlist[1]-self.xlist[0]
        self.Deltay = self.ylist[1]-self.ylist[0]
        self.Nbinsx = len(self.xlist)
        self.Nbinsy = len(self.ylist)

    def w(self,z):
        zR = np.pi*(w0**2)/labda

        return w0*((1+((z/zR)**2))**0.5)

    def R(self,z,guess=False):
        zR = np.pi*(w0**2)/labda

        if guess==False:
            return abs(z*(1+((zR/z)**2)))

        elif guess==True:
            return z

    def z_solve(self,x,y,z,calc_R=False):
        
        R_new = lambda z_tot: self.R(z_tot,guess=False)
        f_solve = lambda z_tot: (R_new(z_tot) - (R_new(z_tot)**2 - (x**2+y**2))**0.5) - (z_tot-z)
        z_sol = scipy.optimize.brentq(f_solve,z,z+1)
        if calc_R==True:
            return R(z_sol)
        else:
            return z_sol

    ### Gaussian beam
    def phi_gauss(self,i,t,dX,dY,side='default',print0 = False):
        dx_list = self.xlist
        dy_list = self.ylist
        Ndata = self.Ndata
        
        [i_self,i_left,i_right] = PAA_LISA.utils.i_slr(i)

        n = Ndata.data.n_func(i,t)
        if print0==True:
            print ('Type of telescope control is: ' + self.tele_control)
            
        if side=='default':
            side = self.side

        if side=='l':
            if self.tele_control=='full control':
                tele = Ndata.tele_l_fc(i,t)
            elif self.tele_control=='no control':
                tele = Ndata.tele_l(i,t)
            elif self.tele_control=='SS':
                r = Ndata.data.r_func(i,t)
                ang_SS = self.tele_SS_l(i,t)
                tele = LA.unit(LA.rotate(r,n,ang_SS))

            beam = Ndata.data.u_l_func_tot(i,t) #...adjust for sending telescope --> PAAM control
        elif side=='r':
            if self.tele_control=='full control':
                tele = Ndata.tele_r_fc(i,t)
            elif self.tele_control=='no control':
                tele = Ndata.tele_r(i,t)
            elif self.tele_control=='SS':
                r = Ndata.data.r_func(i,t)
                ang_SS = self.tele_SS_r(i,t)
                tele = LA.unit(LA.rotate(r,n,ang_SS))
            beam = Ndata.data.u_r_func_tot(i,t)

        beam = self.scale*beam
        [ang_x,ang_y] = LA.beam_ang(beam,tele,n)

        psi = 0 #Gouy angle
        if self.jitter[0]==False:
            Zm = 0 #Exit pupil abberation
        else:
            Zm = self.jitter[0][1](t)

        dX_ac = np.cos(ang_x)*dX
        dY_ac = np.cos(ang_y)*dY
        dZ_ac = np.sin(ang_x)*dX+np.sin(ang_y)*dY
        dR_ac = (dX_ac**2+dY_ac**2)**0.5
        d_ac = self.z_solve(dX_ac,dY_ac,np.linalg.norm(beam)+dZ_ac)
        
        
        Deltax = self.Deltax
        Deltay = self.Deltay
        Nbinsx = self.Nbinsx
        Nbinsy = self.Nbinsy
 
        if self.speed_on==True:
            dx_list = [0]
            dy_list = [0]
            ps = np.empty((1,1),dtype = np.complex64)
        else:
            ps = np.empty((Nbinsx,Nbinsy),dtype = np.complex64)

        for i in range(0,len(dx_list)):
            for j in range(0,len(dy_list)):
                dx = dx_list[i]
                dy = dy_list[j]
                dr = ((dx**2)+(dy**2))**0.5
                if dr<=D:
                    s = (d_ac**2+(dR_ac-dr)**2)**0.5
                    E = E0*np.exp((-(dr**2))/(w0**2))*np.exp(1j*psi)
                    #print(E,E0,dr,psi)
                    R_new = d_ac
                    #R_new = np.linalg.norm(vec)
                    #print(np.exp(1j*((2*np.pi)/labda)*(Zm+s))/s)
                    ps[i,j] = (E*np.exp(1j*((2*np.pi)/labda)*(Zm+s))/s)*Deltax*Deltay*np.exp((-1j*2*np.pi*R_new)/labda)
                else:
                    ps[i,j] = np.nan
        
        diff_inte = np.nansum(ps)
        phi = np.angle(diff_inte)
        diff_inte_mag = (diff_inte.real**2+diff_inte.imag**2)**0.5
        #phi = np.log(diff_inte_new/diff_inte_mag).imag
            
        return diff_inte_mag,phi,ps

    def aperture(self,xlist,ylist,function,dType=np.complex64):
        print('Type of telescope control is: ' + self.tele_control)

        Nbins = len(xlist)
        D = max(xlist)
        ps = np.empty((Nbins,Nbins),dtype=dType)
        for i in range(0,len(xlist)):
            for j in range(0,len(ylist)):
                x = xlist[i]
                y = ylist[j]
                r = (x**2+y**2)**0.5
                if r<=D:
                    ps[i,j] = function(x,y)
                else:
                    ps[i,j] = np.nan
        
        return ps

    def phase_rec(self,i,t,side='default'):
        if self.simple == True:
            print('Simple mode is on, only calculating for center of receiving telescope')
            phi = self.phi_gauss(i,t,0,0,side=side,print0=True)[1]
        else:
            print('Simple mode is off, calculating over whole aperture receiving telescope')
            phi = self.aperture(self.xlist,self.ylist, lambda dx,dy: self.phi_gauss(i,t,dx,dy,side=side)[1])

        return phi
    
    def I_rec(self,i,t):
        if self.simple == True:
            I = self.phi_gauss(i,t,0,0)[0]
        else:
            I = self.aperture(self.xlist,self.ylist, lambda dx,dy: self.phi_gauss(i,t,dx,dy)[0])

        return I

    def plot_aperture(self,ps_list,mode = 'Phi',figs=[1,1],vlim=True):
        def plot_plot(ax_sel,ps):
            if mode == 'Phi' and vlim==True:
                p1 = ax_sel.matshow(ps,cmap='magma',vmin=-np.pi,vmax=np.pi)
            elif mode == 'I' or vlim==False:
                p1 = ax_sel.matshow(ps,cmap='magma')
            ax_sel.set_title(mode)
            f.colorbar(p1,ax=ax_sel)

            return 0 
        
        
        f,ax = plt.subplots(figs[0],figs[1])
        xlist = self.xlist
        ylist = self.ylist
        
        count = 0 
        try:
            for i in len(ax):
                try:
                    for j in len(ax[i]):
                        plot_plot(ax[i,j],ps_list_count)
                        count=count+1
                except:
                    plot_plot(ax[i],ps_list[count])
                    count=count+1
                    break
        except:
            plot_plot(ax,ps_list[count])
            count=count+1
            pass

        return f

    def jitter_tele(self,N,t_end,psd_h,psd_v):
        jitter_h = self.Ndata.Noise_time(0.0001,0.1,4096,psd_h,t_end)
        jitter_v = self.Ndata.Noise_time(0.0001,0.1,4096,psd_v,t_end)

        self.jitter=[jitter_h,jitter_v]

        return 0

    
    def tele_control_noise(self,i,step_max=False,dt=1,side='l'):

        t_vec = self.Ndata.t_all
        d30 = np.radians(30)
        #b0 = self.tele_l_fc(1,t0)
        b0 = d30
        sys = tf([np.radians(b0)], [1,20000,1]) # Define a transfer function C/R = 2/(s^2+2s+2)
        
        if step_max==False:
            t_calc = t_vec
        else:
            t_calc = t_vec[0:step_max] #... adjust for longer/whole time period
        
        T = np.linspace(0,t_vec[1]-t_vec[0],int(t_vec[1]-t_vec[0])/dt)

        yout=[]
        
        if side=='l':
            tele = self.Ndata.tele_l_fc
        elif side=='r':
            tele = self.Ndata.tele_r_fc
       
        T_all = []
        for t in t_calc[0:len(t_calc)-1]: 
            if t==t_calc[0]:
                X0 = LA.angle(tele(i,t),self.Ndata.data.r_func(i,t))/d30
            else:
                X0 = (yout_calc[-1])/d30
            U = LA.angle(tele(i,t+T[-1]),self.Ndata.data.r_func(i,t+T[-1]))/d30
            U = U*1
            
            K = d30*X0
            #sys = tf([K], [1,20000,1]) # Define a transfer function C/R = 2/(s^2+2s+2)
            
            yout_calc = lsim(sys,T=T,U=U,X0 = X0)[0]
            #plt.plot(yout_calc)
            yout.append(np.ndarray.tolist(yout_calc))
            T_all.append(np.ndarray.tolist(T+t))
        
        
        for i in range(0,len(yout)):
            x = T+i*T[-1]
            #plt.plot(x,yout[i])

        yout = LA.flatten(yout)
        T_all = LA.flatten(T_all)

        
        self.tele_control = 'SS'
        f = interp1d(T_all,yout,bounds_error=False)
        self.tele_Colntrol_func = f
        
        return [[T_all,yout],f]

    def tele_control_ss(self,step_max=False,dt=1):
        tele_SS_l = []
        tele_SS_r = []
        
        for i in range(1,4):
            tele_SS_l.append(self.tele_control_noise(i,side='l',step_max=step_max,dt=dt)[1])
            tele_SS_r.append(self.tele_control_noise(i,side='r',step_max=step_max,dt=dt)[1])
        
        self.tele_SS_l = PAA_LISA.utils.func_over_sc(tele_SS_l)
        self.tele_SS_r = PAA_LISA.utils.func_over_sc(tele_SS_r)

        return 0


