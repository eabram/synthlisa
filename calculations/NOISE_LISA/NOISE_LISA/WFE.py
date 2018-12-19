from imports import *
from functions import *
from parameters import *
import PAA_LISA
import NOISE_LISA

class WFE():
    def __init__(self,**kwargs):
        global labda, w0, E0, k, D, LA,L_tele
        
        
        self.Ndata = kwargs.pop('Ndata',False)
        self.tele_control = kwargs.pop('tele_control','full control')
        self.PAAM_control_method = kwargs.pop('PAAM_control','SS')
        self.side = kwargs.pop('side','l')
        self.speed_on = kwargs.pop('speed_on',True)
        self.simple=True
        self.jitter=[False,False]
        self.tele_aim_done = False
        self.jitter_tele_done = False
        self.tilt_send = False
        self.zmn_func = {}
        self.thmn_func = {}
        self.zmn={}
        self.thmn = {}
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
        self.xlist = np.linspace(-D_calc*0.5,D_calc*0.5,Nbins)
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
        z_sol = scipy.optimize.brentq(f_solve,0.5*z,z*1.5)
        if calc_R==True:
            return R(z_sol)
        else:
            return z_sol

    ### Gaussian beam
    def phi_gauss(self,i_self,t,dX,dY,side='l',scale=1,tilt_send = 'Default'):
        if tilt_send=='Default':
            tilt_send = self.tilt_send


        dx_list = self.xlist
        dy_list = self.ylist
        Ndata = self.Ndata
        Deltax = self.Deltax
        Deltay = self.Deltay
        Nbinsx = self.Nbinsx
        Nbinsy = self.Nbinsy

        [i_self,i_left,i_right] = PAA_LISA.utils.i_slr(i_self)

        n = Ndata.data.n_func(i_self,t)
        
        if side=='l':
            tdel = self.Ndata.data.L_rl_func_tot(i_self,t)
            n_left = Ndata.data.n_func(i_left,t-tdel)
            if self.tele_control=='full control':
                tele_rec = Ndata.tele_l_fc(i_self,t)
                tele_send = Ndata.tele_r_fc(i_left,tdel)
            elif self.tele_control=='no control':
                tele_rec = Ndata.tele_l(i_self,t)
                tele_send = Ndata.tele_r(i_left,t-tdel)
            elif self.tele_control=='SS':
                r = Ndata.data.r_func(i_self,t)
                ang_SS = self.tele_SS_l(i_self,t)
                ang_SS_left = self.tele_SS_r(i_left,t-tdel)
                tele_rec = LA.unit(LA.rotate(Ndata.tele_l(i_self,t),n,ang_SS))
                tele_send = LA.unit(LA.rotate(Ndata.tele_r(i_left,t-tdel),n_left,ang_SS_left))

            n_beam = self.Ndata.data.n_func(i_left,(t-tdel))
            beam = self.Ndata.PAA_point_r(i_left,t-tdel)#...adjust for sending telescope --> PAAM control
            v_pos = self.Ndata.data.v_r_func_tot(i_left,t-tdel)
            ps_send = self.phasefront_send(i_self,t,side='l')
            
        elif side=='r':
            tdel = self.Ndata.data.L_rr_func_tot(i,t)
            if self.tele_control=='full control':
                tele_rec = Ndata.tele_r_fc(i_self,t)
                tele_send = Ndata.tele_l_fc(i_right,t-tdel)
            elif self.tele_control=='no control':
                tele_rec = Ndata.tele_r(i_self,t)
                tele_send = Ndata.tele_l(i_right,t-tdel)
            elif self.tele_control=='SS':
                r = Ndata.data.r_func(i_self,t)
                ang_SS = self.tele_SS_r(i_self,t)
                ang_SS_right = self.tele_SS_l(i_right,t-tdel)
                
                tele_rec = LA.unit(LA.rotate(Ndata.tele_r(i_self,t),n,ang_SS))
                tele_send = LA.unit(LA.rotate(Ndata.tele_l(i_right,t-tdel),n_right,ang_SS_right))
            
            n_beam = self.Ndata.data.n_func(i_right,(t-tdel))
            beam = self.Ndata.PAA_point_l(i_right,t-tdel)*scale
            v_pos = self.Ndata.data.v_l_func_tot(i_right,t-tdel)*scale
            ps_send = self.phasefront_send(i_self,t,side='r')
        
        # Adding telescope jitter
        





        [ang_x,ang_y] = LA.beam_ang(beam,tele_rec,n)

        #print(LA.beam_ang(beam,tele_rec,n))
        E0 = 1
        psi=0 #Gouy angle

        # calculate offset
        offset = v_pos - beam
        [xoff,yoff,zoff] = LA.beam_coor(offset,LA.unit(-beam),n_beam)
        R = np.linalg.norm(v_pos)*scale
        
        dX_ac = np.cos(ang_x)*dX
        dY_ac = np.cos(ang_y)*dY
        dZ_ac = np.sin(ang_x)*dX+np.sin(ang_y)*dY
        dR_ac = (dX_ac**2+dY_ac**2)**0.5
        #d_ac = R - dZ_ac
        d_ac = self.z_solve(dX_ac,dY_ac,R+dZ_ac)

        
        ps=np.empty((Nbinsx,Nbinsy),dtype = np.complex64)
        for i in range(0,len(dx_list)):
            for j in range(0,len(dy_list)):
                dx = dx_list[i]
                dy = dy_list[j]
                dr = (dx**2+dy**2)**0.5
                if tilt_send==True:
                    Zn=ps_send[i,j]
                else:
                    Zn=0
                if dr<=0.5*D:
                    E = E0*np.exp(-(dr**2)/(w0**2))*np.exp(1j*psi)              
                    S = (R**2+((dX_ac+xoff)-dx_list[i])**2+((dY_ac+yoff)-dy_list[j])**2)**0.5
                    #print(S)
                    ps[i,j] = np.exp(1j*(2*np.pi/labda)*(Zn+S))*Deltax*Deltay # removed E/S
                else:
                    ps[i,j] = np.nan

            ret = np.nansum(ps)
            I = (np.real(ret)**2 + np.imag(ret**2))**0.5
            #TTL = np.log(np.imag(ret))
            TTL = np.angle(ret)*labda/(2*np.pi)

        return I,TTL,ps

        
        
        #dx_list = self.xlist
        #dy_list = self.ylist
        #Ndata = self.Ndata
        #
        #[i_self,i_left,i_right] = PAA_LISA.utils.i_slr(i)

        #n = Ndata.data.n_func(i,t)
        #if print0==True:
        #    print ('Type of telescope control is: ' + self.tele_control)
        #    
        #if side=='default':
        #    side = self.side

        #if side=='l':
        #    if self.tele_control=='full control':
        #        tele = Ndata.tele_l_fc(i,t)
        #    elif self.tele_control=='no control':
        #        tele = Ndata.tele_l(i,t)
        #    elif self.tele_control=='SS':
        #        r = Ndata.data.r_func(i,t)
        #        ang_SS = self.tele_SS_l(i,t)
        #        tele = LA.unit(LA.rotate(Ndata.tele_l(i,t),n,ang_SS))

        #    beam = Ndata.data.u_l_func_tot(i,t) #...adjust for sending telescope --> PAAM control
        #elif side=='r':
        #    if self.tele_control=='full control':
        #        tele = Ndata.tele_r_fc(i,t)
        #    elif self.tele_control=='no control':
        #        tele = Ndata.tele_r(i,t)
        #    elif self.tele_control=='SS':
        #        r = Ndata.data.r_func(i,t)
        #        ang_SS = self.tele_SS_r(i,t)
        #        tele = LA.unit(LA.rotate(Ndata.tele_l(i,t),n,ang_SS))
        #    beam = Ndata.data.u_r_func_tot(i,t)

        #beam = self.scale*beam
        #[ang_x,ang_y] = LA.beam_ang(beam,tele,n)

        #psi = 0 #Gouy angle
        #self.phasefront_send(i,t,side=side) #...adjust: with delay

        #Zm = self.ttl_s_tot #Exit pupil abberation
        #if self.jitter[0]!=False:
        #    jitter = self.jitter[0][1](t)
        #else:
        #    jitter = 0

        #dX_ac = np.cos(ang_x)*dX
        #dY_ac = np.cos(ang_y)*dY
        #dZ_ac = np.sin(ang_x)*dX+np.sin(ang_y)*dY
        #dR_ac = (dX_ac**2+dY_ac**2)**0.5
        #d_ac = self.z_solve(dX_ac,dY_ac,np.linalg.norm(beam)+dZ_ac)
        #
        #
        #Deltax = self.Deltax
        #Deltay = self.Deltay
        #Nbinsx = self.Nbinsx
        #Nbinsy = self.Nbinsy
 
        #if self.speed_on==True:
        #    dx_list = [0]
        #    dy_list = [0]
        #    ps = np.empty((1,1),dtype = np.complex64)
        #else:
        #    ps = np.empty((Nbinsx,Nbinsy),dtype = np.complex64)

        #for i in range(0,len(dx_list)):
        #    for j in range(0,len(dy_list)):
        #        dx = dx_list[i]
        #        dy = dy_list[j]
        #        dr = ((dx**2)+(dy**2))**0.5
        #        if dr<=D:
        #            Zm_calc = Zm[i,j]+jitter
        #            s = (d_ac**2+(dR_ac-dr)**2)**0.5
        #            #print(s)
        #            E = E0*np.exp((-(dr**2))/(w0**2))*np.exp(1j*psi)
        #            #print(E,E0,dr,psi)
        #            R_new = d_ac
        #            #R_new = np.linalg.norm(vec)
        #            #print(np.exp(1j*((2*np.pi)/labda)*(Zm+s))/s)
        #            ps[i,j] = (E*np.exp(1j*((2*np.pi)/labda)*(Zm[i,j]+s))/s)*Deltax*Deltay*np.exp((-1j*2*np.pi*R_new)/labda)
        #        else:
        #            ps[i,j] = np.nan
        #
        #diff_inte = np.nansum(ps)
        #phi = np.angle(diff_inte)
        #diff_inte_mag = (diff_inte.real**2+diff_inte.imag**2)**0.5
        ##phi = np.log(diff_inte_new/diff_inte_mag).imag
        #    
        #return diff_inte_mag,phi,ps

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

    def TTL_rec(self,i,t,side='default'):
        if self.simple == True:
            print('Simple mode is on, only calculating for center of receiving telescope')
            TTL = self.phi_gauss(i,t,0,0,side=side)[1]
        else:
            print('Simple mode is off, calculating over whole aperture receiving telescope')
            TTL = self.aperture(self.xlist,self.ylist, lambda dx,dy: self.phi_gauss(i,t,dx,dy,side=side)[1])

        return TTL

    #def TTL(self,i,t,side='default'):
    #    if self.simple == True:
    #        print('Simple mode is on, only calculating for center of receiving telescope')
    #        ttl = self.phi_gauss(i,t,0,0,side=side,print0=True)[2]
    #    else:
    #        print('Simple mode is off, calculating over whole aperture receiving telescope')
    #        ttl = self.aperture(self.xlist,self.ylist, lambda dx,dy: self.phi_gauss(i,t,dx,dy,side=side)[2])

    #    return ttl

    
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

# Jitter
    def jitter_tele(self,t_end,psd_in=False,psd_out=False,psd_r=False,psd_in_ang=False,psd_out_ang=False,fmin=0.0001, fmax=0.1,N=4096):
        jitter_in = self.Ndata.Noise_time(fmin,fmax,N,psd_in,t_end,ret='function')
        jitter_out = self.Ndata.Noise_time(fmin,fmax,N,psd_out,t_end,ret='function')
        jitter_r = self.Ndata.Noise_time(fmin,fmax,N,psd_r,t_end,ret='function')
        jitter_in_ang = self.Ndata.Noise_time(fmin,fmax,N,psd_in_ang,t_end,ret='function')
        jitter_out_ang = self.Ndata.Noise_time(fmin,fmax,N,psd_out_ang,t_end,ret='function')


        self.jitter=[jitter_in,jitter_out,jitter_r,jitter_in_ang,jitter_out_ang]

        return 0
    
    def do_jitter_tele(self,psd_in=False,psd_out=False,psd_r=False,psd_in_ang=False,psd_out_ang=False,t_end='all',fmin=0.0001,fmax=0.1,N=4096):

        if t_end=='all':
            t_end = self.Ndata.t_all[-1]

        jitter_l=[]
        jitter_r=[]
        for i in range(1,4):
            self.jitter_tele(t_end,psd_in=psd_in,psd_out=psd_out,psd_r=psd_r,psd_in_ang=psd_in_ang,psd_out_ang=psd_out_ang,fmin=fmin,fmax=fmax,N=N)
            jitter_l.append(self.jitter)
            self.jitter_tele(t_end,psd_in=psd_in,psd_out=psd_out,psd_r=psd_r,psd_in_ang=psd_in_ang,psd_out_ang=psd_out_ang,fmin=fmin,fmax=fmax,N=N)
            jitter_r.append(self.jitter)
        
        self.jitter_tele_done = [jitter_l,jitter_r]

        # Redo telescope pointing vectors
        self.Ndata.tele_l = lambda i,t: self.Ndata.tele_control(i,t,option='no control',side='l',jitter = self.jitter_tele_done[0])
        self.Ndata.tele_r = lambda i,t: self.Ndata.tele_control(i,t,option='no control',side='r',jitter = self.jitter_tele_done[1])
        self.Ndata.tele_l_fc = lambda i,t: self.Ndata.tele_control(i,t,option='full control',side='l',jitter = self.jitter_tele_done[0])
        self.Ndata.tele_r_fc = lambda i,t: self.Ndata.tele_control(i,t,option='full control',side='r',jitter = self.jitter_tele_done[1])


        return 0


# tele and PAA aim

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
        self.tele_control_func = f
        
        return [[T_all,yout],f]

    def step_res(self,t,tarray,yarray): #... has to be adjusted for realistic respons (created beacuase of mal functioning tele_control_noise
        loc=0
        t_l = tarray[loc]
        t_r = tarray[loc+1]
        y = yarray[loc]
        while not (t<t_r and t>=t_l):
            t_l=tarray[loc]
            t_r=tarray[loc+1]
            y = yarray[loc]
            if loc==len(tarray)-2:
                y = np.nan
                break
            loc = loc+1

        return y


    def tele_control_ss(self,step_max=False,dt=1,simple=False):
        tele_SS_l = []
        tele_SS_r = []
        
        for i in range(1,4):
            if simple==False:
                tele_SS_l.append(self.tele_control_noise(i,side='l',step_max=step_max,dt=dt)[1])
                tele_SS_r.append(self.tele_control_noise(i,side='r',step_max=step_max,dt=dt)[1])

            elif simple==True:
                yarrl = []
                yarrr = []
                for t_step in self.Ndata.data.t_all:
                    yarrl.append(LA.angle(self.Ndata.data.v_l_func_tot(i,t_step),self.Ndata.data.r_func(i,t_step)))
                    yarrr.append(LA.angle(self.Ndata.data.v_r_func_tot(i,t_step),self.Ndata.data.r_func(i,t_step)))
                tele_SS_l.append(lambda t: self.step_res(t,self.Ndata.data.t_all,yarrl))
                tele_SS_r.append(lambda t: self.step_res(t,self.Ndata.data.t_all,yarrr))
        
        self.tele_SS_l = PAA_LISA.utils.func_over_sc(tele_SS_l)
        self.tele_SS_r = PAA_LISA.utils.func_over_sc(tele_SS_r)

        return 0

    def tele_aim(self,step_max=False,dt=1,method=False,simple=False):
        tele_aim_l = []
        tele_aim_r = []

        if method == False:
            method = self.tele_control
        print('The telescope control method is: '+method)
        print(' ')
        for i in range(1,4):
            if method =='SS':
                if simple==False:
                    tele_aim_l.append(self.tele_control_noise(i,side='l',step_max=step_max,dt=dt)[1])
                    tele_aim_r.append(self.tele_control_noise(i,side='r',step_max=step_max,dt=dt)[1])

                elif simple==True:
                    yarrl = []
                    yarrr = []
                    for t_step in self.Ndata.data.t_all:
                        yarrl.append(self.Ndata.tele_l_fc(i,t_step))
                        yarrr.append(self.Ndata.tele_r_fc(i,t_step))
                    tele_aim_l.append(lambda t: self.step_res(t,self.Ndata.data.t_all,yarrl))
                    tele_aim_r.append(lambda t: step_res(t,self.Ndata.data.t_all,yarrr))
                    #... check of this method works

            elif method == 'full control':
                self.tele_aim_l = self.Ndata.tele_l_fc
                self.tele_aim_r = self.Ndata.tele_r_fc
                self.tele_aim_done=True
                
            elif method == 'no control': 
                self.tele_aim_l = self.Ndata.tele_l
                self.tele_aim_r = self.Ndata.tele_r
                self.tele_aim_done=True

        if self.tele_aim_done==False:
            self.tele_aim_l = PAA_LISA.utils.func_over_sc(tele_aim_l)
            self.tele_aim_r = PAA_LISA.utils.func_over_sc(tele_aim_r)
        # Note: length of vectors is arbitrary

        return 0

    def PAAM_control_fc_calc(self,i,t,side='l'):
        
        n = self.Ndata.data.n_func(i,t)
        
        if side=='l':
            opt = self.Ndata.data.v_l_func_tot(i,t)
            tele = self.tele_aim_l(i,t)
        elif side=='r':
            opt = self.Ndata.data.v_r_func_tot(i,t)
            tele = self.tele_aim_r(i,t)
        ang_out_opt = LA.ang_out(opt,n)
        ang_out_tele = LA.ang_out(tele,n)
        PAAM = ang_out_opt - ang_out_tele
        
        return PAAM

    def PAAM_control_vec(self,i,t,side='l'):
        n=self.Ndata.data.n_func(i,t)
        if side=='l':
            ang = self.beam_aim_l(i,t)
            tele = self.tele_aim_l(i,t)

        elif side=='r':
            ang = self.beam_aim_r(i,t)
            tele = self.tele_aim_r(i,t)

        tele_mag = np.linalg.norm(tele)
        tele = LA.unit(tele)
        tele_n = np.dot(tele,n)*n
        tele_r = tele - tele_n
        beam_n_new = np.tan(ang)*np.linalg.norm(tele_r)*n
        beam_new = LA.unit(tele_r+beam_n_new)*tele_mag
            
        return beam_new

    def PAAM_control(self,method=False):
        beam_aim_l = []
        beam_aim_r = []
#        beam_aim_vec_l = []
#        beam_aim_vec_r = []

        if method==False:
            method = self.PAAM_control_method
        print('The PAAM control method is: ' +method)
        print(' ')

        for i in range(1,4):
            if method=='SS':
           
                yarrl = []
                yarrr = []

                for t_step in self.Ndata.data.t_all:
                    yarrl.append(self.PAAM_control_fc_calc(i,t_step,side='l'))
                    yarrr.append(self.PAAM_control_fc_calc(i,t_step,side='r'))

                beam_aim_l.append(lambda t: self.step_res(t,self.Ndata.data.t_all,yarrl))
                beam_aim_r.append(lambda t: self.step_res(t,self.Ndata.data.t_all,yarrr))

            elif method=='fc':
                beam_aim_l.append(lambda t: self.PAAM_control_fc_calc(i,t,side='l'))
                beam_aim_r.append(lambda t: self.PAAM_control_fc_calc(i,t,side='r'))

            elif method=='nc':
                beam_aim_l.append(lambda t:0)
                beam_aim_r.append(lambda t:0)
            
            #beam_aim_vec_l.append(lambda t: self.PAAM_control_vec(i,t,side='l'))
            #beam_aim_vec_r.append(lambda t: self.PAAM_control_vec(i,t,side='r'))


        self.beam_aim_l = PAA_LISA.utils.func_over_sc(beam_aim_l)
        self.beam_aim_r = PAA_LISA.utils.func_over_sc(beam_aim_r)
        self.beam_aim_vec_l = lambda i,t: self.PAAM_control_vec(i,t,side='l')
        self.beam_aim_vec_r = lambda i,t: self.PAAM_control_vec(i,t,side='r')

        return 0



### Sending wavefront

    def Rmn_calc(self,m,n,rho,phi):
        R=0
        m_new = m
        m = abs(m)
        for k in range(0,int((n-m)/2)+1):
            R = R+(((-1**k)*math.factorial(n-k))/(math.factorial(k)*math.factorial(int((n+m)/2)-k)*math.factorial(int((n+m)/2)-k)))*(rho**(n-2*k))
        if m_new<0:
            R = R*np.sin(m*phi)
        else:
            R = R*np.cos(m*phi)
        
        return R

    def Cmn_calc(self,m,n,rho,phi,Cmn,thmn):
        index = str(abs(m))+str(n)
        #index_pos = str(abs(m))+str(n)
        if index in Cmn.keys():
            if n==1:
                if m<0:
                    ret = self.Rmn_calc(m,n,rho,phi)*abs(Cmn[index])*np.cos(thmn[index]) #y
                elif m>=0:
                    ret = self.Rmn_calc(m,n,rho,phi)*abs(Cmn[index])*np.sin(thmn[index]) #x
            elif n==0:
                ret=Cmn[index]
            else:
                ret = self.Rmn_calc(m,n,rho,phi)*abs(Cmn[index])*np.exp(1j*thmn[index]*np.sign(m))
                
        else:
            ret=0

        return ret

    def zern(self,m,n,zmn=False,thmn=False,offset=[0,0],mode='ttl'):
        if zmn==False:
            zmn = self.smn
            thmn = self.thmn
        ps = np.empty((self.Nbinsx,self.Nbinsy),dtype=np.float64)
        for i in range(0,len(self.xlist)):
            for j in range(0,len(self.ylist)):
                x = self.xlist[i]
                y = self.ylist[j]
                rho = (x**2+y**2)**0.5
                if x>0:
                    if y>0:
                        phi=np.arctan(abs(y)/abs(x))
                    elif y<0:
                        phi = -np.arctan(abs(y)/abs(x))
                    elif y==0:
                        phi = 0
                elif x<0:
                    if y>0:
                        phi=np.pi - np.arctan(abs(y)/abs(x))
                    elif y<0:
                        phi=np.pi + np.arctan(abs(y)/abs(x))
                    elif y==0:
                        phi = np.pi
                
                elif x==0:
                    phi = np.pi*0.5*np.sign(y)

                if rho<= self.xlist[-1]:
                    if mode =='ttl':
                        ps[i,j] = self.Cmn_calc(m,n,rho,phi,zmn,thmn)
                    elif mode == 'phase':
                        a = self.Cmn_calc(m,n,rho,phi,zmn,thmn)%self.Ndata.data.labda
                        ps[i,j] = 2*np.pi*(a/self.Ndata.data.labda)
                else:
                    ps[i,j] = np.nan
        return ps

    def tilt_calc(self,i,t,side=False):
        if side==False:
            side = self.side
        #...Check for angle conventions (directions)
        if side=='l':
            ksix = 0#self.tele_SS_l(i,t) # ...replace with misalignment
            ksiy = self.Ndata.data.PAA_func['l_out'](i,t)
        elif side =='r':
            ksix = 0#self.tele_SS_r(i,t) # replace with misallignment
            ksiy = self.Ndata.data.PAA_func['r_out'](i,t)

        if ksix==0:
            thmn = np.pi*0.5
        else:
            thmn = np.arctan(ksix/ksiy)

        zmn = np.linalg.norm((ksix**2+ksiy**2)**0.5)

        return [zmn,thmn]

    def get_tilt(self,side=False):
        if side==False:
            side = self.side
        #... add Z00 with laser noise
        self.zmn_func['11'] = lambda i,t: self.tilt_calc(i,t,side=side)[0]
        self.thmn_func['11'] = lambda i,t: self.tilt_calc(i,t,side=side)[1]


    
    def get_zmn_thmn(self,i,t):
        for keys in self.zmn_func.keys():
            self.zmn[keys] = self.zmn_func[keys](i,t)
            self.thmn[keys] = self.thmn_func[keys](i,t)

    def phasefront_send(self,i,t,nmax=4,side = 'l'):
        self.get_tilt(side=side)
        self.get_zmn_thmn(i,t)
        ps = np.empty((self.Nbinsx,self.Nbinsy),dtype=np.float64)
        ps_list = {}
        for n in range(0,nmax+1):
            for m in range(-n,n+1):
                index = str(m)+str(n)
                ps_new = self.zern(m,n,self.zmn,self.thmn)
                ps_list[index] = ps_new

        keys = ps_list.keys()
        ps_tot = ps_new*0
        for ii in range(0,len(ps_list[keys[0]])):
            for j in range(0,len(ps_list[keys[1]][ii])):
                for k in keys:
                    ps_tot[ii,j] = ps_tot[ii,j]+ps_list[k][ii,j]

        self.ttl_s_tot = ps_tot
        self.ttl_s_dict = ps_list

        return 0


    def do_ttl_send(self):

        self.TTL_tele_send_l = lambda i,t: self.phasefront_send(i,t,side='l')
        self.TTL_tele_send_r = lambda i,t: self.phasefront_send(i,t,side='r')
        
        self.tele_point_SS_l = lambda i,t: LA.unit(LA.rotate(self.Ndata.tele_l(i,t),self.Ndata.data.n_func(i,t),self.tele_SS_l(i,t)))
        self.tele_point_SS_r = lambda i,t: LA.unit(LA.rotate(self.Ndata.tele_r(i,t),self.Ndata.data.n_func(i,t),self.tele_SS_r(i,t)))

        return 0

    def TTL_zern_calc(self,i,t,side='l'):
        [i_self,i_left,i_right] = PAA_LISA.utils.i_slr(i)
        n = self.Ndata.data.n_func(i_self,t)
            
        if side=='l':
            i_next = i_left
            tdel = self.Ndata.data.L_rl_func_tot(i_self,t)
            n_left = self.Ndata.data.n_func(i_left,t-tdel)
            if self.tele_control=='full control':
                tele_rec = self.Ndata.tele_l_fc(i_self,t)
                tele_send = self.Ndata.tele_r_fc(i_left,tdel)
            elif self.tele_control=='no control':
                tele_rec = self.Ndata.tele_l(i_self,t)
                tele_send = self.Ndata.tele_r(i_left,t-tdel)
            elif self.tele_control=='SS':
                r = self.Ndata.data.r_func(i_self,t)
                ang_SS = self.tele_SS_l(i_self,t)
                ang_SS_left = self.tele_SS_r(i_left,t-tdel)
                tele_rec = LA.unit(LA.rotate(self.Ndata.tele_l(i_self,t),n,ang_SS))
                tele_send = LA.unit(LA.rotate(self.Ndata.tele_r(i_left,t-tdel),n_left,ang_SS_left))

            n_beam = self.Ndata.data.n_func(i_left,(t-tdel))
            beam = self.Ndata.PAA_point_r(i_left,t-tdel)#...adjust for sending telescope --> PAAM control
            v_pos = self.Ndata.data.v_r_func_tot(i_left,t-tdel)
            ps_send = self.phasefront_send(i_self,t,side='l')
                
        elif side=='r':
            i_next = i_right
            tdel = self.Ndata.data.L_rr_func_tot(i,t)
            if self.tele_control=='full control':
                tele_rec = self.Ndata.tele_r_fc(i_self,t)
                tele_send = self.Ndata.tele_l_fc(i_right,t-tdel)
            elif self.tele_control=='no control':
                tele_rec = self.Ndata.tele_r(i_self,t)
                tele_send = self.Ndata.tele_l(i_right,t-tdel)
            elif self.tele_control=='SS':
                r = Ndata.data.r_func(i_self,t)
                ang_SS = self.tele_SS_r(i_self,t)
                ang_SS_right = self.tele_SS_l(i_right,t-tdel)

                tele_rec = LA.unit(LA.rotate(self.Ndata.tele_r(i_self,t),n,ang_SS))
                tele_send = LA.unit(LA.rotate(self.Ndata.tele_l(i_right,t-tdel),n_right,ang_SS_right))
            
            n_beam = self.Ndata.data.n_func(i_right,(t-tdel))
            beam = self.Ndata.PAA_point_l(i_right,t-tdel)*scale
            v_pos = self.Ndata.data.v_l_func_tot(i_right,t-tdel)*scale
            ps_send = self.phasefront_send(i_self,t,side='r')

        # Calculating tilt
        [xt,yt,zt] = LA.beam_coor(beam,tele_rec,n)
        angx = np.sin(xt/zt)
        angy = np.sin(yt/zt)
        thmn11 = np.arctan(angx/angy)
        zmn11 = (angx**2 + angy**2)**0.5

        # Calculating offset
        pos_send = np.array(self.Ndata.data.LISA.putp(i_next,t-tdel))
        pos_rec = np.array(self.Ndata.data.LISA.putp(i_self,t-tdel)) #... set on Waluschka)

        tele_rel = LA.unit(tele_rec)*L_tele
        O_tele = pos_rec - pos_send +tele_rel # ... Abram vs. Waluschka
        tele_beam = np.dot(O_tele,LA.unit(beam))*LA.unit(beam)

        tele_yoff = LA.outplane(O_tele - tele_beam,n)
        
        tele_xoff = np.linalg.norm(O_tele-tele_beam-tele_yoff)+djitter[0]
        tele_yoff = np.linalg.norm(tele_yoff)+djitter[1]
        tele_zoff = np.linalg.norm(tele_beam) - np.linalg.norm(beam)+djitter[2]


        zmn={}

        zmn['11'] = zmn11
        thmn={}
        thmn['11'] = thmn11
        zxoff = tele_xoff*np.tan(angx)
        zyoff = tele_yoff*np.tan(angy)
        xoff = tele_xoff
        yoff = tele_yoff
        zoff = np.linalg.norm(tele_beam)+zxoff+zyoff
        zmn['00'] = zoff+dr_rec-dr_send
        thmn['00'] = 0
        xoff = xoff/np.cos(angx)
        yoff = yoff/np.cos(angy)
        offset = np.array([xoff,yoff])

        return zmn, thmn, offset, [angx,angy]

    def TTL_zern(self,i,t,side='l'):
        zmn, thmn, offset, [angx,angy] = self.TTL_zern_calc(i,t,side=side)

        ps=np.zeros((self.Nbinsx,self.Nbinsy))
        for n in range(0,2):
            for m in range(-n,n+1):
                if ((m%2) == (n%2)):
                    ps = ps + self.zern(m,n,zmn=zmn,thmn=thmn,offset=offset)
        wave_ttl = np.nanmean(ps)

        return wave_ttl

    def scan_tele(self):
        self.wf_scan_l = lambda i, t: self.TTL_zern(i,t,side='l')
        self.wf_scan_r = lambda i, t: self.TTL_zern(i,t,side='r')

        return 0








#Obtining TTL by pointing

    def zern_aim(self,i_self,t,side='l'):
        [i_self,i_left,i_right] = PAA_LISA.utils.i_slr(i_self)
        n = self.Ndata.data.n_func(i_self,t)

        if side=='l':
            i_next = i_left
            tdel = self.Ndata.data.L_rl_func_tot(i_self,t)
            n_next = self.Ndata.data.n_func(i_left,t-tdel)
            tele_rec = LA.unit(self.tele_aim_l(i_self,t))*L_tele
            beam_send = self.beam_aim_vec_r(i_next,t-tdel)
            dist = tdel*c
            beam_ideal = self.Ndata.data.v_r_func_tot(i_next,t-tdel)

        elif side=='r':
            i_next = i_right
            tdel = self.Ndata.data.L_rr_func_tot(i_self,t)
            n_next = self.Ndata.data.n_func(i_right,t-tdel)
            tele_rec = LA.unit(self.tele_aim_r(i_self,t))*L_tele
            beam_send = self.beam_aim_vec_l(i_next,t-tdel)
            dist = tdel*c
            beam_ideal = self.Ndata.data.v_l_func_tot(i_next,t-tdel)

        if self.jitter_tele_done!=False:
            if side =='l':
                pr = 0
                ps = 1
            elif side == 'r':
                pr = 1
                ps = 0

            dr_send = self.jitter_tele_done[ps][i_next-1][2](t-tdel)
            din_send = self.jitter_tele_done[ps][i_next-1][0](t-tdel)
            dout_send = self.jitter_tele_done[ps][i_next-1][1](t-tdel)
            dr_rec = self.jitter_tele_done[pr][i_self-1][2](t)
            din_rec = self.jitter_tele_done[pr][i_self-1][0](t)
            dout_rec = self.jitter_tele_done[pr][i_self-1][1](t)
            

            r_self = self.Ndata.data.r_func(i_self,t)
            n_self = n
            x_self = LA.unit(np.cross(n_self,r_self))

            dr_rec = LA.unit(r_self)*dr_rec
            dout_rec = LA.unit(n_self)*dout_rec
            din_rec = x_self*din_rec
            #print(dr_rec,dout_rec,din_rec)

            djitter_rec = dr_rec+dout_rec+din_rec
            #print(djitter_rec)
            djitter_rec = LA.tele_coor(djitter_rec,beam_send,n_next)
            djitter = djitter_rec - np.array([din_send,dout_send,dr_send])

            #print(djitter)

        else:
            dr_send = 0
            din_send = 0
            dout_send = 0
            dr_rec = 0
            din_rec = 0
            dout_rec = 0
            
            djitter=np.array([0,0,0])
        
        # Calculating tilt
        tele_beamcoor = LA.tele_coor(-tele_rec,beam_send,n_next)
        angx = -np.arcsin(tele_beamcoor[0]/np.linalg.norm(tele_beamcoor))
        angy = -np.arcsin(tele_beamcoor[1]/np.linalg.norm(tele_beamcoor))

        #Calculating offset
        tele_pos = beam_ideal - tele_rec
        [xoff,yoff,zoff] = LA.tele_coor(tele_pos,beam_send,n_next)
        piston = self.z_solve(xoff,yoff,zoff)
        R = self.R(piston)
        # Adapting it fot jitter
        xoff=xoff+djitter[0]
        yoff=yoff+djitter[1]
        piston=piston+djitter[2]

        # Tilt by offset
        angxoff = np.arcsin(xoff/R)
        angyoff = np.arcsin(yoff/R)

        # Zernike polynomials
        angx_tot = angx+angxoff
        angy_tot = angy+angyoff
        thmn11 = np.arctan(angx_tot/angy_tot)
        zmn11 = (angx_tot**2 + angy_tot**2)**0.5
        zmn00 = piston

        zmn={}
        thmn={}
        zmn['00'] = zmn00
        zmn['11'] = zmn11
        thmn['11'] = thmn11


        return [xoff,yoff,zoff],zmn,thmn

    def obtain_ttl(self,zmn,thmn,offset,piston=True,tilt=True,mode='ttl'):
        ps=np.zeros((self.Nbinsx,self.Nbinsy),dtype=np.float64)
        if piston==False:
            zmn['00']=False
            thmn['00']=False
        if tilt==False:
            zmn['11']=False
            thmn['11']=False
        for n in range(0,2):
            for m in range(-n,n+1):
                if ((m%2) == (n%2)):
                    ps = ps + self.zern(m,n,zmn=zmn,thmn=thmn,offset=offset,mode=mode)
        if mode=='phase':
            ps = ps%(2*np.pi)
        wave_ttl = np.nanmean(ps)

        return ps, wave_ttl

    def ttl_pointing_function_calc(self,i,t,mode='ttl',side='l',piston=True,tilt=True,ret='value',offset=False):
        [xoff,yoff,zoff],zmn,thmn = self.zern_aim(i,t,side=side)
        if offset==False:
            offset=[0,0] # Because of additional tilt in offset, this can be set to 0 
        elif offset==True:
            offset = [xoff,yoff]
        ps,wave_ttl = self.obtain_ttl(zmn,thmn,offset,piston=piston,tilt=tilt,mode=mode)
        
        if ret=='value':
            return wave_ttl

        elif ret=='aperture':
            return ps


    def ttl_pointing_function(self,mode='ttl',ret='value',option='all'):
        
        if option=='piston':
            self.ttl_l = lambda i,t: self.ttl_pointing_function_calc(i,t,mode=mode,side='l',piston=True,tilt=False,ret=ret)
            self.ttl_r = lambda i,t: self.ttl_pointing_function_calc(i,t,mode=mode,side='r',piston=True,tilt=False,ret=ret)
        elif option=='tilt':
            self.ttl_l = lambda i,t: self.ttl_pointing_function_calc(i,t,mode=mode,side='l',piston=False,tilt=True,ret=ret)
            self.ttl_r = lambda i,t: self.ttl_pointing_function_calc(i,t,mode=mode,side='l',piston=False,tilt=True,ret=ret)
        elif option=='all':
            self.ttl_l = lambda i,t: self.ttl_pointing_function_calc(i,t,mode=mode,side='l',piston=True,tilt=True,ret=ret)
            self.ttl_r = lambda i,t: self.ttl_pointing_function_calc(i,t,mode=mode,side='l',piston=True,tilt=True,ret=ret)
        

















