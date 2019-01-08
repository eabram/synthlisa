from imports import *
from functions import *
from parameters import *
import PAA_LISA
import NOISE_LISA

class WFE():
    def __init__(self,**kwargs):
        global labda, w0, E0, k, D, LA,L_tele
        
        
        self.Ndata = kwargs.pop('Ndata',False)
        self.tele_control = kwargs.pop('tele_control','no control')
        self.PAAM_control_method = kwargs.pop('PAAM_control','SS')
        self.side = kwargs.pop('side','l')
        self.speed_on = kwargs.pop('speed_on',0)
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

    def get_pointing(self,tele_method = False,PAAM_method=False): #...add more variables
        aim = NOISE_LISA.AIM(self)
        aim.tele_aim(method=tele_method)
        aim.PAAM_control(method=PAAM_method)

        self.aim = aim


    # Beam properties equations

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


# WFE send
    def WFE_send(self,i,t,side='l',xlist=False,ylist=False): #...adjust wirt real WFE instead of 0
        if self.speed_on==2:
            xlist=np.array([0])
            ylist=np.array([0])
        else:
            if xlist==False:
                xlist = self.xlist
            if ylist == False:
                ylist =self.ylist

        if side=='l':
            angy = self.PAAM_aim_l_ang(i,t)
        elif side=='r':
            angy = self.PAAM_aim_r_ang(i,t)
        angx = 0 #...to do: Add jitter
       
        labda = self.Ndata.data.labda
        function = lambda x,y: 2*np.pi*((x*np.sin(angx)+y*np.sin(angy))/labda)

        w = self.aperture(xlist,ylist,function,dType = np.float64)

        return w

    def w0(self,i,t,side,ksi):
        if side=='l':
            angy = self.PAAM_aim_l_ang(i,t)
        elif side == 'r':
            angy = self.PAAM_aim_r_ang(i,t)
        angx=0#...add jitter

        [x,y] = ksi
        
        labda = self.Ndata.data.labda
        w_error = 2*np.pi*((x*np.sin(angx)+y*np.sin(angy))/labda)

        return w_error

    def u0(self,ksi):#...for gaussian beam
        w = self.w(0)
        [x,y] = ksi

        return np.exp(-((x**2+y**2)/(w**2)))

# WFE receive
    
    def u_rz_calc(self,r,z,SC,t,side,xlist=False,ylist=False):
        if xlist==False:
            xlist = self.xlist
        if ylist==False:
            ylist = self.ylist
        labda = self.Ndata.data.labda
        k = (2*np.pi)/labda
        
        dksi = (xlist[1]-xlist[0])*(ylist[1]-ylist[0])
        ret=0
        for i in range(0,len(xlist)):
            for j in range(0,len(ylist)):
                ksi = np.array([xlist[i],ylist[j]])
                T1 = np.exp((1j*k*np.dot(r,ksi))/z)
                T2 = self.u0(ksi)
                T3 = np.exp(1j*self.w0(SC,t,side,ksi))

                ret = ret+T1*T2*T3
        ret = ret*dksi*(1j*k*np.exp(-(1j*k*(np.linalg.norm(r)**2))/(2*z))/(2*np.pi*z))

        return ret

    def u_rz(self,zmn,thmn,ksi,i,t,side='l',xlist=False,ylist=False):
        [x0,y0] = ksi
        z0 = zmn['00']
        angx = zmn['11']*np.cos(thmn['11'])
        angy = zmn['11']*np.sin(thmn['11'])
        x = x0*np.cos(angx)
        y = y0*np.cos(angy)
        z = z0 + x0*np.sin(angx)+y0*np.sin(angy)
        r = np.array([x,y])

        u = self.u_rz_calc(r,z,i,t,side,xlist=xlist,ylist=ylist)
        w = self.w(z) #...check if proper z is used (z0)
        wac = (self.Ndata.data.D/2)/w
        norm = ((np.pi*((2*np.pi)/self.Ndata.data.labda)*w**2*(1-np.exp(-1/(wac**2))))/(2*np.pi*z))**2
        I = (abs(u)**2)/norm
        phase = np.angle(u)

        return u,I,phase

    def u_rz_aperture(self,zmn,thmn,i,t,side='l',xlist=False,ylist=False,mode='power'):
        if self.speed_on>=1: #Only canculates center (works correctly for only piston and tilt
            xlist = np.array([0])
            ylist = np.array([0])
         
        if mode=='power':
            function = lambda x,y: self.u_rz(zmn,thmn,np.array([x,y]),i,t,side=side,xlist=xlist,ylist=ylist)[1]
            ps = self.aperture(xlist,ylist,function,dType=np.float64)
        elif mode=='u':
            function = lambda x,y: abs(self.u_rz(zmn,thmn,np.array([x,y]),i,t,side=side,xlist=xlist,ylist=ylist)[0])
            ps = self.aperture(xlist,ylist,function,dType=np.float64)
        elif mode=='phase':
            function = lambda x,y: self.u_rz(zmn,thmn,np.array([x,y]),i,t,side=side,xlist=xlist,ylist=ylist)[2]
            ps = self.aperture(xlist,ylist,function,dType=np.float64)

        
        
        return ps

    


    def aperture(self,xlist,ylist,function,dType=np.complex64): # Creates matrix of function over an aperture (circle)
        #print('Type of telescope control is: ' + self.tele_control)
        if type(xlist)==bool:
            if xlist==False:
                xlist = self.xlist
        if type(ylist)==bool:
            if ylist==False:
                ylist = self.ylist

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
    
    def tele_control_ang_fc_calc(self,i,t,side='l'):
        coor = self.coor_SC(i,t)
        if side=='l':
            v = -self.Ndata.data.u_l_func_tot(i,t)
        elif side=='r':
            v = -self.Ndata.data.u_r_func_tot(i,t)

        v_SC = LA.matmul(coor,v)
        
        ang = np.arctan(v_SC[2]/v_SC[0])

        return ang

    def tele_control_ang_fc(self):

        self.tele_ang_l_fc = lambda i,t: self.tele_control_ang_fc_calc(i,t,side='l')
        self.tele_ang_r_fc = lambda i,t: self.tele_control_ang_fc_calc(i,t,side='r')
        
        return 0


    def tele_aim(self,method=False,dt=3600*24,jitter=False):
        self.tele_control_ang_fc()

        if method == False:
            method = self.tele_control

        print('The telescope control method is: '+method)
        print(' ')
        if method=='full control':
            tele_aim_l = self.tele_ang_l_fc
            tele_aim_r = self.tele_ang_r_fc

        elif method=='no control':
            tele_aim_l = lambda i,t: np.radians(-30)
            tele_aim_r = lambda i,t: np.radians(30)

        elif method=='SS':
            tele_aim_l = lambda i,t: self.tele_ang_l_fc(i,t-(t%dt))
            tele_aim_r = lambda i,t: self.tele_ang_r_fc(i,t-(t%dt))

        self.tele_aim_l_0 = tele_aim_l
        self.tele_aim_r_0 = tele_aim_r
        if jitter!=False:
            self.tele_aim_l = lambda i,t: self.add_jitter(self.tele_aim_l_0,i,t,1e-6,1e10,dt=0.1)
            self.tele_aim_r = lambda i,t: self.add_jitter(self.tele_aim_r_0,i,t,1e-6,1e10,dt=0.1)
        else:
            self.tele_aim_l = self.tele_aim_l_0
            self.tele_aim_r = self.tele_aim_r_0

        self.tele_aim_l_vec = lambda i,t: LA.unit(self.coor_tele(i,t,self.tele_aim_l(i,t))[0])*L_tele
        self.tele_aim_r_vec = lambda i,t: LA.unit(self.coor_tele(i,t,self.tele_aim_r(i,t))[0])*L_tele
        
        self.tele_aim_l_coor = lambda i,t: self.coor_tele(i,t,self.tele_aim_l(i,t))
        self.tele_aim_r_coor = lambda i,t: self.coor_tele(i,t,self.tele_aim_r(i,t))

        return 0

    def add_jitter(self,ang_func,i,t,dang,scale_v,dt=0.1):
        # add position jitter
        # add velocity jitter
        v = (ang_func(i,t) - ang_func(i,t-dt))/dt
    
        return np.random.normal(ang_func(i,t),dang*(1+v*scale_v))#...adjust: make correlated errors
        

    
    def PAAM_control(self,method=False,dt=3600*24,jitter=False):
        if method==False:
            method = self.PAAM_control_method
        print('The PAAM control method is: ' +method)
        print(' ')
        
        ang_fc_l = lambda i,t: self.Ndata.data.PAA_func['l_out'](i,t)
        ang_fc_r = lambda i,t: self.Ndata.data.PAA_func['r_out'](i,t)

        if method=='fc':
            ang_l = ang_fc_l
            ang_r = ang_fc_r
        elif method=='nc':
            ang_l = lambda i,t: 0
            ang_r = lambda i,t: 0
        elif method=='SS':
            ang_l = lambda i,t: ang_fc_l(i,t-(t%dt))
            ang_r = lambda i,t: ang_fc_r(i,t-(t%dt))

        self.PAAM_aim_l_ang_0 = ang_l
        self.PAAM_aim_r_ang_0 = ang_r
        
        if jitter!=False:
            self.PAAM_aim_l_ang = lambda i,t: self.add_jitter(self.PAAM_aim_l_ang_0,i,t,1e-8,1e20,dt=3600)
            self.PAAM_aim_r_ang = lambda i,t: self.add_jitter(self.PAAM_aim_r_ang_0,i,t,1e-8,1e20,dt=3600)
        else:
            self.PAAM_aim_l_ang = self.PAAM_aim_l_ang_0
            self.PAAM_aim_r_ang = self.PAAM_aim_r_ang_0

        self.beam_aim_l_coor = lambda i,t: self.beam_tele(i,t,self.tele_aim_l(i,t),self.PAAM_aim_l_ang(i,t))
        self.beam_aim_r_coor = lambda i,t: self.beam_tele(i,t,self.tele_aim_r(i,t),self.PAAM_aim_r_ang(i,t))

        self.beam_aim_l_vec = lambda i,t: self.beam_aim_l_coor(i,t)[0]*self.Ndata.data.L_rl_func_tot(i,t)*c
        self.beam_aim_r_vec = lambda i,t: self.beam_aim_r_coor(i,t)[0]*self.Ndata.data.L_rr_func_tot(i,t)*c
        
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
            zmn = self.zmn
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

    #def TTL_zern_calc(self,i,t,side='l'):
    #    [i_self,i_left,i_right] = PAA_LISA.utils.i_slr(i)
    #        
    #    if side=='l':
    #        i_next = i_left
    #        tdel = self.Ndata.data.L_rl_func_tot(i_self,t)
    #        beam_send = self.beam_aim_r_vec(i_next,t-tdel)
    #        tele_rec_coor = self.tele_aim_l_coor(i_self,t)
    #    
    #    if side=='r':
    #        i_next = i_right
    #        tdel = self.Ndata.data.L_rr_func_tot(i_self,t)
    #        beam_send = self.beam_aim_l_vec(i_next,t-tdel)
    #        tele_rec_coor = self.tele_aim_r_coor(i_self,t)
    #    
    #    [piston,dy,dx] = LA.matmul(-tele_rec_coor,beam_send)
    #    angx = np.arctan(dx/piston)
    #    angy = np.arctan(dy/piston)
    #     
    #    # Calculating tilt
    #    thmn11 = np.arctan(angx/angy)
    #    zmn11 = (angx**2 + angy**2)**0.5

    #    # Calculating offset
    #    pos_send = np.array(self.Ndata.data.LISA.putp(i_next,t-tdel))
    #    pos_rec = np.array(self.Ndata.data.LISA.putp(i_self,t-tdel)) #... set on Waluschka)

    #    tele_rel = LA.unit(tele_rec)*L_tele
    #    O_tele = pos_rec - pos_send +tele_rel # ... Abram vs. Waluschka
    #    tele_beam = np.dot(O_tele,LA.unit(beam))*LA.unit(beam)

    #    tele_yoff = LA.outplane(O_tele - tele_beam,n)
    #    
    #    tele_xoff = np.linalg.norm(O_tele-tele_beam-tele_yoff)+djitter[0]
    #    tele_yoff = np.linalg.norm(tele_yoff)+djitter[1]
    #    tele_zoff = np.linalg.norm(tele_beam) - np.linalg.norm(beam)+djitter[2]


    #    zmn={}

    #    zmn['11'] = zmn11
    #    thmn={}
    #    thmn['11'] = thmn11
    #    zxoff = tele_xoff*np.tan(angx)
    #    zyoff = tele_yoff*np.tan(angy)
    #    xoff = tele_xoff
    #    yoff = tele_yoff
    #    zoff = np.linalg.norm(tele_beam)+zxoff+zyoff
    #    zmn['00'] = zoff+dr_rec-dr_send
    #    thmn['00'] = 0
    #    xoff = xoff/np.cos(angx)
    #    yoff = yoff/np.cos(angy)
    #    offset = np.array([xoff,yoff])

    #    return zmn, thmn, offset, [angx,angy]

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

        if side=='l':
            i_next = i_left
            tdel = self.Ndata.data.L_rl_func_tot(i_self,t)
            beam_send = self.beam_aim_r_vec(i_next,t-tdel)
            beam_send_coor = self.beam_aim_r_coor(i_next,t-tdel)
            tele_rec_coor = self.tele_aim_l_coor(i_self,t)

        if side=='r':
            i_next = i_right
            tdel = self.Ndata.data.L_rr_func_tot(i_self,t)
            beam_send = self.beam_aim_l_vec(i_next,t-tdel)
            beam_send_coor = self.beam_aim_l_coor(i_next,t-tdel)
            tele_rec_coor = self.tele_aim_r_coor(i_self,t)

        tele_rec = L_tele*tele_rec_coor[0]
        [piston,dy,dx] = LA.matmul(-tele_rec_coor,beam_send)
        
        # Calculating tilt
        angx = np.arctan(dx/piston)
        angy = np.arctan(dy/piston)
        
        # Calculating offset
        pos = beam_send - tele_rec # position of aperture
        [dz,dy,dx] = LA.matmul(beam_send_coor,pos)



        if self.jitter_tele_done!=False: #..adjust for new method
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
        
        
        # Calculating offset
        pos = beam_send - tele_rec # position of aperture
        [zoff,yoff,xoff] = LA.matmul(beam_send_coor,pos)
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
        self.zmn = zmn
        self.thmn = thmn

        return [xoff,yoff,zoff],zmn,thmn

# Calulate P
    def P_calc(self,i,t,side='l'): # Only simple calculation (middel pixel)
        if side=='l':
            P = np.cos(self.zern_aim(i,t,side='l')[1]['11'])*self.Ndata.data.P_L
        elif side=='r':
            P = np.cos(self.zern_aim(i,t,side='r')[1]['11'])*self.Ndata.data.P_L
        
        return P

    def get_zern_poly(self,zmn=False,thmn=False,nmax=4):
        if zmn==False:
            zmn = self.zmn
        if thmn==False:
            thmn = self.thmn

        for n in range(0,nmax+1):
            for m in range(-n,n+1):
                if (m%2)==(n%2):
                    key=str(m)+str(n)
                    if key not in zmn.keys():
                        zmn[key]=0
                    if key not in thmn.keys():
                        thmn[key]=0

        self.zmn=zmn
        self.thmn=thmn

        return zmn,thmn


    def zern_para(self,z=False,zmn=False,thmn=False,x=0,y=0):
        if zmn==False:
            zmn = self.zmn
        if thmn==False:
            thmn = self.thmn

        x0 = x
        y0 = y
        angx = zmn['11']*np.cos(thmn['11'])
        angy = zmn['11']*np.sin(thmn['11'])

        x = x0*np.cos(angx)
        y = y0*np.cos(angy)
        if z==False:
            z = zmn['00']
        z = z + x0*np.sin(angx)+y0*np.sin(angy)

        labda = self.Ndata.data.labda
        w = self.w(z)
        r0 = self.Ndata.data.D*0.5
        k = (2*np.pi)/labda #...adjust for laser noise
        wac = r0/w
        q = -1.0/(wac**2) # ...In Sasso paper not with minus sign!!!
        print(wac,q) 
        zmn,thmn = self.get_zern_poly(zmn=zmn,thmn=thmn)

        z02 = zmn['02']
        z04 = zmn['04']
        z22abs = zmn['22']
        z33abs = zmn['33']
        th33 = thmn['33']
        th22 = thmn['22']
        th11 = thmn['11']
        z13abs = zmn['13']
        th13 = thmn['13']

        dzx = zmn['11']*np.cos(thmn['11'])#... adjust to (36a)
        dzy = zmn['11']*np.sin(thmn['11'])#...ajust to (36b)


        A2 = (1+np.exp(q)+2*(1-np.exp(q))*(wac**2))/(1-np.exp(q))
        A4 = (1-np.exp(q)+6*(1+np.exp(q))*(wac**2) + 12*(1-np.exp(q))*(wac**4))/(1-np.exp(q))
        B = (-2*(1+3*(wac**2)+6*(wac**4)+6*(1-np.exp(q))*(wac**6)))/(1-np.exp(q))
        C = (-2*(1+5*(wac**2)+2*(7+2*np.exp(q))*(wac**4)+18*(1-np.exp(q))*(wac**6)))/(1-np.exp(q))
        D = (4*(np.exp(q)+6*np.exp(q)*(wac**2)-2*(2-np.exp(q) - np.exp(2*q))*(wac**4) - 12*((1-np.exp(q))**2)*(wac**6)))/((1-np.exp(q))**2)
        G = (24*((np.exp(q)*(wac**2))-(2-9*np.exp(q)+np.exp(2*q))*(wac**4) -2*(7-2*np.exp(q)-5*np.exp(2*q))*(wac**6) -30*(1-np.exp(q))*(wac**8)))/((1-np.exp(q))**2)

        #G = 0#3(24*(np.exp(q)*(wac**2)-(2-9*np.exp(q) + np.exp(2*q))*(wac**4) - 2*(7-2*np.exp(q) - 5*np.exp(2*q))*(wac**6) -30((1-np.exp(q))**2)*(wac**8)))/((1-np.exp(q))**2)
        E = (2*(np.exp(q) - ((1-np.exp(2*q))**2)*(wac**4)))/((1-np.exp(q))**2)
        F = (-1*(1+2*(wac**2)+2*(1-np.exp(q))*(wac**4)))/(1-np.exp(q))
        H = (6*(2*np.exp(q)*(wac**2) - (1-np.exp(2*q))*(wac**4)-4*((1-np.exp(q))**2)*(wac**6)))/((1-np.exp(q))**2)

        b0 = A2*z02+A4*z04
        b1 = B*z22abs*z33abs*np.cos(th33-th22-th11)+C*z22abs*z13abs*np.cos(th22-th13-th11)+D*z02*z13abs*np.cos(th13-th11)+G*z04*z13abs*np.cos(th13-th11)
        b2 = E*z02+F*z22abs*np.cos(th22-2*th11)+H*z04

        b00 = b0
        b10 = B*np.cos(th33-th22)*z33abs*z22abs +C*np.cos(th22-th13)*z13abs*z22abs+D*np.cos(th13)*z13abs*z02 +G*np.cos(th13)*z13abs*z22abs
        b01 = B*np.sin(th33-th22)*z33abs*z22abs +C*np.sin(th22-th13)*z13abs*z22abs+D*np.sin(th13)*z13abs*z02 +G*np.sin(th13)*z13abs*z22abs
        b20 = E*z02+F*np.cos(th22)*z22abs+H*z04
        b02 = E*z02-F*np.cos(th22)*z22abs+H*z04
        b11 = 2*F*np.sin(th22)*z22abs

        
        # Power density
        c1 = -(4*(1+2*(2+np.exp(q))*(wac**2)+6*(1-np.exp(q))*(wac**4)))/(1-np.exp(-q))
        c2 = (2*(1+(1-np.exp(q))*(wac**2)))/(1-np.exp(-q))
        I = 1+c1*z13abs*(np.cos(th13)*dzx+np.sin(th13)*dzy)-c2*(dzx**2+dzy**2)
        print(c1,c2,I)
       
        u0,[a0,a1,a2,a3] = self.WFE_rec(zmn,thmn,xlist=False,ylist=False,fast_cal=True)

        u_2 = ((np.pi**2)*(w**4)*abs(a0+1j*a1+a2)**2*(k**2))/((2*np.pi*z)**2)
        I_2 = u_2/(((np.pi*k*w**2)*(1-np.exp(-1.0/(wac**2)))/(2*np.pi*z))**2)




        return u_2,I,I_2

#    def power_rec(self,y=0,x=0):
#        u0,I0 = self.zern_para(z=False,zmn=False,thmn=False)
#        
#        x0 = x
#        y  
#        x = 










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
        

















