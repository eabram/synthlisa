from imports import *
from functions import *
from parameters import *
import PAA_LISA
import NOISE_LISA

# tele and PAA aim
class AIM():
    
    def __init__(self,wfe,**kwargs):
        print('Start calculating telescope and PAAM aim')
        self.wfe = wfe
        self.PAAM_method = wfe.PAAM_control_method
        self.tele_method = wfe.tele_control
        self.offset_control = kwargs.pop('offset_control',True)
        global LA
        LA = PAA_LISA.la()

    def static_tele_angle(self,select,i,dt=False,side='l'):
        if select=='PAAM':
            if side=='l':
                func = self.wfe.Ndata.data.PAA_func['l_out']
                #func_y = self.wfe.Ndata.data.PAA_func['l_out']
            elif side=='r':
                func = self.wfe.Ndata.data.PAA_func['r_out']
                #func_y = self.wfe.Ndata.data.PAA_func['r_out']
        
        elif select=='tele':
            if side=='l':                                
                func = self.tele_ang_l_fc    
            elif side=='r':
                func = self.tele_ang_r_fc 

        t_all = self.wfe.Ndata.data.t_all
        if dt==False:
            dt = t_all[1]-t_all[0]
        t_vec = np.linspace(t_all[0],t_all[1],(t_all[1]-t_all[0])/dt)
        val=[]
        for t in t_vec:
            val.append(func(i,t))

        return np.mean(val)

    def do_static_tele_angle(self,select,dt=False):
        tele_ang_off_l = lambda i: self.static_tele_angle(select,i,dt=dt,side='l')
        tele_ang_off_r = lambda i: self.static_tele_angle(select,i,dt=dt,side='r')
        
        if select=='PAAM':
            self.offset_PAAM_l = tele_ang_off_l
            self.offset_PAAM_r = tele_ang_off_r
        elif select=='tele':
            self.offset_tele_l = tele_ang_off_l
            self.offset_tele_r = tele_ang_off_r

        return 0

    def tele_control_ang_fc_calc(self,i,t,side='l'):
        coor = NOISE_LISA.coor_SC(self.wfe,i,t)
        if side=='l':
            v = -self.wfe.Ndata.data.u_l_func_tot(i,t)
        elif side=='r':
            v = -self.wfe.Ndata.data.u_r_func_tot(i,t)

        v_SC = LA.matmul(coor,v)
        
        #print(v_SC)
        ang = np.arcsin(v_SC[2]/np.linalg.norm(v_SC[0])) # Angle betweed x and r component, whih is the optimal telescope pointing angle (inplane)

        return ang

    def tele_control_ang_fc(self):
        # Obtaines functions for optimal telescope pointing vector
        self.tele_ang_l_fc = lambda i,t: self.tele_control_ang_fc_calc(i,t,side='l')
        self.tele_ang_r_fc = lambda i,t: self.tele_control_ang_fc_calc(i,t,side='r')

        return 0


    def tele_aim(self,method=False,dt=3600*24*10,jitter=False,tau=3600*24*5,mode='overdamped'):
        self.tele_control_ang_fc()

        if method == False:
            method = self.tele_method

        print('The telescope control method is: '+method)
        print(' ')

        # Calculating teescope angle for 'full control', 'no control' and 'SS' (step and stair)
        if method=='full control':
            tele_l = self.tele_ang_l_fc
            tele_r = self.tele_ang_r_fc

        elif method=='no control':
            self.do_static_tele_angle('tele')
            if self.offset_control==True:
                tele_l = lambda i,t: self.offset_tele_l(i)
                tele_r = lambda i,t: self.offset_tele_r(i)
            elif self.offset_control==False:
                tele_l = lambda i,t: np.radians(-30)
                tele_r = lambda i,t: np.radians(30)


        elif method=='SS': #After dt, the telescope is pointed again
            tele_l_SS = lambda i,t: self.tele_ang_l_fc(i,t-(t%dt))
            tele_r_SS = lambda i,t: self.tele_ang_r_fc(i,t-(t%dt))
            print('Taken '+mode+' step response for telescope SS control with tau='+str(tau)+'sec')
            tele_l = self.step_response(tele_l_SS,dt,tau=tau,mode=mode)
            tele_r = self.step_response(tele_r_SS,dt,tau=tau,mode=mode)
            self.tele_l_ang_SS = tele_l_SS
            self.tele_r_ang_SS = tele_r_SS

        else:
            raise ValueError('Please select a valid telescope pointing method')

        # Adding jitter
        if jitter!=False:
            self.tele_l_ang = lambda i,t: self.add_jitter(tele_l,i,t,1e-6,1e10,dt=0.1)
            self.tele_r_ang = lambda i,t: self.add_jitter(tele_r,i,t,1e-6,1e10,dt=0.1)
        else:
            self.tele_l_ang = tele_l
            self.tele_r_ang = tele_r

        # Calculating new pointing vectors and coordinate system
        self.tele_l_vec = lambda i,t: LA.unit(NOISE_LISA.coor_tele(self.wfe,i,t,self.tele_l_ang(i,t))[0])*L_tele
        self.tele_r_vec = lambda i,t: LA.unit(NOISE_LISA.coor_tele(self.wfe,i,t,self.tele_r_ang(i,t))[0])*L_tele

        self.tele_l_coor = lambda i,t: NOISE_LISA.coor_tele(self.wfe,i,t,self.tele_l_ang(i,t))
        self.tele_r_coor = lambda i,t: NOISE_LISA.coor_tele(self.wfe,i,t,self.tele_r_ang(i,t))

        return 0



    def add_jitter(self,ang_func,i,t,dang,scale_v,dt=0.1):
        # add position jitter
        # add velocity jitter
        v = (ang_func(i,t) - ang_func(i,t-dt))/dt

        return np.random.normal(ang_func(i,t),dang*(1+v*scale_v))#...adjust: make correlated errors


    def SS_control(self,function,i,t,dt):
        t0 = t-(t%dt)
        if t0==0 or t0==self.wfe.Ndata.data.t_all[-1]:
            ret = np.nan
        else:
            ret = function(i,t-(t%dt))

        return ret

    def PAAM_control(self,method=False,dt=3600*24,jitter=False,tau=3600*12,mode='overdamped'):
        if method==False:
            method = self.PAAM_method
        print('The PAAM control method is: ' +method)
        print(' ')

        ang_fc_l = lambda i,t: self.wfe.Ndata.data.PAA_func['l_out'](i,t)
        ang_fc_r = lambda i,t: self.wfe.Ndata.data.PAA_func['r_out'](i,t)

        # Obtaining PAAM angles for 'fc' (full control), 'nc' (no control) and 'SS' (step and stair)
        
        if method=='full control':
            ang_l = ang_fc_l
            ang_r = ang_fc_r
        elif method=='no control':
            self.do_static_tele_angle('PAAM')
            if self.offset_control==True:
                ang_l = lambda i,t: self.offset_PAAM_l(i)
                ang_r = lambda i,t: self.offset_PAAM_r(i)

            elif self.offset_control==False:
                ang_l = lambda i,t: 0
                ang_r = lambda i,t: 0
        elif method=='SS':
            ang_l_SS = lambda i,t: self.SS_control(ang_fc_l,i,t,dt) # Adjusting the pointing every dt seconds
            ang_r_SS = lambda i,t: self.SS_control(ang_fc_r,i,t,dt)
            print('Taken '+method+' step response for PAAM SS control with tau='+str(tau)+'sec')
            ang_l = self.step_response(ang_l_SS,dt,tau=tau,mode=mode)
            ang_r = self.step_response(ang_r_SS,dt,tau=tau,mode=mode)
            f_noise_l = lambda i,t: (ang_l(i,t)-ang_l_SS(i,t))**2
            f_noise_r = lambda i,t: (ang_r(i,t)-ang_r_SS(i,t))**2
            self.PAAM_ang_l_SS = ang_l_SS
            self.PAAM_ang_r_SS = ang_r_SS
        else:
            raise ValueError('Please select a valid PAAM pointing method')

        # Adding jitter
        if jitter!=False:
            self.PAAM_l_ang = lambda i,t: self.add_jitter(ang_l,i,t,1e-8,1e20,dt=3600)
            self.PAAM_r_ang = lambda i,t: self.add_jitter(ang_r,i,t,1e-8,1e20,dt=3600)
        else:
            self.PAAM_l_ang = ang_l
            self.PAAM_r_ang = ang_r

        # Calculating new pointing vectors and coordinate system
        self.beam_l_coor = lambda i,t: NOISE_LISA.beam_tele(self.wfe,i,t,self.tele_l_ang(i,t),self.PAAM_l_ang(i,t))
        self.beam_r_coor = lambda i,t: NOISE_LISA.beam_tele(self.wfe,i,t,self.tele_r_ang(i,t),self.PAAM_r_ang(i,t))

        self.beam_l_vec = lambda i,t: self.beam_l_coor(i,t)[0]*self.wfe.Ndata.data.L_rl_func_tot(i,t)*c
        self.beam_r_vec = lambda i,t: self.beam_r_coor(i,t)[0]*self.wfe.Ndata.data.L_rr_func_tot(i,t)*c

        return 0

    
    def step_response_calc(self,function,i,t,dt,tau,mode='overdamped'):
        if mode=='overdamped':
            t0 = t-(t%dt)
            t1 = t0+dt
            Y0 = function(i,t0)
            Y1 = function(i,t1)
            if t0==0:
                Y0=Y1
            return Y1+(Y0-Y1)*np.exp(-(t-t0)/tau)
        elif mode==False:
            return function(i,t)

    def step_response(self,function,dt,tau=3600,mode='overdamped'):
        return lambda i,t: self.step_response_calc(function,i,t,dt,tau=tau,mode=mode)


