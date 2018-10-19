from synthlisa import *
import numpy as np
import matplotlib.pyplot as plt
import os
from fractions import Fraction
import math
import datetime
from scipy.interpolate import interp1d
from scipy.interpolate import RegularGridInterpolator
from class_orbit import orbit
import warnings
warnings.filterwarnings('ignore', 'The iteration is not making good progress')
#warnings.filterwarnings("error")
import scipy.optimize

plot_on=True
year2sec=32536000
day2sec=year2sec/365.25
c=300000000
noise_check=False
#noise_check=True

#home='/home/ester/Dropbox/Master_Project/Synthetic_LISA/'
#filename='NGO_1M_10deg_synthlisa.txt'
#directory_imp='LISAorbits/NGO_1M_10deg/'
home='/home/ester/git/synthlisa/'
directory_imp='lisasim/data/'
filename='positions.txt'
#filename='Folkner_orbit.txt'

num_back=0

class la():
    def norm(self,v):
        return np.linalg.norm(v)

    def unit(self,v):
        return v/self.norm(v)

    def angle(self,v1,v2,dot=False):
        norm_v1 = self.norm(v1)
        norm_v2 = self.norm(v2)
        if norm_v1!=0 and norm_v2!=0:
            if dot==False:
                sin = self.norm(np.cross(v1,v2)/(norm_v1*norm_v2))
                return np.arcsin(sin)
            elif dot == True:
                cos = np.dot(v1,v2)/(norm_v1*norm_v2)
                return np.sign(np.dot(v1,v2))*np.arccos(cos)
        else:
            #print('norm v1: '+str(norm_v1))
            #print('norm v2: '+str(norm_v2))

            return np.nan

    def inplane(self,v,n):
        inplane_calc = v - (np.dot(v,n)/(np.linalg.norm(n)**2))*n
        return inplane_calc

    def outplane(self,v,n):
        outplane_calc = (np.dot(v,n)/(np.linalg.norm(n)**2))*n
        return outplane_calc

    def ang_out(self,v,n):
        sign = np.sign(np.dot(self.outplane(v,n),n))
        return sign*self.angle(self.inplane(v,n),v)

    def ang_in(self,v,n,r):
        #ang_out_calc = self.ang_out(v,n)
        inplane_calc = self.inplane(v,n)
        ang_in_calc = self.angle(inplane_calc,r)

        return ang_in_calc

    def ang_in_dot(self,v,v_stat,n,r):
        inplane_calc = self.inplane(v,n)
        costheta_beam = np.dot(inplane_calc,r)/(self.norm(inplane_calc)*self.norm(r))
        inplane_calc = self.inplane(v_stat,n)
        costheta_sc = np.dot(inplane_calc,r)/(self.norm(inplane_calc)*self.norm(r))

        return np.abs(np.arccos(costheta_beam)) - np.abs(np.arccos(costheta_sc))


    def ang_in_direct(self,v,v_stat,n,r):
        inplane_calc = self.inplane(v,n)
        inplane_stat = self.inplane(v_stat,n)
        ang_out_calc = self.angle(inplane_calc,inplane_stat)
        #ang1 = self.angle(inplane_calc,r)
        #ang2 = self.angle(inplane_stat,r)
        #sign = 1#np.sign(ang1 - ang2)

        return ang_out_calc#ang1-ang2

    def print_component(self,v,v_in,v_out,v_arm):
        n = self.norm(v)
        n_in = self.norm(v_in)
        n_out = self.norm(v_out)
        n_arm = self.norm(v_arm)

        print(n_in/n)
        print((n_out**2+n_in**2+n_arm**2)/n**2)
        print('')

        return 0

    def ang_in_out(self,v1,v2,n,v_stat,r):
        n = self.unit(n)
        v1_out = (np.dot(v1,n)*n)/(self.norm(n)**2)
        v1_arm = (np.dot(v1,v_stat)*v_stat)/(self.norm(v_stat)**2)
        v1_in = v1 - v1_out - v1_arm

        v2_out = (np.dot(v2,n)*n)/(self.norm(n)**2)
        v2_arm = (np.dot(v2,v_stat)*v_stat)/(self.norm(v_stat)**2)
        v2_in = v2 - v2_out - v2_arm

        ang_out_1 = np.arcsin(self.norm(v1_out)/self.norm(v1))
        ang_out_1 = ang_out_1 * np.sign(np.dot(v1_out,n))
        ang_out_2 = np.arcsin(self.norm(v2_out)/self.norm(v2))
        ang_out_2 = ang_out_2 * np.sign(np.dot(v2_out,n))
        ang_out = ang_out_1 - ang_out_2
        
        q = np.cross(n,v_stat)
        q = self.unit(q) 
      
        v1_in_calc = v1 - v1_out
        v2_in_calc = v2 - v2_out
        ang_in_1 = np.arcsin(self.norm(np.cross(v1_in_calc,r))/(self.norm(v1_in_calc)*self.norm(r)))
        ang_in_2 = np.arcsin(self.norm(np.cross(v2_in_calc,r))/(self.norm(v2_in_calc)*self.norm(r)))
        
        #ang_in = ang_in_1 - ang_in_2
        #ang_in = 2*np.arctan(self.norm(v1_in)/self.norm(v1_arm))
        ang_in_1 = np.arcsin(self.norm(np.cross(v1_in_calc,v_stat))/(self.norm(v1_in_calc)*self.norm(v_stat)))
        ang_in_2 = np.arcsin(self.norm(np.cross(v2_in_calc,v_stat))/(self.norm(v2_in_calc)*self.norm(v_stat)))
        ang_in = abs(ang_in_1) + abs(ang_in_2)
        #v1_calc = v1 - v1_out
        #v2_calc = v2 - v2_out
        #ang_in = np.linalg.norm(np.cross(v1_calc,v2_calc))/(np.linalg.norm(v1_calc)*np.linalg.norm(v2_calc))
        
        #ang_in_1 = np.arcsin(self.norm(v1_in)/self.norm(v1))
        #ang_in_1 = ang_in_1 * np.sign(np.dot(v1_in,q))
        #ang_in_2 = np.arcsin(self.norm(v2_in)/self.norm(v2))
        #ang_out_2 = ang_in_2 * np.sign(np.dot(v2_in,q))
        #ang_in = ang_in_1# - ang_in_2


        #ang_in = np.arcsin(self.norm(np.cross(v1_arm,v2_arm))/(self.norm(v1_in)*self.norm(v2_in)))

        #ang_in = 2*(self.norm(v1_in)/self.norm(v_stat))

        #self.print_component(v1,v1_in,v1_out,v1_arm) 
        #self.print_component(v2,v2_in,v2_out,v2_arm) 
        #ang_in_1 = np.arcsin(self.norm(np.cross(v1_in,r)) / (self.norm(v1_in)*self.norm(r)))
        #ang_in_2 = np.arcsin(self.norm(np.cross(v2_in,r)) / (self.norm(v1_in)*self.norm(r)))
        #ang_in = ang_in_1 - ang_in_2
        #ang_out = np.arcsin(self.norm(np.cross(v1_out,v2_out)) / (self.norm(v1_out)*self.norm(v2_out)))
        #ang_out = ang_out*np.sign(np.dot(v1_out-v2_out,n))

        #v1_out=self.unit(v1_out)
        #v1_in=self.unit(v1_in)
        #v2_out=self.unit(v2_out)
        #v2_in=self.unit(v2_in)
        
        #ang_out = np.arcsin(self.norm(v1_out-v2_out)/self.norm(v1-v2))
        #ang_out = ang_out*np.sign(np.dot(v1_out-v2_out,n))
        #ang_in = np.arcsin(self.norm(v1_in-v2_in)/self.norm(v1-v2))
        #ang_in = ang_out*np.sign(np.dot(v1_in-v2_in,n))


        #ang_in = np.linalg.norm(np.cross(v1_in,v2_in))/(np.linalg.norm(v1_in)*np.linalg.norm(v2_in))
        #ang_in = np.arcsin(ang_in)
        #ang_out = np.linalg.norm(np.cross(v1_out,v2_out))/(np.linalg.norm(v1_out)*np.linalg.norm(v2_out))
        #ang_out = np.arcsin(ang_out)

        #ang_in = self.angle(v1_in,v2_in)
        #ang_out = self.angle(v1_out,v2_out)


        #v_out = np.cross(v1_out,v2_out)
        #v_out = self.unit(v_out)
        #ang_out = np.arcsin(v_out)
        #v_in = np.cross(v1_in,v2_in)
        #v_in = self.unit(v_in)
        #ang_in = np.arcsin(v_in)


        #ang_out = np.arccos(np.dot(v1_out,v2_out))
        #ang_in = np.arccos(np.dot(v1_in,v2_in))
        
        return [ang_in,ang_out]




class PAA():
    def __init__(self,**kwargs):
        self.home = kwargs.pop('home',os.getcwd())
        self.filename = kwargs.pop('filename','')
        self.directory_imp = kwargs.pop('directory_imp','')
        self.read_max = kwargs.pop('read_max','all')
        self.num_back = kwargs.pop('num_back',0)
        self.plot_on = kwargs.pop('plot_on',True)
        self.scale = kwargs.pop('scale','Default')
        if self.scale=='Default':
            print('Getting scale by filename:')
            a = self.filename
            a1 = a.split('.')[0]
            a1 = a1.split('_')
            for k in range(0,len(a1)):
                if 'scale' == a1[k]:
                    self.scale = float(a1[k+1])
                    print(self.scale)
                    
        else:
            print(self.scale)
        print('')
        
        self.method = kwargs.pop('method','fsolve')        
        self.dir_savefig = kwargs.pop('dir_savefig',False)
        self.dir_extr = kwargs.pop('dir_extr','')
        self.new_folder = kwargs.pop('new_folder',True)
        self.timeunit = kwargs.pop('timeunit','Default')
        if self.timeunit=='Default':
            print('Getting timestep by filename:')
            a = self.filename
            a1 = a.split('.')[0]
            a1 = a1.split('_')
            for k in range(0,len(a1)):
                if 'timestep' == a1[k]:
                    self.timeunit = a1[k+1]
                    print(self.timeunit)
            if self.timeunit!='days' and self.timeunit!='seconds':
                print('Could not obtain proper timestep')
        else:
            print(self.timeunit)
        print('')
            
        self.delay = kwargs.pop('delay',True)
        self.LISA = kwargs.pop('LISA',False)
        self.arm_influence = kwargs.pop('arm_influence',True)
        self.tstep = kwargs.pop('tstep',False)
        
        self.PAA_func()


    def PAA_func(self):
        home = self.home
        filename = self.filename
        directory_imp = self.directory_imp
        read_max = self.read_max
        num_back = self.num_back
        plot_on = self.plot_on
        scale = self.scale
        dir_savefig = self.dir_savefig
        dir_extr = self.dir_extr
        new_folder = self.new_folder 
        timeunit = self.timeunit
        delay = self.delay
        LISA = self.LISA
        arm_influence = self.arm_influence
        tstep = self.tstep

        print('')
        print('Importing Orbit')
        tic=time.clock()
        Orbit=orbit(home=home,filename=filename,directory_imp=directory_imp,num_back=num_back,scale=scale,read_max=read_max,plot_on=False,timeunit=timeunit)
        print(str(Orbit.linecount)+' datapoints')
        print('Done in '+str(time.clock()-tic))


        def i_slr(i):
            i_self = i
            i_left = (i+1)%3
            i_right = (i+2)%3

            i_ret = [i_self,i_left,i_right]
            for j in range(0,len(i_ret)):
                if i_ret[j]==0:
                    i_ret[j]=3

            return i_ret

        def func_error(func,t): 
            y = func(t)
            check=False
            for i in y:
                if i!=0:
                    check=True
                    break
            if check!=False:
                raise ValueError('Timestamp out of bounds')
            return y

        def func_error2(func,t):
            check=True
            try:
                func(t)
            except:
                check=False
                return np.nan
                pass

            if check==True:
                return func(t)
            
        #def r_val(v_l,v_r,m=[2,2,2]):

        def solve_num(func,guess,method='fsolve'):
            if method == 'fsolve':
                guess = np.array(guess)
                ret = scipy.optimize.fsolve(func,guess)
            elif method == 'excitingmixing':
                ret = scipy.optimize.excitingmixing(func,guess)
            elif method == 'linearmixing':
                ret = scipy.optimize.linearmixing(func,guess)
            elif method == 'newton_krylov':
                guess = [guess]
                ret = scipy.optimize.newton_krylov(func,guess)
            elif method == 'anderson':
                ret = scipy.optimize.anderson(func,guess)
            elif method == 'broyden1':
                ret = scipy.optimize.broyden1(func,guess)

            return ret
            
            
            
        def nominal_arm(orbit,i,t):
            def func_arm(i):
                L_vec=[]
                t_vec=orbit.t
                for j in range(0,len(t_vec)):
                    L_vec.append(np.linalg.norm(orbit.L[i-1][j]))

                f=interp1d(t_vec,L_vec)


                return f

            f=func_arm(i)

            return f(t)

        lisa=Orbit.lisa_obj
        func_nominal_arm = lambda i,time: nominal_arm(Orbit,i,time)

        print("")
        print('Obtaining proper LISA object')
        tic=time.clock()
        lisa_orb=PyLISA(lisa,func_nominal_arm)
        lisa_cache=CacheLISA(lisa_orb) # Works with retard() and putp and putn
        if LISA==True:
            #LISA='poep'
            LISA=lisa_cache
        else:
            LISA=False
        #lisa_cache=synthlisa.makeSampledLISA("NGO_1M_10deg_synthlisa.txt")
        print('Done in '+str(time.clock()-tic))

        t=np.array(Orbit.t)


        def func_pos(orbit,i,LISA=False):
            if LISA==False:
                i = (i%3) - 1
                t_inter=orbit.t
                [x,y,z]=[orbit.p[i][:,0],orbit.p[i][:,1],orbit.p[i][:,2]]
                fx=interp1d(t_inter,x)
                fy=interp1d(t_inter,y)
                fz=interp1d(t_inter,z)

                #fx = lambda time: func_error2(fx,time)
                #fy = lambda time: func_error2(fy,time)
                #fz = lambda time: func_error2(fz,time)

                return lambda time: np.array([fx(time),fy(time),fz(time)])
            else:
                L = lambda time: np.array(LISA.putp(i,time))
                f = lambda time: func_error2(L,time)
                return f 

        def func_arm_sec(orbit,i,side,LISA=False): # Returns arm is seconds
            if LISA==False:
                t_inter=orbit.t
                [i_self,i_left,i_right] = i_slr(i)
                pos_self = func_pos(orbit,i_self,LISA=LISA)
                if side=='r':
                    pos_right = func_pos(orbit,i_righti,LISA=LISA)

                    arms = lambda time: np.linalg.norm(pos_right(time) - pos_self(time))/c
                    #arms = orbit.v_r
                elif side=='l':
                    pos_left = func_pos(orbit,i_left,LISA=LISA)
                    arms = lambda time: np.linalg.norm(pos_left(time) - pos_self(time))/c
                    #arms = orbit.v_l
                return arms

            else:
                [i_self,i_l,i_r] = i_slr(i)
                #i_l=((i+1)%3)
                #i_r=((i+2)%3)
                #i_self=i

                if side=='l':
                    arms = lambda time: (np.array(LISA.putp(i_l,time)) - np.array(LISA.putp(i_self,time)))/c
                elif side=='r':
                    arms = lambda time: (np.array(LISA.putp(i_r,time)) - np.array(LISA.putp(i_self,time)))/c
                return arms

        self.delay_func=[] 
        def L_PAA(orbit,pos_self,pos_left,pos_right,LISA=False):
            t_inter = orbit.t

            if LISA!=False:
                tl_guess = LISA.armlength(1,0)/c
            else:
                tl_guess= np.linalg.norm(Orbit.p[0][0,:] - Orbit.p[1][0,:])/c
            tr_guess = tl_guess
            tl_sol=[]
            tr_sol=[]
            trl_sol=[]
            trr_sol=[]
            t_send_l_sol=[]
            t_send_r_sol=[]
            t_rec_l_sol=[]
            t_rec_r_sol=[]
            
            catch_lim=np.array([0,0,0])

            t_min=t_inter[0]
            t_max=t_inter[-1]
            for t in t_inter:
                try:
                    #print(type(t+dtl))
                    #print(type(t))
                    #print('')
                    s1 = lambda x: pos_left(x)
                    s2 = lambda x: pos_self(x)
                    x_0 = t
                    s3 = lambda dt: s1(x_0+dt) - s2(x_0)
                    s4 = lambda dt: np.linalg.norm(s3(dt))
                    s5 = lambda dt: s4(dt) - c*dt
                    s6 = lambda dt: s5(dt)/c

                    #L_send_l = lambda dtl: (np.linalg.norm(pos_left(t+dtl)-pos_self(t)) - c*(dtl))/c # L_rec in seconds
                    #tl_sol.append(fsolve(L_send_l,tl_guess)[0])
                    res = scipy.optimize.brentq(s6,0,tl_guess*4)
                    #if (t_min < res + t) and (res + t< t_max):
                    tl_sol.append(res)
                    t_send_l_sol.append(t)
                    
                    self.delay_func.append(s6)
                except ValueError or TypeError:
                    pass

                try:
                    s1 = lambda x: pos_right(x)
                    s2 = lambda x: pos_self(x)
                    x_0 = t
                    s3 = lambda dt: s1(x_0+dt) - s2(x_0)
                    s4 = lambda dt: np.linalg.norm(s3(dt))
                    s5 = lambda dt: s4(dt) - c*dt
                    s6 = lambda dt: s5(dt)/c

                    #L_send_r = lambda dtr: (np.linalg.norm(pos_right(t+dtr)-pos_self(t)) - c*(dtr))/c # L_rec in seconds
                    res = scipy.optimize.brentq(s6,0,tr_guess*4)
                    #if (t_min < res + t) and (res + t< t_max):
                    tr_sol.append(res)
                    t_send_r_sol.append(t)
                    self.delay_func.append(s6)

                except ValueError or TypeError:
                    pass

                try:
                    s1 = lambda x: pos_self(x)
                    s2 = lambda x: pos_left(x)
                    x_0 = t
                    s3 = lambda dt: s1(x_0) - s2(x_0-dt)
                    s4 = lambda dt: np.linalg.norm(s3(dt))
                    s5 = lambda dt: s4(dt) - c*dt
                    s6 = lambda dt: s5(dt)/c
                    self.delay_func.append(s6)

                    #L_rec_l = lambda dtrl: (np.linalg.norm(pos_self(t) - pos_left(t-dtrl)) - c*(dtrl))/c # L_rec in seconds
                    res = scipy.optimize.brentq(s6,0,tl_guess*4)
                    #if (t_min < -res + t) and (-res + t< t_max):
                    trl_sol.append(res)
                    t_rec_l_sol.append(t)

                except ValueError or TypeError:
                    pass


                try:
                    s1 = lambda x: pos_self(x)
                    s2 = lambda x: pos_right(x)
                    x_0 = t
                    s3 = lambda dt: s1(x_0) - s2(x_0-dt)
                    s4 = lambda dt: np.linalg.norm(s3(dt))
                    s5 = lambda dt: s4(dt) - c*dt
                    s6 = lambda dt: s5(dt)/c

                    #L_rec_r = lambda dtrr: (np.linalg.norm(pos_self(t) - pos_right(t-dtrr)) - c*(dtrr))/c # L_rec in seconds
                    res = scipy.optimize.brentq(s6,0,tr_guess*4)
                    #if (t_min < -res + t) and (-res + t< t_max):
                    trr_sol.append(res)
                    t_rec_r_sol.append(t)
                    self.delay_func.append(s6)

                except ValueError or TypeError:
                    pass

            L_sl=interp1d(t_send_l_sol,tl_sol) #...adjust to better fit
            L_sr=interp1d(t_send_r_sol,tr_sol)
            L_rl=interp1d(t_rec_l_sol,trl_sol)
            L_rr=interp1d(t_rec_r_sol,trr_sol)
           
            dt_l_tot = np.array(tl_sol)+np.array(trl_sol)
            dt_r_tot = np.array(tr_sol)+np.array(trr_sol)

            self.dt_sol_vec = [tl_sol,tr_sol,trl_sol,trr_sol]

            return [[L_sl, L_sr, L_rl, L_rr],[dt_l_tot,dt_r_tot]]

        def n_r_lisa(i,LISA,m=[2,2,2]):
            [i_self,i_left,i_right] = i_slr(i)

            v_l = lambda time: np.array(LISA.putp(i_left,time)) - np.array(LISA.putp(i_self,time))
            v_r = lambda time: np.array(LISA.putp(i_right,time)) - np.array(LISA.putp(i_self,time))
            COM = lambda time: (m[i_left-1]*np.array(LISA.putp(i_left,time)) + m[i_right-1]*np.array(LISA.putp(i_right,time)) + m[i_self-1]*np.array(LISA.putp(i_self,time)))/sum(m)
            
            r = lambda time: COM(time) - np.array(LISA.putp(i_self,time))

            n = lambda time: np.cross(v_l(time),v_r(time))
            #n = lambda time: n(time)/np.linalg.norm(n(time))
            #n = lambda time: n(time)/np.linalg.norm(n)
            #... Not normalized
            return [n,r]
        def r_calc(v_l,v_r,i,m=[2,2,2]):

            [i_self,i_left,i_right] = i_slr(i)
            r =  (v_l*m[i_left-1]+v_r*m[i_right-1])/sum(m)
            
            return r



       

        def send_func(orbit,i,delay=self.delay,LISA=False,arm_influence=True):
            [i_self,i_left,i_right] = i_slr(i)
            #i_left=((i+1)%3) - 1
            #i_self = i-1
            #i_right = ((i+2)%3) - 1

            pos_left = func_pos(orbit,i_left,LISA=LISA)
            pos_self = func_pos(orbit,i_self,LISA=LISA)
            pos_right = func_pos(orbit,i_right,LISA=LISA)

            if delay==True:
                if (LISA==False) or (LISA!=False and arm_influence == False):
                    pos_left_bad = func_pos(orbit,i_left,LISA=False)
                    pos_self_bad = func_pos(orbit,i_self,LISA=False)
                    pos_right_bad = func_pos(orbit,i_right,LISA=False)

                    [[L_sl,L_sr,L_rl,L_rr],[dt_l_tot,dt_r_tot]] = L_PAA(orbit,pos_self_bad,pos_left_bad,pos_right_bad,LISA=LISA)
                else:
                    [[L_sl,L_sr,L_rl,L_rr],[dt_l_tot,dt_r_tot]] = L_PAA(orbit,pos_self,pos_left,pos_right,LISA=LISA)
            elif delay=='Not ahead':
                L_sl = lambda t: np.linalg.norm(pos_left(t) - pos_self(t))/c #func_arm_sec(orbit,i_self,'l',LISA=LISA)
                L_sr = lambda t: np.linalg.norm(pos_right(t) - pos_self(t))/c #func_arm_sec(orbit,i_self,'r',LISA=LISA)
                L_rl=L_sl
                L_rr=L_sr
                dt_l_tot = False
                dt_r_tot = False

            elif delay==False:
                L_sl = lambda t: 0
                L_sr = lambda t: 0
                L_rl=L_sl
                L_rr=L_sr
                dt_l_tot = False
                dt_r_tot = False


            v_send_l = lambda t: pos_left(t+L_sl(t)) - pos_self(t)
            v_send_r = lambda t: pos_right(t+L_sr(t)) - pos_self(t)
            v_rec_l = lambda t: pos_self(t) - pos_left(t - L_rl(t))
            v_rec_r = lambda t: pos_self(t) - pos_right(t - L_rr(t))

            return [[v_send_l,v_send_r,v_rec_l,v_rec_r],[dt_l_tot,dt_r_tot]]
        
        

        LA=la()
        t_calc=[]
        v_l=[]
        v_r=[]
        u_l=[]
        u_r=[]


        t_calc_new=[]
        if tstep!=False:
            for i in range(1,4):
                tmin=Orbit.t[0]
                tmax=Orbit.t[-1]

                t_calc_new.append(np.linspace(tmin,tmax,((tmax-tmin)/tstep)+1))

        self.dt_sol=[]
        self.dt_l_tot=[]
        self.dt_r_tot=[]

        for i in range(1,4):
            t_calc_vec=[]
            v_l_calc=[]
            v_r_calc=[]
            u_l_calc=[]
            u_r_calc=[]
            [[v_l_func,v_r_func,u_l_func,u_r_func],[dt_l_tot_vec,dt_r_tot_vec]] = send_func(Orbit,i,delay=self.delay,LISA=LISA,arm_influence=arm_influence)
            self.dt_l_tot.append(dt_l_tot_vec)
            self.dt_r_tot.append(dt_r_tot_vec)

            if delay == True:
                self.dt_sol.append(self.dt_sol_vec)
            #v_l_func_na,v_r_func_na,u_l_func_na,u_r_func_na = send_func(Orbit,i,delay='Not ahead',LISA=LISA)
            #v_l_func_nd,v_r_func_nd,u_l_func_nd,u_r_func_nd = send_func(Orbit,i,delay=False,LISA=LISA)       
             
            if tstep!=False:
                t_vec = t_calc_new[i-1]
            else:
                t_vec = Orbit.t
            for t in t_vec:
                calc_check=True
                try:
                    v_l_tmp = v_l_func(t)
                    v_r_tmp = v_r_func(t)
                    u_l_tmp = u_l_func(t)
                    u_r_tmp = u_r_func(t)
                except ValueError:
                    calc_check=False
                    print('Not valid timestamp')
                    #if t>0:
                    #    print(t)
                    #    print(v_l_func(t))
                    #    print(v_r_func(t))
                    #    print(u_l_func(t))
                    #    print(u_r_func(t))

                    pass
                if calc_check==True:
                    v_l_calc.append(v_l_tmp)
                    v_r_calc.append(v_r_tmp)
                    u_l_calc.append(u_l_tmp)
                    u_r_calc.append(u_r_tmp)
                    t_calc_vec.append(t)



            v_l.append(np.array(v_l_calc))
            v_r.append(np.array(v_r_calc))
            u_l.append(np.array(u_l_calc))
            u_r.append(np.array(u_r_calc))
            t_calc.append(t_calc_vec)

        t_plot=[]
        PAA_l_in=[]
        PAA_l_out=[]
        PAA_r_in=[]
        PAA_r_out=[]
        ang_sc=[]
        ang_beam=[]
        L_l=[]
        L_r=[]
        diff_L_l=[]
        diff_L_r=[]
        ang_sr_l=[]
        ang_sr_r=[]
        ang_wob=[]
        ang_wob_stat=[]
        ang_wob_diff=[]
        ang_beam_in_l=[]
        ang_beam_out_l=[]
        ang_beam_in_r=[]
        ang_beam_out_r=[]
        pos=[]
        normal_vec=[]

        for i in range(1,4):        
            t_plot_vec=[]
            PAA_l_in_vec=[]
            PAA_l_out_vec=[]
            PAA_r_in_vec=[]
            PAA_r_out_vec=[]
            ang_sc_vec=[]
            ang_beam_vec=[]
            L_l_vec=[]
            L_r_vec=[]
            diff_L_l_vec=[]
            diff_L_r_vec=[]
            ang_sr_l_vec=[]
            ang_sr_r_vec=[]
            ang_wob_vec=[]
            ang_wob_stat_vec=[]
            ang_wob_diff_vec=[]       
            ang_beam_in_l_vec=[]
            ang_beam_out_l_vec=[]
            ang_beam_in_r_vec=[]
            ang_beam_out_r_vec=[]
            pos_vec=[]
            normal_vec_vec=[]

            [i_self,i_left,i_right] = i_slr(i)
            
            pos_left_func = func_pos(Orbit,i_left,LISA=LISA)
            pos_self_func = func_pos(Orbit,i_self,LISA=LISA)
            pos_right_func = func_pos(Orbit,i_right,LISA=LISA)
            if LISA!=False:
                [n_func,r_func]=n_r_lisa(i,LISA)
            i=i-1 
            for j in range(0,len(t_calc[i])):     
                check_good=True
                try:
                    pos_left = pos_left_func(t_calc[i][j])
                    pos_self = pos_self_func(t_calc[i][j])
                    pos_right = pos_right_func(t_calc[i][j])
                    v_l_stat = pos_left - pos_self
                    v_r_stat = pos_right - pos_self
                    
                    pos_vec.append(pos_self)

                    if LISA==False:
                        #n = Orbit.n_func[i](t_calc[i][j])
                        n = LA.unit(np.cross(v_l_stat,v_r_stat))
                        #n = Orbit.n_new_func[i](t_calc[i][j])
                        r = r_calc(v_l_stat,v_r_stat,i+1) #(v_l_stat+v_r_stat)/2.0#Orbit.r_func[i](t_calc[i][j])
                    else:                  
                        n = LA.unit(n_func(t_calc[i][j]))
                        r = r_func(t_calc[i][j])
                


                    normal_vec_vec.append(n)
                    #v_l_stat = Orbit.v_l_func[i](t_calc[i][j]) #...IS NOT WORKING
                    #v_r_stat = Orbit.v_r_func[i](t_calc[i][j])
                    
                    
                    #n = np.cross(v_l_stat,v_r_stat)
                    #n = n/np.linalg.norm(n)
                except ValueError:
                    check_good=False
                    print('not a good value')
                    pass
                if check_good == True:
                    v_l_calc=v_l[i][j,:]
                    v_r_calc=v_r[i][j,:]
                    u_l_calc=u_l[i][j,:]
                    u_r_calc=u_r[i][j,:]

                    t_plot_vec.append(t_calc[i][j])                   
                   
                    ang_in_v_l = LA.ang_in_dot(v_l_calc,v_l_stat,n,r)                   
                    ang_in_u_l = LA.ang_in_dot(-u_l_calc,v_l_stat,n,r)
                    ang_out_v_l = LA.ang_out(v_l_calc,n)
                    ang_out_u_l = LA.ang_out(-u_l_calc,n)
                    
                    ang_in_v_r = LA.ang_in_dot(v_r_calc,v_r_stat,n,r)
                    ang_in_u_r = LA.ang_in_dot(-u_r_calc,v_r_stat,n,r)
                    ang_out_v_r = LA.ang_out(v_r_calc,n)
                    ang_out_u_r = LA.ang_out(-u_r_calc,n)
                     
                    ang_sc_vec.append(LA.angle(v_l_stat,v_r_stat))
                    ang_beam_vec.append(LA.angle(v_l_calc,v_r_calc))

                    #ang_beam_in_l_vec.append(ang_in_v_l - ang_in_u_l)
                    #ang_beam_out_l_vec.append(ang_out_v_l - ang_out_u_l)
                    #ang_beam_in_r_vec.append(ang_in_v_r - ang_in_u_r)
                    #ang_beam_out_r_vec.append(ang_out_v_r - ang_out_u_r)
                    
                    #[calc_angv_l_in,calc_angv_l_out]=LA.ang_in_out(v_l_calc,r,n)
                    #[calc_angv_r_in,calc_angv_r_out]=LA.ang_in_out(v_r_calc,r,n)
 
                    #[calc_angu_l_in,calc_angu_l_out]=LA.ang_in_out(-u_l_calc,v_l_stat,n)
                    #[calc_angu_r_in,calc_angu_r_out]=LA.ang_in_out(-u_r_calc,v_r_stat,n)
                    #
                    #ang_beam_in_l_vec.append(calc_angv_l_in - calc_angu_l_in)
                    #ang_beam_out_l_vec.append(calc_angv_l_out - calc_angu_l_out)
                    #ang_beam_in_r_vec.append(calc_angv_r_in - calc_angu_r_in)
                    #ang_beam_out_r_vec.append(calc_angv_r_out - calc_angu_r_out)
                    
                    [calc_ang_l_in,calc_ang_l_out]=LA.ang_in_out(v_l_calc,-u_l_calc,n,v_l_stat,r)
                    [calc_ang_r_in,calc_ang_r_out]=LA.ang_in_out(v_r_calc,-u_r_calc,n,v_r_stat,r)
                   
                    ang_beam_in_l_vec.append(calc_ang_l_in)
                    ang_beam_out_l_vec.append(calc_ang_l_out)
                    ang_beam_in_r_vec.append(calc_ang_r_in)
                    ang_beam_out_r_vec.append(calc_ang_r_out)



                    L_l_vec.append(LA.norm(v_l_calc))
                    L_r_vec.append(LA.norm(v_r_calc))

                    diff_L_l_vec.append(LA.norm(v_l_calc) - LA.norm(v_l_stat))
                    diff_L_r_vec.append(LA.norm(v_r_calc) - LA.norm(v_r_stat))
                    ang_sr_l_vec.append(LA.angle(v_l_calc,-u_l_calc,dot=True))
                    ang_sr_r_vec.append(LA.angle(v_r_calc,-u_r_calc,dot=True))
                    ang_wob_vec.append(LA.angle((v_l_calc+v_r_calc)/2.0,r))
                    ang_wob_stat_vec.append(LA.angle((v_l_stat+v_r_stat)/2.0,r))
                    ang_wob_diff_vec.append(ang_wob_vec[-1] - ang_wob_stat_vec[-1])
                    

                    PAA_l_in_vec.append(ang_in_v_l)
                    PAA_l_out_vec.append(ang_out_v_l)
                    PAA_r_in_vec.append(ang_in_v_r)
                    PAA_r_out_vec.append(ang_out_v_r)
                 
            PAA_l_in.append(np.array(PAA_l_in_vec))
            PAA_l_out.append(np.array(PAA_l_out_vec))
            PAA_r_in.append(np.array(PAA_r_in_vec))
            PAA_r_out.append(np.array(PAA_r_out_vec))
            ang_sc.append(np.array(ang_sc_vec))
            ang_beam.append(np.array(ang_beam_vec))
            L_l.append(np.array(L_l_vec))
            L_r.append(np.array(L_r_vec))
            diff_L_l.append(np.array(diff_L_l_vec))
            diff_L_r.append(np.array(diff_L_r_vec))
            ang_sr_l.append(np.array(ang_sr_l_vec))
            ang_sr_r.append(np.array(ang_sr_r_vec))
            ang_wob.append(np.array(ang_wob_vec))
            ang_wob_stat.append(np.array(ang_wob_stat_vec))
            ang_wob_diff.append(np.array(ang_wob_diff_vec))
            ang_beam_in_l.append(np.array(ang_beam_in_l_vec))
            ang_beam_out_l.append(np.array(ang_beam_out_l_vec))
            ang_beam_in_r.append(np.array(ang_beam_in_r_vec))
            ang_beam_out_r.append(np.array(ang_beam_out_r_vec))
            pos.append(np.array(pos_vec))
            normal_vec.append(np.array(normal_vec_vec))

            t_plot.append(np.array(t_plot_vec))

        PAA_beam_next_sc = [PAA_l_in,PAA_r_out,PAA_r_in,PAA_r_out]
        PAA_ret = [ang_beam_in_l,ang_beam_out_l,ang_beam_in_r,ang_beam_out_r]
        
        other_ret=[ang_sc,ang_beam,L_l,L_r,diff_L_l,diff_L_r,ang_sr_l,ang_sr_r,ang_wob,ang_wob_stat,ang_wob_diff,PAA_beam_next_sc,t_plot]
       
        self.PAA_ret = PAA_ret
        self.PAA_beam_next_sc = PAA_beam_next_sc
        self.ang_sc = ang_sc
        self.ang_beam = ang_beam
        self.L_l = L_l
        self.L_r = L_r
        self.diff_L_l = diff_L_l
        self.diff_L_r = diff_L_r
        self.ang_sr_l = ang_sr_l
        self.ang_sr_r = ang_sr_r
        self.ang_wob = ang_wob
        self.ang_wob_stat = ang_wob_stat
        self.ang_wob_diff = ang_wob_diff
        self.t_plot = t_plot
        self.Orbit = Orbit
        self.LISA = LISA
        self.pos = pos
        self.normal_vec = normal_vec
        
        # Calculating velocity
        v_inplane=[]
        v_outplane=[]
        v_rel_l = []
        v_rel_r=[]
        v_abs=[]
        for i in range(0,len(t_calc)):
            v_abs_vec = []
            x = t_calc[i]
            y1 = pos[i][:,0]
            y2 = pos[i][:,1]
            y3 = pos[i][:,2]

            df1dt = np.diff(y1)/np.diff(x)
            df2dt = np.diff(y2)/np.diff(x)
            df3dt = np.diff(y3)/np.diff(x)

            x_new = x[0:-1]
            v_inplane_vec=[]
            v_outplane_vec=[]
            v_rel_l_vec=[]
            v_rel_r_vec = []
            for j in range(0,len(x_new)):
                n = normal_vec[i][j]
                v = np.array([df1dt[j],df2dt[j],df3dt[j]])
                inplane_calc = LA.inplane(v,n)
                v_inplane_vec.append(np.linalg.norm(inplane_calc))
                outplane_calc = v - inplane_calc
                v_outplane_vec.append(np.linalg.norm(outplane_calc)*(np.sign(np.dot(outplane_calc,n))))
                v_abs_vec.append(v)
            v_inplane.append(v_inplane_vec)
            v_outplane.append(v_outplane_vec)
            v_abs.append(np.array(v_abs_vec))

        self.v_abs = v_abs
        self.v_inplane = v_inplane
        self.v_outplane = v_outplane

        v_rel_l_inplane = []
        v_rel_r_inplane = []
        v_rel_l_outplane = []
        v_rel_r_outplane = []
        v_rel_l_alongarm= []
        v_rel_r_alongarm = []

        for i in range(0,len(self.v_abs)):
            [i_self,i_left,i_right] = i_slr(i+1)
            i_self = i_self - 1
            i_left = i_left - 1
            i_right = i_right -1
            #print(i_self,i_left,i_right)
            v_rel_l_inplane_vec = []
            v_rel_r_inplane_vec = []
            v_rel_l_outplane_vec = []
            v_rel_r_outplane_vec = []
            v_rel_l_alongarm_vec = []
            v_rel_r_alongarm_vec = []


            for j in range(0,len(self.v_abs[i])):               
                v_rel_l_calc = self.v_abs[i_left][j] - self.v_abs[i_self][j]
                v_rel_r_calc = self.v_abs[i_right][j] - self.v_abs[i_self][j]
                n = normal_vec[i][j]
               
                r_l = LA.unit(v_l[i_self][j])
                r_r = LA.unit(v_r[i_self][j])

                r_in_l = np.cross(n,r_l)
                r_in_r = np.cross(n,r_r)
                #inplane_l = LA.outplane(v_rel_l_calc,r_l)
                #inplane_r = LA.outplane(v_rel_r_calc,r_r) 
                #outplane_l = LA.outplane(v_rel_l_calc,n)*np.sign(np.dot(v_rel_l_calc,n))
                #outplane_r = LA.outplane(v_rel_r_calc,n)*np.sign(np.dot(v_rel_r_calc,n))

                outplane_l = (np.dot(v_rel_l_calc,n)*n) / (np.linalg.norm(n)**2)
                alongarm_l = (np.dot(v_rel_l_calc,r_l)*r_l) / (np.linalg.norm(r_l)**2)
                inplane_l = v_rel_l_calc - outplane_l - alongarm_l
                outplane_r = (np.dot(v_rel_r_calc,n)*n) / (np.linalg.norm(n)**2)
                alongarm_r = (np.dot(v_rel_r_calc,r_r)*r_r) / (np.linalg.norm(r_r)**2)
                inplane_r = v_rel_r_calc - outplane_r - alongarm_r


                
                #outplane_l = self.v_outplane[i_left][j] - self.v_outplane[i_self][j]
                #alongarm_l = np.dot(v_rel_l_calc - outplane_l,LA.unit(r_l))
                #inplane_l = v_rel_l_calc - outplane_l - alongarm_l
                #outplane_r = self.v_outplane[i_right][j] - self.v_outplane[i_self][j]
                #alongarm_r = np.dot(v_rel_r_calc - outplane_r,LA.unit(r_r))
                #inplane_r = v_rel_r_calc - outplane_r - alongarm_r

                #v_rel_l_inplane_vec.append(np.linalg.norm(inplane_l)*np.sign(np.dot(inplane_l,r_l)))
                #v_rel_r_inplane_vec.append(np.linalg.norm(inplane_r)*np.sign(np.dot(inplane_r,r_r)))
                #v_rel_l_outplane_vec.append(np.linalg.norm(outplane_l)*np.sign(np.dot(outplane_l,n)))
                #v_rel_r_outplane_vec.append(np.linalg.norm(outplane_r)*np.sign(np.dot(outplane_r,n)))
                v_rel_l_outplane_vec.append(np.linalg.norm(outplane_l)*np.sign(np.dot(v_rel_l_calc,n)))
                v_rel_r_outplane_vec.append(np.linalg.norm(outplane_r)*np.sign(np.dot(v_rel_r_calc,n)))
                v_rel_l_inplane_vec.append(np.linalg.norm(inplane_l)*np.sign(np.dot(inplane_l,r_in_l)))           
                v_rel_r_inplane_vec.append(np.linalg.norm(inplane_r)*np.sign(np.dot(inplane_r,r_in_r)))
                v_rel_l_alongarm_vec.append(np.linalg.norm(alongarm_l)*np.sign(np.dot(alongarm_l,r_l)))           
                v_rel_r_alongarm_vec.append(np.linalg.norm(alongarm_r)*np.sign(np.dot(alongarm_r,r_r)))
                

            v_rel_l_inplane.append(np.array(v_rel_l_inplane_vec))
            v_rel_r_inplane.append(np.array(v_rel_r_inplane_vec))
            v_rel_l_outplane.append(np.array(v_rel_l_outplane_vec))
            v_rel_r_outplane.append(np.array(v_rel_r_outplane_vec))
            v_rel_l_alongarm.append(np.array(v_rel_l_alongarm_vec))
            v_rel_r_alongarm.append(np.array(v_rel_r_alongarm_vec))
            
        self.v_rel_l_inplane = v_rel_l_inplane
        self.v_rel_r_inplane = v_rel_r_inplane
        self.v_rel_l_outplane = v_rel_l_outplane
        self.v_rel_r_outplane = v_rel_r_outplane
        self.v_rel_l_alongarm = v_rel_l_alongarm
        self.v_rel_r_alongarm = v_rel_r_alongarm
                
        ### Plotting
        if plot_on==True:
            def get_date():
                now = datetime.datetime.now()
                date=str(now.year)+str(now.month)+str(now.day)

                date=date+'-'+dir_extr
                return date

            date_str=get_date()
            filename_save=filename.split('.')[-2]
            filename_save=filename_save.split('/')[-1]
            dir_extr_new=dir_extr
            if dir_savefig==False:
                dir_savefig=os.path.dirname(os.path.realpath(__file__))+'/figures/'+date_str+'/'

            if not os.path.exists(dir_savefig):
                os.makedirs(dir_savefig)
            elif new_folder==True:
                i=1
                while os.path.exists(dir_savefig[0:-1]+'_('+str(i)+')/'):
                    i=i+1
                dir_extr_new = dir_extr+'_('+str(i)+')/'
                dir_savefig = dir_savefig[0:-1]+'_('+str(i)+')/'
                os.makedirs(dir_savefig)



            figs=[]
            if tstep!=False:
                offset = int(round((Orbit.t[1]-Orbit.t[0])/tstep))
            else:
                offset=3
            print('Offset: '+str(offset))
            #PAA=[PAA_in_l,PAA_in_r,PAA_out_l,PAA_out_r]
            f,ax = plt.subplots(4,3,figsize=(15,15))
            f.suptitle('PAA')
            plt.subplots_adjust(hspace=0.6,wspace=0.2)
            titles=['In plane left','Out of plane left','In plane right','Out of plane right']
            for i in range(0,len(ax)):
                for j in range(0,len(PAA_ret[i])):
                    minp=offset
                    maxp=len(t_plot[j])-offset
                    print('maxp: ',maxp)
                    print('')
                    x=t_plot[j][minp:maxp]
                    y=PAA_ret[i][j][minp:maxp]*1000000
                    x_adj=x/day2sec
                    ax[i,j].plot(x_adj,y,label='S'+str(j+1))
                    ax[i,j].axhline(max(y),label='Max='+"{0:.2e}".format(max(y)),color='0',linestyle='--')
                    ax[i,j].axhline(min(y),label='Min='+"{0:.2e}".format(min(y)),color='0',linestyle='--')

                    ax[i,j].set_title(titles[i])
                    ax[i,j].set_xlabel('Time (days)')
                    ax[i,j].set_ylabel('Angle (microrad)')
                    ax[i,j].legend(loc='best')
            figs.append(f)

            PAA_ret_diff=[]
            t_plot_diff=[]
            f,ax = plt.subplots(4,3,figsize=(15,15))
            f.suptitle('Derivative PAA (per time)')
            plt.subplots_adjust(hspace=0.6,wspace=0.2)        
            for i in range(0,len(ax)):
                PAA_ret_diff_vec=[]
                for j in range(0,len(PAA_ret[i])):
                    minp=offset
                    maxp=len(t_plot[j]) -offset
                    x=t_plot[j][minp:maxp]
                    y=PAA_ret[i][j][minp:maxp]
                    dfdx = np.diff(y)/np.diff(x)
                    PAA_ret_diff_vec.append(np.array(dfdx))
                    x_calc=x/day2sec
                    x_adj = x_calc[0:-1]
                    dfdx = dfdx*1000000
                    if i==0:
                        t_plot_diff.append(x_adj)
                    ax[i,j].plot(x_adj,dfdx,label='S'+str(j+1))
                    ax[i,j].axhline(max(dfdx),label='Max='+"{0:.2e}".format(max(dfdx)),color='0',linestyle='--')
                    ax[i,j].axhline(min(dfdx),label='Min='+"{0:.2e}".format(min(dfdx)),color='0',linestyle='--')

                    ax[i,j].set_title(titles[i])
                    ax[i,j].set_xlabel('Time (days)')
                    ax[i,j].set_ylabel('Angle/time (microrad/sec)')
                    ax[i,j].legend(loc='best')
                PAA_ret_diff.append(PAA_ret_diff_vec)
            figs.append(f)

            f,ax = plt.subplots(4,2,figsize=(40,15))
            f.suptitle('PAA - Overview')
            plt.subplots_adjust(hspace=0.6,wspace=0.2)        
            for i in range(0,len(ax)):
                for j in range(0,len(t_plot)):
                    minp=offset
                    maxp=len(t_plot[j])-offset
                    ax[i,0].plot(t_plot[j][minp:maxp],PAA_ret[i][j][minp:maxp]*1000000,label='SC'+str(j+1))
                    ax[i,0].set_xlabel('Time (days)')
                    ax[i,0].set_ylabel('Angle (microrad)')
                    ax[i,0].set_title(titles[i])
                    
                    ax[i,1].plot(t_plot_diff[j],PAA_ret_diff[i][j]*1000000,label='SC'+str(j+1))
                    ax[i,1].set_xlabel('Time (days)')
                    ax[i,1].set_ylabel('Angle/time (microrad/sec)')
                    ax[i,1].set_title(titles[i])
                ax[i,0].legend(loc='best')
                ax[i,1].legend(loc='best')
            figs.append(f)


            f,ax = plt.subplots(3,2,figsize=(15,15))
            f.suptitle('Breathing angles')
            plt.subplots_adjust(hspace=0.6)
            ax[0,0].plot(t_plot[0][minp:maxp],(sum(ang_sc)/len(ang_sc))[minp:maxp],label='Mean')
            ax[0,0].axhline(math.radians(60),color='0',linestyle='--')
            ax[1,0].axhline(math.radians(60),color='0',linestyle='--')
            ax[2,0].axhline(math.radians(0),color='0',linestyle='--')

            for i in range(0,len(ang_sc)):
                x=t_plot[i][minp:maxp]
                y=ang_sc[i][minp:maxp]
                minp=offset
                maxp=len(t_plot[j])-offset
                ax[0,0].plot(x,y,label='SC'+str(i+1))
                dfdx = np.diff(y)/np.diff(x)
                x_diff = x[0:-1]
                ax[0,1].plot(x_diff,dfdx,label='SC'+str(i+1))

                y=ang_beam[i][minp:maxp]
                ax[1,0].plot(x,y,label='SC'+str(i+1))
                dfdx = np.diff(y)/np.diff(x)
                ax[1,1].plot(x_diff,dfdx,label='SC'+str(i+1))

                y=ang_sc[i][minp:maxp]-ang_beam[i][minp:maxp]
                ax[2,0].plot(t_plot[i][minp:maxp],ang_sc[i][minp:maxp]-ang_beam[i][minp:maxp],label='SC'+str(i+1))
                dfdx = np.diff(y)/np.diff(x)
                ax[2,1].plot(x_diff,dfdx,label='SC'+str(i+1))

            for j in range(0,len(ax)):
                ax[j,0].legend(loc='best')
                ax[j,0].set_xlabel('Time (sec)')
                ax[j,0].set_ylabel('Angle (rad)')
                if j<2:
                    ax[j,1].legend(loc='best')
                    ax[j,1].set_xlabel('Time (sec)')
                    ax[j,1].set_ylabel('Angle/time (rad/sec)')

            ax[0,0].set_title('Angle by positions')
            ax[1,0].set_title('Angle by outgoing beams')
            ax[2,0].set_title('Difference of the two above angles')
            ax[0,1].set_title('Angle/time by positions')
            ax[1,1].set_title('Angle/time by outgoing beams')
            ax[2,1].axis('off')

            figs.append(f)

            f,ax = plt.subplots(1,2,figsize=(15,10))
            plt.subplots_adjust(wspace=2)
            f.suptitle('Angle between sending and receiving beam')
            plt.subplots_adjust(hspace=0.6,wspace=0.2)
            titles=['In plane left','OUt of plane left','In plane right','Out of plane right']
            for i in range(0,len(t_plot)):
                ax[0].plot(t_plot[i][minp:maxp]/day2sec,ang_sr_l[i][minp:maxp],label='SC'+str(i+1))        
                ax[1].plot(t_plot[i][minp:maxp]/day2sec,ang_sr_r[i][minp:maxp],label='SC'+str(i+1))        
            ax[0].set_title('Left')
            ax[0].set_xlabel('Time (days)')
            ax[0].set_ylabel('Angle (rad)')
            ax[0].legend(loc='best')
            ax[1].set_title('Right')
            ax[1].set_xlabel('Time (days)')
            ax[1].set_ylabel('Angle (rad)')
            ax[1].legend(loc='best')  

            figs.append(f)

            f,ax = plt.subplots(3,1,figsize=(15,10))
            f.suptitle('Wobbling angle')
            plt.subplots_adjust(hspace=0.6,wspace=0.2)
            for i in range(0,len(t_plot)):
                ax[0].plot(t_plot[i][minp:maxp]/day2sec,ang_wob[i][minp:maxp],label='SC'+str(i+1))
                ax[1].plot(t_plot[i][minp:maxp]/day2sec,ang_wob_stat[i][minp:maxp],label='SC'+str(i+1))
                ax[2].plot(t_plot[i][minp:maxp]/day2sec,ang_wob_diff[i][minp:maxp],label='SC'+str(i+1))
            ax[0].set_title('Dynamic positions while light is traveling')
            ax[0].set_xlabel('Time (days)')
            ax[0].set_ylabel('Angle (rad)')
            ax[0].legend(loc='best')
            ax[1].set_title('Static positions while light is traveling')
            ax[1].set_xlabel('Time (days)')
            ax[1].set_ylabel('Angle (rad)')
            ax[1].legend(loc='best')
            ax[2].set_title('Difference of the two figures above')
            ax[2].set_xlabel('Time (days)')
            ax[2].set_ylabel('Angle (rad)')
            ax[2].legend(loc='best')
            figs.append(f)

            f,ax = plt.subplots(1,2,figsize=(15,10))
            plt.subplots_adjust(wspace=2)
            f.suptitle('Armlengths')
            plt.subplots_adjust(hspace=0.6,wspace=0.2)
            for i in range(0,len(t_plot)):
                ax[0].plot(t_plot[i][minp:maxp]/day2sec,L_l[i][minp:maxp]/1000.0,label='L'+str(i+1))
                ax[1].plot(t_plot[i][minp:maxp]/day2sec,L_r[i][minp:maxp]/1000.0,label='L'+str(i+1))
            ax[0].set_title('Left')
            ax[0].set_xlabel('Time (days)')
            ax[0].set_ylabel('Length (km)')
            ax[0].legend(loc='best')
            ax[1].set_title('Right')
            ax[1].set_xlabel('Time (days)')
            ax[1].set_ylabel('Length (km)')
            ax[1].legend(loc='best')
            figs.append(f)


            f,ax = plt.subplots(2,1,figsize=(15,10))
            plt.subplots_adjust(wspace=2)
            f.suptitle('Velocity')
            plt.subplots_adjust(hspace=0.6,wspace=0.2)
            for i in range(0,len(v_inplane)):
                x = np.array(t_calc)[i][0:-1]/day2sec
                ax[0].plot(x,v_inplane[i],label='SC'+str(i+1))
                ax[1].plot(x,v_outplane[i],label='SC'+str(i+1))
            ax[0].set_title('In plane')
            ax[0].set_xlabel('Time (days)')
            ax[0].set_ylabel('Velocity (m/s)')
            ax[0].legend(loc='best')
            ax[1].set_title('Out of plane')
            ax[1].set_xlabel('Time (days)')
            ax[1].set_ylabel('Velocity (m/s)')
            ax[1].legend(loc='best')
            figs.append(f)

            f,ax = plt.subplots(3,1,figsize=(15,10))
            plt.subplots_adjust(wspace=2)
            f.suptitle('Retarded time')
            plt.subplots_adjust(hspace=0.6,wspace=0.2)
            labels = ['send left','send right','receive left','received right']
            for i in range(0,len(self.dt_sol)):
                for j in range(0,len(self.dt_sol[i])):
                    ax[i].plot(self.dt_sol[i][j],label = labels[j])
                ax[i].set_title('SC'+str(i+1))
                ax[i].set_xlabel('Timestep (AU)')
                ax[i].set_ylabel('Time (sec)')
                ax[i].legend(loc='best')
            figs.append(f)

            f,ax = plt.subplots(3,3,figsize=(15,15))
            plt.subplots_adjust(wspace=2)
            f.suptitle('Relative velocity')
            plt.subplots_adjust(hspace=0.6,wspace=0.2)
            #labels = ['send left','send right','receive left','received right']
            for i in range(0,len(self.v_rel_l_inplane)):
                x = np.array(t_calc[i][0:-1])/day2sec
                ax[i,0].set_title('Angular component of inplane velocity relative to SC '+str(i+1))
                ax[i,0].plot(x,self.v_rel_l_inplane[i],label='left')
                ax[i,0].plot(x,self.v_rel_r_inplane[i],label='right')
                ax[i,1].set_title('Out of plane velocity relative to SC '+str(i+1))
                ax[i,1].plot(x,self.v_rel_l_outplane[i],label='left')
                ax[i,1].plot(x,self.v_rel_r_outplane[i],label='right')
                ax[i,2].set_title('Radial component of inplane velocity relative to SC '+str(i+1))
                ax[i,2].plot(x,self.v_rel_l_alongarm[i],label='left')
                ax[i,2].plot(x,self.v_rel_r_alongarm[i],label='right')
                
                ax[i,0].set_xlabel('Time (days)')
                ax[i,0].set_ylabel('Velocity (m/s)')
                ax[i,0].legend(loc='best')
                ax[i,1].set_xlabel('Time (days)')
                ax[i,1].set_ylabel('Velocity (m/s)')
                ax[i,1].legend(loc='best')           
                ax[i,2].set_xlabel('Time (days)')
                ax[i,2].set_ylabel('Velocity (m/s)')
                ax[i,2].legend(loc='best')           
            figs.append(f)

            if delay==True:
                f,ax = plt.subplots(4,3,figsize=(15,15))
                plt.subplots_adjust(wspace=2)
                f.suptitle('PAA estimate')
                plt.subplots_adjust(hspace=0.6,wspace=0.2)

                for i in range(0,len(t_calc)):
                    x = np.array(t_calc[i][0:-1])/day2sec
                    y_l_in = (self.v_rel_l_inplane[i][0:len(x)]*self.dt_l_tot[i][0:len(x)])/L_l[i][0:len(x)]
                    y_l_in = y_l_in*1000000
                    y_r_in = (self.v_rel_r_inplane[i][0:len(x)]*self.dt_r_tot[i][0:len(x)])/L_r[i][0:len(x)]
                    y_r_in = y_r_in*1000000
                    
                    y_l_out = (self.v_rel_l_outplane[i][0:len(x)]*self.dt_l_tot[i][0:len(x)])/L_l[i][0:len(x)]
                    y_l_out = y_l_out*1000000
                    y_r_out = (self.v_rel_r_outplane[i][0:len(x)]*self.dt_r_tot[i][0:len(x)])/L_r[i][0:len(x)]
                    y_r_out = y_r_out*1000000

                    ax[0,i].plot(x,y_l_in,label = 'SC'+str(i+1))
                    ax[0,i].set_title('In plane left')
                    ax[1,i].plot(x,y_l_out,label = 'SC'+str(i+1))
                    ax[1,i].set_title('Out of plane left')
                    ax[2,i].plot(x,y_r_in,label = 'SC'+str(i+1))
                    ax[2,i].set_title('In plane right')
                    ax[3,i].plot(x,y_r_out,label = 'SC'+str(i+1))
                    ax[3,i].set_title('Out of plane right')
                   
                    for j in range(0,len(ax[:,i])):
                        ax[j,i].set_xlabel('Time (days)')
                        ax[j,i].set_ylabel('Angle (microrad)')
                        ax[j,i].legend(loc='best')
                figs.append(f)





            def save_fig(figs):
                titles=['-PAA','-diffPAA','-PAA_all','-Breathing_angles','-send_receive_angle','-Wobbling_angle','-Armlengths','-Velocity','-Retarded_time','-Relative_velocity']
                if delay==True:
                    titles.append('-PPA_estimate')
                for i in range(0,len(figs)):
                    title=filename_save+titles[i]+'.png'
                    figs[i].savefig(dir_savefig+title)

                    print('Figure '+title+' saved in:')
                    print(dir_savefig)

            save_fig(figs)

            print('')
            print('')


            plt.close()
            self.figs = figs

        
        self.beam = [v_l,v_r,u_l,u_r] 

        return [[lisa_cache,Orbit],[v_l,v_r,u_l,u_r],PAA_ret,other_ret]

dir_orbits='/home/ester/git/synthlisa/orbits/'
LISA_opt = True
delay=True#e'Not ahead'#False

filename_list=[]

for (dirpath, dirnames, filenames) in os.walk(dir_orbits):
    print(filenames)
    for i in filenames:
        if i.split('.')[-1]=='txt':
            a = dirpath+'/'+i
            a = a.replace('//','/')
            filename_list.append(a)

#filename_list=[filename_list[0]]
#timeunit=['days']
#dir_extr='new_1_test'
dir_extr='new_1_synthlisa_armcalc'
#dir_extr='new_4_interp_arminterp'
#timeunit=['seconds','days','days']
timeunit='Default'#['days']
arm_influence=True
length_calc=200#'all'
#tstep=3600
tstep=False
count=0
method = 'fsolve'
select='Hallion'#'Folkner_orbit_timestep_seconds_scale_10'
PAA_res={}
filename_done=[]
for i in filename_list:
    filename_name = i.split('/')[-1]
    if i == filename_list[0]:
        new_folder=False # Adjust if you (don't) want to override
    else:
        new_folder=False
    print('Dir_extr:'+dir_extr)
    if select == 'all':
        if '/try/' in i:
            execute = False
        else:
            execute = True
    else:
        if select in i:
            execute = True
        else:
            execute = False

    if filename_name in filename_done:
        execute = False
    
    if execute == True:
        filename_save = i.split('/')[-1].split('_')[0]
        PAA_res[filename_save]=PAA(home = home,filename = i,directory_imp=False,read_max = length_calc,plot_on=True,dir_extr=dir_extr,new_folder=new_folder,timeunit=timeunit,LISA=LISA_opt,arm_influence=arm_influence,tstep=tstep,delay=delay,method=method)
        PAA_res[str(count+1)] = PAA_res[filename_save]
        filename_done.append(filename_name)
        count=count+1

                  



