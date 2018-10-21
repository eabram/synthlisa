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
#warnings.filterwarnings("error")
from scipy.optimize import fsolve

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

def PAA(home,filename,directory_imp,read_max='all',num_back=0,plot_on=True,scale=1000,dir_savefig=False,dir_extr='',new_folder=True,timeunit='days',delay=True,LISA=False,arm_influence=True,tstep=False):
    
    if timeunit=='Default':
        print('Getting timestep by filename:')
        a = filename
        a1 = a.split('.')[0]
        a1 = a1.split('_')
        for k in range(0,len(a1)):
            if 'timestep' == a1[k]:
                timeunit = a1[k+1]
                print(timeunit)
        if timeunit!='days' and timeunit!='seconds':
            print('Could not obtain proper timestep')
        
    else:
        print(timeunit)
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
            #i = (i%3)-1
            t_inter=orbit.t
            [i_self,i_left,i_right] = i_slr(i)
            if side=='r':
                arms = lambda time: func_pos(orbit,i_right)(time) - func_pos(orbit,i_self)(time)
                #arms = orbit.v_r
            elif side=='l':
                arms = lambda time: func_pos(orbit,i_left)(time) - func_pos(orbit,i_self)(time)
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
    
    def L_PAA(orbit,pos_self,pos_left,pos_right,LISA=LISA):
        t_inter = orbit.t
        tl_guess=3#lisa_cache.armlength(i,t_inter[0])/c
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
                s5 = lambda dt: s4(dt) - c*(dt)
                s6 = lambda dt: s5(dt)/c

                #L_send_l = lambda dtl: (np.linalg.norm(pos_left(t+dtl)-pos_self(t)) - c*(dtl))/c # L_rec in seconds
                #tl_sol.append(fsolve(L_send_l,tl_guess)[0])
                res = fsolve(s6,tl_guess)[0]
                #if (t_min < res + t) and (res + t< t_max):
                tl_sol.append(res)
                t_send_l_sol.append(t)

            except ValueError or TypeError:
                pass

            try:
                s1 = lambda x: pos_right(x)
                s2 = lambda x: pos_self(x)
                x_0 = t
                s3 = lambda dt: s1(x_0+dt) - s2(x_0)
                s4 = lambda dt: np.linalg.norm(s3(dt))
                s5 = lambda dt: s4(dt) - c*(dt)
                s6 = lambda dt: s5(dt)/c

                #L_send_r = lambda dtr: (np.linalg.norm(pos_right(t+dtr)-pos_self(t)) - c*(dtr))/c # L_rec in seconds
                res = fsolve(s6,tr_guess)[0]
                #if (t_min < res + t) and (res + t< t_max):
                tr_sol.append(res)
                t_send_r_sol.append(t)

            except ValueError or TypeError:
                pass

            try:
                s1 = lambda x: pos_self(x)
                s2 = lambda x: pos_left(x)
                x_0 = t
                s3 = lambda dt: s1(x_0) - s2(x_0-dt)
                s4 = lambda dt: np.linalg.norm(s3(dt))
                s5 = lambda dt: s4(dt) - c*(dt)
                s6 = lambda dt: s5(dt)/c

                #L_rec_l = lambda dtrl: (np.linalg.norm(pos_self(t) - pos_left(t-dtrl)) - c*(dtrl))/c # L_rec in seconds
                res = fsolve(s6,tl_guess)[0]
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
                s5 = lambda dt: s4(dt) - c*(dt)
                s6 = lambda dt: s5(dt)/c

                #L_rec_r = lambda dtrr: (np.linalg.norm(pos_self(t) - pos_right(t-dtrr)) - c*(dtrr))/c # L_rec in seconds
                res = fsolve(s6,tr_guess)[0]
                #if (t_min < -res + t) and (-res + t< t_max):
                trr_sol.append(res)
                t_rec_r_sol.append(t)

            except ValueError or TypeError:
                pass

        L_sl=interp1d(t_send_l_sol,tl_sol) #...adjust to better fit
        L_sr=interp1d(t_send_r_sol,tr_sol)
        L_rl=interp1d(t_rec_l_sol,trl_sol)
        L_rr=interp1d(t_rec_r_sol,trr_sol)
        
        return L_sl, L_sr, L_rl, L_rr

    def n_r_lisa(i,LISA,m=[2,2,2]):
        [i_self,i_left,i_right] = i_slr(i)

        v_l = lambda time: np.array(LISA.putp(i_left,time)) - np.array(LISA.putp(i_self,time))
        v_r = lambda time: np.array(LISA.putp(i_right,time)) - np.array(LISA.putp(i_self,time))
        COM = lambda time: (m[i_left-1]*np.array(LISA.putp(i_left,time)) + m[i_right-1]*np.array(LISA.putp(i_right,time)) + m[i_self-1]*np.array(LISA.putp(i_self,time)))/sum(m)
        
        r = lambda time: COM(time) - np.array(LISA.putp(i_self,time))

        n = lambda time: np.cross(v_l(time),v_r(time))
        #n = lambda time: n(time)/np.linalg.norm(n)
        #... Not normalized
        return [n,r]
    def r_calc(v_l,v_r,i,m=[2,2,2]):

        [i_self,i_left,i_right] = i_slr(i)
        r =  (v_l*m[i_left-1]+v_r*m[i_right-1])/sum(m)
        
        return r



   

    def send_func(orbit,i,delay=True,LISA=False,arm_influence=True):
        [i_self,i_left,i_right] = i_slr(i)
        #i_left=((i+1)%3) - 1
        #i_self = i-1
        #i_right = ((i+2)%3) - 1

        pos_left = func_pos(orbit,i_left,LISA=LISA)
        pos_self = func_pos(orbit,i_self,LISA=LISA)
        pos_right = func_pos(orbit,i_right,LISA=LISA)

        if delay==True:
            if LISA!=False and arm_influence == False:
                pos_left_bad = func_pos(orbit,i_left,LISA=False)
                pos_self_bad = func_pos(orbit,i_self,LISA=False)
                pos_right_bad = func_pos(orbit,i_right,LISA=False)

                L_sl,L_sr,L_rl,L_rr = L_PAA(orbit,pos_self_bad,pos_left_bad,pos_right_bad)
            else:
                L_sl,L_sr,L_rl,L_rr = L_PAA(orbit,pos_self,pos_left,pos_right)
        elif delay=='Not ahead':
            L_sl = func_arm_sec(orbit,i_self,'l',LISA=LISA)
            L_sr = func_arm_sec(orbit,i_self,'r',LISA=LISA)
            L_rl=L_sl
            L_rr=L_sr
        elif delay==False:
            L_sl = lambda: 0
            L_sr = lambda: 0
            L_rl=L_sl
            L_rr=L_sr

        v_send_l = lambda t: pos_left(t+L_sl(t)) - pos_self(t)
        v_send_r = lambda t: pos_right(t+L_sr(t)) - pos_self(t)
        v_rec_l = lambda t: pos_self(t) - pos_left(t - L_rl(t))
        v_rec_r = lambda t: pos_self(t) - pos_right(t - L_rr(t))

        return v_send_l,v_send_r,v_rec_l,v_rec_r
    
    

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

    for i in range(1,4):
        t_calc_vec=[]
        v_l_calc=[]
        v_r_calc=[]
        u_l_calc=[]
        u_r_calc=[]
        v_l_func,v_r_func,u_l_func,u_r_func = send_func(Orbit,i,delay=delay,LISA=LISA,arm_influence=arm_influence)
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
#
#            calc_check=True
#            try:
#                v_l_tmp_na = v_l_func_na(t)
#                v_r_tmp_na = v_r_func_na(t)
#                u_l_tmp_na = u_l_func_na(t)
#                u_r_tmp_na = u_r_func_na(t)
#            except ValueError:
#                calc_check=False
#                pass
#            if calc_check==True:
#                v_l_calc_na.append(v_l_tmp_na)
#                v_r_calc_na.append(v_r_tmp_na)
#                u_l_calc_na.append(u_l_tmp_na)
#                u_r_calc_na.append(u_r_tmp_na)
#                t_calc_vec.append(t)
#
#            calc_check=True
#            try:
#                v_l_tmp = v_l_func(t)
#                v_r_tmp = v_r_func(t)
#                u_l_tmp = u_l_func(t)
#                u_r_tmp = u_r_func(t)
#            except ValueError:
#                calc_check=False
#                pass
#            if calc_check==True:
#                v_l_calc.append(v_l_tmp)
#                v_r_calc.append(v_r_tmp)
#                u_l_calc.append(u_l_tmp)
#                u_r_calc.append(u_r_tmp)
#                t_calc_vec.append(t)







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

                if LISA==False:
                    #n = Orbit.n_func[i](t_calc[i][j])
                    n = LA.unit(np.cross(v_l_stat,v_r_stat))
                    #n = Orbit.n_new_func[i](t_calc[i][j])
                    r = r_calc(v_l_stat,v_r_stat,i+1) #(v_l_stat+v_r_stat)/2.0#Orbit.r_func[i](t_calc[i][j])
                else:                  
                    n = LA.unit(n_func(t_calc[i][j]))
                    r = r_func(t_calc[i][j])
            



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
                
                #[ang_in_v_l,ang_out_v_l] = LA.ang_in_out(v_l_calc,n,r)
                #[ang_in_v_r,ang_out_v_r] = LA.ang_in_out(v_r_calc,n,r)
                #[ang_in_v_l_stat,ang_out_v_l_stat] = LA.ang_in_out(v_l_stat,n,r)
                #[ang_in_v_r_stat,ang_out_v_r_stat] = LA.ang_in_out(v_r_stat,n,r)

                #PAA_l_in_vec.append(ang_in_v_l - ang_in_v_l_stat)
                #PAA_l_out_vec.append(ang_out_v_l)
                #PAA_r_in_vec.append(ang_in_v_r - ang_in_v_r_stat)
                #PAA_r_out_vec.append(ang_out_v_r)
                
                ang_in_v_l = LA.ang_in_dot(v_l_calc,v_l_stat,n,r)
                ang_in_u_l = LA.ang_in_dot(-u_l_calc,v_l_stat,n,r)
                
                #ang_in_v_l = LA.ang_in_direct(v_l_calc,v_l_stat,n,r)
                #ang_in_v_l = LA.ang_in(v_l_calc,n,r)
                #ang_in_v_l = LA.angle(v_l_calc,v_l_stat) # Total angle
                
                ang_out_v_l = LA.ang_out(v_l_calc,n)
                ang_out_u_l = LA.ang_out(-u_l_calc,n)
                
                ang_in_v_r = LA.ang_in_dot(v_r_calc,v_r_stat,n,r)
                ang_in_u_r = LA.ang_in_dot(-u_r_calc,v_r_stat,n,r)
                #ang_in_v_r = LA.ang_in_direct(v_r_calc,v_r_stat,n,r)
                #ang_in_v_r = LA.ang_in(v_r_calc,n,r)
                #ang_in_v_r = LA.angle(v_r_calc,v_r_stat) # Total angle

                ang_out_v_r = LA.ang_out(v_r_calc,n)
                ang_out_u_r = LA.ang_out(-u_r_calc,n)
                 
                #ang_in_v_l_stat = LA.ang_in(v_l_stat,n,r)
                #ang_out_v_l_stat = LA.ang_out(v_l_stat,n)
                #ang_in_v_r_stat = LA.ang_in(v_r_stat,n,r)
                #ang_out_v_r_stat = LA.ang_out(v_r_stat,n)
                

                #[ang_in_v_l_stat,ang_out_v_l_stat] = LA.ang_in_out(v_l_stat,n,r)
                #[ang_in_v_r_stat,ang_out_v_r_stat] = LA.ang_in_out(v_r_stat,n,r)
                #print('ang_vl: ',math.degrees(LA.angle(v_l_calc,v_l_stat))) #...IMPORTANT CHECK
                #print('ang_vr: ',math.degrees(LA.angle(v_r_calc,v_r_stat))) #...IMPORTANT CHECK
                #print(v_r_stat)

                #print('')
                #print(ang_in_v_l)
                #print(ang_in_v_l_stat)

                ang_sc_vec.append(LA.angle(v_l_stat,v_r_stat))
                ang_beam_vec.append(LA.angle(v_l_calc,v_r_calc))
                ang_beam_in_l_vec.append(ang_in_v_l - ang_in_u_l)
                ang_beam_out_l_vec.append(ang_out_v_l - ang_out_u_l)
                ang_beam_in_r_vec.append(ang_in_v_r - ang_in_u_r)
                ang_beam_out_r_vec.append(ang_out_v_r - ang_out_u_r)

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

        t_plot.append(np.array(t_plot_vec))

    PAA_beam_next_sc = [PAA_l_in,PAA_r_out,PAA_r_in,PAA_r_out]
    PAA_ret = [ang_beam_in_l,ang_beam_out_l,ang_beam_in_r,ang_beam_out_r]
    
    other_ret=[ang_sc,ang_beam,L_l,L_r,diff_L_l,diff_L_r,ang_sr_l,ang_sr_r,ang_wob,ang_wob_stat,ang_wob_diff,PAA_beam_next_sc,t_plot]

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
        f,ax = plt.subplots(4,3,figsize=(40,15))
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
                ax[i,j].axhline(min(y),label='Min='+"{0:.2e}".format(max(y)),color='0',linestyle='--')

                ax[i,j].set_title(titles[i])
                ax[i,j].set_xlabel('Time (days)')
                ax[i,j].set_ylabel('Angle (microrad)')
                ax[i,j].legend(loc='best')
        figs.append(f)

        PAA_ret_diff=[]
        t_plot_diff=[]
        f,ax = plt.subplots(4,3,figsize=(40,15))
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
            #dfdx = np.diff(y)/np.diff(x)
            #ax[2,1].plot(x_diff,dfdx,label='SC'+str(i+1))

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




        def save_fig(figs):
            titles=['-PAA','-diffPAA','-PAA_all','-Breathing_angles','-send_receive_angle','-Wobbling_angle','-Armlengths']
            for i in range(0,len(figs)):
                title=filename_save+titles[i]+'.png'
                figs[i].savefig(dir_savefig+title)

                print('Figure '+title+' saved in:')
                print(dir_savefig)

        save_fig(figs)

        print('')
        print('')


        plt.close()


    return [[lisa_cache,Orbit],[v_l,v_r,u_l,u_r],PAA_ret,other_ret]



                







dir_orbits='/home/ester/git/synthlisa/orbits/new/'
LISA_opt = True 
#LISA_opt = False

filename_list=[]

for (dirpath, dirnames, filenames) in os.walk(dir_orbits):
    print(filenames)
    for i in filenames:
        if i.split('.')[-1]=='txt':
            filename_list.append(dirpath+i)

#filename_list=[filename_list[0]]
#timeunit=['days']
#dir_extr='new_1_test'
#dir_extr='new_1_synthlisa_armcalc_tstep_1'
dir_extr='new_1_interp_interp'
#timeunit=['seconds','days','days']
timeunit='Default'#['days']
arm_influence=False#True
length_calc=20#'all'
#tstep=3600
tstep=False
count=0
for i in filename_list:
    if i == filename_list[0]:
        new_folder=False # Adjust if you (don't) want to override
    else:
        new_folder=False
    print('Dir_extr:'+dir_extr)
    [[LISA,Orbit],[v_l,v_r,u_l,u_r],PAA_ret,other_ret] = PAA(home,i,False,length_calc,plot_on=True,dir_extr=dir_extr,new_folder=new_folder,timeunit=timeunit,LISA=LISA_opt,arm_influence=arm_influence,tstep=tstep)
    count=count+1

              




