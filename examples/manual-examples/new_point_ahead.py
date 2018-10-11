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
            print('norm v1: '+str(norm_v1))
            print('norm v2: '+str(norm_v2))

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

def PAA(home,filename,directory_imp,read_max='all',num_back=0,plot_on=True,scale=1000,dir_savefig=False,dir_extr='',new_folder=True,timeunit='days',delay=True):

    print('Importing Orbit')
    tic=time.clock()
    Orbit=orbit(home=home,filename=filename,directory_imp=directory_imp,num_back=num_back,scale=scale,read_max=read_max,plot_on=False,timeunit=timeunit)
    print(str(Orbit.linecount)+' datapoints')
    print('Done in '+str(time.clock()-tic))



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
    #lisa_cache=synthlisa.makeSampledLISA("NGO_1M_10deg_synthlisa.txt")
    print('Done in '+str(time.clock()-tic))

    t=np.array(Orbit.t)


    def func_pos(orbit,i):
        i = (i%3) - 1
        t_inter=orbit.t
        [x,y,z]=[orbit.p[i][:,0],orbit.p[i][:,1],orbit.p[i][:,2]]
        fx=interp1d(t_inter,x)
        fy=interp1d(t_inter,y)
        fz=interp1d(t_inter,z)

        return lambda time: np.array([fx(time),fy(time),fz(time)])


    def func_arm_sec(orbit,i,side): # Returns arm is seconds
        i = (i%3)-1
        t_inter=orbit.t
        
        if side=='r':
            arms = orbit.v_r
        elif side=='l':
            arms = orbit.v_l

        L=[]
        for j in range(0,len(t_inter)):
            L.append(np.linalg.norm(arms[i][j,:])/c)
        L_max=max(L)
        return interp1d(t_inter,L),L_max

    def L_PAA(orbit,pos_self,pos_left,pos_right):
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

        for t in t_inter:
            try:
                L_send_l = lambda dtl: (np.linalg.norm(pos_left(t+dtl)-pos_self(t)) - c*(dtl))/c # L_rec in seconds
                tl_sol.append(fsolve(L_send_l,tl_guess)[0])
                t_send_l_sol.append(t)
            except ValueError:
                pass

            try:
                L_send_r = lambda dtr: (np.linalg.norm(pos_right(t+dtr)-pos_self(t)) - c*(dtr))/c # L_rec in seconds
                tr_sol.append(fsolve(L_send_r,tr_guess)[0])
                t_send_r_sol.append(t)

            except ValueError:
                pass

            try:
                L_rec_l = lambda dtrl: (np.linalg.norm(pos_self(t) - pos_left(t-dtrl)) - c*(dtrl))/c # L_rec in seconds
                trl_sol.append(fsolve(L_rec_l,tl_guess)[0])
                t_rec_l_sol.append(t)

            except ValueError:
                pass


            try:
                L_rec_r = lambda dtrr: (np.linalg.norm(pos_self(t) - pos_right(t-dtrr)) - c*(dtrr))/c # L_rec in seconds
                trr_sol.append(fsolve(L_rec_r,tr_guess)[0])
                t_rec_r_sol.append(t)

            except ValueError:
                pass

        L_sl=interp1d(t_send_l_sol,tl_sol)
        L_sr=interp1d(t_send_r_sol,tr_sol)
        L_rl=interp1d(t_rec_l_sol,trl_sol)
        L_rr=interp1d(t_rec_r_sol,trr_sol)
        
        return L_sl, L_sr, L_rl, L_rr


    def send_func(orbit,i,delay=True):
        i_left=((i+1)%3) - 1
        i_self = i-1
        i_right = ((i+2)%3) - 1

        pos_left = func_pos(orbit,i_left)
        pos_self = func_pos(orbit,i_self)
        pos_right = func_pos(orbit,i_right)

        if delay==True:
            L_sl,L_sr,L_rl,L_rr = L_PAA(orbit,pos_self,pos_left,pos_right)
        elif delay=='Not ahead':
            L_sl = func_arm_sec(orbit,i_self,'l')[0]
            L_sr = funv_arm_sec(orbit,i_self,'r')[0]
            L_rl=L_sl
            L_rr=L_sr
        else:
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

    for i in range(1,4):
        t_calc_vec=[]
        v_l_calc=[]
        v_r_calc=[]
        u_l_calc=[]
        u_r_calc=[]
        v_l_func,v_r_func,u_l_func,u_r_func = send_func(Orbit,i,delay=delay)
        for t in Orbit.t:
            calc_check=True
            try:
                v_l_tmp = v_l_func(t)
                v_r_tmp = v_r_func(t)
                u_l_tmp = u_l_func(t)
                u_r_tmp = u_r_func(t)
            except ValueError:
                calc_check=False
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


        i_left = ((i+1)%3) - 1
        i_self = (i%3) - 1
        i_right = ((i+2)%3) - 1
        i=i-1
        for j in range(0,len(t_calc[i])):     
            try:
                n = Orbit.n_func[i](t_calc[i][j])
                #n = Orbit.n_new_func[i](t_calc[i][j])
                r = Orbit.r_func[i](t_calc[i][j])

                pos_left = func_pos(Orbit,i_left)
                pos_self = func_pos(Orbit,i_self)
                pos_right = func_pos(Orbit,i_right)

                v_l_stat = pos_left(t_calc[i][j]) - pos_self(t_calc[i][j])
                v_r_stat = pos_right(t_calc[i][j]) - pos_self(t_calc[i][j])

                #v_l_stat = Orbit.v_l_func[i](t_calc[i][j]) #...IS NOT WORKING
                #v_r_stat = Orbit.v_r_func[i](t_calc[i][j])
                pos_self
                
                #n = np.cross(v_l_stat,v_r_stat)
                #n = n/np.linalg.norm(n)

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
                #ang_in_v_l = LA.ang_in_direct(v_l_calc,v_l_stat,n,r)
                #ang_in_v_l = LA.ang_in(v_l_calc,n,r)
                #ang_in_v_l = LA.angle(v_l_calc,v_l_stat) # Total angle
                
                ang_out_v_l = LA.ang_out(v_l_calc,n)
                
                ang_in_v_r = LA.ang_in_dot(v_r_calc,v_r_stat,n,r)
                #ang_in_v_r = LA.ang_in_direct(v_r_calc,v_r_stat,n,r)
                #ang_in_v_r = LA.ang_in(v_r_calc,n,r)
                #ang_in_v_r = LA.angle(v_r_calc,v_r_stat) # Total angle

                ang_out_v_r = LA.ang_out(v_r_calc,n)
                 
                #ang_in_v_l_stat = LA.ang_in(v_l_stat,n,r)
                #ang_out_v_l_stat = LA.ang_out(v_l_stat,n)
                #ang_in_v_r_stat = LA.ang_in(v_r_stat,n,r)
                #ang_out_v_r_stat = LA.ang_out(v_r_stat,n)
                

                #[ang_in_v_l_stat,ang_out_v_l_stat] = LA.ang_in_out(v_l_stat,n,r)
                #[ang_in_v_r_stat,ang_out_v_r_stat] = LA.ang_in_out(v_r_stat,n,r)
                #print(math.degrees(LA.angle(v_l_calc,v_l_stat))) #...IMPORTANT CHECK
                #print('')
                #print(ang_in_v_l)
                #print(ang_in_v_l_stat)

                ang_sc_vec.append(LA.angle(v_l_stat,v_r_stat))
                ang_beam_vec.append(LA.angle(v_l_calc,v_r_calc))
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
            
            except ValueError:
                pass
        
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
        t_plot.append(np.array(t_plot_vec))

    PAA_ret = [PAA_l_in,PAA_r_out,PAA_r_in,PAA_r_out]

    other_ret=[ang_sc,ang_beam,L_l,L_r,diff_L_l,diff_L_r,ang_sr_l,ang_sr_r,ang_wob,ang_wob_stat,ang_wob_diff]

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

        #PAA=[PAA_in_l,PAA_in_r,PAA_out_l,PAA_out_r]
        f,ax = plt.subplots(4,3,figsize=(40,15))
        f.suptitle('Point ahead angle')
        plt.subplots_adjust(hspace=0.6,wspace=0.2)
        titles=['In plane left','Out of plane left','In plane right','Out of plane right']
        for i in range(0,len(ax)):
            for j in range(0,len(PAA_ret[i])):
                x=t_plot[j]
                y=PAA_ret[i][j]
                x_adj=x/day2sec
                ax[i,j].plot(x_adj,y,label='S'+str(j+1))
                ax[i,j].axhline(max(y),label='Max='+"{0:.2e}".format(max(y)),color='0',linestyle='--')
                ax[i,j].axhline(min(y),label='Min='+"{0:.2e}".format(max(y)),color='0',linestyle='--')

                ax[i,j].set_title(titles[i])
                ax[i,j].set_xlabel('Time (days)')
                ax[i,j].set_ylabel('Angle (rad)')
                ax[i,j].legend(loc='best')
        figs.append(f)

        PAA_ret_diff=[]
        t_plot_diff=[]
        f,ax = plt.subplots(4,3,figsize=(40,15))
        f.suptitle('Point ahead angle per time')
        plt.subplots_adjust(hspace=0.6,wspace=0.2)        
        for i in range(0,len(ax)):
            PAA_ret_diff_vec=[]
            for j in range(0,len(PAA_ret[i])):
                x=t_plot[j]
                y=PAA_ret[i][j]
                dfdx = np.diff(y)/np.diff(x)
                PAA_ret_diff_vec.append(np.array(dfdx))
                x_calc=x/day2sec
                x_adj = x_calc[0:-1]
                if i==0:
                    t_plot_diff.append(x_adj)
                ax[i,j].plot(x_adj,dfdx,label='S'+str(j+1))
                ax[i,j].axhline(max(dfdx),label='Max='+"{0:.2e}".format(max(dfdx)),color='0',linestyle='--')
                ax[i,j].axhline(min(dfdx),label='Min='+"{0:.2e}".format(min(dfdx)),color='0',linestyle='--')

                ax[i,j].set_title(titles[i])
                ax[i,j].set_xlabel('Time (days)')
                ax[i,j].set_ylabel('Angle/time (rad/sec)')
                ax[i,j].legend(loc='best')
            PAA_ret_diff.append(PAA_ret_diff_vec)
        figs.append(f)

        f,ax = plt.subplots(4,2,figsize=(40,15))
        f.suptitle('Point ahead angle - Overview')
        plt.subplots_adjust(hspace=0.6,wspace=0.2)        
        for i in range(0,len(ax)):
            for j in range(0,len(t_plot)):
                ax[i,0].plot(t_plot[j],PAA_ret[i][j],label='SC'+str(j+1))
                ax[i,0].set_xlabel('Time (days)')
                ax[i,0].set_ylabel('Angle (rad)')
                ax[i,0].set_title(titles[i])
                
                ax[i,1].plot(t_plot_diff[j],PAA_ret_diff[i][j],label='SC'+str(j+1))
                ax[i,1].set_xlabel('Time (days)')
                ax[i,1].set_ylabel('Angle/time (rad/sec)')
                ax[i,1].set_title(titles[i])
            ax[i,0].legend(loc='best')
            ax[i,1].legend(loc='best')
        figs.append(f)


        f,ax = plt.subplots(3,2,figsize=(15,15))
        f.suptitle('Breathing angles')
        plt.subplots_adjust(hspace=0.6)
        ax[0,0].axhline(math.radians(60),color='0',linestyle='--')
        ax[1,0].axhline(math.radians(60),color='0',linestyle='--')
        ax[2,0].axhline(math.radians(0),color='0',linestyle='--')

        for i in range(0,len(ang_sc)):
            x=t_plot[i]
            y=ang_sc[i]
            ax[0,0].plot(x,y,label='SC'+str(i+1))
            dfdx = np.diff(y)/np.diff(x)
            x_diff = x[0:-1]
            ax[0,1].plot(x_diff,dfdx,label='SC'+str(i+1))

            y=ang_beam[i]
            ax[1,0].plot(x,y,label='SC'+str(i+1))
            dfdx = np.diff(y)/np.diff(x)
            ax[1,1].plot(x_diff,dfdx,label='SC'+str(i+1))

            y=ang_sc[i]-ang_beam[i]
            ax[2,0].plot(t_plot[i],ang_sc[i]-ang_beam[i],label='SC'+str(i+1))
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
            ax[0].plot(t_plot[i]/day2sec,ang_sr_l[i],label='SC'+str(i+1))        
            ax[1].plot(t_plot[i]/day2sec,ang_sr_r[i],label='SC'+str(i+1))        
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
            ax[0].plot(t_plot[i]/day2sec,ang_wob[i],label='SC'+str(i+1))
            ax[1].plot(t_plot[i]/day2sec,ang_wob_stat[i],label='SC'+str(i+1))
            ax[2].plot(t_plot[i]/day2sec,ang_wob_diff[i],label='SC'+str(i+1))
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







#L_l,L_r,diff_L_l,diff_L_r








        def save_fig(figs):
            titles=['-PAA','-diffPAA','-PAA_all','-Breathing_angles','-send_receive_angle','-Wobbling_angle']
            for i in range(0,len(figs)):
                title=filename_save+titles[i]+'.png'
                figs[i].savefig(dir_savefig+title)

                print('Figure '+title+' saved in:')
                print(dir_savefig)

        save_fig(figs)

        print('')
        print('')


        plt.close()


    return [Orbit,[v_l,v_r,u_l,u_r],PAA_ret,other_ret]



                







dir_orbits='/home/ester/git/synthlisa/orbits/'


filename_list=[]

for (dirpath, dirnames, filenames) in os.walk(dir_orbits):
    print(filenames)
    for i in filenames:
        if i.split('.')[-1]=='txt':
            filename_list.append(dirpath+i)

#filename_list=[filename_list[1]]
#timeunit=['seconds']
dir_extr='new-n_old-'
timeunit=['seconds','days','days']

count=0
for i in filename_list:
    if i == filename_list[0]:
        new_folder=False # Adjust if you (don't) want to override
    else:
        new_folder=False
    print('Dir_extr:'+dir_extr)
    [Orbit,[v_l,v_r,u_l,u_r],PAA_ret,other_ret] = PAA(home,i,False,'all',plot_on=True,dir_extr=dir_extr,new_folder=new_folder,timeunit=timeunit[count])
    count=count+1

              




