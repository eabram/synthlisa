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
warnings.filterwarnings("error")
from scipy.optimize import fsolve

year2sec=32536000
c=300000000
noise_check=False
#noise_check=True

home='/home/ester/Dropbox/Master_Project/Synthetic_LISA/'
filename='NGO_1M_10deg_synthlisa.txt'
directory_imp='LISAorbits/NGO_1M_10deg/'
dir_savefig='/home/ester/git/synthlisa/figures/'
#home='/home/ester/git/synthlisa/'
#directory_imp='lisasim/data/'
#filename='positions.txt'
num_back=0




def PAA(home,filename,directory_imp,read_max='all',num_back=0,plot_on=True,scale=1000,dir_savefig=False):
    
    date_str='zzz'
    filename_save=filename.split('.')[-2]
    filename_save=filename_save.split('/')[-1]
    if dir_savefig==False:
        dir_savefig=os.path.dirname(os.path.realpath(__file__))+'/figures/'+date_str+'/'

    if not os.path.exists(dir_savefig):
        os.makedirs(dir_savefig)

    print('Imprting Orbit: '+filename)
    tic=time.clock()
    Orbit=orbit(home=home,filename=filename,directory_imp=directory_imp,num_back=num_back,scale=scale,read_max=read_max,plot_on=False)
    print('Done in '+str(time.clock()-tic))

    def angle(v1, v2):
        return math.degrees(math.asin(np.cross(v1,v2)/(np.linalg.norm(v1)*np.linalg.norm(v2))))

        #return math.degrees(math.acos(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))))

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
        if i>3:
            i=i-3
        t_inter=orbit.t
        i=i-1 
        [x,y,z]=[orbit.p[i][:,0],orbit.p[i][:,1],orbit.p[i][:,2]]
        fx=interp1d(t_inter,x)
        fy=interp1d(t_inter,y)
        fz=interp1d(t_inter,z)

        return lambda time: np.array([fx(time),fy(time),fz(time)])


    def func_arm_sec(lisa,orbit,i): # Returns arm is seconds
        if i>3:
            i=i-3
        #print(orbit.t)
        t_inter=orbit.t
        L=[]
        for j in range(0,len(t_inter)):
            L.append(lisa.armlength(i,t_inter[j])/c)
        L_max=max(L)
        return interp1d(t_inter,L),L_max

    def send_res2ang(send,receive,orbit):
        s_vec=[]
        r_vec=[]
        t_ret=[]
        ang_inplane=[]
        ang_outplane=[]
        for t in orbit.t:
            try:
                s=send(t)
                r=receive(t)

                s_vec.append(s)        
                r_vec.append(r)
                ang_inplane.append(angle(s_vec[-1][0:2],-r_vec[-1][0:2]))
                t_ret.append(t)
            except ValueError:
                pass

        return [s_vec,r_vec,t_ret],[ang_inplane,ang_outplane]

    def L_PAA(func_prev,func_self,func_next,orbit,lisa,i):
        t_inter = orbit.t
        tp_guess=lisa.armlength(i,t_inter[0])/c
        tp_sol=[]
        for t in t_inter:
            L_rec = lambda tp: np.linalg.norm(func_self(t)-func_prev(t-tp))
            L_rec_val = c*t
            tp_sol.append(fsolve(L_rec - L_rec_val,tp_guess))

        return tp_sol





        


    def send_res(lisa,orbit,i,delay=True):
        
        if delay==True:
            L_out=func_arm_sec(lisa,orbit,i+1)[0]
            L_in=func_arm_sec(lisa,orbit,i+2)[0]
        else:
            L_out=lambda t: 0
            L_in=lambda t: 0
        
        pos_next=func_pos(orbit,i+2)
        pos_prev=func_pos(orbit,i+1)
        pos_self=func_pos(orbit,i)

        tp_sol = L_PAA(pos_prev,pos_self,pos_next,orbit,lisa,i)

        send=lambda t: pos_next(t+L_out(t))-pos_self(t)
        receive= lambda t: pos_self(t)-pos_prev(t-L_in(t))
        rec_dif= lambda t: pos_prev(t) - pos_prev(t-L_in(t))
        return send_res2ang(send,receive,orbit),rec_dif

    send_rec=[]
    for i in range(1,4):
        send_rec.append(send_res(lisa_cache,Orbit,1)[0][0])
        [send,rec,t_calc]=send_res(lisa_cache,Orbit,1)[0][0]


    def in_out_plane(v,n):
        inplane = v - (np.dot(v,n)/(np.linalg.norm(n)**2))*n
        outplane =(np.dot(v,n)/(np.linalg.norm(n)**2))*n

        #ang_out=math.acos(np.dot(v,n)/(np.linalg.norm(v)*np.linalg.norm(n)))
        try:
            #ang_out = np.linalg.norm(np.cross(inplane,outplane))/(np.linalg.norm(inplane)*np.linalg.norm(outplane))
            ang_out = (np.linalg.norm(outplane)/np.linalg.norm(inplane))
            ang_out=math.atan(ang_out)*np.sign(ang_out)
        except RuntimeWarning:
            ang_out = 0
        
        return [inplane,outplane,ang_out]

    def ang_2vectors(v1,v2):        
        try:
            v1=v1/np.linalg.norm(v1)
            v2=v2/np.linalg.norm(v2)
            sintheta=np.linalg.norm(np.cross(v1,v2))
        except RuntimeWarning or ValueError:
            #print(v1)
            #print(v2)
            sintheta=0
        return (math.asin(sintheta))

    def angles_calc(lisa,orbit,delay=True,range_t='all'):
        ang_all=[]
        ang_in=[]
        ang_out=[]
        ang_send_out=[]
        ang_rec_out=[]
        t_all=[]
        for j in range(1,4):
            [send,rec,t_calc]=send_res(lisa,orbit,j,delay=delay)[0][0]
            ang_all_sc=[]
            ang_in_sc=[]
            ang_out_sc=[]
            ang_send_out_sc=[]
            ang_rec_out_sc=[]
            for i in range(0,len(t_calc)):
                out=send[i]
                inc=-rec[i]
                
                normal=np.cross(out,inc)
                normal=normal/np.linalg.norm(normal)

                [out_in,out_out,out_ang]=in_out_plane(out,normal)
                [inc_in,inc_out,inc_ang]=in_out_plane(inc,normal)

                ang_all_sc.append(ang_2vectors(out,inc))
                ang_in_sc.append(ang_2vectors(out_in,inc_in))
                ang_out_sc.append(ang_2vectors(out_out,inc_out))
                ang_send_out_sc.append(inc_ang)
                ang_rec_out_sc.append(out_ang)
            ang_all.append(np.array(ang_all_sc))
            ang_in.append(np.array(ang_in_sc))
            ang_out.append(np.array(ang_out_sc))
            ang_send_out.append(np.array(ang_send_out_sc))
            ang_rec_out.append(np.array(ang_rec_out_sc))
            
            t_all.append(t_calc)

        return ang_all,ang_in,ang_out,ang_send_out,ang_rec_out,t_all

    ang_all,ang_in,ang_out,ang_send_out,ang_rec_out,t_all=angles_calc(lisa_cache,Orbit)
    ang_all_nd,ang_in_nd,ang_out_nd,ang_send_out_nd,ang_rec_out_nd,t_all_nd=angles_calc(lisa_cache,Orbit,delay=False)


    f,ax = plt.subplots(4,1)
    ax[0].axhline(math.radians(60),color='0',linestyle='--')
    ax[1].axhline(math.radians(60),color='0',linestyle='--')

    for i in range(0,len(ang_all)):
        ax[0].plot(t_all[i],ang_all[i],label='Spacecraft '+str(i+1))
        ax[1].plot(t_all[i],ang_in[i],label='Spacecraft '+str(i+1))
        #ax[2].plot(t_all[i],ang_out[i],label='Spacecraft '+str(i+1))
        ax[2].plot(t_all[i],ang_send_out[i],label='Spacecraft '+str(i+1))
        ax[3].plot(t_all[i],ang_rec_out[i],label='Spacecraft '+str(i+1))
    ax[0].legend(loc='best')
    ax[1].legend(loc='best')
    ax[2].legend(loc='best')
    ax[3].legend(loc='best')
    #ax[4].legend(loc='best')
    ax[0].set_title('Total angle')
    ax[1].set_title('Inplane angle')
    ax[2].set_title('Outplane angle send')
    ax[3].set_title('Outplane angle receive angle')
    #ax[4].set_title('Outplane angle')

    f.savefig(dir_savefig+filename_save+'-Angle_with_delay.png')

    f,ax = plt.subplots(3,1)
    for i in range(0,len(ang_all)):
        ax[0].plot(t_all[i],ang_all[i]-ang_all_nd[i][0:len(ang_all[i])],label='Spacecraft '+str(i+1))
        ax[1].plot(t_all[i],ang_in[i]-ang_in_nd[i][0:len(ang_in[i])],label='Spacecraft '+str(i+1))
        ax[2].plot(t_all[i],ang_out[i]-ang_out_nd[i][0:len(ang_out[i])],label='Spacecraft '+str(i+1))
    ax[0].legend(loc='best')
    ax[1].legend(loc='best')
    ax[2].legend(loc='best')
    ax[0].set_title('Total angle, difference delay')
    ax[1].set_title('Inplane angle, difference delay')
    ax[2].set_title('Outplane angle, difference delay')

    f.savefig(dir_savefig+filename_save+'-Angle_difference_delay.png')


    f,ax = plt.subplots(3,1)
    for i in range(0,len(ang_all)):
        ax[0].plot(t_all[i],ang_all_nd[i][0:len(ang_all[i])],label='Spacecraft '+str(i+1))
        ax[1].plot(t_all[i],ang_in_nd[i][0:len(ang_in[i])],label='Spacecraft '+str(i+1))
        ax[2].plot(t_all[i],ang_out_nd[i][0:len(ang_out[i])],label='Spacecraft '+str(i+1))
    ax[0].legend(loc='best')
    ax[1].legend(loc='best')
    ax[2].legend(loc='best')
    ax[0].set_title('Total angle, no delay')
    ax[1].set_title('Inplane angle, no delay')
    ax[2].set_title('Outplane angle, no delay')

    f.savefig(dir_savefig+filename_save+'-Angle_no_delay.png')

    if plot_on==True:
        plt.show()


dir_orbits='/home/ester/git/ynthlisa/orbits/'

filename_list=[]

for (dirpath, dirnames, filenames) in os.walk(dir_orbits):
    if filenames[-1].split('.')[-1]=='txt':
        filename_list.append(dirpath+'/'+filenames[-1]) 

filename_list=[filename_list[-1]]
for filename in filename_list:
    PAA(home,filename,False,200,plot_on=False)


