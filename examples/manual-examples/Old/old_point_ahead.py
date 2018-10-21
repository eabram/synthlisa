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


def PAA(home,filename,directory_imp,read_max='all',num_back=0,plot_on=True,scale=1000,dir_savefig=False,dir_extr='',new_folder=True,timeunit='days'):

#    
#    def get_date():
#        now = datetime.datetime.now()
#        date=str(now.year)+str(now.month)+str(now.day)
#        
#        date=date+'-'+dir_extr
#        return date
#    
#    date_str=get_date()
#    filename_save=filename.split('.')[-2]
#    filename_save=filename_save.split('/')[-1]
#    if dir_savefig==False:
#        dir_savefig=os.path.dirname(os.path.realpath(__file__))+'/figures/'+date_str+'/'
#
#    if not os.path.exists(dir_savefig):
#        os.makedirs(dir_savefig)
#    elif new_folder==True:
#        i=1
#        while os.path.exists(dir_savefig[0:-1]+'_('+str(i)+')/'):
#            i=i+1
#        dir_savefig = dir_savefig[0:-1]+'_('+str(i)+')/'
#        os.makedirs(dir_savefig)


    print('Importing Orbit')
    tic=time.clock()
    Orbit=orbit(home=home,filename=filename,directory_imp=directory_imp,num_back=num_back,scale=scale,read_max=read_max,plot_on=False,timeunit=timeunit)
    print(str(Orbit.linecount)+' datapoints')
    print('Done in '+str(time.clock()-tic))

    def angle(v1, v2):
        return (np.arcsin(np.linalg.norm(np.cross(v1,v2)/(np.linalg.norm(v1)*np.linalg.norm(v2)))))

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
                #ang_inplane.append(angle(s_vec[-1][0:2],-r_vec[-1][0:2]))
                t_ret.append(t)
            except ValueError:
                pass

        return [s_vec,r_vec,t_ret],[ang_inplane,ang_outplane]

    def L_PAA(orbit,pos_prev,pos_self,pos_next):
        t_inter = orbit.t
        tp_guess=3#lisa_cache.armlength(i,t_inter[0])/c
        tp_sol=[]
        tn_sol=[]
        t_rec_sol=[]
        t_send_sol=[]
        
        for t in t_inter:
            try:
                L_rec = lambda dtp: (np.linalg.norm(pos_self(t)-pos_prev(t-dtp)) - c*(dtp))/c # L_rec in seconds
                tp_sol.append(fsolve(L_rec,tp_guess)[0])
                t_rec_sol.append(t)
        
        #print(tp_sol)
            except ValueError:
                pass

            try:
                L_send = lambda dtn: (np.linalg.norm(pos_next(t+dtn)-pos_self(t)) - c*(dtn))/c # L_rec in seconds
                tn_sol.append(fsolve(L_send,tp_guess)[0])
                t_send_sol.append(t)

            except ValueError:
                pass

        #print(len(t_rec_sol),len(tp_sol))
        L_in=interp1d(t_rec_sol,tp_sol)
        L_out=interp1d(t_send_sol,tn_sol)
        
        return L_in, L_out

    def send_res(lisa,orbit,i,delay=True,around_clock=False):
        
        if around_clock==False:
            pos_next=func_pos(orbit,i+2)
            pos_prev=func_pos(orbit,i+1)
            pos_self=func_pos(orbit,i)

            if delay == True:
                L_in,L_out=L_PAA(orbit,pos_prev,pos_self,pos_next)
            elif delay == 'Not ahead':
                L_out=func_arm_sec(lisa,orbit,i+1)[0]
                L_in=func_arm_sec(lisa,orbit,i+2)[0]
            else: 
                L_out=lambda t: 0
                L_in=lambda t: 0

            #tp_sol = L_PAA(pos_prev,pos_self,pos_next,orbit,lisa,i)
            tp_sol=0
            send=lambda t: pos_next(t+L_out(t))-pos_self(t)
            receive= lambda t: pos_self(t)-pos_prev(t-L_in(t))
            rec_dif= lambda t: pos_prev(t) - pos_prev(t-L_in(t))
            return send_res2ang(send,receive,orbit),rec_dif,tp_sol

        else:
            pos_next=func_pos(orbit,i+1)
            pos_prev=func_pos(orbit,i+2)
            pos_self=func_pos(orbit,i)

            if delay == True:
                L_in,L_out=L_PAA(orbit,pos_prev,pos_self,pos_next)
            elif delay == 'Not ahead':
                L_out=func_arm_sec(lisa,orbit,i+1)[0]
                L_in=func_arm_sec(lisa,orbit,i+2)[0]
            else: 
                L_out=lambda t: 0
                L_in=lambda t: 0

            #tp_sol = L_PAA(pos_prev,pos_self,pos_next,orbit,lisa,i)
            tp_sol=0
            send=lambda t: pos_next(t+L_out(t))-pos_self(t)
            receive= lambda t: pos_self(t)-pos_prev(t-L_in(t))
            rec_dif= lambda t: pos_prev(t) - pos_prev(t-L_in(t))
            return send_res2ang(send,receive,orbit),rec_dif,tp_sol




    #send_rec=[]
    #for i in range(1,4):
    #    send_rec.append(send_res(lisa_cache,Orbit,i,around_clock=True)[0][0])
    #    [send,rec,t_calc]=send_res(lisa_cache,Orbit,i)[0][0]


    def in_out_plane(v,n):
        inplane = v - (np.dot(v,n)/(np.linalg.norm(n)**2))*n
        outplane =(np.dot(v,n)/(np.linalg.norm(n)**2))*n

        #ang_out=math.acos(np.dot(v,n)/(np.linalg.norm(v)*np.linalg.norm(n)))
        try:
            #ang_out = np.linalg.norm(np.cross(inplane,outplane))/(np.linalg.norm(inplane)*np.linalg.norm(outplane))
            #ang_out = np.cross(outplane,n)
            #ang_norm = np.linalg.norm(ang_out)
            #ang_out = ang_out/np.linalg.norm(ang_out)
            #ang_out = np.arcsin(ang_out)
            ang_out = (np.linalg.norm(outplane)/np.linalg.norm(inplane))
            ang_out=np.arctan(ang_out)*np.sign(np.dot(outplane,n))
        except RuntimeWarning:
            ang_out = 0
        
        return [inplane,outplane,ang_out]

    def ang_2vectors(v1,v2): 
        v1_norm=np.linalg.norm(v1)
        v2_norm=np.linalg.norm(v2)
        if v1_norm!=0 and v2_norm!=0:
            try:
                #v1=v1/np.linalg.norm(v1)
                #v2=v2/np.linalg.norm(v2)
                sintheta=np.cross(v1,v2)/(np.linalg.norm(v1)*np.linalg.norm(v2))
                sintheta=np.linalg.norm(sintheta)
            except RuntimeWarning or ValueError:
                print(v1)
                print(v2)
                sintheta=0
        else:
            sintheta=0
        return (np.arcsin(sintheta))

    def normal_pos(orbit,i):
        ret=[]
        l_inc = (i+2)%3 - 1
        l_out = (i+1)%3 - 1
        
        for j in range(0,len(orbit.t)):
            L_inc = orbit.p[l_inc][j,:]
            L_out = orbit.p[l_out][j,:]

            n_inc=L_inc/np.linalg.norm(L_inc)
            n_out=L_out/np.linalg.norm(L_out)

            ret.append(np.cross(n_inc,n_out))

        ret=np.array(ret)

        return ret

    def angles_calc(lisa,orbit,delay=True,range_t='all',normal_set='sc',around_clock=False):
        ang_all=[]
        ang_in=[]
        ang_out=[]
        ang_send_out=[]
        ang_rec_out=[]
        t_all=[]
        wobbling_ang_list=[]
        out_list=[]
        inc_list=[]
        for j in range(1,4):
            [send,rec,t_calc]=send_res(lisa,orbit,j,delay=delay,around_clock=around_clock)[0][0]
            ang_all_sc=[]
            ang_in_sc=[]
            ang_out_sc=[]
            ang_send_out_sc=[]
            ang_rec_out_sc=[]
            normal_all = orbit.n[j-1] #normal_pos(orbit,j)
            wobbling_ang=[]
            out_sc=[]
            inc_sc=[]
            for i in range(0,len(t_calc)):
                out=send[i]
                inc=-rec[i]
                out_sc.append(out)
                inc_sc.append(inc)
                
                # Normal respect to laser beams
                normal_beam=np.cross(out,inc)
                normal_beam=normal_beam/np.linalg.norm(normal_beam)
                # Normal with respect to positions spacecrafts
                normal_sc = normal_all[i,:]
                if normal_set=='beam':
                    normal = normal_beam
                elif normal_set == 'sc':
                    normal = normal_sc
                
                wobbling_ang.append(ang_2vectors(normal_beam,normal_sc))


                [out_in,out_out,out_ang]=in_out_plane(out,normal)
                [inc_in,inc_out,inc_ang]=in_out_plane(inc,normal)

                ang_all_sc.append(ang_2vectors(out,inc))
                ang_in_sc.append(ang_2vectors(out_in,inc_in))
                ang_out_sc.append(ang_2vectors(out_out,inc_out))
                ang_send_out_sc.append(inc_ang)
                ang_rec_out_sc.append(out_ang)
            
                # Angle with COM
                #COM=orbit.COM[i]

                #in_out_plane(out
            
            
            
            ang_all.append(np.array(ang_all_sc))
            ang_in.append(np.array(ang_in_sc))
            ang_out.append(np.array(ang_out_sc))
            ang_send_out.append(np.array(ang_send_out_sc))
            ang_rec_out.append(np.array(ang_rec_out_sc))
            wobbling_ang_list.append(np.array(wobbling_ang))
            out_list.append(np.array(out_sc))
            inc_list.append(np.array(inc_sc))
            t_all.append(t_calc)

        return [[ang_all,ang_in,ang_out,ang_send_out,ang_rec_out,t_all],[wobbling_ang_list,inc_list,out_list]]



    [[ang_all,ang_in,ang_out,ang_send_out,ang_rec_out,t_all],[wobbling,inc_l,out_r]]=angles_calc(lisa_cache,Orbit)
    [[ang_all_nd,ang_in_nd,ang_out_nd,ang_send_out_nd,ang_rec_out_nd,t_all_nd],[wobbling_nd,inc_nd_l,out_nd_r]]=angles_calc(lisa_cache,Orbit,delay='Not ahead')

    [[ang_all_c,ang_in_c,ang_out_c,ang_send_out_c,ang_rec_out_c,t_all_c],[wobbling_c,inc_r,out_l]]=angles_calc(lisa_cache,Orbit,around_clock=True)
    [[ang_all_nd_c,ang_in_nd_c,ang_out_nd_c,ang_send_out_nd_c,ang_rec_out_nd_c,t_all_nd_c],[wobbling_nd_c,inc_nd_r,out_nd_l]]=angles_calc(lisa_cache,Orbit,delay='Not ahead',around_clock=True)

    # Obtaining angle between left and right sending beams
    def ang_send_lr(out_l,out_r,t_vec,orbit):
        ang_list=[]
        ang_out_l_list=[]
        ang_out_r_list=[]
        ang_in_l_list=[]
        ang_in_r_list=[]

        t_list = []
        for i in range(0,len(out_l)):
            ang=[]
            ang_out_l_vec=[]
            ang_out_r_vec=[]
            ang_in_l_vec=[]
            ang_in_r_vec=[]
            t_ret=[]
            for j in range(0,len(t_vec[i])):
                try: 
                    normal = orbit.n_func[i](t_vec[i][j]) #orbit.n[i][j,:]
                    r = orbit.r_func[i](t_vec[i][j]) #orbit.COM[j,:] - orbit.p[i][j,:]
                    r = r/np.linalg.norm(r)
                    t_ret.append(t_vec[i][j])
                    v_l = orbit.v_l_func[i](t_vec[i][j])
                    v_r = orbit.v_r_func[i](t_vec[i][j])

                    v1 = out_l[i][j,:]
                    v2 = out_r[i][j,:]

                    v1=v1/np.linalg.norm(v1)
                    v2=v2/np.linalg.norm(v2)
                   
                    ang.append(ang_2vectors(v1,v2))
                    ang_in_l_vec.append(ang_2vectors(v1,r) - ang_2vectors(v_l,r))
                    ang_in_r_vec.append(ang_2vectors(v2,r) - ang_2vectors(v_r,r))

                    #calc=np.cross(v1,v2)
                    #sign_calc = np.sign(np.dot(v1,v2))
                    #calc = calc/(np.linalg.norm(v1)*np.linalg.norm(v2))
                    #ang.append(np.arcsin(sign_calc*(np.linalg.norm(calc))))
                    ##ang.append(np.arccos(np.dot(v1,v2)/(np.linalg.norm(v1)*np.linalg.norm(v2))))
                    #
                    #calc = np.cross(v1,v_l)
                    #sign_calc = 1#np.sign(np.dot(v1,r))
                    #calc = np.linalg.norm(calc)/(np.linalg.norm(v1)*np.linalg.norm(v_l))
                    ##ang_in_l_vec.append(np.arccos(np.dot(v1,r)/(np.linalg.norm(v1)*np.linalg.norm(r))))
                    #ang_in_l_vec.append(np.arcsin(sign_calc*(np.linalg.norm(calc))))
                    #
                    #calc = np.cross(v2,v_r)#v2,r)
                    #sign_calc = 1#np.sign(np.dot(v2,r))
                    #calc = np.linalg.norm(calc)/(np.linalg.norm(v2)*np.linalg.norm(v_r))
                    ##ang_in_r_vec.append(np.arccos(np.dot(v2,r)/(np.linalg.norm(v2)*np.linalg.norm(r))))
                    #ang_in_r_vec.append(np.arcsin(sign_calc*(np.linalg.norm(calc))))
                    
                    [inplane_l,outplane_l,ang_out_l]=in_out_plane(v1,normal) # Left
                    [inplane_r,outplane_r,ang_out_r]=in_out_plane(v2,normal) # Right
                    ang_out_l_vec.append(ang_out_l)
                    ang_out_r_vec.append(ang_out_r)
                except ValueError:
                    pass



            ang_list.append(np.array(ang))
            ang_out_l_list.append(np.array(ang_out_l_vec))
            ang_out_r_list.append(np.array(ang_out_r_vec))
            ang_in_l_list.append(np.array(ang_in_l_vec))
            ang_in_r_list.append(np.array(ang_in_r_vec))
            t_list.append(np.array(t_ret))
        
        PPA_in_l=[]
        PPA_in_r=[]
        for i in range(1,4):           
            PPA_in_l.append(ang_in_l_list[i-1])# - orbit.ang_stat_l[i-1](t_list[i-1]))
            PPA_in_r.append(ang_in_r_list[i-1])# - orbit.ang_stat_r[i-1](t_list[i-1]))



        return [ang_list,ang_in_l_list,ang_in_r_list,ang_out_l_list,ang_out_r_list,PPA_in_l,PPA_in_r,t_list]

    [ang_sendbeams,ang_in_l,ang_in_r,PAA_out_l,PAA_out_r,PAA_in_l,PAA_in_r,t_list]=ang_send_lr(out_l,out_r,t_all,Orbit)


    def plot_angles(t_all,ang_plus,ang_min=False,titles=False):
        [ang_all,ang_in,ang_send_out,ang_rec_out]=ang_plus
        
        if ang_min==False:
            ang_all=ang_plus
        else:
            ang_all=ang_plus-ang_min
        
        f,ax = plt.subplots(4,1)
        plt.subplots_adjust(hspace=1)
        ax[0].axhline(math.radians(60),color='0',linestyle='--')
        ax[1].axhline(math.radians(60),color='0',linestyle='--')

        for i in range(0,len(ang_all)):
            ax[0].plot(t_all[i],ang_all[i],label='Spacecraft '+str(i+1))
            ax[1].plot(t_all[i],ang_in[i],label='Spacecraft '+str(i+1))
            #ax[2].plot(t_all[i],ang_out[i],label='Spacecraft '+str(i+1))
            ax[2].plot(t_all[i],ang_send_out[i],label='Spacecraft '+str(i+1))
            ax[3].plot(t_all[i],ang_rec_out[i],label='Spacecraft '+str(i+1))
        for j in range(0,len(ax)):
            ax[j].legend(loc='best')
            ax[j].set_xlabel('Time (sec)')
            ax[j].set_ylabel('Angle (rad)')
            if titles!=False:
                ax[j].set_title(titles[j])
        plt.close()

        return f,ax

    titles=['Total angle','Inplane angle','Outplane angle send','Outplane angle receive angle']

    #plot_angles(t_all,[np.array(ang_all),np.array(ang_in),np.array(ang_send_out),np.array(ang_rec_out)],titles=titles)
     


    figs=[]
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


        titles=['Wobbling angle of spacecraft - With delay','Wobbling angle of spacecraft - No delay','Difference in wobbling angle']

        f,ax = plt.subplots(3,1,figsize=(15,15))
        for i in range(0,len(wobbling)):
            ax[0].plot(t_all[i],wobbling[i],label='Spacecraft '+str(i+1))
            ax[1].plot(t_all[i],wobbling_nd[i],label='Spacecraft '+str(i+1))
            ax[2].plot(t_all[i],wobbling[i]-wobbling_nd[i],label='Spacecraft '+str(i+1))
        for j in range(0,len(ax)):
            ax[j].set_title(titles[j])
            ax[j].set_xlabel('Time (sec)')
            ax[j].set_ylabel('Angle (rad)')
            ax[j].legend(loc='best')
           
        figs.append(f)
        #f.savefig(dir_savefig+filename_save+'-Wobbling_angles.png')
        plt.close()


        f,ax = plt.subplots(4,1,figsize=(15,15))
        plt.subplots_adjust(hspace=1)
        ax[0].axhline(math.radians(60),color='0',linestyle='--')
        ax[1].axhline(math.radians(60),color='0',linestyle='--')

        for i in range(0,len(ang_all)):
            ax[0].plot(t_all[i],ang_all[i],label='Spacecraft '+str(i+1))
            ax[1].plot(t_all[i],ang_in[i],label='Spacecraft '+str(i+1))
            #ax[2].plot(t_all[i],ang_out[i],label='Spacecraft '+str(i+1))
            ax[2].plot(t_all[i],ang_send_out[i],label='Spacecraft '+str(i+1))
            ax[3].plot(t_all[i],ang_rec_out[i],label='Spacecraft '+str(i+1))
        for j in range(0,len(ax)):
            ax[j].legend(loc='best')
            ax[j].set_xlabel('Time (sec)')
            ax[j].set_ylabel('Angle (rad)')

        ax[0].set_title('Total angle')
        ax[1].set_title('Inplane angle')
        ax[2].set_title('Outplane angle send')
        ax[3].set_title('Outplane angle receive angle')

        figs.append(f)
        #f.savefig(dir_savefig+filename_save+'-Angle_with_delay.png')
        plt.close()

        f,ax = plt.subplots(3,1,figsize=(15,15))

        PAA_all={}
        PAA_in={}
        PAA_out={}
        PAA_t={}
        for i in range(0,len(ang_all)):
            PAA_all[str(i+1)]=ang_all[i]-ang_all_nd[i][0:len(ang_all[i])]
            PAA_in[str(i+1)]=ang_in[i]-ang_in_nd[i][0:len(ang_in[i])]
            PAA_out[str(i+1)]=ang_out[i]-ang_out_nd[i][0:len(ang_out[i])]
            PAA_t[str(i+1)]=t_all[i]

            ax[0].plot(t_all[i],ang_all[i]-ang_all_nd[i][0:len(ang_all[i])],label='Spacecraft '+str(i+1))
            ax[1].plot(t_all[i],ang_in[i]-ang_in_nd[i][0:len(ang_in[i])],label='Spacecraft '+str(i+1))
            ax[2].plot(t_all[i],ang_out[i]-ang_out_nd[i][0:len(ang_out[i])],label='Spacecraft '+str(i+1))
        for j in range(0,len(ax)):
            ax[j].legend(loc='best')
            ax[j].set_xlabel('Time (sec)')
            ax[j].set_ylabel('Angle (rad)')
        ax[0].set_title('Total angle, difference delay')
        ax[1].set_title('Inplane angle, difference delay')
        ax[2].set_title('Outplane angle, difference delay')
        plt.close()
        figs.append(f)


        f,ax = plt.subplots(3,1,figsize=(15,15))
        for i in range(0,len(ang_all)):
            ax[0].plot(t_all[i],ang_all_nd[i][0:len(ang_all[i])],label='Spacecraft '+str(i+1))
            ax[1].plot(t_all[i],ang_in_nd[i][0:len(ang_in[i])],label='Spacecraft '+str(i+1))
            ax[2].plot(t_all[i],ang_out_nd[i][0:len(ang_out[i])],label='Spacecraft '+str(i+1))

        for j in range(0,len(ax)):
            ax[j].legend(loc='best')
            ax[j].set_xlabel('Time (sec)')
            ax[j].set_ylabel('Angle (rad)')

        ax[0].set_title('Total angle, no delay')
        ax[1].set_title('Inplane angle, no delay')
        ax[2].set_title('Outplane angle, no delay')
        plt.close()

        figs.append(f)

    PAA_dict=PAA_in
    keys=PAA_dict.keys()
    f,ax = plt.subplots(2,len(keys)+1,figsize=(15,15))
    for i in range(0,len(keys)):
        y=PAA_dict[keys[i]]
        x=PAA_t[keys[i]]

        sc='SC'+str(i+1)
        #func=interp1d(x,y)

        dfdt = np.diff(y)/np.diff(x)
        #dfdt=np.concatenate(dfdt,np.array([(y[-1] - y[-2])/(x[-1] - x[-2])]))
        AX_tot=ax[0,len(keys)],ax[1,len(keys)]
        AX_tot[0].plot(x,y,label='SC'+str(i+1))
        AX_tot[1].plot(x[0:-1],dfdt,label='SC'+str(i+1))

        AX=ax[0,i],ax[1,i]
        AX[0].plot(x,y)
        AX[1].plot(x[0:-1],dfdt)
        AX[0].set_title('In plane angle of '+sc)
        AX[0].set_ylabel('Angle (rad)')
        AX[1].set_title('In plane derivative angle of '+sc)
        AX[1].set_xlabel('Time (sec)')
        AX[1].set_ylabel('dAngle/dt (rad/sec}')

    AX_tot[0].set_title('PAA all spacecrafts')
    AX_tot[0].set_ylabel('Angle (rad)')
    AX_tot[1].set_xlabel('Time (sec)')
    AX_tot[1].set_ylabel('dAngle/dt (rad/sec)')
    plt.close()
    figs.append(f)   


    PAA=[PAA_in_l,PAA_in_r,PAA_out_l,PAA_out_r]
    f,ax = plt.subplots(4,3,figsize=(40,15))
    f.suptitle('Point ahead angle')
    plt.subplots_adjust(hspace=0.6,wspace=0.6)
    titles=['In plane left','In plane right','Out of plane left','Out of plane right']
    for i in range(0,len(ax)):
        for j in range(0,len(PAA[i])):
            x=t_list[j]
            y=PAA[i][j]
            x_adj=x/day2sec
            ax[i,j].plot(x_adj,PAA[i][j],label='S'+str(j+1))
            ax[i,j].set_title(titles[i])
            ax[i,j].set_xlabel('Time (days)')
            ax[i,j].set_ylabel('Angle (rad)')
            ax[i,j].legend(loc='best')
    figs.append(f)

    f,ax = plt.subplots(4,3,figsize=(40,15))
    f.suptitle('Point ahead angle per time')
    plt.subplots_adjust(hspace=0.6,wspace=0.6)
    titles=['In plane left','In plane right','Out of plane left','Out of plane right']
    for i in range(0,len(ax)):
        for j in range(0,len(PAA[i])):
            x=t_list[j]
            y=PAA[i][j]
            dfdx = np.diff(y)/np.diff(x)
            x_calc=x/day2sec
            x_adj = x_calc[0:-1]
            ax[i,j].plot(x_adj,dfdx,label='S'+str(j+1))
            ax[i,j].set_title(titles[i])
            ax[i,j].set_xlabel('Time (days)')

            ax[i,j].set_ylabel('Angle/time (rad/sec)')
            ax[i,j].legend(loc='best')
    figs.append(f)




    def save_fig(figs):
        titles=['-Wobbling_angles','-Angle_with_delay','-Angle_difference_delay','-Angle_no_delay','-PPA','-newPAA','-diffPPA']
        for i in range(0,len(figs)):
            title=filename_save+titles[i]+'.png'
            figs[i].savefig(dir_savefig+title)
            
            print('Figure '+title+' saved in:')
            print(dir_savefig)

    save_fig(figs)
    
    print('')
    print('')
    
    
    plt.close()

    return [[figs,Orbit,lisa_cache],[inc_l,inc_r,out_l,out_r,dir_extr_new]]

dir_orbits='/home/ester/git/synthlisa/orbits/'


filename_list=[]

for (dirpath, dirnames, filenames) in os.walk(dir_orbits):
    print(filenames)
    for i in filenames:
        if i.split('.')[-1]=='txt':
            filename_list.append(dirpath+i)

#filename_list=[filename_list[1]]
#timeunit=['seconds']
dir_extr='zzz'
timeunit=['seconds','days','days']

count=0
for i in filename_list:
    if i == filename_list[0]:
        new_folder=False # Adjust if you (don't) want to override
    else:
        new_folder=False
    print('Dir_extr:'+dir_extr)
    [[figs,Orbit,Lisa],[inc_l,inc_r,out_l,out_r,dir_extr]] = PAA(home,i,False,200,plot_on=True,dir_extr=dir_extr,new_folder=new_folder,timeunit=timeunit[count])
    count=count+1
#ang_list=[]
#for i in range(0,len(out_l)):
#    ang=[]
#    for j in range(0,len(out_l[i])):
#        v1 = out_l[i][j,:]
#        v2 = out_r[i][j,:]
#        
#        v1=v1/np.linalg.norm(v1)
#        v2=v2/np.linalg.norm(v2)
#        calc = np.cross(v1,v2)/(np.linalg.norm(v1)*np.linalg.norm(v2))
#        ang.append((np.linalg.norm(calc)))
#    ang_list.append(np.array(ang))
#    
#    plt.plot(ang-Orbit.ang_in[i][1:-1])




#x1 = Orbit.p[0][:,0]
#y1 = Orbit.p[0][:,1]
#x2 = Orbit.p[1][:,0]
#y2 = Orbit.p[1][:,1]
#x3 = Orbit.p[2][:,0]
#y3 = Orbit.p[2][:,1]


plt.close()
#plt.figure()
#plt.plot(x1,y1)
#plt.plot(x2,y2)
#plt.plot(x3,y3)
#plt.show()






#plt.show()



