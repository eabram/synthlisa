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

year2sec=32536000
c=300000000
noise_check=False
#noise_check=True

#home='/home/ester/Dropbox/Master_Project/Synthetic_LISA/'
#filename='NGO_1M_10deg_synthlisa.txt'
#directory_imp='LISAorbits/NGO_1M_10deg/'
home='/home/ester/git/synthlisa/'
directory_imp='lisasim/data/'
filename='positions.txt'
num_back=0
print('Importing Orbit')
tic=time.clock()
Orbit=orbit(home=home,filename=filename,directory_imp=directory_imp,num_back=num_back,scale=1000,read_max=200,plot_on=False)
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

    ang_out=math.acos(np.dot(v,n)/(np.linalg.norm(v)*np.linalg.norm(n)))
    
    return [inplane,outplane]

def ang_2vectors(v1,v2): 
    sintheta=np.linalg.norm(np.cross(v1,v2))/(np.linalg.norm(v1)*np.linalg.norm(v2))
    
    return (math.asin(sintheta))

def angles_calc(lisa,orbit,delay=True,range_t='all'):
    ang_all=[]
    ang_in=[]
    ang_out=[]
    t_all=[]
    for j in range(1,4):
        [send,rec,t_calc]=send_res(lisa,orbit,j,delay=delay)[0][0]
        ang_all_sc=[]
        ang_in_sc=[]
        ang_out_sc=[]
        for i in range(0,len(t_calc)):
            out=send[i]
            inc=-rec[i]
            
            normal=np.cross(out,inc)
            normal=normal/np.linalg.norm(normal)

            [out_in,out_out]=in_out_plane(out,normal)
            [inc_in,inc_out]=in_out_plane(inc,normal)

            ang_all_sc.append(ang_2vectors(out,inc))
            ang_in_sc.append(ang_2vectors(out_in,inc_in))
            ang_out_sc.append(ang_2vectors(out_out,inc_out))
        ang_all.append(np.array(ang_all_sc))
        ang_in.append(np.array(ang_in_sc))
        ang_out.append(np.array(ang_out_sc))
        t_all.append(t_calc)

    return ang_all,ang_in,ang_out,t_all

ang_all,ang_in,ang_out,t_all=angles_calc(lisa_cache,Orbit)
ang_all_nd,ang_in_nd,ang_out_nd,t_all_nd=angles_calc(lisa_cache,Orbit,delay=False)


f,ax = plt.subplots(3,1)
ax[0].axhline(math.radians(60),color='0',linestyle='--')
ax[1].axhline(math.radians(60),color='0',linestyle='--')

for i in range(0,len(ang_all)):
    ax[0].plot(t_all[i],ang_all[i],label='Spacecraft '+str(i+1))
    ax[1].plot(t_all[i],ang_in[i],label='Spacecraft '+str(i+1))
    ax[2].plot(t_all[i],ang_out[i],label='Spacecraft '+str(i+1))
ax[0].legend(loc='best')
ax[1].legend(loc='best')
ax[2].legend(loc='best')
ax[0].set_title('Total angle')
ax[1].set_title('Inplane angle')
ax[2].set_title('Outplane angle')


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
plt.show()



























x1=Orbit.p[0][:,0]
y1=Orbit.p[0][:,1]
x2=Orbit.p[1][:,0]
y2=Orbit.p[1][:,0]
x3=Orbit.p[2][:,0]
y3=Orbit.p[2][:,0]

ang_test_z=[]
for i in range(0,len(x1)):
    v1=[x2[i],y2[i]]
    v2=[x3[i],y3[i]]
    ang_test_z.append(math.degrees(math.atan((v2[1]-v1[1])/(v2[0]-v1[0]))))





def point_ahead(lisa,orbit,i,delay=True,single_point=False):
    t_vec=orbit.t

    i=i-1
    i_res=((i+2)%3)+1
    i_send=(i)+1
    i_arm=((i+1)%3)+1
    
    
    f_pos_res=func_pos(orbit,i_res)
    f_pos_send=func_pos(orbit,i_send)
    f_L,L_max=func_arm_sec(lisa,orbit,i_arm)

    t_ret=[]
    for j in t_vec:
        if j+L_max<=t_vec[-1]:
            t_ret.append(j)

    if delay==True: # Vector from
        f_point=lambda t: f_pos_res(t+f_L(t))-f_pos_send(t)
    elif delay==False:
        f_point=lambda t: f_pos_res(t)-f_pos_send(t)

    if single_point==False:
        return f_point,t_ret
    else:
        try:
            f_point(single_point)
            return f_point,single_point
        except:
            print("Please select other timestamp")
            return 0

print("")

L=[]
L_read=[]
t_stop=len(t)
for j in range(0,t_stop):
    if j%20 ==0:
        print(str(j))
    L.append(lisa_orb.armlength(1,t[j]))
    L_read.append(np.linalg.norm(Orbit.L[0][j]))

plot_on=False
if plot_on==True:
    plt.plot(t[0:t_stop],L,label='LISA')
    plt.plot(t[0:t_stop],L_read,label='Orbit')
    plt.legend(loc='best')
    plt.show()


def center_of_mass_vec(orbit,i): #...adjust for difference in feight spacecrafts
    i=i-1
    t_vec=orbit.t
    orbit_center=(orbit.p[0]+orbit.p[1]+orbit.p[2])/3.0
    orbit_center=-orbit.p[i]+orbit_center

    [fx,fy,fz]=[interp1d(t_vec,orbit_center[:,0]),interp1d(t_vec,orbit_center[:,1]),interp1d(t_vec,orbit_center[:,2])]

    return lambda t:np.array([fx(t),fy(t),fz(t)])



pos_func_on=False

if pos_func_on==True:
    print('Obtaining position functions of LISA object')
    tic=time.clock()
    f_del=[]
    f_nodel=[]
    orb_center=[]

    for i in range(1,4):
        f_del.append(point_ahead(lisa_cache,Orbit,i))
        print(str(round(((i*2)/6)*100,0))+'%')
        f_nodel.append(point_ahead(lisa_cache,Orbit,i,delay=False))
        print(str(round(((i*2+1)/6)*100,333))+'%')
        orb_center.append(center_of_mass_vec(Orbit,i))
    print('Done in '+str(time.clock()-tic))

    def angle(v1, v2):
        return math.degrees(math.asin(np.cross(v1,v2)/(np.linalg.norm(v1)*np.linalg.norm(v2))))
        
        #return math.degrees(math.acos(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))))

    def angle_calc(f,i):
        i=i-1
        t_vec=f[i][1]
        inplane=[]
        outplane=[]
        tot_angle=[]
        #orb_center=center_of_mass_vec(Orbit,t_vec,i+1)
        for t in t_vec:
            out=np.array(f[i][0](t))
            #inc=np.array(orb_center[i](t))
            inc=f[(i+1)%3][0](t) # t --> t-L
            

            inplane.append(angle(out[0:2],inc[0:2]))
            outplane.append(math.degrees(math.asin(out[2]/np.linalg.norm(out))))
            tot_angle.append(0)
        return [[inplane,outplane],tot_angle]

    print("")
    print("Obtaining the point ahead angles")
    tic=time.clock()
    [q1_del,tot_1]=angle_calc(f_del,1)
    [q2_del,tot_2]=angle_calc(f_del,2)
    [q3_del,tot_3]=angle_calc(f_del,3)

    [q1_nodel,tot_no1]=angle_calc(f_nodel,1)
    [q2_nodel,tot_no2]=angle_calc(f_nodel,2)
    [q3_nodel,tot_no3]=angle_calc(f_nodel,3)
    print('Done in '+str(time.clock()-tic))

    print("")
    print("Plotting angles")
    q_del=[q1_del,q2_del,q3_del]
    q_nodel=[q1_nodel,q2_nodel,q3_nodel]

    #x=np.array(f_del[0][1])/year2sec
    x=range(0,len(q_del[0][0]))
    f,axarr = plt.subplots(3,2)
    for i in range(0,len(q_del)):
        axarr[i,0].plot(f_del[i][1],q_del[i][0],'b-',label="Delayed")
        axarr[i,0].plot(f_nodel[i][1],q_nodel[i][0],'r-',label="Not delayed")
        axarr[i,0].set_title("Inplane angle spacecraft "+str(i+1))
        axarr[i,0].legend(loc='best')

        axarr[i,1].plot(f_del[i][1],q_del[i][1],'b-',label="Delayed")
        axarr[i,1].plot(f_nodel[i][1],q_nodel[i][1],'r-',label="Not delayed")
        axarr[i,1].set_title("Outplane angle spacecraft "+str(i+1))
        axarr[i,1].legend(loc='best')

    plt.show()

#Orbit.
    



