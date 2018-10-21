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

noise_check=False
#noise_check=True

home='/home/ester/Dropbox/Master_Project/Synthetic_LISA/'
filename='NGO_1M_10deg_synthlisa.txt'
directory_imp='LISAorbits/NGO_1M_10deg/'
num_back=0
Orbit=orbit(home=home,filename=filename,directory_imp=directory_imp,num_back=4)

lisa=Orbit.lisa_obj
t=Orbit.t

def nominal_arm(i,t):
    def func_arm(i):
        L_vec=[]
        t_vec=Orbit.t
        for j in range(0,len(t_vec)):
            L_vec.append(np.linalg.norm(Orbit.L[i-1][j]))
    
        f=interp1d(t_vec,L_vec)

        return f

    f=func_arm(i)

    return f(t)

func_nominal_arm = lambda i,time: nominal_arm(i,time)

lisa_orb=PyLISA(lisa,func_nominal_arm)
lisa_cache=CacheLISA(lisa_orb) # Works with retard() and putp and putn
t=np.array(Orbit.t)*365.25*24*3600

def point_ahead(lisa,orbit,i,delay=True):
    t_vec=np.array(orbit.t)*365.25*24*3600 # Check which value is used)
    
    L=[]
    #for j in range(0,len(t_vec)):
    #    L.append(lisa.armlength(i,orbit.t[j]))
    #L=np.array(L)
    #L=interp1d(t_vec,L)
    #L=lambda t: lisa.armlength(i,t/(365.25*360*24)
    #print(L)
    #print(t_vec[loc])
    #print(t_vec[loc]-L)
    i_res=(i+1)%3
    i_send=(i+2)%3
    f_pos_res=func_pos(orbit,i_res,t_vec)
    f_pos_send=func_pos(orbit,i_send,t_vec)

    f_point=[]
    t_vec_arr=[]
    for loc in range(0,len(t_vec)):
        L=lisa.armlength(i,orbit.t[loc])    
        if delay==True:
            try:
                f_point.append(f_pos_res(t_vec[loc])-f_pos_send(t_vec[loc]-L))
                t_vec_arr.append(t_vec[loc])
            except:
                pass
        elif delay==False:
            f_point.append(f_pos_res(t_vec[loc])-f_pos_send(t_vec[loc]))
            t_vec_arr.append(t_vec[loc])

    L_max=[0,0]
    for j in range(0,len(t)):
        L_new=lisa.armlength(i,orbit.t[j])
        if L_new>L_max[0]:
            L_max[0]=L_new
            L_max[1]=j


    return np.array(f_point),L_max,np.array(t_vec_arr)

def func_pos(orbit,i,t_inter):
    i=i-1
    [x,y,z]=[orbit.p[i,:,0],orbit.p[i,:,1],orbit.p[i,:,2]]
    fx=interp1d(t_inter,x)
    fy=interp1d(t_inter,y)
    fz=interp1d(t_inter,z)

    return lambda time: np.array([fx(time),fy(time),fz(time)])

def point_ahead_func(lisa,orbit,i,loc,opt1='value',delay=True):
    [ret,L_max,t_ret]=point_ahead(lisa,orbit,i,delay=delay)
    return [ret,t_ret]

    print('Done for arm '+str(i))
#pos_func_res=func_pos(Orbit,1,t)
#pos_func(

q_delay=[]
q_nodelay=[]

for i in range(1,4):
    q_delay.append(point_ahead_func(lisa_cache,Orbit,i,500,opt1='func',delay=True))
    q_nodelay.append(point_ahead_func(lisa_cache,Orbit,i,500,opt1='func',delay=False))

#qqq=point_ahead(lisa_cache,Orbit,1,500)

def angle(v1, v2):
  return math.degrees(math.acos(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))))



def ang_plane(q,i):
    time=q[0][1]
    ang_inplane=[]
    ang_outplane=[]
    for loc in range(0,len(time)):
        i=i-1
        out=q[(i+1)%3][0][loc,:]
        inc=q[(i+2)%3][0][loc,:]
        ang_inplane.append(angle(-inc[0:2],out[0:2]))
        ang_outplane.append(math.degrees(math.acos(np.linalg.norm(out[1:3])/np.linalg.norm(out))))
    
    return [ang_inplane,ang_outplane]

f2, axarr = plt.subplots(3,2)
for i in range(1,4):
    [inplane_del,outplane_del]=ang_plane(q_delay,i)
    [inplane_nodel,outplane_nodel]=ang_plane(q_nodelay,i)
    axarr[i-1,0].set_title('Point ahead angle spacecraft '+str(i))
    axarr[i-1,0].plot(inplane_del,'b-',label='Inplane delayed')
    axarr[i-1,1].plot(inplane_nodel,'r-',label='Inplane not delayed')
    axarr[i-1,0].plot(outplane_del,'b--',label='Outplane delayed')
    axarr[i-1,1].plot(outplane_nodel,'r--',label='Outplane not delayed')
    axarr[i-1,0].legend(loc='best')
    axarr[i-1,1].legend(loc='best')

plt.show()








skip=True
if skip==False:
    plt.figure()
    L_plot=[]
    L_lisa=[]
    L_lisa_orb=[]
    for i in range(0,len(t)):
        L_lisa.append(lisa.armlength(1,t[i]))
        L_lisa_orb.append(lisa_orb.armlength(1,t[i]))
        #L_plot.append(np.linalg.norm(Orbit.L[0][i]))
    plt.plot(t,L_lisa,label='SampledLISA')
    plt.plot(t,L_lisa_orb,label='PyLISA')
    plt.plot(t,L_plot,label='baseLISA')
    plt.legend(loc='best')
    plt.show()



if noise_check==True:
    proofnoise = PowerLawNoise(1.0,256.0,2.5e-48,-2.0,interplen=1)
    samples = 2**20  # 2**20 takes 13s on a 1.25GHz G4; 2*22 used for plot
    sstep = 0.1
    patches = 64
    noises_pr = getobsc(samples,sstep,proofnoise)

    optnoise = PowerLawNoise(1.0,256.0,1.8e-37,2.0,interplen=1)
    samples = 2**20 # 2**20 takes 13s on a 1.25GHz G4; 2*22 used for plot
    sstep = 0.1
    patches = 128
    noises_opt = getobsc(samples,sstep,optnoise)

    lasernoise = PowerLawNoise(1.0,256.0,1.1e-26,0.0,interplen=1)
    samples = 2**20 # 2**20 takes 13s on a 1.25GHz G4; 2*22 used for plot
    sstep = 0.1
    patches = 64
    noises_ls = getobsc(samples,sstep,lasernoise)

    #InterpolateNoise(noises_ls,len(noises_ls),sstep,Orbit.t[0],1
    #InterpolateNoise(1.0,256.0,1.1e-26,0.0,1)
    #NoisyLISA(lisa,Orbit.Dt,optnoise)
    #noise=TDInoise(lisa,noises_pr,noises_opt,noises_ls)

    #Simple Binary
    f=0
    phi0=0
    i=0
    A=0
    beta=0
    labda=0
    psi=0

    wave=SimpleBinary(f,phi0,i,A,beta,labda,psi)

    #TDI=TDIsignal(lisa,wave)

    stime=8
    TDI_noise = TDInoise(lisa_orb,
                           stime, 2.5e-48, # proof-mass noise parameters
                           stime, 1.8e-37, # optical-path noise parameters
                           stime, 1.1e-26) # laser frequency noise parameters
    TDI = TDInoise(lisa_orb,
                           stime, 0, # proof-mass noise parameters
                           stime, 0, # optical-path noise parameters
                           stime, 0) # laser frequency noise parameters
    print("Obtained TDI")


    def TDI_obs(lisa,TDI,stime,name=False):
        if name==False:
            name='NEWFILE'
        samples = 2**12 / stime #2**25 / stime
    #samples = 2**18 / stime

        patches = 25
        [noiseXm, noiseX1,noiseX2,noiseX3] = np.transpose(getobs(samples,stime,[TDI.Xm,TDI.X1,TDI.X2,TDI.X3]))

        myspecX=spect(noiseXm, stime,patches)
        myspecX1=spect(noiseX1, stime,patches)
        myspecX2=spect(noiseX2, stime,patches)
        myspecX3=spect(noiseX3, stime,patches)

        writearray('data/'+name+'.txt', myspecX[1:])

        return myspecX
    q1=TDI_obs(lisa_orb,TDI,stime,name='TDI_nonoise')
    q2=TDI_obs(lisa_orb,TDI_noise,stime,name='TDI_noise')

    os.system("gnuplot> load 'plotfile.plt'")
    #freqs = myspecX[1:,0]

    plt.figure()
    plt.plot(q1[:,0],q1[:,1],label="No Noise")
    plt.plot(q2[:,0],q2[:,1],label="Noise")
    plt.show()
