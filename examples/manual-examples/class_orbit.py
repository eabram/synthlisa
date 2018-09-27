##!/usr/bin/env python

from synthlisa import *
import numpy as np
import matplotlib.pyplot as plt
import os
from fractions import Fraction
import math
import datetime


class orbit():
    def __init__(self,home=os.getcwd(),**kwargs):
        self.home=home
        self.num_back=kwargs.pop('num_back',0)
        self.filename=kwargs.pop('filename','None')
        self.directory_imp=home+kwargs.pop('directory_imp','orbit_files/')
        directory_plot=kwargs.pop('directory_plot','/home/ester/git/synthlisa/figures')
        if self.filename=='None':
            print('Please select filename')
        else:
            self.import_file()

        self.calculations()
        self.plot_func(directory_plot)

    def import_file(self):
        directory=self.directory_imp
        file_orb=self.filename
        num_back=self.num_back
        
        
        par=['t','p1x','p1y','p1z','p2x','p2y','p2z','p3x','p3y','p3z']
        p=[[],[],[]]
        t=[]
        #direc=self.fold_dir()
        direc=self.directory_imp
        file=open(direc+file_orb,'r')
        line_num=1
        for line in file.readlines():
            a=line.split(' ')
            b=[]
            for j in a:
                b.append(np.float64(float(j)))
            a=b
            p[0].append(np.array([a[1],a[2],a[3]]))
            p[1].append(np.array([a[4],a[5],a[6]]))
            p[2].append(np.array([a[7],a[8],a[9]]))
            t.append(a[0])
        p=np.array([p[0],p[1],p[2]])
        Dt=t[1]-t[0] # Assuming Dt s constant
        lisa_obj=SampledLISA(p[0],p[1],p[2],Dt,t[0],2)
        self.p=p
        self.t=t
        self.Dt=Dt
        self.pos=[t,p]
        self.par=par
        self.lisa_obj=lisa_obj

    def fold_dir(self):
        '''returns num_back folders before direct'''
        
        direct=self.directory_imp
        num_back=self.num_back
        a=direct.split('/')
        ret=''

        i=0
        while i<len(a)-num_back:
            ret=ret+a[i]+'/'
            i=i+1
        print(ret)
        return ret

    def calculations(self):
        [t,[p1,p2,p3]]=self.pos

        ang1=[]
        ang2=[]
        ang3=[]
        def angle_e(v1,v2):
            ang=numpy.degrees(numpy.arccos(np.dot(v1,v2)/(np.linalg.norm(v1) * np.linalg.norm(v2))))

            return ang
        s1=[]
        s2=[]
        s3=[]
        for i in range(0,len(t)):
            s1.append(p2[i,:]-p3[i,:])
            s2.append(p3[i,:]-p1[i,:])
            s3.append(p1[i,:]-p2[i,:])

            ang1.append(angle_e(s3[-1],-s2[-1]))
            ang2.append(angle_e(s1[-1],-s3[-1]))
            ang3.append(angle_e(s2[-1],-s1[-1]))

        self.ang=[np.array(ang1),np.array(ang2),np.array(ang3)]
        self.L=[np.array(s1),np.array(s2),np.array(s3)]
   
    def plot_func(self,directory):
        print("Figures will be saved in: "+directory)
        def save_fig(f,directory,title=False,ext='.png'):
            if not os.path.exists(directory):
               os.makedirs(directory)
            if title==False:
                title="New_figure_"+datetime.datetime.now().isoformat()

            f.savefig(directory+title+ext)

        par=self.par

  
        
        f1, axarr = plt.subplots(9,1)
        count=1
        for i in range(0,len(self.p[:,:,:])):
            for j in range(0,len(self.p[i,0,:])):
                axarr[count-1].plot(self.t,self.p[i,:,j])
                axarr[count-1].set_title(par[count])
                count=count+1

        save_fig(f1,directory,title='positions')

        f2, axarr = plt.subplots(len(self.ang)+1,1)
        for i in range(0,len(self.ang)):
            axarr[i].plot(self.t-self.t[0],self.ang[i])
            axarr[i].set_xlabel('Time [days]')
            axarr[i].set_ylabel('Angle [deg]')
            axarr[i].set_title('Angle at spacecraft '+str(i+1))
        axarr[3].plot(self.t-self.t[0],sum(self.ang))
        save_fig(f2,directory,title='angles')


        plt.show()



