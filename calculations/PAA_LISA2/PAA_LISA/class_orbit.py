##!/usr/bin/env python

from synthlisa import *
import numpy as np
import matplotlib.pyplot as plt
import os
from fractions import Fraction
import math
import datetime
from scipy.interpolate import interp1d
year2sec=32536000
day2sec = year2sec/365.25

class orbit():
    def __init__(self,home=os.getcwd(),**kwargs):
        self.home=home
        self.num_back=kwargs.pop('num_back',0)
        self.filename=kwargs.pop('filename','None')
        self.directory_imp=kwargs.pop('directory_imp',False)
        if self.directory_imp != False:
            self.directory_imp=home+self.directory_imp
        directory_plot=kwargs.pop('directory_plot','/home/ester/git/synthlisa/figures')
        plot_on=kwargs.pop('plot_on',False)
        self.scale=kwargs.pop('scale',1)
        self.read_max=kwargs.pop('read_max','all')
        self.timeunit = kwargs.pop('timeunit','days')
        
        if self.filename=='None':
            print('Please select filename')
        else:
            self.import_file(read_max=self.read_max)

        self.calculations()
        if plot_on==True:
            self.plot_func(directory_plot)

    
    def ang_2vectors(self,v1,v2):
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

    def import_file(self,read_max='all'):
        directory=self.directory_imp
        file_orb=self.filename
        num_back=self.num_back
        
        
        par=['t','p1x','p1y','p1z','p2x','p2y','p2z','p3x','p3y','p3z']
        p=[[],[],[]]
        t=[]
        #direc=self.fold_dir()
        direc=self.directory_imp
        if directory==False:
            file=open(file_orb,'r')
        else:
            file=open(direc+file_orb,'r')
        line_num=1
        scale=self.scale
        line_count=0
        for line in file.readlines():
            if read_max!='all':
                if line_count==read_max:
                    break
            if line[0]!= '#':
                a=line.split(' ')
                cleanup=False
                while cleanup==False:
                    try:
                        a.remove('')
                    except ValueError:
                        cleanup=True

                b=[]
                read_check=True
                try:
                    for j in a:
                        b.append(np.float64(float(j)))
                    a=b
                    #scale=1000*2.5004844051425046)       
                    p[0].append(np.array([a[1]*scale,a[2]*scale,a[3]*scale]))
                    p[1].append(np.array([a[4]*scale,a[5]*scale,a[6]*scale]))
                    p[2].append(np.array([a[7]*scale,a[8]*scale,a[9]*scale]))
                    t.append(a[0]) # ... [s]
                except ValueError:
                    read_check=False
                    pass
                if read_check==True:
                    line_count=line_count+1
        p=np.array([p[0],p[1],p[2]],np.float64)
        if self.timeunit == 'days':
            self.t=(np.array(t) - np.array([t[0]]*len(t)))*day2sec #... in sec, fist point at t=0
        elif self.timeunit == 'years':
            self.t=(np.array(t) - np.array([t[0]]*len(t)))*years2sec #... in sec, fist point at t=0
        else: #Already in seconds
            self.t=np.array(t) - np.array([t[0]]*len(t)) #... in sec, fist point at t=0
        #print('t:')
        #print(self.t)
        #print('')
        Dt=self.t[1]-self.t[0] # Assuming Dt s constant
        self.lisa_obj=SampledLISA(p[0],p[1],p[2],Dt,self.t[0],2)
        self.p=p
        self.Dt=Dt
        self.pos=[self.t,p]
        self.par=par
        lisa_obj=self.lisa_obj
        self.linecount=line_count

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
        def COM(m1=2,m2=2,m3=2): #... adjust to mass proofmass
            p=self.p
            m=[m1,m2,m3]
            COM_list=[]
            for j in range(0,len(p[0])):
                com_vec=[]
                for i in range(0,len(p)):
                    com_vec.append(m[i]*p[i][j,:])
                COM = sum(com_vec)/(sum(m))

                COM_list.append([COM[0],COM[1],COM[2]])

            self.COM=np.array(COM_list)
            #self.COM_func = interp1d(self.t,self.COM)

            r_list=[]
            r_list_func=[]
            for i in range(1,4):
                i=i-1
                r=[]
                for j in range(0,len(self.t)):
                    r.append(self.COM[j,:] - self.p[i][j,:])
                r_list.append(np.array(r))
                r=np.array(r)
                
                [x,y,z] =[interp1d(self.t,r[:,0]),interp1d(self.t,r[:,1]),interp1d(self.t,r[:,2])]
                func_calc = lambda t: np.array([x(t),y(t),z(t)])
                r_list_func.append(func_calc)
            
            self.r = r_list
            self.r_func = r_list_func


            
            return np.array(COM_list)


        def normal_pos():
            normal_vec=[]
            normal_vec_func=[]
            for i in range(1,4):
                ret=[]
                #l_inc = (i+2)%3 - 1
                #l_out = (i+1)%3 - 1
                i_left = ((i+1)%3) - 1
                i_right = ((i+2)%3) - 1
                i_self = i-1

                for j in range(0,len(self.t)):
                    L_left = self.p[i_left][j,:] - self.p[i_self][j,:]
                    L_right = self.p[i_right][j,:] - self.p[i_self][j,:]

                    #L_inc = -self.p[l_inc][j,:] #... Now both vectors ar pointed from de spacecraft away
                    #L_out = self.p[l_out][j,:]
                    
                    ret_calc=np.cross(L_left,L_right)/(np.linalg.norm(L_left)*np.linalg.norm(L_right))
                    ret.append(ret_calc/np.linalg.norm(ret_calc))

                ret=np.array(ret)
                normal_vec.append(ret)
                [x,y,z] = [interp1d(self.t,ret[:,0]),interp1d(self.t,ret[:,1]),interp1d(self.t,ret[:,2])]
                func_calc = lambda t: np.array([x(t),y(t),z(t)])
                normal_vec_func.append(func_calc)

            self.n=normal_vec
            self.n_func=normal_vec_func
            return ret

        def angle_pos():
            ang_list_in = []
            ang_list_out=[]
            ang_stat_l_list=[]
            ang_stat_r_list=[]
            v_l_list=[]
            v_r_list=[]
            n_new_list=[]
            v_l_func=[]
            v_r_func=[]
            n_new_func=[]
            for i in range(1,4):
                i_left = ((i+1)%3) - 1
                i_right = ((i+2)%3) - 1
                i_self = i-1
                ang_in=[]
                ang_out=[]
                ang_stat_l=[]
                ang_stat_r=[]
                v_l_vec=[]
                v_r_vec=[]
                n_new_vec=[]
                for j in range(0,len(self.t)):
                    p_left=self.p[i_left][j,:]
                    p_right=self.p[i_right][j,:]
                    p_self=self.p[i_self][j,:]

                    v_l = p_left - p_self
                    v_r = p_right - p_self

                    ang_in.append(np.linalg.norm(np.cross(v_l,v_r))/(np.linalg.norm(v_l)*np.linalg.norm(v_r)))
                    ang_out.append(0)

                    ang_stat_l.append(self.ang_2vectors(v_l,self.r[i-1][j,:]))
                    ang_stat_r.append(self.ang_2vectors(v_r,self.r[i-1][j,:]))
                    v_l_vec.append(v_l)
                    v_r_vec.append(v_r)
                    n_calc=np.cross(v_l,v_r)
                    n_new_vec.append(n_calc/np.linalg.norm(n_calc))
                ang_list_in.append(np.array(ang_in))
                ang_list_out.append(np.array(ang_out))
                ang_stat_l_list.append(interp1d(self.t,np.array(ang_stat_l)))
                ang_stat_r_list.append(interp1d(self.t,np.array(ang_stat_r)))
                v_l_vec = np.array(v_l_vec)
                v_r_vec = np.array(v_r_vec)
                n_new_vec = np.array(n_new_vec)
                v_l_list.append(v_l_vec)
                v_r_list.append(v_r_vec)
                n_new_list.append(n_new_vec)

                [x,y,z] = [interp1d(self.t,v_l_vec[:,0]),interp1d(self.t,v_l_vec[:,1]),interp1d(self.t,v_l_vec[:,2])]
                func_calc = lambda t: np.array([x(t),y(t),z(t)])
                v_l_func.append(func_calc)
                [x,y,z] = [interp1d(self.t,v_r_vec[:,0]),interp1d(self.t,v_r_vec[:,1]),interp1d(self.t,v_r_vec[:,2])]
                func_calc = lambda t: np.array([x(t),y(t),z(t)])
                v_r_func.append(func_calc)
                [x,y,z] = [interp1d(self.t,n_new_vec[:,0]),interp1d(self.t,n_new_vec[:,1]),interp1d(self.t,n_new_vec[:,2])]
                func_calc = lambda t: np.array([x(t),y(t),z(t)])
                n_new_func.append(func_calc)



            self.ang_in = ang_list_in
            self.ang_out = ang_list_out
            self.ang_stat_l = ang_stat_l_list
            self.ang_stat_r = ang_stat_r_list
            self.v_l = v_l_list
            self.v_r = v_r_list
            self.n_new = n_new_list
            self.v_l_func = v_l_func
            self.v_r_func = v_r_func
            self.n_new_func = n_new_func
            return 0
   

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
        self.COM = COM()
        normal_pos()
        angle_pos()



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



