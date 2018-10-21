##!/usr/bin/env python

from synthlisa import *
import numpy as np
import matplotlib.pyplot as plt
import os
from fractions import Fraction
import math
import datetime
from class_orbit import orbit

class laser_aim():
    # Assumed gaussian beam
    def __init__(self,orbit):
        # Import orbit and create LISA object
        home='/home/ester/Dropbox/Master_Project/Synthetic_LISA/'
        filename='NGO_1M_10deg_synthlisa.txt'
        directory_imp='LISAorbits/NGO_1M_10deg/'
        num_back=0
        self.orbit=orbit(home=home,filename=filename,directory_imp=directory_imp,num_back=4)
        self.lisa=Orbit.lisa_obj

        # Parameters
        self.I0=2 # Laser power [W]
        self.w0=0.01 # Waist radius (w(0)) [m]
        self.A=0.01 # 
        self.labda= #laser frequency

        pos=self.lisa.pos
        L=self.lisa.L
        ang=self.lisa.ang
        t=self.lisa.t

        front()

    def angle_spacecraft(self,pos):
                

    def point_ahead(self):

        return

    def waist(self,z):
        zR=(np.pi*(self.w0**2))/self.labda
        return w0*((1+((z/zR)**2))**0.5)

    def Intensity(self,x,y,z):
        r=((x**2)+(y**2))**0.5
        I=self.I0*((self.w0/self.w(z))**2)*exp((-2*(r**2))/w(z)**2)

        return I

    def Enery(self,I,angle=90,unit_ang='deg'):
        # Has to be adjusted for proper mirror orientation (integrate over surface)
        if unit_ang=='deg':
            angle=np.dag2rad(angle)
        
        laser_energy = I*np.sin(angle)*self.A
        return laser_enery

    
        


