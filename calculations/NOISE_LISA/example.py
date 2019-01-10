import PAA_LISA
import NOISE_LISA
import os

import matplotlib.pyplot as plt
import numpy as np
import random
import os
from fractions import Fraction
import math
import datetime
from scipy.interpolate import interp1d
from scipy.interpolate import RegularGridInterpolator
import warnings
warnings.filterwarnings('ignore', 'The iteration is not making good progress')
#warnings.filterwarnings("error")
import scipy.optimize

year2sec=32536000
day2sec=year2sec/365.25
c=300000000

input_param = {
        'calc_method': 'Waluschka',
        'plot_on':False, #If plots will be made       'dir_savefig': os.getcwd(), # The directory where the figures will be saved. If False, it will be in the current working directory
        'noise_check':False,
        'home':'/home/ester/git/synthlisa/', # Home directory
        'directory_imp': False,
        'num_back': 0,
        'dir_orbits': '/home/ester/git/synthlisa/orbits/', # Folder with orbit files
        'length_calc': 40, # Length of number of imported datapoints of orbit files. 'all' is also possible
        'dir_extr': 'zzzWaluschka_no_abberation', # This will be added to the folder name of the figures
        'timeunit':'Default', # The timeunit of the plots (['minutes'],['days']['years'])
        'LISA_opt':True, # If a LISA object from syntheticLISA will be used for further calculations (not sure if it works properly if this False)
        'arm_influence': True, # Set True to consider the travel time of the photons when calculating the nominal armlengths
        'tstep':False,
        'delay':True, #'Not ahead' or False
        'method':'fsolve', # Method used to solve the equation for the photon traveling time
        'valorfunc':'Function', #
        'select':'Hallion', # Select which orbit files will be imported ('all' is all)
        'test_calc':False,
        'abberation':False,
        'delay': True
        }

# Make PAA_LISA object
data_all = PAA_LISA.runfile.do_run(input_param) #Dictionary with all PAA_LISA objects

# Make NOISE_LISA object
Ndata = NOISE_LISA.Noise(data=data_all['1']) # Using data_all['1'] as PAA_LISA argument

# Make WFE object
wfe = NOISE_LISA.WFE(Ndata=Ndata)

# Obtain telescope and PAAM pointing functions (WFE.AIM() object)
wfe.get_pointing(PAAM_method='SS',tele_method='no control')

# Define aperture (resolution)
self.pupil(Nbins=10)
wfe.speed_on=0 #=0 consider transmitting and receiving aperture surface (slow),=1 only transmitting surface and receiving point,=2 both transmitting and receiving point

# Calculations
i=1 # SC number
t = 12321 # Time (s)
side='l'

zmn,thmn = wfe.zern_aim(i,t,side=side)[1:3] # Zernike polynomials
power = wfe.u_rz_aperture(wfe.zmn,wfe.thmn,i,t,side=side,mode='power')
u = wfe.u_rz_aperture(wfe.zmn,wfe.thmn,i,t,side=side,mode='u')
phase = wfe.u_rz_aperture(wfe.zmn,wfe.thmn,i,t,side=side,mode='phase')

# TTL (incomplete)
NOISE_LISA.calc_values.ttl(wfe,PAAM_control_method='full control')
NOISE_LISA.plot_func(wfe).plot_ttl(1,side)







