from synthlisa import *
import numpy as np
import matplotlib.pyplot as plt
import os
from fractions import Fraction
import math
import datetime
from class_orbit import orbit

home='/home/ester/Dropbox/Master_Project/Synthetic_LISA/'
filename='NGO_1M_10deg_synthlisa.txt'
directory_imp='LISAorbits/NGO_1M_10deg/'
num_back=0
Orbit=orbit(home=home,filename=filename,directory_imp=directory_imp,num_back=4)

lisa=Orbit.lisa_obj
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
NoisyLISA(lisa,Orbit.Dt,optnoise)
#noise=TDInoise(lisa,noises_pr,noises_opt,noises_ls)

