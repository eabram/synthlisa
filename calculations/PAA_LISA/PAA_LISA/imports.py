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
warnings.filterwarnings('ignore', 'The iteration is not making good progress')
#warnings.filterwarnings("error")
import scipy.optimize
from sympy import *
#from calc import *
from calc2 import *
import runfile
import plotfile2
import save_fig
import writefile
import functions as utils

year2sec=32536000
day2sec=year2sec/365.25
c=300000000

