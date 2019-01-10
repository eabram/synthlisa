from synthlisa import *
import numpy as np
import matplotlib.pyplot as plt
import os
from fractions import Fraction
import math
import datetime
from scipy.interpolate import interp1d
from scipy.interpolate import RegularGridInterpolator
from control import * # package for control theory calculations
#from control.matlab import * # package for control theory calculations

import warnings
warnings.filterwarnings('ignore', 'The iteration is not making good progress')
#warnings.filterwarnings("error")
import scipy.optimize
import sympy as sp
from calc import *
import runfile
from WFE import *



year2sec=32536000
day2sec=year2sec/365.25
c=300000000

