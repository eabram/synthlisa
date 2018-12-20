from imports import *
from functions import *
from parameters import *
import PAA_LISA
import NOISE_LISA

def ttl(wfe,tele_control='nocontrol',PAAM_control_method='SS',simple=True):
    
    wfe.tele_control = tele_control
    wfe.PAAM_control_method = PAAM_control_method
    wfe.tele_aim()
    wfe.PAAM_control()
    
    ttl = {}
    wfe.ttl_pointing_function(option='all')
    ttl['pointing_all'] = wfe.ttl_l,wfe.ttl_r
    wfe.ttl_pointing_function(option='tilt')
    ttl['pointing_tilt'] = wfe.ttl_l,wfe.ttl_r
    wfe.ttl_pointing_function(option='piston')
    ttl['pointing_piston'] = wfe.ttl_l,wfe.ttl_r

    ttl_PAAM_l = []
    ttl_PAAM_r = []
    ttl_PAAM = wfe.Ndata.PAAMnoise(wfe=wfe)
    for i in range(1,4):
        [i_self,i_left,i_right] = PAA_LISA.utils.i_slr(i)
        keyl = str(i_self)+str(i_left)
        keyr = str(i_self)+str(i_right)
        ttl_PAAM_l.append(ttl_PAAM[keyl])
        ttl_PAAM_r.append(ttl_PAAM[keyr])

    ttl['PAAM'] = PAA_LISA.utils.func_over_sc(ttl_PAAM_l),PAA_LISA.utils.func_over_sc(ttl_PAAM_r)

    wfe.ttl_val = ttl
    
    return ttl


