from imports import *
from functions import *
from parameters import *
import calc_values

import PAA_LISA
import NOISE_LISA

def do_plot(wfe,i,side,title='',folder='Figures/TTL/',dt=36000):
    import os
    title=' of '+title 
    ttl = wfe.ttl_val
    t_vec = wfe.Ndata.t_all

    t0 = t_vec[1]
    tend = t_vec[-2]
    (t0-tend)/dt

    t_plot = np.linspace(t0,tend,round((tend-t0)/dt,1)+1)
    
    if side=='l':
        side=0
    elif side =='r':
        side=1

    ttl_sample={}
    for t in t_plot:
        for k in ttl.keys():
            if k not in ttl_sample.keys():
                ttl_sample[k]=[]
            ttl_sample[k].append(ttl[k][side](i,t))
    
    for k in ttl_sample.keys():
        ttl_sample[k] = np.array(ttl_sample[k])

    f,ax = plt.subplots(len(ttl_sample.keys()),1,figsize=(20,10))
    plt.subplots_adjust(hspace=0.6,wspace=0.2)
    for k in range(0,len(ttl_sample.keys())):
        key = ttl_sample.keys()[k]
        y = np.array(ttl_sample[key])
        ax[k].plot(t_plot,y,label='SC'+str(i))
        ax[k].set_title(key)
        ax[k].set_xlabel('time (s)')
        ax[k].set_ylabel('Length (m)')

    f.suptitle('TTL'+title+', telescope control='+wfe.tele_control+' , PAAM control = '+wfe.PAAM_control_method)
    
    directory = os.getcwd()+'/'+folder
    if not os.path.exists(directory):
        os.makedirs(directory)

    fig_title = wfe.tele_control+'_'+wfe.PAAM_control_method
    f.savefig(directory+fig_title+'.png')

    return [[t_plot,ttl_sample],[f,ax]]




