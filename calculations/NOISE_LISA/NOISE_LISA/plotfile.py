from imports import *
from functions import *
from parameters import *
import calc_values

import PAA_LISA
import NOISE_LISA

class plot_func():
    def __init__(self,wfe,**kwargs):
        self.dt = kwargs.pop('dt',3600)
        self.make_t_plot(wfe,dt=self.dt)
        self.ttl_sample_all={}
        self.save_fig()
        self.wfe = wfe

    def save_fig(self,directory=False,extra_folder=False):
        if directory==False:
            directory = 'Figures/TTL/'
        directory = os.getcwd()+'/'+directory
        if not os.path.exists(directory):
            os.makedirs(directory)

        self.directory = directory



    def make_t_plot(self,wfe=False,t0=False,tend=False,dt=3600):
        if wfe==False:
            wfe=self.wfe
        if t0==False:
            t0 = wfe.Ndata.t_all[1]
        if tend==False:
            tend = wfe.Ndata.t_all[-2]
        N = int(np.round((tend-t0)/dt))+1
        self.t_plot = np.linspace(t0,tend,N)
    
    def plot_ttl(self,i,side,wfe=False,title='',dt=False):
        if wfe==False:
            wfe = self.wfe
        title=' of '+title 
        ttl = wfe.ttl_val
        t_vec = wfe.Ndata.t_all

        if dt!=False:
            self.make_t_plot(dt=dt)
        t_plot = self.t_plot
            
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
        
        fig_title = wfe.tele_control+'_'+wfe.PAAM_control_method
        f.savefig(self.directory+fig_title+'.png')
        
        self.f_ttl = f
        self.ttl_sample = ttl_sample
        key_tele = wfe.tele_control
        key_PAAM = wfe.PAAM_control_method
        if key_tele not in self.ttl_sample_all.keys():
            self.ttl_sample_all[key_tele]={}
        self.ttl_sample_all[key_tele][key_PAAM] = ttl_sample

        return [[t_plot,ttl_sample],[f,ax]]

    def plot_ttl_overview(self,ttl_sample_all=False):
        t_plot = self.t_plot
        if ttl_sample_all==False:
            ttl_sample_all = self.ttl_sample_all

        f,ax = plt.subplots(4,1,figsize=(20,10))
        plt.subplots_adjust(hspace=0.6,wspace=0.2)
        ttl0 = ttl_sample_all['no control']['nc']
        for tele_key in ttl_sample_all.keys():
            for PAAM_key in ttl_sample_all[tele_key].keys():
                ttl = ttl_sample_all[tele_key][PAAM_key]
                
                for k in range(0,len(ttl.keys())):
                    key = ttl.keys()[k]
                    y = ttl[key]-ttl0[key]
                    ax[k].semilogy(t_plot,y,label='tele='+tele_key+', PAAM='+PAAM_key)
        
        for i in range(0,len(ax)):
            ax[i].legend(loc='best')
            ax[i].set_title(ttl.keys()[i])
            ax[i].set_xlabel('time (s)')
            ax[i].set_ylabel('TTL (m)')
        f.suptitle('TTL difference compared to no telescope and PAAM control')

        f.savefig('Figures/TTL/Overview.png')

    def plot_P(self,i,side='l'):
        
        y = []
        for t in self.t_plot:
            y.append(self.wfe.P_calc(i,t,side=side))
        plt.plot(self.t_plot,y)
        plt.title('Power for SC'+str(i)+', side='+side)
        plt.xlabel('Time (s)')
        plt.ylabel('Relative power')
        plt.show()

