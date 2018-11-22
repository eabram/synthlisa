from imports import *
import functions as utils

### Plotting
def do_plot(self,dir_extr,filename,new_folder,tstep,plot_on=True):
    LA=utils.la()    
    #PAA_ret = self.PAA_ret
    
    #t_plot = self.t_plot 
    
 
    
    def string_length(l,string):
        while len(string)<l:
            string = '0'+string

        return string

    def get_date():
        now = datetime.datetime.now() 
        date=string_length(2,str(now.year))+string_length(2,str(now.month))+string_length(2,str(now.day))
        date=date+'-'+dir_extr
        return date

    def obtain_values(self,i,t_vec,func,key=False):
        y=[]
        if i!=False:
            for t in t_vec:
                y.append(func(i,t))
        else:
            for t in t_vec:
                y.append(func[key](t))

        return np.array(y)

    date_str=get_date()
    filename_save=filename.split('.')[-2]
    filename_save=filename_save.split('/')[-1]
    self.filename_save = filename_save
    dir_extr_new=dir_extr
    if self.dir_savefig==False:
        self.dir_savefig=os.path.dirname(os.path.realpath(__file__))+'/figures/'+date_str+'/'
    else:
        self.dir_savefig=self.dir_savefig+'/figures/'+date_str+'/'

    if not os.path.exists(self.dir_savefig):
        os.makedirs(self.dir_savefig)
    elif new_folder==True:
        i=1
        while os.path.exists(self.dir_savefig[0:-1]+'_('+str(i)+')/'):
            i=i+1
        dir_extr_new = dir_extr+'_('+str(i)+')/'
        self.dir_savefig = self.dir_savefig[0:-1]+'_('+str(i)+')/'
        os.makedirs(self.dir_savefig)


    self.write_info={}

    figs=[]
    
    if plot_on==True:     
        if tstep!=False:
            offset = int(round((Orbit.t[1]-Orbit.t[0])/tstep))
        else:
            offset=3
        self.offset = offset
        self.t_plot = self.t_all[offset:-offset]
        t_plot = self.t_plot
        self.calc={}

        # PAA PLOTS
        #--- PAA ---
        print('Offset: '+str(offset))
        f,ax = plt.subplots(4,3,figsize=(40,15))
        f.suptitle('PAA')
        plt.subplots_adjust(hspace=0.6,wspace=0.2)
        titles=['In plane left','Out of plane left','In plane right','Out of plane right']
        select=['l_in','l_out','r_in','r_out']
        y_all = self.PAA_func
        self.calc['PAA']={}
        for i in range(0,len(ax)):
            self.calc['PAA'][select[i]]=[]
            for j in range(0,len(ax[0])):
                x = t_plot
                func = y_all[select[i]]
                y = obtain_values(self,j+1,x,func)*1000000.0
                x_adj=x/day2sec
                ax[i,j].plot(x_adj,y,label='PAA')
                ax[i,j].axhline(max(y),label='Max='+"{0:.2e}".format(max(y)),color='0',linestyle='--')
                ax[i,j].axhline(min(y),label='Min='+"{0:.2e}".format(min(y)),color='0',linestyle='--')
                ax[i,j].axhline(sum(y)/len(y),label='Mean='+"{0:.2e}".format(sum(y)/len(y)),color='0',linestyle='-')

                ax[i,j].set_title(titles[i]+' SC'+str(j+1))
                ax[i,j].set_xlabel('Time (days)')
                ax[i,j].set_ylabel('Angle (microrad)')
                ax[i,j].legend(loc='best')

                self.calc['PAA'][select[i]].append(y)
                self.write_info[f._suptitle.get_text().replace('-','')+'-'+ax[i,j].get_title()+'_SC'+str(i+1)] = np.array([x_adj,y/1000000.0])
        figs.append(f)

        #--- PAA derivative
        f,ax = plt.subplots(4,3,figsize=(40,15))
        f.suptitle('Derivative PAA (per time)')
        plt.subplots_adjust(hspace=0.6,wspace=0.2)
        self.calc['PAAdiff'] = {}
        for i in range(0,len(ax)):
            self.calc['PAAdiff'][select[i]]=[]
            PAA_ret_diff_vec=[]
            for j in range(0,len(ax[i])):
                x=t_plot
                func = y_all[select[i]]
                y=obtain_values(self,j+1,x,func)
                x_calc=x/day2sec
                dfdx = np.diff(y)/np.diff(x_calc)
                x_adj = x_calc[0:-1]
                dfdx = dfdx*1000000
                ax[i,j].plot(x_adj,dfdx,label='PAA')
                ax[i,j].axhline(max(dfdx),label='Max='+"{0:.2e}".format(max(dfdx)),color='0',linestyle='--')
                ax[i,j].axhline(min(dfdx),label='Min='+"{0:.2e}".format(min(dfdx)),color='0',linestyle='--')
                ax[i,j].axhline(sum(dfdx)/len(dfdx),label='Mean='+"{0:.2e}".format(sum(dfdx)/len(dfdx)),color='0',linestyle='--')

                ax[i,j].set_title(titles[i]+' SC'+str(j+1))
                ax[i,j].set_xlabel('Time (days)')
                ax[i,j].set_ylabel('Angle/time (microrad/day)')
                ax[i,j].legend(loc='best')
                self.calc['PAAdiff'][select[i]].append(dfdx) 
                self.write_info[f._suptitle.get_text().replace('-','')+'-'+ax[i,j].get_title()+'_SC'+str(i+1)] = np.array([x_adj,dfdx])      
        self.calc['xPAAdiff'] = x_adj
        figs.append(f)

        #--- PAA Overview ---
        f,ax = plt.subplots(4,2,figsize=(40,15))
        f.suptitle('PAA - Overview')
        plt.subplots_adjust(hspace=0.6,wspace=0.2)
        for i in range(0,len(ax)):
            for j in range(0,3):
                
                ax[i,0].plot(t_plot,self.calc['PAA'][select[i]][j],label='SC'+str(j+1))
                ax[i,0].set_xlabel('Time (days)')
                ax[i,0].set_ylabel('Angle (microrad)')
                ax[i,0].set_title(titles[i])

                ax[i,1].plot(self.calc['xPAAdiff'],self.calc['PAAdiff'][select[i]][j],label='SC'+str(j+1))
                ax[i,1].set_xlabel('Time (days)')
                ax[i,1].set_ylabel('Angle/time (microrad/day)')
                ax[i,1].set_title(titles[i])
            ax[i,0].legend(loc='best')
            ax[i,1].legend(loc='best')
        figs.append(f)

        
        # Geometry plots
        #--- Breathing angle
        ang_breathing_din = lambda i, time: LA.angle(self.v_l_func_tot(i,time),self.v_r_func_tot(i,time))
        ang_breathing_stat = lambda i, time: LA.angle(self.v_l_stat_func_tot(i,time),self.v_r_stat_func_tot(i,time))
        f,ax = plt.subplots(3,2,figsize=(15,15))
        f.suptitle('Breathing angles')
        plt.subplots_adjust(hspace=0.6)
        ax[0,0].axhline(math.radians(60),color='0',linestyle='--')
        ax[1,0].axhline(math.radians(60),color='0',linestyle='--')
        ax[2,0].axhline(math.radians(0),color='0',linestyle='--')

        for j in range(0,len(ax)):
            ax[j,0].legend(loc='best')
            ax[j,0].set_xlabel('Time (sec)')
            ax[j,0].set_ylabel('Angle (rad)')
            if j<2:
                ax[j,1].legend(loc='best')
                ax[j,1].set_xlabel('Time (sec)')
                ax[j,1].set_ylabel('Angle/time (rad/sec)')

        ax[0,0].set_title('Angle by positions')
        ax[1,0].set_title('Angle by outgoing beams')
        ax[2,0].set_title('Difference of the two above angles')
        ax[0,1].set_title('Angle/time by positions')
        ax[1,1].set_title('Angle/time by outgoing beams')
        ax[2,1].axis('off')

        w1 = []
        w2=[]
        w3=[]
        w4=[]
        w5=[]

        for i in range(0,3):
            func_stat=obtain_values(self,i+1,t_plot,ang_breathing_stat)
            func_din=obtain_values(self,i+1,t_plot,ang_breathing_din)
            x=t_plot
            y=func_stat
            ax[0,0].plot(x,y,label='SC'+str(i+1))
            w1.append(x)
            w1.append(y)

            dfdx = np.diff(y)/np.diff(x)
            x_diff = x[0:-1]
            ax[0,1].plot(x_diff,dfdx,label='SC'+str(i+1))
            w2.append(x_diff)
            w2.append(dfdx)

            y=func_din
            ax[1,0].plot(x,y,label='SC'+str(i+1))
            w3.append(x)
            w3.append(y)

            dfdx = np.diff(y)/np.diff(x)
            ax[1,1].plot(x_diff,dfdx,label='SC'+str(i+1))
            w4.append(x_diff)
            w4.append(dfdx)

            y=func_stat-func_din
            ax[2,0].plot(t_plot,y,label='SC'+str(i+1))
            w5.append(t_plot)
            w5.append(y)

        w=[w1,w2,w3,w4,w5]
        count=0
        for i in ax:
            for j in i:
                count = count +1
                if count == 6:
                    break
                self.write_info[f._suptitle.get_text().replace('-','')+'-'+j.get_title()] = w[count-1] 
        self.w = w1
        figs.append(f)

        #--- Armlenghts ---
        f,ax = plt.subplots(1,2,figsize=(15,10))
        plt.subplots_adjust(wspace=2)
        f.suptitle('Armlengths')
        plt.subplots_adjust(hspace=0.6,wspace=0.2)
        w1=[]
        w2=[]
        func_l = lambda i,time: np.linalg.norm(self.v_l_stat_func_tot(i,time))
        func_r = lambda i,time: np.linalg.norm(self.v_r_stat_func_tot(i,time))
            
        for i in range(0,3):
            y_l = obtain_values(self,i+1,t_plot,func_l)
            y_r = obtain_values(self,i+1,t_plot,func_r)
    

            ax[0].plot(t_plot/day2sec,y_l/1000.0,label='L'+str(i+1))
            w1.append(t_plot)
            w1.append(y_l)

            ax[1].plot(t_plot/day2sec,y_r/1000.0,label='L'+str(i+1))
            w2.append(t_plot)
            w2.append(y_r)
        ax[0].set_title('Left')
        ax[0].set_xlabel('Time (days)')
        ax[0].set_ylabel('Length (km)')
        ax[0].legend(loc='best')
        ax[1].set_title('Right')
        ax[1].set_xlabel('Time (days)')
        ax[1].set_ylabel('Length (km)')
        ax[1].legend(loc='best')
        w=[w1,w2]
        for k in range(0,len(w)):
            self.write_info[f._suptitle.get_text().replace('-','')+'-'+ax[k].get_title()] = np.array(w[k])
        figs.append(f)
        

        #--- Relative Velocity ---
        f,ax = plt.subplots(3,3,figsize=(15,15))
        plt.subplots_adjust(wspace=2)
        f.suptitle('Relative velocity')
        plt.subplots_adjust(hspace=0.6,wspace=0.2)
        for i in range(0,3):
            x = t_plot/day2sec
            [i_self,i_left,i_right] = utils.i_slr(i+1)
            keyl = str(i_left)+str(i_self)
            keyr = str(i_right)+str(i_right)

            y_in_l = obtain_values(self,i+1,t_plot,self.v_in_mag_l)
            y_in_r = obtain_values(self,i+1,t_plot,self.v_in_mag_r)
            y_out_l = obtain_values(self,i+1,t_plot,self.v_out_mag_l)
            y_out_r = obtain_values(self,i+1,t_plot,self.v_out_mag_r)
            y_arm_l = obtain_values(self,i+1,t_plot,self.v_arm_mag_l)
            y_arm_r = obtain_values(self,i+1,t_plot,self.v_arm_mag_r)
            
            ax[i,0].set_title('Angular component of inplane velocity relative to SC '+str(i+1))
            ax[i,0].plot(x,y_in_l,label='left')
            ax[i,0].plot(x,y_in_r,label='right')
            ax[i,1].set_title('Out of plane velocity relative to SC '+str(i+1))
            ax[i,1].plot(x,y_out_l,label='left')
            ax[i,1].plot(x,y_out_r,label='right')
            ax[i,2].set_title('Radial component of inplane velocity relative to SC '+str(i+1))
            ax[i,2].plot(x,y_arm_l,label='left')
            ax[i,2].plot(x,y_arm_r,label='right')

            ax[i,0].set_xlabel('Time (days)')
            ax[i,0].set_ylabel('Velocity (m/s)')
            ax[i,0].legend(loc='best')
            ax[i,1].set_xlabel('Time (days)')
            ax[i,1].set_ylabel('Velocity (m/s)')
            ax[i,1].legend(loc='best')
            ax[i,2].set_xlabel('Time (days)')
            ax[i,2].set_ylabel('Velocity (m/s)')
            ax[i,2].legend(loc='best')
        figs.append(f)

        #--- PAA Estimate ---
        show_PAA_estimate=False
        if show_PAA_estimate:
            f,ax = plt.subplots(4,3,figsize=(40,15))
            plt.subplots_adjust(wspace=2)
            f.suptitle('PAA estimate')
            plt.subplots_adjust(hspace=0.6,wspace=0.2)

            f_lin = lambda i,time: np.arcsin(self.v_in_mag_l(i,time)*(self.L_rl_func_tot(1,time)+self.L_sl_func_tot(i,time))/np.linalg.norm(self.v_l_stat_func_tot(i,time)))*1000000
            f_lout = lambda i,time: np.arcsin(self.v_out_mag_l(i,time)*(self.L_rl_func_tot(1,time)+self.L_sl_func_tot(i,time))/np.linalg.norm(self.v_l_stat_func_tot(i,time)))*-1000000
            f_rin = lambda i,time: np.arcsin(self.v_in_mag_r(i,time)*(self.L_rr_func_tot(1,time)+self.L_sr_func_tot(i,time))/np.linalg.norm(self.v_r_stat_func_tot(i,time)))*1000000
            f_rout = lambda i,time: np.arcsin(self.v_out_mag_r(i,time)*(self.L_rr_func_tot(1,time)+self.L_sr_func_tot(i,time))/np.linalg.norm(self.v_r_stat_func_tot(i,time)))*-1000000
            t_vec=self.t_all
            for j in range(0,len(ax[0])):
                y_lin = obtain_values(self,j+1,t_plot,f_lin)
                y_lout = obtain_values(self,j+1,t_plot,f_lout)
                y_rin = obtain_values(self,j+1,t_plot,f_rin)
                y_rout = obtain_values(self,j+1,t_plot,f_rout)
                y_all = [y_lin,y_lout,y_rin,y_rout]

                for k in range(0,len(y_all)):
                    y = y_all[k]
                    ax[k,j].plot(t_plot,y,label='SC'+str(j+1))
                    ax[k,j].set_title('In plane left')
                    #ax[k,j].axhline(max(y),label='Max='+"{0:.2e}".format(max(y)),color='0',linestyle='--')
                    #ax[k,j].axhline(min(y),label='Min='+"{0:.2e}".format(min(y)),color='0',linestyle='--')
                    #ax[k,j].axhline(sum(y)/len(y),label='Mean='+"{0:.2e}".format(sum(y)/len(y)),color='0',linestyle='-')
                
                for i in range(0,len(ax)):
                    ax[i,j].set_xlabel('Time (days)')
                    ax[i,j].set_ylabel('Angle (microrad)')
                    ax[i,j].legend(loc='best')
                
            figs.append(f)

   
        
    self.figs = figs
    

    return self

