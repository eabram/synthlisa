from imports import *


### Plotting
def do_plot(self,dir_extr,filename,new_folder,tstep):
    
    PAA_ret = self.PAA_ret
    t_plot = self.t_plot 
    
    
    
    def string_length(l,string):
        while len(string)<l:
            string = '0'+string

        return string

    def get_date():
        now = datetime.datetime.now() 
        date=string_length(2,str(now.year))+string_length(2,str(now.month))+string_length(2,str(now.day))
        date=date+'-'+dir_extr
        return date

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
    if tstep!=False:
        offset = int(round((Orbit.t[1]-Orbit.t[0])/tstep))
    else:
        offset=3
    print('Offset: '+str(offset))
    #PAA=[PAA_in_l,PAA_in_r,PAA_out_l,PAA_out_r]
    f,ax = plt.subplots(4,3,figsize=(40,15))
    f.suptitle('PAA')
    plt.subplots_adjust(hspace=0.6,wspace=0.2)
    titles=['In plane left','Out of plane left','In plane right','Out of plane right']
    for i in range(0,len(ax)):
        for j in range(0,len(PAA_ret[i])):
            minp=offset
            maxp=len(t_plot[j])-offset
            #print('maxp: ',maxp)
            #print('')
            x=t_plot[j][minp:maxp]
            y=PAA_ret[i][j][minp:maxp]*1000000
            x_adj=x/day2sec
            ax[i,j].plot(x_adj,y,label='S'+str(j+1))
            ax[i,j].axhline(max(y),label='Max='+"{0:.2e}".format(max(y)),color='0',linestyle='--')
            ax[i,j].axhline(min(y),label='Min='+"{0:.2e}".format(min(y)),color='0',linestyle='--')
            ax[i,j].axhline(sum(y)/len(y),label='Mean='+"{0:.2e}".format(sum(y)/len(y)),color='0',linestyle='-')

            ax[i,j].set_title(titles[i])
            ax[i,j].set_xlabel('Time (days)')
            ax[i,j].set_ylabel('Angle (microrad)')
            ax[i,j].legend(loc='best')

            self.write_info[f._suptitle.get_text().replace('-','')+'-'+ax[i,j].get_title()+'_SC'+str(i+1)] = np.array([x_adj,y/1000000.0])
    figs.append(f)

    PAA_ret_diff=[]
    t_plot_diff=[]
    f,ax = plt.subplots(4,3,figsize=(40,15))
    f.suptitle('Derivative PAA (per time)')
    plt.subplots_adjust(hspace=0.6,wspace=0.2)        
    for i in range(0,len(ax)):
        PAA_ret_diff_vec=[]
        for j in range(0,len(PAA_ret[i])):
            minp=offset
            maxp=len(t_plot[j]) -offset
            x=t_plot[j][minp:maxp]
            y=PAA_ret[i][j][minp:maxp]
            dfdx = np.diff(y)/np.diff(x)
            PAA_ret_diff_vec.append(np.array(dfdx))
            x_calc=x/day2sec
            x_adj = x_calc[0:-1]
            dfdx = dfdx*1000000
            if i==0:
                t_plot_diff.append(x_adj)
            ax[i,j].plot(x_adj,dfdx,label='S'+str(j+1))
            ax[i,j].axhline(max(dfdx),label='Max='+"{0:.2e}".format(max(dfdx)),color='0',linestyle='--')
            ax[i,j].axhline(min(dfdx),label='Min='+"{0:.2e}".format(min(dfdx)),color='0',linestyle='--')
            ax[i,j].axhline(sum(dfdx)/len(dfdx),label='Mean='+"{0:.2e}".format(sum(dfdx)/len(dfdx)),color='0',linestyle='--')

            ax[i,j].set_title(titles[i])
            ax[i,j].set_xlabel('Time (days)')
            ax[i,j].set_ylabel('Angle/time (microrad/sec)')
            ax[i,j].legend(loc='best')
            self.write_info[f._suptitle.get_text().replace('-','')+'-'+ax[i,j].get_title()+'_SC'+str(i+1)] = np.array([x_adj,dfdx])
        PAA_ret_diff.append(PAA_ret_diff_vec)
    figs.append(f)
    self.PAA_ret_diff = PAA_ret_diff

    f,ax = plt.subplots(4,2,figsize=(40,15))
    f.suptitle('PAA - Overview')
    plt.subplots_adjust(hspace=0.6,wspace=0.2)        
    for i in range(0,len(ax)):
        for j in range(0,len(t_plot)):
            minp=offset
            maxp=len(t_plot[j])-offset
            ax[i,0].plot(t_plot[j][minp:maxp],PAA_ret[i][j][minp:maxp]*1000000,label='SC'+str(j+1))
            ax[i,0].set_xlabel('Time (days)')
            ax[i,0].set_ylabel('Angle (microrad)')
            ax[i,0].set_title(titles[i])
            
            ax[i,1].plot(t_plot_diff[j],PAA_ret_diff[i][j]*1000000,label='SC'+str(j+1))
            ax[i,1].set_xlabel('Time (days)')
            ax[i,1].set_ylabel('Angle/time (microrad/sec)')
            ax[i,1].set_title(titles[i])
        ax[i,0].legend(loc='best')
        ax[i,1].legend(loc='best')
    figs.append(f)


    f,ax = plt.subplots(3,2,figsize=(15,15))
    f.suptitle('Breathing angles')
    plt.subplots_adjust(hspace=0.6)
    ax[0,0].plot(t_plot[0][minp:maxp],(sum(self.ang_sc)/len(self.ang_sc))[minp:maxp],label='Mean')
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

    for i in range(0,len(self.ang_sc)):
        x=t_plot[i][minp:maxp]
        y=self.ang_sc[i][minp:maxp]
        minp=offset
        maxp=len(t_plot[j])-offset
        ax[0,0].plot(x,y,label='SC'+str(i+1))
        w1.append(x)
        w1.append(y)

        dfdx = np.diff(y)/np.diff(x)
        x_diff = x[0:-1]
        ax[0,1].plot(x_diff,dfdx,label='SC'+str(i+1))
        w2.append(x_diff)
        w2.append(dfdx)

        y=self.ang_beam[i][minp:maxp]
        ax[1,0].plot(x,y,label='SC'+str(i+1))
        w3.append(x)
        w3.append(y)

        dfdx = np.diff(y)/np.diff(x)
        ax[1,1].plot(x_diff,dfdx,label='SC'+str(i+1))
        w4.append(x_diff)
        w4.append(dfdx)

        y=self.ang_sc[i][minp:maxp]-self.ang_beam[i][minp:maxp]
        ax[2,0].plot(t_plot[i][minp:maxp],self.ang_sc[i][minp:maxp]-self.ang_beam[i][minp:maxp],label='SC'+str(i+1))
        w5.append(t_plot[i][minp:maxp])
        w5.append(self.ang_sc[i][minp:maxp] - self.ang_beam[i][minp:maxp])

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

    f,ax = plt.subplots(1,2,figsize=(15,10))
    plt.subplots_adjust(wspace=2)
    f.suptitle('Angle between sending and receiving beam')
    plt.subplots_adjust(hspace=0.6,wspace=0.2)
    titles=['In plane left','Out of plane left','In plane right','Out of plane right']
    w1=[]
    w2=[]
    for i in range(0,len(t_plot)):
        ax[0].plot(t_plot[i][minp:maxp]/day2sec,self.ang_sr_l[i][minp:maxp],label='SC'+str(i+1))
        w1.append(t_plot[i][minp:maxp])
        w1.append(self.ang_sr_l[i][minp:maxp])
        ax[1].plot(t_plot[i][minp:maxp]/day2sec,self.ang_sr_r[i][minp:maxp],label='SC'+str(i+1))
        w2.append(t_plot[i][minp:maxp])
        w2.append(self.ang_sr_r[i][minp:maxp])

    ax[0].set_title('Left')
    ax[0].set_xlabel('Time (days)')
    ax[0].set_ylabel('Angle (rad)')
    ax[0].legend(loc='best')
    ax[1].set_title('Right')
    ax[1].set_xlabel('Time (days)')
    ax[1].set_ylabel('Angle (rad)')
    ax[1].legend(loc='best')  

    w=[w1,w2]
    for k in range(0,len(w)):
        self.write_info[f._suptitle.get_text().replace('-','')+'-'+ax[k].get_title()] = np.array(w[k])

    figs.append(f)

    f,ax = plt.subplots(3,1,figsize=(15,10))
    f.suptitle('Wobbling angle')
    plt.subplots_adjust(hspace=0.6,wspace=0.2)
    w1=[]
    w2=[]
    w3=[]
    for i in range(0,len(t_plot)):
        ax[0].plot(t_plot[i][minp:maxp]/day2sec,self.ang_wob[i][minp:maxp],label='SC'+str(i+1))
        w1.append(t_plot[i][minp:maxp])
        w1.append(self.ang_wob[i][minp:maxp])
        ax[1].plot(t_plot[i][minp:maxp]/day2sec,self.ang_wob_stat[i][minp:maxp],label='SC'+str(i+1))
        w2.append(t_plot[i][minp:maxp])
        w2.append(self.ang_wob_stat[i][minp:maxp])
        ax[2].plot(t_plot[i][minp:maxp]/day2sec,self.ang_wob_diff[i][minp:maxp],label='SC'+str(i+1))
        w3.append(t_plot[i][minp:maxp])
        w3.append(self.ang_wob_diff[i][minp:maxp])

    ax[0].set_title('Dynamic positions while light is traveling')
    ax[0].set_xlabel('Time (days)')
    ax[0].set_ylabel('Angle (rad)')
    ax[0].legend(loc='best')
    ax[1].set_title('Static positions while light is traveling')
    ax[1].set_xlabel('Time (days)')
    ax[1].set_ylabel('Angle (rad)')
    ax[1].legend(loc='best')
    ax[2].set_title('Difference of the two figures above')
    ax[2].set_xlabel('Time (days)')
    ax[2].set_ylabel('Angle (rad)')
    ax[2].legend(loc='best')
    w=[w1,w2,w3]
    for k in range(0,len(w)):
        self.write_info[f._suptitle.get_text().replace('-','')+'-'+ax[k].get_title()] = np.array(w[k])

    figs.append(f)

    f,ax = plt.subplots(1,2,figsize=(15,10))
    plt.subplots_adjust(wspace=2)
    f.suptitle('Armlengths')
    plt.subplots_adjust(hspace=0.6,wspace=0.2)
    w1=[]
    w2=[]
    for i in range(0,len(t_plot)):
        ax[0].plot(t_plot[i][minp:maxp]/day2sec,self.L_l[i][minp:maxp]/1000.0,label='L'+str(i+1))
        w1.append(t_plot[i][minp:maxp])
        w1.append(self.L_l[i][minp:maxp])

        ax[1].plot(t_plot[i][minp:maxp]/day2sec,self.L_r[i][minp:maxp]/1000.0,label='L'+str(i+1))
        w2.append(t_plot[i][minp:maxp])
        w2.append(self.L_r[i][minp:maxp])
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


    f,ax = plt.subplots(2,1,figsize=(15,10))
    plt.subplots_adjust(wspace=2)
    f.suptitle('Velocity')
    plt.subplots_adjust(hspace=0.6,wspace=0.2)
    w1=[]
    w2=[]
    for i in range(0,len(self.v_inplane)):
        x = np.array(self.t_calc)[i][0:-1]/day2sec
        ax[0].plot(x,self.v_inplane[i],label='SC'+str(i+1))
        w1.append(x*day2sec)
        w1.append(self.v_inplane[i])
        ax[1].plot(x,self.v_outplane[i],label='SC'+str(i+1))
        w2.append(x*day2sec)
        w2.append(self.v_outplane[i])

    ax[0].set_title('In plane')
    ax[0].set_xlabel('Time (days)')
    ax[0].set_ylabel('Velocity (m/s)')
    ax[0].legend(loc='best')
    ax[1].set_title('Out of plane')
    ax[1].set_xlabel('Time (days)')
    ax[1].set_ylabel('Velocity (m/s)')
    ax[1].legend(loc='best')
    w=[w1,w2]
    for k in range(0,len(w)):
        self.write_info[f._suptitle.get_text().replace('-','')+'-'+ax[k].get_title()] = np.array(w[k])
    figs.append(f)

    f,ax = plt.subplots(3,1,figsize=(15,10))
    plt.subplots_adjust(wspace=2)
    f.suptitle('Retarded time')
    plt.subplots_adjust(hspace=0.6,wspace=0.2)
    labels = ['send left','send right','receive left','received right']
    for i in range(0,len(self.dt_sol)):
        for j in range(0,len(self.dt_sol[i])):
            ax[i].plot(self.dt_sol[i][j],label = labels[j])
        ax[i].set_title('SC'+str(i+1))
        ax[i].set_xlabel('Timestep (AU)')
        ax[i].set_ylabel('Time (sec)')
        ax[i].legend(loc='best')
    figs.append(f)

    f,ax = plt.subplots(3,3,figsize=(15,15))
    plt.subplots_adjust(wspace=2)
    f.suptitle('Relative velocity')
    plt.subplots_adjust(hspace=0.6,wspace=0.2)
    #labels = ['send left','send right','receive left','received right']
    for i in range(0,len(self.v_rel_l_inplane)):
        x = np.array(self.t_calc[i][0:-1])/day2sec
        ax[i,0].set_title('Angular component of inplane velocity relative to SC '+str(i+1))
        ax[i,0].plot(x,self.v_rel_l_inplane[i],label='left')
        ax[i,0].plot(x,self.v_rel_r_inplane[i],label='right')
        ax[i,1].set_title('Out of plane velocity relative to SC '+str(i+1))
        ax[i,1].plot(x,self.v_rel_l_outplane[i],label='left')
        ax[i,1].plot(x,self.v_rel_r_outplane[i],label='right')
        ax[i,2].set_title('Radial component of inplane velocity relative to SC '+str(i+1))
        ax[i,2].plot(x,self.v_rel_l_alongarm[i],label='left')
        ax[i,2].plot(x,self.v_rel_r_alongarm[i],label='right')
        
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

    if self.delay==True:
        f,ax = plt.subplots(4,3,figsize=(40,15))
        plt.subplots_adjust(wspace=2)
        f.suptitle('PAA estimate')
        plt.subplots_adjust(hspace=0.6,wspace=0.2)

        for i in range(0,len(self.t_calc)):
            x = np.array(self.t_calc[i][0:-1])/day2sec
            y_l_in = (self.v_rel_l_inplane[i][0:len(x)]*self.dt_l_tot[i][0:len(x)])/self.L_l[i][0:len(x)]
            y_l_in = y_l_in*1000000
            y_r_in = (self.v_rel_r_inplane[i][0:len(x)]*self.dt_r_tot[i][0:len(x)])/self.L_r[i][0:len(x)]
            y_r_in = y_r_in*1000000
            
            y_l_out = (self.v_rel_l_outplane[i][0:len(x)]*self.dt_l_tot[i][0:len(x)])/self.L_l[i][0:len(x)]
            y_l_out = y_l_out*1000000
            y_r_out = (self.v_rel_r_outplane[i][0:len(x)]*self.dt_r_tot[i][0:len(x)])/self.L_r[i][0:len(x)]
            y_r_out = y_r_out*1000000

            ax[0,i].plot(x,y_l_in,label = 'SC'+str(i+1))
            ax[0,i].set_title('In plane left')
            ax[0,i].axhline(max(y_l_in),label='Max='+"{0:.2e}".format(max(y_l_in)),color='0',linestyle='--')
            ax[0,i].axhline(min(y_l_in),label='Min='+"{0:.2e}".format(min(y_l_in)),color='0',linestyle='--')
            ax[0,i].axhline(sum(y_l_in)/len(y_l_in),label='Mean='+"{0:.2e}".format(sum(y_l_in)/len(y_l_in)),color='0',linestyle='-')
            
            ax[1,i].plot(x,y_l_out,label = 'SC'+str(i+1))
            ax[1,i].set_title('Out of plane left')
            ax[1,i].axhline(max(y_l_out),label='Max='+"{0:.2e}".format(max(y_l_out)),color='0',linestyle='--')
            ax[1,i].axhline(min(y_l_out),label='Min='+"{0:.2e}".format(min(y_l_out)),color='0',linestyle='--')
            ax[1,i].axhline(sum(y_l_out)/len(y_l_out),label='Mean='+"{0:.2e}".format(sum(y_l_out)/len(y_l_out)),color='0',linestyle='-')
            
            ax[2,i].plot(x,y_r_in,label = 'SC'+str(i+1))
            ax[2,i].set_title('In plane right')
            ax[2,i].axhline(max(y_r_in),label='Max='+"{0:.2e}".format(max(y_r_in)),color='0',linestyle='--')
            ax[2,i].axhline(min(y_r_in),label='Min='+"{0:.2e}".format(min(y_r_in)),color='0',linestyle='--')
            ax[2,i].axhline(sum(y_r_in)/len(y_r_in),label='Mean='+"{0:.2e}".format(sum(y_r_in)/len(y_r_in)),color='0',linestyle='-')
            
            ax[3,i].plot(x,y_r_out,label = 'SC'+str(i+1))
            ax[3,i].set_title('Out of plane right')
            ax[3,i].axhline(max(y_r_out),label='Max='+"{0:.2e}".format(max(y_r_out)),color='0',linestyle='--')
            ax[3,i].axhline(min(y_r_out),label='Min='+"{0:.2e}".format(min(y_r_out)),color='0',linestyle='--')
            ax[3,i].axhline(sum(y_r_out)/len(y_r_out),label='Mean='+"{0:.2e}".format(sum(y_r_out)/len(y_r_out)),color='0',linestyle='-')
            
            for j in range(0,len(ax[:,i])):
                ax[j,i].set_xlabel('Time (days)')
                ax[j,i].set_ylabel('Angle (microrad)')
                ax[j,i].legend(loc='best')
        figs.append(f)
    self.figs = figs
    self.offset = offset

    return self

