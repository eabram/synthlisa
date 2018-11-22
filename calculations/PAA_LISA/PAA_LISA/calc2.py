from imports import *

year2sec=32536000
day2sec=year2sec/365.25
c=300000000

class PAA():
    def __init__(self,**kwargs):
        self.home = kwargs.pop('home',os.getcwd())
        self.filename = kwargs.pop('filename','')
        self.directory_imp = kwargs.pop('directory_imp','')
        self.read_max = kwargs.pop('read_max','all')
        self.num_back = kwargs.pop('num_back',0)
        self.plot_on = kwargs.pop('plot_on',True)
        self.scale = kwargs.pop('scale','Default')
        if self.scale=='Default':
            print('Getting scale by filename:')
            a = self.filename
            a1 = a.split('.')[0]
            a1 = a1.split('_')
            for k in range(0,len(a1)):
                if 'scale' == a1[k]:
                    self.scale = float(a1[k+1])
                    print(self.scale)
                    
        else:
            print(self.scale)
        print('')
        
        self.method = kwargs.pop('method','fsolve')        
        self.dir_savefig = kwargs.pop('dir_savefig',os.getcwd())
        self.dir_extr = kwargs.pop('dir_extr','')
        self.new_folder = kwargs.pop('new_folder',True)
        self.timeunit = kwargs.pop('timeunit','Default')
        if self.timeunit=='Default':
            print('Getting timestep by filename:')
            a = self.filename
            a1 = a.split('.')[0]
            a1 = a1.split('_')
            for k in range(0,len(a1)):
                if 'timestep' == a1[k]:
                    self.timeunit = a1[k+1]
                    print(self.timeunit)
            if self.timeunit!='days' and self.timeunit!='seconds':
                print('Could not obtain proper timestep')
        else:
            print(self.timeunit)
        print('')
            
        self.delay = kwargs.pop('delay',True)
        self.LISA = kwargs.pop('LISA',False)
        self.arm_influence = kwargs.pop('arm_influence',True)
        self.tstep = kwargs.pop('tstep',False)
        self.valorfunc = kwargs.pop('valorfunc','Value')
        self.calc_method = kwargs.pop('calc_method','Waluschka')
        print(self.calc_method)
        self.abb = kwargs.pop('abberation',False)
        import parameters
        parameters.do_obtain_var(parameters.do_para(),output=self)
    
    
    def PAA_func(self):
        import functions as utils
        print('')
        print('Importing Orbit')
        tic=time.clock()
        Orbit=orbit(home=self.home,filename=self.filename,directory_imp=self.directory_imp,num_back=self.num_back,scale=self.scale,read_max=self.read_max,plot_on=False,timeunit=self.timeunit)
        print(str(Orbit.linecount)+' datapoints')
        self.orbit = Orbit
        utils.LISA_obj(self,type_select='cache')
        print('Done in '+str(time.clock()-tic))

#        def velocity2(self,i,t,tdel=False,side='l'):
#            if tdel==False:
#                tdel = lambda time: 0
#
#            [i_self,i_left,i_right]=self.i_slr(i)
#            if side=='l':
#                i_next = i_left
#            elif side == 'r':
#                i_next=i_right
#            pos_next = func_pos(self,i_next,LISA=self.LISA)
#            pos_self = func_pos(self,i_self,LISA=self.LISA)
#            
#            if self.LISA==False:
#                pos = np.array([0,0,0])
#            
#            else:
#                pos_next = pos_next(t-tdel(t))
#                if self.calc_method=='Abram':
#                    pos_self = pos_self(t)
#                    #pos = np.array(self.LISA.putp(i_self,t)) - np.array(self.LISA.putp(i_left,t-tdel(t)))
#                elif self.calc_method=='Waluschka':
#                    pos_self = pos_self(t-tdel(t))
#                    #pos = np.array(self.LISA.putp(i_self,t-tdel(t))) - np.array(self.LISA.putp(i_left,t-tdel(t)))
#                #if side=='r':
#                #    if self.calc_method=='Abram':
#                #        
#                #        #pos = np.array(self.LISA.putp(i_self,t)) - np.array(self.LISA.putp(i_right,t-tdel(t)))
#                #    elif self.calc_method=='Waluschka':
#                #        pos = np.array(self.LISA.putp(i_self,t-tdel(t))) - np.array(self.LISA.putp(i_right,t-tdel(t)))
#
#            return [pos_self,pos_next]
#
#        def velocity_calc2(self,i,tdel,side):
#            def differentiate(func,t):
#                h=1e-5
#                ret=[]
#                for j in range(0,3):
#                    ret.append((func(t+h)[j]-func(t)[j])/h)
#                
#                return np.array(ret)
#
#            func_self = lambda t: velocity2(self,i,t,tdel,side=side)[0]
#            func_next = lambda t: velocity2(self,i,t,tdel,side=side)[1]
#            #v_self = differentiate(func_self)
#            #v_next = differentiate(func_next)
#            func_integ = lambda time: func_self(time) - func_next(time)
#            #v_ret = lambda time: v_self(time) - v_next(time)
#            v_ret = lambda time: differentiate(func_integ,time)
#
#            #time = Symbol("t")
#            #f = func(time)
#            #deriv = Derivative(f,time)
#            #ret = lambda t: deriv.doit().subs({time:t})
#
#            return v_ret
#
#        def abberation(self,i,t,u_func,tdel,side):
#            #u_func = np.array(u_func(t))
#            v = velocity_calc(self,i,tdel,side)
#            u_calc = u_func(i,t)
#            v_calc = np.array(v(t))
#            print('u_calc',u_calc,np.linalg.norm(u_calc))
#            print('v_calc',v_calc,np.linalg.norm(v_calc))
#            print(c)
#            q1 = v_calc
#            q2 = c*(u_calc/np.linalg.norm(u_calc))
#            print('q1q2',q1,q2)
#            u_new = ((q1+q2)/np.linalg.norm(q1+q2))*np.linalg.norm(u_calc)
#            print('u_new',u_new,np.linalg.norm(u_new))
#            print('')
#            #u_mag = np.linalg.norm(u_calc)
#            #u_calc = v_calc+c*(u_calc/np.linalg.norm(u_calc))
#            #u_new = u_mag*(u_calc/np.linalg.norm(u_calc))
#            #u_func_mag = np.linalg.norm(u_func(t))
#            #print((u_func(t)/np.linalg.norm(u_func(t))*c))
#            #print(v(t))
#            #u_func_new = v(t) + c*(u_func(t)/u_func_mag)
#            #print(u_func_new)
#            #print('')
#            #print(np.linalg.norm(v(t)),u_func_mag,c)
#            #u_func_new = u_func_mag*(u_func_new/np.linalg.norm(u_func_new))
#            #print(u_func_new)
#            return u_new
#
#        def abberation_calc(self,u_func_l,u_func_r):
#            u_func_new_l=[]
#            u_func_new_r=[]
#            
#            for i in range(1,4):
#                #tdel_l = lambda time: self.L_rl_func_tot(i,time)
#                #u_func_l = lambda time: self.u_l_func_tot_old(i,time)
#                u_func_new_l.append(lambda time: abberation(self,i,time,u_func_l,self.L_rl_func_tot_old[i-1],'l'))
#                #tdel_r = lambda time: self.L_rr_func_tor(i,time)
#                #u_func_r = lambda time: self.u_r_func_tot_old(i,time)
#                u_func_new_r.append(lambda time: abberation(self,i,time,u_func_r,self.L_rr_func_tot_old[i-1],'r'))
#            
#            return u_func_new_l,u_func_new_r

        # Calculations
        LA=utils.la()
        v_l_func_tot=[]
        v_r_func_tot=[]
        u_l_func_tot=[]
        u_r_func_tot=[]
        L_sl_func_tot_old=[]
        L_sr_func_tot_old=[]
        L_rl_func_tot_old=[]
        L_rr_func_tot_old=[]
        v_l_stat_func_tot=[]
        v_r_stat_func_tot=[]
        
        for i in range(1,4):
            #t_calc_vec=[]
            #v_l_calc=[]
            #v_r_calc=[]
            #u_l_calc=[]
            #u_r_calc=[]
            [[v_l_func,v_r_func,u_l_func,u_r_func],[L_sl_func,L_sr_func,L_rl_func,L_rr_func]] = utils.send_func(self,i,calc_method = self.calc_method)

            v_l_func_tot.append(v_l_func)
            v_r_func_tot.append(v_r_func)
            u_l_func_tot.append(u_l_func)
            u_r_func_tot.append(u_r_func)
            L_sl_func_tot_old.append(L_sl_func)
            L_sr_func_tot_old.append(L_sr_func)
            L_rl_func_tot_old.append(L_rl_func)
            L_rr_func_tot_old.append(L_rr_func)
            [i_self,i_left,i_right] = utils.i_slr(i)
            
            #pos_left_func = func_pos(self,i_left)
            #pos_self_func = func_pos(self,i_self)
            #pos_right_func = func_pos(self,i_right)
            v_l_stat_func_tot.append(utils.get_armvec_func(self,i_self,'l'))
            v_r_stat_func_tot.append(utils.get_armvec_func(self,i_self,'r'))


            #pos_rel_l = lambda time: self.LISA(i_left,time) - self.LISA(i_self,time)
            #pos_rel_r = lambda time: self.LISA(i_right,time) - self.LISA(i_self,time)

            #v_l_stat_func_tot.append(lambda time: pos_left_func(time) - pos_self_func(time))
            #v_r_stat_func_tot.append(lambda time: pos_right_func(time) - pos_self_func(time))
            #self.n_func.append(lambda time: LA.unit(np.cross(self.v_l_stat_func_tot[-1](time),self.v_r_stat_func_tot[-1](time))))
            #self.r_func.append(lambda time: r_calc(self.v_l_stat_func_tot[-1](time),self.v_r_stat_func_tot[-1](time),i))

        

        self.v_l_func_tot = utils.func_over_sc(v_l_func_tot)
        self.v_l_func_list = v_l_func_tot
        self.v_r_func_tot = utils.func_over_sc(v_r_func_tot)
        self.v_r_func_list = v_r_func_tot
        self.u_l_func_tot_old = utils.func_over_sc(u_l_func_tot)
        self.u_l_func_tot_list = u_l_func_tot
        self.u_r_func_tot_old = utils.func_over_sc(u_r_func_tot)
        self.u_r_func_tot_list = u_r_func_tot
        self.L_sl_func_tot = utils.func_over_sc(L_sl_func_tot_old)
        self.L_sr_func_tot = utils.func_over_sc(L_sr_func_tot_old)
        self.L_rl_func_tot = utils.func_over_sc(L_rl_func_tot_old)
        self.L_rr_func_tot = utils.func_over_sc(L_rr_func_tot_old)
        #self.u_l_func_tot_new_new,self.u_r_func_tot_new_new = abberation_calc(self,self.u_l_func_tot_old,self.u_r_func_tot_old)
        #self.u_l_func_tot_new = func_over_sc(self.u_l_func_tot_new_new)
        #self.u_r_func_tot_new = func_over_sc(self.u_r_func_tot_new_new)

        self.v_l_stat_func_tot = utils.func_over_sc(v_l_stat_func_tot)
        self.v_l_stat_func_list = v_l_stat_func_tot
        self.v_r_stat_func_tot = utils.func_over_sc(v_r_stat_func_tot)
        self.v_r_stat_func_list = v_r_stat_func_tot
        self.n_func = lambda i,t: LA.unit(np.cross(self.v_l_stat_func_tot(i,t),self.v_r_stat_func_tot(i,t)))
        self.r_func = lambda i,t: utils.r_calc(self.v_l_stat_func_tot(i,t),self.v_r_stat_func_tot(i,t),i)
             

        self.u_l_func_tot_calc = self.u_l_func_tot_old
        self.u_r_func_tot_calc = self.u_r_func_tot_old

        #--- Obtaining Velocity
        utils.velocity_func(self,hstep=100)
        utils.velocity_abs(self,hstep=100)

        #--- Obtaining PAA --- 
        print('Abberation: '+str(self.abb))
        selections=['l_in','l_out','r_in','r_out']
        PAA_func_val={}
        PAA_func_val[selections[0]] = lambda i,t: utils.calc_PAA_lin(self,i,t)
        PAA_func_val[selections[1]] = lambda i,t: utils.calc_PAA_lout(self,i,t)
        PAA_func_val[selections[2]] = lambda i,t: utils.calc_PAA_rin(self,i,t)
        PAA_func_val[selections[3]] = lambda i,t: utils.calc_PAA_rout(self,i,t)

        self.PAA_func = PAA_func_val 
       
        self.ang_breathing_din = lambda i, time: LA.angle(self.v_l_func_tot(i,time),self.v_r_func_tot(i,time))
        self.ang_breathing_stat = lambda i, time: LA.angle(self.v_l_stat_func_tot(i,time),self.v_r_stat_func_tot(i,time))

        return self #...adjust


