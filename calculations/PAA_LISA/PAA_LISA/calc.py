from imports import *

year2sec=32536000
day2sec=year2sec/365.25
c=300000000

class la():
    def norm(self,v):
        return np.linalg.norm(v)

    def unit(self,v):
        return v/self.norm(v)

    def angle(self,v1,v2,dot=False):
        norm_v1 = self.norm(v1)
        norm_v2 = self.norm(v2)
        if norm_v1!=0 and norm_v2!=0:
            if dot==False:
                sin = self.norm(np.cross(v1,v2)/(norm_v1*norm_v2))
                return np.arcsin(sin)
            elif dot == True:
                cos = np.dot(v1,v2)/(norm_v1*norm_v2)
                return np.sign(np.dot(v1,v2))*np.arccos(cos)
        else:
            #print('norm v1: '+str(norm_v1))
            #print('norm v2: '+str(norm_v2))

            return np.nan

    def inplane(self,v,n):
        inplane_calc = v - (np.dot(v,n)/(np.linalg.norm(n)**2))*n
        return inplane_calc

    def outplane(self,v,n):
        outplane_calc = (np.dot(v,n)/(np.linalg.norm(n)**2))*n
        return outplane_calc

    def ang_out(self,v,n):
        sign = np.sign(np.dot(self.outplane(v,n),n))
        return sign*self.angle(self.inplane(v,n),v)

    def ang_in(self,v,n,r):
        #ang_out_calc = self.ang_out(v,n)
        inplane_calc = self.inplane(v,n)
        ang_in_calc = self.angle(inplane_calc,r)

        return ang_in_calc

    def ang_in_dot(self,v,v_stat,n,r):
        inplane_calc = self.inplane(v,n)
        costheta_beam = np.dot(inplane_calc,r)/(self.norm(inplane_calc)*self.norm(r))
        inplane_calc = self.inplane(v_stat,n)
        costheta_sc = np.dot(inplane_calc,r)/(self.norm(inplane_calc)*self.norm(r))

        return np.abs(np.arccos(costheta_beam)) - np.abs(np.arccos(costheta_sc))


    def ang_in_direct(self,v,v_stat,n,r):
        inplane_calc = self.inplane(v,n)
        inplane_stat = self.inplane(v_stat,n)
        ang_out_calc = self.angle(inplane_calc,inplane_stat)
        #ang1 = self.angle(inplane_calc,r)
        #ang2 = self.angle(inplane_stat,r)
        #sign = 1#np.sign(ang1 - ang2)

        return ang_out_calc#ang1-ang2

    def print_component(self,v,v_in,v_out,v_arm):
        n = self.norm(v)
        n_in = self.norm(v_in)
        n_out = self.norm(v_out)
        n_arm = self.norm(v_arm)

        print(n_in/n)
        print((n_out**2+n_in**2+n_arm**2)/n**2)
        print('')

        return 0

    def ang_in_out(self,v1,v2,n,v_stat,r):
        n = self.unit(n)
        v1_out = (np.dot(v1,n)*n)/(self.norm(n)**2)
        v1_arm = (np.dot(v1,v_stat)*v_stat)/(self.norm(v_stat)**2)
        v1_in = v1 - v1_out - v1_arm

        v2_out = (np.dot(v2,n)*n)/(self.norm(n)**2)
        v2_arm = (np.dot(v2,v_stat)*v_stat)/(self.norm(v_stat)**2)
        v2_in = v2 - v2_out - v2_arm

        ang_out_1 = np.arcsin(self.norm(v1_out)/self.norm(v1))
        ang_out_1 = ang_out_1 * np.sign(np.dot(v1_out,n))
        ang_out_2 = np.arcsin(self.norm(v2_out)/self.norm(v2))
        ang_out_2 = ang_out_2 * np.sign(np.dot(v2_out,n))
        ang_out = ang_out_1 - ang_out_2
        
        q = np.cross(n,v_stat)
        q = self.unit(q) 
      
        v1_in_calc = v1 - v1_out
        v2_in_calc = v2 - v2_out
        ang_in_1 = np.arcsin(self.norm(np.cross(v1_in_calc,r))/(self.norm(v1_in_calc)*self.norm(r)))
        ang_in_2 = np.arcsin(self.norm(np.cross(v2_in_calc,r))/(self.norm(v2_in_calc)*self.norm(r)))
        

        ang_in = self.norm(np.cross(v1_in_calc,v2_in_calc))/(self.norm(v1_in_calc)*self.norm(v2_in_calc)) 
        #ang_in_1 = self.norm(np.cross(v1_in_calc,v_stat))/(self.norm(v1_in_calc)*self.norm(v_stat)) 
        #ang_in_2 = self.norm(np.cross(v2_in_calc,v_stat))/(self.norm(v2_in_calc)*self.norm(v_stat)) 
        
        #ang_in = abs(ang_in_1)+abs(ang_in_2)
       


        return [ang_in,ang_out]




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
        self.dir_savefig = kwargs.pop('dir_savefig',False)
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
        #self.PAA_func()
         

    def PAA_func(self): 
        print('')
        print('Importing Orbit')
        tic=time.clock()
        Orbit=orbit(home=self.home,filename=self.filename,directory_imp=self.directory_imp,num_back=self.num_back,scale=self.scale,read_max=self.read_max,plot_on=False,timeunit=self.timeunit)
        print(str(Orbit.linecount)+' datapoints')
        print('Done in '+str(time.clock()-tic))


        def i_slr(i):
            '''Obtaining the correct spacecraft numbers'''

            i_self = i
            i_left = (i+1)%3
            i_right = (i+2)%3

            i_ret = [i_self,i_left,i_right]
            for j in range(0,len(i_ret)):
                if i_ret[j]==0:
                    i_ret[j]=3

            return i_ret
        
        
        def solve_num(func,guess,method='fsolve'):
            '''Select method for numerically solving an equation'''
            if method == 'fsolve':
                guess = np.array(guess)
                ret = scipy.optimize.fsolve(func,guess)
            elif method == 'excitingmixing':
                ret = scipy.optimize.excitingmixing(func,guess)
            elif method == 'linearmixing':
                ret = scipy.optimize.linearmixing(func,guess)
            elif method == 'newton_krylov':
                guess = [guess]
                ret = scipy.optimize.newton_krylov(func,guess)
            elif method == 'anderson':
                ret = scipy.optimize.anderson(func,guess)
            elif method == 'broyden1':
                ret = scipy.optimize.broyden1(func,guess)

            return ret
            
            
            
        def nominal_arm(orbit,i,t):
            def func_arm(i):
                L_vec=[]
                t_vec=orbit.t
                for j in range(0,len(t_vec)):
                    L_vec.append(np.linalg.norm(orbit.L[i-1][j]))

                f=interp1d(t_vec,L_vec,bounds_error)


                return f

            f=func_arm(i)

            return f(t)

        # Obtaining LISA object (or False)
        lisa=Orbit.lisa_obj
        func_nominal_arm = lambda i,time: nominal_arm(Orbit,i,time)
        print("")
        print('Obtaining proper LISA object')
        tic=time.clock()
        lisa_orb=PyLISA(lisa,func_nominal_arm)
        lisa_cache=CacheLISA(lisa_orb) # Works with retard() and putp and putn
        
        if self.LISA==True:
            LISA=lisa_cache
        else:
            LISA=False
        print('Done in '+str(time.clock()-tic))
        t=np.array(Orbit.t)
        self.lisasampled = Orbit.lisa_obj


        def func_pos(orbit,i,LISA=False):
            '''Generate functions of the positions'''
            if LISA==False:
                t_inter=orbit.t
                [x,y,z]=[orbit.p[i-1][:,0],orbit.p[i-1][:,1],orbit.p[i-1][:,2]]
                fx=interp1d(t_inter,x,bounds_error=False)
                fy=interp1d(t_inter,y,bounds_error=False)
                fz=interp1d(t_inter,z,bounds_error=False)
                return lambda time: np.array([fx(time),fy(time),fz(time)])
            else:
                L = lambda time: np.array(LISA.putp(i,time))
                return L

        def func_arm_sec(orbit,i,side,LISA=False): # Returns arm is seconds
            '''Obtains functions for nominal armlengths in seconds'''
            [i_self,i_left,i_right] = i_slr(i)
            pos_self = func_pos(orbit,i_self,LISA=LISA)
            if side=='r':
                pos_right = func_pos(orbit,i_righti,LISA=LISA)
                arms = lambda time: np.linalg.norm(pos_right(time) - pos_self(time))/c
            elif side=='l':
                pos_left = func_pos(orbit,i_left,LISA=LISA)
                arms = lambda time: np.linalg.norm(pos_left(time) - pos_self(time))/c
            return arms


        self.delay_func=[] #...remove (and verwijzingen) 

        def solve_L_PAA(t,orbit,pos_self,pos_left,pos_right,LISA=False,select='sl',method='Value',calc_method='Waluschka'):
            if LISA!=False:
                tl_guess = LISA.armlength(1,0)/c
            else:
                tl_guess= np.linalg.norm(Orbit.p[0][0,:] - Orbit.p[1][0,:])/c
            
            zzz = 'not_solve'
        
            
            s6 = lambda dt: np.linalg.norm(pos_left(t+dt) - pos_self(t+dt)) - c*dt
            res = scipy.optimize.brentq(s6,0,tl_guess*4)       
            self.zzzz = res 
            
            
            if zzz!='solve':
                try:
                    if select=='sl' or select=='rl':
                        s1 = lambda x: pos_left(x)
                    elif select=='sr' or select=='rr':
                        s1 = lambda x: pos_right(x)

                    s2 = lambda x: pos_self(x)
                    x_0 = t
                    if select=='sl' or select=='sr':
                        if calc_method=='Abram':
                            s3 = lambda dt: s1(x_0+dt) - s2(x_0)
                        else:
                            s3 = lambda dt: s1(x_0+dt) - s2(x_0+dt)
                    elif select=='rl' or select=='rr':
                        if calc_method=='Abram':
                            s3 = lambda dt: s1(x_0-dt) - s2(x_0)
                        else:
                            s3 = lambda dt: s1(x_0-dt) - s2(x_0-dt)
                    s4 = lambda dt: np.linalg.norm(s3(dt))
                    s5 = lambda dt: s4(dt) - c*dt
                    s6 = lambda dt: s5(dt)/c

                    res = scipy.optimize.brentq(s6,0,tl_guess*4)       
                    
                except ValueError:
                    res = np.nan
                    print('ValueError')
                    pass
                except TypeError:
                    print(type(x_0))
                    res = np.nan
                    pass

                return res

        def L_PAA(orbit,pos_self,pos_left,pos_right,LISA=False,method='Value',calc_method='Walushka'):
            '''Obtain time of flight of beam between spacecrafts'''
            t_inter = orbit.t
            
            tl_sol=[]
            tr_sol=[]
            trl_sol=[]
            trr_sol=[]
            t_send_l_sol=[]
            t_send_r_sol=[]
            t_rec_l_sol=[]
            t_rec_r_sol=[]
            
            catch_lim=np.array([0,0,0])

            t_min=t_inter[0]
            t_max=t_inter[-1]
           
            selections = ['sl','sr','rl','rr']

            for t in t_inter:
                
                tl_sol.append(solve_L_PAA(t,orbit,pos_self,pos_left,pos_right,LISA=False,select = selections[0]))
                tr_sol.append(solve_L_PAA(t,orbit,pos_self,pos_left,pos_right,LISA=False,select = selections[1]))
                trl_sol.append(solve_L_PAA(t,orbit,pos_self,pos_left,pos_right,LISA=False,select = selections[2]))
                trr_sol.append(solve_L_PAA(t,orbit,pos_self,pos_left,pos_right,LISA=False,select = selections[3]))

            L_sl=interp1d(t_inter,tl_sol,bounds_error=False)
            L_sr=interp1d(t_inter,tr_sol,bounds_error=False)
            L_rl=interp1d(t_inter,trl_sol,bounds_error=False)
            L_rr=interp1d(t_inter,trr_sol,bounds_error=False)
            dt_l_tot = np.array(tl_sol)+np.array(trl_sol)
            dt_r_tot = np.array(tr_sol)+np.array(trr_sol)
            self.dt_sol_vec = [tl_sol,tr_sol,trl_sol,trr_sol] 
            
            L_sl_func =  lambda time: solve_L_PAA(time,orbit,pos_self,pos_left,pos_right,LISA=False,select=selections[0],calc_method=calc_method)
            L_sr_func =  lambda time: solve_L_PAA(time,orbit,pos_self,pos_left,pos_right,LISA=False,select=selections[1],calc_method=calc_method)
            L_rl_func =  lambda time: solve_L_PAA(time,orbit,pos_self,pos_left,pos_right,LISA=False,select=selections[2],calc_method=calc_method)
            L_rr_func =  lambda time: solve_L_PAA(time,orbit,pos_self,pos_left,pos_right,LISA=False,select=selections[3],calc_method=calc_method)
            
            self.zzz = L_sl_func
            if method=='Value':
                return [[L_sl, L_sr, L_rl, L_rr],[dt_l_tot,dt_r_tot]]
            if method=='Function':
                return [[L_sl_func, L_sr_func, L_rl_func, L_rr_func],[dt_l_tot,dt_r_tot]]
                
        def n_r_lisa(i,LISA,m=[2,2,2]):
            '''Obtaining normal, r and COM vectors'''
            [i_self,i_left,i_right] = i_slr(i)

            v_l = lambda time: np.array(LISA.putp(i_left,time)) - np.array(LISA.putp(i_self,time))
            v_r = lambda time: np.array(LISA.putp(i_right,time)) - np.array(LISA.putp(i_self,time))
            COM = lambda time: (m[i_left-1]*np.array(LISA.putp(i_left,time)) + m[i_right-1]*np.array(LISA.putp(i_right,time)) + m[i_self-1]*np.array(LISA.putp(i_self,time)))/sum(m)
            
            r = lambda time: COM(time) - np.array(LISA.putp(i_self,time))

            n = lambda time: np.cross(v_l(time),v_r(time))
            
            return [n,r]
        
        def r_calc(v_l,v_r,i,m=[2,2,2]):

            [i_self,i_left,i_right] = i_slr(i)
            r =  (v_l*m[i_left-1]+v_r*m[i_right-1])/sum(m)
            
            return r



        def send_func(orbit,i,delay=self.delay,LISA=False,arm_influence=True,method=self.valorfunc,calc_method='Waluschka'):
            [i_self,i_left,i_right] = i_slr(i)
            
            pos_left = func_pos(orbit,i_left,LISA=LISA)
            pos_self = func_pos(orbit,i_self,LISA=LISA)
            pos_right = func_pos(orbit,i_right,LISA=LISA)

            if delay==True:
                if (LISA==False) or (LISA!=False and arm_influence == False):
                    pos_left_bad = func_pos(orbit,i_left,LISA=False)
                    pos_self_bad = func_pos(orbit,i_self,LISA=False)
                    pos_right_bad = func_pos(orbit,i_right,LISA=False)

                    [[L_sl,L_sr,L_rl,L_rr],[dt_l_tot,dt_r_tot]] = L_PAA(orbit,pos_self_bad,pos_left_bad,pos_right_bad,LISA=LISA,method=method,calc_method=calc_method)
                else:
                    [[L_sl,L_sr,L_rl,L_rr],[dt_l_tot,dt_r_tot]] = L_PAA(orbit,pos_self,pos_left,pos_right,LISA=LISA,method=method,calc_method=calc_method)
            elif delay=='Not ahead':
                L_sl = lambda t: np.linalg.norm(pos_left(t) - pos_self(t))/c #func_arm_sec(orbit,i_self,'l',LISA=LISA)
                L_sr = lambda t: np.linalg.norm(pos_right(t) - pos_self(t))/c #func_arm_sec(orbit,i_self,'r',LISA=LISA)
                L_rl=L_sl
                L_rr=L_sr
                dt_l_tot = False
                dt_r_tot = False

            elif delay==False:
                L_sl = lambda t: 0
                L_sr = lambda t: 0
                L_rl=L_sl
                L_rr=L_sr
                dt_l_tot = False
                dt_r_tot = False

            
            if calc_method=='Abram':
                #Abram2018
                v_send_l = lambda t: pos_left(t+L_sl(t)) - pos_self(t)
                v_send_r = lambda t: pos_right(t+L_sr(t)) - pos_self(t)
                v_rec_l = lambda t: pos_self(t) - pos_left(t - L_rl(t))
                v_rec_r = lambda t: pos_self(t) - pos_right(t - L_rr(t))
                
            else:
                #Waluschka2003
                v_send_l = lambda t: pos_left(t+L_sl(t)) - pos_self(t+L_sl(t))
                v_send_r = lambda t: pos_right(t+L_sr(t)) - pos_self(t+L_sr(t))
                v_rec_l = lambda t: pos_self(t-L_rl(t)) - pos_left(t - L_rl(t))
                v_rec_r = lambda t: pos_self(t-L_rl(t)) - pos_right(t - L_rr(t))


            return [[v_send_l,v_send_r,v_rec_l,v_rec_r],[dt_l_tot,dt_r_tot],[L_sl,L_sr,L_rl,L_rr]]
       
       
 
        LA=la()
        t_calc=[]
        v_l=[]
        v_r=[]
        u_l=[]
        u_r=[]


        t_calc_new=[]
        if self.tstep!=False:
            for i in range(1,4):
                tmin=Orbit.t[0]
                tmax=Orbit.t[-1]

                t_calc_new.append(np.linspace(tmin,tmax,((tmax-tmin)/tstep)+1))

        self.dt_sol=[]
        self.dt_l_tot=[]
        self.dt_r_tot=[]

        self.v_l_func_tot=[]
        self.v_r_func_tot=[]
        self.u_l_func_tot=[]
        self.u_r_func_tot=[]
        self.L_sl_func_tot=[]
        self.L_sr_func_tot=[]
        self.L_rl_func_tot=[]
        self.L_rr_func_tot=[]
        for i in range(1,4):
            t_calc_vec=[]
            v_l_calc=[]
            v_r_calc=[]
            u_l_calc=[]
            u_r_calc=[]
            [[v_l_func,v_r_func,u_l_func,u_r_func],[dt_l_tot_vec,dt_r_tot_vec],[L_sl_func,L_sr_func,L_rl_func,L_rr_func]] = send_func(Orbit,i,delay=self.delay,LISA=LISA,arm_influence=self.arm_influence)
            self.dt_l_tot.append(dt_l_tot_vec)
            self.dt_r_tot.append(dt_r_tot_vec)

            self.v_l_func_tot.append(v_l_func)
            self.v_r_func_tot.append(v_r_func)
            self.u_l_func_tot.append(u_l_func)
            self.u_r_func_tot.append(u_r_func)
            self.L_sl_func_tot.append(L_sl_func)
            self.L_sr_func_tot.append(L_sr_func)
            self.L_rl_func_tot.append(L_rl_func)
            self.L_rr_func_tot.append(L_rr_func)

            if self.delay == True:
                self.dt_sol.append(self.dt_sol_vec)
             
            if self.tstep!=False:
                t_vec = t_calc_new[i-1]
            else:
                t_vec = Orbit.t
            for t in t_vec:
                calc_check=True
                try:
                    v_l_tmp = v_l_func(t)
                    v_r_tmp = v_r_func(t)
                    u_l_tmp = u_l_func(t)
                    u_r_tmp = u_r_func(t)
                except ValueError:
                    calc_check=False
                    print('Not valid timestamp')
                    pass

                if calc_check==True:
                    v_l_calc.append(v_l_tmp)
                    v_r_calc.append(v_r_tmp)
                    u_l_calc.append(u_l_tmp)
                    u_r_calc.append(u_r_tmp)
                    t_calc_vec.append(t)

            v_l.append(np.array(v_l_calc))
            v_r.append(np.array(v_r_calc))
            u_l.append(np.array(u_l_calc))
            u_r.append(np.array(u_r_calc))
            t_calc.append(t_calc_vec)
        
        self.t_calc = t_calc
        t_plot=[]
        PAA_l_in=[]
        PAA_l_out=[]
        PAA_r_in=[]
        PAA_r_out=[]
        ang_sc=[]
        ang_beam=[]
        L_l=[]
        L_r=[]
        diff_L_l=[]
        diff_L_r=[]
        ang_sr_l=[]
        ang_sr_r=[]
        ang_wob=[]
        ang_wob_stat=[]
        ang_wob_diff=[]
        ang_beam_in_l=[]
        ang_beam_out_l=[]
        ang_beam_in_r=[]
        ang_beam_out_r=[]
        pos=[]
        normal_vec=[]

        self.PAA_func={}
        self.PAA_func['l_in']=[]
        self.PAA_func['l_out']=[]
        self.PAA_func['r_in']=[]
        self.PAA_func['r_out']=[]

        for i in range(1,4):        
            t_plot_vec=[]
            PAA_l_in_vec=[]
            PAA_l_out_vec=[]
            PAA_r_in_vec=[]
            PAA_r_out_vec=[]
            ang_sc_vec=[]
            ang_beam_vec=[]
            L_l_vec=[]
            L_r_vec=[]
            diff_L_l_vec=[]
            diff_L_r_vec=[]
            ang_sr_l_vec=[]
            ang_sr_r_vec=[]
            ang_wob_vec=[]
            ang_wob_stat_vec=[]
            ang_wob_diff_vec=[]       
            ang_beam_in_l_vec=[]
            ang_beam_out_l_vec=[]
            ang_beam_in_r_vec=[]
            ang_beam_out_r_vec=[]
            pos_vec=[]
            normal_vec_vec=[]

            [i_self,i_left,i_right] = i_slr(i)
            
            pos_left_func = func_pos(Orbit,i_left,LISA=LISA)
            pos_self_func = func_pos(Orbit,i_self,LISA=LISA)
            pos_right_func = func_pos(Orbit,i_right,LISA=LISA)
            if LISA!=False:
                [n_func,r_func]=n_r_lisa(i,LISA)
             

            
            # Obtain functions:
            v_l_stat_func = lambda time: pos_left_func(time) - pos_self_func(time)
            v_r_stat_func = lambda time: pos_right_func(time) - pos_self_func(time)
            v_l_func = self.v_l_func_tot[i-1]
            v_r_func = self.v_r_func_tot[i-1]
            u_l_func = self.u_l_func_tot[i-1]
            u_r_func = self.u_r_func_tot[i-1]           

            n_func = lambda time: LA.unit(np.cross(v_l_stat_func(time),v_r_stat_func(time)))
            r_func_func = lambda time: r_calc(v_l_stat_func(time),v_r_stat_func(time),i)

            # PAA
            calc_ang_l_in_func=lambda time: LA.ang_in_out(v_l_func(time),-u_l_func(time),n_func(time),v_l_stat_func(time),r_func_func(time))[0]
            calc_ang_l_out_func=lambda time: LA.ang_in_out(v_l_func(time),-u_l_func(time),n_func(time),v_l_stat_func(time),r_func_func(time))[1]
            calc_ang_r_in_func=lambda time: LA.ang_in_out(v_r_func(time),-u_r_func(time),n_func(time),v_r_stat_func(time),r_func_func(time))[0]
            calc_ang_r_out_func=lambda time: LA.ang_in_out(v_r_func(time),-u_r_func(time),n_func(time),v_r_stat_func(time),r_func_func(time))[1]
            
            self.PAA_func['l_in'].append(calc_ang_l_in_func)
            self.PAA_func['l_out'].append(calc_ang_l_out_func)
            self.PAA_func['r_in'].append(calc_ang_r_in_func)
            self.PAA_func['r_out'].append(calc_ang_r_out_func)


            #------------------
            

            i=i-1
            for j in range(0,len(t_calc[i])):     
                check_good=True
                                
                try:
                    pos_left = pos_left_func(t_calc[i][j])
                    pos_self = pos_self_func(t_calc[i][j])
                    pos_right = pos_right_func(t_calc[i][j])
                    v_l_stat = pos_left - pos_self
                    v_r_stat = pos_right - pos_self
                    
                    pos_vec.append(pos_self)

                    if LISA==False:
                        #n = Orbit.n_func[i](t_calc[i][j])
                        n = LA.unit(np.cross(v_l_stat,v_r_stat))
                        #n = Orbit.n_new_func[i](t_calc[i][j])
                        r = r_calc(v_l_stat,v_r_stat,i+1) #(v_l_stat+v_r_stat)/2.0#Orbit.r_func[i](t_calc[i][j])
                    else:                  
                        n = LA.unit(n_func(t_calc[i][j]))
                        r = r_func(t_calc[i][j])
                


                    normal_vec_vec.append(n)
                    #v_l_stat = Orbit.v_l_func[i](t_calc[i][j]) #...IS NOT WORKING
                    #v_r_stat = Orbit.v_r_func[i](t_calc[i][j])
                    
                    
                    #n = np.cross(v_l_stat,v_r_stat)
                    #n = n/np.linalg.norm(n)
                except ValueError:
                    check_good=False
                    print('not a good value')
                    pass
                if check_good == True:
                    v_l_calc=v_l[i][j,:]
                    v_r_calc=v_r[i][j,:]
                    u_l_calc=u_l[i][j,:]
                    u_r_calc=u_r[i][j,:]

                    t_plot_vec.append(t_calc[i][j])                   

                    # PAA
                    [ang_in_v_l,ang_out_v_l]=LA.ang_in_out(v_l_calc,r,n,v_l_stat,r)
                    [ang_in_v_r,ang_out_v_r]=LA.ang_in_out(v_r_calc,r,n,v_r_stat,r)
                    [ang_in_u_l,ang_out_u_l]=LA.ang_in_out(-u_l_calc,r,n,v_l_stat,r)
                    [ang_in_u_r,ang_out_u_r]=LA.ang_in_out(-u_r_calc,r,n,v_l_stat,r)
                    #------------------
                    
                    ang_sc_vec.append(LA.angle(v_l_stat,v_r_stat))
                    ang_beam_vec.append(LA.angle(v_l_calc,v_r_calc))
                   
                    [calc_ang_l_in,calc_ang_l_out]=LA.ang_in_out(v_l_calc,-u_l_calc,n,v_l_stat,r)
                    [calc_ang_r_in,calc_ang_r_out]=LA.ang_in_out(v_r_calc,-u_r_calc,n,v_r_stat,r)
                   
                    ang_beam_in_l_vec.append(calc_ang_l_in)
                    ang_beam_out_l_vec.append(calc_ang_l_out)
                    ang_beam_in_r_vec.append(calc_ang_r_in)
                    ang_beam_out_r_vec.append(calc_ang_r_out)

                    L_l_vec.append(LA.norm(v_l_calc))
                    L_r_vec.append(LA.norm(v_r_calc))

                    diff_L_l_vec.append(LA.norm(v_l_calc) - LA.norm(v_l_stat))
                    diff_L_r_vec.append(LA.norm(v_r_calc) - LA.norm(v_r_stat))
                    ang_sr_l_vec.append(LA.angle(v_l_calc,-u_l_calc,dot=False))
                    ang_sr_r_vec.append(LA.angle(v_r_calc,-u_r_calc,dot=False))
                                    
                    ang_wob_vec.append(LA.angle((v_l_calc+v_r_calc)/2.0,r))
                    ang_wob_stat_vec.append(LA.angle((v_l_stat+v_r_stat)/2.0,r))
                    ang_wob_diff_vec.append(ang_wob_vec[-1] - ang_wob_stat_vec[-1])
                    
                    PAA_l_in_vec.append(ang_in_v_l)
                    PAA_l_out_vec.append(ang_out_v_l)
                    PAA_r_in_vec.append(ang_in_v_r)
                    PAA_r_out_vec.append(ang_out_v_r)


            PAA_l_in.append(np.array(PAA_l_in_vec))
            PAA_l_out.append(np.array(PAA_l_out_vec))
            PAA_r_in.append(np.array(PAA_r_in_vec))
            PAA_r_out.append(np.array(PAA_r_out_vec))
            ang_sc.append(np.array(ang_sc_vec))
            ang_beam.append(np.array(ang_beam_vec))
            L_l.append(np.array(L_l_vec))
            L_r.append(np.array(L_r_vec))
            diff_L_l.append(np.array(diff_L_l_vec))
            diff_L_r.append(np.array(diff_L_r_vec))
            ang_sr_l.append(np.array(ang_sr_l_vec))
            ang_sr_r.append(np.array(ang_sr_r_vec))
            ang_wob.append(np.array(ang_wob_vec))
            ang_wob_stat.append(np.array(ang_wob_stat_vec))
            ang_wob_diff.append(np.array(ang_wob_diff_vec))
            ang_beam_in_l.append(np.array(ang_beam_in_l_vec))
            ang_beam_out_l.append(np.array(ang_beam_out_l_vec))
            ang_beam_in_r.append(np.array(ang_beam_in_r_vec))
            ang_beam_out_r.append(np.array(ang_beam_out_r_vec))
            pos.append(np.array(pos_vec))
            normal_vec.append(np.array(normal_vec_vec))

            t_plot.append(np.array(t_plot_vec))

        PAA_beam_next_sc = [PAA_l_in,PAA_l_out,PAA_r_in,PAA_r_out]
        PAA_ret = [ang_beam_in_l,ang_beam_out_l,ang_beam_in_r,ang_beam_out_r]
        
        other_ret=[ang_sc,ang_beam,L_l,L_r,diff_L_l,diff_L_r,ang_sr_l,ang_sr_r,ang_wob,ang_wob_stat,ang_wob_diff,PAA_beam_next_sc,t_plot]
       
        self.PAA_ret = PAA_ret
        self.PAA_beam_next_sc = PAA_beam_next_sc
        self.ang_sc = ang_sc
        self.ang_beam = ang_beam
        self.L_l = L_l
        self.L_r = L_r
        self.diff_L_l = diff_L_l
        self.diff_L_r = diff_L_r
        self.ang_sr_l = ang_sr_l
        self.ang_sr_r = ang_sr_r
        self.ang_wob = ang_wob
        self.ang_wob_stat = ang_wob_stat
        self.ang_wob_diff = ang_wob_diff
        self.t_plot = t_plot
        self.Orbit = Orbit
        self.LISA = LISA
        self.pos = pos
        self.normal_vec = normal_vec
        self.ang_sc = ang_sc        
        def func_over_sc(func_tot):
            f = lambda i,t: func_tot[i-1](t)

            return f

        for keys in self.PAA_func.keys():
            self.PAA_func[keys] = func_over_sc(self.PAA_func[keys])

        # Calculating velocity
        v_inplane=[]
        v_outplane=[]
        v_rel_l = []
        v_rel_r=[]
        v_abs=[]
        for i in range(0,len(t_calc)):
            v_abs_vec = []
            x = t_calc[i]
            y1 = pos[i][:,0]
            y2 = pos[i][:,1]
            y3 = pos[i][:,2]

            df1dt = np.diff(y1)/np.diff(x)
            df2dt = np.diff(y2)/np.diff(x)
            df3dt = np.diff(y3)/np.diff(x)

            x_new = x[0:-1]
            v_inplane_vec=[]
            v_outplane_vec=[]
            v_rel_l_vec=[]
            v_rel_r_vec = []
            for j in range(0,len(x_new)):
                n = normal_vec[i][j]
                v = np.array([df1dt[j],df2dt[j],df3dt[j]])
                inplane_calc = LA.inplane(v,n)
                v_inplane_vec.append(np.linalg.norm(inplane_calc))
                outplane_calc = v - inplane_calc
                v_outplane_vec.append(np.linalg.norm(outplane_calc)*(np.sign(np.dot(outplane_calc,n))))
                v_abs_vec.append(v)
            v_inplane.append(v_inplane_vec)
            v_outplane.append(v_outplane_vec)
            v_abs.append(np.array(v_abs_vec))

        self.v_abs = v_abs
        self.v_inplane = v_inplane
        self.v_outplane = v_outplane

        v_rel_l_inplane = []
        v_rel_r_inplane = []
        v_rel_l_outplane = []
        v_rel_r_outplane = []
        v_rel_l_alongarm= []
        v_rel_r_alongarm = []

        for i in range(0,len(self.v_abs)):
            [i_self,i_left,i_right] = i_slr(i+1)
            i_self = i_self - 1
            i_left = i_left - 1
            i_right = i_right -1
            #print(i_self,i_left,i_right)
            v_rel_l_inplane_vec = []
            v_rel_r_inplane_vec = []
            v_rel_l_outplane_vec = []
            v_rel_r_outplane_vec = []
            v_rel_l_alongarm_vec = []
            v_rel_r_alongarm_vec = []


            for j in range(0,len(self.v_abs[i])):               
                v_rel_l_calc = self.v_abs[i_left][j] - self.v_abs[i_self][j]
                v_rel_r_calc = self.v_abs[i_right][j] - self.v_abs[i_self][j]
                n = normal_vec[i][j]
               
                r_l = LA.unit(v_l[i_self][j])
                r_r = LA.unit(v_r[i_self][j])

                r_in_l = np.cross(n,r_l)
                r_in_r = np.cross(n,r_r)

                outplane_l = (np.dot(v_rel_l_calc,n)*n) / (np.linalg.norm(n)**2)
                alongarm_l = (np.dot(v_rel_l_calc,r_l)*r_l) / (np.linalg.norm(r_l)**2)
                inplane_l = v_rel_l_calc - outplane_l - alongarm_l
                outplane_r = (np.dot(v_rel_r_calc,n)*n) / (np.linalg.norm(n)**2)
                alongarm_r = (np.dot(v_rel_r_calc,r_r)*r_r) / (np.linalg.norm(r_r)**2)
                inplane_r = v_rel_r_calc - outplane_r - alongarm_r
                
                v_rel_l_outplane_vec.append(np.linalg.norm(outplane_l)*np.sign(np.dot(v_rel_l_calc,n)))
                v_rel_r_outplane_vec.append(np.linalg.norm(outplane_r)*np.sign(np.dot(v_rel_r_calc,n)))
                v_rel_l_inplane_vec.append(np.linalg.norm(inplane_l)*np.sign(np.dot(inplane_l,r_in_l)))
                v_rel_r_inplane_vec.append(np.linalg.norm(inplane_r)*np.sign(np.dot(inplane_r,r_in_r)))
                v_rel_l_alongarm_vec.append(np.linalg.norm(alongarm_l)*np.sign(np.dot(alongarm_l,r_l)))
                v_rel_r_alongarm_vec.append(np.linalg.norm(alongarm_r)*np.sign(np.dot(alongarm_r,r_r)))
                
            v_rel_l_inplane.append(np.array(v_rel_l_inplane_vec))
            v_rel_r_inplane.append(np.array(v_rel_r_inplane_vec))
            v_rel_l_outplane.append(np.array(v_rel_l_outplane_vec))
            v_rel_r_outplane.append(np.array(v_rel_r_outplane_vec))
            v_rel_l_alongarm.append(np.array(v_rel_l_alongarm_vec))
            v_rel_r_alongarm.append(np.array(v_rel_r_alongarm_vec))
            
        self.v_rel_l_inplane = v_rel_l_inplane
        self.v_rel_r_inplane = v_rel_r_inplane
        self.v_rel_l_outplane = v_rel_l_outplane
        self.v_rel_r_outplane = v_rel_r_outplane
        self.v_rel_l_alongarm = v_rel_l_alongarm
        self.v_rel_r_alongarm = v_rel_r_alongarm
        self.beam = [v_l,v_r,u_l,u_r]
        self.lisa_orb = lisa_orb
        self.lisa_cache = lisa_cache

        return [self,[[lisa_cache,lisa_orb,Orbit],[v_l,v_r,u_l,u_r],PAA_ret,other_ret]]

