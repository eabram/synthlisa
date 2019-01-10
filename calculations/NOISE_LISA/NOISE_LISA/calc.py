from imports import *
from functions import *
from parameters import *
import PAA_LISA
class Noise():
    def __init__(self,**kwargs):
        self.data = kwargs.pop('data',False)
        self.t_all = self.data.t_all 

        if self.data==False:
            print('Please import a PAA_LISA object')
            ready=False
        else:
            print('Obtaining noise')
            ready=True

        if ready==True:
            y_all={}
            z_all={}

            PSD_laser = lambda f: 400 #...source
            f0 = 1e-6
            f_max = 1e-3
            print('Obaining lasernoise')
            self.lasernoise(PSD_laser,f0=f0,f_max=f_max,N=4096)
            y_all['laser'] = self.y_lasernoise
            z_all['laser'] = self.z_lasernoise

            print('Obtaining shotnoise')
            self.shotnoise()
            y_all['shot'] = self.y_shotnoise
            z_all['shot'] = self.z_shotnoise

            self.noise_y = y_all
            self.noise_z = z_all
            self.keys = self.y_lasernoise.keys()

            self.PAAMnoise()
            y_all['PAAM'] = self.y_PAAM
            z_all['PAAM'] = self.z_PAAM
            
            self.tele_l = lambda i,t: self.tele_control(i,t,option='no control',side='l') #...adjust for telescope length
            self.tele_r = lambda i,t: self.tele_control(i,t,option='no control',side='r')
            self.tele_l_fc = lambda i,t: self.tele_control(i,t,option='full control',side='l')
            self.tele_r_fc = lambda i,t: self.tele_control(i,t,option='full control',side='r')
            




    def Noise(self,f0,f_max,N,psd,unit='freq'):
        '''Obtains IFFT of psd between f0 and f_max'''
        M = N/2.0 +1
        df = (f_max-f0)/M

        f_lin = np.linspace(f0,f_max,M)
        psd_lin = psd(f_lin)

        ASD=[]
        phi = []
        Z_neg = []
        Z_pos = []
        for f in f_lin:
            if unit =='phasepercycle':
                ASD.append((((psd(f)**2)*2*np.pi)*2)**0.5)
            else:
                ASD.append((psd(f)*2)**0.5)
            phi.append(random.random()*2*np.pi)
            Z_pos.append(ASD[-1]*np.exp(1j*phi[-1]))
            Z_neg.append(ASD[-1]*np.exp(1j*random.random()*2*np.pi))


        Z1 = [0]*(N/2)
        Z_tot=[0]
        for i in range(0,len(Z_pos)):
            Z_tot.append(Z_pos[i])
        for i in range(0,len(Z_neg)):
            Z_tot.append(Z_neg[i])

        IFFT = np.fft.ifft(Z_tot)
        Dt = 1/(2.0*f_max)
        t0 = 0
        t_max = (N-1)*Dt
        t_IFFT = np.linspace(t0,t_max,len(IFFT))

        return [t_IFFT,IFFT]

    def Noise_time(self,f0,f_max,N,psd,t_stop,unit='freq',t=False,ret='all'):
        '''Returns the noise in time domain untill t_stop (or further)'''
        if psd!=False:
            t_max = 0
            count=0
            while t_max< t_stop:
                #print(count)
                [t_IFFT,IFFT] = self.Noise(f0,f_max,N,psd,unit=unit)
                if count!=0:
                    t_tot = np.concatenate((t_tot,t_IFFT+t_tot[-1]))
                    noise_tot = np.concatenate((noise_tot,IFFT))
                elif count==0:
                    t_tot = t_IFFT
                    noise_tot = IFFT
                t_max = t_tot[-1]
                #print(t_max)
                count=count+1

            t_ret = []
            noise_ret = []
            for i in range(0,len(t_tot)):
                if t_tot[i]<=t_stop:
                    t_ret.append(t_tot[i])
                    noise_ret.append(noise_tot[i])
            t_ret = np.array(t_ret)
            noise_ret = np.array(noise_ret)

            func_noise = interpolate(t_ret,np.real(noise_ret))
        
        if psd==False:
            func_noise = lambda y:0
            ret='function'
            
        if ret=='all':
            return [t_ret,noise_ret],func_noise

        elif ret=='function':
            return func_noise

    def lasernoise(self,PSD,f0=1e-6,f_max=1e-3,N=4096):
        '''Obatains lasernoise''' #..adjust for phase locking
                    
        def y_laser(data,C,C_star,t):
            y={}
            for s in range(1,4):
                for r in range(1,4):
                    check=True
                    if (s - r)%3 ==1:
                        sign = 'pos'
                        t_del = data.L_rl_func_tot(r,t)
            
                    elif (s - r)%3 ==2:
                        sign = 'neg'
                        t_del = data.L_rr_func_tot(r,t)
                    
                    elif s==r:
                        check=False
                    if check==True:
                        if sign=='pos':
                            y[str(s)+str(r)] = C_star[str(s)](t-t_del) - C[str(r)](t)
                        elif sign=='neg':
                            y[str(s)+str(r)] = C_star[str(s)](t-t_del) - C[str(r)](t)
            return y

        def z_laser(data,C,C_star,t):
            z={}
            for s in range(1,4):
                for r in range(1,4):
                    check=True
                    
                    if (s - r)%3 ==1:
                        sign = 'pos'
            
                    elif (s - r)%3 ==2:
                        sign = 'neg'
                    
                    elif s==r:
                        check=False
                    if check==True:
                        if sign=='pos':
                            z[str(s)+str(r)] = C_star[str(r)](t) - C[str(r)](t)
                        elif sign=='neg':
                            z[str(s)+str(r)] = C[str(r)](t) - C_star[str(r)](t)
                            
            return z

        def yz_laser_calc(data,C,C_star,yorz):
            t_all = self.data.t_all
            y_ret={}
            y_ret_interp={}
            for t in t_all:
                if yorz=='y':
                    y = y_laser(data,C,C_star,t) 
                elif yorz=='z':
                    y = z_laser(data,C,C_star,t)

                for keys in y.keys():
                    if keys not in y_ret.keys():
                        y_ret[keys] = []
                        t_ret[keys]=[]
                    y_ret[keys].append(y[keys])
                    t_ret[keys].append(t)
            
            for keys in y_ret:
                y_ret_interp[keys] = interpolate(t_all,y_ret[keys])

            return y_ret_interp
        
        C_func = {}
        C_func_star = {}
        t_ret = {}
        for j in range(0,6):
            i=j%3
            t_stop = self.data.t_plot[-1]
            [t_ret_calc,noise_ret],func_noise = self.Noise_time(f0,f_max,N,PSD,t_stop)
            

            if j<3:
                C_func[str(i+1)] = interp1d(t_ret_calc,np.real(noise_ret/nu_0),bounds_error=False)
            else:
                C_func_star[str(i+1)] = interp1d(t_ret_calc,np.real(noise_ret/nu_0),bounds_error=False)
            t_ret[str(i+1)] = t_ret_calc

        self.C_func = C_func
        self.C_func_star = C_func_star

        self.y_lasernoise = yz_laser_calc(self.data,C_func,C_func_star,'y')
        self.z_lasernoise = yz_laser_calc(self.data,C_func,C_func_star,'z')
    
    def shotnoise(self):
        def shotnoise_calc(self,i,t,side='l',retsel=0):
            data = self.data
            LA = PAA_LISA.la()
            if side=='l':
                v = data.v_l_func_tot(i,t)
                v_in = data.v_l_in_func_tot(i,t)
                v_out = data.v_l_out_func_tot(i,t)
            if side=='r':
                v = data.v_r_func_tot(i,t)
                v_in = data.v_r_in_func_tot(i,t)
                v_out = data.v_r_out_func_tot(i,t)
            
            L = np.linalg.norm(v) #...adjust    
            r = data.r_func(i,t)
            ang_v_in = LA.angle(v_in,r) - np.radians(30) #...adjust for telescope control
            ang_v_out = np.arcsin(np.linalg.norm(v_out)/L)
            P_rec0 = ((np.pi**2)*(D**4)*(eta_opt*P_L))/(16*(labda**2)*(L**2))
            P_rec = P_rec0*np.cos(ang_v_in)*np.cos(ang_v_out)
            phi_sn_q_vec = ((h/(eta_pd*(P_rec/4.0)))**0.5) # Shot noise per quadrant (upper limit)
            
            ret = [P_rec0,P_rec,phi_sn_q_vec]

            return ret[retsel]

        self.P0_l = lambda i,time: shotnoise_calc(self,i,time,side='l',retsel=0)
        self.P_l = lambda i,time: shotnoise_calc(self,i,time,side='l',retsel=1)
        self.phi_sn_l = lambda i,time: shotnoise_calc(self,i,time,side='l',retsel=2)
        self.P0_r = lambda i,time: shotnoise_calc(self,i,time,side='r',retsel=0)
        self.P_r = lambda i,time: shotnoise_calc(self,i,time,side='r',retsel=1)
        self.phi_sn_r = lambda i,time: shotnoise_calc(self,i,time,side='r',retsel=2)
        

        def y_shot(self,s,r):
            y={}
            # ... Neclecting telescope control
            check=False
            if (s - r)%3 ==1:
                P = self.P0_l(r,0) # .. at t=0 neglectig time dependance
                check=True 
            elif (s - r)%3 ==2:
                P = self.P0_r(r,0) # ... at t=0 neglecting time dependance
                check=True

            if check==True:
                hcbar = 3.16152649e-26
                psd_y_shot = lambda f: ((hcbar)/(2*np.pi*labda*P))*1e-18 #...adjust
                    
                [t_ret,noise_ret],func_noise = self.Noise_time(1e-8,1e-5,256,psd_y_shot,self.data.t_all[-1])
                return func_noise

        def yz_shot_calc(self,yorz):
            t_all = self.t_all
            y_ret={}
            y_ret_interp = {}

            for s in range(1,4):
                for r in range(1,4):
                    if s!=r:
                        key = str(s)+str(r)
                        y_ret[key] = y_shot(self,s,r)
            
            if yorz=='y':
                return y_ret
            elif yorz=='z':
                return False
        
        self.y_shotnoise = yz_shot_calc(self,'y')
        self.z_shotnoise = yz_shot_calc(self,'z')

    def tele_control(self,i,t,option='full control',side='l',L_tele=2,jitter=False):
        LA = PAA_LISA.la()
        r = self.data.r_func(i,t)
        n = self.data.n_func(i,t)
        
        if side=='l':
            v = self.data.v_l_func_tot(i,t)
        if side=='r':
            v = self.data.v_r_func_tot(i,t)
        
        x = LA.unit(np.cross(n,r))
        sign = np.sign(np.dot(x,v))
        if option=='no control':
            r_unit = LA.unit(r)
            ret_x = np.tan(np.radians(30))*np.linalg.norm(r_unit)*x*sign
            ret = LA.unit(ret_x+r_unit)*L_tele
        
        elif option=='full control':
            v_n = np.dot(v,LA.unit(n))*LA.unit(n)
            v_tele = v - v_n
            ret = LA.unit(v_tele)*L_tele
        
        if jitter!=False:
            jitter_ang_in = jitter[i-1][-2](t)
            jitter_ang_out = jitter[i-1][-1](t)
            ret = LA.rotate(ret,n,jitter_ang_in) #...check. Check orientation (when time dependent)
            ret = LA.rotate(ret,x,jitter_ang_out)
        
        return ret 

    def PAA_control(self,wfe=False):
        #Remove, only using pointing in WFE class

        '''Returns the point ahead angle''' #..adjust with controller
        magnification = self.data.MAGNIFICATION
        t_vec = self.t_all
        alpha={}
        alpha_func={}
        if wfe==False:
            print('Not using WFE')
            for keys in self.keys:
                alpha_vec=[]
                sc_s = int(keys[0])
                sc_r = int(keys[1])
                keys_diff = sc_r - sc_s
                
                if keys_diff == -1 or keys_diff == 2:
                    alpha_func[keys] = lambda time: self.data.PAA_func['l_out'](int(keys[1]),time)*0.5
                elif keys_diff == 1 or keys_diff == -2:
                    alpha_func[keys] = lambda time: self.data.PAA_func['r_out'](int(keys[1]),time)*0.5
                else:
                    raise ValueError('Wrong key')

            self.alpha_func = alpha_func

            return alpha_func
        
        else:
            print('Using PAAM pointing from WFE class')
            for i in range(1,4):
                [i_self,i_left,i_right] = PAA_LISA.utils.i_slr(i)
                keyl = str(i_self)+str(i_left)
                keyr = str(i_self)+str(i_right)
                alpha_func[keyl] = lambda time: (wfe.aim.PAAM_l_ang(i_self,time)*0.5)/magnification
                alpha_func[keyr] = lambda time: (wfe.aim.PAAM_r_ang(i_self,time)*0.5)/magnification

            self.alpha_func = alpha_func

            return alpha_func
    
    def create_noise(self,mean,sigma,t0,tend,dt=100):
        t_vec = np.linspace(t0,tend,int((tend-t0)/dt))
        ret = []
        for t in t_vec:
            ret.append(np.random.normal(loc=mean,scale=sigma))
        ret = np.array(ret)

        return interpolate(t_vec,ret)

    def PAAMnoise(self,C_func='default',C_func_star='default',wfe=False):
        if wfe==False:
            alpha_func = self.PAA_control()
            if C_func=='default':
                C_func = self.C_func
                C_func_star = self.C_func_star
        else:
            alpha_func = self.PAA_control(wfe=wfe)
            if C_func=='default':
                C_func = wfe.Ndata.C_func
                C_func_star = wfe.Ndata.C_func_star
        t_vec = self.t_all
        sigma_Dx = (1000e-6)/3.0
        sigma_Dy = (50e-6)/3.0

        def OPD_PAAM(alpha_func,Dx_func,Dy_func,r,t,dt=0.001,scale_din=5e20):
            alpha = alpha_func(t)
            alpha_dot = (alpha_func(t)-alpha_func(t-dt))/dt

            Dx = Dx_func(t)*(1+scale_din*alpha_dot)
            Dy = Dy_func(t)*(1+scale_din*alpha_dot)

            OPD = (Dy-Dx*np.tan(alpha))*(np.sin(alpha)/np.sin(np.radians(135)))*(1-np.cos(np.radians(90)-2*alpha))
            return OPD
        
        y_PAAM={}
        
        ttl={}
        for keys in alpha_func.keys():
            Dx = self.create_noise(0,sigma_Dx,t_vec[0],t_vec[-1]) # ... creating discrete uncorrelated noise may cause signal in frequency range
            Dy = self.create_noise(0,sigma_Dy,t_vec[0],t_vec[-1]) # ... creating discrete uncorrelated noise may cause signal in frequency range
            
            sc_s = int(keys[0])
            sc_r = int(keys[1])
            diff_sc = sc_s - sc_r
            if diff_sc == 1 or diff_sc == -2:
                side = 'l'
                C = lambda time: C_func[str(sc_r)](time-self.data.L_sl_func_tot(sc_r,time))
            elif diff_sc == -1 or diff_sc == 2:
                side = 'r'
                C = lambda time: C_func_star[str(sc_r)](time-self.data.L_sr_func_tot(sc_r,time))
            OPD = lambda time: OPD_PAAM(alpha_func[keys],Dx,Dy,sc_r,time) 
            ttl[keys] = OPD
            OPD_phase = lambda time: (OPD(time)*2*np.pi*nu_0*(C(time))/c)
            
            y_PAAM[keys] = OPD_phase


#            for i in range(0,len(alpha_func[keys])):
#                Dx = np.random.normal(loc=0,scale=sigma_Dx)
#                Dy = np.random.normal(loc=0,scale=sigma_Dy)
#                if i==0:
#                    delta_alpha.append(0)
#                    delta_Dx.append(0)
#                else:
#                    delta_alpha.append(alpha[keys][i] - alpha[keys][i-1])
#                    delta_Dx.append(Dx[i]-Dx[i-1])
#
#                angular_jitter.append((2**0.5)*(Dx[i]*alpha[keys][i] +Dy[i])*delta_alpha[i])
#                long_jitter.append((2**0.5)*(1+0.5*(alpha[keys][i]**2))*delta_Dx[i])
#                rot_long_jitter.append((2**0.5)*delta_alpha[i])
#
#                
#                OPD = (Dy[i]-Dx[i]*np.tan(alpha[keys][i]))*(np.sin(alpha[keys][i])/np.sin(np.radians(135)+alpha[keys][i]))*(1-np.cos(np.radians(90) - 2*alpha[keys][i]))
#                OPD_phase.append((OPD*2*np.pi*nu_0*(C_func_calc(self.t_all[i]))/c))
#                
#                
#                #OPD_PAAM.append((Dy[i]-Dx[i]*np.tan(alpha[i]))*(np.sin(alpha[i])/np.sin(np.radians(135)+alpha[i]))*(1-np.cos(np.radians(90) - 2*alpha[i])))
#
                
#            y_PAAM[keys] = interpolate(self.t_all,np.array(OPD_phase))

        self.y_PAAM = y_PAAM
        self.z_PAAM = False

        return ttl

    def PAA_point_calc(self,i,t,side='l',noise=False):#..adjust with noise
        LA = PAA_LISA.la()
        data = self.data
        
        [i_self,i_left,i_right] = PAA_LISA.utils.i_slr(i)
        keyl = str(i_self)+str(i_left)
        keyr = str(i_self)+str(i_right)
        n = data.n_func(i_self,t)
        r = data.r_func(i_self,t)
        rotas = LA.unit(np.cross(r,n))
        if side=='l':
            beam_send = LA.rotate(data.v_l_func_tot(i_self,t),rotas,self.alpha_func[keyl](t))
        elif side=='r':
            beam_send = LA.rotate(data.v_r_func_tot(i_self,t),rotas,self.alpha_func[keyr](t))

        return beam_send

    def PAA_point(self):
        self.PAA_point_l = lambda i,t: self.PAA_point_calc(i,t,side='l')
        self.PAA_point_r = lambda i,t: self.PAA_point_calc(i,t,side='r')
        
        return 0









class TDI():
    def __init__(self,**kwargs):
        print('TDI')
        noise = kwargs.pop('noise',False) 
        data = noise.data
        t_vec = noise.t_all
        y_all = noise.noise_y
        z_all = noise.noise_z
        
        para = ['X','Y','Z','alpha','beta','gamma','X1','Y1','Z1','P','Q','R']
        tdi_all = {}
        tdi_all_func = {}
        tdi_all['all']={}
        tdi_all_func['all']={}

        count=0
        for keys in y_all.keys():
            count = count+1
            tdi_all[keys]={}
            tdi_all_func[keys]={}
            for i in range(0,len(para)):
                X_calc,func= self.tdi_val(data,t_vec,y_all,z_all,y_all_include=[keys],para=para[i])
                #func = lambda time: interpolate(t,X_calc)
                tdi_all[keys][para[i]],tdi_all_func[keys][para[i]] = X_calc,func

                if count==1:
                    X_calc,func= self.tdi_val(data,t_vec,y_all,z_all,y_all_include='all',para=para[i])
                    #X_func= lambda time: self.tdi_val(data,time,y_all,z_all,y_all_include='all',para=para[i],function=True)
                    
                    #func = interpolate(t,X_calc)
                    tdi_all['all'][para[i]],tdi_all_func['all'][para[i]] = X_calc,func

        self.tdi_all = tdi_all
        self.tdi_all_func = tdi_all_func
        #self.X_func = X_func


    def tdi_val(self,data,t_vec,y_all,z_all,y_all_include='all',side='l',generation=0,para='X',print_on=False,delay_on=True,function=False):
        
        if para=='X':
            shift = 0
            method='XYZ'
        elif para=='Y':
            shift = 1
            method='XYZ'
        elif para == 'Z':
            shift = 2
            method='XYZ'
        elif para=='alpha':
            shift = 0
            method='aby'
        elif para=='beta':
            shift = 1
            method='aby'
        elif para == 'gamma':
            shift = 2
            method='aby'
        elif para=='X1':
            shift = 0
            method='X1Y1Z1'
        elif para=='Y1':
            shift = 1
            method='X1Y1Z1'
        elif para == 'Z1':
            shift = 2
            method='X1Y1Z1'
        elif para=='P':
            shift = 0
            method='PQR'
        elif para=='Q':
            shift = 1
            method='PQR'
        elif para == 'R':
            shift = 2
            method='PQR'
        
        
        
        
        
        id={}
        for kk in range(-3,4):
            if kk!=0:
                sign = np.sign(kk)
                abs_kk = abs(kk) - 1
                
                id[str(kk)] = (((abs_kk+shift)%3)+1)*sign
        if print_on==True:
            print(id)
        
        def id_convert(s,r):
            return str(id[str(s)])+str(id[str(r)])
            
        X = []
        
        def obtain_tdi_aby(data,id,para,y,t):
            a=[]
            
            t_del = delay(data,[],t,para=para,delay_on=delay_on)
            a.append(y[id_convert(3,1)](t-t_del))
            
            t_del = delay(data,[],t,para=para,delay_on=delay_on)
            a.append(-y[id_convert(2,1)](t-t_del))
            
            t_del = delay(data,[id['2']],t,para=para,delay_on=delay_on)
            a.append(y[id_convert(2,3)](t-t_del))
            
            t_del = delay(data,[id['3']],t,para=para,delay_on=delay_on)
            a.append(-y[id_convert(3,2)](t-t_del))
            
            t_del = delay(data,[id['1'],id['2']],t,para=para,delay_on=delay_on)
            a.append(y[id_convert(1,2)](t-t_del))

            t_del = delay(data,[id['1'],id['3']],t,para=para,delay_on=delay_on)
            a.append(-y[id_convert(1,3)](t-t_del))

            a = np.array(a)
            
            return sum(np.real(a))
        
        def obtain_tdi_P(data,id,para,y,z,t):
            P = []
            
            
            t_del = delay(data,[id['-2']],t,para=para,delay_on=delay_on)
            P.append(y[id_convert(1,2)](t+t_del))
            
            t_del = delay(data,[id['3']],t,para=para,delay_on=delay_on)
            P.append(-y[id_convert(1,3)](t+t_del))
            
            t_del = delay(data,[id['-2']],t,para=para,delay_on=delay_on)
            P.append(-y[id_convert(3,2)](t+t_del))
            
            t_del = delay(data,[id['3']],t,para=para,delay_on=delay_on)
            P.append(y[id_convert(2,3)](t+t_del))
            
            t_del = delay(data,[id['1'],id['3']],t,para=para,delay_on=delay_on)
            P.append(y[id_convert(3,2)](t+t_del))
            
            t_del = delay(data,[id['-1'],id['-2']],t,para=para,delay_on=delay_on)
            P.append(-y[id_convert(2,3)](t+t_del))
            
            t_del = delay(data,[id['3'],id['1'],id['-1']],t,para=para,delay_on=delay_on)
            P.append(y[id_convert(1,3)](t+t_del))
            
            t_del = delay(data,[id['-2'],id['1'],id['-1']],t,para=para,delay_on=delay_on)
            P.append(-y[id_convert(1,2)](t+t_del))
            
            if z!=False:
                
                t_del = delay(data,[id['-2'],id['3']],t,para=para,delay_on=delay_on)
                P.append(0.5*z[id_convert(2,1)](t+t_del))
                
                t_del = delay(data,[id['-2'],id['3']],t,para=para,delay_on=delay_on)
                P.append(-0.5*z[id_convert(3,1)](t+t_del))
                
                t_del = delay(data,[id['1'],id['-1'],id['-2'],id['3']],t,para=para,delay_on=delay_on)
                P.append(-0.5*z[id_convert(2,1)](t+t_del))
                
                t_del = delay(data,[id['1'],id['-1'],id['-2'],id['3']],t,para=para,delay_on=delay_on)
                P.append(0.5*z[id_convert(3,1)](t+t_del))
                
                
                t_del = delay(data,[id['-2']],t,para=para,delay_on=delay_on)
                P.append(0.5*z[id_convert(3,2)](t+t_del))
                
                t_del = delay(data,[id['-2']],t,para=para,delay_on=delay_on)
                P.append(-0.5*z[id_convert(1,2)](t+t_del))
                
                t_del = delay(data,[id['1'],id['-1'],id['-2']],t,para=para,delay_on=delay_on)
                P.append(-0.5*z[id_convert(3,2)](t+t_del))
                
                t_del = delay(data,[id['1'],id['-1'],id['-2']],t,para=para,delay_on=delay_on)
                P.append(0.5*z[id_convert(1,2)](t+t_del))
                
                
                
                t_del = delay(data,[id['3']],t,para=para,delay_on=delay_on)
                P.append(0.5*z[id_convert(1,3)](t+t_del))
                
                t_del = delay(data,[id['3']],t,para=para,delay_on=delay_on)
                P.append(-0.5*z[id_convert(2,3)](t+t_del))
                
                t_del = delay(data,[id['1'],id['-1'],id['3']],t,para=para,delay_on=delay_on)
                P.append(-0.5*z[id_convert(1,3)](t+t_del))
                
                t_del = delay(data,[id['1'],id['-1'],id['3']],t,para=para,delay_on=delay_on)
                P.append(0.5*z[id_convert(2,3)](t+t_del))
            
            
            P = np.array(P)
            
            return sum(np.real(P))
            
            
            
        def obtain_tdi_X_adjusted(data,id,para,y,z,t):
            X=[]
            
            t_del = delay(data,[],t,para=para,delay_on=delay_on)
            X.append(y[id_convert(3,1)](t+t_del))
            
            t_del = delay(data,[id['2']],t,para=para,delay_on=delay_on)
            X.append(y[id_convert(1,3)](t+t_del))
            
            t_del = delay(data,[id['-2'],id['2']],t,para=para,delay_on=delay_on)
            X.append(y[id_convert(2,1)](t+t_del))
            
            t_del = delay(data,[id['2'],id['-2'],id['-3']],t,para=para,delay_on=delay_on)
            X.append(y[id_convert(1,2)](t+t_del))
            
            t_del = delay(data,[id['3'],id['-3'],id['-2'],id['2']],t,para=para,delay_on=delay_on)
            X.append(y[id_convert(2,1)](t+t_del))
            
            t_del = delay(data,[id['3'],id['-3'],id['-2'],id['2'],id['-3']],t,para=para,delay_on=delay_on)
            X.append(y[id_convert(1,2)](t+t_del))
            
            t_del = delay(data,[id['3'],id['-3'],id['3'],id['-3'],id['-2'],id['2']],t,para=para,delay_on=delay_on)
            X.append(y[id_convert(3,1)](t+t_del))
            
            t_del = delay(data,[id['3'],id['-3'],id['3'],id['-3'],id['-2'],id['2'],id['2']],t,para=para,delay_on=delay_on)
            X.append(y[id_convert(1,3)](t+t_del))
            
            t_del = delay(data,[],t,para=para,delay_on=delay_on)
            X.append(-y[id_convert(2,1)](t+t_del))
            
            t_del = delay(data,[id['-3']],t,para=para,delay_on=delay_on)
            X.append(-y[id_convert(1,2)](t+t_del))
            
            t_del = delay(data,[id['3'],id['-3']],t,para=para,delay_on=delay_on)
            X.append(-y[id_convert(3,1)](t+t_del))
            
            t_del = delay(data,[id['3'],id['-3'],id['2']],t,para=para,delay_on=delay_on)
            X.append(-y[id_convert(1,3)](t+t_del))
            
            t_del = delay(data,[id['-2'],id['2'],id['3'],id['-3']],t,para=para,delay_on=delay_on)
            X.append(-y[id_convert(3,1)](t+t_del))
            
            t_del = delay(data,[id['-2'],id['2'],id['3'],id['-3'],id['2']],t,para=para,delay_on=delay_on)
            X.append(-y[id_convert(1,3)](t+t_del))
            
            t_del = delay(data,[id['-2'],id['2'],id['-2'],id['2'],id['2'],id['3'],id['-3']],t,para=para,delay_on=delay_on)
            X.append(-y[id_convert(2,1)](t+t_del))
            
            t_del = delay(data,[id['-2'],id['2'],id['-2'],id['2'],id['2'],id['3'],id['-3'],id['-3']],t,para=para,delay_on=delay_on)
            X.append(-y[id_convert(1,2)](t+t_del))
            
            if z!=False:
                t_del = delay(data,[],t,para=para,delay_on=delay_on)
                X.append(0.5*z[id_convert(2,1)](t-t_del))

                t_del = delay(data,[],t,para=para,delay_on=delay_on)
                X.append(-0.5*z[id_convert(3,1)](t-t_del))

                t_del = delay(data,[id['3'],id['-3']],t,para=para,delay_on=delay_on)
                X.append(-0.5*z[id_convert(2,1)](t-t_del))

                t_del = delay(data,[id['3'],id['-3']],t,para=para,delay_on=delay_on)
                X.append(0.5*z[id_convert(3,1)](t-t_del))

                t_del = delay(data,[id['-2'],id['2']],t,para=para,delay_on=delay_on)
                X.append(-0.5*z[id_convert(2,1)](t-t_del))

                t_del = delay(data,[id['-2'],id['2']],t,para=para,delay_on=delay_on)
                X.append(0.5*z[id_convert(3,1)](t-t_del))

                t_del = delay(data,[id['3'],id['-3'],id['3'],id['-3'],id['-2'],id['2']],t,para=para,delay_on=delay_on)
                X.append(0.5*z[id_convert(2,1)](t-t_del))

                t_del = delay(data,[id['3'],id['-3'],id['3'],id['-3'],id['-2'],id['2']],t,para=para,delay_on=delay_on)
                X.append(-0.5*z[id_convert(3,1)](t-t_del))

                t_del = delay(data,[id['-2'],id['2'],id['-2'],id['2'],id['3'],id['-3']],t,para=para,delay_on=delay_on)
                X.append(0.5*z[id_convert(2,1)](t-t_del))

                t_del = delay(data,[id['-2'],id['2'],id['-2'],id['2'],id['3'],id['-3']],t,para=para,delay_on=delay_on)
                X.append(-0.5*z[id_convert(3,1)](t-t_del))

                t_del = delay(data,[id['-2'],id['2'],id['3'],id['-3'],id['3'],id['-3'],id['-2'],id['2']],t,para=para,delay_on=delay_on)
                X.append(-0.5*z[id_convert(2,1)](t-t_del))

                t_del = delay(data,[id['-2'],id['2'],id['3'],id['-3'],id['3'],id['-3'],id['-2'],id['2']],t,para=para,delay_on=delay_on)
                X.append(0.5*z[id_convert(3,1)](t-t_del))

            X = np.array(X)
            
            return sum(np.real(X))
            
            
        def obtain_tdi(data,id,para,y,z,t):
            X=[]
            
            t_del = delay(data,[id['3'],id['2'],id['-2']],t,para=para,delay_on=delay_on)
            X.append(y[id_convert(1,2)](t+t_del))
            
            t_del = delay(data,[id['2'],id['-2']],t,para=para,delay_on=delay_on)
            X.append(y[id_convert(2,1)](t+t_del))
            
            t_del = delay(data,[id['-2']],t,para=para,delay_on=delay_on)
            X.append(y[id_convert(1,3)](t+t_del))
            
            t_del = delay(data,[],t,para=para,delay_on=delay_on)
            X.append(y[id_convert(3,1)](t+t_del))
            
            t_del = delay(data,[id['-2'],id['-3'],id['3']],t,para=para,delay_on=delay_on)
            X.append(-y[id_convert(1,3)](t+t_del))
            
            t_del = delay(data,[id['-3'],id['3']],t,para=para,delay_on=delay_on)
            X.append(-y[id_convert(3,1)](t+t_del))
            
            t_del = delay(data,[id['3']],t,para=para,delay_on=delay_on)
            X.append(-y[id_convert(1,2)](t+t_del))
            
            t_del = delay(data,[],t,para=para,delay_on=delay_on)
            X.append(-y[id_convert(2,1)](t+t_del))
            
            
            if z!=False:
                t_del = delay(data,[id['2'],id['-2'],id['-3'],id['3']],t,para=para,delay_on=delay_on)
                X.append(-0.5*z[id_convert(3,1)](t-t_del))

                t_del = delay(data,[id['-3'],id['3'],id['2'],id['-2']],t,para=para,delay_on=delay_on)
                X.append(0.5*z[id_convert(2,1)](t-t_del))

                t_del = delay(data,[],t,para=para,delay_on=delay_on)
                X.append(-0.5*z[id_convert(3,1)](t-t_del))

                t_del = delay(data,[],t,para=para,delay_on=delay_on)
                X.append(0.5*z[id_convert(2,1)](t-t_del))

                t_del = delay(data,[id['2'],id['-2']],t,para=para,delay_on=delay_on)
                X.append(0.5*z[id_convert(3,1)](t-t_del))

                t_del = delay(data,[id['2'],id['-2']],t,para=para,delay_on=delay_on)
                X.append(-0.5*z[id_convert(2,1)](t-t_del))

                t_del = delay(data,[id['-3'],id['3']],t,para=para,delay_on=delay_on)
                X.append(-0.5*z[id_convert(3,1)](t-t_del))

                t_del = delay(data,[id['-3'],id['3']],t,para=para,delay_on=delay_on)
                X.append(0.5*z[id_convert(2,1)](t-t_del))
        
            X = np.array(X)
            
            return sum(np.real(X))
               
        def obtain_X_val(self,data,t,y_all,z_all,y_all_include,para,method='XYZ'):
                   
            X_val = 0 
            if y_all_include=='all':
                y_all_include = y_all.keys()
            for source in y_all_include:
                if method=='XYZ':
                   X_val = X_val + obtain_tdi(data,id,para,y_all[source],z_all[source],t)
                elif method=='aby':
                   X_val = X_val + obtain_tdi_aby(data,id,para,y_all[source],t)
                elif method=='X1Y1Z1':
                   X_val = X_val + obtain_tdi_X_adjusted(data,id,para,y_all[source],z_all[source],t)
                elif method=='PQR':
                   X_val = X_val + obtain_tdi_P(data,id,para,y_all[source],z_all[source],t)
            
            
            return X_val
        
        X_ret=[]
        X_func = [] 
        for t in t_vec:
            X_ret.append(obtain_X_val(self,data,t,y_all,z_all,y_all_include,para,method=method))
        X_ret =  np.array(X_ret)
        X_func = lambda time: obtain_X_val(self,data,time,y_all,z_all,y_all_include,para,method=method)
        
        return X_ret,X_func
  
        
        
        
        




            


