from imports import *
from functions import *
from parameters import *

class Noise():
    def __init__(self,**kwargs):
        self.data = kwargs.pop('data',False)
        self.t_all = self.data.t_plot[0]

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

    def Noise_time(self,f0,f_max,N,psd,t_stop,unit='freq',t=False):
        '''Returns the noise in time domain untill t_stop (or further)'''
        
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

        return [t_ret,noise_ret],func_noise

    def lasernoise(self,PSD,f0=1e-6,f_max=1e-3,N=4096):
        '''Obatains lasernoise''' #..adjust for phase locking
                    
        def y_laser(data,C,C_star,t):
            y={}
            for s in range(1,4):
                for r in range(1,4):
                    check=True
                    if (s - r)%3 ==1:
                        sign = 'pos'
                        t_del = data.L_rl_func_tot[r-1](t)
            
                    elif (s - r)%3 ==2:
                        sign = 'neg'
                        t_del = data.L_rr_func_tot[r-1](t)
                    
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
            t_all = self.t_all
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
            t_stop = self.data.t_plot[i][-1]
            [t_ret_calc,noise_ret],func_noise = self.Noise_time(f0,f_max,N,PSD,t_stop)
            nu_0 = c/labda

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
        data = self.data

        P_v_l=[]
        P_v_l0=[]
        P_v_r=[]
        P_v_r0=[]
        phi_l_sn_q = []
        phi_r_sn_q = []

        for i in range(0,len(data.t_plot)):
            i_self = i+1
            t = data.t_plot[i]
            v_l_func = data.v_l_func_tot[i]
            v_r_func = data.v_r_func_tot[i]
            u_l_func = data.u_l_func_tot[i]
            u_r_func = data.u_r_func_tot[i]
            P_v_l_vec=[]
            P_v_l0_vec=[]
            P_v_r_vec=[]
            P_v_r0_vec=[]
            phi_l_sn_q_vec = []
            phi_r_sn_q_vec = []
            
            for j in range(0,len(t)):
                v_l = v_l_func(t[i])
                L = np.linalg.norm(v_l)
                ang_v_l_in = data.PAA_beam_next_sc[0][i][j] - np.radians(30)
                ang_v_l_out = data.PAA_beam_next_sc[1][i][j] - np.radians(30)
                P_rec0 = ((np.pi**2)*(D**4)*(eta_opt*P_L))/(16*(labda**2)*(L**2))
                P_rec = P_rec0*np.cos(ang_v_l_in)*np.cos(ang_v_l_out)
                P_v_l_vec.append(P_rec)
                P_v_l0_vec.append(P_rec0)

                phi_l_sn_q_vec.append((h/(eta_pd*(P_rec/4.0)))**0.5) # Shot noise per quadrant (upper limit)
                    
                v_r = v_r_func(t[i])
                L = np.linalg.norm(v_r)
                ang_v_r_in = data.PAA_beam_next_sc[2][i][j] - np.radians(30)
                ang_v_r_out = data.PAA_beam_next_sc[3][i][j] - np.radians(30)
                P_rec0 = ((np.pi**2)*(D**4)*(eta_opt*P_L))/(16*(labda**2)*(L**2))
                P_rec = P_rec0*np.cos(ang_v_r_in)*np.cos(ang_v_r_out)
                P_v_r_vec.append(P_rec)
                P_v_r0_vec.append(P_rec0)

                phi_r_sn_q_vec.append((h/(eta_pd*(P_rec/4.0)))**0.5) # Shot noise per quadrant (upper limit)
            
            P_v_l.append(np.array(P_v_l_vec))
            P_v_l0.append(np.array(P_v_l0_vec))
            P_v_r.append(np.array(P_v_r_vec))
            P_v_r0.append(np.array(P_v_r0_vec))
            phi_l_sn_q.append(np.array(phi_l_sn_q_vec))
            phi_r_sn_q.append(np.array(phi_r_sn_q_vec))
            
        P_r_l_func = {}
        P_r_r_func={}
        P_r_l_max = {}
        P_r_r_max = {}

        for i in range(0,3):
            P_r_l_func[str(i+1)] = interp1d(t,P_v_l0[i],bounds_error=False)
            P_r_r_func[str(i+1)] = interp1d(t,P_v_r0[i],bounds_error=False)
            P_r_l_max[str(i+1)] = max(P_v_l0[i])
            P_r_r_max[str(i+1)] = max(P_v_r0[i])
                

        def y_shot(data,P_r_l_max,P_r_r_max,t):
            y={}
            # ... Using the maximum value of P_r_r(/P_r_l) to cancel time dependance
            for s in range(1,4):
                for r in range(1,4):
                    check=True
                    if (s - r)%3 ==1:
                        P_r = P_r_l_max[str(r)]
            
                    elif (s - r)%3 ==2:
                        P_r = P_r_r_max[str(r)]
                    
                    elif s==r:
                        check=False
                    hcbar = 3.16152649e-26
                    if check==True:
                        psd = lambda f: ((hcbar)/(2*np.pi*labda*P_r))*1e-18 #...adjust
                        
                        [t_ret,noise_ret],func_noise = self.Noise_time(1e-8,1e-5,256,psd,data.t_plot[r-1][-1])
                    
                        y[str(s)+str(r)] = func_noise(t)
            
            return y,psd

        def yz_shot_calc(data,P_r_l_max,P_r_r_max,yorz):
            t_all = self.t_all
            y_ret={}
            y_ret_interp = {}
            for t in t_all:
            
                if yorz=='y':
                    y = y_shot(data,P_r_l_max,P_r_r_max,t)[0]
                elif yorz=='z':
                    y = False
                    break

                for keys in y.keys():
                    if keys not in y_ret.keys():
                        y_ret[keys] = []
                    y_ret[keys].append(y[keys])
            
            for keys in y_ret:
                y_ret_interp[keys] = interpolate(t_all,y_ret[keys])
            
            if yorz=='y':
                return y_ret_interp
            elif yorz=='z':
                return False
        
        self.y_shotnoise = yz_shot_calc(self.data,P_r_l_max,P_r_r_max,'y')
        self.z_shotnoise = yz_shot_calc(self.data,P_r_l_max,P_r_r_max,'z')











class TDI():
    def __init__(self,**kwargs):
        print('TDI')
        noise = kwargs.pop('noise',False)
        y_all_select = kwargs.pop('y_all_select','all')
        data = noise.data
        t_vec = noise.t_all
        y_all = noise.noise_y
        z_all = noise.noise_z
        
        para = ['X','Y','Z','alpha','beta','gamma','X1','Y1','Z1','P','Q','R']
        tdi_all = {}
        tdi_all['all']={}
        count=0
        for keys in y_all.keys():
            count = count+1
            tdi_all[keys]={}
            for i in range(0,len(para)):
                t,X_calc= self.tdi_val(data,t_vec,y_all,z_all,y_all_include=[keys],para=para[i])
                func = lambda time: interpolate(t,X_calc)
                tdi_all[keys][para[i]] = [np.array(t),np.array(X_calc),func]

            if count==1:
                t,X_calc= self.tdi_val(data,t_vec,y_all,z_all,y_all_include='all',para=para[i])
                func = interpolate(t,X_calc)
                tdi_all['all'][para[i]] = [np.array(t),np.array(X_calc),func]

        self.tdi_all = tdi_all


    def tdi_val(self,data,t_vec,y_all,z_all,y_all_include='all',side='l',generation=0,para='X',print_on=False,delay_on=True):
        
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
        
        X_ret=[]
        t_ret=[]
        for t in t_vec:
            try:
                X = 0
                if y_all_include=='all':
                    y_all_include = y_all.keys()
                for source in y_all_include:
                    if method=='XYZ':
                        X = X + obtain_tdi(data,id,para,y_all[source],z_all[source],t)
                    elif method=='aby':
                        X = X + obtain_tdi_aby(data,id,para,y_all[source],t)
                    elif method=='X1Y1Z1':
                        X = X + obtain_tdi_X_adjusted(data,id,para,y_all[source],z_all[source],t)
                    elif method=='PQR':
                        X = X + obtain_tdi_P(data,id,para,y_all[source],z_all[source],t)
                
                X_ret.append(X)
                t_ret.append(t)
                #print(X_ret[-1])
            except ValueError:
                #print "Unexpected error:", sys.exc_info()[0]
                #print(t)
                #raise ValueError('can not calculate at t='+str(t))
                pass
            
        return np.array(t_ret),np.array(X_ret)
        
        
        




            


