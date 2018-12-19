def phi_gauss(self,i,t,dX,dY,side='default',print0 = False):
    dx_list = self.xlist
    dy_list = self.ylist
    Ndata = self.Ndata

    [i_self,i_left,i_right] = PAA_LISA.utils.i_slr(i)

    n = Ndata.data.n_func(i,t)
    if print0==True:
        print ('Type of telescope control is: ' + self.tele_control)

    if side=='default':
        side = self.side

    if side=='l':
        if self.tele_control=='full control':
            tele = Ndata.tele_l_fc(i,t)
        elif self.tele_control=='no control':
            tele = Ndata.tele_l(i,t)
        elif self.tele_control=='SS':
            r = Ndata.data.r_func(i,t)
            ang_SS = self.tele_SS_l(i,t)
            tele = LA.unit(LA.rotate(r,n,ang_SS))

        beam = Ndata.data.u_l_func_tot(i,t) #...adjust for sending telescope --> PAAM control
    elif side=='r':
        if self.tele_control=='full control':
            tele = Ndata.tele_r_fc(i,t)
        elif self.tele_control=='no control':
            tele = Ndata.tele_r(i,t)
        elif self.tele_control=='SS':
            r = Ndata.data.r_func(i,t)
            ang_SS = self.tele_SS_r(i,t)
            tele = LA.unit(LA.rotate(r,n,ang_SS))
        beam = Ndata.data.u_r_func_tot(i,t)

    beam = self.scale*beam
    [ang_x,ang_y] = LA.beam_ang(beam,tele,n)

    psi = 0 #Gouy angle
    self.phasefront_send(i,t,side=side) #...adjust: with delay

    Zm = self.ttl_s_tot #Exit pupil abberation
    if self.jitter[0]!=False:
        jitter = self.jitter[0][1](t)
    else:
        jitter = 0

    dX_ac = np.cos(ang_x)*dX
    dY_ac = np.cos(ang_y)*dY
    dZ_ac = np.sin(ang_x)*dX+np.sin(ang_y)*dY
    dR_ac = (dX_ac**2+dY_ac**2)**0.5
    d_ac = self.z_solve(dX_ac,dY_ac,np.linalg.norm(beam)+dZ_ac)


    Deltax = self.Deltax
    Deltay = self.Deltay
    Nbinsx = self.Nbinsx
    Nbinsy = self.Nbinsy

    if self.speed_on==True:
        dx_list = [0]
        dy_list = [0]
        ps = np.empty((1,1),dtype = np.complex64)
    else:
        ps = np.empty((Nbinsx,Nbinsy),dtype = np.complex64)

    for i in range(0,len(dx_list)):
        for j in range(0,len(dy_list)):
            dx = dx_list[i]
            dy = dy_list[j]
            dr = ((dx**2)+(dy**2))**0.5
            if dr<=D:
                Zm_calc = Zm[i,j]+jitter
                s = (d_ac**2+(dR_ac-dr)**2)**0.5
                print(s)
                E = E0*np.exp((-(dr**2))/(w0**2))*np.exp(1j*psi)
                #print(E,E0,dr,psi)
                R_new = d_ac
                #R_new = np.linalg.norm(vec)
                #print(np.exp(1j*((2*np.pi)/labda)*(Zm+s))/s)
                ps[i,j] = (E*np.exp(1j*((2*np.pi)/labda)*(Zm[i,j]+s))/s)*Deltax*Deltay*np.exp((-1j*2*np.pi*R_new)/labda)
            else:
                ps[i,j] = np.nan

    diff_inte = np.nansum(ps)
    phi = np.angle(diff_inte)
    diff_inte_mag = (diff_inte.real**2+diff_inte.imag**2)**0.5
    #phi = np.log(diff_inte_new/diff_inte_mag).imag

    return diff_inte_mag,phi,ps

