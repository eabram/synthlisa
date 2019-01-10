from imports import *
import PAA_LISA
import NOISE_LISA

global LA
LA = PAA_LISA.la()

# Changes of coordinate system
def coor_SC(wfe,i,t):
    # r,n,x (inplane) format

    r = LA.unit(wfe.Ndata.data.r_func(i,t))
    n = LA.unit(wfe.Ndata.data.n_func(i,t))
    x = np.cross(n,r)

    return np.array([r,n,x])

def coor_tele(wfe,i,t,ang_tele,L_tele=2):
    # Retunrs the coordinate system of telescope (same as SC but rotated over ang_tele inplane)

    [r,n,x] = coor_SC(wfe,i,t)
    tele = r*L_tele
    tele = LA.rotate(tele,n,ang_tele)
    r = LA.unit(tele)
    x = np.cross(n,r)

    return np.array([r,n,x])

def beam_tele(wfe,i,t,ang_tele,ang_paam):
    # Retunrs the coordinate system of the transmitted beam (same as SC but rotated over ang_tele inplane and ang_tele outplane)
    [r,n,x] = coor_tele(wfe,i,t,ang_tele) #Telescope coordinate system

    r = LA.unit(LA.rotate(r,x,ang_paam)) # Rotate r in out of plane over ang_paam
    n = np.cross(r,x)

    return np.array([r,n,x])




def i_slr(i):
    i_self = i
    i_left = (i+1)%3
    i_right = (i+2)%3

    i_ret = [i_self,i_left,i_right]
    for j in range(0,len(i_ret)):
        if i_ret[j]==0:
            i_ret[j]=3

    return i_ret

def delay(data,l_array,t,para='X',delay_on=True):
    t_del = 0
    if delay_on==True:
        for k in range(0,len(l_array)):
            j = -1-k
            i_r = (abs(l_array[j])+1)%3
            try:
                if l_array[j]>0:
                    t_del = t_del - data.L_rl[i_r-1](t - t_del)
                elif l_array[j]<0:
                    t_del = t_del - data.L_rr[i_r-1](t - t_del)
            except:
                pass

    return t_del

def PSD(f_list,SD_list):
    return interp1d(f_list,SD_list,bounds_error=False,fill_value=0)

def PowerLaw(SD_val,f0,exp=1):
    return lambda f: (SD_val)*((f/f0)**exp)

def add_func_help(func_list,f):

    func_ret = func_list[0]

    if len(func_list)>1:
        for i in range(1,len(func_list)):
            func_ret = func_ret(f)+func_list[i](f)

    return func_ret

def add_func(func_list):

    return lambda f: add_func_help(func_list,f)

def interpolate(x,y,method='interp1d'):
    if method=='interp1d':
        return interp1d(x,y,bounds_error=False)
    else:
        print('Please select proper interpolation method (interp1d)')


