import numpy as np
import matplotlib.pyplot as plt
from sympy import *
from matplotlib.ticker import FormatStrFormatter
import inspect

###############################################################################################
pc_kpc = 1e3		#number of pc in one kpc
cm_km = 1e5		#number of cm in one km
s_day = 24*3600		#number of seconds in one day
s_min = 60		#number of seconds in one hour
s_hr = 3600		#number of seconds in one hour
cm_Rsun = 6.957e10	#solar radius in cm
g_Msun = 1.989e33	#solar mass in g
cgs_G = 6.674e-8 
cms_c = 2.998e10
g_mH = 1.6736e-24
g_me = 9.10938e-28
cgs_h = 6.626e-27
deg_rad = 180e0/np.pi
arcmin_deg = 60e0
arcsec_deg = 3600e0
kpc_cm = 3.086e+21  #number of ccentimeters in one parsec
Myr_s = 1e+6*(365*24*60*60) #megayears to seconds

############################################################################################################################
# Defining the Observables
q = Symbol('q')
omega = Symbol('\Omega')
sigma = Symbol('\Sigma')
sigmatot = Symbol('Sigma_tot')
sigmasfr = Symbol('Sigma_SFR')
T = Symbol('T')


# Defining the Constants
calpha = Symbol('C_alpha')
gamma = Symbol('gamma')
boltz = Symbol('k_B')
mu = Symbol('mu')
mh = Symbol('m_H')
G = Symbol('G')
xio = Symbol('xi_0')
delta = Symbol('\delta')
mstar = Symbol('m_*')
cl = Symbol('C_l')
kappa = Symbol('kappa')
mach = Symbol('M')
E51 = Symbol('E_51')
Rk = Symbol('R_k')
zet = Symbol('zeta')
psi = Symbol('psi')
kalpha = Symbol('K_alpha')
bet = Symbol('beta')

# Defining the general parameters
u = Symbol('u')
tau = Symbol('tau')
l = Symbol('l')
h = Symbol('h')

cs = (gamma*boltz*T/(mu*mh))**Rational(1/2)

def list_transpose(x):
    array = np.array(x)
    transposed_array = array.T
    return transposed_array.tolist()
######################################################################################################################
def retrieve_name(var):
    callers_local_vars = inspect.currentframe().f_back.f_locals.items()
    return [var_name for var_name, var_val in callers_local_vars if var_val is var][-1]

###############################################################################################################

def exp_analytical_data(express, data_pass):

    const = [(gamma, 1.5), (boltz, 1.3807e-16), (mh, 1.67e-24), (mu, 14/11), (mach, 1), (G, 6.67e-8),
            (cl, 3/4), (xio, 0.4), (mstar, 0.85*g_Msun), (delta, 8e-3), (E51, 1), (kalpha, 1)]
  
    express = express.subs(const).simplify(force=True) 

    exp = np.array([express.evalf(subs={ sigmatot:sigt, sigma: sig, sigmasfr: sigsfr, q:qs, omega: oms, zet:zets, T:t, psi:ps, bet :b, calpha: ca, Rk: rk}) for sigt,sig, sigsfr,qs, oms, zets, t, ps, b, ca, rk in data_pass])
    
    return exp

############################################################################################################################

def datamaker(quan, data_pass, h_f, tau_f = None):
    quan_val = exp_analytical_data(quan, data_pass)
    if tau_f is None: 
        tau_f = np.ones(len(h_f))
    if retrieve_name(quan) == 'Bbar':
        Bbar_in = np.array([quan_val[i].evalf(subs= {h : hf, tau : tauf} ) for i, (hf, tauf) in enumerate(zip(h_f, tau_f))])
        return np.float64(Bbar_in*(np.float64(Bbar_in*Bbar_in>0)))
    else: 
        return  np.array([np.float64( quan_val[i].evalf(subs= {h : hf, tau : tauf} )) for i, (hf, tauf) in enumerate(zip(h_f, tau_f))])
    
##############################################################################################################################################

def root_finder(h_exp, data_pass):
    h_val = exp_analytical_data(h_exp, data_pass)

    h_f = []
    for hv in h_val:
        func = lambda x : np.array([np.float64((h-hv).evalf(subs={h : i})) for i in x])
        from scipy.optimize import fsolve
        h_initial_guess = 7e+25
        h_solution = fsolve(func, h_initial_guess)
        h_f.append(h_solution[0])
    h_f = np.array(h_f)

    return h_f

#####################################################################################################
def getmean(i,x,y):
    good = np.where((x>kpc_xbnd_M31[i-1]) & (x<=kpc_xbnd_M31[i]))
    xtemp = x[good]
    ytemp = y[good]
    meanval = np.mean(ytemp)
    return xtemp,meanval

def makefig(fn,x,y,xr,yr,xt,yt,l_leg,x2=[],y2=[],lab1='',lab2='',y1e=[],y2e=[],y3=[],nticks=30):
    fig, ax1 = plt.subplots(figsize=(length, breadth))
    ax1.set_xlim(xr)
    ax1.set_ylim(yr)
    plt.errorbar(x2,y2,yerr=y2e,ecolor='k',elinewidth=1,capsize=1,c='tab:blue',mfc='k',mec='k',barsabove=True,label=lab2,marker='D',markersize=1.2)
    plt.plot(x,y,marker='o',markersize=1.2,c='tab:orange',mfc='k',mec='k',label=lab1)
    for i in range(len(kpc_xbnd_M31)):
        plt.plot([kpc_xbnd_M31[i],kpc_xbnd_M31[i]],[0,10*yr[1]],linestyle='dotted',linewidth=1,c='k')
        #plot average values within each bin
        if i>0:
            xtemp,meanval = getmean(i,x,y)
            plt.plot([xtemp[0],xtemp[-1]],[meanval,meanval],c='tab:orange',alpha=0.5)
            if(y2!=[]):
                x2temp,meanval2 = getmean(i,x2,y2)
                plt.plot([x2temp[0],x2temp[-1]],[meanval2,meanval2],c='tab:blue',alpha=0.5)
    if(y3!=[]):
        plt.plot(kpc_xmid_M31,y3,marker='*',mfc='yellow',mec='tab:green',mew=1,linewidth=0,label='Van Eck et al. (2015)')  
    if l_leg:
        handles,labels = ax1.get_legend_handles_labels()
        #change order for legend
        if(y3!=[]):
            if(y2!=[]):
                handles = [handles[0],handles[2],handles[1]]
                labels = [labels[0],labels[2],labels[1]]
        leg= ax1.legend(handles, labels,loc='best',handlelength=4,ncol=1,prop={'size':leg_textsize,'family':'Times New Roman'},fancybox=True,framealpha=0.9,handletextpad=0.7,columnspacing=0.7)
    #
    minor_ticks_x=  np.arange(xr[0],xr[1]+1, 1)	
    minor_ticks_y=  np.arange(yr[0],yr[1]+yr[1]/nticks, yr[1]/nticks)	
    ax1.tick_params(axis='both', which='minor', labelsize=axis_textsize, colors='k', length=3 , width=1   )
    ax1.tick_params(axis='both', which='major', labelsize=axis_textsize, colors='k', length=5 , width=1.25)
    ax1.set_xticks(minor_ticks_x,minor=True)
    ax1.set_yticks(minor_ticks_y,minor=True)
    #
    ax1.set_xlabel(xtit   ,size=axis_textsize)
    ax1.set_xticklabels(ax1.get_xticks(),fontname = "Times New Roman", fontsize=axis_textsize)
    plt.gca().xaxis.set_major_formatter(FormatStrFormatter('%g'))	
        
    ax1.set_ylabel(ytit,size=axis_textsize)
    ax1.set_yticklabels(ax1.get_yticks(),fontname = "Times New Roman", fontsize=axis_textsize)
    plt.gca().yaxis.set_major_formatter(FormatStrFormatter('%g'))	
    
    plt.savefig(fn, format='png', bbox_inches='tight', dpi= 300)
