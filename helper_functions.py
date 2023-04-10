import numpy as np
from sympy import *
import inspect
from scipy.optimize import curve_fit, fsolve

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
alphak = Symbol('alpha_k')

# Defining the general parameters
u = Symbol('u')
tau = Symbol('tau')
l = Symbol('l')
h = Symbol('h')

cs = (gamma*boltz*T/(mu*mh))**Rational(1/2)


def power_law(x, a, b):
    return a*np.power(x, b)

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

    const = [(gamma, 1.5), (boltz, 1.3807e-16), (mh, 1.67e-24), (mach, 1), (G, 6.67e-8),
            (cl, 3/4), (xio, 0.4), (mstar, 0.85*g_Msun), (delta, 8e-3), (E51, 1), (kalpha, 1)]
  
    express = express.subs(const).simplify(force=True) 

    exp = np.array([express.evalf(subs={ sigmatot:sigt, sigma: sig, sigmasfr: sigsfr, q:qs, omega: oms, zet:zets, T:t, psi:ps, bet :b, calpha: ca, Rk: rk, mu: m}) for sigt,sig, sigsfr,qs, oms, zets, t, ps, b, ca, rk, m in data_pass])
    
    return exp

############################################################################################################################

def datamaker(quan, data_pass, h_f, tau_f = None, alphak_f = None, scal_rel = False):
    quan_val = exp_analytical_data(quan, data_pass)
    if tau_f is None: 
        tau_f = np.ones(len(h_f))
    if alphak_f is None:
        return  np.array([np.float64( quan_val[i].evalf(subs= {h : hf, tau : tauf} )) for i, (hf, tauf) in enumerate(zip(h_f, tau_f))])
    else: 
        Bbar_in = np.array([quan_val[i].evalf(subs= {h : hf, tau : tauf, alphak : alphakf} ) for i, (hf, tauf, alphakf) in enumerate(zip(h_f, tau_f, alphak_f))])
        return np.float64(Bbar_in*(np.float64(Bbar_in*Bbar_in>0)))
        
    
##############################################################################################################################################

def root_finder(h_val):
    h_f = []
    for hv in h_val:
        func = lambda x : np.array([np.float64((h-hv).evalf(subs={h : i})) for i in x])
        h_initial_guess = 7e+25
        h_solution = fsolve(func, h_initial_guess)
        h_f.append(h_solution[0])
    h_f = np.array(h_f)

    return h_f

#####################################################################################################
def scal_helper(express, data_pass, observable = zet, _range = np.linspace(1,5000,50)):

    const = [(gamma, 1.5), (boltz, 1.3807e-16), (mh, 1.67e-24), (mach, 1), (G, 6.67e-8),
            (cl, 3/4), (xio, 0.4), (mstar, 0.85*g_Msun), (delta, 8e-3), (E51, 1), (kalpha, 1)]
  
    express = express.subs(const).simplify(force=True)
    sigt,sig, sigsfr,qs, oms, zets, t, ps, b, ca, rk, m = data_pass[0]
    val_subs = { sigmatot:sigt, sigma: sig, sigmasfr: sigsfr, q:qs, omega: oms, zet:zets, T:t, psi:ps, bet :b, calpha: ca, Rk: rk, mu: m}
    try:
        val_subs.pop(observable)
        obs_val = _range#(val_subs.pop(observable))*_range
    except ValueError:
        print('Observable does not exist!')
    exp = express.evalf(subs = val_subs)
    return obs_val, np.array([exp.evalf(subs = {observable:o}) for o in obs_val])


def scal_finder(h_exp, tau_exp, quan_exp, observable, data_pass, _range = np.linspace(1,5000,50)):
    obs_val, h_val = scal_helper(h_exp, data_pass, observable, _range)
    h_scal = root_finder(h_val)
    obs_val, tau_val = scal_helper(tau_exp, data_pass, observable, _range)
    tau_scal = np.array([np.float64( tau_val[i].evalf(subs= {h : hf} )) for i, hf in enumerate(h_scal)]) 

    obs_val, quan_val = scal_helper(quan_exp, data_pass, observable, _range)
    quan_f = np.array([np.float64( quan_val[i].evalf(subs= {h : hf, tau : tauf} )) for i, (hf, tauf) in enumerate(zip(h_scal, tau_scal))])

    pg, cov = curve_fit(f=power_law, xdata=obs_val, ydata=quan_f, p0=[0, 0], bounds=(-np.inf, np.inf))

    return quan_f, round(pg[1],3)

