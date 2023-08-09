import numpy as np
from sympy import *
import inspect
from scipy.optimize import curve_fit, fsolve, root


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
Gamma = Symbol('Gamma')


# Defining the general parameters
u = Symbol('u')
tau = Symbol('tau')
l = Symbol('l')
h = Symbol('h')

cs = (gamma*boltz*T/(mu*mh))**Rational(1/2)

###############################################################################################

g_Msun = 1.989e33  # solar mass in g
cgs_G = 6.674e-8  # gravitational constant in cgs units
g_mH = 1.6736e-24  # mass of hydrogen atom in grams
cgs_kB = 1.3807e-16  # boltzmann constant in cgs units

gval, clval, xioval, mstarval, deltaval, e51val, kaval, Gammaval = tuple(
    np.genfromtxt('constants.in', delimiter='=', dtype=np.float64)[:, -1])

const = [(boltz, cgs_kB), (mh, g_mH), (G, cgs_G), (gamma, gval),
         (cl, clval), (xio, xioval), (mstar, mstarval*g_Msun), (delta, deltaval), (E51, e51val), (kalpha, kaval), (Gamma, Gammaval)]

######################################################################################################################


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

# Function which takes the expression and the data to be substituted
def exp_analytical_data(express, data_pass):
    # Substitute the constants in the given expression
    express = express.subs(const).simplify(force=True)
    # Substitute the data for the observables as well as the parameters for each radii
    exp = np.array([express.evalf(subs={sigmatot: sigt, sigma: sig, sigmasfr: sigsfr, q: qs, omega: oms, zet: zets, T: t,
                   psi: ps, bet: b, calpha: ca, Rk: rk, mu: m}) for sigt, sig, sigsfr, qs, oms, zets, t, ps, b, ca, rk, m in data_pass])

    return exp

############################################################################################################################


def datamaker(quan, data_pass, h_f, tau_f=None, alphak_f=None):
    quan_val = exp_analytical_data(quan, data_pass)
    if tau_f is None:
        tau_f = np.ones(len(h_f))
    if alphak_f is None:
        return np.array([np.float64(quan_val[i].evalf(subs={h: hf, tau: tauf})) for i, (hf, tauf) in enumerate(zip(h_f, tau_f))])
    else:
        Bbar_in = np.array([quan_val[i].evalf(subs={h: hf, tau: tauf, alphak: alphakf}) for i, (
            hf, tauf, alphakf) in enumerate(zip(h_f, tau_f, alphak_f))])
        return np.float64(Bbar_in*(np.float64(Bbar_in*Bbar_in > 0)))


##############################################################################################################################################
# Function which takes the RHS of the expression of h as input (h_val). h_init is the initial guess for the solution in cgs units
def root_finder(h_val, h_init=7e+25):
    #define an empty array
    h_f = []
    for hv in h_val:
        # Function to find the root for
        def func(x): 
            # h is an expression. hv is the RHS of the expression for h for a particular radii
            return np.array(
            [np.float64((h-hv).evalf(subs={h: i})) for i in x])
        # Derivative of the function
        def dfunc(x): 
            return np.array(
            [np.float64(diff((h-hv), h).evalf(subs={h: i})) for i in x])# diff is a symbolic derivative
        # solve for the function using the fsolve routine. First element of this array is the solution
        h_solution = fsolve(func, h_init, fprime=dfunc )
        # append the solution in an array.
        h_f.append(h_solution[0])
    # Convert array to numpy
    h_f = np.array(h_f)
    return h_f

#####################################################################################################


def scal_helper(express, data_pass, observable=zet, _range=np.linspace(1, 5000, 50)):
    express = express.subs(const).simplify(force=True)
    sigt, sig, sigsfr, qs, oms, zets, t, ps, b, ca, rk, m = data_pass[0]
    val_subs = {sigmatot: sigt, sigma: sig, sigmasfr: sigsfr, q: qs,
                omega: oms, zet: zets, T: t, psi: ps, bet: b, calpha: ca, Rk: rk, mu: m}
    try:
        #val_subs.pop(observable)
        obs_val = (val_subs.pop(observable))*_range
    except ValueError:
        print('Observable does not exist!')
    exp = express.evalf(subs=val_subs)
    return obs_val, np.array([exp.evalf(subs={observable: o}) for o in obs_val])


def scal_finder(h_exp, quan_exp, observable, data_pass, tau_exp=None, alpha_exp=None, _range=np.linspace(1, 5000, 200), init_h=7e+20):
    def scal_dat(quan, data_pass, h_f, tau_f=None, alphak_f=None):
        quan_val = scal_helper(quan, data_pass, observable, _range)[1]
        if tau_f is None:
            tau_f = np.ones(len(h_f))
        if alphak_f is None:
            return np.array([np.float64(quan_val[i].evalf(subs={h: hf, tau: tauf})) for i, (hf, tauf) in enumerate(zip(h_f, tau_f))])
        else:
            Bbar_in = np.array([quan_val[i].evalf(subs={h: hf, tau: tauf, alphak: alphakf}) for i, (
            hf, tauf, alphakf) in enumerate(zip(h_f, tau_f, alphak_f))])
            return np.float64(np.abs(Bbar_in))
    obs_val, h_val = scal_helper(h_exp, data_pass, observable, _range)
    h_scal = root_finder(h_val, init_h)
    if tau_exp is not None:
        tau_scal = scal_dat(tau_exp, data_pass, h_scal)
    else:
        tau_scal = None
    if alpha_exp is not None:
        alphak_scal = scal_dat(alpha_exp, data_pass, h_scal, tau_scal)
    else:
        alphak_scal = None
    quan_f = scal_dat(quan_exp, data_pass, h_scal, tau_scal, alphak_scal)
    # pg, cov = curve_fit(f=power_law, xdata=obs_val, ydata=quan_f, p0=np.asarray([10**5,-1]))
    # perr = np.sqrt(np.diag(cov))
    pg = np.mean(
        np.gradient(np.log(np.abs(quan_f)))/np.gradient(np.log(obs_val)))
    return obs_val, quan_f, pg #[1], perr[1]
