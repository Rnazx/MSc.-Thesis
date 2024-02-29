import numpy as np
from sympy import *
import inspect
from scipy.optimize import curve_fit, fsolve, root
from scipy.integrate import quad
import os

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
A = Symbol('A')
K = Symbol('K')


# Defining the general parameters
u = Symbol('u')
tau = Symbol('tau')
l = Symbol('l')
h = Symbol('h')

cs = (gamma*boltz*T/(mu*mh))**Rational(1/2)

###############################################################################################
base_path = os.environ.get('MY_PATH')

g_Msun = 1.989e33  # solar mass in g
cgs_G = 6.674e-8  # gravitational constant in cgs units
g_mH = 1.6736e-24  # mass of hydrogen atom in grams
cgs_kB = 1.3807e-16  # boltzmann constant in cgs units

gval, clval, xioval, mstarval, deltaval, e51val, kaval, Gammaval, Rkval = tuple(
    np.genfromtxt(os.path.join(base_path,'inputs','constants.in'), delimiter='=', dtype=np.float64)[:, -1])

const = [(boltz, cgs_kB), (mh, g_mH), (G, cgs_G), (gamma, gval),
         (cl, clval), (xio, xioval), (mstar, mstarval*g_Msun), (delta, deltaval), (E51, e51val), (kalpha, kaval), (Gamma, Gammaval),(Rk, Rkval)]

######################################################################################################################

def parameter_read(filepath):
#opening these files and making them into dictionaries
    params = {}
    with open(filepath, 'r') as FH:
        for file in FH.readlines():
            line = file.strip()
            try:
                par_name, value = line.split('= ')
            except ValueError:
                print("Record: ", line)
                raise Exception(
                    "Failed while unpacking. Not enough arguments to supply.")
            try:
                params[par_name] = np.float64(value)
            except ValueError: #required cz of 14/11 in parameter.in file
                try:
                    num, denom = value.split('/')
                    params[par_name] = np.float64(num) / np.float64(denom)
                except ValueError:
                    params[par_name] = value
            
    return params

###################################################################################################################################
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
                   psi: ps, bet: b, calpha: ca, K: k, mu: m, A:a}) for sigt, sig, qs, oms, sigsfr, t, zets, ps, b, ca, k, m, a in data_pass])

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


####################################################################################################################################################################################

def pitch_angle_integrator(kpc_r, tanpB_f, tanpb_f,Bbar_f, bani_f, tanpB_err,tanpb_err, Bbar_err, bani_err):
    pB = np.arctan(-tanpB_f)
    pB_err = -tanpB_err/(1+tanpB_f**2)
    pb = np.arctan(tanpb_f)
    pb_err = tanpb_err/(1+tanpb_f**2)
    # Old Expression
    # pbo = (1/2)*((1+(2*Bbar_f*bani_f*np.cos(pbb-pB))/
    #             (bani_f**2+Bbar_f**2))*np.arctan((Bbar_f*np.sin(pB) + bani_f*np.sin(pbb))/
    #                                             ((Bbar_f*np.cos(pB)) + bani_f*np.cos(pbb)))
    #                                                         + (1-(2*Bbar_f*bani_f*np.cos(pbb-pB))/
    #                                                         (bani_f**2+Bbar_f**2))*np.arctan((Bbar_f*np.sin(pB) - bani_f*np.sin(pbb))/
    #                                                                                             ((Bbar_f*np.cos(pB)) - bani_f*np.cos(pbb))))

    def pogen(b, B, pb, pB, s):
        return (np.exp(-b**2/(2*s**2))/
                (np.sqrt(2*(np.pi))*s))*(1+(2*B*b*np.cos(pb-pB))/
                            (b**2 + B**2))*np.arctan((B*np.sin(pB) + b*np.sin(pb))/
                                                        ((B*np.cos(pB)) + b*np.cos(pb)))
    brms = np.sqrt(np.average(bani_f**2))

    h = 1e-8
    def dpodbani(b, B, pb, pB, s):
        return (pogen(b, B, pb, pB, s+h)-pogen(b, B, pb, pB, s-h))/(2*h)
    def dpodBbar(b, B, pb, pB, s):
        return (pogen(b, B+h, pb, pB, s)-pogen(b, B-h, pb, pB, s))/(2*h)
    h = 0.01
    def dpodpB(b, B, pb, pB, s):
        return (pogen(b, B, pb, pB+h, s)-pogen(b, B, pb, pB-h, s))/(2*h)
    def dpodpb(b, B, pb, pB, s):
        return (pogen(b, B, pb+h, pB, s)-pogen(b, B, pb-h, pB, s))/(2*h)

    def integrator(fn, interval = 1e+2): #1e+3
        for i in range(len(kpc_r)):
            print(i,Bbar_f[i], pb[i], pB[i], bani_f[i])
        return np.array([quad(fn, -interval, interval, args=(Bbar_f[i], pb[i], pB[i], bani_f[i]),
                points=[-interval*brms, interval*brms])[0] for i in range(len(kpc_r))])
    po = integrator(pogen)

    inte = 1e+3
    po_err = np.array([quad(pogen, -inte, inte, args=(Bbar_f[i], pb[i], pB[i], bani_f[i]),
                points=[-inte*brms, inte*brms])[1] for i in range(len(kpc_r))]) 
    po_err += np.sqrt((integrator(dpodbani,inte)*bani_err)**2 +(integrator(dpodBbar,inte)*Bbar_err)**2
                    +(integrator(dpodpB,inte)*pB_err)**2+(integrator(dpodpb,inte)*pb_err)**2)
    
    return pB, po, pb, pB_err, po_err, pb_err