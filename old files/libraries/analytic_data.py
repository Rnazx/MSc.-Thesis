import numpy as np
import matplotlib.pyplot as plt
from sympy import *
from fractions import Fraction
import pickle
import os

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

def analytic_data(kms_dat_v, dat_sigma, dat_omega, dat_q, dat_sigmasfr, model_no =3, let = 'a'):
    kms = 1e+5
    kpcm = 3.086e+21
    pcm = kpcm/1e+3
    Msun = 1.989e+33
    cs = (gamma*boltz*T/(mu*mh))**Rational(1/2)
    cms_cs = cs.evalf(subs={gamma: 1.5, boltz: 1.38e-16, T:1e+4, mu:14/11, mh:1.67e-24})
    print('Value of cs is {} and the max value of u is {}'.format(cms_cs,(kms_dat_v*kms).max()))
    let = 'a'
    if np.sqrt(2)*cms_cs < (kms_dat_v*kms).max(): 
        print('$M>>\sqrt(2)$')
        let = 'b'
    else : 
        print('$M<<\sqrt(2)$')
    script_dir = os.getcwd()  # <-- absolute dir the script is in
    a = [(gamma, 1.5), (boltz, 1.3807e-16), (mh, 1.67e-24), (mu, 14/11), (mach, 1), (G, 6.67e-8),
            (cl, 3/4), (xio, 0.4), (calpha, 1), (Rk, 0.3), (mstar, 0.85*Msun), (delta, 8e-3), (E51, 1),
            (sigmatot, 10*Msun/(pcm)**2), (sigma, dat_sigma.mean()),
                (sigmasfr, dat_sigma.mean()), (omega, dat_omega.mean()), (q, dat_q.mean()), (T, 1e+4)]
    def model(let, model_no, a):
        model_name = "\model"+str(model_no)+let+".txt"
        rel_path = "model_scripts" + model_name
        abs_file_path = os.path.join(script_dir, rel_path)
        with open(abs_file_path, "rb") as inf:
            quantities = pickle.load(inf)
        if quantities[3].subs(a) < simplify(quantities[1].subs(a)/quantities[2].subs(a)):
            return quantities, let
        else:
            
            if let =='a':
                let = 'c'
                print('$tau^e>tau^r$. Therefore model changed to '+str(model_no) +let)
                return model(let, model_no, a)
            elif let =='b':
                let = 'd'
                print('$tau^e>tau^r$. Therefore model changed to '+str(model_no) +let)
                return model(let, model_no, a) 
            

    quantities, let = model(let, model_no, a)
    vel = quantities[2]
    const = [(gamma, 1.5), (boltz, 1.3807e-16), (mh, 1.67e-24), (mu, 14/11), (mach, 1), (G, 6.67e-8),
         (cl, 3/4), (xio, 0.4), (calpha, 1), (Rk, 0.3), (mstar, 0.85*Msun), (delta, 8e-3), (E51, 1)]
    variables = [(sigmatot, 10*Msun/(pcm)**2), (sigma, 1),
                (sigmasfr, 1), (omega, 1), (q, 1), (T, 1e+4)]

    # plotting the scaling relations
    observ = [sigma, sigmasfr, q, omega]

    for obs in observ:
        variables.remove((obs, 1))
    final = const + variables
    z = vel.subs(final).simplify(force=True)
    scalreldata = []
    for sig, sigsfr in zip(dat_sigma, dat_sigmasfr):
        scalreldata.append(z.evalf(subs={sigma: sig, sigmasfr: sigsfr}))
    scalreldata = np.array(np.float64(scalreldata))/kms
    return scalreldata


