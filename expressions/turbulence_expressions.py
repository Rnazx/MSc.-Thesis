import numpy as np
import matplotlib.pyplot as plt
from sympy import *
from fractions import Fraction
import pickle
import sys

pc_kpc = 1e3  # number of pc in one kpc
cm_km = 1e5  # number of cm in one km
s_day = 24*3600  # number of seconds in one day
s_min = 60  # number of seconds in one hour
s_hr = 3600  # number of seconds in one hour
cm_Rsun = 6.957e10  # solar radius in cm
g_Msun = 1.989e33  # solar mass in g
cgs_G = 6.674e-8
cms_c = 2.998e10
g_mH = 1.6736e-24
g_me = 9.10938e-28
cgs_h = 6.626e-27
deg_rad = 180e0/np.pi
arcmin_deg = 60e0
arcsec_deg = 3600e0
cm_kpc = 3.086e+21  # number of ccentimeters in one parsec
s_Myr = 1e+6*(365*24*60*60)  # megayears to seconds


################## defining symbols#######################################
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

##############################################################################################################

# Defining the expressions
cs = (gamma*boltz*T/(mu*mh))**Rational(1/2) #sound speed, eq 36
rho = sigma/(2*h) #gas density, eq 35
n = rho/(mu*mh) #converting mass density rho to number density n, eq 22

#model 2
# lsn = psi*cl*h #lsn= driving scale of isolated SNe, psi=fixed parameter used since u isnt same as velocity dispersion
# l = lsn

#model 3
lsn = psi*0.14*cm_kpc*(E51)**Fraction(16, 51) * \
    (n/0.1)**Fraction(-19, 51)*(cs/(cm_km*10))**Fraction(-1, 3)
l = ((Gamma-1)/Gamma)*cl*lsn

l = simplify(l)

nu = (delta*sigmasfr)/(2*h*mstar)
u = simplify(((4*pi/3)*l*lsn**3*cs**2*nu)**Fraction(1, 3))
#u = cs #assumption in models 1 and 2

#scale height expressions 
hg = zet*(u**2 + cs**2)/(3*pi*G*sigmatot)
hsub = zet*(cs**2)/(3*pi*G*sigmatot)
hsup = zet*(u**2)/(3*pi*G*sigmatot)


rho = sigma/(2*h)
n = rho/((14/11)*mh)
taue = simplify(l/u) #eq 30
taur = simplify(6.8*s_Myr*(1/4)*(nu*cm_kpc**3*s_Myr/50)**(-1)*(E51)
                ** Fraction(-16, 17) * (n/0.1)**Fraction(19, 17)*(cs/(cm_km*10))) #eq 32

#eq 11
alphak1 = calpha*tau**2*u**2*omega/h
alphak2 = calpha*tau*u**2/h
alphak3 = kalpha*u


turb_expr = hg, rho, nu, u, l, taue, taur, alphak1, alphak2, alphak3


with open('turb_exp.pickle', 'wb') as f:
    pickle.dump(turb_expr, f)

print('Solved the turbulence expressions')