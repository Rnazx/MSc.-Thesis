import numpy as np
import matplotlib.pyplot as plt
from sympy import *
from fractions import Fraction
import pickle
import sys

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


##################defining symbols#######################################
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

##############################################################################################################

with open('turb_exp.pickle', 'rb') as f:
     hg, rho, nu, u, l, taue, taur = pickle.load(f)


alphak1 = calpha*tau**2*u**2*omega/h
alphak2 = calpha*tau*u**2/h
alphak3 = kalpha*u

alphareg = int(sys.argv[1])

if alphareg == 1:
    alphak = alphak1
elif alphareg == 2:
    alphak = alphak2
else :
    alphak = alphak3

Beq = bet*u*(4*pi*rho)**Rational(1/2)
biso = (Beq*(xio**(1/2)))/mach
biso = simplify(biso)
biso = biso.powsimp(force=True)


bani = biso*(Rational(1/3)*2*q*omega*tau*(1+(q*omega*tau)/2))**Rational(1/2)  #+ (Uo*tau/l)*(1+1/(1+q*omega*tau)**2)
bani = simplify(bani)
bani = bani.powsimp(force=True)

Rk = Symbol('R_k')
eta = (1/3)*tau*u**2

Ralpha = alphak*h/eta
Romega = -q*omega*h**2/eta
Dk = Ralpha*Romega
Dc = -(pi**5)/32
Bbar = (pi*Beq*l*(Rk*(Dk/Dc-1))**(0.5))/h
Bbar = simplify(Bbar)
#Bbar = Bbar.powsimp(force=True)

tanpB = -((pi**2)*tau*(u**2))/(12*q*omega*(h**2))
tanpB = simplify(tanpB)
tanpB = tanpB.subs([(tau, tau), (l, l)])
tanpB = simplify(tanpB)

tanpb = 1/(1+q*omega*tau)

mag_expr = biso, bani, Bbar, tanpb, tanpB, Beq, eta, alphak


with open('mag_exp.pickle', 'wb') as f:
    pickle.dump(mag_expr, f)
