import numpy as np
import matplotlib.pyplot as plt
from sympy import *
from fractions import Fraction
import pickle
import os


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
Rk = Symbol('R_k')

##############################################################################################################
#U_0=0, f_SB=0 assumed
current_directory = str(os.getcwd())
with open(current_directory+r'\expressions\turb_exp.pickle', 'rb') as f:
    hg, rho, nu, u, l, taue, taur, alphak1, alphak2, alphak3 = pickle.load(f)

# 'Rational' class ensures that the number remains in exact fractional form and is not approximated by a decimal.
cs = (gamma*boltz*T/(mu*mh))**Rational(1/2) #eq 36

Beq = bet*u*(4*pi*rho)**Rational(1/2) #equipartition, eq 2

biso = (Beq*(xio**(1/2)))/Max(1,u/cs) #xio= (b_iso/B_eq)^2, eq 1
biso = simplify(biso)
biso = biso.powsimp(force=True) #simplify expressions involving powers

bani = biso*(Rational(1/3)*2*q*omega*tau*(1+(q*omega*tau)/2))**Rational(1/2)  # eq 46 # + (Uo*tau/l)*(1+1/(1+q*omega*tau)**2)
bani = simplify(bani)
bani = bani.powsimp(force=True)

eta = (1/3)*tau*u**2 #eq 9, turbulent diffusivity of vec(B)

Ralpha = alphak*h/eta #eq 8
Romega = -q*omega*h**2/eta #eq 8
Dk = Ralpha*Romega #eq 10
Dc = -(pi**5)/32 #eq 12
Bbar = (pi*Beq*l*(Rk*(Dk/Dc-1))**(0.5))/h #eq 7, 40. #where is xi0?
Bbar = simplify(Bbar)
# Bbar = Bbar.powsimp(force=True)

tanpB = -((pi**2)*tau*(u**2))/(12*q*omega*(h**2)) #eq 41
tanpB = simplify(tanpB)
tanpB = tanpB.subs([(tau, tau), (l, l)])
tanpB = simplify(tanpB)

tanpb = 1/(1+q*omega*tau)

mag_expr = biso, bani, Bbar, tanpb, tanpB, Beq, eta, cs

with open(current_directory+ r'\expressions\mag_exp.pickle', 'wb') as f:
    pickle.dump(mag_expr, f)
print('Solved the magnetic expressions')