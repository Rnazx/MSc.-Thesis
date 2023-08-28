import numpy as np
import matplotlib.pyplot as plt
from sympy import *
from fractions import Fraction
import pickle
import sys


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
Nsb = Symbol('N_SB')#number of SNR in a SB
eta = Symbol('eta')#fractionof SB that is mechanical
xi = Symbol('xi')#SB horizontal radius at blowout as fraction of H
fsb = Symbol('f_SB')
Uo = Symbol('U_0') #outflow speed

# Defining the general parameters
u = Symbol('u')
tau = Symbol('tau')
l = Symbol('l')
h = Symbol('h')

##############################################################################################################
with open('turb_exp.pickle', 'rb') as f:
    hg, rho, nu, u, l, lsn, lsb,_Esn_Esb, taue, taur, alphak1, alphak2, alphak3, Rsb = pickle.load(f)
cs = (gamma*boltz*T/(mu*mh))**Rational(1/2)


Beq = bet*u*(4*pi*rho)**Rational(1/2)
biso = (Beq*(xio**(1/2)))/Max(1,u/cs)
biso = simplify(biso)
biso = biso.powsimp(force=True)


bani = biso*(Rational(1/3)*2*q*omega*tau*(1+(q*omega*tau)/2)
             )**Rational(1/2)   #+ (Uo*tau/l)*(1+(1/(1+q*omega*tau)**2))
bani = simplify(bani)
bani = bani.powsimp(force=True)

Rk = Symbol('R_k')
eta_t = (1/3)*tau*u**2

Ralpha = alphak*h/eta_t
Romega = -q*omega*h**2/eta_t
Ru = Uo*h/eta_t
Dk = Ralpha*Romega
Dc = -((pi**5)/32)*(1+Ru/(pi**2))
Bbar = (Beq*l*((Ru+pi**2*Rk)*(Dk/Dc-1))**(0.5))/h
Bbar = simplify(Bbar)
# Bbar = Bbar.powsimp(force=True)

tanpB = (pi**2+Ru)/(4*Romega)
tanpB = simplify(tanpB)
tanpB = tanpB.subs([(tau, tau), (l, l)])
tanpB = simplify(tanpB)

tanpb = 1/(1+q*omega*tau)

mag_expr = biso, bani, Bbar, tanpb, tanpB, Beq, eta_t, cs


with open('mag_exp.pickle', 'wb') as f:
    pickle.dump(mag_expr, f)
print('Solved the magnetic expressions')