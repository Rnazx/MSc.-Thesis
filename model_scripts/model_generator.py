import numpy as np
import matplotlib.pyplot as plt
from sympy import *
from fractions import Fraction

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
zet = Symbol('zeta')
E51 = Symbol('E_51')

# Defining the general parameters
u = Symbol('u')
tau = Symbol('tau')
l = Symbol('l')
mach = Symbol('M')
n = Symbol('n')

#conversion factors
cm = 1
kpc = 3.086e+21*cm
kms = 1e+5*cm
Myr = 1e+6*(365*24*60*60)
# Defining the expressions
cs = (gamma*boltz*T/(mu*mh))**Rational(1/2)
# Ralpha = alphak*h/eta
# Romega = -q*omega*h**2/eta

###############################################################################################
###########################Function to generate the model###########################################
#######################################################################################
def model_gen(model_no, let, tbool = True):
    def find_hlut(model_no, let):
        def choose_lreg(h, model_no):
            rho = sigma/(2*h)
            n = rho/((14/11)*mh)
            if [3, 4].count(model_no)>0 :
                lsn = 0.14*kpc*(E51)**Fraction(16, 51) * (n/0.1)**Fraction(-19, 51)*(cs/(kms*10))**Fraction(-1, 3)
                l = (3/10)*lsn
                l = simplify(l)
            else:
                l = simplify(cl*h)
                lsn = simplify(cl*h)
            return l, lsn, n
        if let == 'a':
            if model_no == 4:
                h = zet*cs/omega
            else: 
                h = zet*(cs**2)/(3*pi*G*sigma)
            nu = (delta*sigmasfr)/(2*h*mstar)
            l, lsn, n = choose_lreg(h, model_no)
            usn = ((4*pi/3)*l*lsn**3*cs**2*nu)**Fraction(1, 3)
            u = simplify(usn)
            h = simplify(h)
            
        else:
            h = Symbol('h')
            nu = (delta*sigmasfr)/(2*h*mstar)
            l, lsn, n = choose_lreg(h, model_no)
            u = ((4*pi/3)*l*lsn**3*cs**2*(nu))**Fraction(1, 3)
            usn = u.subs(h, 1)
            if model_no == 4:
                h = simplify((zet*(usn)/(omega))**Fraction(153, 128))
            elif model_no == 3:
                h = simplify((zet*(usn**2)/(3*pi*G*sigma))**Fraction(153, 103))
            else :
                h = simplify((zet*(usn**2)/(3*pi*G*sigmatot))**Fraction(-1, 1))
            l, lsn, n = choose_lreg(h, model_no)
            nu = (delta*sigmasfr)/(2*h*mstar)
            u = ((4*pi/3)*l*lsn**3*cs**2*(nu))**Fraction(1, 3)

        if tbool:
            taue = simplify(l/u)
            tau = taue
        else :
            taur = simplify(6.8*Myr*(1/4)*(nu*kpc**3*Myr/50)**(-1)*(E51)**Fraction(-16, 51) * (n/0.1)**Fraction(19, 17)*(cs/(kms*10)))
            tau = taur
        return [h, l ,u, tau]
    h, l, u, tau = find_hlut(model_no, let)
    rho = sigma/(2*h)
    Beq = u*(4*pi*rho)**Rational(1/2)
    biso = (Beq*(xio**(1/2)))
    biso = simplify(biso)
    biso = biso.powsimp(force=True)

    bani = biso*(Rational(2/3)*q*omega)**Rational(1/2)*(tau**Rational(1/2))
    bani = simplify(bani)
    bani = bani.powsimp(force=True)

    Rk = Symbol('R_k')
    Dk = -(9*calpha*q*(h**2)*(omega**2))/u**2
    Dc = -(pi**5)/32
    rho = sigma/(2*h)
    Beq = (4*pi*rho)**Rational(1/2)*u
    Bbar = (pi*Beq*l*(Rk*(Dk/Dc))**Rational(1/2))/h
    Bbar = simplify(Bbar)
    Bbar = Bbar.powsimp(force=True)




    tanpb = -((pi**2)*tau*(u**2))/(12*q*omega*(h**2))
    tanpb = simplify(tanpb)
    tanpb = tanpb.subs([(tau, tau), (l, l)])
    tanpb = simplify(tanpb)

    quantities = [ h, l, u, tau, biso, bani, Bbar, tanpb ]

    return quantities

#########function for bining the data#####################################
def bin_data(x, y, start, stop, step ):
     bdata = []
     for i in range(start, stop, step):
         idx = np.where((x>i)*(x<i+2))
         bdata.append(np.take(y,idx).mean())
     return np.array(bdata)
###############################################################################################
###########################Function to derive the model in terms of observables###########################################
#######################################################################################