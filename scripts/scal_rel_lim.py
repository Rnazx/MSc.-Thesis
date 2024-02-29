# import libraries
import numpy as np
import matplotlib.pyplot as plt
from sympy import *
from fractions import Fraction
import os


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
mu0 = Symbol('Mu_0')
mh = Symbol('m_H')
G = Symbol('G')
xio = Symbol('xi_0')
delta = Symbol('\delta')
mstar = Symbol('m_*')
cl = Symbol('C_l')
kappa = Symbol('kappa')
E51 = Symbol('E_51')
Rk = Symbol('R_k')
kalpha = Symbol('K_alpha')
Gamma = Symbol('Gamma')
eta = Symbol('eta')
Nsb = Symbol('N_sb')


# Defining the general parameters
u = Symbol('u')
tau = Symbol('tau')
l = Symbol('l')
max_mach = Symbol('M')
n = Symbol('n')
nu = Symbol('nu')

# conversion factors
pc_kpc = 1e3  # number of pc in one kpc
cm_kpc = 3.086e+21  # number of centimeters in one parsec
cm_km = 1e5  # number of cm in one km
s_Myr = 1e+6*(365*24*60*60)  # number of seconds in one megayear

base_path = os.environ.get('MY_PATH')

#Model generator based on the regime

def model_gen_regime(hreg, lreg = 'supernovae-driven', ureg = 'supernovae-driven', taureg = 'eddy turnover time', alphareg='regime 1'):
    cs = (gamma*boltz*T/(mu*mh))**Rational(1/2)

    # Function to return the expression of l given the regime
    def choose_lreg(h, lreg):
        nu = (delta*sigmasfr)/(2*h*mstar)
        # dependence of gas density on scale height(h)
        rho = sigma/(2*h)
        # number density
        n = rho/(mu*mh)
        # If we include correlation length from the supernovae model from Chamandy and Shukurov (2020)
        #Model 3
        #fsb = 0
        if lreg == 'supernovae-driven':
            ls = 0.14*cm_kpc*(E51)**Fraction(16, 51) * \
                (n/0.1)**Fraction(-19, 51)*(cs/(cm_km*10))**Fraction(-1, 3)
            nu = nu
            l = ((Gamma-1)/Gamma)*cl*ls
            l = simplify(l)
        #fsb = 1, Nsb is set to 1 later on as it does not affect the scaling relations
        elif lreg == 'superbubble-driven':
            #Eqn 10 Chamandy and Sukurov (2020)
            Rsb = 0.53*cm_kpc*(eta/0.1)**Fraction(1, 3)*(Nsb/100)**Fraction(1, 3)*(E51)**Fraction(1, 3)*(n/0.1)**Fraction(-1, 3)*(cs/(cm_km*10))**Fraction(-2, 3)
            ls = Rsb
            nu = nu/Nsb
            l = ((Gamma-1)/Gamma)*cl*ls
            l = simplify(l)
        # Minimalistic model for l in model 1 and 2. Turbulence is driven at the maximum scale (h)
        else:
            l = simplify(cl*h)
            ls = simplify(cl*h)
        return l, ls, n, nu
    # Function to return the expression of u given the regime.
    def choose_ureg(h, ureg):
        # We first choose l according to the regime
        l, ls, n, nu = choose_lreg(h, lreg)
        # If u is set equal to the sound speed cs
        #Model 1
        if ureg == 'sound speed' :
            u = mu0*cs
        #If the regime is supernovae driven, the expression for u is taken from Chamandy and Shukurov (2020)
        #Models 2 and 3
        else:
            u = ((4*pi/3)*l*ls**3*cs**2*(nu))**Fraction(1, 3)            
        return l,ls,u,n,nu
    # Edge case where in the expression for h, we assume u<<cs (w = u^2 +cs^2)
    if hreg == 'subsonic':
        h = (cs**2)/(3*pi*G*sigmatot)
    # Edge case where in the expression for h, we assume u>>cs
    # This case is slightly more complicate since u indirectly depends on h through l
    # Hence we need to algebraically solve the equation 
    elif hreg == 'supersonic':
        #First define h as a symbol
        h = Symbol('h')
        #Solve for l and u as a function of h which is just a symbol
        l,ls,u,n,nu = choose_ureg(h, ureg)
        # define a new quantity us where h is set to 1 in the expression for u
        us = u.subs(h, 1)
        # The standard expression for h raised to the power of 1/(1-2*(exponent of h in u))
        h = simplify(((us**2)/(3*pi*G*sigmatot))
                            ** (1/(1-2*diff(log(u), h)*h)))#exponent is found by differentiating the log
    #alternate expression for h
    else :
        h = cs/omega
    l,ls,u,n,nu = choose_ureg(h, ureg)
    # Define max_mach as max(1,u/cs) which will be used later on in the expression for b_iso
    if hreg == 'supersonic':
        max_mach = u/cs
    else:
        max_mach = 1
    # Regime for the turbulence correlation time
    # If the regime for tau is eddy turnover time
    taue = simplify(l/u)
    taur = simplify(((4/3)*pi*nu*(ls**3))**(-1))
    if taureg == 'eddy turnover time':
        tau = taue
    # If the regime for tau is supernovae renovation time
    # The quantities are all converted to cgs units for standardization
    else:
        tau = taur # does not converge if model no = 2
    # after solving for h we put it into the other expressions
    #define rho again for the supersonic case as we have to input the final expression for h
    rho = sigma/(2*h)
    
    # Magnetic field considering equipartition
    Beq = u*(4*pi*rho)**Rational(1/2)
    
    #Expression for the isotropic magnetic field from Federrath et. al. (2011)
    #max_mach is max(1, u/cs)
    biso = (Beq*(xio**(1/2)))/max_mach# change the name
    biso = simplify(biso)
    biso = biso.powsimp(force=True)

    #Expression for the anisotropic magnetic field considering differential rotation
    bani = biso*(Rational(2/3)*q*omega*tau)**(1/2)# mention the approximations
    bani = simplify(bani)
    bani = bani.powsimp(force=True)
    # Expression for the mean magnetic field from Dynamo theory
    Rk = Symbol('R_k')
    eta_t = (1/3)*tau*u**2
    #Three regimes for alpha_k are chosen and scaling relations can be found for each regime
    alphak1 = calpha*tau**2*u**2*omega/h
    alphak2 = calpha*tau*u**2/h
    alphak3 = kalpha*u
    if alphareg == 'regime 1':
        alphak = alphak1
    elif alphareg == 'regime 2':
        alphak = alphak2
    else:
        alphak = alphak3
    # Substitute alpha_k in the reynolds numbers
    Ralpha = alphak*h/eta_t
    Romega = -q*omega*h**2/eta_t
    #Dynamo numbers
    Dk = Ralpha*Romega
    Dc = -(pi**5)/32
    # Final expression for mean magnetic field after some approximations
    Bbar = simplify((pi*Beq*l*(Rk*(Dk/Dc))**(1/2))/h)

    # Expression for the pitch angle of the mean magnetic field
    tanpB = -((pi**2)*tau*(u**2))/(12*q*omega*(h**2))
    tanpB = simplify(tanpB)
    # Substitute tau and l in tanpB
    tanpB = tanpB.subs([(tau, tau), (l, l)])# change the names
    tanpB = simplify(tanpB)

    # Expression for the pitch angle of the random magnetic field
    tanpb = 1/(1+q*omega*tau)

    #Put all the expression in a single list
    quantities = [h, l, u, cs, alphak, taue, taur, biso, bani, Bbar, tanpB, tanpb, Dk/Dc]
#hg, l, u, cs, alphak1, taue, taur, biso, bani, Bbar, tanpB, tanpb, Dk/Dc
    return quantities
g_Msun = 1.989e33  # solar mass in g
cgs_G = 6.674e-8  # gravitational constant in cgs units
g_mH = 1.6736e-24  # mass of hydrogen atom in grams
cgs_kB = 1.3807e-16  # boltzmann constant in cgs units

# Reading the Constant values
gval, clval, xioval, mstarval, deltaval, e51val, kaval, Gammaval, Caval, Rkval, muval, mu0val, g_Msun, cgs_G, g_mH, cgs_kB = tuple(
    np.genfromtxt(os.path.join(base_path,'scripts','constants.in'), delimiter='=', dtype=np.float64)[:, -1])

# List of tuples for substituting the values in the symbol. 
# The firt element of each tuple is the symbol for which the value needs to be substituted
# The second element is the numerical value which is stored in constants.in file
const = [(boltz, cgs_kB), (mh, g_mH), (G, cgs_G), (gamma, gval),
         (calpha, Caval), (Rk, Rkval), (mu, muval), (cl,
                                               clval), (xio, xioval), (mstar, mstarval*g_Msun),
         (delta, deltaval), (E51, e51val), (kalpha, kaval), (Gamma, Gammaval), (mu0, mu0val)]
def scal_finder(quantities):
    #q ,omega ,sigma,sigmatot,sigmasfr,T
    variables = [(q, 1),(omega, 1),(sigma, 1),(sigmatot, 1),  (sigmasfr, 1),
                  (T, 1)]
    # observable to be varied
    observ = [variables[i][0] for i in range(len(variables))]
    # plotting the scaling relations
    exps=[]
    for quan in quantities:
        power = []
        for obs in observ:
            variables = [(sigmatot, 1), (sigma, 1), (sigmasfr, 1),
                        (omega, 1), (q, 1), (T, 1)]
            variables.remove((obs, 1))
            final = const + variables
            z = quan.subs(final)
            po = np.float64((diff(log(z), obs)*obs).subs(obs, 1))
            power.append(po)
        exps.append(np.array(power))
    return np.array(exps)

hregs = ['subsonic', 'supersonic']
for hreg in hregs:
    quantities = model_gen_regime(hreg)
    exps = scal_finder(quantities)
    os.chdir(os.path.join(base_path,'inputs'))
    np.save('scal_exponents_'+hreg,np.array(exps))
