
from helper_functions import scal_finder
import numpy as np
import matplotlib.pyplot as plt
from sympy import *
import pickle
import os
import sys
current_directory = str(os.getcwd())

# Defining the Observables
q = Symbol('q')
omega = Symbol('\Omega')
sigma = Symbol('\Sigma')
sigmatot = Symbol('Sigma_tot')
sigmasfr = Symbol('Sigma_SFR')
T = Symbol('T')
zet = Symbol('zeta')

# os.chdir(current_directory + '\data')
# sys.path.append(current_directory + '\expressions')
# sys.path.append(current_directory + '\data')

os.chdir(current_directory + '\expressions')
with open('turb_exp.pickle', 'rb') as f:
    hg, rho, nu, u, l, taue, taur, alphak1, alphak2, alphak3, Rsb = pickle.load(f)

from expressions.turbulence_expressions import hsup, hsub

with open('mag_exp.pickle', 'rb') as f:
    biso, bani, Bbar, tanpb, tanpB, Beq, eta_t, cs = pickle.load(f)

os.chdir(current_directory + '\data')
with open('zip_data.pickle', 'rb') as f:
    kpc_r, data_pass = pickle.load(f)
r = kpc_r.size

# with open('zip_data.pickle', 'rb') as f:
#     kpc_r, data_pass = pickle.load(f)

# os.chdir(current_directory + '\expressions')
# import magnetic_expressions as m
# import turbulence_expressions as t

observable = sigma
quantity = Rsb/hsup

def scal_plotter(h_exp, quantity, observable, h_regime):
    zr, quan_f, coeff= scal_finder(h_exp, quantity, observable, data_pass, taue, alphak1, np.linspace(1,5000,25), 1e+25)
    # coeff = np.mean(
    #         zr[1:]*(np.gradient(np.log(np.abs(quan_f[1:])))/np.gradient(zr[1:])))
    #print(coeff)
    plt.plot(zr, quan_f, marker='^', linestyle='-', mfc='k', mec='k', markersize=1.2,
            label=str(h_regime)+r'$'+str(latex(observable))+r'$ exponent = ' + str(round(coeff, 3)))

scal_plotter(hg, quantity, observable, 'general')
scal_plotter(hsub, quantity, observable, 'subsonic')
scal_plotter(hsup, quantity, observable, 'supersonic')


plt.xlabel(r'$'+str(latex(observable))+r'$')
plt.ylabel(r'$u(km/s)$', size='large')
plt.title('Scaling relation convergence') #at T = '+ "{:e}".format((round(t[0],2))))
# plt.yscale('log')
# plt.xscale('log')
plt.grid(True, which="both", ls="--")
plt.legend()
# plt.savefig('t_'+str(t[0])+'.png')

plt.show()
