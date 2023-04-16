
from helper_functions import scal_finder
import numpy as np
import matplotlib.pyplot as plt
from sympy import *
import pickle
import os
from matplotlib.ticker import FormatStrFormatter
import sys
current_directory = str(os.getcwd())

# Defining the Observables
q = Symbol('q')
omega = Symbol('\Omega')
sigma = Symbol('\Sigma')
sigmatot = Symbol('Sigma_tot')
sigmasfr = Symbol('Sigma_SFR')
T = Symbol('T')

os.chdir(current_directory + '\data')
sys.path.append(current_directory + '\expressions')
sys.path.append(current_directory + '\data')



with open('zip_data.pickle', 'rb') as f:
    data_pass = pickle.load(f)

os.chdir(current_directory + '\expressions')
import magnetic_expressions as m
import turbulence_expressions as t

observable = sigma
quantity = m.biso

zr, quan_f, coeff = scal_finder(t.hsup, quantity, observable, data_pass, t.taue, t.alphak1, np.linspace(1,5000,100))
# coeff = np.mean(
#         zr[1:]*(np.gradient(np.log(np.abs(quan_f[1:])))/np.gradient(zr[1:])))
print(round(coeff, 3))
zr, quan_f2, coeff2 = scal_finder(
    t.hsup, t.u/m.cs, observable, data_pass, t.taue, t.alphak1, np.linspace(1, 5000, 100))


# plt.plot(zr, quan_f, color='r', marker='^', linestyle='-', mfc='k', mec='k', markersize=1.2,
#          label=r'$'+str(latex(observable))+r'$ exponent = ' + str(round(coeff, 3)))
plt.plot(zr, quan_f2, color='g', marker='^', linestyle='-', mfc='k', mec='k', markersize=1.2,
         label=r'$'+str(latex(observable))+r'$ exponent = ' + str(round(coeff2, 2)))

# plt.yscale('log')
# plt.xscale('log')
plt.xlabel(r'$'+str(latex(observable))+r'$')
plt.ylabel(r'$u(km/s)$', size='large')
# plt.title('Scaling relation convergence at T = '+ "{:e}".format((round(t[0],2))))
plt.legend()
# plt.savefig('t_'+str(t[0])+'.png')

plt.show()
