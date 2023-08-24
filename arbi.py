import numpy as np
from sympy import *
import pickle
import os
import subprocess
import sys

current_directory = str(os.getcwd())  # store the path of the current directory

params = {}
switch = {}

with open('parameter_file.in', 'r') as FH:
    for file in FH.readlines():
        line = file.strip()
        try:
            par_name, value = line.split('=')
        except ValueError:
            print("Record: ", line)
            raise Exception(
                "Failed while unpacking. Not enough arguments to supply.")
        try:
            params[par_name] = np.float64(value)
        except ValueError:
            num, denom = value.split('/')
            params[par_name] = np.float64(num) / np.float64(denom)

with open('switches.in', 'r') as FH:
    for file in FH.readlines():
        line = file.strip()
        try:
            sw_name, value = line.split('=')
        except ValueError:
            print("Record: ", line)
            raise Exception(
                "Failed while unpacking. Not enough arguments to supply.")
        switch[sw_name] = value

print('Succesfully read the parameters and switches')


# data extraction for galaxy
os.chdir(current_directory + '\data')
if switch['chem_or_claude'] is 'Chemin':
    c_or_cl = 1
else:
    c_or_cl = 0
#args = sys.argv[1:]
args = np.float64(np.array([float(s) for s in sys.argv[1:]]))
print('The value for the Observables are sigmatot, sigmaHI, q, omega, sigmasfr\n', args)
obs_array = [np.array([arg]) for arg in args]
dat_sigmatot, dat_sigmaHI, dat_q, dat_omega, dat_sigmasfr, molfrac = tuple(obs_array)
kpc_r = np.ones(1)
r = kpc_r.size  # common radius of the interpolated data

dat_sigmaH2 = dat_sigmaHI*(1/(1-molfrac))

T_tb = (0.017*kpc_r + 0.5)*1e+4

if switch['incl_moldat'] is 'Yes':
    dat_sigma = dat_sigmaHI + dat_sigmaH2
else:
    dat_sigma = dat_sigmaHI

dat_sigma = params['mu']*dat_sigma

ks_exp = params['ks_exp']
ks_const = (dat_sigmasfr/(dat_sigma)**(ks_exp)).mean()

ks_split = switch['force_kennicut_scmidt'].split()
if ks_split[0] is 'Yes':
    if ks_split[1] is 'sigmasfrdata':
        dat_sigma = (dat_sigmasfr/ks_const)**(1/ks_exp)
    else:
        dat_sigmasfr = ks_const*(dat_sigma)**(ks_exp)

os.chdir(current_directory)


zet = params['\zeta']*np.ones(r)
psi = params['\psi']*np.ones(r)
bet = params[r'\beta']*np.ones(r)

T = T_tb  # 1e+4*np.ones(r)
ca = params[r'C_\alpha']*np.ones(r)
rk = params['R_\kappa']*np.ones(r)
mu = params['mu']*np.ones(r)
data_pass = list(zip(dat_sigmatot, dat_sigma, dat_sigmasfr,
                 dat_q, dat_omega, zet, T, psi, bet, ca, rk, mu))


os.chdir(current_directory + '\data')
print('Succesfully zipped the data and the parameters')

with open('arb_data.pickle', 'wb') as f:
    pickle.dump(data_pass, f)
