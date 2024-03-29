import numpy as np
import matplotlib.pyplot as plt
from sympy import *
from fractions import Fraction
import pickle
import os
from matplotlib.ticker import FormatStrFormatter
import sys
from scipy.interpolate import griddata
from helper_functions import scal_finder
current_directory = str(os.getcwd())

q = Symbol('q')
omega = Symbol('\Omega')
sigma = Symbol('\Sigma')
sigmatot = Symbol('Sigma_tot')
sigmasfr = Symbol('Sigma_SFR')
T = Symbol('T')

os.chdir(current_directory + '\expressions')
with open('turb_exp.pickle', 'rb') as f:
    hg, rho, nu, u, l, taue, taur, alphak1, alphak2, alphak3 = pickle.load(f)

with open('mag_exp.pickle', 'rb') as f:
    biso, bani, Bbar, tanpb, tanpB, Beq, eta, cs = pickle.load(f)

os.chdir(current_directory + '\data')
with open('zip_data.pickle', 'rb') as f:
    kpc_r, data_pass = pickle.load(f)
r = kpc_r.size

observables = [q ,omega ,sigma,sigmatot,sigmasfr,T]
quantities = [hg, l, u, cs, alphak1, taue, biso, bani, Bbar, tanpB, tanpb]
#kpc_r, h_f, l_f, u_f, cs_f, alphak_f, tau_f, taue_f, taur_f, biso_f, bani_f, Bbar_f, tanpB_f, tanpb_f, dkdc_f, alpham_f, omt, kah
exps = []
for i,quan in enumerate(quantities):
    exps_quan = np.array([scal_finder(hg, quan, obs, data_pass, taue, alphak1, np.linspace(1,5000,100), 1e+15)[2] for obs in observables])
    exps.append(exps_quan)

np.save('scal_exponents_model1',np.array(exps))
print('The scaling relations are calculated')
