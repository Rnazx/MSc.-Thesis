import numpy as np
import matplotlib.pyplot as plt
from sympy import *
from fractions import Fraction
import pickle
import os
from matplotlib.ticker import FormatStrFormatter
import sys
import subprocess

current_directory = str(os.getcwd())

#conversion factors
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

#data extraction for galaxy
subprocess.run(["python", "zipped_data.py"])

os.chdir(current_directory +'\data')

with open('zip_data.pickle', 'rb') as f:
     data_pass = pickle.load(f)

#extracting the expressions
os.chdir(current_directory + '\expressions')
subprocess.run(["python", "turbulence_expressions.py"])
subprocess.run(["python", "magnetic_expressions.py"])


with open('turb_exp.pickle', 'rb') as f:
     hg, rho, nu, u, l, taue, taur, alphak1, alphak2, alphak3 = pickle.load(f)

with open('mag_exp.pickle', 'rb') as f:
     biso, bani, Bbar, tanpb, tanpB, Beq, eta = pickle.load(f)

os.chdir(current_directory)
from helper_functions import datamaker, root_finder, exp_analytical_data

h_f = root_finder(exp_analytical_data(hg, data_pass))

l_f = datamaker(l, data_pass, h_f)
u_f = datamaker(u, data_pass, h_f)
taue_f = datamaker(taue, data_pass, h_f)
taur_f = datamaker(taur, data_pass, h_f)
tau_f = taue_f

for i in range(len(tau_f)):
     if taue_f[i]<taur_f[i]: tau_f[i] = taue_f[i]
     else: tau_f[i] = taur_f[i]

omega = Symbol('\Omega')
kalpha = Symbol('K_alpha')
calpha = Symbol('C_alpha')

omt = datamaker(omega, data_pass, h_f, tau_f)*tau_f
kah = datamaker(kalpha/calpha, data_pass, h_f, tau_f)*(h_f/(tau_f*u_f))

alphak_f = []

for i in range(len(omt)):
     if min(1, kah[i]) >= omt[i]:
          alpha_k = alphak1         
     elif min(omt[i], kah[i]) >= 1:
          alpha_k = alphak2
     else:
          alpha_k = alphak3
     alphak_f.append(datamaker(alpha_k, [data_pass[i]], np.array([h_f[i]]), np.array([tau_f[i]]))[0])

alphak_f = np.array(alphak_f)




biso_f = datamaker(biso, data_pass, h_f, tau_f)
bani_f = datamaker(bani, data_pass, h_f, tau_f)

Bbar_f = datamaker(Bbar, data_pass, h_f, tau_f, alphak_f)

tanpB_f = datamaker(tanpB, data_pass, h_f, tau_f)
tanpb_f = datamaker(tanpb, data_pass, h_f, tau_f)

mag_obs = h_f, l_f, u_f, alphak_f, tau_f, biso_f, bani_f, Bbar_f, tanpB_f, tanpb_f

os.chdir(current_directory)

with open('mag_observables.pickle', 'wb') as f:
    pickle.dump(mag_obs, f)