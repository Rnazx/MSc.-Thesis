from helper_functions import datamaker, root_finder, exp_analytical_data
import numpy as np
from sympy import *
import pickle
import os

import subprocess

current_directory = str(os.getcwd())

# data extraction for galaxy
# subprocess.run(["python", "zipped_data.py"])

os.chdir(current_directory + '\data')

with open('zip_data.pickle', 'rb') as f:
    kpc_r, data_pass = pickle.load(f)

# extracting the expressions
os.chdir(current_directory + '\expressions')
# subprocess.run(["python", "turbulence_expressions.py"])
# subprocess.run(["python", "magnetic_expressions.py"])


with open('turb_exp.pickle', 'rb') as f:
    hg, rho, nu, u, l, taue, taur, alphak1, alphak2, alphak3 = pickle.load(f)

with open('mag_exp.pickle', 'rb') as f:
    biso, bani, Bbar, tanpb, tanpB, Beq, eta, cs, Dk, Dc = pickle.load(f)

os.chdir(current_directory)

cs_f = exp_analytical_data(cs, data_pass).astype(np.float64)
# h = Symbol('h')
#print(exp_analytical_data(hg, data_pass)[0].subs(h,1), cs_f)
h_init_trys = [1e+15, 1e+25, 1e+35]
for i,hi in enumerate(h_init_trys):
    try:
        print('Try {} for initial guess of h as {:e} cm'.format(i,np.round(hi)))
        h_f = root_finder(exp_analytical_data(hg, data_pass), hi)
        print('Root found succesfully')
        l_f = datamaker(l, data_pass, h_f)
        u_f = datamaker(u, data_pass, h_f)
        taue_f = datamaker(taue, data_pass, h_f)
        taur_f = datamaker(taur, data_pass, h_f)
        tau_f = np.minimum(taue_f, taur_f)

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
            alphak_f.append(datamaker(alpha_k, [data_pass[i]], np.array(
                [h_f[i]]), np.array([tau_f[i]]))[0])

        alphak_f = np.array(alphak_f)


        biso_f = datamaker(biso, data_pass, h_f, tau_f)
        bani_f = datamaker(bani, data_pass, h_f, tau_f)

        dkdc_f = datamaker((Dk/Dc), data_pass, h_f, tau_f, alphak_f)
        alpham_f = alphak_f*((1/dkdc_f)-1)

        Bbar_f = datamaker(Bbar, data_pass, h_f, tau_f, alphak_f)

        tanpB_f = datamaker(tanpB, data_pass, h_f, tau_f)
        tanpb_f = datamaker(tanpb, data_pass, h_f, tau_f)

        mag_obs = kpc_r, h_f, l_f, u_f, cs_f, alphak_f, tau_f, taue_f, taur_f, biso_f, bani_f, Bbar_f, tanpB_f, tanpb_f, dkdc_f, alpham_f, omt, kah


        os.chdir(current_directory)

        with open('mag_observables.pickle', 'wb') as f:
            pickle.dump(mag_obs, f)
        break

    except:
        continue

else: 
    print('*************************************************************************************')
    print('Please change the values of the initial guesses')
    print('*************************************************************************************')
