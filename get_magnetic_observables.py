from helper_functions import datamaker, root_finder, exp_analytical_data
import numpy as np
from sympy import *
import pickle
import os
import time
import subprocess
from multiprocessing import Pool

start_time = time.time()


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
    hg, rho, nu, u, l, taue, taur, alphak1, alphak2, alphak3, Rsb = pickle.load(f)

with open('mag_exp.pickle', 'rb') as f:
    biso, bani, Bbar, tanpb, tanpB, Beq, eta_t, cs = pickle.load(f)

os.chdir(current_directory)

def paralel_find(data_pass):
        cs_f = exp_analytical_data(cs, data_pass).astype(np.float64)
        #print(exp_analytical_data(hg, data_pass))
        try:
            h_f = root_finder(exp_analytical_data(hg, data_pass), 1e+25)
            print('Root found succesfully')
        except:
            print('*************************************************************************************')
            print('Please change the value of the initial guess')
            print('*************************************************************************************')
        #print(h_f)
        l_f = datamaker(l, data_pass, h_f)
        Rsb_f = datamaker(Rsb, data_pass, h_f)
        #print(Rsb_f/h_f)

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

        Bbar_f = datamaker(Bbar, data_pass, h_f, tau_f, alphak_f)

        tanpB_f = datamaker(tanpB, data_pass, h_f, tau_f)
        tanpb_f = datamaker(tanpb, data_pass, h_f, tau_f)
        return h_f, l_f, u_f, cs_f, alphak_f, tau_f, biso_f, bani_f, Bbar_f, tanpB_f, tanpb_f

if __name__ == "__main__":  
    def fn(data_pass):
        with Pool(5) as p:
            return p.map(paralel_find, [[data_pass[i]] for i in range(len(data_pass))])
    final_mag = (np.array(fn(data_pass))).flatten()
    final_mag = (np.reshape(final_mag, (-1,len(kpc_r)))).T
    print(final_mag)
    mag_obs = np.concatenate((np.array([kpc_r]), final_mag),  axis=0)
    mag_obs = tuple(mag_obs)
    print(mag_obs)
    os.chdir(current_directory)
    with open('mag_observables.pickle', 'wb') as f:
        pickle.dump(mag_obs, f)
    print("--- %s seconds ---" % (time.time() - start_time))
