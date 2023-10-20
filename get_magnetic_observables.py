from helper_functions import datamaker, root_finder, exp_analytical_data
import numpy as np
from sympy import *
import pickle
import os

from zipped_data_edited import save_files_dir, gal

########################################################################################################
current_directory = str(os.getcwd())
os.chdir(current_directory + '\data')

with open('zip_data.pickle', 'rb') as f: # there is only 1 zip_data file. dont change name of pickle file to be used here
    kpc_r, data_pass = pickle.load(f) 
print(data_pass[-1][:])
########################################################################################################

# extracting the expressions
os.chdir(current_directory + '\expressions')

with open('turb_exp.pickle', 'rb') as f:
    hg, rho, nu, u, l, taue, taur, alphak1, alphak2, alphak3 = pickle.load(f)

with open('mag_exp.pickle', 'rb') as f:
    biso, bani, Bbar, tanpb, tanpB, Beq, eta, cs = pickle.load(f)
########################################################################################################

os.chdir(current_directory)

cs_f = exp_analytical_data(cs, data_pass).astype(np.float64)

try:
    h_f = root_finder(exp_analytical_data(hg, data_pass), 1e+27) #1e+15 is initial guess, can be changed if its not converging
    print('Root found succesfully')
except:
    print('*************************************************************************************')
    print('Please change the value of the initial guess')
    print('*************************************************************************************')
l_f = datamaker(l, data_pass, h_f) #subscript f stand for final
u_f = datamaker(u, data_pass, h_f)
taue_f = datamaker(taue, data_pass, h_f)
taur_f = datamaker(taur, data_pass, h_f)

#17/10/2023
# include tau_r
# tau_f = np.minimum(taue_f, taur_f) #to see if results vary significantly for M51 and 6946
tau_f=taue_f #ignoring tau_r in tau calculation 

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

mag_obs = kpc_r, h_f, l_f, u_f, cs_f, alphak_f, tau_f, biso_f, bani_f, Bbar_f, tanpB_f, tanpb_f

with open('mag_observables_{}.pickle'.format(gal), 'wb') as f:
    pickle.dump(mag_obs, f)

########################################################################################################
# params to go to folder 

#opening these files and making them into dictionaries
params = {}

with open(current_directory+'\parameter_file.in', 'r') as FH:
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
        except ValueError: #required cz of 14/11 in parameter.in file
            num, denom = value.split('/')
            params[par_name] = np.float64(num) / np.float64(denom)

os.chdir(save_files_dir)
with open(save_files_dir+r'\mag_observables_{}.pickle'.format(gal), 'wb') as f:
    pickle.dump(mag_obs, f)
########################################################################################################
