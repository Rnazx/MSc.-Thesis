import numpy as np
from sympy import *
import pickle
import os
import subprocess

current_directory=r'D:\Documents\Gayathri_college\MSc project\codes\MSc.-Thesis' # store the path of the current directory

params = {}
switch = {}

#opening these files and making them into dictionaries
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

with open(current_directory+'\switches.in', 'r') as FH:
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

#choose "data_conv_M33.py" or "data_conv_M31.py" in the next two lines to choose the galaxy
subprocess.run(["python", current_directory+"\data\data_m51.py"]) #, str(int(c_or_cl)
with open('data_M51.pickle', 'rb') as f:
    data = pickle.load(f)

kpc_r,dat_sigmatot,dat_sigmaHI,dat_sigmaH2,dat_sigmagas, dat_q, dat_omega, dat_sigmasfr,dat_sigmasfr_fuv,temp,vel_disp = data

r = kpc_r.size  # common radius of the interpolated data

#this is if total gas density is isnt available
#not used for M51 data
#adds HI and H2 data based on switch 
if switch['incl_moldat'] == 'Yes': #this is set to 'YES' in switches.in file
    dat_sigma = dat_sigmaHI + dat_sigmaH2
else:
    dat_sigma = dat_sigmaHI

#finding total gas density
dat_sigmagas = params['mu']*dat_sigma #mu= 14/11 set in parameters.in file

#to apply kennicut-schmidt relation
ks_exp = params['ks_exp']
ks_const = (dat_sigmasfr/(dat_sigma)**(ks_exp)).mean()

ks_split = switch['force_kennicut_scmidt'].split() #currently this is set to 'No, sigmadata'
if ks_split[0] == 'Yes':
    if ks_split[1] == 'sigmasfrdata':
        dat_sigma = (dat_sigmasfr/ks_const)**(1/ks_exp)
    else:
        dat_sigmasfr = ks_const*(dat_sigma)**(ks_exp)

#to switch between fuv and h_alpha data
#sfr_list={'h_alpha':0,'FUV':1}
sfr_val=0
if sfr_val==0: #choosing H_alpha data
    del list(data)[8]
else: #chose FUV data
    del list(data)[7]


#changes based on galaxy and quality of fit
zet = params['\zeta']*np.ones(r)
psi = params['\psi']*np.ones(r)
bet = params[r'\beta']*np.ones(r)

T = temp  #obtained from paper
ca = params[r'C_\alpha']*np.ones(r)
rk = params['R_\kappa']*np.ones(r)
mu = params['mu']*np.ones(r)

# zip function makes array of corresponding elements of each array passed into it
data_pass = kpc_r, list(zip(dat_sigmatot, dat_sigma, dat_sigmasfr,
                 dat_q, dat_omega, zet, T, psi, bet, ca, rk, mu))


with open(current_directory+'\data\zip_data.pickle', 'wb') as f:
    pickle.dump(data_pass, f)

with open('zip_data_M51.pickle', 'wb') as f:
    pickle.dump(data_pass, f)

print('Succesfully zipped the data and the parameters, and made pickle file')

