import numpy as np
from sympy import *
import pickle
import os

current_directory = str(os.getcwd())

params = {}
switch = {}

#change galaxy name here only
#this is imported to get_magnetic_observables, error_approximations and plots.ipynb
#nothing to be changed in the other files now except plots.ipynb
#'6946' for ngc6946, everything lowercase
gal='m33'

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

#choose galaxy name in the next line to choose the galaxy
with open(current_directory+ '\data\data_{}.pickle'.format(gal), 'rb') as f:
    data = pickle.load(f)

# #use this for m51, m31 and m33
kpc_r, dat_sigmatot, dat_sigmaHI,dat_sigmaH2, dat_q, dat_omega, dat_sigmasfr, T= data
# finding total gas density
# adds HI and H2 data based on switch 

if switch['incl_moldat'] == 'Yes': #this is set to 'Yes' in switches.in file
    print('switch for moldat',switch['incl_moldat'])
    dat_sigma = dat_sigmaHI + dat_sigmaH2
else:
    dat_sigma = dat_sigmaHI
    print('switch for moldat',switch['incl_moldat'])


#use this for 6946
#difference as here we have to mandatorily consider molecular gas density
# kpc_r, dat_sigmatot, dat_sigma, dat_q, dat_omega, dat_sigmasfr, T= data

r = kpc_r.size  # common radius of the interpolated data
dat_sigmagas = params['mu']*dat_sigma #mu= 14/11 set in parameters.in file

#####################################################################################################################

#####################################################################################################################
#to apply kennicut-schmidt relation
ks_exp = params['ks_exp']
ks_const = (dat_sigmasfr/(dat_sigma)**(ks_exp)).mean()

ks_split = switch['force_kennicut_scmidt'].split() #currently this is set to 'No, sigmadata'
if ks_split[0] == 'Yes':
    if ks_split[1] == 'sigmasfrdata':
        dat_sigma = (dat_sigmasfr/ks_const)**(1/ks_exp)
    else:
        dat_sigmasfr = ks_const*(dat_sigma)**(ks_exp)
#####################################################################################################################

#changes based on galaxy and quality of fit
zet = params['\zeta']*np.ones(r)
psi = params['\psi']*np.ones(r)
bet = params[r'\beta']*np.ones(r)

ca = params[r'C_\alpha']*np.ones(r)
rk = params['R_\kappa']*np.ones(r)
mu = params['mu']*np.ones(r)
A = params['A']*np.ones(r)

# zip function makes array of corresponding elements of each array passed into it
data_pass = kpc_r, list(zip(dat_sigmatot, dat_sigmagas, dat_sigmasfr,
                 dat_q, dat_omega, zet, T, psi, bet, ca, rk, mu, A))

#change folder name here when changing galaxies
save_files_dir=current_directory+r'\plots\M33\taue,moldat_{},z_{},psi_{},ca_{},beta_{},rk_{},A_{}'.format(switch['incl_moldat'],params[r'\zeta'],params[r'\psi'],
                                params[r'C_\alpha'],params[r'\beta'],params[r'R_\kappa'],params['A'])

#dont change this pickle file name
with open(current_directory+ '\data\zip_data.pickle', 'wb') as f:
    pickle.dump(data_pass, f)

#so that new directories wont be made when the directory name is imported from this file
if __name__ == '__main__':
    print('gal=', gal) #to ensure that we are working with correct galaxy
    try:
        os.makedirs(save_files_dir)
        os.chdir(save_files_dir)
    except FileExistsError:
        # Handle the case where the directory already exists
        print(f"The directory '{save_files_dir}' already exists, going there.")
        os.chdir(save_files_dir)
        #anything in this folder before will be re-written
    except OSError as e:
        # Handle other OSError exceptions if they occur
        print(f"An error occurred while creating the directory: {e}")


with open(save_files_dir+ '\zip_data.pickle', 'wb') as f:
    pickle.dump(data_pass, f)
print('Succesfully zipped the data and the parameters, and made pickle file')
