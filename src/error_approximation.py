
####################################################################################################################

import numpy as np
import pickle
import os
from helper_functions import parameter_read

base_path = os.environ.get('MY_PATH')
galaxy_name = os.environ.get('galaxy_name')
params = parameter_read(os.path.join(base_path,'inputs','parameter_file.in'))

current_directory = str(os.getcwd())
#kpc_r, h_f, l_f, u_f, cs_f, alphak_f, tau_f, taue_f, taur_f, biso_f, bani_f, Bbar_f, tanpB_f, tanpb_f, dkdc_f
os.chdir(os.path.join(base_path,'outputs'))

with open(f'{galaxy_name}output_ca_'+str(params[r'C_\alpha'])+'K_'+str(params[r'K'])+'z_'+
          str(params[r'\zeta'])+'psi_'+str(params[r'\psi'])+'b_'+str(params[r'\beta'])+'.out', 'rb') as f:
                model_f = pickle.load(f)

os.chdir(os.path.join(base_path,'inputs'))

with open('zip_data.in', 'rb') as f:
    kpc_r, data_pass = pickle.load(f)
hregs = ['subsonic', 'supersonic']
for hreg in hregs:
    os.chdir(os.path.join(base_path,'inputs'))
    exps = np.load(f'scal_exponents_{hreg}.npy')
    r = kpc_r.size

    # dat_sigmatot, dat_sigma, dat_sigmasfr, dat_q, dat_omega, zet, T_tb, psi, bet, ca, rk, mu = (
    #     np.array([data_pass[i][j] for i in range(r)]) for j in range(len(data_pass[0])))

    # galaxy_name = os.environ.get('galaxy_name')
    # os.chdir(os.path.join(base_path,'data', 'supplementary_data'))
    # os.chdir(galaxy_name)

    # from data_conv_M31 import kpc_r_Chemin, kms_vcirc_error_Chemin, cm_km, cm_kpc, pc_kpc, kms_vcirc_Chemin
    # #make a seperate file for all the errors
    # kmskpc_Om = kms_vcirc_Chemin/kpc_r_Chemin
    # err_kmskpc_Om = kms_vcirc_error_Chemin/kpc_r_Chemin
    # err_omega = err_kmskpc_Om*cm_km/cm_kpc
    # err_q = (kpc_r_Chemin/(kmskpc_Om*np.gradient(kpc_r_Chemin)))*((np.gradient(kmskpc_Om)/kmskpc_Om)-1)*err_kmskpc_Om
    # err_q = griddata(kpc_r_Chemin, err_q, kpc_r, method='linear',
    #                 fill_value=nan, rescale=False)
    # err_omega = griddata(kpc_r_Chemin, err_omega, kpc_r,
    #                     method='linear', fill_value=nan, rescale=False)

    percent = 0.1
    relerr_q = percent*np.ones(r)#err_q/np.abs(dat_q)
    relerr_omega = percent*np.ones(r)#err_omega/np.abs(dat_omega)
    relerr_sigma = percent*np.ones(r)
    relerr_sigmatot = percent*np.ones(r)
    relerr_sigmasfr = percent*np.ones(r)

    err_T = (0.005*kpc_r + 0.1)*1e+4 #from Tabatabaei+13b equation ??
    relerr_T = percent*np.ones(r) #err_T/np.abs(T_tb)
    rel_err = np.array([relerr_q, relerr_omega, relerr_sigma, relerr_sigmatot,relerr_sigmasfr, relerr_T])

    relerr_quan = np.sqrt(np.matmul(exps**2,rel_err**2))
    err_quantities = model_f[1:]*relerr_quan

    # os.chdir(current_directory)
    #print(err_quantities)
    os.chdir(os.path.join(base_path,'outputs'))
    with open(f'errors_{hreg}.out', 'wb') as f:
        pickle.dump(err_quantities, f)
print('Found the errors from the scaling relations')