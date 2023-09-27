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
#kpc_r, h_f, l_f, u_f, cs_f, alphak_f, tau_f, taue_f, taur_f, biso_f, bani_f, Bbar_f, tanpB_f, tanpb_f, dkdc_f, alpham_f, omt, kah
with open('mag_observables.pickle', 'rb') as f:
    model_f = pickle.load(
        f)

os.chdir(current_directory + '\data')
with open('zip_data.pickle', 'rb') as f:
    kpc_r, data_pass = pickle.load(f)

r = kpc_r.size

dat_sigmatot, dat_sigma, dat_sigmasfr, dat_q, dat_omega, zet, T_tb, psi, bet, ca, rk, mu = (
    np.array([data_pass[i][j] for i in range(r)]) for j in range(len(data_pass[0])))


from data.data_conv_M31 import kpc_r_Chemin, kms_vcirc_error_Chemin, cm_km, cm_kpc, pc_kpc, kms_vcirc_Chemin

kmskpc_Om = kms_vcirc_Chemin/kpc_r_Chemin
err_kmskpc_Om = kms_vcirc_error_Chemin/kpc_r_Chemin
err_omega = err_kmskpc_Om*cm_km/cm_kpc
err_q = (kpc_r_Chemin/(kmskpc_Om*np.gradient(kpc_r_Chemin)))*((np.gradient(kmskpc_Om)/kmskpc_Om)-1)*err_kmskpc_Om
err_q = griddata(kpc_r_Chemin, err_q, kpc_r, method='linear',
                fill_value=nan, rescale=False)
err_omega = griddata(kpc_r_Chemin, err_omega, kpc_r,
                    method='linear', fill_value=nan, rescale=False)


relerr_q = err_q/np.abs(dat_q)
relerr_omega = err_omega/np.abs(dat_omega)
relerr_sigma = 0.01*np.ones(r)
relerr_sigmatot = 0.01*np.ones(r)
relerr_sigmasfr = 0.01*np.ones(r)

err_T = (0.005*kpc_r + 0.1)*1e+4 #from Tabatabaei+13b equation ??
relerr_T = err_T/np.abs(T_tb)
rel_err = np.array([relerr_q, relerr_omega, relerr_sigma, relerr_sigmatot,relerr_sigmasfr, relerr_T])
exps = np.load('scal_exponents.npy')
relerr_quan = np.sqrt(np.matmul(exps**2,rel_err**2))
err_quantities = model_f[1:]*relerr_quan

os.chdir(current_directory)
#print(err_quantities)
with open('errors_quan.pickle', 'wb') as f:
    pickle.dump(err_quantities, f)
print('Found the errors from the scaling relations')
# plt.errorbar(kpc_r, u_f/cm_km, err_h/cm_km, c='r', linestyle='-',ms=2, mew=2, capsize=2,
#                   ecolor = 'y', marker='o', mfc='k', mec='k', label=r' $h$(pc)')

# plt.xlabel(r'Radius (kpc)', fontsize=15)
# plt.ylabel(r'Length scale (pc)', fontsize=15)
# plt.legend(fontsize = 10)
# plt.show()