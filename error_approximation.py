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
#kpc_r, h_f, l_f, u_f, cs_f, alphak_f, tau_f, biso_f, bani_f, Bbar_f, tanpB_f, tanpb_f
with open(current_directory+'\mag_observables_m31.pickle', 'rb') as f: #change file name here
    model_f = pickle.load(f)

os.chdir(current_directory + '\data')
with open('zip_data.pickle', 'rb') as f:
    kpc_r, data_pass = pickle.load(f)

r = kpc_r.size

dat_sigmatot, dat_sigma, dat_sigmasfr, dat_q, dat_omega, zet, T_tb, psi, bet, ca, rk, mu = (
    np.array([data_pass[i][j] for i in range(r)]) for j in range(len(data_pass[0])))

os.chdir(current_directory)
#################################################################################################################################
#galaxy specific data 

# for m31
from data.data_conv_M31 import kpc_r_Chemin, kms_vcirc_error_Chemin, cm_km, cm_kpc, pc_kpc, kms_vcirc_Chemin
kmskpc_Om = kms_vcirc_Chemin/kpc_r_Chemin

err_kmskpc_Om = kms_vcirc_error_Chemin/kpc_r_Chemin


err_omega = err_kmskpc_Om*cm_km/cm_kpc


err_q = (kpc_r_Chemin/(kmskpc_Om*np.gradient(kpc_r_Chemin)))*((np.gradient(kmskpc_Om)/kmskpc_Om)-1)*err_kmskpc_Om


err_q = griddata(kpc_r_Chemin, err_q, kpc_r, method='linear',
                fill_value=nan, rescale=False)


err_omega = griddata(kpc_r_Chemin, err_omega, kpc_r,
                    method='linear', fill_value=nan, rescale=False)



#################################################################################################################################

#for m51
# from data.data_m51 import error_list, quant_list,radius_list,nandeleted_data, cm_km,cm_kpc

# kpc_r, dat_sigmatot, dat_sigmaHI,dat_sigmaH2, dat_q, dat_omega, dat_sigmasfr, T= nandeleted_data

# omega=quant_list[4]/(cm_km/cm_kpc) # undo unit corrected done in data file
# radius_omega=radius_list[4] # in kpc
# err_RC=error_list[1] # in km/s

# radius_temp=radius_list[-2] #in kpc
# err_temp=error_list[0]

# err_Om = err_RC/radius_omega
# err_q = (radius_omega/(omega*np.gradient(radius_omega)))*((np.gradient(omega)/omega)-1)*err_Om
# err_q = griddata(radius_omega, err_q, kpc_r, method='linear',
#                 fill_value=nan, rescale=False)
# err_omega = griddata(radius_omega, err_Om, kpc_r,
#                     method='linear', fill_value=nan, rescale=False)
#################################################################################################################################

#for m33

# from data.data_conv_M33 import kpc_r_kam, kms_vcirc_error_Kam, cm_km, cm_kpc, pc_kpc, kms_vcirc_Kam,T
# kmskpc_Om = kms_vcirc_Kam/kpc_r_kam
# err_kmskpc_Om = kms_vcirc_error_Kam/kpc_r_kam
# err_omega = err_kmskpc_Om*cm_km/cm_kpc
# err_q = (kpc_r_kam/(kmskpc_Om*np.gradient(kpc_r_kam)))*((np.gradient(kmskpc_Om)/kmskpc_Om)-1)*err_kmskpc_Om
# err_q = griddata(kpc_r_kam, err_q, kpc_r, method='linear',
#                 fill_value=nan, rescale=False)
# err_omega = griddata(kpc_r_kam, err_omega, kpc_r,
#                     method='linear', fill_value=nan, rescale=False)

#################################################################################################################################
#for 6946
# from data.data_6946 import error_list, quant_list,radius_list,nandeleted_data, cm_km,cm_kpc

# kpc_r, dat_sigmatot, dat_sigmaHI,dat_sigmaH2, dat_q, dat_omega, dat_sigmasfr, T= nandeleted_data

# omega=quant_list[4]/(cm_km/cm_kpc) # undo unit corrected done in data file
# radius_omega=radius_list[4] # in kpc
# err_RC=error_list[1] # in km/s

# radius_temp=radius_list[-2] #in kpc
# err_temp=error_list[0]

# err_Om = err_RC/radius_omega
# err_q = (radius_omega/(omega*np.gradient(radius_omega)))*((np.gradient(omega)/omega)-1)*err_Om
# err_q = griddata(radius_omega, err_q, kpc_r, method='linear',
#                 fill_value=nan, rescale=False)
# err_omega = griddata(radius_omega, err_Om, kpc_r,
#                     method='linear', fill_value=nan, rescale=False)
##############################################

# r=len(dat_sigmatot) # length of interpolated data
r=len(kpc_r)
relerr_q = err_q/np.abs(dat_q)
relerr_omega = err_omega/np.abs(dat_omega)
relerr_sigma = 0.01*np.ones(r)
relerr_sigmatot = 0.01*np.ones(r)
relerr_sigmasfr = 0.01*np.ones(r)

####################################################################################################################################
#for m31
err_T = (0.005*kpc_r + 0.1)*1e+4 #from Tabatabaei+13b equation ??
relerr_T = err_T/np.abs(T_tb)

#for m51
# os.chdir(current_directory+r'\data')
# err_T= griddata(radius_temp, err_temp, kpc_r, method='linear',fill_value=nan, rescale=False)
# relerr_T = err_T/np.abs(T)

#for 6946
# os.chdir(current_directory+r'\data')
# err_T= griddata(radius_temp, err_temp, kpc_r, method='linear',fill_value=nan, rescale=False)
# relerr_T = err_T/np.abs(T)

#for m33
# # err_T_OIII=(669*kpc_r + 470) #from Lin+17 eq 1&2
# err_T_OIII=(669*kpc_r + 470)*0 #when OIII data is excluded
# err_T_NII=(936*kpc_r + 316)
# err_T = (err_T_OIII+err_T_NII)/2 
# relerr_T = err_T/np.abs(T_tb)
####################################################################################################################################

#change scal_exponents file name depending on the model to be used
rel_err = np.array([relerr_q, relerr_omega, relerr_sigma, relerr_sigmatot,relerr_sigmasfr, relerr_T])
exps = np.load(current_directory+r'\data\scal_exponents.npy')
relerr_quan = np.sqrt(np.matmul(exps**2,rel_err**2))

err_quantities = model_f[1:]*relerr_quan

os.chdir(current_directory)

with open('errors_quan.pickle', 'wb') as f:
    pickle.dump(err_quantities, f)
print('Found the errors from the scaling relations')   
