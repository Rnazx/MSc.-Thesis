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

with open('mag_observables.pickle', 'rb') as f:
    kpc_r, h_f, l_f, u_f, cs_f, alphak_f, tau_f, biso_f, bani_f, Bbar_f, tanpB_f, tanpb_f = pickle.load(
        f)

os.chdir(current_directory + '\expressions')
import expressions.magnetic_expressions as m
import expressions.turbulence_expressions as t

q = Symbol('q')
omega = Symbol('\Omega')
sigma = Symbol('\Sigma')
sigmatot = Symbol('Sigma_tot')
sigmasfr = Symbol('Sigma_SFR')
T = Symbol('T')

os.chdir(current_directory + '\data')
with open('zip_data.pickle', 'rb') as f:
    kpc_r, data_pass = pickle.load(f)
r = kpc_r.size

from data.data_conv_M31 import kpc_r_Chemin, kms_vcirc_error_Chemin, cm_km, cm_kpc, pc_kpc

err_kmskpc_Om = kms_vcirc_error_Chemin/kpc_r_Chemin
err_omega = err_kmskpc_Om*cm_km/cm_kpc
err_q = -1 * kpc_r_Chemin / err_kmskpc_Om * np.gradient(err_kmskpc_Om)/np.gradient(kpc_r_Chemin)
err_q = griddata(kpc_r_Chemin, err_q, kpc_r, method='linear',
                fill_value=nan, rescale=False)
err_omega = griddata(kpc_r_Chemin, err_omega, kpc_r,
                    method='linear', fill_value=nan, rescale=False)



dat_sigmatot, dat_sigma, dat_sigmasfr, dat_q, dat_omega, zet, T_tb, psi, bet, ca, rk, mu = (
    np.array([data_pass[i][j] for i in range(r)]) for j in range(len(data_pass[0])))


relerr_q = err_q/dat_q
relerr_omega = err_omega/dat_omega
relerr_sigma = 0.01*np.ones(r)
relerr_sigmatot = 0.01*np.ones(r)
#relerr_sigmasfr = 0.01*np.ones(r)

err_T = (0.005*kpc_r + 0.1)*1e+4 #from Tabatabaei+13b equation ??
relerr_T = err_T/T_tb
observables = [q ,omega ,sigma,sigmatot ,T]
quantity = t.hg

exps_quan = np.array([scal_finder(t.hg, quantity, obs, data_pass, t.taue, t.alphak1, np.linspace(1,5000,100))[2]*np.ones(r) for obs in observables])
rel_err = np.array([relerr_q, relerr_omega, relerr_sigma, relerr_sigmatot, relerr_T])
relerr_h = np.sum(exps_quan*rel_err, axis = 0)
err_h = relerr_h*h_f
plt.errorbar(kpc_r, h_f*pc_kpc/cm_kpc, err_h*pc_kpc/cm_kpc, c='r', linestyle='-',ms=2, mew=2, capsize=2,
                  ecolor = 'y', marker='o', mfc='k', mec='k', label=r' $h$(pc)')

plt.xlabel(r'Radius (kpc)', fontsize=15)
plt.ylabel(r'Length scale (pc)', fontsize=15)
plt.legend(fontsize = 10)
plt.show()