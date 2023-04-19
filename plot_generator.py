import matplotlib
from helper_functions import datamaker, exp_analytical_data
from scipy.integrate import quad
from scipy.interpolate import griddata
import numpy as np
import matplotlib.pyplot as plt
from sympy import *
from fractions import Fraction
import pickle
import os
from matplotlib.ticker import FormatStrFormatter
import subprocess

# Defining the Observables
q = Symbol('q')
omega = Symbol('\Omega')
sigma = Symbol('\Sigma')
sigmatot = Symbol('Sigma_tot')
sigmasfr = Symbol('Sigma_SFR')
T = Symbol('T')


# Defining the Constants
gamma = Symbol('gamma')
boltz = Symbol('k_B')
mu = Symbol('mu')
mh = Symbol('m_H')


# Defining the general parameters
u = Symbol('u')
tau = Symbol('tau')
l = Symbol('l')
h = Symbol('h')
##############################################################################################################################
#reading the parameters
params = {}
with open('parameter_file.in', 'r') as FH:
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
        except ValueError:
            num, denom = value.split('/')
            params[par_name] = np.float64(num) / np.float64(denom)
########################################################################################################################

#subprocess.run(["python", "get_magnetic_observables.py"])


#####################################################################################################################
# Magnetic field data from Beck et. al. 2019 Tables 3 and 4
G_dat_Bord = np.array([4.9, 5.2, 4.9, 4.6])
G_dat_Breg = np.array([1.8, 2.1, 2.6, 2.7])
G_dat_Btot = np.array([7.3, 7.5, 7.1, 6.3])
RM_dat_po = np.array([30, 29, 26, 27, 27]) * np.pi/180
err_RM_dat_po = np.array([5, 4, 3, 2, 3]) * np.pi/180
rmdat_tanpo = np.tan(RM_dat_po)
rm_errdat_tanpo = 1/(np.cos(err_RM_dat_po))**2
mrange_endps = np.array([6.8, 9.0, 11.3, 13.6, 15.8])
mrange = (mrange_endps[1:] + mrange_endps[:-1])/2

M_dat_pb = np.array([13, 19, 11, 8]) * np.pi/180
err_M_dat_pb = np.array([4, 3, 3, 3]) * np.pi/180
RM_dat_pb = np.array([4, 9, 7, 7, 5]) * np.pi/180
err_RM_dat_pb = np.array([5, 3, 3, 2, 3]) * np.pi/180
rmdat_tanpb = np.tan(RM_dat_pb)
rm_errdat_tanpb = 1/(np.cos(err_RM_dat_pb))**2
m_errdat_tanpb = 1/(np.cos(err_M_dat_pb))**2

rmrange = np.arange(7.5, 12, 1)
kpc_dat_r = np.array([7, 9, 11, 13])

pc_dat_h = np.array([316.4, 371.9, 437.1, 513.7])
#######################################################################################################################################
# conversion factors
pc_kpc = 1e3  # number of pc in one kpc
cm_km = 1e5  # number of cm in one km
cm_kpc = 3.086e+21  # number of centimeters in one parsec
cm_kpc = 3.086e+21  # number of ccentimeters in one parsec
s_Myr = 1e+6*(365*24*60*60)  # megayears to seconds
deg_rad = 180e0/np.pi
arcmin_deg = 60e0
arcsec_deg = 3600e0


########################################################################################################
current_directory = str(os.getcwd())



with open('mag_observables.pickle', 'rb') as f:
    kpc_r, h_f, l_f, u_f, cs_f, alphak_f, tau_f, biso_f, bani_f, Bbar_f, tanpB_f, tanpb_f = pickle.load(
        f)
with open('errors_quan.pickle', 'rb') as f:
    h_err, l_err, u_err, cs_err, alphak_err, tau_err, biso_err, bani_err, Bbar_err, tanpB_err, tanpb_err = pickle.load(
        f)
os.chdir(current_directory + '\data')

with open('zip_data.pickle', 'rb') as f:
    kpc_r, data_pass = pickle.load(f)

# with open('data_m31.pickle', 'rb') as f:
#     data_m31 = pickle.load(f)
# kpc_r, dat_sigmatot, dat_sigma, dat_q, dat_omega, dat_sigmasfr, molfrac = data_m31
#######################################################################################################################################

os.chdir(current_directory + '\expressions')
from expressions.magnetic_expressions import Dk, Dc


dkdc_f = datamaker((Dk/Dc), data_pass, h_f, tau_f, alphak_f)
alpham_f = alphak_f*((1/dkdc_f)-1)

alphasat_f = alphak_f + alpham_f


#####################################################################################################################################

kms_sigmaLOS = np.array([17.10, 8.58, 12.94, 17.12, 12.99, 10.03, 13.13, 15.38, 16.60, 14.95, 15.63, 14.50, 12.63, 12.25, 12.57, 12.51, 13.53, 13.79, 13.93, 13.21, 11.56, 10.87, 10.81, 11.06, 11.77, 11.12, 12.21, 12.12, 12.12, 11.70, 11.41, 11.70, 12.54, 12.15, 11.65, 11.18, 11.21, 11.36, 11.72, 10.62, 10.23, 10.47, 10.62, 10.65, 9.76, 10.82,
                        11.45, 10.92, 10.91, 10.46, 10.01, 10.04, 9.86, 9.94, 10.17, 9.86, 9.55, 9.72, 9.79, 10.15, 9.82, 9.76, 9.86, 9.37, 9.73, 9.58, 9.66, 10.21, 10.21, 9.20, 8.95, 9.59, 10.20, 8.94, 9.01, 9.12, 9.29, 9.46, 9.80, 8.34, 9.01, 8.46, 8.98, 7.70, 7.99, 7.04, 6.43, 5.80, 8.92, 8.76, 6.95, 6.81, 6.60, 9.30, 13.69, 11.38, 6.81, 4.58, 9.05, 10.89])
kms_sigmaLOS_warp = np.array([12.05, 15.64, 15.70, 10.75, 8.71, 5.82, 6.08, 8.23, 9.48, 11.61, 13.61, 12.37, 7.77, 6.22, 5.98, 5.38, 5.47, 6.62, 8.16, 8.31, 8.67, 9.30, 9.67, 8.51, 7.51, 6.32, 6.13, 6.52, 9.29, 10.66, 11.48, 10.77, 10.61, 11.23, 12.05, 11.80, 11.42, 10.91, 8.56, 6.32, 5.45, 4.81, 5.43, 6.68, 10.43, 10.71,
                             10.95, 10.25, 9.16, 8.87, 9.05, 8.72, 8.98, 8.49, 9.61, 6.84, 6.16, 6.88, 9.87, 9.87, 9.85, 8.62, 7.97, 7.54, 7.19, 7.99, 8.30, 8.56, 7.14, 7.23, 7.85, 7.05, 7.49, 6.64, 8.34, 7.07, 5.53, 5.61, 4.06, 4.75, 8.15, 7.09, 5.89, 6.24, 5.59, 5.59, 4.49, 3.81, 2.61, 3.57, 4.67, 5.37, 6.67, 4.19, 5.09, 4.45, 4.40, 4.53, 5.62, 5.20])

arcsec_r = np.array([100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0, 1100.0, 1200.0, 1300.0, 1400.0, 1500.0, 1600.0, 1700.0, 1800.0, 1900.0, 2000.0, 2100.0, 2200.0, 2300.0, 2400.0, 2500.0, 2600.0, 2700.0, 2800.0, 2900.0, 3000.0, 3100.0, 3200.0, 3300.0, 3400.0, 3500.0, 3600.0, 3700.0, 3800.0, 3900.0, 4000.0, 4100.0, 4200.0, 4300.0, 4400.0, 4500.0, 4600.0, 4700.0, 4800.0, 4900.0, 5000.0,
                    5100.0, 5200.0, 5300.0, 5400.0, 5500.0, 5600.0, 5700.0, 5800.0, 5900.0, 6000.0, 6100.0, 6200.0, 6300.0, 6400.0, 6500.0, 6600.0, 6700.0, 6800.0, 6900.0, 7000.0, 7100.0, 7200.0, 7300.0, 7400.0, 7500.0, 7600.0, 7700.0, 7800.0, 7900.0, 8000.0, 8100.0, 8200.0, 8300.0, 8400.0, 8500.0, 8600.0, 8700.0, 8800.0, 8900.0, 9000.0, 9100.0, 9200.0, 9300.0, 9400.0, 9500.0, 9600.0, 9700.0, 9800.0, 9900.0, 10000.0])

kpc_gal_dist = 780e0  # from Beck et al
kpc_radius = kpc_gal_dist*arcsec_r/(arcsec_deg*deg_rad)

dat_u = griddata(kpc_radius, np.sqrt(3)*kms_sigmaLOS, kpc_r, method='linear',
                 fill_value=nan, rescale=False)*1e+5
dat_u_warp = griddata(kpc_radius, np.sqrt(3)*kms_sigmaLOS_warp, kpc_r, method='linear',
                 fill_value=nan, rescale=False)*1e+5

#########################################################################################################################################################


os.chdir(current_directory)

#from get_magnetic_observables import omt, kah, taue_f, taur_f

pB = np.arctan(-tanpB_f)
pB_err = -tanpB_err/(1+tanpB_f**2)
pbb = np.arctan(tanpb_f)
pbb_err = tanpb_err/(1+tanpb_f**2)
pbo = (1/2)*((1+(2*Bbar_f*bani_f*np.cos(pbb-pB))/
              (bani_f**2+Bbar_f**2))*np.arctan((Bbar_f*np.sin(pB) + bani_f*np.sin(pbb))/
                                               ((Bbar_f*np.cos(pB)) + bani_f*np.cos(pbb)))
                                                        + (1-(2*Bbar_f*bani_f*np.cos(pbb-pB))/
                                                           (bani_f**2+Bbar_f**2))*np.arctan((Bbar_f*np.sin(pB) - bani_f*np.sin(pbb))/
                                                                                            ((Bbar_f*np.cos(pB)) - bani_f*np.cos(pbb))))

def pogen(b, B, pb, pB, s):
    return (np.exp(-b**2/(2*s**2))/
            (np.sqrt(2*(np.pi))*s))*(1+(2*B*b*np.cos(pb-pB))/
                                     (b**2 + B**2))*np.arctan((B*np.sin(pB) + b*np.sin(pb))/
                                                                                            ((B*np.cos(pB)) + b*np.cos(pb)))

brms = np.sqrt(np.average(bani_f**2))
brms
interval = 1e+3
pog = np.array([quad(pogen, -interval, interval, args=(Bbar_f[i], pbb[i], pB[i], brms),
               points=[-interval*brms, interval*brms])[0] for i in range(len(kpc_r))])

G_scal_Bbartot = np.sqrt(biso_f**2 + bani_f**2 + Bbar_f**2)
G_scal_Bbarreg = Bbar_f
G_scal_Bbarord = np.sqrt(bani_f**2 + Bbar_f**2)
print(h_err, l_err, u_err, cs_err, alphak_err, tau_err, biso_err, bani_err, Bbar_err, tanpB_err, tanpb_err)
G_scal_Bbartot_err = (biso_err*biso_f + bani_err*bani_f + Bbar_err*Bbar_f)/G_scal_Bbartot
G_scal_Bbarreg_err = Bbar_err
G_scal_Bbarord_err = (bani_err*bani_f + Bbar_err*Bbar_f)/G_scal_Bbarord

####################################################################################################################################################
m = 2
fs = 15
lfs = 10
matplotlib.rc('xtick', labelsize=fs)
matplotlib.rc('ytick', labelsize=fs)
matplotlib.ticker.AutoMinorLocator(n=None)
plt.rcParams["xtick.minor.visible"] =  True
plt.rcParams["ytick.minor.visible"] =  True
plt.rcParams["legend.loc"] = 'upper right'
plt.rcParams["errorbar.capsize"] = 2

os.chdir(current_directory+'\plots')

from matplotlib.backends.backend_pdf import PdfPages
# PDF = PdfPages('ca_'+str(params[r'C_\alpha'])+'rk_'+str(params[r'R_\kappa'])+'z_'+str(
#      params[r'\zeta'])+'psi_'+str(params[r'\psi'])+'f_'+str(params[r'\beta'])+'.pdf')#('plots_model'+str(model_no)+let+'t_vary_'+'ca_'+str(ca)+'rk_'+str(rk)+'z_'+str(z.mean())+'.pdf')
fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(10, 10), tight_layout=True)
i = 0
j = 0
#print(l_err*pc_kpc/cm_kpc)
ax[i][j].errorbar(kpc_r, h_f*pc_kpc/cm_kpc, h_err*pc_kpc/cm_kpc,ecolor = 'y', c='r', linestyle='-', mfc='k',
              mec='k', markersize=m, marker='o', label=r' $h$(pc)')
ax[i][j].plot(kpc_dat_r, pc_dat_h, c='b', linestyle='--', mfc='k', mec='k',
              markersize=m, marker='o', label=r'Data from Chamandy et.al.(16) $h(pc)$')
ax[i][j].errorbar(kpc_r, l_f*pc_kpc/cm_kpc,l_err*pc_kpc/cm_kpc, c='g',
              linestyle='-', mfc='k', mec='k', markersize=m, marker='o', label=r'Correlation length l(pc)')
# ax[i][j].plot(kpc_r, datamaker(lsn , data_pass, h_f, tau_f)*pc_kpc/cm_kpc,c = 'y',linestyle='--',mfc='k',mec='k', marker='o')
ax[i][j].xaxis.set_ticks(np.arange(6, 20, 2))
ax[i][j].xaxis.set_major_formatter(FormatStrFormatter('%g'))
ax[i][j].set_xlabel(r'Radius (kpc)', fontsize=fs)
ax[i][j].set_ylabel(r'Length scale (pc)', fontsize=fs)
ax[i][j].legend(fontsize = lfs)

j = 1
ax[i][j].errorbar(kpc_r, u_f/cm_km, u_err/cm_km, color='y', marker='o', mfc='k',
              mec='k', markersize=m, label=r'$u$')
ax[i][j].plot(kpc_r, alphak_f/cm_km, color='b', marker='o',
              mfc='k', mec='k', markersize=m, label=r'$\alpha_k$')
# ax[i][j].plot(kpc_r, alpham_f/cm_km, color='r', marker='o',
#               mfc='k', mec='k', markersize=m, label=r'$\alpha_m$')
ax[i][j].plot(kpc_r, alphasat_f/cm_km, color='m', marker='o',
              mfc='k', mec='k', markersize=m, label=r'$\alpha_{sat}$')


# ax[i][j].plot(kpc_r, np.sqrt(u_f**2 + cs_f**2) /
#               cm_km, color='r', marker='o', mfc='k', mec='k', markersize=m, label=r'$\sqrt{u^2+c_s^2}$')
ax[i][j].errorbar(kpc_r, cs_f /
              cm_km, cs_err /cm_km, color='g', linestyle='--', label=r'$c_s$')
ax[i][j].plot(kpc_r, dat_u/cm_km, marker='o', markersize=m,
              c='tab:orange', mfc='k', mec='k', linestyle='--', label='data')
ax[i][j].plot(kpc_r, dat_u_warp/cm_km, marker='o', markersize=m,
              c='tab:cyan', mfc='k', mec='k', linestyle='--', label='data (warp)', alpha = 0.5)
ax[i][j].set_ylim(0)
ax[i][j].xaxis.set_ticks(np.arange(6, 20, 2))
ax[i][j].xaxis.set_major_formatter(FormatStrFormatter('%g'))
ax[i][j].set_xlabel(r'Radius ($kpc$)', fontsize=fs)
ax[i][j].set_ylabel(r'Turbulence parameters (km/s)',  fontsize=fs)
ax[i][j].legend(fontsize=8)

i = 1

ax[i][j].plot(mrange, G_dat_Btot, c='b', linestyle='--', marker='o', mfc='k',
              mec='k', markersize=m)#, label='Average Binned data $B_{tot}$ ($\mu G$)')
ax[i][j].errorbar(kpc_r, G_scal_Bbartot*1e+6,G_scal_Bbartot_err*1e+6, c='b', linestyle='-', marker='o', mfc='k', mec='k',
              markersize=m, label=r' $B_{tot}=\sqrt{\bar{B}^2+b_{iso}^2+b_{ani}^2}$')

ax[i][j].plot(mrange, G_dat_Breg, c='r', linestyle='--', marker='o', mfc='k',
              mec='k', markersize=m)#, label='Average Binned data $B_{reg}$ ($\mu G$)')
ax[i][j].errorbar(kpc_r, G_scal_Bbarreg*1e+6,G_scal_Bbarreg_err*1e+6, c='r', linestyle='-', marker='o',
              mfc='k', mec='k', markersize=m, label=r' $B_{reg} = \bar{B}$')

ax[i][j].plot(mrange, G_dat_Bord, c='y', linestyle='--', marker='o', mfc='k',
              mec='k', markersize=m)#, label='Average Binned data $B_{ord}$ ($\mu G$)')
ax[i][j].errorbar(kpc_r, G_scal_Bbarord*1e+6,G_scal_Bbarord_err*1e+6, c='y', linestyle='-', marker='o', mfc='k',
              mec='k', markersize=m, label=r' $B_{ord} = \sqrt{\bar{B}^2+b_{ani}^2}$')

ax[i][j].set_xlabel(r'Radius ($kpc$)', fontsize=fs)
ax[i][j].xaxis.set_ticks(np.arange(6, 20, 2))
ax[i][j].xaxis.set_major_formatter(FormatStrFormatter('%g'))
ax[i][j].set_ylabel('Magnetic field strength ($\mu G$)', fontsize=fs)
ax[i][j].legend(fontsize = lfs)
# ax[i][j].plot(kpc_r, biso_f*1e+6,c = 'b',linestyle='-', marker='o',markersize=m,label = r' $b_{iso}$')

# ax[i][j].plot(mrange, G_dat_Breg, c = 'r',linestyle='--', marker='o',mfc='k',mec='k',markersize=m, label = r'Data $B_{reg}$ ($\mu G$)')
# ax[i][j].plot(mrange, np.sqrt(G_dat_Btot**2 - G_dat_Bord**2), c = 'b',mfc='k',mec='k',markersize=m,linestyle='--', marker='o', label = r'Data $b_{iso}$ ($\mu G$)')
# ax[i][j].plot(mrange, np.sqrt(G_dat_Bord**2 - G_dat_Breg**2), c = 'y',mfc='k',mec='k',markersize=m,linestyle='--', marker='o', label = r'Data $b_{ani}$ ($\mu G$)')
# ax[i][j].plot(kpc_r, Bbar_f*1e+6,c = 'r',linestyle='-', marker='o',mfc='k',mec='k',markersize=m,label = r' $B_{reg} = \bar{B}$')
# ax[i][j].plot(kpc_r, bani_f*1e+6,c = 'y',linestyle='-', marker='o',mfc='k',mec='k',markersize=m,label = r' $b_{ani}$')
# ax[i][j].set_xlabel(r'Radius ($kpc$)', fontsize = fs)
# ax[i][j].set_ylabel('Magnetic field ($\mu G$)')
# ax[i][j].legend(fontsize = lfs)

j = 0
# ax[i][j].errorbar(mrange, 180*M_dat_pb/np.pi, yerr = 180*err_M_dat_pb/np.pi, ms=4, mew=2, capsize=2, c = 'b',linestyle='--', marker='o',mfc='k',mec='k', label = r'Data $p_{B}$ (mean field)(M)', alpha = 0.5)
ax[i][j].errorbar(rmrange, 180*RM_dat_pb/np.pi, yerr=180*err_RM_dat_pb/np.pi, ms=2, mew=2, capsize=2,
                  c='r', linestyle='--', marker='o', mfc='k', mec='k')#, label=r'Data $p_{B}$ (mean field)(RM)')
ax[i][j].errorbar(rmrange, 180*RM_dat_po/np.pi, yerr=180*err_RM_dat_po/np.pi, ms=2, mew=2, capsize=2,
                  c='y', linestyle='--', marker='o', mfc='k', mec='k')#, label=r'Data $p_{o}$ (ordered field)(RM)')

ax[i][j].errorbar(kpc_r, 180*pB/np.pi,180*pB_err/np.pi, c='r', linestyle='-', marker='o',
              markersize=m, mfc='k', mec='k', label=r' $p_{B}$ (mean field)')
# ax[i][j].plot(kpc_r, 180*pbb/np.pi,c = 'r',linestyle='--', marker='o',label = r' $p_{b}$ (anisotropic field)')
ax[i][j].plot(kpc_r, 180*pog/np.pi, c='y', linestyle='-', mfc='k', markersize=m,
              mec='k', marker='o', label=r' $p_{o}$ (ordered field)')

ax[i][j].set_xlabel(r'Radius ($kpc$)', fontsize=fs)
ax[i][j].xaxis.set_ticks(np.arange(6, 20, 2))
ax[i][j].xaxis.set_major_formatter(FormatStrFormatter('%g'))
ax[i][j].set_ylabel(r'Pitch angle (deg)', fontsize=fs)
ax[i][j].legend(fontsize = lfs)
# PDF.savefig(fig)
# fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(10, 10), tight_layout=True)
# i = 0
# j=0
# ax[i][j].axhline(y=1, color='g', linestyle='-', label=r'1')
# ax[i][j].plot(kpc_r, omt, marker='o', markersize=m,
#               c='tab:orange', mfc='k', mec='k', label=r'$\Omega\tau$')
# ax[i][j].plot(kpc_r, kah, marker='o',
#               markersize=m, c='tab:blue', mfc='k', mec='k', label=r'$\frac{K_\alpha h}{C_\alpha \tau u}$')

# ax[i][j].set_xlabel('Radius(kpc)', fontsize=fs)
# ax[i][j].set_ylabel(r'Condition for $\alpha_k$ (Myr)', fontsize=fs)
# ax[i][j].legend(fontsize = lfs)

# j = 1
# ax[i][j].plot(kpc_r, taue_f/s_Myr, c='b', markersize=m,
#               linestyle='-', marker='o', mfc='k', mec='k', label=r'$\tau^e$')
# ax[i][j].plot(kpc_r, taur_f/s_Myr, c='g',
#               markersize=m, linestyle='-', marker='o', mfc='k', mec='k', label=r'$\tau^r$')
# ax[i][j].plot(kpc_r, h_f/(u_f*s_Myr), c='y', markersize=m,
#               linestyle='-', marker='o', mfc='k', mec='k', label=r'$h/u$')

# ax[i][j].set_xlabel('Radius(kpc)', fontsize=fs)
# ax[i][j].set_ylabel(r'Correlation Time $\tau$ (Myr)', fontsize=fs)
# ax[i][j].legend(fontsize = lfs)

# i = 1
# ax[i][j].plot(kpc_r, dkdc_f, markersize=m, linestyle='-',
#               marker='o', mfc='k', mec='k', label=r'$D_k/D_c$')
# ax[i][j].plot(kpc_r, 1*np.ones(len(kpc_r)))
# ax[i][j].set_xlabel('Radius(kpc)', fontsize=fs)
# ax[i][j].set_ylabel(r'$D_k/D_c$',  fontsize=fs)
# ax[i][j].legend(fontsize = lfs)

# j = 0
# ax[i][j].plot(kpc_r, (((np.pi**2)*(tau_f*(u_f**2))/3*(np.sqrt(dkdc_f)-1)/(4*h_f**2))**(-1))/
#               (s_Myr*1e+3), c='g', markersize=m, linestyle='-', marker='o', mfc='k', mec='k', label=r'$\gamma$')
# ax[i][j].set_xlabel('Radius(kpc)', fontsize=fs)
# ax[i][j].set_ylabel(r'local e-folding time $\gamma$ (Gyr)', fontsize=fs)
# ax[i][j].legend(fontsize = lfs)
# #fig.suptitle(r'$T_{avg} = $ '+str(T_tb.mean()), fontsize=15)
# PDF.savefig(fig)
# PDF.close()

# plt.savefig('ca_'+str(params[r'C_\alpha'])+'rk_'+str(params[r'R_\kappa'])+'z_'+str(
#     params[r'\zeta'])+'psi_'+str(params[r'\psi'])+'f_'+str(params[r'\beta'])+'.png')
plt.show()
