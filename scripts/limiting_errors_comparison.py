import matplotlib

import numpy as np
import matplotlib.pyplot as plt
from sympy import *
import pickle
import os
from matplotlib.ticker import FormatStrFormatter
from scipy.interpolate import griddata
import sys

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
base_path = os.environ.get('MY_PATH')
galaxy_name = os.environ.get('galaxy_name')
sys.path.append(os.path.join(base_path,'src'))
from helper_functions import datamaker, parameter_read

params = parameter_read(os.path.join(base_path,'inputs','parameter_file.in'))
sys.path.append(os.path.join(base_path,'data','supplementary_data', galaxy_name))
from observables import *

########################################################################################################################
# subprocess.run(["python", "zipped_data.py"])
# subprocess.run(["python", "get_magnetic_observables.py"])

# conversion factors
pc_kpc = 1e3  # number of pc in one kpc
cm_km = 1e5  # number of cm in one km
cm_kpc = 3.086e+21  # number of centimeters in one parsec
s_Myr = 1e+6*(365*24*60*60)  # megayears to seconds
deg_rad = 180e0/np.pi
arcmin_deg = 60e0
arcsec_deg = 3600e0

########################################################################################################
current_directory = str(os.getcwd())


os.chdir(os.path.join(base_path,'inputs'))

with open('zip_data.in', 'rb') as f:
    kpc_r, data_pass = pickle.load(f)

#######################################################################################################################################

dat_u = griddata(kpc_radius, np.sqrt(3)*kms_sigmaLOS, kpc_r, method='linear',
                 fill_value=nan, rescale=False)*1e+5
try:
    dat_u_warp = griddata(kpc_radius, np.sqrt(3)*kms_sigmaLOS_warp, kpc_r, method='linear',
                    fill_value=nan, rescale=False)*1e+5
except NameError:
    pass
os.chdir(current_directory)

m = 2
dm = 2.5
fs = 15
lfs = 10
leg_textsize = 10
axis_textsize = 10
rc = {"font.family" : "serif", 
      "mathtext.fontset" : "stix"}
plt.rcParams.update(rc)
plt.rcParams["font.serif"] = ["Times New Roman"] + plt.rcParams["font.serif"]
matplotlib.rc('xtick', labelsize=fs)
matplotlib.rc('ytick', labelsize=fs)
matplotlib.ticker.AutoMinorLocator(n=None)
plt.rcParams["xtick.minor.visible"] =  True
plt.rcParams["ytick.minor.visible"] =  True
plt.rcParams["legend.loc"] = 'upper right'
plt.rcParams["errorbar.capsize"] = 2

def axis_pars(ax):
    #ax.xaxis.set_ticks(np.arange(6, 20, 2))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%g'))
    ax.tick_params(axis='both', which='minor',
                   labelsize=axis_textsize, colors='k', length=3, width=1)
    ax.tick_params(axis='both', which='major',
                   labelsize=axis_textsize, colors='k', length=5, width=1.25)
    # ax.spines['top'].set_visible(False)
    # ax.spines['right'].set_visible(False)
    # ax.spines['bottom'].set_linewidth(2)
    # ax.spines['left'].set_linewidth(2)
    ax.legend(fontsize=lfs, frameon=False, handlelength=4, ncol=1, prop={
            'size': leg_textsize, 'family': 'Times New Roman'}, fancybox=True, framealpha=0.9, handletextpad=0.7, columnspacing=0.7)
    
def fill_error(ax, quan_f, quan_err, color = 'red', alpha = 0.2, error_exists = True):
    if error_exists:
        ax.fill_between(kpc_r, (quan_f+quan_err), (quan_f-quan_err)
                        , alpha=alpha, edgecolor='k', facecolor=color, where = None, interpolate=True)
    else:
        return


os.chdir(os.path.join(base_path,'outputs'))

with open(f'{galaxy_name}output_ca_'+str(params[r'C_\alpha'])+'K_'+str(params[r'K'])+'z_'+str(params[r'\zeta'])+'psi_'+str(params[r'\psi'])+'b_'+str(params[r'\beta'])+'.out', 'rb') as f:
    kpc_r, h_f, l_f, u_f, cs_f, alphak_f, taue_f, taur_f, biso_f, bani_f, Bbar_f, tanpB_f, tanpb_f , dkdc_f = pickle.load(
        f)


fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7, 10), tight_layout=True)
fig.suptitle(r'$C_\alpha$ = '+str(params[r'C_\alpha'])+r'    $K$ = '+str(params[r'K'])+
             r'    $\zeta$ = '+str(params[r'\zeta'])+r'    $\psi$ = '+str(params[r'\psi'])+r'    $\beta$ = '+str(params[r'\beta']), weight = 15)
i = 0
ax.plot(kpc_r, h_f*pc_kpc/cm_kpc, c='r', linestyle='-', mfc='k',
              mec='k', markersize=m, marker='o', label=r' $h$(pc)')
try:
    ax.plot(kpc_dat_r, pc_dat_h, c='b', linestyle='dotted', 
                marker='*',mfc='y',mec='b',mew=1, markersize = 7, label=r'Fiducial values from Chamandy et.al.(2016) $h(pc)$')
except NameError:
    pass
ax.plot(kpc_r, l_f*pc_kpc/cm_kpc, c='g',
              linestyle='-', mfc='k', mec='k', markersize=m, marker='o', label=r'Correlation length l(pc)')
# ax.plot(kpc_r, datamaker(lsn , data_pass, h_f, tau_f)*pc_kpc/cm_kpc,c = 'y',linestyle='--',mfc='k',mec='k', marker='o')
ax.axhline(y=100, color='black', linestyle='--', alpha = 0.2)
#ax.set_yticks(list(plt.yticks()[0])+[100])
with open('errors_subsonic.out', 'rb') as f:
    h_err, l_err, u_err, cs_err, alphak_err, tau_err, taur_err, biso_err, bani_err, Bbar_err, tanpB_err, tanpb_err, dkdc_err = pickle.load(
        f)
try:
    fill_error(ax, h_f*pc_kpc/cm_kpc, h_err*pc_kpc/cm_kpc,'blue')
except NameError:
    pass
with open('errors_supersonic.out', 'rb') as f:
    h_err, l_err, u_err, cs_err, alphak_err, tau_err, taur_err, biso_err, bani_err, Bbar_err, tanpB_err, tanpb_err, dkdc_err = pickle.load(
        f)
try:
    fill_error(ax, h_f*pc_kpc/cm_kpc, h_err*pc_kpc/cm_kpc,'green')
except NameError:
    pass
axis_pars(ax)

    
#ax.set_xlabel(r'Radius (kpc)', fontsize=fs)
ax.set_ylabel(r'Length scale (pc)', fontsize=fs)

plt.show()
os.chdir(current_directory)
