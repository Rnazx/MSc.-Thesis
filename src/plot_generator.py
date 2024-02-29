import matplotlib
from helper_functions import datamaker, parameter_read
import numpy as np
import matplotlib.pyplot as plt
from sympy import *
import pickle
import os
from matplotlib.ticker import FormatStrFormatter
from scipy.interpolate import griddata
import sys
from datetime import date

today = date.today()

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
g=galaxy_name
params = parameter_read(os.path.join(base_path,'inputs','parameter_file.in'))
switch = parameter_read(os.path.join(base_path,'inputs','switches.in'))

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

os.chdir(os.path.join(base_path,'outputs'))

with open(f'{galaxy_name}output_ca_'+str(params[r'C_\alpha'])+'K_'+str(params[r'K'])+'z_'+str(params[r'\zeta'])+'psi_'+str(params[r'\psi'])+'b_'+str(params[r'\beta'])+'.out', 'rb') as f:
    kpc_r, h_f, l_f, u_f, cs_f, alphak_f, taue_f, taur_f, biso_f, bani_f, Bbar_f, tanpB_f, tanpb_f , dkdc_f = pickle.load(
        f)
with open('errors_subsonic.out', 'rb') as f:
        subsonic_errors= pickle.load(f)
with open('errors_supersonic.out', 'rb') as f:
        supersonic_errors= pickle.load(f)
h_err, l_err, u_err, cs_err, alphak_err, tau_err, \
        taur_err, biso_err, bani_err, Bbar_err, \
                tanpB_err, tanpb_err, dkdc_err = [np.maximum(sub, sup) for sub,sup in zip(subsonic_errors, supersonic_errors)]
# h_f=(10**(-3))*h_f

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

from helper_functions import pitch_angle_integrator

pB, po, pb, pB_err, po_err, pb_err = pitch_angle_integrator(kpc_r, tanpB_f,tanpb_f, \
                                   Bbar_f, bani_f, tanpB_err,tanpb_err, Bbar_err, bani_err)


G_scal_Bbartot = np.sqrt(biso_f**2 + bani_f**2 + Bbar_f**2)
G_scal_Bbarreg = Bbar_f
G_scal_Bbarord = np.sqrt(bani_f**2 + Bbar_f**2)


G_scal_Bbartot_err = np.sqrt((biso_err*biso_f )**2+ (bani_err*bani_f)**2 + (Bbar_err*Bbar_f)**2)/G_scal_Bbartot
G_scal_Bbarreg_err = Bbar_err
G_scal_Bbarord_err = np.sqrt((bani_err*bani_f)**2 + (Bbar_err*Bbar_f)**2)/G_scal_Bbarord

for i in range(len(kpc_r)):
    print('#####################################')
    print(i,G_scal_Bbarreg_err[i],G_scal_Bbarord_err[i],G_scal_Bbartot_err[i])
m = 9 #marker size
lw = 3
dm = 2.5
fs = 20
lfs = 10
leg_textsize = 20
axis_textsize = 20
hd=1.6 #handlelength: changes length of line in legend
legend_labelspace= 0.17 #handletextpad: gap between label and symbol in legend


rc = {"font.family" : "serif", "mathtext.fontset" : "stix"}
plt.rcParams.update(rc)
plt.rcParams["font.serif"] = ["Times New Roman"] + plt.rcParams["font.serif"]
matplotlib.rc('xtick', labelsize=fs)
matplotlib.rc('ytick', labelsize=fs)
matplotlib.ticker.AutoMinorLocator(n=None)
plt.rcParams['axes.xmargin'] = 0
plt.rcParams['axes.ymargin'] = 0
plt.rcParams["xtick.minor.visible"] =  True
plt.rcParams["ytick.minor.visible"] =  True
plt.rcParams["legend.loc"] = 'upper right'
plt.rcParams["errorbar.capsize"] = 2

def axis_pars(ax):
    if g=='m31':
        ax.xaxis.set_ticks(np.arange(6, 18, 1)) #for m31
    elif g=='m33':
        ax.xaxis.set_ticks(np.arange(0, 10, 1)) #for m33
    elif g=='m51':
        ax.xaxis.set_ticks(np.arange(0, 11, 1)) #for m51
    else:
        ax.xaxis.set_ticks(np.arange(0, 21, 2)) #for ngc

    ax.xaxis.set_major_formatter(FormatStrFormatter('%g'))
    ax.tick_params(axis='both', which='minor',
                   labelsize=axis_textsize, colors='k', length=3, width=1)

    ax.tick_params(axis='both', which='major',
                   labelsize=axis_textsize, colors='k', length=5, width=1.25)
    # ax.spines['top'].set_visible(False)
    # ax.spines['right'].set_visible(False)
    # ax.spines['bottom'].set_linewidth(2)
    # ax.spines['left'].set_linewidth(2)

    #next line is commented as legend location needs to be customised for each plot
    # ax.legend(fontsize=lfs, frameon=False, handlelength=4, ncol=1, prop={
    #         'size': leg_textsize, 'family': 'Times New Roman'}, fancybox=True, framealpha=0.9, handletextpad=0.7, columnspacing=0.7)
    
def fill_error(ax, quan_f, quan_err, color = 'red', alpha = 0.2, error_exists = True):
    if error_exists:
        ax.fill_between(kpc_r, (quan_f+quan_err), (quan_f-quan_err)
                        , alpha=alpha, edgecolor='k', facecolor=color, where = None, interpolate=True)
    else:
        return

save_files_dir=current_directory+r'\{},moldat_{},taue,z_{},psi_{},ca_{},beta_{},A_{}'.format(str(today),switch['incl_moldat'],params[r'\zeta'],params[r'\psi'],
                                params[r'C_\alpha'],params[r'\beta'],params['A'])

#so that new directories wont be made when the directory name is imported from this file
if __name__ == '__main__':
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

##################################################################################################################

fig,ax = plt.subplots(nrows=1, ncols=1, figsize=(7, 5), tight_layout=True)

h=h_f*pc_kpc/cm_kpc
l=l_f*pc_kpc/cm_kpc

# ax.plot(kpc_r, h, c='r', linestyle='-', mfc='k',mec='k', markersize=m, marker='o', label=r'$h$ (pc)')
ax.plot(kpc_r, h_f*pc_kpc/cm_kpc, c='b', linestyle='-', linewidth=lw, label=r'Scale height')


# ax.plot(kpc_r, l, c='g',linestyle='-.', mfc='k', mec='k', markersize=m, marker='o', label=r'$l$ (pc)')
ax.plot(kpc_r, l_f*pc_kpc/cm_kpc, c='g',linestyle='-.', linewidth=lw, label=r'Correlation length')

#mfc- marker fill color, mec- marker edge color
ax.plot(kpc_dat_r, pc_dat_h, zorder=2,linestyle=' ',marker='*',c='b',mfc='b',mec='k',mew=1, markersize = 13, label=r'Milky Way scaling (C16)')

# ax.plot(kpc_r, datamaker(lsn , data_pass, h_f, tau_f)*pc_kpc/cm_kpc,c = 'y',linestyle='--',mfc='k',mec='k', marker='o')
# ax.axhline(y=100, color='black', linestyle=':', alpha = 1)
# ax.set_yticks(list(plt.yticks()[0])+[100])
axis_pars(ax)

h_err_corr_units= h_err*pc_kpc/cm_kpc
try:
    fill_error(ax, h_f*pc_kpc/cm_kpc, h_err_corr_units,'b', 0.2)
except NameError:
    pass

if g=='m31':
    ax.xaxis.set_ticks(np.arange(6, 18, 1)) #for m31
elif g=='m33':
    ax.xaxis.set_ticks(np.arange(0, 10, 1)) #for m33
elif g=='m51':
    ax.xaxis.set_ticks(np.arange(0, 11, 1)) #for m51
else:
    ax.xaxis.set_ticks(np.arange(0, 21, 2)) #for ngc


if g=='ngc6946':
    ax.yaxis.set_ticks(np.arange(0,max(h_err_corr_units+h)+400,200))
else:
    ax.yaxis.set_ticks(np.arange(0,max(h_err_corr_units+h)+400,200))

ax.set_xlabel(r'Radius (kpc)', fontsize=fs)
ax.set_ylabel(r'Length scales (pc)', fontsize=fs)
# ax.axhline(y=0, color='black', linestyle=':', alpha = 1)

# if g=='m33':
#     ax.legend(fontsize=lfs, frameon=True, handlelength=hd, ncol=1, bbox_to_anchor=(0.6, 1),prop={
#             'size': leg_textsize, 'family': 'Times New Roman'}, fancybox=True, framealpha=0.9, handletextpad=0.7, columnspacing=0.7)
if g=='m31':
    ax.legend(fontsize=lfs, frameon=False, handlelength=hd, ncol=1, bbox_to_anchor=(0.6, 1), prop={
            'size': leg_textsize,  'family': 'Times New Roman'}, fancybox=True, framealpha=0.9, handletextpad=legend_labelspace, columnspacing=0.7)
if g=='m51':
    ax.legend(fontsize=lfs, frameon=False, handlelength=hd, ncol=1, bbox_to_anchor=(0.65, 1), prop={
            'size': leg_textsize,  'family': 'Times New Roman'}, fancybox=True, framealpha=0.9, handletextpad=legend_labelspace, columnspacing=0.7)
else:
    ax.legend(fontsize=lfs, frameon=False, handlelength=hd, ncol=1, bbox_to_anchor=(0.7, 1),prop={
            'size': leg_textsize, 'family': 'Times New Roman'}, fancybox=True, framealpha=0.9,handletextpad=legend_labelspace, columnspacing=0.7)
# else:
#     ax.legend(fontsize=lfs, frameon=True, handlelength=hd, ncol=1, bbox_to_anchor=(0.6, 1),prop={
#             'size': leg_textsize, 'family': 'Times New Roman'}, fancybox=True, framealpha=0.9, handletextpad=0.7, columnspacing=0.7)

plt.savefig(save_files_dir+r'\1 h,l')

fig,ax = plt.subplots(nrows=1, ncols=1, figsize=(7, 5), tight_layout=True)

##################################################################################################################
u=u_f/cm_km
cs=cs_f /cm_km
sig = (np.sqrt(u_f**2 + (cs_f)**2))/cm_km
dat_u=dat_u/cm_km

if g=='m31':
    ref_cs= ' (T13)'
    ref_data= ' (no warp)'
elif g=='m33':
    ref_cs= ' (L17)'
    ref_data= ' (K17)'
elif g=='m51':
    ref_cs= ' (B04)'
    ref_data= ' (S07)'
else:
    ref_cs= ' (G13)'
    ref_data= ' (B08)'

# ax.plot(kpc_r, u, color='tab:orange', marker='o', mfc='k',linestyle=':',mec='k', markersize=m, label=r'$u$')
ax.plot(kpc_r, u, color='tab:orange',linestyle='-.',linewidth=lw,label=r'$u$')
try:
    fill_error(ax, u_f/cm_km,u_err/cm_km, 'tab:orange', 0.5)
except NameError:
    pass
percent_err_u= u_err/u_f
print('percent error u',percent_err_u)

ax.plot(kpc_r, sig, color='r', linewidth=lw, label=r'$\sqrt{u^2+c_\mathrm{s}^2}$',linestyle='-')
sig_err=np.sqrt((u_f*u_err)**2 + (cs_f*cs_err)**2)/(sig*cm_km**2)
try:
    fill_error(ax, sig,sig_err, 'r')
except NameError:
    pass

#B04 for M51 (Bresolin), G13 for ngc (Gusev), L17 for M33 (Lin), T13 for M31 (tabatabei)
# ax.plot(kpc_r, cs, color='g', linewidth=lw, linestyle='--', label=r'$c_\mathrm{s}$ (G13)', alpha = 0.5)
ax.plot(kpc_r, cs, color='g', linewidth=lw, linestyle='--', label=r'$c_\mathrm{s}$'+ '{}'.format(ref_cs), alpha = 0.5)
try:
    fill_error(ax, cs_f/cm_km,cs_err/cm_km, 'green')
except NameError:
    pass
percent_err_cs= cs_err/cs_f
print('percent error cs',percent_err_cs)

percent_err_tot= sig_err/sig
print('percent error sig',percent_err_tot)
# next line is for warped u data for M31
if g=='m31':
    dat_u_warp=dat_u_warp/cm_km
    ax.plot(kpc_r, dat_u_warp, zorder=2,
              c='tab:cyan',  linestyle=' ', label='$\sigma_\mathrm{g}$ (with warp)', alpha = 1,marker='D',mfc='b',mec='k',mew=1, markersize = 8.75)

#data S07-M51 B08-ngc K17-M33
ax.plot(kpc_r, dat_u, zorder=2,c='r', label='$\sigma_\mathrm{g}$' + '{}'.format(ref_data),linestyle=' ',alpha = 1,marker='*',mfc='r',mec='k',mew=1, markersize = 13)
# ax.axhline(y=0, color='black', linestyle=':', alpha = 1) 


axis_pars(ax)
if g=='m31':
    ax.yaxis.set_ticks(np.arange(4,max(dat_u)+18,4))
    ax.legend(fontsize=lfs, frameon=False, handlelength=hd, ncol=2, bbox_to_anchor=(1, 1),prop={
            'size': leg_textsize, 'family': 'Times New Roman'}, fancybox=True, framealpha=0.9, handletextpad=legend_labelspace, columnspacing=0.7)
elif g=='m33':
    ax.yaxis.set_ticks(np.arange(2,max(dat_u)+8,2))
    ax.legend(fontsize=lfs, frameon=False, handlelength=hd, ncol=2, bbox_to_anchor=(1, 1),prop={
            'size': leg_textsize, 'family': 'Times New Roman'}, fancybox=True, framealpha=0.9, handletextpad=legend_labelspace, columnspacing=0.7)
elif g=='m51':
    ax.yaxis.set_ticks(np.arange(0,max(dat_u)+3,2))
    ax.legend(fontsize=lfs, frameon=False, handlelength=hd, ncol=2, bbox_to_anchor=(1, 0.3),prop={
            'size': leg_textsize, 'family': 'Times New Roman'}, fancybox=True, framealpha=0.9, handletextpad=legend_labelspace, columnspacing=0.7)
else:
    ax.yaxis.set_ticks(np.arange(4,max(sig)+8,4))
    ax.legend(fontsize=lfs, frameon=False, handlelength=hd, ncol=2, bbox_to_anchor=(1, 1),prop={
            'size': leg_textsize, 'family': 'Times New Roman'}, fancybox=True, framealpha=0.9, handletextpad=legend_labelspace, columnspacing=0.7)

ax.set_xlabel(r'Radius (kpc)', fontsize=fs)
ax.set_ylabel(r'Speed (km s$^{-1}$)',  fontsize=fs)

plt.savefig(save_files_dir+r'\2 speeds')

##################################################################################################################

Btot = G_scal_Bbartot*1e+6
Breg = G_scal_Bbarreg*1e+6
Bord = G_scal_Bbarord*1e+6

fig,ax = plt.subplots(nrows=1, ncols=1, figsize=(7, 5), tight_layout=True)

ax.plot(kpc_r, Btot , c='b', linewidth=lw, label=r' $B_{\mathrm{tot}}$')
ax.plot(kpc_r, Breg , c='r', linewidth=lw, label=r' $B_{\mathrm{reg}}$')
ax.plot(kpc_r, Bord , c='g', linewidth=lw, label=r' $B_{\mathrm{ord}}$')

try:
    fill_error(ax, G_scal_Bbartot*1e+6,G_scal_Bbartot_err*1e+6, 'b', 0.2)
    fill_error(ax, G_scal_Bbarreg*1e+6,G_scal_Bbarreg_err*1e+6, 'r', 0.2)
    fill_error(ax, G_scal_Bbarord*1e+6,G_scal_Bbarord_err*1e+6, 'g', 0.2)
except NameError:
    pass

#plotting
ax.plot(mrange, G_dat_Btot, zorder=2,c='b', linestyle=' ',  marker='*',mfc='b',mec='k',mew=1,markersize = 13,label=r' $B_{\mathrm{tot}}$ (B19)')#, label='Average Binned data $B_{tot}$ ($\mu G$)')
ax.plot(mrange, G_dat_Breg, zorder=2,c='r', linestyle=' ', marker='s',mfc='r',mec='k',mew=1,markersize = 8.75,label=r' $B_{\mathrm{reg}}$ (B19)')#, label='Average Binned data $B_{reg}$ ($\mu G$)')
ax.plot(mrange, G_dat_Bord, zorder=2,c='g', linestyle=' ', marker='D',mfc='g',mec='k',mew=1,markersize = 8.75,label=r' $B_{\mathrm{ord}}$ (B19)')#, label='Average Binned data $B_{ord}$ ($\mu G$)')

ax.set_xlabel(r'Radius (kpc)', fontsize=fs)
ax.yaxis.set_ticks(np.arange(0,max(Btot)+max(G_scal_Bbartot_err*1e+6)+10,5))

ax.xaxis.set_major_formatter(FormatStrFormatter('%g'))
ax.set_ylabel('Magnetic field strength ($\mathrm{\mu G}$)', fontsize=fs)
axis_pars(ax)
ax.legend(fontsize=lfs, frameon=False, handlelength=hd, ncol=2, bbox_to_anchor=(1,1),prop={
            'size': leg_textsize, 'family': 'Times New Roman'}, fancybox=True, framealpha=0.9, handletextpad=legend_labelspace, columnspacing=0.7)


plt.savefig(save_files_dir+r'\3 B')

fig,ax = plt.subplots(nrows=1, ncols=1, figsize=(7, 5), tight_layout=True)

#plotting

#markers: sq/circle for pB/po, solid lines for model
# model
ax.plot(kpc_r, 180*pB/np.pi , c='r', linestyle='-', linewidth=lw, label=r' $p_{\mathrm{B}}$')
try:
    fill_error(ax, 180*pB/np.pi,180*pB_err/np.pi, 'r')
except NameError:
    pass

# ax.plot(kpc_r, pb,c = 'b',linestyle='-', label = r' $p_{b}$')
# fill_error(ax,pb,dpb_dr, 'b')

ax.plot(kpc_r,180*po/np.pi , c='g', linestyle='-', linewidth=lw,label=r' $p_{\mathrm{ord}}$')
try:
    fill_error(ax, 180*po/np.pi,180*po_err/np.pi, 'g')
except NameError:
    pass

#plotting observations

if g=='m31':
###################################     M31    ######################################################################
    ax.errorbar(mrange, 180*M_pb_beck19/np.pi,zorder=2,elinewidth=1, yerr=180*err_M_pb_beck19/np.pi,ecolor='k', mew=1, capsize=2,
                c='k', linestyle=' ', mfc='orange', mec='k',label=r' $p_{\mathrm{B}}$ (B19)',barsabove=True,marker='D',markersize=11)
    ax.errorbar(range_po_beck19, 180*RM_pb_beck19/np.pi,zorder=2,elinewidth=1, yerr=180*err_RM_pb_beck19/np.pi,ecolor='k',  mew=1, capsize=2,
                c='k', linestyle=' ', mfc='orange', mec='k',label=r' $p_{\mathrm{B}}$ (B19)',barsabove=True,marker='^',markersize=11)
    ax.errorbar(range_po_beck19, 180*po_beck19/np.pi,zorder=2,elinewidth=1, yerr=180*err_po_beck19/np.pi, mew=1, capsize=2,
                c='k', linestyle=' ', marker='P', mfc='g', mec='k',label=r'$p_{\mathrm{ord}}$ (B19)',ecolor='k',markersize=11)
###################################     M31     ######################################################################

elif g=='m33':
###################################     M33     ######################################################################
    ax.errorbar(mrange, 180*pb_beck19/np.pi,zorder=2,elinewidth=1, yerr=180*err_pb_beck19/np.pi,ecolor='k', ms=11, mew=1, capsize=2,
                  c='k', linestyle=' ', mfc='orange', mec='k',marker='D',label=r' $p_{\mathrm{B}}$ (B19)',barsabove=True)
    ax.errorbar(po_range_beck19, 180*po_beck19/np.pi,zorder=2,elinewidth=1, yerr=180*err_po_beck19/np.pi, ms=11, mew=1, capsize=2,
                  c='k', linestyle=' ', marker='^', mfc='g', mec='k',label=r' $p_{\mathrm{ord}}$ (B19)',ecolor='k')
###################################     M33     ######################################################################

elif g=='m51':
###################################     M51     ######################################################################
# #Beck+19 data
    ax.errorbar(range_po_beck19, po_beck19/(np.pi/180),zorder=2,elinewidth=1, yerr=err_po/(np.pi/180), markersize=11, mew=1, capsize=2,
                  c='k', linestyle=' ', marker='o', mfc='y', mec='k',ecolor='k', label=r'$p_{\mathrm{ord}}$ (B19)')

#Borlaff+23 data
    ax.plot(range_po_borlaff23,po_borlaff23, zorder=2,c='k', linestyle=' ', marker='*',
              markersize=13, mfc='b', mec='k', label=r' $p_{\mathrm{ord}}$ (B23)')

#Borlaff+21 data
    ax.plot(range_po_borlaff21,po_borlaff21,zorder=2, c='k', linestyle=' ',
              markersize=m, mfc='orange', mec='k', label=r' $p_{\mathrm{ord}}$ (B21)', marker='D')

#Surgent+23 data
    ax.plot(range_po_surgent23,po_surgent23,zorder=2, c='k', linestyle=' ', marker='P',
              markersize=11, mfc='pink', mec='k', label=r' $p_{\mathrm{ord}}$ (S23)')

# #Beck+19 data p_B
    ax.errorbar(range_pb, dat_pb/(np.pi/180),zorder=2,elinewidth=1, yerr=180*err_pb/np.pi,ecolor='k', markersize=11, mew=1, capsize=2,
                  c='k', linestyle=' ', mfc='r', mec='k',label=r' $p_{\mathrm{B}}$ (B19)',barsabove=True,marker='v')

###################################     M51     ######################################################################

else:
###################################     NGC 6946     ######################################################################
    print('no pB data for 6946 currently')

#Borlaff+23
    ax.plot(range_po_borlaff23, po_borlaff23, zorder=2, marker='s',linestyle=' ',
              markersize=m, mfc='pink', mec='k', label=r' $p_{\mathrm{ord}}$ (B23)')

#Surgent+23
    ax.plot(range_po_Surgent23, po_Surgent23, zorder=2, marker='o',linestyle=' ',
              markersize=9, mfc='orange', mec='k', label=r' $p_{\mathrm{ord}}$ (S23)')

#Beck+19
    ax.errorbar(range1_beck19, 180*po_range1_beck19/np.pi, yerr=180*err_po_range1_beck19/np.pi, zorder=2,ms=9, mew=1, capsize=2,linestyle=' ',
                   marker='v', mfc='yellow',label=r' $p_{\mathrm{ord}}$ (B19 range 1)', mec='k',ecolor='k')#, label=r'Data $p_{o}$ (ordered field)(RM)')
    ax.errorbar(range2_beck19, 180*po_range2_beck19/np.pi, yerr=180*err_po_range2_beck19/np.pi, zorder=2,ms=9, mew=1, capsize=2,linestyle=' ',
                   marker='^', mfc='yellow',label=r' $p_{\mathrm{ord}}$ (B19 range 2)', mec='k',ecolor='k')#, label=r'Data $p_{o}$ (ordered field)(RM)')
###################################     NGC 6946     ######################################################################

axis_pars(ax)
ax.set_xlabel(r'Radius (kpc)', fontsize=fs)
ax.set_ylabel(r'Pitch angle (degrees)', fontsize=fs)

ax.yaxis.set_ticks(np.arange(0,100,10))


if g=='m51':
    ax.legend(fontsize=lfs, frameon=False, handlelength=hd, ncol=3, bbox_to_anchor=(1,1),prop={
            'size': leg_textsize-2, 'family': 'Times New Roman'}, fancybox=True, framealpha=0.9, handletextpad=legend_labelspace, columnspacing=0.7)
elif g=='ngc6946':
    ax.legend(fontsize=lfs, frameon=False, handlelength=hd, ncol=3, bbox_to_anchor=(1,1.05),prop={
            'size': leg_textsize-2, 'family': 'Times New Roman'}, fancybox=True, framealpha=0.9, handletextpad=legend_labelspace, columnspacing=0.7)
else:
    ax.legend(fontsize=lfs, frameon=False, handlelength=hd, ncol=2, bbox_to_anchor=(1,1),prop={
            'size': leg_textsize, 'family': 'Times New Roman'}, fancybox=True, framealpha=0.9, handletextpad=legend_labelspace, columnspacing=0.7)

plt.savefig(save_files_dir+r'\4 pitch angles')
##################################################################################################################

# alpha plotting
fig,ax = plt.subplots(nrows=1, ncols=1, figsize=(7, 5), tight_layout=True)

ax.plot(kpc_r, alphak_f/cm_km, linewidth=lw, color='b', label=r'$\mathrm{\alpha_k}$')
# try:
#   fill_error(ax, alphak_f/cm_km,alphak_err/cm_km, 'blue')
# except NameError:
#   pass
    
alpham_f = alphak_f*((1/dkdc_f)-1)
alphasat_f = alphak_f + alpham_f

a=alpham_f/cm_km
b=alphak_f/cm_km


ax.plot(kpc_r, alpham_f/cm_km, linewidth=lw, color='r', label=r'$\mathrm{\alpha_m}$')
ax.plot(kpc_r, alphasat_f/cm_km, linewidth=lw, color='m',label=r'$\mathrm{\alpha_{\mathrm{sat}}}$')

# if g=='m31':
#     ax.yaxis.set_ticks(np.arange(-(round(min(abs(a)),1)+2.2),round(max(b),2)+0.4,0.8))
# elif g=='m33':
#     ax.yaxis.set_ticks(np.arange(-0.5,0.7,0.2))
# elif g=='m51':
#     ax.yaxis.set_ticks(np.arange(-(round(min(abs(a)),1)+2.5),round(max(b),1)+0.4,0.8))
# else:
#     ax.yaxis.set_ticks(np.arange(-(round(min(abs(a)),1)+2.1),round(max(b),1)+0.8,0.8))

ax.axhline(y=0, color='black', linestyle=':', alpha = 1)


axis_pars(ax)
ax.set_xlabel(r'Radius (kpc)', fontsize=fs)
ax.set_ylabel(r'$\alpha$ (km s$^{-1}$)', fontsize=fs)
ax.legend(fontsize=lfs, frameon=False, handlelength=hd, ncol=3, bbox_to_anchor=(1,1),prop={
            'size': leg_textsize+5, 'family': 'Times New Roman'}, fancybox=True, framealpha=0.9, handletextpad=legend_labelspace, columnspacing=0.7)

plt.savefig(save_files_dir+r'\5 alphas')

##################################################################################################################
fig,ax = plt.subplots(nrows=1, ncols=1, figsize=(7, 5), tight_layout=True)


omega = Symbol('\Omega')
kalpha = Symbol('K_alpha')
calpha = Symbol('C_alpha')
tau_f = taue_f
omt = datamaker(omega, data_pass, h_f, tau_f)*tau_f
kah = datamaker(kalpha/calpha, data_pass, h_f, tau_f)*(h_f/(tau_f*u_f))
ax.axhline(y=0, color='k', linestyle=':')

ax.axhline(y=1, color='g', linewidth=lw,linestyle=':',alpha=1)
ax.plot(kpc_r, omt, c='tab:orange', linewidth=lw,linestyle='-.', label=r'$\Omega\tau$')
ax.plot(kpc_r, kah, c='tab:blue',linewidth=lw,linestyle='-',  label=r'$\frac{K_\mathrm{\alpha} h}{C_\mathrm{\alpha} \tau u}$')

if g=='ngc6946':
    ax.yaxis.set_ticks(np.arange(0,max(kah)+2,1))
else:
    ax.yaxis.set_ticks(np.arange(0,max(kah)+1,1))

if galaxy_name=='m31':
    ax.legend(fontsize=lfs, frameon=False, handlelength=hd, ncol=1, bbox_to_anchor=(0.3,1),prop={
            'size': leg_textsize+5, 'family': 'Times New Roman'}, fancybox=True, framealpha=0.9, handletextpad=legend_labelspace, columnspacing=0.7)
elif galaxy_name=='m33':
    ax.legend(fontsize=lfs, frameon=False, handlelength=hd, ncol=1, bbox_to_anchor=(1,0.6),prop={
            'size': leg_textsize+5, 'family': 'Times New Roman'}, fancybox=True, framealpha=0.9, handletextpad=legend_labelspace, columnspacing=0.7)
elif galaxy_name=='m51':
    ax.legend(fontsize=lfs, frameon=False, handlelength=hd, ncol=1, bbox_to_anchor=(1,0.6),prop={
            'size': leg_textsize+5, 'family': 'Times New Roman'}, fancybox=True, framealpha=0.9, handletextpad=legend_labelspace, columnspacing=0.7)
else:
    ax.legend(fontsize=lfs, frameon=False, handlelength=hd, ncol=1, bbox_to_anchor=(1,0.6),prop={
            'size': leg_textsize+5, 'family': 'Times New Roman'}, fancybox=True, framealpha=0.9, handletextpad=legend_labelspace, columnspacing=0.7)
axis_pars(ax)
ax.set_xlabel('Radius (kpc)', fontsize=fs)
ax.set_ylabel(r'Condition for $\mathrm{\alpha_k}$ (Myr)', fontsize=fs)
plt.savefig(save_files_dir+r'\6 omega')

##################################################################################################################

fig,ax = plt.subplots(nrows=1, ncols=1, figsize=(7, 5), tight_layout=True)
ax.set_yscale('log')
ax.plot(kpc_r, taue_f/s_Myr,linewidth=lw, c='b',linestyle='-', label=r'$\tau^\mathrm{e}$')
ax.plot(kpc_r, taur_f/s_Myr, linewidth=lw,c='g',linestyle='--',label=r'$\tau^\mathrm{r}$')
ax.plot(kpc_r, h_f/(u_f*s_Myr), linewidth=lw,c='y',linestyle='-.', label=r'$h/u$')
ax.plot(kpc_r, 1/omt, c='r',linewidth=lw,linestyle=':', label=r'$\Omega^{-1}$')

t=h_f/(u_f*s_Myr)
# ax.yaxis.set_ticks(np.arange(1,round(max(taur_f/s_Myr),1)+0.4,0.8))
if g=='m33':
    plt.ylim(1, 500)
elif g=='m31':
    plt.ylim(2,400)
elif g=='m51':
    plt.ylim(0.1,100)
else:
    plt.ylim(0.08, 90)

if g=='m31':
    ax.legend(fontsize=lfs, frameon=False, handlelength=hd, ncol=2, bbox_to_anchor=(0.6,1),prop={
            'size': leg_textsize+5, 'family': 'Times New Roman'}, fancybox=True, framealpha=0.9, handletextpad=legend_labelspace, columnspacing=0.7)
elif g=='m33':
    ax.legend(fontsize=lfs, frameon=False, handlelength=hd, ncol=2, bbox_to_anchor=(0.5,1),prop={
            'size': leg_textsize+5, 'family': 'Times New Roman'}, fancybox=True, framealpha=0.9, handletextpad=legend_labelspace, columnspacing=0.7)
elif g=='m51':
    ax.legend(fontsize=lfs, frameon=False, handlelength=hd, ncol=2, bbox_to_anchor=(0.5,1.05),prop={
            'size': leg_textsize+5, 'family': 'Times New Roman'}, fancybox=True, framealpha=0.9, handletextpad=legend_labelspace, columnspacing=0.7)
else:
    ax.legend(fontsize=lfs, frameon=False, handlelength=hd, ncol=2, bbox_to_anchor=(1,0.3),prop={
            'size': leg_textsize+5, 'family': 'Times New Roman'}, fancybox=True, framealpha=0.9, handletextpad=legend_labelspace, columnspacing=0.7)
axis_pars(ax)
ax.set_xlabel('Radius (kpc)', fontsize=fs)
ax.set_ylabel(r'Correlation time $\tau$ (Myr)', fontsize=fs)
ax.axhline(y=0, color='k', linestyle=':')

plt.savefig(save_files_dir+r'\7 times')

##################################################################################################################

#plotting dynamo number
fig,ax = plt.subplots(nrows=1, ncols=1, figsize=(7, 5), tight_layout=True)
# ax.xaxis.set_ticks(np.arange(0, 10, 1)) #for m51

gamma= ((((np.pi**2)*(tau_f*(u_f**2))/3*(np.sqrt(dkdc_f)-1)/(4*h_f**2))**(-1))/
              (s_Myr*1e+3))**(-1)
p1=ax.plot(kpc_r,gamma , c='b', linewidth=lw,linestyle='-', label=r'$\gamma$')

ax.set_ylabel(r'Local growth rate $\gamma$ ($\mathrm{Gyr^{-1}}$)', fontsize=fs)

ax2 = ax.twinx()
p2=ax2.plot(kpc_r, dkdc_f,c='r',linestyle='-',linewidth=lw,label=r'$D_\mathrm{k}/D_\mathrm{c}$')
# ax2.plot(kpc_r, 1*np.ones(len(kpc_r)),linestyle=':',)
axis_pars(ax)
ax2.set_ylabel(r'$D_\mathrm{k}/D_\mathrm{c}$',  fontsize=fs)
ax.set_xlabel('Radius (kpc)', fontsize=fs)
ax2.axhline(y=1, color='black', linestyle=':', alpha = 1)

# try:
#     fill_error(ax2, dkdc_f,dkdc_err, 'r', 0.5)
# except NameError:
#     pass

if g=='m31':
    ax.xaxis.set_ticks(np.arange(6, 18, 1)) #for m31
elif g=='m33':
    ax.xaxis.set_ticks(np.arange(0, 10, 1)) #for m33
elif g=='m51':
    ax.xaxis.set_ticks(np.arange(0, 11, 1)) #for m51
else:
    ax.xaxis.set_ticks(np.arange(0, 21, 2)) #for ngc

if galaxy_name=='m33':
    ax.yaxis.set_ticks(np.arange(0,max(gamma)+1,1))
else:
    ax.yaxis.set_ticks(np.arange(0,max(gamma)+2,2))
    
ax2.yaxis.set_ticks(np.arange(0,max(dkdc_f)+6,2))
# ax2.yaxis.set_ticks(np.arange(0,max(dkdc_f+dkdc_err)+6,10))


leg = p1 + p2
labs = [l.get_label() for l in leg]
ax.legend(leg, labs,fontsize=lfs, frameon=False, handlelength=hd, ncol=1, bbox_to_anchor=(0.8,1),prop={
            'size': leg_textsize+5, 'family': 'Times New Roman'}, fancybox=True, framealpha=0.9, handletextpad=legend_labelspace, columnspacing=0.7)

plt.savefig(save_files_dir+r'\8 dynamo number_efolding')


##################################################################################################################

# for plotting in pdf
from datetime import date
os.chdir(os.path.join(base_path,'plots'))

from matplotlib.backends.backend_pdf import PdfPages
PDF = PdfPages(f'{galaxy_name}_ca_'+str(params[r'C_\alpha'])+'K_'+str(params[r'K'])+'z_'+str(
      params[r'\zeta'])+'psi_'+str(params[r'\psi'])+'b_'+str(params[r'\beta'])+'.pdf')#('plots_model'+str(model_no)+let+'t_vary_'+'ca_'+str(ca)+'rk_'+str(rk)+'z_'+str(z.mean())+'.pdf')

fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(7, 10), tight_layout=True)
fig.suptitle(r'$C_\alpha$ = '+str(params[r'C_\alpha'])+r'    $K$ = '+str(params[r'K'])+
             r'    $\zeta$ = '+str(params[r'\zeta'])+r'    $\psi$ = '+str(params[r'\psi'])+r'    $\beta$ = '+str(params[r'\beta']), weight = 15)


#scale height and correlation length
i = 0


ax[i].plot(kpc_r, h_f*pc_kpc/cm_kpc, c='r', linestyle='-', linewidth=lw, label=r'$h$ (pc)')
try:
    ax[i].plot(kpc_dat_r, pc_dat_h, linestyle='--',marker='*',mfc='b',mec='k',mew=1, markersize = 9, label=r'$h$ (pc) (Chamandy+16)')
except NameError:
    pass
ax[i].plot(kpc_r, l_f*pc_kpc/cm_kpc, c='g',linestyle='-.', linewidth=lw, label=r'$l$ (pc)')
# ax[i].plot(kpc_r, datamaker(lsn , data_pass, h_f, tau_f)*pc_kpc/cm_kpc,c = 'y',linestyle='--',mfc='k',mec='k', marker='o')
ax[i].axhline(y=100, color='black', linestyle='--', alpha = 0.2)
#ax[i].set_yticks(list(plt.yticks()[0])+[100])
axis_pars(ax[i])
try:
    fill_error(ax[i], h_f*pc_kpc/cm_kpc, h_err*pc_kpc/cm_kpc)
except NameError:
    pass
    
ax[i].set_xlabel(r'Radius (kpc)', fontsize=fs)
ax[i].set_ylabel(r'Length scale (pc)', fontsize=fs)
ax[i].legend(fontsize=lfs, frameon=True, handlelength=hd, ncol=1, bbox_to_anchor=(0.45, 1),prop={
            'size': leg_textsize, 'family': 'Times New Roman'}, fancybox=True, framealpha=0.9, handletextpad=0.7, columnspacing=0.7)


##################################################################################################################

#velociy plots
i = 1

u=u_f/cm_km
cs=cs_f /cm_km
sig = (np.sqrt(u_f**2 + (cs_f)**2))/cm_km

ax[i].plot(kpc_r, u, color='tab:orange',linestyle='-.',linewidth=lw,label=r'$u$')
try:
    fill_error(ax[i], u,u_err/cm_km, 'tab:orange', 0.5)
except NameError:
    pass
    

# ax[i].plot(kpc_r, alphak_f/cm_km, color='b', label=r'$\alpha_k$')
# try:
#   fill_error(ax[i], alphak_f/cm_km,alphak_err/cm_km, 'blue')
# except NameError:
#   pass
    

# ax[i].plot(kpc_r, alpham_f/cm_km, color='r', label=r'$\alpha_m$')
# ax[i].plot(kpc_r, alphasat_f/cm_km, color='m',label=r'$\alpha_{sat}$')

ax[i].plot(kpc_r, sig , color='r', linewidth=lw, label=r'$\sqrt{u^2+c_s^2}$',linestyle='-')
try:
    fill_error(ax[i], sig,np.sqrt((u*u_err/cm_km)**2 + (cs*cs_err/cm_km)**2)/(sig), 'r')
except NameError:
    pass
    

ax[i].plot(kpc_r, cs, color='g', linewidth=lw, linestyle='--', label=r'$c_s$', alpha = 0.5)
try:
    fill_error(ax[i], cs,cs_err/cm_km, 'green')
except NameError:
    pass
ax[i].axhline(y=0, color='black', linestyle=':', alpha = 0.2)

#################################################################################################################################################################
#velocity dispersion of each galaxy  
try:
    ax[i].plot(kpc_r, dat_u/cm_km, c='k', label='$\sigma_\star$',alpha = 1,marker='*',mfc='y',mec='k',mew=1, markersize = 9)
    ax[i].plot(kpc_r, dat_u_warp/cm_km, c='k',linestyle='dashdot', label='$\sigma_\star$ with warp', alpha = 1,marker='D',mfc='tab:cyan',mec='k',mew=1, markersize = m)
except NameError:
    pass
#################################################################################################################################################################

axis_pars(ax[i])

ax[i].set_ylim(0)
ax[i].set_xlabel(r'Radius ($kpc$)', fontsize=fs)
ax[i].set_ylabel(r'Speed (km/s)',  fontsize=fs)
ax[i].legend(fontsize=lfs, frameon=True, handlelength=hd, ncol=1, bbox_to_anchor=(1, 1),prop={
            'size': leg_textsize, 'family': 'Times New Roman'}, fancybox=True, framealpha=0.9, handletextpad=0.7, columnspacing=0.7)
# plt.savefig('D:\Documents\Gayathri_college\MSc project\codes\MSc.-Thesis\final_plots'+r'\3 speeds')
PDF.savefig(fig)
###############################################################################################################################################################

#pitch angle plots
fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(7, 10), tight_layout=True)

i = 1

#model plots
ax[i].axhline(y=0, color='k', linestyle=':')

ax[i].plot(kpc_r, 180*po/np.pi, c='g', linestyle='-.',linewidth=lw, label=r'$p_{o}$')
try:
    fill_error(ax[i], 180*po/np.pi,180*po_err/np.pi, 'g')
except NameError:
    pass

ax[i].plot(kpc_r, 180*pB/np.pi, c='r', linestyle='-', linewidth=lw, label=r'$p_{B}$')
try:
    fill_error(ax[i], 180*pB/np.pi,180*pB_err/np.pi, 'r')
except NameError:
    pass

#this part varies according to galaxy
# m31
# try:
#     ax[i].errorbar(range_po_beck19, RM_pb_beck19/(np.pi/180),elinewidth=1, yerr=err_RM_pb_beck19/(np.pi/180),ecolor='k', ms=m, mew=1, capsize=2,
#                     c='r', mfc='r', mec='k',barsabove=True,marker='D',label=r' $p_{B}$ (Beck+19, RM)')
#     ax[i].errorbar(range_po_beck19, po_beck19/(np.pi/180),elinewidth=1, yerr=err_po_beck19/(np.pi/180), ms=m, mew=1, capsize=2,
#                     c='g', marker='s', mfc='g', mec='k',barsabove=True,ecolor='k',label=r'$p_{o}$ (Beck+19)')
#     ax[i].errorbar(range_pb_beck19, M_pb_beck19/(np.pi/180),elinewidth=1, yerr=err_M_pb_beck19/(np.pi/180), ms=m, mew=1, capsize=2,
#                     c='r', marker='o', mfc='r', mec='g',barsabove=True,ecolor='k',label=r' $p_{B}$ (Beck+19, M)')
# except NameError:
#     pass

#m33
# try:
#     ax[i].errorbar(pb_range_beck19, pb_beck19/(np.pi/180),elinewidth=1, yerr=err_pb_beck19/(np.pi/180),ecolor='k', ms=m, mew=1, capsize=2,
#                     c='r', mfc='r', mec='k',barsabove=True,marker='D',label=r' $p_{B}$ (Beck+19)')
#     ax[i].errorbar(po_range_beck19, po_beck19/(np.pi/180),elinewidth=1, yerr=err_po_beck19/(np.pi/180), ms=m, mew=1, capsize=2,
#                     c='g', marker='s', mfc='g', mec='k',barsabove=True,ecolor='k',label=r'$p_{o}$ (Beck+19)')
# except NameError:
#     pass

#m51
try:
    ax[i].errorbar(range_po_beck19, po_beck19/(np.pi/180),elinewidth=1,linestyle=' ', yerr=err_po/(np.pi/180), ms=m, mew=1, capsize=2,
                     marker='s', mfc='g', mec='k',barsabove=True,ecolor='k',label=r' $p_{o}$ Beck+19 (i=$20^\circ$)')
    ax[i].plot(range_po_surgent23,po_surgent23, linestyle=' ', marker='D',markersize=m, mfc='orange', mec='k', label=r' $p_{o}$ Surgent+23 (i=$22.5^\circ$)')
    ax[i].plot(range_po_borlaff23,po_borlaff23, linestyle=' ', marker='*',
              markersize=9, mfc='pink', mec='k', label=r' $p_{o}$ Borlaff+23 (i=$22^\circ$)')
    ax[i].plot(range_po_borlaff21,po_borlaff21, linestyle=' ', marker='P',
              markersize=m, mfc='b', mec='k', label=r' $p_{o}$ Borlaff+21 (i=$22.5^\circ$)')
    ax[i].errorbar(range_pb, dat_pb/(np.pi/180),elinewidth=1,linestyle=' ', yerr=err_pb/(np.pi/180), ms=m, mew=1, capsize=2,
                     marker='o', mfc='r', mec='k',barsabove=True,ecolor='k',label=r' $p_{B}$ Beck+19')
except NameError:
    pass

#NGC 6946
# try:
#     ax[i].errorbar(range1_beck19, po_range1_beck19/(np.pi/180),yerr=180*err_po_range1_beck19/np.pi,ms=9, mew=1, capsize=2,linestyle=' ',marker='s', mfc='g', mec='k',ecolor='k',label=r' $p_{o}$ (Beck+19 range1 i= $30^\circ$)')
#     ax[i].errorbar(range2_beck19, po_range2_beck19/(np.pi/180),yerr=180*err_po_range2_beck19/np.pi,ms=9, mew=1, capsize=2,linestyle=' ',marker='D', mfc='b',mec='k',ecolor='k',label=r' $p_{o}$ (Beck+19 range2 i= $30^\circ$)')
#     ax[i].plot(range_po_Surgent23,po_Surgent23,marker='*',linestyle=' ',markersize=9, mfc='orange', mec='k', label=r' $p_{o}$ (Surgent+23, i= $38.4^\circ)$')
#     ax[i].plot(range_po_borlaff23,po_borlaff23,marker='P',linestyle=' ',markersize=9, mfc='pink', mec='k', label=r' $p_{o}$ (Borlaff+23, i= $38.4^\circ)$')
# except NameError:
#     pass

axis_pars(ax[i])
ax[i].set_xlabel(r'Radius ($kpc$)', fontsize=fs)
ax[i].set_ylabel(r'Pitch angle (deg)', fontsize=fs)
ax[i].legend(fontsize=lfs, frameon=True, handlelength=hd, ncol=1, bbox_to_anchor=(1,1),prop={
            'size': 9.5, 'family': 'Times New Roman'}, fancybox=True, framealpha=0.9, handletextpad=0.7, columnspacing=0.7)

##################################################################################################################
#field strengths
i = 0

ax[i].plot(kpc_r, G_scal_Bbartot*1e+6, c='b', linestyle='-',label=r' $B_{tot}=\sqrt{\bar{B}^2+b_{iso}^2+b_{ani}^2}$')
ax[i].plot(kpc_r, G_scal_Bbarreg*1e+6, c='r', linestyle='-',label=r' $B_{reg} = \bar{B}$')
ax[i].plot(kpc_r, G_scal_Bbarord*1e+6, c='green', linestyle='-', label=r' $B_{ord} = \sqrt{\bar{B}^2+b_{ani}^2}$')
ax[i].axhline(y=0, color='k', linestyle=':')

try:
    fill_error(ax[i], G_scal_Bbartot*1e+6,G_scal_Bbartot_err*1e+6, 'b')
    fill_error(ax[i], G_scal_Bbarreg*1e+6,G_scal_Bbarreg_err*1e+6, 'r', 0.2)
    fill_error(ax[i], G_scal_Bbarord*1e+6,G_scal_Bbarord_err*1e+6, 'g', 0.2)
except NameError:
    pass

#this part changes according to galaxy
#M31
# try:
#     ax[i].plot(range_pb_beck19, G_dat_Btot, c='b', linestyle='-.',  marker='*',mfc='b',mec='k',mew=1,markersize = 9,label=r' $B_{tot}$ (Beck+19)')
#     ax[i].plot(range_pb_beck19, G_dat_Breg, c='r', linestyle='-.', marker='s',mfc='r',mec='k',mew=1,markersize = 7,label=r' $B_{reg}$ (Beck+19)')
#     ax[i].plot(range_pb_beck19, G_dat_Bord, c='g', linestyle='-.', marker='D',mfc='g',mec='k',mew=1,markersize = 7,label=r' $B_{ord}$ (Beck+19)')
# except NameError:
#     pass

# #M33
# try:
#     ax[i].plot(pb_range_beck19, G_dat_Btot, c='b', linestyle='-.',  marker='*',mfc='b',mec='k',mew=1,markersize = 9,label=r' $B_{tot}$ (Beck+19)')
#     ax[i].plot(pb_range_beck19, G_dat_Breg, c='r', linestyle='-.', marker='s',mfc='r',mec='k',mew=1,markersize = 7,label=r' $B_{reg}$ (Beck+19)')
#     ax[i].plot(pb_range_beck19, G_dat_Bord, c='g', linestyle='-.', marker='D',mfc='g',mec='k',mew=1,markersize = 7,label=r' $B_{ord}$ (Beck+19)')
# except NameError:
#     pass

#M51
try:
    ax[i].plot(range_pb, G_dat_Btot, c='b', linestyle='-.',  marker='*',mfc='b',mec='k',mew=1,markersize = 9,label=r' $B_{tot}$ (Beck+19)')
    ax[i].plot(range_pb, G_dat_Breg, c='r', linestyle='-.', marker='s',mfc='r',mec='k',mew=1,markersize = 7,label=r' $B_{reg}$ (Beck+19)')
    ax[i].plot(range_pb, G_dat_Bord, c='g', linestyle='-.', marker='D',mfc='g',mec='k',mew=1,markersize = 7,label=r' $B_{ord}$ (Beck+19)')
except NameError:
    pass

#NGC 6946
# try:
#     ax[i].plot(mrange, G_dat_Btot, c='b', linestyle='-.',  marker='*',mfc='b',mec='k',mew=1,markersize = 9,label=r' $B_{tot}$ (Beck+19)')
#     ax[i].plot(mrange, G_dat_Breg, c='r', linestyle='-.', marker='s',mfc='r',mec='k',mew=1,markersize = 7,label=r' $B_{reg}$ (Beck+19)')
#     ax[i].plot(mrange, G_dat_Bord, c='g', linestyle='-.', marker='D',mfc='g',mec='k',mew=1,markersize = 7,label=r' $B_{ord}$ (Beck+19)')
# except NameError:
#     pass

ax[i].xaxis.set_major_formatter(FormatStrFormatter('%g'))
ax[i].legend(fontsize=lfs, frameon=True, handlelength=hd, ncol=1, bbox_to_anchor=(1,1),prop={'size': 9.5, 'family': 'Times New Roman'}, fancybox=True, framealpha=0.9, handletextpad=0.7, columnspacing=0.7)

ax[i].set_xlabel(r'Radius ($kpc$)', fontsize=fs)
ax[i].set_ylabel('Magnetic field strength ($\mu G$)', fontsize=fs)


axis_pars(ax[i])

PDF.savefig(fig)
##################################################################################################################

#common for all galaxies
fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(10, 10), tight_layout=True)

i = 0
j = 0
axis_pars(ax[i][j])

omega = Symbol('\Omega')
kalpha = Symbol('K_alpha')
calpha = Symbol('C_alpha')
tau_f = taue_f
omt = datamaker(omega, data_pass, h_f, tau_f)*tau_f
kah = datamaker(kalpha/calpha, data_pass, h_f, tau_f)*(h_f/(tau_f*u_f))
ax[i][j].axhline(y=0, color='k', linestyle=':')

ax[i][j].axhline(y=1, color='g', linestyle=':',linewidth=lw, label=r'1')
ax[i][j].plot(kpc_r, omt, c='tab:orange', linestyle='-.',linewidth=lw, label=r'$\Omega\tau$')
ax[i][j].plot(kpc_r, kah, c='tab:blue',linestyle='-',linewidth=lw,  label=r'$\frac{K_\alpha h}{C_\alpha \tau u}$')
ax[i][j].legend(fontsize=lfs, frameon=True, handlelength=hd, ncol=1, bbox_to_anchor=(0.3,1),prop={
            'size': leg_textsize, 'family': 'Times New Roman'}, fancybox=True, framealpha=0.9, handletextpad=0.7, columnspacing=0.7)
ax[i][j].set_xlabel('Radius(kpc)', fontsize=fs)
ax[i][j].set_ylabel(r'Condition for $\alpha_k$ (Myr)', fontsize=fs)


j = 1
axis_pars(ax[i][j])

ax[i][j].plot(kpc_r, taue_f/s_Myr, c='b',linestyle='-.', label=r'$\tau^e$')
ax[i][j].plot(kpc_r, taur_f/s_Myr, c='g',linestyle='--',label=r'$\tau^r$')
ax[i][j].plot(kpc_r, h_f/(u_f*s_Myr), c='y',linestyle='-', label=r'$h/u$')
ax[i][j].legend(fontsize=lfs, frameon=True, handlelength=hd, ncol=1, bbox_to_anchor=(0.3,1),prop={
            'size': leg_textsize, 'family': 'Times New Roman'}, fancybox=True, framealpha=0.9, handletextpad=0.7, columnspacing=0.7)
ax[i][j].set_xlabel('Radius(kpc)', fontsize=fs)
ax[i][j].set_ylabel(r'Correlation Time $\tau$ (Myr)', fontsize=fs)
ax[i][j].axhline(y=0, color='k', linestyle=':')


#plotting dynamo number
i = 1
axis_pars(ax[i][j])

ax[i][j].plot(kpc_r, dkdc_f,c='r',linestyle='-',linewidth=lw,label=r'$D_k/D_c$')
ax[i][j].plot(kpc_r, 1*np.ones(len(kpc_r)),linestyle=':',)
ax[i][j].set_xlabel('Radius(kpc)', fontsize=fs)
ax[i][j].set_ylabel(r'$D_k/D_c$',  fontsize=fs)

try:
    fill_error(ax[i][j], dkdc_f,dkdc_err, 'r', 0.5)
except NameError:
    pass
ax[i][j].legend(fontsize=lfs, frameon=True, handlelength=hd, ncol=1, bbox_to_anchor=(1,1),prop={
            'size': leg_textsize, 'family': 'Times New Roman'}, fancybox=True, framealpha=0.9, handletextpad=0.7, columnspacing=0.7)
ax[i][j].axhline(y=1, color='k', linestyle=':')

#e-folding time
j = 0
axis_pars(ax[i][j])

ax[i][j].plot(kpc_r, (((np.pi**2)*(tau_f*(u_f**2))/3*(np.sqrt(dkdc_f)-1)/(4*h_f**2))**(-1))/
              (s_Myr*1e+3), c='g', linestyle='-', label=r'$\frac{1}{\gamma}$')
ax[i][j].set_xlabel('Radius(kpc)', fontsize=fs)
ax[i][j].set_ylabel(r'local e-folding time $\frac{1}{\gamma}$ (Gyr)', fontsize=fs)
# ax[i][j].legend(fontsize = lfs)
ax[i][j].legend(fontsize=lfs, frameon=True, handlelength=hd, ncol=1, bbox_to_anchor=(1,1),prop={
            'size': leg_textsize, 'family': 'Times New Roman'}, fancybox=True, framealpha=0.9, handletextpad=0.7, columnspacing=0.7)

PDF.savefig(fig)
PDF.close()
os.chdir(current_directory)
