import numpy as np
import matplotlib.pyplot as plt
from sympy import *
from fractions import Fraction
import pickle
import os
import sys
from matplotlib.ticker import FormatStrFormatter


from scipy.interpolate import griddata



kpc = 1e+3
kpcm = 3.086e+21
pcm = kpcm/1e+3
pc_kpc = 1e3		#number of pc in one kpc
cm_km = 1e5		#number of cm in one km
s_day = 24*3600		#number of seconds in one day
s_min = 60		#number of seconds in one hour
s_hr = 3600		#number of seconds in one hour
cm_Rsun = 6.957e10	#solar radius in cm
g_Msun = 1.989e33	#solar mass in g
cgs_G = 6.674e-8 
cms_c = 2.998e10
g_mH = 1.6736e-24
g_me = 9.10938e-28
cgs_h = 6.626e-27
deg_rad = 180e0/np.pi
arcmin_deg = 60e0
arcsec_deg = 3600e0

kpc_D_M33_Plot = 840e0 

kpc_D_M33_Tabatabaei = 840e0  #used by Tabatabaei+08 (magnetic field data)
kpc_D_M33_Kam = 840e0  #used by Kam+2015,2017 -- quote an uncertainty of +/- 40 kpc (HI data)
kpc_D_M33_Koch = 840e0  #used by Koch+2018a (HI data)
kpc_D_M33_Verley = 840e0  #used by Verley+09 (SFR data)
kpc_D_M33_Gratier = 840e0  #used by Gratier+10 (H2 data)
kpc_rmax_M33 = 20e0
kpc_xbnd_M33_Tabatabaei = np.array([1,3,5])	#radial bins for magnetic field data (Tabatabaei+08)
kpc_xmid_M33_Tabatabaei = np.array([2,4])	#midpoints of radial bins for magnetic field data (Tabatabaei+08)

kpc_xbnd_M33 = kpc_xbnd_M33_Tabatabaei * kpc_D_M33_Plot / kpc_D_M33_Tabatabaei
kpc_xmid_M33 = kpc_xmid_M33_Tabatabaei * kpc_D_M33_Plot / kpc_D_M33_Tabatabaei

arcmin_xbnd_M33 = kpc_xbnd_M33 / kpc_D_M33_Plot * deg_rad * arcmin_deg
arcmin_xmid_M33 = kpc_xmid_M33 / kpc_D_M33_Plot * deg_rad * arcmin_deg

print('arcmin_xbnd_M33',arcmin_xbnd_M33)
print('arcmin_xmid_M33',arcmin_xmid_M33)

length, breadth = [5, 2.5]

#Specifiy text size in legend
leg_textsize = 10
axis_textsize = 10
rc = {"font.family" : "serif", 
      "mathtext.fontset" : "stix"}
plt.rcParams.update(rc)
plt.rcParams["font.serif"] = ["Times New Roman"] + plt.rcParams["font.serif"]


### M33

##HI data from Koch+2018a (MNRAS 479, 2505) and Kam+2017 (ApJ 154:41)

# Radius for circular velocity data from Kam+17, Tab. 4 (Tilted-ring model)
arcmin_r_Kam = np.array([2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62, 64, 66, 68, 70, 72, 74, 76, 78, 80, 82, 84, 86, 88, 90, 92, 94, 96])

# Vcirc from Kam+17, Tab. 4 (Tilted-ring model)
kms_vcirc_Kam = np.array([42.0, 58.8, 69.4, 79.3, 86.7, 91.4, 94.2, 96.5, 99.8, 102.1, 103.6, 105.9, 105.7, 106.8, 107.3, 108.3, 109.7, 112.0, 116.1, 117.2, 116.5, 115.7, 117.4, 116.8, 115.7, 115.1, 117.1, 118.2, 118.4, 118.2, 117.5, 119.6, 118.6, 122.6, 124.1, 125.0, 125.5, 125.2, 122.0, 120.4, 114.0, 110.0, 98.7, 100.1, 104.3, 101.2, 123.5, 115.3])

# Vcirc_error from Kam+17, Tab. 4 (Tilted-ring model)
kms_vcirc_error_Kam = np.array([2.4, 1.5, 0.4, 4.0, 1.8, 3.1, 4.8, 5.5, 3.9, 1.7, 0.4, 0.7, 1.7, 2.2, 3.0, 4.0, 4.0, 4.8, 2.2, 2.5, 6.5, 8.1, 8.2, 8.9, 9.6, 7.7, 5.1, 3.2, 1.4, 1.8, 2.4, 0.8, 1.5, 0.5, 2.9, 2.2, 2.5, 8.1, 9.8, 8.5, 26.6, 34.6, 27.4, 33.4, 35.2, 27.4, 39.1, 26.7])

# Radius for circular velocity data from Koch+18a, Tab. C1 (DISKFIT) (correction in Koch+18b table 1 gives identical radius values -- ??)
arcsec_r_Koch = np.array([9, 27, 45, 63, 81, 99, 117, 135, 153, 171, 189, 207, 225, 243, 261, 279, 297, 315, 333, 351, 369, 387, 405, 423, 441, 459, 477, 495, 513, 531, 549, 567, 585, 603, 621, 639, 657, 675, 693, 711, 729, 747, 765, 783, 801, 819, 837, 855, 873, 891, 909, 927, 945, 963, 981, 999, 1017, 1035, 1053, 1071, 1089, 1107, 1125, 1143, 1161, 1179, 1197, 1215, 1233, 1251, 1269, 1287, 1305, 1323, 1341, 1359, 1377, 1395, 1413, 1431, 1449, 1467, 1485, 1503, 1521, 1539, 1557, 1575, 1593, 1611, 1629, 1647, 1665, 1683, 1701, 1719, 1737, 1755, 1773, 1791, 1809, 1827])


# Vcirc from Koch+18a, Tab. C1 (DISKFIT) -- identical to corrected version in Tab. 1 of Koch+18b  (??)
kms_vcirc_Koch = np.array([1.58, 11.88, 22.71, 30.53, 29.95, 34.20, 36.76, 44.61, 47.80, 51.67, 52.73, 53.24, 55.30, 57.58, 58.70, 59.47, 60.30, 64.32, 67.66, 69.07, 72.05, 75.89, 74.09, 73.44, 77.06, 77.06, 77.82, 79.01, 80.96, 80.34, 80.35, 82.90, 86.23, 86.92, 86.13, 87.43, 87.89, 87.82, 91.08, 90.37, 88.22, 90.67, 92.32, 93.31, 94.15, 93.40, 94.02, 94.91, 95.57, 95.15, 93.19, 94.52, 95.00, 95.72, 96.89, 98.58, 98.09, 99.87, 99.10, 99.01, 97.58, 98.60, 99.61, 99.61, 102.02, 101.84, 102.76, 102.66, 102.57, 104.08, 103.24, 103.17, 103.33, 103.41, 103.91, 102.91, 105.50, 104.05, 106.03, 104.89, 106.68, 105.53, 105.80, 104.03, 105.58, 105.69, 106.32, 105.97, 105.58, 106.18, 106.08, 106.24, 106.39, 106.47, 107.19, 105.88, 106.58, 106.09, 106.30, 107.23, 104.01, 104.50])

# Vcirc_error from Koch+18a, Tab. C1 (DISKFIT) -- identical to corrected version in Tab. 1 of Koch+18b (??)
kms_vcirc_error_Koch = np.array([6.24, 5.26, 5.25, 5.12, 4.47, 4.03, 3.99, 3.79, 3.49, 3.08, 2.96, 3.06, 3.22, 3.30, 3.42, 3.50, 3.59, 3.78, 3.75, 4.03, 4.49, 4.82, 4.91, 3.29, 2.76, 2.80, 2.60, 2.48, 2.59, 2.50, 2.65, 2.66, 2.66, 2.59, 2.63, 2.57, 2.63, 2.35, 2.44, 2.25, 2.09, 2.36, 2.48, 2.26, 2.15, 2.37, 2.44, 2.50, 2.48, 2.22, 2.46, 2.63, 2.77, 2.90, 2.81, 2.73, 2.61, 2.53, 2.41, 2.40, 2.21, 2.26, 2.22, 2.49, 2.40, 2.52, 2.38, 2.62, 2.68, 2.58, 2.39, 2.42, 2.50, 2.64, 2.56, 2.63, 2.47, 2.47, 2.47, 2.57, 2.48, 2.36, 2.33, 2.23, 2.23, 2.26, 2.17, 2.32, 2.20, 2.26, 2.28, 2.36, 2.37, 2.43, 2.44, 2.58, 2.52, 2.56, 2.35, 2.61, 2.82, 2.59])

## HI velocity dispersion

# HI velocity dispersion from Kam+17 Tab. 8
# Entries that are listed as 0._ in the original work are probably typos. They should read 10._ ; this has been corrected in the following data.
kms_sigmaLOS_Kam = np.array([10.7, 9.5, 9.3, 9.6, 10.0, 10.1, 9.1, 8.2, 8.0, 8.4, 7.9, 8.0, 7.6, 7.6, 8.0, 8.9, 9.8, 10.2, 10.3, 10.5, 9.5, 8.8, 9.1, 9.0, 8.5, 8.2, 8.3, 7.6, 7.7, 6.8, 7.4, 7.5, 6.7, 6.6, 6.4, 8.0, 8.2, 6.5, 6.5, 7.4, 6.5, 5.6, 6.1, 5.7, 5.7, 5.6, 6.4, 7.5])

# HI velocity dispersion from Koch+18a Fig 16 
# kms_sigmaLOS = np.array([])

#kms_sigmaLOS = np.array([])

## HI Surface density

# HI surface density of disc (from Kam+17 Tab. 8)
Msunpc2_SigmaHI_Kam = np.array([7.89, 7.73, 8.26, 8.52, 8.23, 7.66, 8.00, 8.61, 8.14, 8.19, 7.61, 7.57, 7.54, 6.50, 4.99, 3.50, 2.51, 1.75, 1.34, 1.00, 0.75, 0.67, 0.56, 0.63, 0.49, 0.36, 0.23, 0.18, 0.12, 0.08, 0.04, 0.03, 0.03, 0.00, 0.03, 0.07, 0.03, 0.04, 0.04, 0.04, 0.08, 0.09, -1, -1, 0.09, -1, -1, -1])

# HI surface density of disc (from Van Eck+2015, Tab. 2 -- also used in Chamandy+2016)
Msunpc2_SigmaHI_Vaneck = np.array([11.3, 9.43])

# HI surface density from Koch+18a Fig. 6
#Msunpc2_SigmaHI_Koch = np.array([])

# Total surface density of disc (from Kam+15/Kam+17): See script data/M33/Sigma_Star.py

## SFR surface density data from Verley+2009, A&A 493,453, Fig. 13 (SFR S.D. according to 5 different tracers: 24\mu m, H\alpha, FUV, FIR, TIR)

# Radius, Verley+09 assumes a distance of 840 kpc
kpc_r_SFR = np.array([ 0.3715, 0.6153, 0.8592, 1.103, 1.346, 1.590, 1.834, 2.089, 2.333, 2.577, 2.821, 3.065, 3.309, 3.552, 3.796, 4.040, 4.284, 4.528, 4.772, 5.027, 5.271, 5.515, 5.759, 6.002, 6.246, 6.490, 6.734, 6.978, 7.222, 7.465, 7.709, 7.953])

# Star formation rate surface density (extinction-corrected FUV -- Fig. 13a, use Dexter)
Msunpc2Gyr_Sigma_SFR = np.array([21.11, 15.96, 15.01, 13.25, 12.46, 13.27, 12.47, 10.19, 8.079, 6.500, 5.231, 5.315, 7.606, 5.574, 4.848, 4.283, 3.500, 2.997, 2.862, 2.232, 1.853, 2.132, 1.663, 1.338, 0.9077, 0.7191, 0.4801, 0.3205, 0.2106, 0.1721, 0.1281, 0.06569])

# Values used in Van Eck+15/Chamandy+16
Msunpc2Gyr_Sigma_SFR_Vaneck = np.array([9.64, 3.99])


## Molecular gas surface density (using Gratier+10, A&A 522,A3, Fig. 8 -- from Dexter)

# Radius for molecular gas surface density (assumes distance of 840 kpc, Gratier+10, Tab. 1)
kpc_r_SigmaH2 = np.array([0.2671, 0.7504, 1.257, 1.752, 2.259, 2.754, 3.249, 3.756, 4.251, 4.758, 5.253, 5.748, 6.243, 6.750, 7.245, 7.752, 8.247])	
# Molecular gas surface density
Msunpc2_SigmaH2 = np.array([9.557, 8.997, 4.844, 4.844, 4.700, 2.569, 1.871, 2.346, 0.9629, 1.022, 0.4735, 0.2791, 0.4197, 0.3554, 0.3056, 0.1436, 0.1094])

# H2 surface density of disc (from Van Eck+2015, Tab. 2)
Msunpc2_SigmaH2_Vaneck = np.array([1.90, 1.28])

def s(kpc_x,Ups):
    mu = mu0 + 1.10857 * kpc_x / kpc_Rd
    return Ups * 10**(-0.4 * (mu - C)) 

MsunLsun_Ups = [0.72,0.52]  #see Kam+17 Sec 5.1
C = 24.8
mu0 = 18.01
kpc_Rd = 1.82

nR = 10001
nU = len(MsunLsun_Ups)

#Specify range 
xr=[0,20]
yr=[0,400]
 
kpc_xbnd_M33 = np.array([1,3,5])	#radial bins for magnetic field data (Tabatabaei+08)
kpc_xmid_M33 = np.array([2,4])	#midpoints of radial bins for magnetic field data (Tabatabaei+08)

Delta_kpc_R = xr[1]
shift_kpc_R = 0

kpc_R = np.arange(nR) / nR * Delta_kpc_R + shift_kpc_R

Msunpc2_Sigma_star = np.zeros((nR,nU))

for j in range(nU):
    Msunpc2_Sigma_star[:,j] = s(kpc_R,MsunLsun_Ups[j])

kpc_r_kam = arcmin_r_Kam / arcmin_deg / deg_rad * kpc_D_M33_Kam

rad_data = [kpc_r_kam, kpc_r_SFR, kpc_r_SigmaH2]
kpc_r = rad_data[np.argmin(np.array([d.size for d in rad_data]))]

dat_sigmastar = griddata(kpc_R, Msunpc2_Sigma_star[:,0]*g_Msun/(pcm)**2, kpc_r, method='linear', fill_value=nan, rescale=False)
dat_sigma = griddata(kpc_r_kam, Msunpc2_SigmaHI_Kam*g_Msun/(pcm)**2, kpc_r, method='linear', fill_value=nan, rescale=False)
dat_sigmah2 = griddata(kpc_r_SigmaH2, Msunpc2_SigmaH2 *g_Msun/(pcm)**2, kpc_r, method='linear', fill_value=nan, rescale=False)

dat_sigmasfr = Msunpc2Gyr_Sigma_SFR*g_Msun/((10**9*365*24*60*60)*(pcm)**2)
kmskpc_Om = kms_vcirc_Kam /kpc_r_kam
dat_omega = kmskpc_Om*1e+5/kpcm
dat_q = -1 * kpc_r_kam/ kmskpc_Om * np.gradient(kmskpc_Om)/np.gradient(kpc_r_kam)
dat_omega = griddata(kpc_r_kam, dat_omega, kpc_r, method='linear', fill_value=nan, rescale=False)
dat_q = griddata(kpc_r_kam, dat_q, kpc_r, method='linear', fill_value=nan, rescale=False)


dat_sigmasfr = griddata(kpc_r_SFR, dat_sigmasfr, kpc_r, method='linear', fill_value=nan, rescale=False)

molbool = True
if molbool:
    dat_sigma += dat_sigmah2

dat_sigmatot = dat_sigma + dat_sigmastar


data  = kpc_r, dat_sigmatot, dat_sigma, dat_q, dat_omega, dat_sigmasfr
with open('data_m33.pickle', 'wb') as f:
    pickle.dump(data, f)
print(kpc_r.size)