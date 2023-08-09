import numpy as np
import matplotlib.pyplot as plt
from sympy import *
import pickle
import sys


from scipy.interpolate import griddata



pc_kpc = 1e3  # number of pc in one kpc
cm_km = 1e5  # number of cm in one km
s_day = 24*3600  # number of seconds in one day
s_min = 60  # number of seconds in one hour
s_hr = 3600  # number of seconds in one hour
cm_Rsun = 6.957e10  # solar radius in cm
g_Msun = 1.989e33  # solar mass in g
cgs_G = 6.674e-8
cms_c = 2.998e10
g_mH = 1.6736e-24
g_me = 9.10938e-28
cgs_h = 6.626e-27
deg_rad = 180e0/np.pi
arcmin_deg = 60e0
arcsec_deg = 3600e0
cm_kpc = 3.086e+21  # number of centimeters in one parsec
cm_pc = cm_kpc/1e+3
s_Myr = 1e+6*(365*24*60*60)  # megayears to seconds

kpc_D_M31_Plot = 780e0  # distance to M31 used for plots (see Beck+2019)

kpc_D_M31_Chemin = 785e0  # distance to M31 used by Chemin+2009
kpc_D_M31_TB10 = 780e0  # distance to M31 used by Tabatabaei+Berkhuijsen 2010
kpc_D_M31_Fletcher = 690e0  # distance to M31 used by Fletcher+2004
kpc_rmax_M31 = 46
# radial bins for magnetic field data (Fletcher+04)
kpc_xbnd_M31_Fletcher = np.array([6, 8, 10, 12, 14])
# midpoints of radial bins for magnetic field data (Fletcher+04)
kpc_xmid_M31_Fletcher = np.array([7, 9, 11, 13])

kpc_xbnd_M31 = kpc_xbnd_M31_Fletcher * kpc_D_M31_Plot / kpc_D_M31_Fletcher
kpc_xmid_M31 = kpc_xmid_M31_Fletcher * kpc_D_M31_Plot / kpc_D_M31_Fletcher

arcmin_xbnd_M31 = kpc_xbnd_M31 / kpc_D_M31_Plot * deg_rad * arcmin_deg
arcmin_xmid_M31 = kpc_xmid_M31 / kpc_D_M31_Plot * deg_rad * arcmin_deg


length, breadth = [5, 2.5]

# Specifiy text size in legend
leg_textsize = 10
axis_textsize = 10
rc = {"font.family": "serif",
      "mathtext.fontset": "stix"}
plt.rcParams.update(rc)
plt.rcParams["font.serif"] = ["Times New Roman"] + plt.rcParams["font.serif"]


# M31

# HI data from Claude Carignan (private communication--see data.pdf for additional info.

# Reference: Chemin+2009 (ApJ, 705, 1395)
# Distance to M31 assumed in that paper: 785 +/- 25 kpc
# velocity dispersion data copied from m31_HIdata_velocdisp.txt
# Circular velocity data copied from model2c copie.dat

# Radius for velocity dispersion data (from m31_HIdata_velocdisp.txt)
arcsec_r = np.array([100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0, 1100.0, 1200.0, 1300.0, 1400.0, 1500.0, 1600.0, 1700.0, 1800.0, 1900.0, 2000.0, 2100.0, 2200.0, 2300.0, 2400.0, 2500.0, 2600.0, 2700.0, 2800.0, 2900.0, 3000.0, 3100.0, 3200.0, 3300.0, 3400.0, 3500.0, 3600.0, 3700.0, 3800.0, 3900.0, 4000.0, 4100.0, 4200.0, 4300.0, 4400.0, 4500.0, 4600.0, 4700.0, 4800.0, 4900.0, 5000.0,
                    5100.0, 5200.0, 5300.0, 5400.0, 5500.0, 5600.0, 5700.0, 5800.0, 5900.0, 6000.0, 6100.0, 6200.0, 6300.0, 6400.0, 6500.0, 6600.0, 6700.0, 6800.0, 6900.0, 7000.0, 7100.0, 7200.0, 7300.0, 7400.0, 7500.0, 7600.0, 7700.0, 7800.0, 7900.0, 8000.0, 8100.0, 8200.0, 8300.0, 8400.0, 8500.0, 8600.0, 8700.0, 8800.0, 8900.0, 9000.0, 9100.0, 9200.0, 9300.0, 9400.0, 9500.0, 9600.0, 9700.0, 9800.0, 9900.0, 10000.0])

# Radius from Chemin+2009 (Tab. 4)
arcmin_r_Chemin = np.array([1.67, 3.33, 5.00, 6.67, 8.33, 10.00, 11.67, 13.33, 15.00, 16.67, 18.33, 20.00, 21.67, 23.33, 25.00, 26.67, 28.33, 30.00, 31.67, 33.33, 35.00, 36.67, 38.33, 40.00, 41.67, 43.33, 45.00, 46.67, 48.33, 50.00, 51.67, 53.33, 55.00, 56.67, 58.33, 60.00, 61.67, 63.33, 65.00, 66.67, 68.33, 70.00, 71.67, 73.33, 75.00, 76.67, 78.33, 80.00, 81.67, 83.33, 85.00, 86.67,
                           88.33, 90.00, 91.67, 93.33, 95.00, 96.67, 98.33, 100.00, 101.67, 103.33, 105.00, 106.67, 108.33, 110.00, 111.67, 113.33, 115.00, 116.67, 118.33, 120.00, 121.67, 123.33, 125.00, 126.67, 128.33, 130.00, 131.67, 133.33, 135.00, 136.67, 138.33, 140.00, 141.67, 143.33, 145.00, 146.67, 148.33, 150.00, 151.67, 153.33, 155.00, 156.67, 158.33, 160.00, 161.67, 163.33, 165.00, 166.67])

# Vcirc from Chemin+2009 (Tab. 4)
kms_vcirc_Chemin = np.array([0, 0, 336.2, 324.6, 339.0, 243.6, 235.2, 238.9, 239.3, 226.3, 202.6, 207.3, 202.5, 208.9, 221.6, 232.2, 237.6, 239.8, 235.6, 241.7, 244.3, 248.8, 251.8, 253.0, 258.8, 259.0, 262.2, 266.8, 266.8, 265.9, 264.4, 264.7, 265.3, 265.2, 262.0, 260.8, 259.2, 258.1, 258.4, 259.2, 262.7, 266.1, 270.0, 269.8, 269.1, 268.5, 263.0, 257.1, 254.1,
                            251.9, 249.5, 245.7, 243.7, 245.9, 242.3, 239.2, 239.5, 236.1, 233.8, 233.1, 230.1, 232.1, 228.7, 229.1, 227.9, 226.9, 225.1, 225.4, 230.3, 229.0, 229.9, 230.1, 229.8, 230.4, 230.9, 229.8, 228.8, 238.3, 243.6, 247.3, 247.8, 248.4, 248.1, 244.5, 244.4, 241.7, 237.7, 237.6, 244.9, 247.9, 256.3, 253.5, 244.3, 249.3, 255.7, 255.0, 271.1, 269.8, 258.2, 275.1])

# Vcirc_error from Chemin+2009 (Tab. 4)
kms_vcirc_error_Chemin = np.array([0, 0, 171.7, 125.1, 52.8, 25.8, 17.0, 5.7, 18.3, 16.1, 4.7, 10.7, 21.7, 15.6, 13.4, 13.7, 8.3, 2.2, 6.1, 3.3, 6.4, 6.4, 5.5, 9.2, 9.6, 9.5, 10.8, 13.0, 11.7, 9.9, 7.2, 5.3, 5.4, 3.9, 2.4, 4.5, 5.4, 6.8, 6.1, 4.6, 4.2, 4.0, 2.3, 1.0, 3.9, 7.1,
                                  15.0, 13.5, 10.6, 11.1, 8.1, 7.4, 6.6, 7.5, 6.1, 6.3, 4.7, 1.8, 1.7, 3.3, 5.5, 5.0, 1.7, 1.8, 5.0, 2.0, 1.6, 1.8, 2.0, 2.3, 4.8, 6.6, 3.0, 5.2, 2.9, 2.1, 1.8, 3.3, 1.4, 3.1, 1.3, 2.0, 1.5, 1.6, 3.0, 4.3, 1.9, 6.1, 3.6, 3.2, 3.1, 4.1, 4.7, 5.8, 4.5, 5.8, 7.8, 4.7, 10.7, 4.8])

# HI surface density from Chemin+2009 (Tab. 4)
Msunpc2_SigmaHI_Chemin = np.array([2.77, 3.19, 4.21, 3.94, 4.49, 2.14, 1.52, 1.51, 1.57, 2.21, 3.00, 3.54, 2.75, 2.31, 2.09, 1.86, 1.41, 1.43, 1.60, 2.07, 2.84, 3.85, 4.61, 4.39, 3.78, 3.34, 3.16, 3.37, 4.34, 5.02, 4.95, 4.45, 4.03, 3.98, 3.91, 3.67, 3.39, 3.37, 3.01, 2.61, 2.39, 2.00, 1.92, 1.79, 2.21, 1.91, 1.74,
                                  1.60, 1.44, 1.35, 1.24, 1.15, 1.18, 1.10, 1.08, 0.83, 0.75, 0.76, 0.88, 0.82, 0.75, 0.67, 0.58, 0.49, 0.43, 0.42, 0.45, 0.46, 0.42, 0.39, 0.38, 0.29, 0.23, 0.20, 0.21, 0.17, 0.10, 0.08, 0.07, 0.07, 0.10, 0.08, 0.07, 0.06, 0.06, 0.06, 0.02, 0.01, 0.01, 0.01, 0.02, 0.03, 0.02, 0.02, 0.02, 0.02, 0.03, 0.03, 0.03, 0.02])

# HI velocity dispersion: sigma_los using model WITHOUT warp (C. Carignan, private communication; Claude suggested I should use this model)
kms_sigmaLOS = np.array([17.10, 8.58, 12.94, 17.12, 12.99, 10.03, 13.13, 15.38, 16.60, 14.95, 15.63, 14.50, 12.63, 12.25, 12.57, 12.51, 13.53, 13.79, 13.93, 13.21, 11.56, 10.87, 10.81, 11.06, 11.77, 11.12, 12.21, 12.12, 12.12, 11.70, 11.41, 11.70, 12.54, 12.15, 11.65, 11.18, 11.21, 11.36, 11.72, 10.62, 10.23, 10.47, 10.62, 10.65, 9.76, 10.82,
                        11.45, 10.92, 10.91, 10.46, 10.01, 10.04, 9.86, 9.94, 10.17, 9.86, 9.55, 9.72, 9.79, 10.15, 9.82, 9.76, 9.86, 9.37, 9.73, 9.58, 9.66, 10.21, 10.21, 9.20, 8.95, 9.59, 10.20, 8.94, 9.01, 9.12, 9.29, 9.46, 9.80, 8.34, 9.01, 8.46, 8.98, 7.70, 7.99, 7.04, 6.43, 5.80, 8.92, 8.76, 6.95, 6.81, 6.60, 9.30, 13.69, 11.38, 6.81, 4.58, 9.05, 10.89])

# HI velocity dispersion: sigma_los using model WITH warp (C. Carignan, private communication)
kms_sigmaLOS_warp = np.array([12.05, 15.64, 15.70, 10.75, 8.71, 5.82, 6.08, 8.23, 9.48, 11.61, 13.61, 12.37, 7.77, 6.22, 5.98, 5.38, 5.47, 6.62, 8.16, 8.31, 8.67, 9.30, 9.67, 8.51, 7.51, 6.32, 6.13, 6.52, 9.29, 10.66, 11.48, 10.77, 10.61, 11.23, 12.05, 11.80, 11.42, 10.91, 8.56, 6.32, 5.45, 4.81, 5.43, 6.68, 10.43, 10.71,
                             10.95, 10.25, 9.16, 8.87, 9.05, 8.72, 8.98, 8.49, 9.61, 6.84, 6.16, 6.88, 9.87, 9.87, 9.85, 8.62, 7.97, 7.54, 7.19, 7.99, 8.30, 8.56, 7.14, 7.23, 7.85, 7.05, 7.49, 6.64, 8.34, 7.07, 5.53, 5.61, 4.06, 4.75, 8.15, 7.09, 5.89, 6.24, 5.59, 5.59, 4.49, 3.81, 2.61, 3.57, 4.67, 5.37, 6.67, 4.19, 5.09, 4.45, 4.40, 4.53, 5.62, 5.20])

# Radius for circular velocity and surface density data (from model2c copie.dat) -- assumes distance of 785 kpc
kpc_r_Claude = np.array([0.00, 0.48, 0.96, 1.44, 1.92, 2.40, 2.88, 3.36, 3.84, 4.32, 4.80, 5.28, 5.76, 6.24, 6.72, 7.20, 7.68, 8.16, 8.64, 9.12, 9.60, 10.08, 10.56, 11.04, 11.52, 12.00, 12.48, 12.96, 13.44, 13.92, 14.40, 14.88, 15.36, 15.84, 16.32, 16.80, 17.28, 17.76, 18.24, 18.72, 19.20, 19.68, 20.16, 20.64, 21.12, 21.60, 22.08, 22.56,
                        23.04, 23.52, 24.00, 24.48, 24.96, 25.44, 25.92, 26.40, 26.88, 27.36, 27.84, 28.32, 28.80, 29.28, 29.76, 30.24, 30.72, 31.20, 31.68, 32.16, 32.64, 33.12, 33.60, 34.08, 34.56, 35.04, 35.52, 36.00, 36.48, 36.96, 37.44, 37.92, 38.40, 38.88, 39.36, 39.84, 40.32, 40.80, 41.28, 41.76, 42.24, 42.72, 43.20, 43.68, 44.16, 44.64, 45.12, 45.60])

# Circular velocity (from the VEL. TOTAL column of model2c copie.dat, which seems to be a smoothed version of the entries in Tab. 4 of Chemin+09?)
kms_vcirc = np.array([0.00, 150.97, 221.65, 247.13, 252.51, 253.89, 252.25, 251.54, 251.64, 252.07, 250.78, 246.00, 240.36, 236.65, 235.64, 235.92, 235.76, 233.72, 230.73, 229.63, 228.34, 228.24, 230.80, 233.99, 236.79, 238.75, 239.80, 240.51, 241.17, 241.85, 242.65, 244.03, 245.74, 247.20, 247.87, 247.59, 246.79, 245.90, 245.12, 244.22, 243.14, 242.09, 241.13, 240.15, 239.12, 238.18, 237.41,
                     236.67, 235.80, 234.78, 233.72, 232.66, 231.61, 230.55, 229.56, 228.63, 227.70, 226.84, 226.10, 225.41, 224.72, 224.02, 223.31, 222.62, 221.94, 221.28, 220.64, 220.02, 219.42, 218.84, 218.28, 217.74, 217.22, 216.72, 216.24, 215.79, 215.36, 214.95, 214.55, 214.19, 213.84, 213.51, 213.22, 212.93, 212.67, 212.43, 212.20, 212.00, 211.82, 211.66, 211.53, 211.42, 211.32, 211.19, 211.05, 210.85])

# HI surface density of disc (from model2c copie.dat column HI S.D.)
Msunpc2_SigmaHI = np.array([0.1212E-02, 0.1975E+00, 0.3949E+00, 0.5924E+00, 0.7899E+00, 0.9873E+00, 0.1185E+01, 0.1382E+01, 0.1596E+01, 0.1816E+01, 0.2046E+01, 0.2282E+01, 0.2491E+01, 0.2678E+01, 0.2797E+01, 0.2836E+01, 0.2960E+01, 0.3228E+01, 0.3525E+01, 0.3895E+01, 0.4232E+01, 0.4438E+01, 0.4621E+01, 0.4700E+01, 0.4771E+01, 0.4829E+01, 0.4886E+01, 0.4926E+01, 0.4960E+01, 0.4721E+01, 0.4424E+01, 0.3954E+01, 0.3421E+01, 0.2932E+01, 0.2467E+01, 0.2110E+01, 0.1831E+01, 0.1661E+01, 0.1583E+01, 0.1537E+01, 0.1526E+01, 0.1506E+01, 0.1468E+01, 0.1429E+01, 0.1391E+01, 0.1334E+01, 0.1240E+01,
                           0.1134E+01, 0.1012E+01, 0.8889E+00, 0.7868E+00, 0.6936E+00, 0.6259E+00, 0.5683E+00, 0.5107E+00, 0.4532E+00, 0.3956E+00, 0.3411E+00, 0.2969E+00, 0.2634E+00, 0.2298E+00, 0.1962E+00, 0.1736E+00, 0.1581E+00, 0.1437E+00, 0.1293E+00, 0.1149E+00, 0.1006E+00, 0.8617E-01, 0.7176E-01, 0.5736E-01, 0.4317E-01, 0.2953E-01, 0.2473E-01, 0.1993E-01, 0.1513E-01, 0.1034E-01, 0.3429E-02, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00])

# HI surface density of disc (from Van Eck+2015, Tab. 2 -- also used in Chamandy+2016)
Msunpc2_SigmaHI_Vaneck = np.array([1.47, 2.17, 3.64, 4.05])

# Total surface density of disc (from model2c copie.dat column DISK S.D.)
Msunpc2_SigmaTot = np.array([0.7652E+04, 0.4847E+04, 0.3022E+04, 0.1854E+04, 0.1275E+04, 0.9529E+03, 0.7519E+03, 0.6224E+03, 0.5232E+03, 0.4383E+03, 0.3604E+03, 0.3059E+03, 0.2761E+03, 0.2597E+03, 0.2437E+03, 0.2238E+03, 0.2018E+03, 0.1831E+03, 0.1741E+03, 0.1668E+03, 0.1589E+03, 0.1589E+03, 0.1530E+03, 0.1444E+03, 0.1340E+03, 0.1230E+03, 0.1138E+03, 0.1059E+03, 0.9915E+02, 0.9326E+02, 0.8844E+02, 0.8397E+02, 0.7813E+02, 0.7083E+02, 0.6349E+02, 0.5694E+02, 0.5197E+02, 0.4787E+02, 0.4386E+02, 0.3959E+02, 0.3629E+02, 0.3349E+02, 0.3085E+02, 0.2815E+02, 0.2599E+02, 0.2420E+02, 0.2252E+02,
                            0.2042E+02, 0.1849E+02, 0.1691E+02, 0.1558E+02, 0.1430E+02, 0.1311E+02, 0.1205E+02, 0.1141E+02, 0.1040E+02, 0.9532E+01, 0.9259E+01, 0.8529E+01, 0.7806E+01, 0.7141E+01, 0.6533E+01, 0.5976E+01, 0.5466E+01, 0.4998E+01, 0.4568E+01, 0.4180E+01, 0.3826E+01, 0.3500E+01, 0.3201E+01, 0.2925E+01, 0.2677E+01, 0.2450E+01, 0.2241E+01, 0.2048E+01, 0.1874E+01, 0.1715E+01, 0.1568E+01, 0.1433E+01, 0.1313E+01, 0.1201E+01, 0.1096E+01, 0.1004E+01, 0.9188E+00, 0.8388E+00, 0.7689E+00, 0.7030E+00, 0.6424E+00, 0.5885E+00, 0.5375E+00, 0.4922E+00, 0.4503E+00, 0.4112E+00, 0.1963E+00, 0.0000E+00, 0.0000E+00])

# SFR surface density data from Tabatabaei+Berkhuijsen 2010

# Radius assuming distance to M31 of D = 780 +/- 40 kpc
kpc_r_SFR_TB10 = np.array([6.25, 6.75, 7.25, 7.75, 8.25, 8.75, 9.25, 9.75, 10.25, 10.75,
                          11.25, 11.75, 12.25, 12.75, 13.25, 13.75, 14.25, 14.75, 15.25, 15.75, 16.25, 16.75])

# Star formation rate surface density
Msunpc2Gyr_Sigma_SFR = np.array([0.40, 0.43, 0.52, 0.41, 0.43, 0.56, 0.68, 0.75, 0.86,
                                1.00, 0.79, 0.54, 0.33, 0.22, 0.19, 0.16, 0.12, 0.07, 0.08, 0.10, 0.09, 0.09])

# Molecular gas fraction (from Nieten+06 Fig. 10, N+S data points, obtained using Dexter software)

# Radius for molecular gas fraction data
arcmin_r_MolFrac = np.array([5.710, 6.744, 7.711, 8.698, 9.688, 10.69, 11.70, 12.70, 13.69, 14.68, 15.66, 16.65, 17.65, 18.64, 19.71, 20.70, 21.70, 22.69, 23.68, 24.69, 25.69, 26.68, 27.65, 28.64, 29.69, 30.68, 30.68, 31.67, 32.66, 33.66, 34.65, 35.64, 36.63, 37.70, 38.62, 39.61,
                            40.68, 41.68, 42.67, 43.67, 44.66, 45.65, 46.65, 47.64, 48.63, 49.62, 50.68, 51.68, 52.67, 53.66, 54.65, 55.64, 56.63, 57.62, 58.68, 59.60, 60.59, 61.58, 62.65, 63.57, 64.63, 65.69, 66.61, 67.61, 68.53, 69.52, 70.59, 71.58, 72.57, 73.56, 74.55, 75.62, 76.61, 77.68, 78.60, 79.59])

# Molecular gas fraction data (error bars in Nieten+06 vary, but in the range we are interested, are between ~0.003 and 0.014)
MolFrac = np.array([0.1084, 0.3376, 0.1946, 0.1561, 0.1423, 0.2400, 0.3307, 0.3500, 0.3576, 0.3107, 0.2238, 0.2415, 0.2569, 0.2723, 0.3061, 0.2830, 0.2984, 0.2761, 0.2846, 0.3692, 0.4200, 0.3869, 0.2853, 0.2292, 0.1953, 0.1453, 0.1461, 0.1307, 0.1576, 0.1615, 0.1630, 0.1446, 0.1361, 0.1653, 0.1569, 0.1530, 0.1715, 0.1892, 0.2061,
                   0.2169, 0.2107, 0.2046, 0.2084, 0.2107, 0.2023, 0.1930, 0.1846, 0.1838, 0.1746, 0.1684, 0.1507, 0.1292, 0.1069, 0.08923, 0.07384, 0.06538, 0.06307, 0.06000, 0.05000, 0.04307, 0.03615, 0.03461, 0.03384, 0.03615, 0.03384, 0.03230, 0.03230, 0.02846, 0.02384, 0.02153, 0.01230, 0.009230, 0.008461, 0.01846, 0.01615, 0.01307])

# Values used in Van Eck+15/Chamandy+16
Msunpc2Gyr_Sigma_SFR_Vaneck = np.array([0.443, 0.621, 0.794, 0.227])

kpc_r_Chemin_orig = arcmin_r_Chemin / arcmin_deg / deg_rad * kpc_D_M31_Chemin

kpc_r_molfrac = arcmin_r_MolFrac / arcmin_deg / deg_rad * kpc_D_M31_Plot

# correct to distance used for our plots
kpc_r_cl = kpc_r_Claude * kpc_D_M31_Plot / kpc_D_M31_Chemin

kpc_r_Chemin_orig = arcmin_r_Chemin / arcmin_deg / deg_rad * kpc_D_M31_Chemin
kpc_r_Chemin = kpc_r_Chemin_orig * kpc_D_M31_Plot / \
    kpc_D_M31_Chemin  # correct to distance used for our plots

kpc_r_SFR = kpc_r_SFR_TB10 * kpc_D_M31_Plot / \
    kpc_D_M31_TB10  # correct to distance used for our plots

if __name__ == '__main__':
    # Select which data to use for sigma, q and omega
    chemin_flag = bool(int(sys.argv[1]))

    if chemin_flag:
        print(chemin_flag)
        kpc_r_f = kpc_r_Chemin
        Msunpc2_SigmaHI_f = Msunpc2_SigmaHI_Chemin
        kms_vcirc_f = kms_vcirc_Chemin
    else:
        kpc_r_f = kpc_r_cl
        Msunpc2_SigmaHI_f = Msunpc2_SigmaHI
        kms_vcirc_f = kms_vcirc

    argr = np.argwhere(kpc_r_f>0)
    kpc_r_f = kpc_r_f[argr].flatten()
    Msunpc2_SigmaHI_f = Msunpc2_SigmaHI_f[argr].flatten()
    kms_vcirc_f = kms_vcirc_f[argr].flatten()
    rad_data = [kpc_r_cl, kpc_r_Chemin, kpc_r_SFR, kpc_r_molfrac]
    kpc_r = rad_data[np.argmin(np.array([d.size for d in rad_data]))]

    dat_sigmatot = griddata(kpc_r_cl, Msunpc2_SigmaTot*g_Msun/(cm_pc)
                            ** 2, kpc_r, method='linear', fill_value=nan, rescale=False)
    dat_sigma = griddata(kpc_r_f, Msunpc2_SigmaHI_f*g_Msun/(cm_pc)
                        ** 2, kpc_r, method='linear', fill_value=nan, rescale=False)
    kmskpc_Om = kms_vcirc_f/kpc_r_f
    dat_omega = kmskpc_Om*cm_km/cm_kpc
    dat_q = -1 * kpc_r_f / kmskpc_Om * np.gradient(kmskpc_Om)/np.gradient(kpc_r_f)
    dat_q = griddata(kpc_r_f, dat_q, kpc_r, method='linear',
                    fill_value=nan, rescale=False)
    dat_omega = griddata(kpc_r_f, dat_omega, kpc_r,
                        method='linear', fill_value=nan, rescale=False)

    dat_sigmasfr = griddata(kpc_r_SFR, Msunpc2Gyr_Sigma_SFR*g_Msun/((1e+3*s_Myr)
                            * (cm_pc)**2), kpc_r, method='linear', fill_value=nan, rescale=False)
    molfrac = griddata(kpc_r_molfrac, MolFrac, kpc_r,
                    method='linear', fill_value=nan, rescale=False)
    # plt.plot(kpc_r, dat_sigmasfr)
    # plt.show()

    data = kpc_r, dat_sigmatot, dat_sigma, dat_q, dat_omega, dat_sigmasfr, molfrac
    # kpc_r, dat_sigmatot, dat_sigma, dat_q, dat_omega, dat_sigmasfr, molfrac
    with open('data_m31.pickle', 'wb') as f:
        pickle.dump(data, f)
