import numpy as np
from sympy import *
from fractions import Fraction
import pickle
from scipy.interpolate import griddata
import pandas as pd
import os

current_directory = str(os.getcwd())

#converting from 2nd unit to 1st
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

#REQUIRED FUNCTIONS
###########################################################################################################################################
#importing data from csv file to np array
def file_reader(list_of_file_names):
    dataframe_list=[]
    for i in list_of_file_names:
        df = pd.read_csv(current_directory+r'\data\M51 data\{}.dat'.format(i),delimiter=',')
        dataframe_list.append(df)
    return dataframe_list

#extrapolation
def interpolation(list1,list2,standard):
    interpolated_data = griddata(list1, list2, standard, method='linear', fill_value=nan, rescale=False)
    return interpolated_data
###########################################################################################################################################

file_names=['smdf','sigma_HI_bigiel','sigma_H2_bigiel','q_valuesM51sofue+18','omega_sofue+18',
            'sigma_sfr_bigiel','HI vel dispersion Hitschfeld']
# file_names=['smdf','HI m51..','H2 m51.','q_valuesM51sofue+18','omega_sofue+18',
#             'SFR_Halpha24 m51.','SFR_FUV24 m51.','CO vel dispersion schuster.']

dataframe_list=file_reader(file_names)

#to obtain radius data from every df
distance_m51= 8.5 #Mpc
distances_Mpc=[9.6,8,8,9.6,9.6,8,8.4]
for i in distances_Mpc:
    print(i)
radius_list=[np.array(dataframe_list[i]['r']) for i in range(len(dataframe_list))]
radius_list=[radius_list[i]*(distance_m51/distances_Mpc[i]) for i in range(len(radius_list))] #applying distance corrections

#obtain arrays of quantities
col_names=['smdf','sigma_HI','sigma_H2','q','omega',
           'sigma_sfr','vel disp']
conv_factors=[(g_Msun/(cm_pc**2) ), g_Msun/(cm_pc**2),g_Msun/(cm_pc**2), 1,
              cm_km/cm_kpc, g_Msun/((s_Myr*10**(-6))*(cm_kpc**2)), (3**0.5)] #last term to make 3D vdisp

# to switch between fuv and h_alpha data for sfr
# sfr_val=0
# tbd=[col_names,dataframe_list,conv_factors,radius_list]
# if sfr_val==0: #choosing H_alpha data
#     for i in range(len(tbd)):
#         del tbd[i][6]
#         # del tbd[i][5]
# else: #chose FUV data
#     for i in range(len(tbd)):
#         del tbd[i][5]
#         # del tbd[i][4]

# col_names=tbd[0]
# dataframe_list=tbd[1]
# conv_factors=tbd[2]
# radius_list=tbd[3]

# get the quantities and convert to cgs units
quant_list=[np.array(dataframe_list[i][col_names[i]]) for i in range(len(dataframe_list))]
quant_list=[quant_list[i]*conv_factors[i] for i in range(len(quant_list))] #list with 8 elements including vel disp 

#find array with least length and correct it for radius
#if vel disp data is the smallest one, remove that and repeat the process
kpc_r = min(radius_list, key=len) 
#temp data #fit obtained using scipy.stats.linregress in temperature fitting.py
T=np.array([((280.1606531283593*r)+4794.905648749489) for r in kpc_r])
error_temp=np.std(T)

#interpolating the data and appending the common radius list 
#interpolated arrays has an _ip ending
quant_list_ip=[interpolation(radius_list[i],quant_list[i],kpc_r) for i in range(len(col_names))]
quant_list_ip.insert(len(quant_list_ip)-1,T)
quant_list_ip.insert(0,kpc_r)

#removing nan values for points whr interpolation is impossible
# nan_max = np.argmax([np.sum(np.isnan(d)) for d in quant_list_ip])
# nan_max_data = quant_list_ip[nan_max]
# nan_mask = ~np.isnan(nan_max_data)

# nandeleted_data = []
# for i,d in enumerate(quant_list_ip):
#     nandeleted_data.append(d[nan_mask])

#removing nan values for points whr interpolation is impossible
nan_mask = np.zeros(kpc_r.size)
for d in quant_list_ip:
    nan_mask += np.isnan(d)
nan_mask = ~(nan_mask>0)
nandeleted_data = []
for i,d in enumerate(quant_list_ip):
    nandeleted_data.append(d[nan_mask])

#separating vel disp data to maintain uniformity
vel_disp=nandeleted_data[-1]
del nandeleted_data[-1]
nandeleted_data = tuple(nandeleted_data)

#storing that in a pickle file
#nandeleted_data follows same order of arrays as M31 and M33 data
#kpc_r, dat_sigmatot, dat_sigmaHI,dat_sigmaH2, dat_q, dat_omega, dat_sigmasfr, T= data
with open(current_directory+'\data\data_{}.pickle'.format('m51'), 'wb') as f:
    pickle.dump(nandeleted_data, f)
########################################################################################################################################
                                                      #ERROR DATA#
########################################################################################################################################

#to load errors if any given in the dat files
#this assumes that in all files, the error columns are named as 'quant_error'
error_list=[]

temp_error=np.array([error_temp]*len(radius_list[-2])) # chose len(radius_list[4]) to match err_Om calc in error_approximations.py
error_list.append(temp_error)


RC_error=np.array([15*(cm_km/cm_kpc)]*len(radius_list[4]))
error_list.append(RC_error)







