import numpy as np
from sympy import *
from fractions import Fraction
import pickle
from scipy.interpolate import griddata
import pandas as pd
import math as m #for inclination angle correction
import os
import math as m

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
        df = pd.read_csv(current_directory+r'\data\6946 data\{}.dat'.format(i),delimiter=',')
        dataframe_list.append(df)
    return dataframe_list

#extrapolation
def interpolation(list1,list2,standard):
    extrapolated_data = griddata(list1, list2, standard, method='linear', fill_value=nan, rescale=False)
    return extrapolated_data
###########################################################################################################################################
    
# file_names=['smdf','HI 6946','H2 6946','q_values6946sofue+18','omega_sofue+18',
#             'SFR_Halpha24 6946','SFR_FUV24 6946','velocity disp.']
file_names=['smdf','sigma_HI_bigiel','sigma_H2_bigiel','q_values6946sofue+18','omega_sofue+18',
            'sigma_sfr_bigiel','velocity disp.']
dataframe_list=file_reader(file_names)

#to obtain radius data from every df
distance_6946= 7.72 #Mpc
distances_Mpc=[5.5,5.9,5.9,5.5,5.5,5.9,6]

radius_list=[np.array(dataframe_list[i]['r']) for i in range(len(dataframe_list))]
radius_list=[radius_list[i]*(distance_6946/distances_Mpc[i]) for i in range(len(radius_list))] #applying distance corrections

#obtain arrays of quantities
col_names=['smdf','sigma_HI','sigma_H2','q','omega','sigma_sfr','vel disp']
conv_factors=[(g_Msun/(cm_pc**2) ),g_Msun/(cm_pc**2),g_Msun/(cm_pc**2),1,
              cm_km/cm_kpc,g_Msun/((s_Myr*10**(-6))*(cm_kpc**2)),(3**0.5)]

#to switch between fuv and h_alpha data for sfr
# sfr_val=0
# tbd=[col_names,dataframe_list,conv_factors,radius_list]
# if sfr_val==0: #choosing H_alpha data
#     for i in range(len(tbd)):
#         del tbd[i][5]
# else: #chose FUV data
#     for i in range(len(tbd)):
#         del tbd[i][4]

#inclination correction for quantities
inclinations=[30,33,33,30,30,33,38]
i_6946=m.radians(38) #in radians
quant_list=[np.array(dataframe_list[i][col_names[i]])*(m.cos(i_6946)/m.cos(m.radians(inclinations[i]))) 
                     for i in range(len(dataframe_list))] #corrected for different inclination angles
quant_list=[quant_list[i]*conv_factors[i] for i in range(len(quant_list))] 

#find array with least length and correct it for radius
#if vel disp data is the smallest one, remove that and repeat the process
kpc_r = min(radius_list, key=len) #this shows that temp data has least number of points

#temperature fit #fit obtained using scipy.stats.linregress in temperature fitting.py
T=np.array([(266.33525528933166*r)+5749.778318143174 for r in kpc_r]) 
error_temp=np.std(T)

#interpolating the data and appending the common radius list 
#interpolated arrays has an _ip ending
quant_list_ip=[interpolation(radius_list[i],quant_list[i],kpc_r) for i in range(len(dataframe_list))]
quant_list_ip.insert(len(quant_list_ip)-1,T)
quant_list_ip.insert(0,(kpc_r))

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
with open(current_directory+'\data\data_{}.pickle'.format('6946'), 'wb') as f:
    pickle.dump(nandeleted_data, f)

########################################################################################################################################
                                                      #ERROR DATA#
########################################################################################################################################

#to load errors if any given in the dat files
#this assumes that in all files, the error columns are named as 'quant_error'
error_list=[]

#temp error
temp_error=np.array([error_temp]*len(radius_list[4]))
error_list.append(temp_error)

#RC error
RC_error=np.array([15*(cm_km/cm_kpc)]*len(radius_list[4]))
error_list.append(RC_error)






