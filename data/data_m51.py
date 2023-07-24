import numpy as np
from sympy import *
from fractions import Fraction
import pickle
from scipy.interpolate import griddata
import pandas as pd
# current_directory=r'D:\Documents\Gayathri_college\MSc project\codes\data'
current_directory=r'D:\Documents\Gayathri_college\MSc project\codes\MSc.-Thesis\data\M51 data'

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
        df = pd.read_csv(current_directory+r'\{}.dat'.format(i),delimiter=',')
        dataframe_list.append(df)
    return dataframe_list

#extrapolation
def interpolation(list1,list2,standard):
    extrapolated_data = griddata(list1, list2, standard, method='linear', fill_value=nan, rescale=False)
    return extrapolated_data
###########################################################################################################################################
    
file_names=['smdf','HI m51..','H2 m51.','q_valuesM51sofue+18','omega_sofue+18',
            'SFR_Halpha24 m51.','SFR_FUV24 m51.','temperature','CO vel dispersion schuster.']
dataframe_list=file_reader(file_names)

#to obtain radius data from every df
radius_list=[np.array(dataframe_list[i]['r']) for i in range(len(dataframe_list))]

#obtain arrays of quantities
col_names=['smdf','sigma_HI','sigma_H2','q','omega',
           'sigma_sfr','sigma_sfr_fuv','temp','vel disp']
conv_factors=[(g_Msun/(cm_pc**2) ),g_Msun/(cm_pc**2),g_Msun/(cm_pc**2),1,
              cm_km/cm_kpc,g_Msun/((s_Myr*10**(-6))*(cm_pc**2)),g_Msun/((s_Myr*10**(-6))*(cm_pc**2)),1,1]
quant_list=[np.array(dataframe_list[i][col_names[i]])*conv_factors[i] for i in range(len(dataframe_list))]

#to load errors if any given in the dat files
#this assumes that in all files, the error columns are named as 'quant_error'
error_list=[]
for i in dataframe_list:
    if 'quant_error' in i.columns:
        error_list.append(np.array(i['quant_error']))
    else:
        continue

#find array with least length and correct it for radius
array_with_least_length = min(radius_list, key=len) #this shows that temp data has least number of points
corrected_radius=array_with_least_length*(8.5/7.6) #using known correction for temperature data

#interpolating the data and appending the common radius list 
#interpolated arrays has an _ip ending
quant_list_ip=[np.round(interpolation(radius_list[i],quant_list[i],corrected_radius),3) for i in range(len(dataframe_list))]
quant_list_ip.insert(0,np.round(corrected_radius,3))

#removing nan values for points whr interpolation is impossible
nan_max = np.argmax([np.sum(np.isnan(d)) for d in quant_list_ip])
nan_max_data = quant_list_ip[nan_max]
nan_mask = ~np.isnan(nan_max_data)

nandeleted_data = []
for i,d in enumerate(quant_list_ip):
    nandeleted_data.append(d[nan_mask])
nandeleted_data = tuple(nandeleted_data)

#storing that in a pickle file
#nandeleted_data follows same order of arrays as M31 and M33 data
code_directory=r'D:\Documents\Gayathri_college\MSc project\codes\MSc.-Thesis\data'
with open(code_directory+'\data_{}.pickle'.format('m51'), 'wb') as f:
    pickle.dump(nandeleted_data, f)





