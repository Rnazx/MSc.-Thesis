import numpy as np
from sympy import *
from fractions import Fraction
import pickle
from scipy.interpolate import griddata
import pandas as pd
import os
import math as m

current_directory = str(os.getcwd())

#converting from 2nd unit to 1st
pc_kpc = 1e3  # number of pc in one kpc
cm_km = 1e5  # number of cm in one km
s_day = 24*3600  # number of seconds in one day
s_min = 60  # number of seconds in one minute
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
        df = pd.read_csv(current_directory+r'\data\M31 data\{}.dat'.format(i),delimiter=',')
        dataframe_list.append(df)
    return dataframe_list

#extrapolation
def interpolation(list1,list2,standard):
    extrapolated_data = griddata(list1, list2, standard, method='linear', fill_value=nan, rescale=False)
    return extrapolated_data
###########################################################################################################################################

file_names=['sigmatot_Claude_M31','sigmaHI_Chemin_M31','sigmaHI_Claude_M31','vcirc_chemin_M31','vcirc_chemin_M31','vcirc_claude_M31',
            'sigma_sfr_M31','molfrac_M31','vdisp_M31','vdisp_M31']  #'vdisp_M31' repeated since warp and no-warp data are in same file
                                                                    # 'vcirc_chemin_M31' repeated to get error for vcirc
dataframe_list=file_reader(file_names)

#to obtain radius data from every df
distance_m31= 0.78 #Mpc
distances_Mpc=[0.785,0.785,0.785,0.785,0.785,0.785,0.780,0.780,0.785,0.785]
r_col_names=['r kpc','r arcmin','r kpc','r arcmin','r arcmin','r kpc','r kpc','r arcmin','r arcsec','r arcsec']
radius_conv_corr=[distance_m31/distances_Mpc[0],(1/(arcmin_deg*deg_rad))*(distance_m31*1000),distance_m31/distances_Mpc[2],
                  (1/(arcmin_deg*deg_rad))*(distance_m31*1000),(1/(arcmin_deg*deg_rad))*(distance_m31*1000),distance_m31/distances_Mpc[5],
                  distance_m31/distances_Mpc[6],(1/(arcmin_deg*deg_rad))*(distance_m31*1000),(1/(60*arcmin_deg*deg_rad))*(distance_m31*1000),
                  (1/(60*arcmin_deg*deg_rad))*(distance_m31*1000)] 

#obtain arrays of quantities
col_names=['sigma_tot','sigma_HI','sigma_HI','vcirc kms','error kms','vcirc kms',
           'sigma_sfr','molfrac','vdisp kms no_warp','vdisp kms warp']
conv_factors=[(g_Msun/(cm_pc**2) ), g_Msun/(cm_pc**2), g_Msun/(cm_pc**2), 1,1,1,
              g_Msun/((s_Myr*10**(3))*(cm_pc**2)),1, (3**0.5),(3**0.5)] #last term to make 3D vdisp

#to switch between Claude and Chemin data for sigmaHI and vcirc
Claude_Chemin_switch=1 #Claude data chosen
tbd=[col_names,dataframe_list,conv_factors,distances_Mpc,r_col_names,radius_conv_corr]
if Claude_Chemin_switch==0: #choosing Chemin data
    for i in range(len(tbd)):
        del tbd[i][5] #remove vcirc_Claude
        del tbd[i][2] #remove sigmaHI_Claude
else: #chose Claude data
    for i in range(len(tbd)):
        del tbd[i][3] #remove vcirc_Chemin
        del tbd[i][1] #remove sigmaHI_Chemin

#redefining the lists without Chemin/Claude data
col_names        =tbd[0]
dataframe_list   =tbd[1]
conv_factors     =tbd[2]
distances_Mpc    =tbd[3]
r_col_names      =tbd[4]
radius_conv_corr =tbd[5]

# get the radii
radius_list=[np.array(dataframe_list[i][r_col_names[i]]) for i in range(len(dataframe_list))]
# convert all radii to kpc with distance correction
radius_list_kpc=[radius_list[i]*radius_conv_corr[i] for i in range(len(radius_list))]

#inclination angle correction
#need to get inclination used in Claude data
i_m31= 75 #deg
inclinations=[i_m31,77,i_m31,77,77,i_m31,75,77.5,i_m31,i_m31] #used i_m31 as no inclination correction is needed for Claude data

# get the quantities and convert to cgs units
quant_list=[np.array(dataframe_list[i][col_names[i]])*(m.cos(m.radians(i_m31))/m.cos(m.radians(inclinations[i]))) 
            for i in range(len(dataframe_list))]
quant_list=[quant_list[i]*conv_factors[i] for i in range(len(quant_list))]  

#find array with least length
#if vel disp data is the smallest one, remove that and repeat the process
kpc_r = min(radius_list_kpc, key=len) 

#temp data
T = (0.017*kpc_r + 0.5)*1e+4  #obtained from paper
error_temp=np.std(T)

#calculate omega
kpc_r_om_q=radius_list[3]
kmskpc_Om = quant_list[3]/radius_list[3]
dat_omega = kmskpc_Om*cm_km/cm_kpc
#calculate q
dat_q = -1 * radius_list[3]/ kmskpc_Om * np.gradient(kmskpc_Om)/np.gradient(radius_list[3])
#add omega, q and radius
quant_list.insert(3,dat_q)
radius_list_kpc.insert(3,radius_list[3])
quant_list.insert(4,dat_omega)
radius_list_kpc.insert(4,radius_list[3])

#remove vcirc data #now sfr in 5th pos
del quant_list[5]
del radius_list_kpc[5]

#interpolating the data and appending the common radius list 
#interpolated arrays has an _ip ending
quant_list_ip=[interpolation(radius_list_kpc[i],quant_list[i],kpc_r) for i in range(len(quant_list))]

#calculating sigmaH2 
# chose quantities after interpolation so that sigmaHI and molfrac are at same radii
dat_sigmaH2 = quant_list_ip[1]*(1/(1-quant_list_ip[-3]))
quant_list_ip.insert(2,dat_sigmaH2)
# removing molfrac data
molfrac=quant_list_ip[-3]
del quant_list_ip[-3]

quant_list_ip.insert(len(quant_list_ip)-2,T) # -2 because i need v_disp to be at end to remove easily in next step
quant_list_ip.insert(0,kpc_r)

#removing nan values for points whr interpolation is impossible
nan_max = np.argmax([np.sum(np.isnan(d)) for d in quant_list_ip])
nan_max_data = quant_list_ip[nan_max]
nan_mask = ~np.isnan(nan_max_data)
nandeleted_data = []
for i,d in enumerate(quant_list_ip):
    nandeleted_data.append(d[nan_mask])

#separating vel disp data since to maintain uniformity with other galaxies. To be used while plotting
dat_u=nandeleted_data[-2] #no warp data
dat_u_warp= nandeleted_data[-1]

del nandeleted_data[-1]
del nandeleted_data[-1]

#separating vcirc error since to maintain uniformity with other galaxies 
RC_error= nandeleted_data[4]
del nandeleted_data[4] #removing the vcirc error from nandeleted data
nandeleted_data = tuple(nandeleted_data)

#storing that in a pickle file
#nandeleted_data follows same order of arrays as M31 and M33 data
with open(current_directory+'\data\data_m31.pickle', 'wb') as f:
    pickle.dump(nandeleted_data, f)

########################################################################################################################################
                                                      #ERROR DATA#
########################################################################################################################################

#to load errors if any given in the dat files
#this assumes that in all files, the error columns are named as 'quant_error'
error_list=[]

temp_error=np.array([error_temp]*len(radius_list[4]))
error_list.append(temp_error)
error_list.append(RC_error)







