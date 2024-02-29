# import numpy as np
# import pickle
# import os
# import pandas as pd
# from helper_functions import parameter_read
# from data_helpers import *


# #converting from 2nd unit to 1st
# pc_kpc = 1e3  # number of pc in one kpc
# cm_km = 1e5  # number of cm in one km
# s_day = 24*3600  # number of seconds in one day
# s_min = 60  # number of seconds in one hour
# s_hr = 3600  # number of seconds in one hour
# cm_Rsun = 6.957e10  # solar radius in cm
# g_Msun = 1.989e33  # solar mass in g
# cgs_G = 6.674e-8
# cms_c = 2.998e10
# g_mH = 1.6736e-24
# g_me = 9.10938e-28
# cgs_h = 6.626e-27
# deg_rad = 180e0/np.pi
# arcmin_deg = 60e0
# arcsec_deg = 3600e0
# cm_kpc = 3.086e+21  # number of centimeters in one parsec
# cm_pc = cm_kpc/1e+3
# s_Myr = 1e+6*(365*24*60*60)  # megayears to seconds

# #########################################################################################
# #error in distance and inclination

# base_path = os.environ.get('MY_PATH')
# galaxy_name = os.environ.get('galaxy_name')

# #get parameter values
# params = parameter_read(os.path.join(base_path,'inputs','parameter_file.in'))
# switch = parameter_read(os.path.join(base_path,'inputs','switches.in'))

# current_directory = str(os.getcwd())
# os.chdir(os.path.join(base_path,'outputs'))

# with open(f'{galaxy_name}output_ca_'+str(params[r'C_\alpha'])+'K_'+str(params[r'K'])+'z_'+
#           str(params[r'\zeta'])+'psi_'+str(params[r'\psi'])+'b_'+str(params[r'\beta'])+'.out', 'rb') as f:
#                 model_f = pickle.load(f)

# os.chdir(os.path.join(base_path,'inputs'))

# #getting available data and error
# os.chdir(os.path.join(base_path, 'data','model_data', f'{galaxy_name}_data'))
# raw_data = pd.read_csv(f'combined_data_{galaxy_name}.csv', skiprows=1) #needed to obtain vcirc values which isnt present in interpolated_data.csv
# err_data= pd.read_csv(f'error_combined_{galaxy_name}.csv')
# corrections = pd.read_csv(f'correction_data_{galaxy_name}.csv', skiprows=1, index_col=0)


# os.chdir(os.path.join(base_path,'data'))

# # kpc_r, dat_sigmatot, dat_sigmaHI, dat_sigmaH2, dat_q, dat_omega, dat_sigmasfr, T= data
# data_frame = pd.read_csv("data_interpolated.csv")
# columns_as_arrays = [np.array(data_frame[col]) for col in data_frame.columns]
# kpc_r = np.array(data_frame.iloc[:, 0])

# #distance correction
# new_dist= corrections.iloc[-1,0]
# err_new_dist=corrections.iloc[-1,2]
# old_dist = corrections.iloc[:-1,0].values
# err_old_dist=corrections.iloc[:-1,2].values

# #inclination correction
# new_i= corrections.iloc[-1,1] #deg
# err_new_i= corrections.iloc[-1,3] #deg
# old_i= corrections.iloc[:-1,1].values #used new_i as no inclination correction is needed for Claude data
# err_old_i=corrections.iloc[:-1,3].values

# #################################################################################################################
# #raw_data from combined_data.csv is converted to correct units, but no D/i correction done here
# raw_data = incl_distance_correction(raw_data, distance_new=old_dist, distance_old=old_dist,\
#                         i_new=np.radians(old_i), i_old=np.radians(old_i))

# os.chdir(os.path.join(base_path, 'data','model_data', f'{galaxy_name}_data'))

# #unwanted cols are removed from raw_data
# index_of_removed_data=[]
# try:
#     data_rem = pd.read_csv(f'removed_data_{galaxy_name}.csv', dtype=str)
# except:
#     data_rem = []
# for d in data_rem:
#     #change this
#     names=[i for i in raw_data.columns if d in i]
#     index_of_removed_data=[raw_data.columns.get_loc(i) for i in names]
#     remove_data(raw_data, d)
# #################################################################################################################

# #################################################################################################################
# columns_as_arrays_no_radius_temp = columns_as_arrays[1:-1]
# #reverting the correction for interpolated data to be used for error calculation


# if galaxy_name=='m31':
#         ip_list=[0,1,2,-1]
#         if switch['chem_or_claude']=='Claude':
#                 indices = [0,2,-1,-2] #these are indices for accessing old_inc values
#         else:
#                 indices = [0,1,-1,-2]
#         columns_as_arrays_reverted=[columns_as_arrays_no_radius_temp[ip_list[indices.index(i)]]*(np.cos(np.radians(old_i[i]))/np.cos(np.radians(new_i))) for i in indices]

# elif galaxy_name=='m33':
#         columns_as_arrays = [columns_as_arrays[i] for i in range(len(columns_as_arrays)) if i not in index_of_removed_data]
#         # columns_as_arrays_no_radius_temp = [columns_as_arrays_no_radius_temp[i] for i in range(len(columns_as_arrays_no_radius_temp)) if i not in index_of_removed_data]
#         # columns_as_arrays = [value for index, value in enumerate(columns_as_arrays) if index not in (index_of_removed_data+1)]
#         # old= [value for index, value in enumerate(old_i) if index not in index_of_removed_data]
#         # indices = [1,2,3,-2]
#         indices = [0,2,3,-1]
#         ip_list=[0,1,2,-1]
#         columns_as_arrays_reverted=[columns_as_arrays_no_radius_temp[ip_list[indices.index(i)]]*(np.cos(np.radians(old_i[i]))/np.cos(np.radians(new_i))) for i in indices]

# else:
#         indices = [0,1,2,-1]
#         columns_as_arrays_reverted=[columns_as_arrays_no_radius_temp[i]*np.cos(np.radians(old_i[i]))/np.cos(np.radians(new_i)) for i in indices]

# # for i in range(len(columns_as_arrays_reverted)):
# #         print('############################################')
# #         print(i,columns_as_arrays_reverted[i])
# #         print('############################################')

# # print(columns_as_arrays[-2])
# #columns as arrays reverted contains \sigma_tot, \sigma_HI, \sigma_H2, \sigma_SFR

# # print('e',columns_as_arrays_reverted)
# # print(len(columns_as_arrays_reverted[-1]),len(columns_as_arrays_reverted))
# #################################################################################################################

# #################################################################################################################
# #new version vcirc_q_omega calculator that returns the removed vcirc col also

# #need to change in data_helpers.py
# def vcirc_to_qomega_new(df,remove_vcirc=False):
#     df = df.copy()
#     vcirc_data = keep_substring_columns(df, 'vcirc')
#     if vcirc_data[0].empty:
#         return df
#     else:
#         r = df[get_adjacent_column(df, vcirc_data[1][0], False)].to_numpy().flatten()
#         Om = vcirc_data[0].to_numpy().flatten()/r
#         q = -1 * r/Om* np.gradient(Om)/np.gradient(r)
#         index_of_vcirc = df.columns.get_loc(vcirc_data[1][0])
#         df.insert(index_of_vcirc, '\Omega', Om)
#         df.insert(index_of_vcirc, 'r omega', r)
#         df.insert(index_of_vcirc, 'q', q)
#         if remove_vcirc:
#                 #get df without vcirc
#                 df.drop(columns=vcirc_data[1], inplace=True)
#         return df
# #################################################################################################################

# #################################################################################################################

# raw_data= vcirc_to_qomega_new(raw_data) #raw_data contain v_circ also
# raw_data_radii_df = keep_substring_columns(raw_data, 'r ')[0]
# # print(raw_data.iloc[:,-1][3]*g_Msun/((s_Myr*1e3)*(cm_pc**2)))


# err_data = incl_distance_correction(err_data, distance_new=new_dist, distance_old=old_dist,\
#             i_new=np.radians(new_i), i_old=np.radians(old_i)) #correcting errors for D and i
# err_radii_df = keep_substring_columns(err_data, 'R')[0]

# #making kpc_r and err_radii_df same length by adding 'nan' to kpc_r 
# #otherwise interpolation cant work
# list_kpcr=list(kpc_r)
# if len(kpc_r)>len(err_radii_df.iloc[:,0]):
#         print('code wont work')
# for i in range(len(err_radii_df.iloc[:,0])):
#         if len(kpc_r)<=i:
#                 list_kpcr.append(np.nan)
# kpc_r_long=np.array(list_kpcr)

#  #commented on 25/02/24 
# # list_kpcr=list(kpc_r)
# # if galaxy_name=='m33':
# #         h=list(raw_data.columns).index('r kpc.4')
# # elif galaxy_name=='m31':
# #         h=list(raw_data.columns).index('r kpc.2')
# # else:
# #         h=list(raw_data.columns).index('r kpc.3')

# #making length of kpc_r same as the largest element in raw_data

#  #commented on 25/02/24 
# # array_with_max_elements = min(columns_as_arrays, key=len)

# #cant figure out how this works
# array_with_min_elements = min(columns_as_arrays, key=len)
# y = columns_as_arrays.index(array_with_min_elements)

# for i in range(len(raw_data.iloc[:,y])):
#         if len(kpc_r)<=i:
#                 list_kpcr.append(np.nan)
# kpc_r_long2=np.array(list_kpcr)
# #################################################################################################################

# #################################################################################################################
# #interpolation happens here
# err_interpolated_df = df_interpolation(err_data,err_radii_df, kpc_r_long)
# # for i in range(len(raw_data)):
# #         print(i,len(raw_data.iloc[:,i-1]))
# print(raw_data)
# raw_data = df_interpolation(raw_data,raw_data_radii_df,kpc_r_long2)

# if galaxy_name=='m31':
#         vcirc_interpolated=raw_data['vcirc_claude kms']
# elif galaxy_name=='m33':
#         vcirc_interpolated=raw_data['kms_vcirc_Kam']
# else:
#         vcirc_interpolated=raw_data['vcirc']
# #################################################################################################################

# #################################################################################################################
# #remove nans
# nan_mask = np.isnan(raw_data)
# raw_data = raw_data[~(nan_mask.sum(axis=1)>0)]

# nan_mask = np.isnan(err_interpolated_df)
# err_interpolated_df = err_interpolated_df[~(nan_mask.sum(axis=1)>0)]
# #################################################################################################################

# #################################################################################################################
# #conversion of units for errors

# #distance unit conversion already done during dist and i correction
# if galaxy_name=='m31' or galaxy_name=='m33':
#     conv_factors=np.array([1,cm_km,1])
# else:
#     conv_factors=np.array([1,cm_km])

# err_interpolated_df = err_interpolated_df*conv_factors
# # switch = parameter_read(os.path.join(base_path,'inputs','switches.in'))
# #################################################################################################################

# #################################################################################################################
# #functions for error calculations start here
# def err_omega(omega_interpolated,vcirc,kpc_r,sigma_v,oldD,oldi,newdist,newi): #omega=columns_as_arrays[4]
#         oldD[4]=oldD[4]*cm_kpc*(10**3)
#         newdist=newdist*cm_kpc*(10**3)
        
#         oldi=np.radians(oldi)
#         newi=np.radians(newi)

#         if galaxy_name=='m31' or galaxy_name=='m33':
#                 u=4
#         else:
#                 u=3
#         # vcirc=[15*cm_km for i in range(len(kpc_r))]
#         vcirc=vcirc*cm_km
#         # print('vcirc',vcirc)
#         cm_r=kpc_r*cm_kpc
#         cm_r_NoDistCorr=cm_r #undoing the distance correction
#         sigma_v=[sigma_v[i]/np.sin(newi) for i in range(len(kpc_r))]
#         term1 = [(oldD[u]*np.sin(oldi[u])*sigma_v[i])/(newdist*cm_r_NoDistCorr[i]*np.sin(newi)) for i in range(len(kpc_r))]

#         term2 = [(vcirc[i]*oldD[u]*np.cos(oldi[u])*np.radians(err_old_i[u]))/((np.sin(newi))*cm_r_NoDistCorr[i]*newdist) for i in range(len(kpc_r))]       
#         # print('yyyyyyyyyyyyyyyyyyyyyyyy',cm_r_NoDistCorr)
#         # print('2',term2)
#         term3 = [(vcirc[i]*oldD[u]*np.sin(oldi[u])*np.cos(newi)*np.radians(err_new_i))/(((np.sin(newi))**2)*cm_r_NoDistCorr[i]*newdist) for i in range(len(kpc_r))]
#         # print('3',term3)
#         term4 = [(vcirc[i]*np.sin(oldi[u])*err_old_dist[u]*cm_kpc*(10**3))/((np.sin(newi))*cm_r_NoDistCorr[i]*newdist) for i in range(len(kpc_r))]
#         # print('4',term4)
#         term5 = [(vcirc[i]*np.sin(oldi[u])*oldD[u]*err_new_dist*cm_kpc*(10**3))/((np.sin(newi))*cm_r_NoDistCorr[i]*(newdist**2)) for i in range(len(kpc_r))]
#         # print('5',term5)
#         err_omega     = [np.sqrt(term1[i]**2 +term2[i]**2 + term3[i]**2 + term4[i]**2 + term5[i]**2) for i in range(len(kpc_r))]
#         rel_err_omega = err_omega/omega_interpolated
#         print('omega err',rel_err_omega)
#         return rel_err_omega


# def err_sigmaHI(sigmaHI_corrected,sigmaHI_not_corrected,percent_sigmaHI_err,newi,oldi): #sigmaHI_err is set to 6 percent, sigmaHI=columns_as_arrays[2]
#         # print('sigmahi',sigmaHI)
#         newi=np.radians(newi)
#         oldi=np.radians(oldi)

#         if galaxy_name=='m31' or galaxy_name=='m33':
#                 u=2
#         else:
#                 u=1
#         # sigmaHI_not_corrected=(g_Msun/(cm_pc**2) )*sigmaHI_not_corrected
#         sigmaHI_err=sigmaHI_not_corrected*percent_sigmaHI_err
#         term1 = [(sigmaHI_err[i]*np.cos(newi))/np.cos(oldi[u]) for i in range(len(kpc_r))]
#         # print('1',term1)
#         term2 = [(sigmaHI_not_corrected[i]*np.sin(newi)*np.radians(err_new_i))/np.cos(oldi[u]) for i in range(len(kpc_r))]
#         # print('2',term2)
#         term3 = [(sigmaHI_not_corrected[i]*np.cos(newi)*np.radians(err_old_i[u])*np.sin(oldi[u]))/((np.cos(oldi[2]))**2) for i in range(len(kpc_r))]
#         # print('3',term3)
#         err_sigmaHI= [np.sqrt(term1[i]**2 + term2[i]**2 + term3[i]**2) for i in range(len(kpc_r))]
#         # print('err',err_sigmaHI)
#         rel_err_sigmaHI = [err_sigmaHI[i]/sigmaHI_corrected[i] for i in range(len(kpc_r))]
#         print('rel_errHI',rel_err_sigmaHI)
#         print('errHI',err_sigmaHI)
        
#         return err_sigmaHI,rel_err_sigmaHI

# # rel_err_sigmaHI=err_sigmaHI(columns_as_arrays[2],raw_data['sigma_HI_claude'],0.06,new_i,old_i)
# # print(rel_err_sigmaHI)

# def err_sigmaH2(sigmaHI_not_corrected,molfrac_or_sigmaH2,percent_sigmaHI_err,molfrac_err,newi,oldi,percent_sigmaH2_err=0):
#         newi=np.radians(newi)
#         oldi=np.radians(oldi)
#         # sigmaHI_not_corrected=(g_Msun/(cm_pc**2) )*sigmaHI_not_corrected
#         if switch['incl_moldat']=='No':
#                 print('considering ONLY HI surface density')
#                 err_sigmaH2=0
#                 return err_sigmaH2
#         else:
#                 if galaxy_name=='m31':
#                         molfrac=molfrac_or_sigmaH2
#                         sigmaH2 = [(molfrac[i]/(1-molfrac[i]))*sigmaHI_not_corrected[i] for i in range(len(kpc_r))]
#                         molfrac_err_tot = np.sqrt(((raw_data['molfrac'][i]*np.sin(newi)*np.radians(err_new_i))/(np.cos(oldi[-1])))**2 + ((raw_data['molfrac'][i]*np.cos(newi)*np.radians(err_old_i[-1])*np.sin(oldi[-1]))/(np.cos(oldi[-1]))**2)**2 + (molfrac_err[i]*np.cos(newi))/np.cos(oldi[-1])**2)
#                         sigmaHI_err=percent_sigmaHI_err*sigmaHI_not_corrected
#                         term1= [(sigmaHI_err[i]*molfrac[i])/(1-molfrac[i]) for i in range(len(kpc_r))]
#                         term2= [(sigmaHI_not_corrected[i]*molfrac_err_tot[i])/((1-molfrac[i])**2) for i in range(len(kpc_r))]
#                         err_sigmaH2 = np.sqrt(term1**2+ term2**2)
#                         rel_err_sigmaH2 = [err_sigmaH2[i]/sigmaH2]
#                         # print('sigmah2err',rel_err_sigmaH2)
#                         return err_sigmaH2
#                 else:
#                         sigmaH2_not_corrected=molfrac_or_sigmaH2
#                         percent_sigmaH2_err=0.06
#                         sigmaH2_err=sigmaH2_not_corrected*percent_sigmaH2_err
#                         index_for_sigmaH2=((list(columns_as_arrays)).index('sigma_H2'))-1
#                         if galaxy_name=='m33':
#                                 index_for_sigmaH2+=1
#                         term1 = [(sigmaH2_err[i]*np.cos(newi))/np.cos(oldi[index_for_sigmaH2]) for i in range(len(kpc_r))]
#                         # print('1',term1)
#                         term2 = [(sigmaH2_not_corrected[i]*np.sin(newi)*np.radians(err_new_i))/np.cos(oldi[index_for_sigmaH2]) for i in range(len(kpc_r))]
#                         # print('2',term2)
#                         term3 = [(sigmaH2_not_corrected[i]*np.cos(newi)*np.radians(err_old_i[index_for_sigmaH2])*np.sin(oldi[index_for_sigmaH2]))/(np.cos(oldi[index_for_sigmaH2]))**2 for i in range(len(kpc_r))]
#                         # print('3',term3)
#                         err_sigmaH2= [np.sqrt(term1[i]**2 + term2[i]**2 + term3[i]**2) for i in range(len(kpc_r))]
#                         rel_err_sigmaH2 = [err_sigmaH2[i]/sigmaH2_not_corrected]
#                         # print('sigmah2err',rel_err_sigmaH2)
#                         return err_sigmaH2
# # err_sigmaH2=err_sigmaH2()

# def err_sigmaSFR(sigmaSFR_corrected,sigmaSFR_not_corrected,percent_sigmaSFR_err,newi,oldi):
#         newi=np.radians(newi)
#         oldi=np.radians(oldi)

#         # sigmaSFR_not_corrected=g_Msun/((s_Myr*m_to_gconv)*(cm_pc**2))*sigmaSFR_not_corrected
#         sigmaSFR_not_corrected=sigmaSFR_not_corrected
#         sigmaSFR_err=sigmaSFR_not_corrected*percent_sigmaSFR_err
#         if galaxy_name=='m31':
#                 u=-2
#         else:
#                 u=-1
#         term1 = [(sigmaSFR_err[i]*np.cos(newi))/np.cos(oldi[u]) for i in range(len(kpc_r))]
#         # print('1',term1)
#         term2 = [(sigmaSFR_not_corrected[i]*np.sin(newi)*np.radians(err_new_i))/np.cos(oldi[u]) for i in range(len(kpc_r))]
#         # print('2',term2)
#         term3 = [(sigmaSFR_not_corrected[i]*np.cos(newi)*np.radians(err_old_i[u])*np.sin(oldi[u]))/(np.cos(oldi[u]))**2 for i in range(len(kpc_r))]
#         # print('3',term3)
#         err_sigmaSFR= [np.sqrt(term1[i]**2 + term2[i]**2 + term3[i]**2) for i in range(len(kpc_r))]
#         # print('err',err_sigmaSFR)
#         rel_err_sigmaSFR = [err_sigmaSFR[i]/sigmaSFR_corrected[i] for i in range(len(kpc_r))]
#         print('relerrSFR',rel_err_sigmaSFR)
#         return rel_err_sigmaSFR
#         # err_sigmaSFR= np.sqrt(((sigmaSFR_err*np.cos(np.radians(i))/np.cos(np.radians(i_old)))**2 + (sigmaSFR*np.sin(np.radians(i))*err_new_i/np.cos(np.radians(i_old)))**2))
#         # return err_sigmaSFR

# # print(rel_err_sigmaSFR)

# # def err_sigmagas(sigmaHI_corr,sigmaH2_corr,sigmaHI_no_corr,mu_no_corr,err_mu=0):
# #         # print('#################################################')
# #         # print(sigmaHI_corr)
# #         # print(sigmaH2_corr)
# #         # print(sigmaHI_no_corr)
# #         # sigmaHI_no_corr=(g_Msun/(cm_pc**2) )*sigmaHI_no_corr

# #         #change this. its not mu but molfrac thats needed in the equation
# #         sigmaH2_no_corr=sigmaHI_no_corr*(mu_no_corr/(1-mu_no_corr))
# #         #to get corrected sigmagas
# #         if switch['incl_moldat']=='Yes':
# #                 sigmagas_corr=sigmaHI_corr+ sigmaH2_corr #columns_as_arrays[2]+columns_as_arrays[3]
# #                 sigmagas_no_corr=sigmaHI_no_corr+sigmaH2_no_corr
# #         else:
# #                 sigmagas_corr=sigmaHI_corr
# #                 sigmagas_no_corr=sigmaHI_no_corr
# #         mu_err=[err_mu*mu_no_corr[i] for i in range(len(kpc_r))]
# #         sigmaHI_err=err_sigmaHI(sigmaHI_corr,sigmaHI_no_corr,0.06,new_i,old_i)[0]
# #         #might have to changed in case of moldat==yes
# #         sigmaH2_err=[0 for i in range(len(kpc_r))]
# #         sigmagas_err=[np.sqrt(sigmaHI_err[i]**2 + sigmaH2_err[i]**2) for i in range(len(kpc_r))]
# #         # print('########################################################################')
# #         print('step1',sigmagas_err)
# #         term1= [(sigmagas_err[i]*(3*mu_no_corr[i]/(4-mu_no_corr[i]))) for i in range(len(kpc_r))]
# #         print('term1',term1)
# #         term2= [(sigmagas_no_corr[i]*(12*mu_err[i]/(4-mu_no_corr[i])**2)) for i in range(len(kpc_r))]
# #         # print('term2',term2)
# #         err_sigmagas= [np.sqrt((term1[i]**2 + term2[i]**2)) for i in range(len(kpc_r))]
# #         print('err',err_sigmagas)
# #         rel_err_sigmagas= [err_sigmagas[i]/sigmagas_corr[i] for i in range(len(kpc_r))]
# #         print('rel_err_sigmagas',rel_err_sigmagas)
# #         return rel_err_sigmagas

# def err_sigmagas(sigmaHI_corr,sigmaH2_corr,sigmaHI_no_corr,sigmaH2_no_corr,err_mu=0):
#         # print('#################################################')
#         # print(sigmaHI_corr)
#         # print(sigmaH2_corr)
#         # print(sigmaHI_no_corr)
#         # sigmaHI_no_corr=(g_Msun/(cm_pc**2) )*sigmaHI_no_corr
#         mu= params['mu']
#         #change this. its not mu but molfrac thats needed in the equation
#         # sigmaH2_no_corr=sigmaHI_no_corr*(mu_no_corr/(1-mu_no_corr))
#         #to get corrected sigmagas
#         if switch['incl_moldat']=='Yes':
#                 sigmagas_corr=sigmaHI_corr+ sigmaH2_corr #columns_as_arrays[2]+columns_as_arrays[3]
#                 sigmagas_no_corr=sigmaHI_no_corr+sigmaH2_no_corr
#                 sigmaH2_err=err_sigmaH2(sigmaH2_no_corr,sigmaH2_corr,0.06,)
#         else:
#                 sigmagas_corr=sigmaHI_corr
#                 sigmagas_no_corr=sigmaHI_no_corr
#                 sigmaH2_err=[0 for i in range(len(kpc_r))]

#         mu_err=[err_mu*mu for i in range(len(kpc_r))]
#         sigmaHI_err=err_sigmaHI(sigmaHI_corr,sigmaHI_no_corr,0.06,new_i,old_i)[0]
#         #might have to changed in case of moldat==yes
#         sigmagas_err=[np.sqrt(sigmaHI_err[i]**2 + sigmaH2_err[i]**2) for i in range(len(kpc_r))]
#         # print('########################################################################')
#         # print('step1',sigmagas_err)
#         term1= [(sigmagas_err[i]*(3*mu/(4-mu))) for i in range(len(kpc_r))]
#         # print('term1',term1)
#         term2= [(sigmagas_no_corr[i]*(12*mu_err[i]/(4-mu)**2)) for i in range(len(kpc_r))]
#         # print('term2',term2)
#         err_sigmagas= [np.sqrt((term1[i]**2 + term2[i]**2)) for i in range(len(kpc_r))]
#         # print('err',err_sigmagas)
#         rel_err_sigmagas= [err_sigmagas[i]/sigmagas_corr[i] for i in range(len(kpc_r))]
#         print('rel_err_sigmagas',rel_err_sigmagas)
#         return rel_err_sigmagas

# # print(rel_err_sigmagas)

# def err_sigmatot(galaxy_name,sigmatot_corrected,sigmatot_not_corrected,percent_sigmatot_err,newi,oldi,new_dist=0,old_dist=0):
#         newi=np.radians(newi)
#         oldi=np.radians(oldi)
#         # sigmatot_not_corrected=(g_Msun/(cm_pc**2))*sigmatot_not_corrected
#         if galaxy_name=='m33':
#                 gamma=0.52*g_Msun/((cm_pc)**2)
#                 err_gamma=0.1*g_Msun/((cm_pc)**2)
#                 mu0=18.01
#                 kpc_Rd=1.82*cm_kpc 
#                 c36=24.8
#                 new_dist=new_dist*cm_kpc*(10**3)
#                 old_dist=old_dist*cm_kpc*(10**3)
#                 cm_r=list(cm_kpc*(raw_data.iloc[:,0]))
                
#                 mu=[(mu0+(1.10857*(cm_r[i]/kpc_Rd))) for i in range(len(kpc_r))]
#                 # print('mu',mu)
#                 sigma_mu=[((1.10857)/kpc_Rd)*cm_r[i]*np.sqrt(((err_new_dist*cm_kpc*(10**3))/(new_dist))**2 + ((new_dist*err_old_dist[0]*cm_kpc*(10**3))/(old_dist[0])**2)**2) for i in range(len(kpc_r))]
#                 # print('sigmamu',sigma_mu)
#                 a=[10**(-0.4*(mu[i]-c36)) for i in range(len(kpc_r))]
#                 # print('a',a)
#                 da_dt=[(np.log(10))*a[i] for i in range(len(kpc_r))]
#                 # print('dadt',da_dt)
#                 sigma_a= [-0.4*da_dt[i]*sigma_mu[i] for i in range(len(kpc_r))]
#                 # print('sigmaa',sigma_a)
#                 term1= [err_gamma*a[i] for i in range(len(kpc_r))]
#                 # print('term1',term1)
#                 term2= [gamma*sigma_a[i] for i in range(len(kpc_r))]
#                 # print('term2',term2)
#                 sigmatotal=[gamma*a[i] for i in range(len(kpc_r))]
#                 # print('totcalc', sigmatotal)
#                 # print('actual',sigmatot_corrected)
#                 err_sigmatot=[np.sqrt(term1[i]**2 + term2[i]**2) for i in range(len(kpc_r))]
#                 # print('err',err_sigmatot)
#                 # rel_err_sigmatot = [err_sigmatot[i]/sigmatot_not_corrected[i] for in range(len(kpc_r))]
                
#         else:
#                 sigmatot_err=percent_sigmatot_err*sigmatot_not_corrected
#                 term1=[(sigmatot_err[i]*np.cos(newi)/np.cos(oldi[0])) for i in range(len(kpc_r))]
#                 # print('1',term1)
#                 term2=[(sigmatot_not_corrected[i]*np.sin(newi)*np.radians(err_new_i)/np.cos(oldi[0])) for i in range(len(kpc_r))]
#                 # print('2',term2)
#                 term3=[(sigmatot_not_corrected[i]*np.cos(newi)*np.radians(err_old_i[0])*np.sin(oldi[0]))/((np.cos(oldi[0]))**2) for i in range(len(kpc_r))]
#                 # print('3',term3)
#                 err_sigmatot= [np.sqrt(term1[i]**2 + term2[i]**2 + term3[i]**2) for i in range(len(kpc_r))]
        
#         rel_err_sigmatot = [err_sigmatot[i]/sigmatot_corrected[i] for i in range(len(kpc_r))]
#         # print('rel err tot',rel_err_sigmatot)
#         # rel_err_sigmatot=[rel_err_sigmatot[i]*100 for i in range(len(kpc_r))]
#         print('rel_err_sigmatot',rel_err_sigmatot)
#         return rel_err_sigmatot
#         #do the other derivation with upsilon for m33

# # print(rel_err_sigmatot)

# def err_T_q(T_corrected,q_corrected):
#         # print(T_corrected,q_corrected)
#         err_T=np.std(T_corrected)
#         err_q=np.std(q_corrected)
#         # print('*******************',err_T,err_q)
#         rel_err_T=[err_T/T_corrected[i] for i in range(len(kpc_r))]
#         rel_err_q=[err_q/q_corrected[i] for i in range(len(kpc_r))]
#         print('err_T',rel_err_T)
#         print('err_q',rel_err_q)
#         return rel_err_T,rel_err_q
# #################################################################################################################

# #################################################################################################################
# #calculating relative errors 
# hregs = ['subsonic', 'supersonic']
# for hreg in hregs:
#         os.chdir(os.path.join(base_path,'inputs'))
#         exps = np.load(f'scal_exponents_{hreg}.npy')
#         r = kpc_r.size
#         if hreg=='supersonic':
#                 rel_err_omega =err_omega(columns_as_arrays[5],vcirc_interpolated,kpc_r,err_interpolated_df.iloc[:,1],old_dist/(cm_kpc*(10**3)),old_i,new_dist,new_i)
#         else:
#                 rel_err_omega =err_omega(columns_as_arrays[5],vcirc_interpolated,kpc_r,err_interpolated_df.iloc[:,1],old_dist,old_i,new_dist,new_i)
        
#         if galaxy_name=='m33':
#                 rel_err_sigmatot=err_sigmatot(galaxy_name,columns_as_arrays[1],columns_as_arrays_reverted[0],0.1,new_i,old_i,new_dist,old_dist)
#         else:
#                 rel_err_sigmatot=err_sigmatot(galaxy_name,columns_as_arrays[1],columns_as_arrays_reverted[0],0.1,new_i,old_i,new_dist,old_dist)

#         if galaxy_name=='m31':
#                 rel_err_sigmagas=err_sigmagas(columns_as_arrays[2],columns_as_arrays[3],columns_as_arrays_reverted[1],columns_as_arrays_reverted[2],0)
#         else:
#                 rel_err_sigmagas=err_sigmagas(columns_as_arrays[2],columns_as_arrays[3],columns_as_arrays_reverted[1],columns_as_arrays_reverted[2],0)

#         rel_err_sigmasfr=err_sigmaSFR(columns_as_arrays[-2],columns_as_arrays_reverted[-1],0.1,new_i,old_i)
#         rel_err_T,rel_err_q=err_T_q(columns_as_arrays[-1],columns_as_arrays[4])
#         rel_err = np.array([rel_err_q, rel_err_omega, rel_err_sigmagas, rel_err_sigmatot,rel_err_sigmasfr, rel_err_T])

#         relerr_quan = np.sqrt(np.matmul(exps**2,rel_err**2))
#         err_quantities = model_f[1:]*relerr_quan
#         # print('P',rel_err_omega)
# #################################################################################################################
# #################################################################################################################
# #inputting errors into a pickle file
#         os.chdir(os.path.join(base_path,'outputs'))
#         with open(f'errors_{hreg}.out', 'wb') as f:
#                 pickle.dump(err_quantities, f)
# print('Found the errors from the scaling relations')
 
import numpy as np
import pickle
import os
import pandas as pd
from helper_functions import parameter_read
from data_helpers import *


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

#########################################################################################
#error in distance and inclination

base_path = os.environ.get('MY_PATH')
galaxy_name = os.environ.get('galaxy_name')

#get parameter values
params = parameter_read(os.path.join(base_path,'inputs','parameter_file.in'))
switch = parameter_read(os.path.join(base_path,'inputs','switches.in'))

current_directory = str(os.getcwd())
os.chdir(os.path.join(base_path,'outputs'))

with open(f'{galaxy_name}output_ca_'+str(params[r'C_\alpha'])+'K_'+str(params[r'K'])+'z_'+
          str(params[r'\zeta'])+'psi_'+str(params[r'\psi'])+'b_'+str(params[r'\beta'])+'.out', 'rb') as f:
                model_f = pickle.load(f)

os.chdir(os.path.join(base_path,'inputs'))

#getting available data and error
os.chdir(os.path.join(base_path, 'data','model_data', f'{galaxy_name}_data'))
raw_data = pd.read_csv(f'combined_data_{galaxy_name}.csv', skiprows=1) #needed to obtain vcirc values which isnt present in interpolated_data.csv
err_data= pd.read_csv(f'error_combined_{galaxy_name}.csv')
corrections = pd.read_csv(f'correction_data_{galaxy_name}.csv', skiprows=1, index_col=0)

print(galaxy_name)
os.chdir(os.path.join(base_path,'data'))

# kpc_r, dat_sigmatot, dat_sigmaHI, dat_sigmaH2, dat_q, dat_omega, dat_sigmasfr, T= data
data_frame = pd.read_csv("data_interpolated.csv")
columns_as_arrays = [np.array(data_frame[col]) for col in data_frame.columns]
kpc_r = np.array(data_frame.iloc[:, 0])

#distance correction
new_dist= corrections.iloc[-1,0]
err_new_dist=corrections.iloc[-1,2]
old_dist = corrections.iloc[:-1,0].values
err_old_dist=corrections.iloc[:-1,2].values

#inclination correction
new_i= corrections.iloc[-1,1] #deg
err_new_i= corrections.iloc[-1,3] #deg
old_i= corrections.iloc[:-1,1].values #used new_i as no inclination correction is needed for Claude data
err_old_i=corrections.iloc[:-1,3].values

#################################################################################################################
#raw_data from combined_data.csv is converted to correct units, but no D/i correction done here
raw_data = incl_distance_correction(raw_data, distance_new=old_dist, distance_old=old_dist,\
                        i_new=np.radians(old_i), i_old=np.radians(old_i))

os.chdir(os.path.join(base_path, 'data','model_data', f'{galaxy_name}_data'))

#unwanted cols are removed from raw_data
index_of_removed_data=[]
try:
    data_rem = pd.read_csv(f'removed_data_{galaxy_name}.csv', dtype=str)
except:
    data_rem = []
for d in data_rem:
    #change this
    names=[i for i in raw_data.columns if d in i]
    index_of_removed_data=[raw_data.columns.get_loc(i) for i in names]
    remove_data(raw_data, d)
#################################################################################################################

#################################################################################################################
columns_as_arrays_no_radius_temp = columns_as_arrays[1:-1]
#reverting the correction for interpolated data to be used for error calculation


if galaxy_name=='m31':
        ip_list=[0,1,2,-1]
        if switch['chem_or_claude']=='Claude':
                indices = [0,2,-1,-2] #these are indices for accessing old_inc values
        else:
                indices = [0,1,-1,-2]
        columns_as_arrays_reverted=[columns_as_arrays_no_radius_temp[ip_list[indices.index(i)]]*(np.cos(np.radians(old_i[i]))/np.cos(np.radians(new_i))) for i in indices]

elif galaxy_name=='m33':
        columns_as_arrays = [columns_as_arrays[i] for i in range(len(columns_as_arrays)) if i not in index_of_removed_data]
        # columns_as_arrays_no_radius_temp = [columns_as_arrays_no_radius_temp[i] for i in range(len(columns_as_arrays_no_radius_temp)) if i not in index_of_removed_data]
        # columns_as_arrays = [value for index, value in enumerate(columns_as_arrays) if index not in (index_of_removed_data+1)]
        # old= [value for index, value in enumerate(old_i) if index not in index_of_removed_data]
        # indices = [1,2,3,-2]
        indices = [0,2,3,-1]
        ip_list=[0,1,2,-1]
        columns_as_arrays_reverted=[columns_as_arrays_no_radius_temp[ip_list[indices.index(i)]]*(np.cos(np.radians(old_i[i]))/np.cos(np.radians(new_i))) for i in indices]

else:
        indices = [0,1,2,-1]
        columns_as_arrays_reverted=[columns_as_arrays_no_radius_temp[i]*np.cos(np.radians(old_i[i]))/np.cos(np.radians(new_i)) for i in indices]

# for i in range(len(columns_as_arrays_reverted)):
#         print('############################################')
#         print(i,columns_as_arrays_reverted[i])
#         print('############################################')

# print(columns_as_arrays[-2])
#columns as arrays reverted contains \sigma_tot, \sigma_HI, \sigma_H2, \sigma_SFR

# print('e',columns_as_arrays_reverted)
# print(len(columns_as_arrays_reverted[-1]),len(columns_as_arrays_reverted))
#################################################################################################################

#################################################################################################################
#new version vcirc_q_omega calculator that returns the removed vcirc col also

#need to change in data_helpers.py
def vcirc_to_qomega_new(df,remove_vcirc=False):
    df = df.copy()
    vcirc_data = keep_substring_columns(df, 'vcirc')
    if vcirc_data[0].empty:
        return df
    else:
        r = df[get_adjacent_column(df, vcirc_data[1][0], False)].to_numpy().flatten()
        Om = vcirc_data[0].to_numpy().flatten()/r
        q = -1 * r/Om* np.gradient(Om)/np.gradient(r)
        index_of_vcirc = df.columns.get_loc(vcirc_data[1][0])
        df.insert(index_of_vcirc, '\Omega', Om)
        df.insert(index_of_vcirc, 'r omega', r)
        df.insert(index_of_vcirc, 'q', q)
        if remove_vcirc:
                #get df without vcirc
                df.drop(columns=vcirc_data[1], inplace=True)
        return df
#################################################################################################################

#################################################################################################################

raw_data= vcirc_to_qomega_new(raw_data) #raw_data contain v_circ also
raw_data_radii_df = keep_substring_columns(raw_data, 'r ')[0]
# print(raw_data.iloc[:,-1][3]*g_Msun/((s_Myr*1e3)*(cm_pc**2)))


err_data = incl_distance_correction(err_data, distance_new=new_dist, distance_old=old_dist,\
            i_new=np.radians(new_i), i_old=np.radians(old_i)) #correcting errors for D and i
err_radii_df = keep_substring_columns(err_data, 'R')[0]

#making kpc_r and err_radii_df same length by adding 'nan' to kpc_r 
#otherwise interpolation cant work
list_kpcr=list(kpc_r)
if len(kpc_r)>len(err_radii_df.iloc[:,0]):
        print('code wont work')
for i in range(len(err_radii_df.iloc[:,0])):
        if len(list_kpcr)<=i:
                list_kpcr.append(np.nan)
kpc_r_long=np.array(list_kpcr)

 #commented on 25/02/24 
# list_kpcr=list(kpc_r)
# if galaxy_name=='m33':
#         h=list(raw_data.columns).index('r kpc.4')
# elif galaxy_name=='m31':
#         h=list(raw_data.columns).index('r kpc.2')
# else:
#         h=list(raw_data.columns).index('r kpc.3')

#making length of kpc_r same as the largest element in raw_data

 #commented on 25/02/24 
# array_with_max_elements = min(columns_as_arrays, key=len)

#cant figure out how this works
array_with_min_elements = min(columns_as_arrays, key=len)
y = columns_as_arrays.index(array_with_min_elements)
# print(columns_as_arrays)
for i in range(len(raw_data.iloc[:,y])):
        if len(list_kpcr)<=i:
                list_kpcr.append(np.nan)
kpc_r_long2=np.array(list_kpcr)
#################################################################################################################

#################################################################################################################
#interpolation happens here
err_interpolated_df = df_interpolation(err_data,err_radii_df, kpc_r_long)
# for i in range(len(raw_data)):
#         print(i,len(raw_data.iloc[:,i-1]))

raw_data = df_interpolation(raw_data,raw_data_radii_df,kpc_r_long2)

if galaxy_name=='m31':
        vcirc_interpolated=raw_data['vcirc_claude kms']
elif galaxy_name=='m33':
        vcirc_interpolated=raw_data['kms_vcirc_Kam']
else:
        vcirc_interpolated=raw_data['vcirc']
#################################################################################################################

#################################################################################################################
#remove nans
nan_mask = np.isnan(raw_data)
raw_data = raw_data[~(nan_mask.sum(axis=1)>0)]

nan_mask = np.isnan(err_interpolated_df)
err_interpolated_df = err_interpolated_df[~(nan_mask.sum(axis=1)>0)]
#################################################################################################################

#################################################################################################################
#conversion of units for errors

#distance unit conversion already done during dist and i correction
if galaxy_name=='m31' or galaxy_name=='m33':
    conv_factors=np.array([1,cm_km,1])
else:
    conv_factors=np.array([1,cm_km])

err_interpolated_df = err_interpolated_df*conv_factors
switch = parameter_read(os.path.join(base_path,'inputs','switches.in'))
#################################################################################################################

#################################################################################################################
#functions for error calculations start here
def err_omega(omega_interpolated,vcirc,kpc_r,sigma_v,oldD,oldi,newdist,newi): #omega=columns_as_arrays[4]
        oldD[4]=oldD[4]*cm_kpc*(10**3)
        newdist=newdist*cm_kpc*(10**3)
        
        oldi=np.radians(oldi)
        newi=np.radians(newi)

        if galaxy_name=='m31' or galaxy_name=='m33':
                u=4
        else:
                u=3
        # vcirc=[15*cm_km for i in range(len(kpc_r))]
        vcirc=vcirc*cm_km
        # print('vcirc',vcirc)
        cm_r=kpc_r*cm_kpc
        cm_r_NoDistCorr=cm_r #undoing the distance correction
        sigma_v=[sigma_v[i]/np.sin(newi) for i in range(len(kpc_r))]
        term1 = [(oldD[u]*np.sin(oldi[u])*sigma_v[i])/(newdist*cm_r_NoDistCorr[i]*np.sin(newi)) for i in range(len(kpc_r))]

        term2 = [(vcirc[i]*oldD[u]*np.cos(oldi[u])*np.radians(err_old_i[u]))/((np.sin(newi))*cm_r_NoDistCorr[i]*newdist) for i in range(len(kpc_r))]       
        # print('yyyyyyyyyyyyyyyyyyyyyyyy',cm_r_NoDistCorr)
        # print('2',term2)
        term3 = [(vcirc[i]*oldD[u]*np.sin(oldi[u])*np.cos(newi)*np.radians(err_new_i))/(((np.sin(newi))**2)*cm_r_NoDistCorr[i]*newdist) for i in range(len(kpc_r))]
        # print('3',term3)
        term4 = [(vcirc[i]*np.sin(oldi[u])*err_old_dist[u]*cm_kpc*(10**3))/((np.sin(newi))*cm_r_NoDistCorr[i]*newdist) for i in range(len(kpc_r))]
        # print('4',term4)
        term5 = [(vcirc[i]*np.sin(oldi[u])*oldD[u]*err_new_dist*cm_kpc*(10**3))/((np.sin(newi))*cm_r_NoDistCorr[i]*(newdist**2)) for i in range(len(kpc_r))]
        # print('5',term5)
        err_omega     = [np.sqrt(term1[i]**2 +term2[i]**2 + term3[i]**2 + term4[i]**2 + term5[i]**2) for i in range(len(kpc_r))]
        rel_err_omega = err_omega/omega_interpolated
        print('omega err',rel_err_omega)
        return rel_err_omega


def err_sigmaHI(sigmaHI_corrected,sigmaHI_not_corrected,percent_sigmaHI_err,newi,oldi): #sigmaHI_err is set to 6 percent, sigmaHI=columns_as_arrays[2]
        # print('sigmahi',sigmaHI)
        newi=np.radians(newi)
        oldi=np.radians(oldi)

        if galaxy_name=='m31' or galaxy_name=='m33':
                u=2
        else:
                u=1
        # sigmaHI_not_corrected=(g_Msun/(cm_pc**2) )*sigmaHI_not_corrected
        sigmaHI_err=sigmaHI_not_corrected*percent_sigmaHI_err
        term1 = [(sigmaHI_err[i]*np.cos(newi))/np.cos(oldi[u]) for i in range(len(kpc_r))]
        # print('1',term1)
        term2 = [(sigmaHI_not_corrected[i]*np.sin(newi)*np.radians(err_new_i))/np.cos(oldi[u]) for i in range(len(kpc_r))]
        # print('2',term2)
        term3 = [(sigmaHI_not_corrected[i]*np.cos(newi)*np.radians(err_old_i[u])*np.sin(oldi[u]))/((np.cos(oldi[2]))**2) for i in range(len(kpc_r))]
        # print('3',term3)
        err_sigmaHI= [np.sqrt(term1[i]**2 + term2[i]**2 + term3[i]**2) for i in range(len(kpc_r))]
        # print('err',err_sigmaHI)
        rel_err_sigmaHI = [err_sigmaHI[i]/sigmaHI_corrected[i] for i in range(len(kpc_r))]
        print('rel_errHI',rel_err_sigmaHI)
        print('errHI',err_sigmaHI)
        
        return err_sigmaHI,rel_err_sigmaHI

# rel_err_sigmaHI=err_sigmaHI(columns_as_arrays[2],raw_data['sigma_HI_claude'],0.06,new_i,old_i)
# print(rel_err_sigmaHI)

def err_sigmaH2(sigmaHI_not_corrected,molfrac_or_sigmaH2,percent_sigmaHI_err,molfrac_err,newi,oldi,percent_sigmaH2_err=0):
        newi=np.radians(newi)
        oldi=np.radians(oldi)
        # sigmaHI_not_corrected=(g_Msun/(cm_pc**2) )*sigmaHI_not_corrected
        if switch['incl_moldat']=='No':
                print('considering ONLY HI surface density')
                err_sigmaH2=0
                return err_sigmaH2
        else:
                if galaxy_name=='m31':
                        molfrac=molfrac_or_sigmaH2
                        sigmaH2 = [(molfrac[i]/(1-molfrac[i]))*sigmaHI_not_corrected[i] for i in range(len(kpc_r))]
                        molfrac_err_tot = np.sqrt(((raw_data['molfrac'][i]*np.sin(newi)*np.radians(err_new_i))/(np.cos(oldi[-1])))**2 + ((raw_data['molfrac'][i]*np.cos(newi)*np.radians(err_old_i[-1])*np.sin(oldi[-1]))/(np.cos(oldi[-1]))**2)**2 + (molfrac_err[i]*np.cos(newi))/np.cos(oldi[-1])**2)
                        sigmaHI_err=percent_sigmaHI_err*sigmaHI_not_corrected
                        term1= [(sigmaHI_err[i]*molfrac[i])/(1-molfrac[i]) for i in range(len(kpc_r))]
                        term2= [(sigmaHI_not_corrected[i]*molfrac_err_tot[i])/((1-molfrac[i])**2) for i in range(len(kpc_r))]
                        err_sigmaH2 = np.sqrt(term1**2+ term2**2)
                        rel_err_sigmaH2 = [err_sigmaH2[i]/sigmaH2]
                        # print('sigmah2err',rel_err_sigmaH2)
                        return err_sigmaH2
                else:
                        sigmaH2_not_corrected=molfrac_or_sigmaH2
                        percent_sigmaH2_err=0.06
                        sigmaH2_err=sigmaH2_not_corrected*percent_sigmaH2_err
                        index_for_sigmaH2=((list(columns_as_arrays)).index('sigma_H2'))-1
                        if galaxy_name=='m33':
                                index_for_sigmaH2+=1
                        term1 = [(sigmaH2_err[i]*np.cos(newi))/np.cos(oldi[index_for_sigmaH2]) for i in range(len(kpc_r))]
                        # print('1',term1)
                        term2 = [(sigmaH2_not_corrected[i]*np.sin(newi)*np.radians(err_new_i))/np.cos(oldi[index_for_sigmaH2]) for i in range(len(kpc_r))]
                        # print('2',term2)
                        term3 = [(sigmaH2_not_corrected[i]*np.cos(newi)*np.radians(err_old_i[index_for_sigmaH2])*np.sin(oldi[index_for_sigmaH2]))/(np.cos(oldi[index_for_sigmaH2]))**2 for i in range(len(kpc_r))]
                        # print('3',term3)
                        err_sigmaH2= [np.sqrt(term1[i]**2 + term2[i]**2 + term3[i]**2) for i in range(len(kpc_r))]
                        rel_err_sigmaH2 = [err_sigmaH2[i]/sigmaH2_not_corrected]
                        # print('sigmah2err',rel_err_sigmaH2)
                        return err_sigmaH2
# err_sigmaH2=err_sigmaH2()

def err_sigmaSFR(sigmaSFR_corrected,sigmaSFR_not_corrected,percent_sigmaSFR_err,newi,oldi):
        newi=np.radians(newi)
        oldi=np.radians(oldi)

        # sigmaSFR_not_corrected=g_Msun/((s_Myr*m_to_gconv)*(cm_pc**2))*sigmaSFR_not_corrected
        sigmaSFR_not_corrected=sigmaSFR_not_corrected
        sigmaSFR_err=sigmaSFR_not_corrected*percent_sigmaSFR_err
        if galaxy_name=='m31':
                u=-2
        else:
                u=-1
        term1 = [(sigmaSFR_err[i]*np.cos(newi))/np.cos(oldi[u]) for i in range(len(kpc_r))]
        # print('1',term1)
        term2 = [(sigmaSFR_not_corrected[i]*np.sin(newi)*np.radians(err_new_i))/np.cos(oldi[u]) for i in range(len(kpc_r))]
        # print('2',term2)
        term3 = [(sigmaSFR_not_corrected[i]*np.cos(newi)*np.radians(err_old_i[u])*np.sin(oldi[u]))/(np.cos(oldi[u]))**2 for i in range(len(kpc_r))]
        # print('3',term3)
        err_sigmaSFR= [np.sqrt(term1[i]**2 + term2[i]**2 + term3[i]**2) for i in range(len(kpc_r))]
        # print('err',err_sigmaSFR)
        rel_err_sigmaSFR = [err_sigmaSFR[i]/sigmaSFR_corrected[i] for i in range(len(kpc_r))]
        print('relerrSFR',rel_err_sigmaSFR)
        return rel_err_sigmaSFR
        # err_sigmaSFR= np.sqrt(((sigmaSFR_err*np.cos(np.radians(i))/np.cos(np.radians(i_old)))**2 + (sigmaSFR*np.sin(np.radians(i))*err_new_i/np.cos(np.radians(i_old)))**2))
        # return err_sigmaSFR

# print(rel_err_sigmaSFR)

# def err_sigmagas(sigmaHI_corr,sigmaH2_corr,sigmaHI_no_corr,mu_no_corr,err_mu=0):
#         # print('#################################################')
#         # print(sigmaHI_corr)
#         # print(sigmaH2_corr)
#         # print(sigmaHI_no_corr)
#         # sigmaHI_no_corr=(g_Msun/(cm_pc**2) )*sigmaHI_no_corr

#         #change this. its not mu but molfrac thats needed in the equation
#         sigmaH2_no_corr=sigmaHI_no_corr*(mu_no_corr/(1-mu_no_corr))
#         #to get corrected sigmagas
#         if switch['incl_moldat']=='Yes':
#                 sigmagas_corr=sigmaHI_corr+ sigmaH2_corr #columns_as_arrays[2]+columns_as_arrays[3]
#                 sigmagas_no_corr=sigmaHI_no_corr+sigmaH2_no_corr
#         else:
#                 sigmagas_corr=sigmaHI_corr
#                 sigmagas_no_corr=sigmaHI_no_corr
#         mu_err=[err_mu*mu_no_corr[i] for i in range(len(kpc_r))]
#         sigmaHI_err=err_sigmaHI(sigmaHI_corr,sigmaHI_no_corr,0.06,new_i,old_i)[0]
#         #might have to changed in case of moldat==yes
#         sigmaH2_err=[0 for i in range(len(kpc_r))]
#         sigmagas_err=[np.sqrt(sigmaHI_err[i]**2 + sigmaH2_err[i]**2) for i in range(len(kpc_r))]
#         # print('########################################################################')
#         print('step1',sigmagas_err)
#         term1= [(sigmagas_err[i]*(3*mu_no_corr[i]/(4-mu_no_corr[i]))) for i in range(len(kpc_r))]
#         print('term1',term1)
#         term2= [(sigmagas_no_corr[i]*(12*mu_err[i]/(4-mu_no_corr[i])**2)) for i in range(len(kpc_r))]
#         # print('term2',term2)
#         err_sigmagas= [np.sqrt((term1[i]**2 + term2[i]**2)) for i in range(len(kpc_r))]
#         print('err',err_sigmagas)
#         rel_err_sigmagas= [err_sigmagas[i]/sigmagas_corr[i] for i in range(len(kpc_r))]
#         print('rel_err_sigmagas',rel_err_sigmagas)
#         return rel_err_sigmagas

def err_sigmagas(sigmaHI_corr,sigmaH2_corr,sigmaHI_no_corr,sigmaH2_no_corr,err_mu=0):
        # print('#################################################')
        # print(sigmaHI_corr)
        # print(sigmaH2_corr)
        # print(sigmaHI_no_corr)
        # sigmaHI_no_corr=(g_Msun/(cm_pc**2) )*sigmaHI_no_corr
        mu= params['mu']
        #change this. its not mu but molfrac thats needed in the equation
        # sigmaH2_no_corr=sigmaHI_no_corr*(mu_no_corr/(1-mu_no_corr))
        #to get corrected sigmagas
        if switch['incl_moldat']=='Yes':
                sigmagas_corr=sigmaHI_corr+ sigmaH2_corr #columns_as_arrays[2]+columns_as_arrays[3]
                sigmagas_no_corr=sigmaHI_no_corr+sigmaH2_no_corr
                sigmaH2_err=err_sigmaH2(sigmaH2_no_corr,sigmaH2_corr,0.06,)
        else:
                sigmagas_corr=sigmaHI_corr
                sigmagas_no_corr=sigmaHI_no_corr
                sigmaH2_err=[0 for i in range(len(kpc_r))]

        mu_err=[err_mu*mu for i in range(len(kpc_r))]
        sigmaHI_err=err_sigmaHI(sigmaHI_corr,sigmaHI_no_corr,0.06,new_i,old_i)[0]
        #might have to changed in case of moldat==yes
        sigmagas_err=[np.sqrt(sigmaHI_err[i]**2 + sigmaH2_err[i]**2) for i in range(len(kpc_r))]
        # print('########################################################################')
        # print('step1',sigmagas_err)
        term1= [(sigmagas_err[i]*(3*mu/(4-mu))) for i in range(len(kpc_r))]
        # print('term1',term1)
        term2= [(sigmagas_no_corr[i]*(12*mu_err[i]/(4-mu)**2)) for i in range(len(kpc_r))]
        # print('term2',term2)
        err_sigmagas= [np.sqrt((term1[i]**2 + term2[i]**2)) for i in range(len(kpc_r))]
        # print('err',err_sigmagas)
        rel_err_sigmagas= [err_sigmagas[i]/sigmagas_corr[i] for i in range(len(kpc_r))]
        print('rel_err_sigmagas',rel_err_sigmagas)
        return rel_err_sigmagas

# print(rel_err_sigmagas)

def err_sigmatot(galaxy_name,sigmatot_corrected,sigmatot_not_corrected,percent_sigmatot_err,newi,oldi,new_dist=0,old_dist=0):
        newi=np.radians(newi)
        oldi=np.radians(oldi)
        # sigmatot_not_corrected=(g_Msun/(cm_pc**2))*sigmatot_not_corrected
        if galaxy_name=='m33':
                gamma=0.52*g_Msun/((cm_pc)**2)
                err_gamma=0.1*g_Msun/((cm_pc)**2)
                mu0=18.01
                kpc_Rd=1.82*cm_kpc 
                c36=24.8
                new_dist=new_dist*cm_kpc*(10**3)
                old_dist=old_dist*cm_kpc*(10**3)
                cm_r=list(cm_kpc*(raw_data.iloc[:,0]))
                
                mu=[(mu0+(1.10857*(cm_r[i]/kpc_Rd))) for i in range(len(kpc_r))]
                # print('mu',mu)
                sigma_mu=[((1.10857)/kpc_Rd)*cm_r[i]*np.sqrt(((err_new_dist*cm_kpc*(10**3))/(new_dist))**2 + ((new_dist*err_old_dist[0]*cm_kpc*(10**3))/(old_dist[0])**2)**2) for i in range(len(kpc_r))]
                # print('sigmamu',sigma_mu)
                a=[10**(-0.4*(mu[i]-c36)) for i in range(len(kpc_r))]
                # print('a',a)
                da_dt=[(np.log(10))*a[i] for i in range(len(kpc_r))]
                # print('dadt',da_dt)
                sigma_a= [-0.4*da_dt[i]*sigma_mu[i] for i in range(len(kpc_r))]
                # print('sigmaa',sigma_a)
                term1= [err_gamma*a[i] for i in range(len(kpc_r))]
                # print('term1',term1)
                term2= [gamma*sigma_a[i] for i in range(len(kpc_r))]
                # print('term2',term2)
                sigmatotal=[gamma*a[i] for i in range(len(kpc_r))]
                # print('totcalc', sigmatotal)
                # print('actual',sigmatot_corrected)
                err_sigmatot=[np.sqrt(term1[i]**2 + term2[i]**2) for i in range(len(kpc_r))]
                # print('err',err_sigmatot)
                # rel_err_sigmatot = [err_sigmatot[i]/sigmatot_not_corrected[i] for in range(len(kpc_r))]
                
        else:
                sigmatot_err=percent_sigmatot_err*sigmatot_not_corrected
                term1=[(sigmatot_err[i]*np.cos(newi)/np.cos(oldi[0])) for i in range(len(kpc_r))]
                # print('1',term1)
                term2=[(sigmatot_not_corrected[i]*np.sin(newi)*np.radians(err_new_i)/np.cos(oldi[0])) for i in range(len(kpc_r))]
                # print('2',term2)
                term3=[(sigmatot_not_corrected[i]*np.cos(newi)*np.radians(err_old_i[0])*np.sin(oldi[0]))/((np.cos(oldi[0]))**2) for i in range(len(kpc_r))]
                # print('3',term3)
                err_sigmatot= [np.sqrt(term1[i]**2 + term2[i]**2 + term3[i]**2) for i in range(len(kpc_r))]
        
        rel_err_sigmatot = [err_sigmatot[i]/sigmatot_corrected[i] for i in range(len(kpc_r))]
        # print('rel err tot',rel_err_sigmatot)
        # rel_err_sigmatot=[rel_err_sigmatot[i]*100 for i in range(len(kpc_r))]
        print('rel_err_sigmatot',rel_err_sigmatot)
        return rel_err_sigmatot
        #do the other derivation with upsilon for m33

# print(rel_err_sigmatot)

def err_T_q(T_corrected,q_corrected):
        # print(T_corrected,q_corrected)
        err_T=np.std(T_corrected)
        err_q=np.std(q_corrected)
        # print('*******************',err_T,err_q)
        rel_err_T=[err_T/T_corrected[i] for i in range(len(kpc_r))]
        rel_err_q=[err_q/q_corrected[i] for i in range(len(kpc_r))]
        print('err_T',rel_err_T)
        print('err_q',rel_err_q)
        return rel_err_T,rel_err_q
#################################################################################################################

#################################################################################################################
#calculating relative errors 
hregs = ['subsonic', 'supersonic']
for hreg in hregs:
        os.chdir(os.path.join(base_path,'inputs'))
        exps = np.load(f'scal_exponents_{hreg}.npy')
        r = kpc_r.size
        if hreg=='supersonic':
                rel_err_omega =err_omega(columns_as_arrays[5],vcirc_interpolated,kpc_r,err_interpolated_df.iloc[:,1],old_dist/(cm_kpc*(10**3)),old_i,new_dist,new_i)
        else:
                rel_err_omega =err_omega(columns_as_arrays[5],vcirc_interpolated,kpc_r,err_interpolated_df.iloc[:,1],old_dist,old_i,new_dist,new_i)
        
        if galaxy_name=='m33':
                rel_err_sigmatot=err_sigmatot(galaxy_name,columns_as_arrays[1],columns_as_arrays_reverted[0],0.1,new_i,old_i,new_dist,old_dist)
        else:
                rel_err_sigmatot=err_sigmatot(galaxy_name,columns_as_arrays[1],columns_as_arrays_reverted[0],0.1,new_i,old_i,new_dist,old_dist)

        if galaxy_name=='m31':
                rel_err_sigmagas=err_sigmagas(columns_as_arrays[2],columns_as_arrays[3],columns_as_arrays_reverted[1],columns_as_arrays_reverted[2],0)
        else:
                rel_err_sigmagas=err_sigmagas(columns_as_arrays[2],columns_as_arrays[3],columns_as_arrays_reverted[1],columns_as_arrays_reverted[2],0)

        rel_err_sigmasfr=err_sigmaSFR(columns_as_arrays[-2],columns_as_arrays_reverted[-1],0.1,new_i,old_i)
        rel_err_T,rel_err_q=err_T_q(columns_as_arrays[-1],columns_as_arrays[4])
        rel_err = np.array([rel_err_q, rel_err_omega, rel_err_sigmagas, rel_err_sigmatot,rel_err_sigmasfr, rel_err_T])

        relerr_quan = np.sqrt(np.matmul(exps**2,rel_err**2))
        err_quantities = model_f[1:]*relerr_quan
        # print('P',rel_err_omega)
#################################################################################################################
#################################################################################################################
#inputting errors into a pickle file
        os.chdir(os.path.join(base_path,'outputs'))
        with open(f'errors_{hreg}.out', 'wb') as f:
                pickle.dump(err_quantities, f)
print('Found the errors from the scaling relations')
 
