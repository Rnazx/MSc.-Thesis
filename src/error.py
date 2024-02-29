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

# oDIST,oINC,nDIST,nINC
#distance correction
nDIST= corrections.iloc[-1,0]
err_nDIST=corrections.iloc[-1,2]
oDIST = corrections.iloc[:-1,0].values
err_oDIST=corrections.iloc[:-1,2].values

#inclination correction
nINC= corrections.iloc[-1,1] #deg
err_nINC= corrections.iloc[-1,3] #deg
oINC= corrections.iloc[:-1,1].values #used new_i as no inclination correction is needed for Claude data
err_oINC=corrections.iloc[:-1,3].values

#################################################################################################################
#raw_data from combined_data.csv is converted to correct units, but no D/i correction done here
raw_data = incl_distance_correction(raw_data, distance_new=oDIST, distance_old=oDIST,\
                        i_new=np.radians(oINC), i_old=np.radians(oINC))
raw_data_corrected = incl_distance_correction(raw_data, distance_new=nDIST, distance_old=oDIST,\
                        i_new=np.radians(nINC), i_old=np.radians(oINC))
# os.chdir(os.path.join(base_path, 'data','model_data', f'{galaxy_name}_data'))

#unwanted cols are removed from raw_data
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

try:
    data_rem = pd.read_csv(f'removed_data_{galaxy_name}.csv', dtype=str)
except:
    data_rem = []
for d in data_rem:
    remove_data(raw_data, d)
    remove_data(raw_data_corrected, d)

temp_fit = np.genfromtxt(f'temp_{galaxy_name}.csv', skip_header = 1,delimiter=',')

#################################################################################################################

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
    
def molfrac_to_H2_new(df,remove_molfrac=False):
    df = df.copy()
    molfrac_data = keep_substring_columns(df, 'molfrac')
    if molfrac_data[0].empty:
        return df
    else:
        HI_data = keep_substring_columns(df, 'HI')
        sigma_H2 = HI_data[0].multiply((molfrac_data[0]/(1-molfrac_data[0])).values, axis = 0)
        index_of_HI = df.columns.get_loc(HI_data[1][0])
        df.insert(index_of_HI+1, 'sigma_H2', sigma_H2)
        if remove_molfrac:
            df.drop(columns=molfrac_data[1], inplace=True)
        return df
#################################################################################################################

#################################################################################################################

raw_data= vcirc_to_qomega_new(raw_data) #raw_data contain v_circ also
raw_data_radii_df = keep_substring_columns(raw_data, 'r ')[0]

raw_data_corrected= vcirc_to_qomega_new(raw_data_corrected) #raw_data contain v_circ also
raw_data_radii_df_corrected = keep_substring_columns(raw_data_corrected, 'r ')[0]
# print(raw_data.iloc[:,-1][3]*g_Msun/((s_Myr*1e3)*(cm_pc**2)))


err_data = incl_distance_correction(err_data, distance_new=nDIST, distance_old=oDIST,\
            i_new=np.radians(nINC), i_old=np.radians(oINC)) #correcting errors for D and i
err_radii_df = keep_substring_columns(err_data, 'R')[0]

coarsest_radii_mask = raw_data_radii_df.isnull().sum().idxmax()

print("Coarsest radii is {} and the data it corresponds to is {}:".format(coarsest_radii_mask,get_adjacent_column(raw_data,coarsest_radii_mask)))
kpc_r = raw_data_radii_df[coarsest_radii_mask].to_numpy()
# kpc_r_not_interpolated=kpc_r

#print(raw_data)
interpolated_df = df_interpolation(raw_data,raw_data_radii_df, kpc_r)
interpolated_df_corrected = df_interpolation(raw_data_corrected,raw_data_radii_df_corrected, kpc_r)

print(interpolated_df)
print(interpolated_df_corrected)
print(err_data)

##################################################################
# max_rows_A = max(len(err_data[col]) for col in err_data.columns)
# max_rows_B = max(len(interpolated_df[col]) for col in interpolated_df.columns)
# # print(max_rows_B)
# # print(len(interpolated_df.columns))
# # print(len(err_data.columns))

# err_data_long = pd.DataFrame(0, index=range(max(max_rows_A,max_rows_B)), columns=err_data.columns)
# if max_rows_A>max_rows_B: #err_data has greater length. nan added to interpolated_df
#        for i in range(len(interpolated_df.columns)):
#                 interpolated_df.iloc[abs(max_rows_A-max_rows_B):, i] = np.nan
# else:
#        print('this part is working now')
#        for i in range(len(err_data.columns)):
#                 err_data_col_list=np.array(err_data.iloc[:,i])
#                 # print('########################################3')
#                 # print(i,err_data_col_list)
#                 # print('########################################3')
#                 # for j in range(max_rows_A,max_rows_B):
#                 err_add_nan = np.concatenate((err_data_col_list, np.full(abs(max_rows_A-max_rows_B), np.nan)))
#                 #        err_data_col_list[j]=np.nan
#                 err_data_long.iloc[:,i]=err_add_nan
#                 # err_data.iloc[abs(max_rows_A-max_rows_B):, i] = np.nan

# # print(err_data_long)
# err_radii_df = keep_substring_columns(err_data_long, 'R')[0]
# err_interpolated_df = df_interpolation(err_data_long,err_radii_df, kpc_r)

##################################################################
err_interpolated_df = df_interpolation(err_data,err_radii_df, kpc_r)

interpolated_df = molfrac_to_H2_new(interpolated_df)
interpolated_df_corrected = molfrac_to_H2_new(interpolated_df_corrected)

add_temp(temp_fit,interpolated_df)
add_temp(temp_fit,interpolated_df_corrected)

nan_mask                  = np.isnan(interpolated_df)
nan_mask_corrected        = np.isnan(interpolated_df_corrected)
interpolated_df           = interpolated_df[~(nan_mask.sum(axis=1)>0)]
interpolated_df_corrected = interpolated_df_corrected[~(nan_mask_corrected.sum(axis=1)>0)]

#remove nans
nan_mask = np.isnan(raw_data)
raw_data = raw_data[~(nan_mask.sum(axis=1)>0)]

nan_mask_corrected = np.isnan(raw_data_corrected)
raw_data_corrected = raw_data_corrected[~(nan_mask_corrected.sum(axis=1)>0)]

nan_mask = np.isnan(err_interpolated_df)
err_interpolated_df = err_interpolated_df[~(nan_mask.sum(axis=1)>0)]

#interpolated_df.dropna()
# Changed for m51 and ngc6946
if galaxy_name == 'm31' or galaxy_name == 'm33':
    m_to_gconv = 1e3
else:
    m_to_gconv = 1

if galaxy_name=='m31': #since there is extra molfrac list here at 2nd last element
    conv_factors=np.array([1, (g_Msun/(cm_pc**2) ), g_Msun/(cm_pc**2), g_Msun/(cm_pc**2), 1,cm_km/cm_kpc,cm_km,
            g_Msun/((s_Myr*m_to_gconv)*(cm_pc**2)),1,1])
else:
    conv_factors=np.array([1, (g_Msun/(cm_pc**2) ), g_Msun/(cm_pc**2), g_Msun/(cm_pc**2), cm_km,1,cm_km/cm_kpc,
            g_Msun/((s_Myr*m_to_gconv)*(cm_pc**2)),1])
# print('from data common',interpolated_df)

interpolated_df           = interpolated_df*conv_factors
interpolated_df_corrected = interpolated_df_corrected*conv_factors

if galaxy_name=='m31' or galaxy_name=='m33':
    conv_factors=np.array([1,cm_km,1])
else:
    conv_factors=np.array([1,cm_km])

err_interpolated_df = err_interpolated_df*conv_factors

print(interpolated_df)
print(interpolated_df_corrected)
print(err_interpolated_df)

# print(old_dist,old_i)
#nDIST, oDIST, nINC, oINC
#functions for error calculations start here

#unit conversions for all functions
oDIST=oDIST*cm_kpc*(10**3)
nDIST=nDIST*cm_kpc*(10**3)
        
oINC=np.radians(oINC)
nINC=np.radians(nINC)

err_nDIST=err_nDIST*cm_kpc*(10**3)
err_oDIST=err_oDIST*cm_kpc*(10**3)

err_nINC=np.radians(err_nINC)
err_oINC=np.radians(err_oINC)

kpc_r=interpolated_df['kpc_r']
cm_r=kpc_r*cm_kpc

def err_omega(omega,vcirc,sigma_v): #omega=columns_as_arrays[4]
        if galaxy_name=='m31' or galaxy_name=='m33':
                u=4
        else:
                u=3
        # sigma_v=[sigma_v[i]/np.sin(nINC) for i in range(len(kpc_r))]
        # print('inc and d',oDIST[u],nDIST,oINC[u],nINC)
        term1 = [(oDIST[u]*np.sin(oINC[u])*sigma_v[i])/(nDIST*cm_r[i]*np.sin(nINC)) for i in range(len(kpc_r))]
        # print('1',term1)

        term2 = [(vcirc[i]*oDIST[u]*np.cos(oINC[u])*err_oINC[u])/((np.sin(nINC))*cm_r[i]*nDIST) for i in range(len(kpc_r))]       
        # print('2',term2)
        term3 = [(vcirc[i]*oDIST[u]*np.sin(oINC[u])*np.cos(nINC)*err_nINC)/(((np.sin(nINC))**2)*cm_r[i]*nDIST) for i in range(len(kpc_r))]
        # print('3',term3)
        term4 = [(vcirc[i]*np.sin(oINC[u])*err_oDIST[u])/((np.sin(nINC))*cm_r[i]*nDIST) for i in range(len(kpc_r))]
        # print('4',term4)
        term5 = [(vcirc[i]*np.sin(oINC[u])*oDIST[u]*err_nDIST)/((np.sin(nINC))*cm_r[i]*(nDIST**2)) for i in range(len(kpc_r))]
        # print('5',term5)
        err_omega     = [np.sqrt(term1[i]**2 +term2[i]**2 + term3[i]**2 + term4[i]**2 + term5[i]**2) for i in range(len(kpc_r))]
        rel_err_omega = err_omega/omega
        # print('omega err',rel_err_omega)
        return rel_err_omega

# err_omega(interpolated_df_corrected['\Omega'],interpolated_df['vcirc_claude kms'],err_interpolated_df['error vcirc kms'])

def err_sigmaHI(sigmaHI_corrected,sigmaHI_not_corrected,percent_sigmaHI_err): #sigmaHI_err is set to 6 percent, sigmaHI=columns_as_arrays[2]
        if galaxy_name=='m31' or galaxy_name=='m33':
                u=2
        else:
                u=1
        sigmaHI_err=sigmaHI_not_corrected*percent_sigmaHI_err
        term1 = [(sigmaHI_err[i]*np.cos(nINC))/np.cos(oINC[u]) for i in range(len(kpc_r))]
        # print('1',term1)
        term2 = [(sigmaHI_not_corrected[i]*np.sin(nINC)*err_nINC)/np.cos(oINC[u]) for i in range(len(kpc_r))]
        # print('2',term2)
        term3 = [(sigmaHI_not_corrected[i]*np.cos(nINC)*err_oINC[u]*np.sin(oINC[u]))/((np.cos(oINC[2]))**2) for i in range(len(kpc_r))]
        # print('3',term3)
        err_sigmaHI= [np.sqrt(term1[i]**2 + term2[i]**2 + term3[i]**2) for i in range(len(kpc_r))]
        rel_err_sigmaHI = [err_sigmaHI[i]/sigmaHI_corrected[i] for i in range(len(kpc_r))]
        # print('rel_errHI',rel_err_sigmaHI)
        # print('errHI',err_sigmaHI)
        
        return err_sigmaHI,rel_err_sigmaHI

# err_sigmaHI(interpolated_df_corrected['sigma_HI_claude'],interpolated_df['sigma_HI_claude'],0.06)

#no need of corrected sigmas here as relative error isnt required
def err_sigmaH2(sigmaHI_not_corrected,molfrac_or_sigmaH2,percent_sigmaHI_err,molfrac_err,percent_sigmaH2_err):
        # sigmaHI_not_corrected=(g_Msun/(cm_pc**2) )*sigmaHI_not_corrected
        
        if switch['incl_moldat']=='No':
                print('considering ONLY HI surface density')
                err_sigmaH2=[0 for i in range(len(kpc_r))]
                
        else:
                if galaxy_name=='m31':
                        molfrac=molfrac_or_sigmaH2
                        # sigmaH2 = [(molfrac[i]/(1-molfrac[i]))*sigmaHI_not_corrected[i] for i in range(len(kpc_r))]
                        molfrac_err_tot = [np.sqrt(((raw_data['molfrac'][i]*np.sin(nINC)*err_nINC)/(np.cos(oINC[-1])))**2 + ((raw_data['molfrac'][i]*np.cos(nINC)*err_oINC[-1]*np.sin(oINC[-1]))/(np.cos(oINC[-1]))**2)**2 + (molfrac_err[i]*np.cos(nINC))/np.cos(oINC[-1])**2) for i in range(len(kpc_r))]
                        sigmaHI_err=sigmaHI_not_corrected*percent_sigmaHI_err
                        
                        term1= [(sigmaHI_err[i]*molfrac[i])/(1-molfrac[i]) for i in range(len(kpc_r))]
                        term2= [(sigmaHI_not_corrected[i]*molfrac_err_tot[i])/((1-molfrac[i])**2) for i in range(len(kpc_r))]
                        err_sigmaH2 = np.sqrt(term1**2+ term2**2)
                        # rel_err_sigmaH2 = [err_sigmaH2[i]/sigmaH2 for i in range(len(kpc_r))]
                        # print('sigmah2err',rel_err_sigmaH2)

                else:
                        if galaxy_name=='m33':
                            u=3
                        else:
                            u=2
                        sigmaH2_not_corrected=molfrac_or_sigmaH2
                        percent_sigmaH2_err=0.06
                        sigmaH2_err=sigmaH2_not_corrected*percent_sigmaH2_err

                        # index_for_sigmaH2=((list(columns_as_arrays)).index('sigma_H2'))-1
                        term1 = [(sigmaH2_err[i]*np.cos(nINC))/np.cos(oINC[u]) for i in range(len(kpc_r))]
                        # print('1',term1)
                        term2 = [(sigmaH2_not_corrected[i]*np.sin(nINC)*err_nINC)/np.cos(oINC[u]) for i in range(len(kpc_r))]
                        # print('2',term2)
                        term3 = [(sigmaH2_not_corrected[i]*np.cos(nINC)*err_oINC[u]*np.sin(oINC[u]))/(np.cos(oINC[u]))**2 for i in range(len(kpc_r))]
                        # print('3',term3)
                        err_sigmaH2= [np.sqrt(term1[i]**2 + term2[i]**2 + term3[i]**2) for i in range(len(kpc_r))]
                        # rel_err_sigmaH2 = [err_sigmaH2[i]/sigmaH2_not_corrected for i in range(len(kpc_r))]
                        # print('sigmah2err',rel_err_sigmaH2)
        return err_sigmaH2

# errsigmah2=err_sigmaH2(interpolated_df['sigma_HI_claude'],interpolated_df['sigma_H2'],0.06,err_interpolated_df['error molfrac'],0.06)
# print(errsigmah2)

def err_sigmagas(sigmaHI_corr,sigmaH2_corr,sigmaHI_no_corr,sigmaH2_no_corr,err_mu):
        mu= params['mu']
        #change this. its not mu but molfrac thats needed in the equation
        # sigmaH2_no_corr=sigmaHI_no_corr*(mu_no_corr/(1-mu_no_corr))
        #to get corrected sigmagas
        sigmaH2_err=err_sigmaH2(sigmaHI_no_corr,sigmaH2_no_corr,0.06,err_interpolated_df['error molfrac'],0.06)
        sigmaHI_err=err_sigmaHI(sigmaHI_corr,sigmaHI_no_corr,0.06)[0]

        if switch['incl_moldat']=='Yes':
                sigmagas_corr=sigmaHI_corr+ sigmaH2_corr #columns_as_arrays[2]+columns_as_arrays[3]
                sigmagas_no_corr=sigmaHI_no_corr+sigmaH2_no_corr
                # sigmaH2_err=err_sigmaH2(sigmaH2_no_corr,sigmaH2_corr,0.06)
        else:
                sigmagas_corr=sigmaHI_corr
                sigmagas_no_corr=sigmaHI_no_corr
                # sigmaH2_err=[0 for i in range(len(kpc_r))]

        mu_err=[err_mu*mu for i in range(len(kpc_r))]
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
        # print('rel_err_sigmagas',rel_err_sigmagas)
        return rel_err_sigmagas

# err_sigmagas(interpolated_df_corrected['sigma_HI_claude'],interpolated_df_corrected['sigma_H2'],interpolated_df['sigma_HI_claude'],interpolated_df['sigma_H2'],0.1)

def err_sigmaSFR(sigmaSFR_corrected,sigmaSFR_not_corrected,percent_sigmaSFR_err):
        # sigmaSFR_not_corrected=sigmaSFR_not_corrected
        sigmaSFR_err=sigmaSFR_not_corrected*percent_sigmaSFR_err
        if galaxy_name=='m31':
                u=-2
        else:
                u=-1
        term1 = [(sigmaSFR_err[i]*np.cos(nINC))/np.cos(oINC[u]) for i in range(len(kpc_r))]
        # print('1',term1)
        term2 = [(sigmaSFR_not_corrected[i]*np.sin(nINC)*err_nINC)/np.cos(oINC[u]) for i in range(len(kpc_r))]
        # print('2',term2)
        term3 = [(sigmaSFR_not_corrected[i]*np.cos(nINC)*err_oINC[u]*np.sin(oINC[u]))/(np.cos(oINC[u]))**2 for i in range(len(kpc_r))]
        # print('3',term3)
        err_sigmaSFR= [np.sqrt(term1[i]**2 + term2[i]**2 + term3[i]**2) for i in range(len(kpc_r))]
        # print('err',err_sigmaSFR)
        rel_err_sigmaSFR = [err_sigmaSFR[i]/sigmaSFR_corrected[i] for i in range(len(kpc_r))]
        # print('relerrSFR',rel_err_sigmaSFR)
        return rel_err_sigmaSFR

# err_sigmaSFR(interpolated_df_corrected['sigma_sfr'],interpolated_df['sigma_sfr'],0.1)

def err_sigmatot(sigmatot_corrected,sigmatot_not_corrected,percent_sigmatot_err):
        # sigmatot_not_corrected=(g_Msun/(cm_pc**2))*sigmatot_not_corrected
        if galaxy_name=='m33':
                gamma=0.52*g_Msun/((cm_pc)**2)
                err_gamma=0.1*g_Msun/((cm_pc)**2)
                mu0=18.01
                kpc_Rd=1.82*cm_kpc 
                c36=24.8
                
                mu=[(mu0+(1.10857*(cm_r[i]/kpc_Rd))) for i in range(len(kpc_r))]
                # print('mu',mu)
                sigma_mu=[((1.10857)/kpc_Rd)*cm_r[i]*np.sqrt(((err_nDIST*cm_kpc*(10**3))/(nDIST))**2 + ((nDIST*err_oDIST[0]*cm_kpc*(10**3))/(oDIST[0])**2)**2) for i in range(len(kpc_r))]
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
                # sigmatotal=[gamma*a[i] for i in range(len(kpc_r))]
                # print('totcalc', sigmatotal)
                # print('actual',sigmatot_corrected)
                err_sigmatot=[np.sqrt(term1[i]**2 + term2[i]**2) for i in range(len(kpc_r))]
                # print('err',err_sigmatot)
                # rel_err_sigmatot = [err_sigmatot[i]/sigmatot_not_corrected[i] for in range(len(kpc_r))]
                
        else:
                sigmatot_err=percent_sigmatot_err*sigmatot_not_corrected
                term1=[(sigmatot_err[i]*np.cos(nINC)/np.cos(oINC[0])) for i in range(len(kpc_r))]
                # print('1',term1)
                term2=[(sigmatot_not_corrected[i]*np.sin(nINC)*err_nINC/np.cos(oINC[0])) for i in range(len(kpc_r))]
                # print('2',term2)
                term3=[(sigmatot_not_corrected[i]*np.cos(nINC)*np.radians(err_oINC[0])*np.sin(oINC[0]))/((np.cos(oINC[0]))**2) for i in range(len(kpc_r))]
                # print('3',term3)
                err_sigmatot= [np.sqrt(term1[i]**2 + term2[i]**2 + term3[i]**2) for i in range(len(kpc_r))]
        
        rel_err_sigmatot = [err_sigmatot[i]/sigmatot_corrected[i] for i in range(len(kpc_r))]
        # print('rel_err_sigmatot',rel_err_sigmatot)
        return rel_err_sigmatot

# err_sigmatot(interpolated_df_corrected['sigma_tot'],interpolated_df['sigma_tot'],0.1)

def err_T_q(T_corrected,q_corrected):
        # print(T_corrected,q_corrected)
        err_T=np.std(T_corrected)
        err_q=np.std(q_corrected)
        # print('*******************',err_T,err_q)
        rel_err_T=[err_T/T_corrected[i] for i in range(len(kpc_r))]
        rel_err_q=[err_q/q_corrected[i] for i in range(len(kpc_r))]
        return rel_err_T,rel_err_q
#################################################################################################################
#################################################################################################################
# interpolated_df.rename(columns={'kms_vcirc_Kam': 'vcirc'}, inplace=True)
# interpolated_df_corrected.rename(columns={'kms_vcirc_Kam': 'vcirc'}, inplace=True)
# hregs = ['subsonic', 'supersonic']
# for hreg in hregs:
#         os.chdir(os.path.join(base_path,'inputs'))
#         exps = np.load(f'scal_exponents_{hreg}.npy')
#         r = kpc_r.size
#         if hreg=='supersonic':
#                 rel_err_omega = err_omega(interpolated_df_corrected['\Omega'],interpolated_df['vcirc'],err_interpolated_df['error vcirc kms'])

#         else:
#                 rel_err_omega = err_omega(interpolated_df_corrected['\Omega'],interpolated_df['vcirc'],err_interpolated_df['error vcirc kms'])

#         rel_err_sigmatot = err_sigmatot(interpolated_df_corrected['sigma_tot'],interpolated_df['sigma_tot'],0.1)
#         rel_err_sigmagas = err_sigmagas(interpolated_df_corrected['sigma_HI'],interpolated_df_corrected['sigma_H2'],interpolated_df['sigma_HI'],interpolated_df['sigma_H2'],0.1)
#         rel_err_sigmasfr = err_sigmaSFR(interpolated_df_corrected['sigma_sfr'],interpolated_df['sigma_sfr'],0.1)

#         rel_err_T,rel_err_q=err_T_q(interpolated_df_corrected['T'],interpolated_df_corrected['q'])
#         rel_err = np.array([rel_err_q, rel_err_omega, rel_err_sigmagas, rel_err_sigmatot,rel_err_sigmasfr, rel_err_T])

#         relerr_quan = np.sqrt(np.matmul(exps**2,rel_err**2))
#         err_quantities = model_f[1:]*relerr_quan

hregs = ['subsonic', 'supersonic']
for hreg in hregs:
        os.chdir(os.path.join(base_path,'inputs'))
        exps = np.load(f'scal_exponents_{hreg}.npy')
        r = kpc_r.size
        if galaxy_name=='m31':
                rel_err_omega=err_omega(interpolated_df_corrected['\Omega'],interpolated_df['vcirc_claude kms'],err_interpolated_df['error vcirc kms'])
        elif galaxy_name=='m33':
                print('try',interpolated_df['kms_vcirc_Kam'])
                rel_err_omega=err_omega(interpolated_df_corrected['\Omega'],interpolated_df['kms_vcirc_Kam'],err_interpolated_df['error vcirc kms'])
               
        print('ran till here')
        rel_err = np.array([rel_err_q, rel_err_omega, rel_err_sigmagas, rel_err_sigmatot,rel_err_sigmasfr, rel_err_T])

        relerr_quan = np.sqrt(np.matmul(exps**2,rel_err**2))
        err_quantities = model_f[1:]*relerr_quan
#################################################################################################################
#################################################################################################################
#inputting errors into a pickle file
        os.chdir(os.path.join(base_path,'outputs'))
        with open(f'errors_{hreg}.out', 'wb') as f:
                pickle.dump(err_quantities, f)
print('Found the errors from the scaling relations')
 