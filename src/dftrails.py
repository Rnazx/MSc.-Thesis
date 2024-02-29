import numpy as np
import pickle
import os
import pandas as pd
from helper_functions import *
from data_helpers import *

# galaxy_name='m31'
base_path = os.environ.get('MY_PATH')
print('basepath',base_path)
galaxy_name = os.environ.get('galaxy_name')

#get parameter values
params = parameter_read(os.path.join(base_path,'inputs','parameter_file.in'))

#getting available data and error
os.chdir(os.path.join(base_path, 'data','model_data', f'{galaxy_name}_data'))
raw_data = pd.read_csv(f'combined_data_{galaxy_name}.csv', skiprows=1)
err_data= pd.read_csv(f'error_combined_{galaxy_name}.csv')
corrections = pd.read_csv(f'correction_data_{galaxy_name}.csv', skiprows=1, index_col=0)

#distance correction
new_dist= corrections.iloc[-2,0]
err_new_dist=corrections.iloc[-1,0]
old_dist = corrections.iloc[:-1,0].values

#inclination correction
new_i= corrections.iloc[-2,1] #deg
err_new_i= corrections.iloc[-1,1] #deg
old_i= corrections.iloc[:-1,1].values #used new_i as no inclination correction is needed for Claude data

raw_data = incl_distance_correction(raw_data, distance_new=new_dist, distance_old=old_dist,\
                        i_new=np.radians(new_i), i_old=np.radians(old_i))

#correcting error data based on availability
# for M31- vcirc and molfrac data
err_m31_old_dist=[0.785,0.78]
err_m31_old_i=[75,77.5]
err_data = incl_distance_correction(err_data, distance_new=new_dist, distance_old=err_m31_old_dist,\
            i_new=np.radians(new_i), i_old=np.radians(err_m31_old_i))

err_data = vcirc_to_qomega(err_data)
err_radii_df = keep_substring_columns(err_data, 'r ')[0]

print(err_data)
print(err_radii_df)
os.chdir(os.path.join(base_path,'inputs'))
with open('zip_data.in', 'rb') as f:
    kpc_r, data_pass = pickle.load(f)
