import numpy as np
from sympy import *
from fractions import Fraction
import pickle
from scipy.interpolate import griddata
import pandas as pd

current_directory=r'D:\Documents\Gayathri_college\MSc project\data\m51'

#REQUIRED FUNCTIONS
###########################################################################################################################################
#importing data from csv file to np array
def file_reader(list_of_file_names):
    dataframe_list=[]
    for i in list_of_file_names:
        df = pd.read_csv(current_directory+r'\{}.csv'.format(i))
        dataframe_list.append(df)
    return dataframe_list

#extrapolation
def extrapolation(list1,list2,standard):
    extrapolated_data = griddata(list1, list2, standard, method='linear', fill_value=nan, rescale=False)
    return extrapolated_data
###########################################################################################################################################

file_names=['RC_sofue+18','omega_sofue+18','q_valuesM51sofue+18','HI+H2 m51.','H2 m51.',
            'HI m51..','SFR_Halpha24 m51.','SFR_FUV24 m51.','smdf','CO vel dispersion schuster.','temperature']
dataframe_list=file_reader(file_names)

#to obtain radius data from every df
radius_list=[np.array(dataframe_list[i]['r']) for i in range(len(dataframe_list))]

#obtain observables from every df and append to a list named quant_list
col_names=['v','omega','q','sigma_gas','sigma_H2','sigma_HI','sigma_sfr','sigma_sfr','smdf','vel disp','temp']
quant_list=[np.array(dataframe_list[i][col_names[i]]) for i in range(len(dataframe_list))]

#find array with least length
array_with_least_length = min(radius_list, key=len) #this shows that temp data has least number of points

#extrapolating the data and appending the common radius list 
quant_list_ep=[np.round(extrapolation(radius_list[i],quant_list[i],radius_list[-1]),3) for i in range(len(dataframe_list))]
quant_list_ep.append(np.round(radius_list[-1],3))

#removing nan values for points whr interpolation is impossible
nan_max = np.argmax([np.sum(np.isnan(d)) for d in quant_list_ep])
nan_max_data = quant_list_ep[nan_max]
nan_mask = ~np.isnan(nan_max_data)

nandeleted_data = []
for i,d in enumerate(quant_list_ep):
    nandeleted_data.append(d[nan_mask])
nandeleted_data = tuple(nandeleted_data)
print(nandeleted_data)

#storing that in a pickle file
code_directory=r'D:\Documents\Gayathri_college\MSc project\codes'
with open(code_directory+'\data_{}.pickle'.format('M51'), 'wb') as f:
    pickle.dump(nandeleted_data, f)





