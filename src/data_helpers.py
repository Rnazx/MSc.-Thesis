import numpy as np
from scipy.interpolate import griddata
from sympy import nan

deg_rad = 180e0/np.pi
arcmin_deg = 60e0
arcsec_deg = 3600e0


#REQUIRED FUNCTIONS
###########################################################################################################################################
#extrapolation
def interpolation(x,y,standard):
    interpolated_data = griddata(x, y, standard, method='linear', fill_value=nan, rescale=False)
    return interpolated_data

def replace_conversion(df, substring_to_replace, replacement_string):
    # Create a dictionary for column renaming
    rename_dict = {col: col.replace(substring_to_replace, replacement_string) for col in df.columns}
    # Rename the columns using the dictionary
    updated_df = df.rename(columns=rename_dict)

    return updated_df

def keep_substring_columns(dataframe, substring):
    # Get the columns that contain the specified substring
    filtered_columns = [col for col in dataframe.columns if substring in col]
    
    # Create a new DataFrame with only the filtered columns
    result_df = dataframe[filtered_columns].copy()
    
    return result_df, filtered_columns

# get the next or previous column of the dataframe
def get_adjacent_column(df, col, next = True): 
    index_of_target = df.columns.get_loc(col)
    if index_of_target < len(df.columns) - 1:
        if next: #get the column right side of col
            return df.columns[index_of_target+1]
        else: #get the column left side of col
            return df.columns[index_of_target-1]
    else:
        print("The target column is the last column.")

def remove_data(raw_data, data_rem):
    data_removed = keep_substring_columns(raw_data, data_rem)[1]
    radii_removed = [get_adjacent_column(raw_data, dr, False) for dr in data_removed]
    raw_data.drop(columns=radii_removed+data_removed,inplace=True)
    return raw_data

def df_interpolation(df, radii_df, standard):
    result_df = df.copy()
    for cols in radii_df.columns:
        result_df[get_adjacent_column(df, cols)] = interpolation(df[cols],df[get_adjacent_column(df, cols)],standard)
        result_df.drop(columns=[cols], inplace=True)
    result_df.insert(0, 'kpc_r', standard)
    return result_df

# def inclination_correction(df, i_new, i_old):
#     radii_df = keep_substring_columns(df, 'r ')[0]
#     #df.iloc[:,1:].multiply(np.cos(i_new)/np.cos(i_old), axis = 1)
#     return df

def incl_distance_correction(df, distance_new, distance_old, i_new, i_old):
    def find_and_multiply_column(dataframe, substring, dist_multiplier = 1, incl_multiplier = 1):
    # Create a copy of the DataFrame to avoid modifying the original
        result_df = dataframe.copy()
        iter_dist = 0
        iter_incl = 0
        for col in result_df.columns:
            if substring in col:
                try:
                    result_df[col] = result_df[col] * dist_multiplier[iter_dist]
                except:
                    result_df[col] = result_df[col] * dist_multiplier
                iter_dist+=1
            else:
                try:
                    result_df[col] = result_df[col] * incl_multiplier[iter_incl]
                except:
                    result_df[col] = result_df[col] * incl_multiplier
                iter_incl+=1
                
        return result_df
    #radii_df = keep_substring_columns(df, 'r ')[0]
    #radii_df = radii_df.drop(columns='error kms')
    #convert arcmin to kpc
    df = find_and_multiply_column(df, 'arcmin', ((distance_new*1000)/(arcmin_deg*deg_rad)))
    #convert arcsec to kpc
    df = find_and_multiply_column(df, 'arcsec', ((distance_new*1000)/(arcsec_deg*deg_rad)))
    #change the names
    df = replace_conversion(df, 'arcmin', 'kpc;')
    df = replace_conversion(df, 'arcsec', 'kpc;')
    #x = df.copy()
    #distance and inclination correction
    df = find_and_multiply_column(df, 'kpc', distance_new/distance_old, np.cos(i_new)/np.cos(i_old))
    
    #print(x.compare(df))
    return df

def molfrac_to_H2(df):
    df = df.copy()
    molfrac_data = keep_substring_columns(df, 'molfrac')
    if molfrac_data[0].empty:
        return df
    else:
        HI_data = keep_substring_columns(df, 'HI')
        sigma_H2 = HI_data[0].multiply((molfrac_data[0]/(1-molfrac_data[0])).values, axis = 0)
        index_of_HI = df.columns.get_loc(HI_data[1][0])
        df.insert(index_of_HI+1, 'sigma_H2', sigma_H2)
        df.drop(columns=molfrac_data[1], inplace=True)
        return df
    
def add_temp(temp_fit, df):
    m = temp_fit[0]
    c = temp_fit[1]
    r = df.iloc[:,0].to_numpy().flatten()
    try:
        T = m*r +c
    except:
        T = interpolation(m, c, r)
    df.insert(len(df.columns), 'T', T)
    return

def vcirc_to_qomega(df):
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
        df.drop(columns=vcirc_data[1], inplace=True)
        return df

