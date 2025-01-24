import os
from scipy.io import loadmat, savemat
import pandas as pd
import re
from scipy.stats import skew
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import gc
import seaborn as sns
import math
from math import asin, degrees, sqrt
import shutil
import zipfile
from datetime import datetime
import time
import csv



def match_freq_pol(file_name):
    pattern = r'(\d+(?:\.\d+)?)MHz_([XY])pol'
    match  = re.search(pattern, file_name)
    if match:
        frequency = match.group(1)
        pol = match.group(2)  
        # pol = 0 if pol == 'X' else 1
    return frequency, pol

def count_e_tot(dim1, dim2, dim3, input_mat_file):
    E_total= np.zeros((dim1, dim2, dim3), dtype=float)
    E_phi = input_mat_file['Ephi']
    E_theta = input_mat_file['Etheta']
    E_total= np.sqrt(np.abs(E_phi)**2 + np.abs(E_theta)**2)
    return E_total

def count_e_norm(dim1, dim2, dim3,e_tot):
    E_total_max_matrix = np.zeros((dim1,dim2, dim3), dtype = float)
    E_norm = np.zeros((dim1, dim2, dim3), dtype=float)
    for antenna in range(dim3):
        E_total_vector_of_antenna = e_tot[:,:,antenna]
        e_tot_m = np.max(np.max(E_total_vector_of_antenna))
        E_total_max_matrix[:,:,antenna] = e_tot_m
    E_norm = 20 * np.log10(e_tot/ E_total_max_matrix)
    return E_norm

# detect problamtatic regions


def process_xpol(antenna_data):
    """Processes data for Xpol files."""
    df_0_180 = pd.DataFrame(antenna_data[:360, :theta_range + 1])  # phi 0-179.5 degrees
    df_180_360 = pd.DataFrame(antenna_data[360:, :theta_range + 1])  # phi 180-360 degrees
    return pd.concat([df_180_360, df_0_180]).reset_index(drop=True)

def process_ypol(antenna_data):
    """Processes data for Ypol files."""
    df_0_180 = pd.DataFrame(antenna_data[180:540, :theta_range + 1])  #phi  0-179.5 degrees
    df_180_270 = pd.DataFrame(antenna_data[540:, :theta_range + 1])  # phi 180-270 degrees
    df_270_360 = pd.DataFrame(antenna_data[0:180, :theta_range + 1])  # phi 270-360 degrees
    df_180_360 = pd.concat([df_180_270, df_270_360])
    return pd.concat([df_180_360, df_0_180]).reset_index(drop=True)

    
def process_e_nrom(filename, theta_max, e_norm_path = '/data/curtin_eepnn/P_vogel_FEKO/enorms/E'):
    """
    Loads and processes e_norm data from a .mat file.

    This function reads an e_norm .mat file, extracts and rearranges the data for 
    either X-pol or Y-pol polarization, and returns a list of DataFrames, where 
    each DataFrame corresponds to the data for a specific antenna.

    Parameters:
    ----------
    filename : str
        The name of the .mat file containing the normalized electric field data.
        It should end with either 'Xpol_enorm.mat' or 'Ypol_enorm.mat'.
    
    theta_max : float
        The maximum theta value (in degrees) for which data is considered. 
        The function processes data within the range [0, theta_max].

    path : str, optional
        The directory path where the .mat file is located. 
        Default is '/data/curtin_eepnn/P_vogel_FEKO/enorms/E'.

    Returns:
    -------
    list of pandas.DataFrame
        A list of DataFrames where each DataFrame contains the processed e_norm 
        data for one antenna. The data is reorganized to match the desired angular 
        range and order. Row index is indicates the ϕ angle, wheras column name is θ angle.
    """
    x_range = np.arange(0, theta_max +0.5, 0.5)
    y_range = np.arange(-180, 180.5, 0.5)
    X, Y = np.meshgrid(x_range, y_range)
    theta_range = np.where(x_range == theta_max)[0][0]
    
    if e_norm_path:
        print(f"fetching e_norm data from {e_norm_path}...")
        e_norm = loadmat(f'{e_norm_path}/{filename}')['e_norm']
    else:
        print(f"loading {filename}...")
        e_norm = loadmat({filename})['e_norm']
        
    dfs = []
    if filename.endswith('Xpol_enorm.mat'):
        for antenna in range(256):
            df_0_180 = pd.DataFrame(e_norm[:360,:theta_range+1,antenna]) # get phi 0- 179.5 degree
            df_180_360 = pd.DataFrame(e_norm[360:,:theta_range+1,antenna]) # get phi 180- 360 degree
            df = pd.concat([df_180_360,df_0_180]).reset_index(drop = True) # 
            # df=process_xpol(e_norm)
            dfs.append(df)
    elif filename.endswith('Ypol_enorm.mat'):
        for antenna in range(256):
            df_0_180 = pd.DataFrame(e_norm[180:540,:theta_range+1,antenna]) # get phi 0- 179.5 degree
            df_180_270 = pd.DataFrame(e_norm[540:,:theta_range+1,antenna]) # get phi 180-270 degree
            df_270_360 = pd.DataFrame(e_norm[0:180,:theta_range+1,antenna]) # get phi 270-360 degree
            df_180_360 = pd.concat([df_180_270, df_270_360])
            df = pd.concat([df_180_360, df_0_180]).reset_index(drop = True)
            df=process_ypol(e_norm)
            # dfs.append(df)
    print("e_norm file/files were loaded and fliped to ")
    return dfs

def detect_bad_pattern(dfs):
    """
    Detects problematic patterns in a list of DataFrames based on threshold values.

    This function iterates through a list of DataFrames, identifying entries that fall below
    specific thresholds (-6dB < e_norm <= -3dB, -9dB < e_norm <= -6dB, e_norm<=-9). For each threshold, the function collects the
    row index, column index, and corresponding value for problematic entries.
    
    Parameters:
    ----------
    dfs : list of pandas.DataFrame
        A list of DataFrames representing data from 256 antennas. Each DataFrame should have
        numerical values and a consistent shape (721*181).
    
    Returns:
    -------
    tuple of lists
        - problematic_negative_3_all : list of lists
          Contains tuples of (ϕ,θ,e_norm) for entries ≤ -3.
        - problematic_negative_6_all : list of lists
          Contains tuples of (rϕ,θ,e_norm) for entries ≤ -6.
        - problematic_negative_9_all : list of lists
          Contains tuples of (ϕ,θ,e_norm) for entries ≤ -9.

    """

    problematic_negative_3 = -3
    problematic_negative_6 = -6
    problematic_negative_9 = -9
    problematic_negative_3_all = []
    problematic_negative_6_all = []
    problematic_negative_9_all = []

    for antenna in range(len(dfs)):
        # print(dfs[antenna].shape)
        df = dfs[antenna]
        container_problematic_negative_3 = []
        container_problematic_negative_6 = []
        container_problematic_negative_9 = []

        for index, row in df.iterrows():
              for theta in range(df.shape[1]):
                    if row[theta] <= problematic_negative_9:
                          container_problematic_negative_9.append((index,theta, df.loc[index, theta]))
                    elif row[theta] <= problematic_negative_6:
                          container_problematic_negative_6.append((index, theta, df.loc[index, theta]))
                    elif row[theta] <= problematic_negative_3:
                          container_problematic_negative_3.append((index, theta, df.loc[index, theta]))

        problematic_negative_3_all.append(container_problematic_negative_3)
        problematic_negative_6_all.append(container_problematic_negative_6)
        problematic_negative_9_all.append(container_problematic_negative_9)
    
    return problematic_negative_3_all, problematic_negative_6_all, problematic_negative_9_all



def cut_box_based_on_phi(pro_file):
    """
    identify the phi edge of boxes
    """
    phi_edges = []
    for antenna in range(len(pro_file)):
        antenna_problematic_list = pro_file[antenna].copy()
        if antenna_problematic_list:
            temp = []
            temp.append([antenna_problematic_list[0]])
            for index in range(len(antenna_problematic_list)-1):
                # get a tuple (351, 59, -3.9822506724878615), phi, theta,and power. angles need to ddevied by 2
                current_problematic_location = (antenna_problematic_list[index]) 
                next_problematic_location = (antenna_problematic_list[index+1])

                if  next_problematic_location[0] - current_problematic_location[0] >  1 :
                    # print(current_problematic_location, next_problematic_location)
                    temp[-1].append(current_problematic_location)
                    temp.append([next_problematic_location])

            temp[-1].append(antenna_problematic_list[-1])
            phi_edges.append(temp)
        else:
            phi_edges.append([])
    return phi_edges


def get_phi_range(file, pro_file, max_theta):
    all_box = []
    for i in range(len(file)):
        phi = (file[i][0][0], file[i][1][0])
        temp_dbi = []
        temp_theta = []
     
        for k in range(len(pro_file)):
            if pro_file[k][0] in range(phi[0],phi[1]+1):
                temp_theta.append(pro_file[k][1])
                temp_dbi.append (pro_file[k][2])

        theta = (np.min(temp_theta), np.max(temp_theta))

        all_box.append(([(i/2-180) for i in phi], [i/2 for i in theta]))

    return all_box


#get lowest dB of a problematic region
def get_lowest_dB(file, pro_file, max_theta):

    lowesest_dB_list = []
    for i in range(len(file)):
        phi = (file[i][0][0], file[i][1][0])

        temp_dbi = []
        temp_theta = []
        # print(len(pro_file))
        for k in range(len(pro_file)):
            if pro_file[k][0] in range(phi[0],phi[1]+1):
                temp_theta.append(pro_file[k][1])
                temp_dbi.append (pro_file[k][2])
        lowest_dB = np.min(temp_dbi)
        lowesest_dB_list.append(lowest_dB)
    return lowesest_dB_list


def transform_str_to_float(string):
    numbers = string.strip('[]').split(',')
    float_numbers = [float(num.strip()) for num in numbers]
    return float_numbers

def get_region_of_non_problematic(file):
    df = pd.read_csv(file)
    df['theta_range'] = df['theta_range'].apply(lambda x: transform_str_to_float(x))
    df['phi_range'] = df['phi_range'].apply(lambda x: transform_str_to_float(x))
    df['lowest_dB'] = df['lowest_dB'].apply(lambda x: float(x))
    df['theta_left'] = df['theta_range'].apply(lambda x: x[0])
    print(df.shape)
    print(df['theta_left'])
    return df

def detect_problematic_regions(fov, SAVEPATH):
    start_time = time.time()
    output_path = f'{SAVEPATH}/math_model_output.csv'
    with open(output_path, 'w') as f:
        writer = csv.writer(f)
        writer.writerow(['class','theta_range', 'phi_range', 'antenna', 'freq.', 'pol.', 'FOV','lowest_dB'])
        dirs = [os.listdir(SAVEPATH)]
        for obj in dirs[0]:
         
            if obj.endswith('enorm.mat'):
                filename = obj
                print(filename)
                freq, pol = match_freq_pol(filename)
                e_norm = process_e_nrom(filename, fov)
                pro_files = detect_bad_pattern(e_norm)

                for c in range(len(pro_files)):
                    all_problematic_regions = cut_box_based_on_phi(pro_files[c])
                    for antenna in range(len(all_problematic_regions)):
                        test = get_phi_range(all_problematic_regions[antenna], pro_files[c][antenna],fov)
                        lowest_dB_list = get_lowest_dB(all_problematic_regions[antenna], pro_files[c][antenna],fov)
                        for box in range(len(test)):
                            c = c
                            antenna_code = antenna
                            theta_range= test[box][1]
                            phi_range= test[box][0]
                            # print(f"test[box]: {test[box]}, len(test[box]): {len(test[box])},  phi_range: {phi_range}")
                            lowest_dB = lowest_dB_list[box]
                            writer.writerow([c,theta_range, phi_range, antenna_code, freq, pol, fov, lowest_dB])
   
    end_time = time.time()
    running_time = end_time - start_time
    df = get_region_of_non_problematic(output_path)
    df.to_excel(f'{SAVEPATH}/problematic_regions.xlsx',index=False)
    print("Total running time:", running_time/60, "minutes")
    
def plot_2d_contour(filename, theta_max, output_path=None, enorm_path= None):
    x_range = np.arange(0, theta_max +0.5, 0.5)
    y_range = np.arange(-180, 180.5, 0.5)
    X, Y = np.meshgrid(x_range, y_range)
    theta_range = np.where(x_range == theta_max)[0][0]
    e_norm = loadmat(f'{enorm_path}/{filename}')['e_norm']
    freq, pol = match_freq_pol(filename)

    for antenna in range(256):
        if filename.endswith('Xpol_enorm.mat'):
            df_0_180 = pd.DataFrame(e_norm[:360,:theta_range+1,antenna]) # get 0- 179.5 degree
            df_180_360 = pd.DataFrame(e_norm[360:,:theta_range+1,antenna]) # get 180- 360 degree
            df = pd.concat([df_180_360,df_0_180]).reset_index(drop = True) # df_b
            # df=process_xpol(e_norm)
            
        elif filename.endswith('Ypol_enorm.mat'):
            df_right = pd.DataFrame(e_norm[180:540,:theta_range+1,antenna]) # get 0- 179.5 degree
            df_left_1 = pd.DataFrame(e_norm[540:,:theta_range+1,antenna]) # get 180- 360 degree
            df_left_2 = pd.DataFrame(e_norm[0:180,:theta_range+1,antenna]) # get 180- 360 degree
            df_left = pd.concat([df_left_1, df_left_2])
            df = pd.concat([df_left, df_right])
            # df=process_ypol(e_norm)

        # Plot the contour
        plt.figure(figsize = (6.4,5))
        levels = [min(df.values.flatten()),-9, -6, -3]
        contour_plot = plt.contourf(X, Y, df.values,levels=levels,cmap= 'coolwarm',vmin =-12)
        c = plt.contour(X, Y, df.values,levels, linestyles='dashed', colors = 'k', alpha = .8)
        plt.clabel(c, inline=True, fontsize=8, colors = 'k')
        colorbar = plt.colorbar(contour_plot, label='Normalised EEPs, dB')

        # Set labels and title
        title = f"Problematic regions at {freq}MHz in {pol}pol, antenna #{antenna+1}"
        plt.title(title)
        plt.xlabel('(θ deg)')
        plt.ylabel('(φ deg)')
        # plt.ylim(-180,180)
        plt.grid()
        plt.tight_layout()
        if output_path is not None:
            plt.savefig(f"{output_path}/{freq}MHz_{pol}pol_#{antenna+1}.png", dpi=100)
        plt.close()