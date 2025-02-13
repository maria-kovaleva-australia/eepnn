import os
from scipy.io import loadmat, savemat
import pandas as pd
import re
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import math
from math import asin, degrees, sqrt
import csv
import warnings
from tqdm import tqdm
warnings.filterwarnings("ignore", message="The input coordinates to pcolormesh.*")


def match_freq_pol(file_name):
    pattern = r'(\d+(?:\.\d+)?)MHz_([XY])pol'
    match  = re.search(pattern, file_name)
    if match:
        frequency = match.group(1)
        pol = match.group(2)  
        # pol = 0 if pol == 'X' else 1
    return frequency, pol

# def count_e_tot(dim1, dim2, dim3, input_mat_file):
#     E_total= np.zeros((dim1, dim2, dim3), dtype=float)
#     E_phi = input_mat_file['Ephi']
#     E_theta = input_mat_file['Etheta']
#     E_total= np.sqrt(np.abs(E_phi)**2 + np.abs(E_theta)**2)
#     return E_total

def count_e_tot(dim1, dim2, dim3, E_phi, E_theta):
    return np.sqrt(np.abs(E_phi)**2 + np.abs(E_theta)**2)

# def count_e_norm(dim1, dim2, dim3,e_tot):
#     E_total_max_matrix = np.zeros((dim1,dim2, dim3), dtype = float)
#     E_norm = np.zeros((dim1, dim2, dim3), dtype=float)
#     for antenna in range(dim3):
#         E_total_vector_of_antenna = e_tot[:,:,antenna]
#         e_tot_m = np.max(np.max(E_total_vector_of_antenna))
#         E_total_max_matrix[:,:,antenna] = e_tot_m
#     E_norm = 20 * np.log10(e_tot/ E_total_max_matrix)
#     return E_norm

def count_e_norm(dim1, dim2, dim3, e_tot):
    E_total_max_matrix = np.max(e_tot, axis=(0, 1))  # Max across the first two dimensions (dim1, dim2)
    E_norm = 20 * np.log10(e_tot / E_total_max_matrix)  # Normalize in-place
    return E_norm

def detect_shape_of_data(source_files, starstwith='FEKO', ant_start=1, ant_end=256):
    
    files = os.listdir(source_files)
    mat_files = [file for file in files if file.endswith('.mat') and file.startswith(starstwith)]
    data = loadmat(f'{source_files}/{mat_files[0]}', variable_names=['Ephi'])
    dim1, dim2, dim3 = data['Ephi'].shape
    
    print(f"    -. Number of ϕ angles: {dim1}, Number of θ angles: {dim2}, Total antennas: {dim3}")
    
    if ant_end > dim3:
        tot_ant = dim3-ant_start+1
        ant_end = dim3
    
    tot_ant =ant_end-ant_start+1
    print(f"    -. Normalised power calculated for {dim3} antennas")
    return dim1, dim2, dim3
    
    
def cal_and_save_e_norm(dim1, dim2, dim3, source_files, save_path,starstwith='FEKO'):
    os.makedirs(save_path, exist_ok=True)
    enorm_list =[]
    files = os.listdir(source_files)
    mat_files = [file for file in files if file.endswith('.mat') and file.startswith(starstwith)]
    os.makedirs(f'{save_path}/e_norms', exist_ok=True)
    
    for file in tqdm(mat_files, desc="    -. Progress: ", unit="file"):
    # for file in mat_files:
        try:
            mat_file = loadmat(f'{source_files}/{file}', variable_names=['Ephi','Etheta'])
            e_tot = count_e_tot(dim1, dim2, dim3, mat_file['Ephi'], mat_file['Etheta'])
            e_norm = count_e_norm(dim1, dim2, dim3, e_tot)

            pattern = r'\d+(\.\d+)?MHz_[XY]pol'
            file_name = re.search(pattern, file).group()
            save_to = f'{save_path}/e_norms/{file_name}_enorm.mat'

            savemat(save_to, {'e_norm': e_norm})
            enorm_list.append(e_norm)
        except ValueError as err:
            print(f'{err} at {file}')
    
    return enorm_list

def process_xpol(antenna_data, theta_range, ant):
    """Processes data for Xpol files."""
    # Using numpy arrays directly instead of pandas DataFrames
    data_0_180 = antenna_data[:360, :theta_range + 1, ant]  # phi 0-179.5 degrees
    data_180_360 = antenna_data[360:, :theta_range + 1, ant]  # phi 180-360 degrees
    
    # Concatenate the arrays (flip order to ensure 180 is on top)
    return np.concatenate((data_180_360, data_0_180), axis=0)

def process_ypol(antenna_data, theta_range, ant):
    """Processes data for Ypol files."""
    # Using numpy arrays directly instead of pandas DataFrames
    data_0_180 = antenna_data[180:540, :theta_range + 1, ant]  # phi 0-179.5 degrees
    data_180_270 = antenna_data[540:, :theta_range + 1, ant]  # phi 180-270 degrees
    data_270_360 = antenna_data[0:180, :theta_range + 1, ant]  # phi 270-360 degrees
    
    # Concatenate and reorder the arrays for proper contour plot representation
    data_180_360 = np.concatenate((data_180_270, data_270_360), axis=0)
    return np.concatenate((data_180_360, data_0_180), axis=0)


def process_e_norm(e_norm_path, e_norm_filename, fov, ant_start=1, ant_end=256):
    """
    Loads and processes e_norm data from a .mat file.

    This function reads an e_norm .mat file, extracts and rearranges the data for 
    either X-pol or Y-pol polarization, and returns a list of arrays, where 
    each array corresponds to the data for a specific antenna.

    Parameters:
    ----------
    e_norm_filename : str
        The name of the .mat file containing the normalized electric field data.
        It should end with either 'Xpol_enorm.mat' or 'Ypol_enorm.mat'.
    
    fov : float
        The maximum theta value (in degrees) for which data is considered. 
        The function processes data within the range [0, fov].
        
    ant_start : int, optional
        The starting antenna index (1-based). Default is 1.
    
    ant_end : int, optional
        The ending antenna index (1-based). Default is 256.
        
    e_norm_path : str, optional
        The directory path where the .mat file is located. Default is 
        '/data/curtin_eepnn/P_vogel_FEKO/enorms/E'.

    Returns:
    -------
    list of np.array
        A list of arrays where each array contains the processed e_norm 
        data for one antenna. The data is reorganized to match the desired angular 
        range and order. The rows represent ϕ (phi) angles, and the columns represent 
        θ (theta) angles.
    """
    # Define the angular ranges
    x_range = np.arange(0, fov + 0.5, 0.5)
    theta_range = np.where(x_range == fov)[0][0]
    
    # Load e_norm data
    if e_norm_path:
        e_norm = loadmat(f'{e_norm_path}/{e_norm_filename}')['e_norm']
    else:
        e_norm = loadmat(e_norm_filename)['e_norm']
    
    # Initialize the list of arrays to store results
    results = []
    ant_start_idx = ant_start - 1  # Adjust for 0-based indexing
    
    # Scenario 1: When user wants to select a single antenna
    if ant_end == ant_start:
        antennas = [ant_start_idx]
    
    # Scenario 2: When user wants to select a range of antennas
    elif ant_end > ant_start:
        antennas = range(ant_start_idx, ant_end)
    
    # Process data for each antenna in the selected range
    for antenna in antennas:
        if e_norm_filename.endswith('Xpol_enorm.mat'):
            result = process_xpol(e_norm, theta_range, antenna)
        elif e_norm_filename.endswith('Ypol_enorm.mat'):
            result = process_ypol(e_norm, theta_range, antenna)
        results.append(result)
    
    return results


def detect_bad_pattern(results, problematic_threshold=-3):
    """
    Detects problematic patterns in a list of arrays based on threshold values.

    This function iterates through a list of arrays, identifying entries that fall below
    specific thresholds (-6dB < e_norm <= problematic_threshold). For each threshold, the function collects the
    row index, column index, and corresponding value for problematic entries.
    
    Parameters:
    ----------
    results : list of np.array
        A list of arrays representing data from 256 antennas. Each array should have
        numerical values and a consistent shape (721*181).
    
    Returns:
    -------
    tuple of lists
        - PRC_ALL : list of lists
          Contains tuples of (ϕ,θ,e_norm) for entries ≤ -3.
    """
    PRC_ALL = []  # Initialize the list to hold all problematic regions

    for antenna in range(len(results)):
        result = results[antenna]
        pcr_ant = []
        # Iterate through the array for each φ and θ
        for phi in range(result.shape[0]):  # Rows represent φ (phi)
            for theta in range(result.shape[1]):  # Columns represent θ (theta)
                if result[phi, theta] <= problematic_threshold:
                    pcr_ant.append((phi, theta, result[phi, theta]))  # Store (φ, θ, e_norm)
        PRC_ALL.append(pcr_ant)
    
    return PRC_ALL


def identify_phi_edges(problematic_data):
    """
    Identifies the φ (phi) edges of problematic regions for each antenna.

    This function processes a list of problematic regions, represented as a list of tuples
    containing φ, θ, and power values, and groups them into contiguous regions
    based on the φ angle. Non-contiguous regions are separated into distinct groups.

    Parameters:
    ----------
    problematic_data : list of list of tuples
        A list where each element represents a single antenna's problematic regions.
        Each problematic region is defined as a list of tuples (φ, θ, power).

    Returns:
    -------
    list of list of list of tuples
        A nested list where the first level corresponds to antennas, the second level
        contains groups of contiguous φ regions, and each group contains tuples
        representing (φ, θ, power).
    """
    phi_edges = []

    for antenna in problematic_data:
        if antenna:
            temp = [[antenna[0]]]  # Start a new group with the first problematic point
          
            for index in range(len(antenna) - 1):
              
                current = antenna[index]
                next_ = antenna[index + 1]
                
                # Check if the φ difference exceeds 1 (indicating a new group)
                if next_[0] - current[0] > 1:
                    temp[-1].append(current)  # Close the current group
                    temp.append([next_])     # Start a new group

            temp[-1].append(antenna[-1])  # Add the last point to the last group
            phi_edges.append(temp)
        else:
            phi_edges.append([])

    return phi_edges

def calculate_phi_theta_ranges(phi_groups, problematic_data):
    """
    Calculates φ (phi) and θ (theta) ranges for problematic regions.

    This function determines the angular boundaries (φ and θ) for each identified
    problematic region and converts the angles from integer indices to degrees.

    Parameters:
    ----------
    phi_groups : list of list of tuples
        A list of grouped φ regions for each antenna. Each group contains tuples
        of (φ, θ, power).

    problematic_data : list of tuples
        A flat list of tuples representing all problematic points for a given antenna,
        where each tuple contains (φ, θ, power).

    Returns:
    -------
    list of tuples
        A list where each element is a tuple:
        - First element: φ range as (φ_min, φ_max) in degrees.
        - Second element: θ range as (θ_min, θ_max) in degrees.
    """
    all_ranges = []

    for group in phi_groups:
        # Extract φ range
        phi_min, phi_max = group[0][0], group[-1][0]

        # Find corresponding θ and power values within the φ range
        temp_theta = [p[1] for p in problematic_data if phi_min <= p[0] <= phi_max]

        # Calculate θ range
        theta_min, theta_max = np.min(temp_theta), np.max(temp_theta)

        # Convert ranges to degrees
        phi_range = [(phi_min / 2) - 180, (phi_max / 2) - 180]
        theta_range = [theta_min / 2, theta_max / 2]

        all_ranges.append((phi_range, theta_range))
    return all_ranges



#get lowest dB of a problematic region
def get_minimum_power_dB(file, problematic_file):
    """
    this is to get the minimum power of the problematic region.
    """

    min_dB_list = []
    for i in range(len(file)):
        phi = (file[i][0][0], file[i][1][0])

        temp_dbi = []
        temp_theta = []
        # print(len(problematic_file))
        for k in range(len(problematic_file)):
            if problematic_file[k][0] in range(phi[0],phi[1]+1):
                temp_theta.append(problematic_file[k][1])
                temp_dbi.append (problematic_file[k][2])
        min_dB = np.min(temp_dbi)
        min_dB_list.append(min_dB)
    return min_dB_list

def plot_2d_eep(filename, fov, output_path=None, enorm_path= None, problematic_threshold=None):
    x_range = np.arange(0, fov +0.5, 0.5)
    y_range = np.arange(-180, 180.5, 0.5)
    X, Y = np.meshgrid(x_range, y_range)
    theta_range = np.where(x_range == fov)[0][0]
    e_norm = loadmat(f'{enorm_path}/{filename}')['e_norm']
    freq, pol = match_freq_pol(filename)
    
    ant_num = e_norm.shape[2]
    # print(f"shape of enorm, {e_norm.shape}, total antenna number is {ant_num}")
    for antenna in range(ant_num):
        if filename.endswith('Xpol_enorm.mat'):
            arr=process_xpol(e_norm, theta_range,antenna)
        elif filename.endswith('Ypol_enorm.mat'):
            arr=process_ypol(e_norm, theta_range,antenna)

        locat = np.unravel_index(np.argmax(arr), arr.shape)
        plt.scatter(locat[1]/2, (locat[0] / 2)-180, c="gray", marker='+', s=80, linewidths=1.5)
        # print((locat[0] / 2)-180,locat[1]/2 )
        plt.imshow(arr, aspect = 'auto', extent=[0,90,-180,180], alpha = 1, origin = "lower", cmap= 'viridis')
        plt.colorbar()
        
        # Set labels and title
        title = f"Polar coordinate: {freq}MHz in {pol}pol, antenna #{antenna+1}"
        plt.title(title)
        plt.xlabel('(θ deg)')
        plt.ylabel('(φ deg)')
        plt.grid()
        plt.tight_layout()
        if output_path is not None:
            plt.savefig(f"{output_path}/{freq}MHz_{pol}pol_#{antenna+1}.png", dpi=100)
        plt.show()
        plt.close()



def get_kx_ky(FEKO_data_path):
    # all_mat= os.listdir(FEKO_data_path)
    all_mat = [file for file in os.listdir(FEKO_data_path) if file.endswith('.mat')]
    # Load the .mat file once
    data = loadmat(f"{FEKO_data_path}/{all_mat[0]}", variable_names=['kx', 'ky'])
    return data['kx'], data['ky'] 

def plot_uv_plane(eep, kx, ky, freq, pol, antenna, problematic_threshold, locat):
    mesh = plt.pcolormesh(kx, ky, eep, shading='nearest', cmap='viridis', vmin=np.min(eep), vmax=problematic_threshold)  # Heatmap shading='auto'
    mesh.set_clim(np.min(eep), problematic_threshold)
    mesh.cmap.set_over('w')
    cbar = plt.colorbar(mesh, label='Normalised power (dB)')
    cbar.set_ticks(np.linspace(np.min(eep), problematic_threshold, num=5))  # Set proper tick marks
    
    # Add contour line where eep == problematic_threshold
    contour = plt.contour(kx, ky, eep, levels=[problematic_threshold], colors='k', linewidths=1.2)
    plt.clabel(contour, inline=True, fontsize=8, fmt=f"{problematic_threshold} dB", colors='k')

    plt.xlabel('u (kx)')
    plt.ylabel('v (ky)')
    plt.xticks(np.arange(-1, 1.2, 0.25))
    plt.yticks(np.arange(-1, 1.2, 0.25))
    plt.axis('equal')
    title = f"{freq}MHz in {pol}pol, antenna #{antenna}"
    plt.title(title)
    
    
    # Find the maximum power location and Plot the marker for maximum power location
    max_power_idx = np.unravel_index(np.argmax(eep), eep.shape)
    max_kx = kx[max_power_idx[0]][max_power_idx[1]]  # Index corresponds to ky row, kx column
    max_ky = ky[max_power_idx[0]][max_power_idx[1]]  # Index corresponds to ky row, kx column
    
    plt.scatter(max_kx, max_ky, color='k', marker='+', s=50, zorder=5)
    coords_text = f"({locat[1] / 2:.1f}\u00B0, {(locat[0] / 2) - 180:.1f}\u00B0)"
    plt.annotate(coords_text,         #f"(u: {max_kx:.3f}, v: {max_ky:.3f})", 
                 xy=(max_kx, max_ky), 
                 xytext=(3, 3),       # Slight offset for readability
                 textcoords='offset points', #axes fraction
                 fontsize=6, 
                 fontweight='bold',
                 color='k')
    
    # add phi=0
    add_phi0_coor(pol)
    plt.grid()
    
def plot_max_power_location(locat):
    plt.scatter(locat[1]/2, (locat[0] / 2)-180, c="k", marker='+', s=50, linewidths=1.2) # plot location of max power
    # Annotate the exact coordinates
    coords_text = f"({locat[1] / 2:.1f}\u00B0, {(locat[0] / 2) - 180:.1f}\u00B0)"
    plt.annotate(coords_text, 
                 xy=(locat[1] / 2, (locat[0] / 2) - 180), 
                 xytext=(locat[1] / 2 + 2, (locat[0] / 2) - 180 + 10),  # Slight offset for readability
                 fontsize=6, 
                 fontweight='bold',
                 color='k',
                 ha='left',
                 va='top')
    
def plot_polar_coor(X, Y, locat, arr, problematic_threshold, freq, pol, antenna):
    levels = sorted([min(arr.flatten()), problematic_threshold - 10, problematic_threshold - 5, problematic_threshold])
    contour_plot = plt.contourf(X, Y, arr, levels=levels, cmap='viridis', vmin=problematic_threshold - 20)  # Use arr directly as it's a numpy array
    c = plt.contour(X, Y, arr, levels, linestyles='dashed', colors='k', alpha=.8)
    plt.clabel(c, inline=True, fontsize=8, colors='k')
    colorbar = plt.colorbar(contour_plot, label='Normalised power (dB)')
    plot_max_power_location(locat)

    # Set title and labels
    title = f"{freq}MHz in {pol}pol, antenna #{antenna}"
    plt.title(title)
    plt.xlabel('(θ deg)')
    plt.ylabel('(φ deg)')
    plt.grid()

        

def add_phi0_coor(pol):
    # Add phi = 0° annotation
    if pol.lower() == 'x':  # X-pol: φ = 0° at (u=1, v=0)
        phi0_u, phi0_v = 1, 0
    elif pol.lower() == 'y':  # Y-pol: φ = 0° at (u=0, v=1)
        phi0_u, phi0_v = 0, 1
    else:
        phi0_u, phi0_v = None, None  # Unknown polarization

    if phi0_u is not None and phi0_v is not None:
        # plt.scatter(phi0_u, phi0_v, color='r', marker='o', s=50, zorder=5)
        plt.annotate(r"$φ=0^\circ$", 
                     xy=(phi0_u, phi0_v), 
                     xytext=(1, 2), 
                     textcoords='offset points', 
                     fontsize=8,  
                     fontweight='bold',
                     color='k',
                     ha='center',
                     va='center')


def plot_it(enorm_folder, e_norm_filename, kx, ky, output_path, problematic_threshold, ant_start=1, ant_end=256):
    ################ process data #########################
    fov=90
    x_range = np.arange(0, fov +0.5, 0.5)
    y_range = np.arange(-180, 180.5, 0.5)
    X, Y = np.meshgrid(x_range, y_range)
    theta_range = np.where(x_range == fov)[0][0]
    
    e_norm = loadmat(f'{enorm_folder}/{e_norm_filename}')['e_norm']
    freq, pol = match_freq_pol(e_norm_filename)
    
    ant_start_idx = ant_start - 1  # Adjust for 0-based indexing
    # Scenario 1: When user wants to select a single antenna
    if ant_end == ant_start:
        antennas = [ant_start_idx]
    # Scenario 2: When user wants to select a range of antennas
    elif ant_end > ant_start:
        antennas = range(ant_start_idx, ant_end)
    
    ########################### plotting ####################
    for antenna in antennas:
        if e_norm_filename.endswith('Xpol_enorm.mat'):
            arr=process_xpol(e_norm, theta_range,antenna)
        elif e_norm_filename.endswith('Ypol_enorm.mat'):
            arr=process_ypol(e_norm, theta_range,antenna)
   
        # Plot the contour
        locat = np.unravel_index(np.argmax(arr), arr.shape)
        plt.figure(figsize = (13,5))
        plt.subplot(121)
        plot_polar_coor(X,Y, locat, arr, problematic_threshold,freq, pol, antenna+1)
        
        # plot uv
        plt.subplot(122)
        plot_uv_plane(e_norm[:, :, antenna], kx, ky, freq, pol, antenna + 1, problematic_threshold, locat)
        
        plt.tight_layout()
    
        ######################## save and show plots ######################
        if output_path is not None:
            os.makedirs(f"{output_path}/plots", exist_ok = True)
            plt.savefig(f"{output_path}/plots/{freq}MHz_{pol}pol_#{antenna+1}_{problematic_threshold}dB.png", dpi=100)
        # plt.show()
        plt.close()