import functions as F
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
import os
import csv
import time
import argparse


def main(fov, FEKO_SOURCE_PATH, SAVE_PATH, problematic_threshold = -3,ant_start=1, ant_end=256):
    start_time = time.time()
    output_path_csv = f'{SAVE_PATH}/problematic_regions.csv'
    output_path_excel=f'{SAVE_PATH}/problematic_regions.xlsx'
    # F.cal_and_save_e_nrom(721,181,256, FEKO_SOURCE_PATH, SAVE_PATH)
    
    with open(output_path_csv, 'w') as f:
        writer = csv.writer(f)
        writer.writerow(['class','theta_range', 'phi_range', 'antenna', 'freq.', 'pol.', 'FOV',
                         'minimum_dB_in_region', 'ant_max_power', 'locat_max_power'])
        e_norm_files = [os.listdir(f"{SAVE_PATH}/e_norms")][0]
        print(e_norm_files)
        for file_ in e_norm_files:
            if file_.endswith('enorm.mat'):
                print(f"processing {file_}")
                freq, pol = F.match_freq_pol(file_)
                e_norm = F.process_e_norm(f"{SAVE_PATH}/e_norms", file_, fov, ant_start=ant_start, ant_end=ant_end)
                print(f"selected {len(e_norm)} EEPs")
                print(f"type of enorm: {type(e_norm[0])}")
               
           
                PR_container = F.detect_bad_pattern(e_norm, problematic_threshold=problematic_threshold) # output a list(len is antenna number) of list(len is number of PR) of tuples (ϕ,θ,e_norm)
                print(f"Detected {len(PR_container)} problematic region {problematic_threshold}")
  
                print("________________________________tested above ____________________________________________")
    
                # for c in range(len(PR_container)):
                PRs = F.identify_phi_edges(PR_container)
                for antenna in range(len(PRs)):
                    ant_num = antenna+ant_start
                    print(f"ant_num: {ant_num}")
                    
                    print(f"testing F.calculate_phi_theta_ranges....")
                    boxes = F.calculate_phi_theta_ranges(PRs[antenna], PR_container[antenna])
                    
                    print(f"testing F.get_minimum_power_dB....")
                    lowest_dB_list = F.get_minimum_power_dB(PRs[antenna], PR_container[antenna])
                    ant_max_power = e_norm[antenna].max().max()
                    
                    print(f"printing e_norm for antenna: {antenna}, type: {type(e_norm[antenna])}")
                    print( e_norm[antenna])
                    p, t= divmod(e_norm[antenna].values.argmax(), e_norm[antenna].shape[1])
                    location_max_power_phi = (p / 2) - 180
                    location_max_power_the = t/2
                    for box in range(len(boxes)):
                        theta_range= boxes[box][1]
                        phi_range= boxes[box][0]
                        # print(f"boxes[box]: {boxes[box]}, len(boxes[box]): {len(boxes[box])},  phi_range: {phi_range}")
                        minimum_dB_in_region = lowest_dB_list[box]
                        writer.writerow([problematic_threshold,
                                         theta_range, 
                                         phi_range, 
                                         ant_num, 
                                         freq,
                                         pol,
                                         fov,
                                         minimum_dB_in_region,
                                         ant_max_power,
                                         (location_max_power_phi, location_max_power_the)])

                # for c in range(len(PR_container)):
                #     PRs = F.identify_phi_edges(PR_container[c])
                #     pr_class = (c+1)*-3 
                #     for antenna in range(len(PRs)):
                #         ant_num = antenna+1
                #         boxes = F.calculate_phi_theta_ranges(PRs[antenna], PR_container[c][antenna])
                #         lowest_dB_list = F.get_minimum_power_dB(PRs[antenna], PR_container[c][antenna])
                #         ant_max_power = e_norm[antenna].max().max()
                #         location_max_power_phi, location_max_power_the = divmod(e_norm[antenna].values.argmax(), e_norm[antenna].shape[1])
                #         for box in range(len(boxes)):
                #             theta_range= boxes[box][1]
                #             phi_range= boxes[box][0]
                #             # print(f"boxes[box]: {boxes[box]}, len(boxes[box]): {len(boxes[box])},  phi_range: {phi_range}")
                #             minimum_dB_in_region = lowest_dB_list[box]
                #             writer.writerow([pr_class,
                #                              theta_range, 
                #                              phi_range, 
                #                              ant_num, 
                #                              freq,
                #                              pol,
                #                              fov,
                #                              minimum_dB_in_region,
                #                              ant_max_power,
                #                              (location_max_power_phi, location_max_power_the)])
    
    end_time = time.time()
    running_time = end_time - start_time
    df = F.process_problematic_region_data(output_path_csv)
    df.to_excel(output_path_excel,index=False)
    print("Total running time:", running_time/60, "minutes")  


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = 'Detecting problematic regions for antenna within an antenna array layout')
    parser.add_argument('--fov', type = int, required = True, help = 'field of view')
    parser.add_argument('--FEKO_SOURCE_PATH', type = str, required = True, help = 'Path of FEKO .mat files')
    parser.add_argument('--SAVE_PATH', type = str, required = True, help = 'The path for saving calculated e_norm files and problematic regions result.')
    parser.add_argument('--problematic_threshold', type = float, required = True, help = 'The threshold of defining as problematic in logarithmic scale (dB), eg: -3')
    parser.add_argument('--ant_start', type = int, required = True, help = 'The start antenna number of interest')
    parser.add_argument('--ant_end', type = int, required = True, help = 'The end antenna number of interest')
    args = parser.parse_args()   
    main(args.fov, args.FEKO_SOURCE_PATH, args.SAVE_PATH, args.problematic_threshold , args.ant_start , args.ant_end)
