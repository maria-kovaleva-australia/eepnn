import functions as F
import pandas as pd
import time
import csv
import os
import argparse


def main(fov, FEKO_SOURCE_PATH, SAVE_PATH, problematic_threshold,ant_start, ant_end):
    start_time = time.time()
    os.makedirs(f'{SAVE_PATH}/result', exist_ok=True)
    output_path_csv = f'{SAVE_PATH}/result/problematic_regions.csv'
    output_path_excel=f'{SAVE_PATH}/result/problematic_regions.xlsx'
    result_path = f"{SAVE_PATH}/result"
    F.cal_and_save_e_nrom(721,181,256, FEKO_SOURCE_PATH,  result_path)
    
    with open(output_path_csv, 'w') as f:
        writer = csv.writer(f)
        writer.writerow(['class','theta_range', 'phi_range', 'antenna', 'freq.', 'pol.', 'FOV',
                         'minimum_dB_in_region', 'ant_max_power', 'locat_max_power'])
        e_norm_files = [os.listdir(f"{result_path}/e_norms")][0]
        print(e_norm_files)
        for file_ in e_norm_files:
            if file_.endswith('enorm.mat'):
                print(f"processing {file_}")
                freq, pol = F.match_freq_pol(file_)
                e_norm = F.process_e_norm(f"{result_path}/e_norms", file_, fov, ant_start=ant_start, ant_end=ant_end)
                print(f"selected {len(e_norm)} EEPs")
               
                # output a list(len is antenna number) of list(len is number of PR) of tuples (ϕ,θ,e_norm)
                PR_container = F.detect_bad_pattern(e_norm, problematic_threshold=problematic_threshold) 
             
                PRs = F.identify_phi_edges(PR_container)
                for antenna in range(len(PRs)):
                    ant_num = antenna+ant_start
                    boxes = F.calculate_phi_theta_ranges(PRs[antenna], PR_container[antenna])
                    lowest_dB_list = F.get_minimum_power_dB(PRs[antenna], PR_container[antenna])
                    ant_max_power = e_norm[antenna].max().max()

                    p, t= divmod(e_norm[antenna].values.argmax(), e_norm[antenna].shape[1])
                    p_adj = (p / 2) - 180
                    t_adj = t/2
                    for box in range(len(boxes)):
                        theta_range= boxes[box][1]
                        phi_range= boxes[box][0]
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
                                         (p_adj, t_adj)])

    
    end_time = time.time()
    running_time = end_time - start_time
    df = F.process_problematic_region_data(output_path_csv)
    df.to_excel(output_path_excel,index=False)
    print(f"saved result to {result_path} folder")
    print("Total running time:", running_time/60, "minutes")  



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Detect problematic regions for antennas within an array layout')
    parser.add_argument('--fov', type=int, default=45, help='Field of view (0-90 degrees, default: 45)')
    parser.add_argument('--FEKO_SOURCE_PATH', type=str, help='Path to FEKO .mat files')
    parser.add_argument('--SAVE_PATH', type=str, help='Path to an existing directory where results will be saved')
    parser.add_argument('--problematic_threshold', type=float, default=-3, help='Problematic threshold in dB (default: -3)')
    parser.add_argument('--ant_start', type=int, default=1, help='Start antenna number (default: 1)')
    parser.add_argument('--ant_end', type=int, default=256, help='End antenna number (default: 256)')
    
    args = parser.parse_args()
    main(args.fov, args.FEKO_SOURCE_PATH, args.SAVE_PATH, args.problematic_threshold, args.ant_start, args.ant_end)