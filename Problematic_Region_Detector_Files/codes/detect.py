import functions as F
import pandas as pd
import time
import csv
import os
import argparse
from tqdm import tqdm


################ Cal and save e norm #################

class EEPProblematicRegionProcessor:
    """
    A class for processing, analyzing, and detecting problematic regions in EEP (Embedded Element Pattern) data.

    This class provides functionality to:
    - Compute and save normalized EEP values.
    - Generate plots for EEP data in polar and UV planes.
    - Detect problematic regions in EEP data and save the results.

    Attributes:
        fov (float): Field of view in degrees.
        problematic_threshold (float): Threshold to define problematic regions.
        ant_start (int): Starting index of the antenna.
        ant_end (int): Ending index of the antenna.
        FEKO_SOURCE_PATH (str): Path to the FEKO simulation data.
        SAVE_PATH (str): Directory where results will be saved.
        result_path (str): Path where processed results are stored.
        output_path_csv (str): File path for saving detected problematic regions.
    """
    
    def __init__(self, fov, FEKO_SOURCE_PATH, SAVE_PATH, problematic_threshold,ant_start, ant_end):
        """
        Initializes the EEPProblematicRegionProcessor class with the given parameters.

        Args:
            fov (float): Field of view in degrees.
            FEKO_SOURCE_PATH (str): Path to the FEKO simulation data.
            SAVE_PATH (str): Directory where results will be saved.
            problematic_threshold (float): Threshold to define problematic regions.
            ant_start (int): Starting index of the antenna.
            ant_end (int): Ending index of the antenna.
        """
        # Shared attributes for processing
        self.fov=fov
        self.problematic_threshold=problematic_threshold
        self.ant_start=ant_start
        self.ant_end=ant_end
        
        # Shared attributes for saving
        self.FEKO_SOURCE_PATH=FEKO_SOURCE_PATH
        self.SAVE_PATH=SAVE_PATH
        self.result_path = f"{self.SAVE_PATH}/result"
        self.output_path_csv = f"{self.result_path}/problematic_regions{self.problematic_threshold}_fov{self.fov}_{time.time()}.csv"
        os.makedirs(self.result_path, exist_ok=True)

    def compute_and_store_e_norms(self):
        """
        Computes and saves normalized EEP values.

        This function:
        - Detects the shape of the FEKO dataset.
        - Computes and stores the EEP norm values.
        - Ensures the ending antenna index is within the dataset bounds.

        Prints:
            - Progress updates and total computation time.
        """
   
        # get shape of FEKO data
        s_time = time.time()
        print(f"Calculating EEPs in logarithmic scale...\n")
        dim1, dim2, dim3 = F.detect_shape_of_data(self.FEKO_SOURCE_PATH, ant_start=self.ant_start, ant_end=self.ant_end)
        F.cal_and_save_e_norm(dim1, dim2, dim3, self.FEKO_SOURCE_PATH, self.result_path)
        print(f'    -. Saved e_norm to {self.result_path}/e_norms')
  
         # Adapt to the actual last antenna
        if self.ant_end > dim3:
            print(f"    -. Adjusted ant_end from {self.ant_end} to {dim3}")
            self.ant_end = dim3  # Correctly updates the instance variable
        print(f"    -. Total time spent for this task: {(time.time()-s_time)/60:.2f} minutes.")

    
    def plot_(self):
        """
        Generates and saves plots of EEP data in polar and UV coordinate planes.

        This function:
        - Loads computed EEP norms.
        - Plots EEPs using `F.plot_it`.
        - Saves the generated plots to the results directory.

        Prints:
            - Progress updates and total time taken.
        """
        s_time = time.time()
        print(f"\nPlotting EEPs in polar coordinate and uv plane...\n")
        
        kx, ky = F.get_kx_ky(self.FEKO_SOURCE_PATH)
        e_norm_files = [f for f in os.listdir(f"{self.result_path}/e_norms") if f.endswith('enorm.mat')]

        for file_ in tqdm(e_norm_files, desc="    -. Progress: ", unit="file"):
            F.plot_it(f"{self.result_path}/e_norms", 
                      file_, 
                      kx, ky, 
                      output_path=self.result_path, 
                      problematic_threshold=self.problematic_threshold, 
                      ant_start=self.ant_start, ant_end=self.ant_end)

        print(f"    -. Plots include {self.ant_end - self.ant_start + 1} antenna(s), covering antennas {self.ant_start} to {self.ant_end}.")
        print(f'    -. Saved plots to {self.result_path}/plots') 
        print(f"    -. Total time spent for this task: {(time.time()-s_time)/60:.2f} minutes.")
    
    def detect_(self):
        """
        Detects problematic regions in EEP data and saves the results to a CSV file.

        This function:
        - Scans processed EEP data.
        - Identifies regions with EEP values below the `problematic_threshold`.
        - Saves detected problematic regions to a CSV file.

        Prints:
            - Progress updates and total time taken.
        """
        s_time = time.time()
        
        print(f"\nDetecting problematic regions <= {self.problematic_threshold} in the FOV {self.fov}\u00B0")
        print(f"    -. Result included {self.ant_end - self.ant_start +1} antenna/antennas, including antenna {self.ant_start} to {self.ant_end}")
        
        with open(self.output_path_csv, 'w') as f:
            writer = csv.writer(f)
            writer.writerow(['threshold','Theta_range', 'phi_range', 'antenna', 'freq.', 'pol.', 'FOV',
                             'minimum_dB_in_region', 'ant_max_power', 'max_power_coords_in_fov'])
             
            e_norm_files = [f for f in os.listdir(f"{self.result_path}/e_norms") if f.endswith('enorm.mat')]
            
            for file_ in tqdm(e_norm_files, desc="    -. Progress: ", unit="file"):
                    freq, pol = F.match_freq_pol(file_)
                    e_norm = F.process_e_norm(f"{self.result_path}/e_norms", file_, self.fov, ant_start=self.ant_start, ant_end=self.ant_end)
                    PR_container = F.detect_bad_pattern(e_norm, problematic_threshold=self.problematic_threshold) 

                    PRs = F.identify_phi_edges(PR_container)
                    for antenna in range(len(PRs)):
                        ant_num = antenna+self.ant_start
                        boxes = F.calculate_phi_theta_ranges(PRs[antenna], PR_container[antenna])
                        lowest_dB_list = F.get_minimum_power_dB(PRs[antenna], PR_container[antenna])
                        ant_max_power = e_norm[antenna].max().max()

                        p, t= divmod(e_norm[antenna].argmax(), e_norm[antenna].shape[1]) #divmod(e_norm[antenna].values.argmax(), e_norm[antenna].shape[1])
                        p_adj = (p / 2) - 180
                        t_adj = t/2
                        for box in range(len(boxes)):
                            theta_range= boxes[box][1]
                            phi_range= boxes[box][0]
                            minimum_dB_in_region = lowest_dB_list[box]
                            writer.writerow([self.problematic_threshold,
                                             theta_range, 
                                             phi_range, 
                                             ant_num, 
                                             freq,
                                             pol,
                                             self.fov,
                                             minimum_dB_in_region,
                                             ant_max_power,
                                             (t_adj, p_adj)])

        print(f"    -. Saved result to {self.output_path_csv}")
        print(f"    -. Total time spent for this task: {(time.time()-s_time)/60:.2f} minutes.")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Detect problematic regions for antennas within an array layout")
    parser.add_argument("--fov", type=int, default=45, help="Field of view (0-90 degrees, default: 45)")
    parser.add_argument("--FEKO_SOURCE_PATH", type=str, required=True, help="Path to FEKO .mat files")
    parser.add_argument("--SAVE_PATH", type=str, required=True, help="Path to an existing directory where results will be saved")
    parser.add_argument("--problematic_threshold", type=float, default=-3, help="Problematic threshold in dB (default: -3)")
    parser.add_argument("--ant_start", type=int, default=1, help="Start antenna number (default: 1)")
    parser.add_argument("--ant_end", type=int, default=256, help="End antenna number (default: 256) of interest")
    parser.add_argument("--compute_enorms", type=int, choices=[0, 1], default=0, 
                        help="Set to 1 to compute and store E-norms, 0 to disable (default: 0)")
    parser.add_argument("--plot_EEPs", type=int, choices=[0, 1], default=0, 
                        help="Set to 1 to plot and store EEPs, 0 to disable (default: 0)")

    args = parser.parse_args()

    # Create an instance of the processor
    processor = EEPProblematicRegionProcessor(fov=args.fov,
                                              FEKO_SOURCE_PATH=args.FEKO_SOURCE_PATH,
                                              SAVE_PATH=args.SAVE_PATH,
                                              problematic_threshold=args.problematic_threshold,
                                              ant_start=args.ant_start,
                                              ant_end=args.ant_end)

    # Conditionally compute and store E-norms and plot EEPs
    if args.compute_enorms:
        processor.compute_and_store_e_norms()
        
    if args.plot_EEPs:
        processor.plot_()
    
    processor.detect_()
