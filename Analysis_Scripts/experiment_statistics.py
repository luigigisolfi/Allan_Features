# %%
"""
Script to perform statistical analysis on PRIDE Doppler FDETS files.

This script:
- Loads standard SPICE kernels.
- Initializes PRIDE characterization tools.
- Extracts Doppler noise, SNR, and elevation data from FDETS input files.
- Generates user-defined parameter plots (SNR, Doppler noise, raw FDETS).
- Generates elevation plots for each station and the whole experiment.
- Computes statistical summaries (mean, std) of SNR and Doppler noise.
- Combines SNR and elevation plots into single images for each station/date.

This specific demo is configured to analyze:
- Mission: 'juice'
- Experiment: 'ec094b'

Note:
This code has to be run before being able to run dataset_statistics.py .
"""
from pride_characterization_library import PrideDopplerCharacterization
import os
from collections import defaultdict
from datetime import timezone
import datetime
import shutil

pride = PrideDopplerCharacterization() # create PRIDE Object
process_fdets = pride.ProcessFdets() # create Process Fdets Object
utilities = pride.Utilities() # create Utilities Object
analysis = pride.Analysis(process_fdets, utilities) # create Analysis Object

TWO_STEP_FILTER_FLAG = False # if TRUE, allows to save ONLY the filtered data to the user_defined_parameters output files. If FALSE, all original data points are saved.
ALLAN_DEVIATIONS_FLAG = True
PLOT_GAUSSIAN_FLAG = False
COMPARE_FILTERS_FLAG = False
BAD_OBSERVATIONS_FLAG = False # if set to true, it 1) plots the observations as flagged and 2) removes them from the final statistics for mean FoM computation

start_date = datetime.datetime(2000, 1, 1, tzinfo=timezone.utc)
end_date =  datetime.datetime(2024, 12, 31, tzinfo=timezone.utc)
yymm_folders_to_consider = utilities.list_yymm(start_date, end_date)

months_list = list(yymm_folders_to_consider.keys())
days_list = [item for sublist in yymm_folders_to_consider.values() for item in sublist]

root_dir = f'/Users/lgisolfi/Desktop/PRIDE_DATA_NEW/analysed_pride_data/'
experiments_to_analyze = defaultdict(list)
missions_to_analyse = ['juice']
yymmdd_folders_per_mission = defaultdict(list)

for mission_name in missions_to_analyse:
    mission_root = os.path.join(root_dir, mission_name)
    if not os.path.exists(mission_root):
        print(f"⚠️  Warning: Mission root path '{mission_root}' does not exist.")
        continue

    for yymm_folder in months_list:
        month_folder_name = f"{mission_name}_{yymm_folder}"
        month_folder_path = os.path.join(mission_root, month_folder_name)

        if not os.path.exists(month_folder_path):
            continue  # continue to next yymmdd_folder

        for yymmdd in days_list:
            day_folder_name = f"{mission_name}_{yymmdd}"
            day_folder_path = os.path.join(month_folder_path, day_folder_name)

            if not os.path.exists(day_folder_path):
                continue  # continue to next yymmdd_folder
            yymmdd_folders_per_mission[mission_name].append(yymmdd)

if not yymmdd_folders_per_mission:
    print(f'No files found between {start_date} and {end_date} for mission: {mission_name.upper()}.\nAborting...\n')
    exit()
experiments_to_analyze = yymmdd_folders_per_mission
for mission_name, yymmdds in yymmdd_folders_per_mission.items():
    yymmdd_folders = [mission_name + '_' + yymmdd for yymmdd in yymmdds]
    yymm_folders = [mission_name + '_' + yymmdd[:4] for yymmdd in yymmdds]
    experiment_names = [utilities.find_experiment_from_yymmdd(yymmdd) for yymmdd in yymmdds]
    for yymm_folder, yymmdd_folder, experiment_name in zip(yymm_folders, yymmdd_folders, experiment_names):
        fdets_folder_path = f'/Users/lgisolfi/Desktop/PRIDE_DATA_NEW/analysed_pride_data/{mission_name}/{yymm_folder}/{yymmdd_folder}/input/complete' #or insert your path
        output_dir = f'/Users/lgisolfi/Desktop/PRIDE_DATA_NEW/analysed_pride_data/{mission_name}/{yymm_folder}/{yymmdd_folder}/output/' #or insert your path

        fdets_folder_path = f'/Users/lgisolfi/Desktop/PRIDE_DATA_NEW/analysed_pride_data/{mission_name}/{yymm_folder}/{yymmdd_folder}/input/complete' #or insert your path
        output_dir = f'/Users/lgisolfi/Desktop/PRIDE_DATA_NEW/analysed_pride_data/{mission_name}/{yymm_folder}/{yymmdd_folder}/output/' #or insert your path

        if os.path.exists(output_dir):
            print(f'Deleting directory {output_dir}')
            shutil.rmtree(output_dir)

        horizons_target = utilities.mission_name_to_horizons_target(mission_name)
        print(f'Performing Statistical Analysis for mission: {mission_name} (Horizons Code: {horizons_target})...')

        # Get list of FDETS files
        dir_path = fdets_folder_path
        files_list = []
        station_ids = []
        for file in os.listdir(dir_path):
            if file.startswith('Fdets') and file.endswith('r2i.txt'):
                files_list.append(os.path.join(dir_path, file))

        # Extract data
        extracted_data_list = process_fdets.extract_folder_data(dir_path)

        #filtered_extracted_data_list = analysis.two_step_filter(extracted_data_list)

        for extracted_data in extracted_data_list:
            station_id = extracted_data['receiving_station_name']
            for file_name in files_list:
                if station_id != process_fdets.get_station_name_from_file(file_name):
                    continue
                else:
                    # Plot user-defined parameters (doppler noise, SNR, fdets)
                    analysis.plot_user_defined_parameters(
                        extracted_data,
                        save_dir = os.path.join(output_dir, 'user_defined_parameters'),
                        plot_snr = True,
                        plot_doppler_noise = True,
                        plot_fdets= True,
                        suppress = True
                    )

                    # Plot elevation for each file (station)
                    analysis.get_elevation_plot(
                        [file_name],
                        horizons_target,
                        [station_id],
                        mission_name =mission_name ,
                        suppress = True,
                        save_dir = os.path.join(output_dir, 'elevation')
                    )

                    station_ids.append(station_id) # Append station_id to list station_ids

        # Plot SNR and Doppler Noise statistics
        analysis.get_all_stations_statistics(
            fdets_folder_path = fdets_folder_path,
            mission_name = mission_name,
            extracted_parameters_list= extracted_data_list,
            doppler_noise_statistics = True,
            snr_statistics= True,
            save_dir = os.path.join(output_dir, 'statistics')
        )

        # Plot combined elevation plot for all stations
        analysis.get_elevation_plot(
            files_list,
            horizons_target,
            station_ids,
            mission_name,
            suppress = True,
            save_dir = os.path.join(output_dir, 'statistics')
        )


        #Combine Images
        snr_noise_folder = os.path.join(output_dir, 'user_defined_parameters/snr_noise_fdets')
        elevation_folder = os.path.join(output_dir, 'elevation')


        # Get lists of image filenames from both folders
        snr_images = sorted([file for file in os.listdir(snr_noise_folder) if file.endswith('.png')])  # Ensure sorted order
        elevation_images = sorted([file for file in os.listdir(elevation_folder) if file.endswith('.png')])  # Ensure sorted order

        # Zip them into tuples (assuming each folder contains matching files)
        images_to_combine = zip(snr_images, elevation_images)

        # Iterate through image pairs and combine them
        for snr_image, elevation_image in images_to_combine:
            station_id = snr_image.split('_')[0]
            date = snr_image.split('_')[1]

            # Create full paths
            snr_image_path = os.path.join(snr_noise_folder, snr_image)
            elevation_image_path = os.path.join(elevation_folder, elevation_image)

            # Call the function with a list of tuples, maintaining correct filename
            utilities.combine_plots(
               image_paths =[snr_image_path, elevation_image_path],  # Pass a list of image paths
               output_dir = os.path.join(output_dir, 'elevation_snr_noise'),
               output_file_name=f'{station_id}_{date}_combined.png',
               direction='vertical'
            )

    print('Done.')

