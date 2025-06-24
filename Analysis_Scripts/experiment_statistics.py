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
import shutil

# %%
from pride_characterization_library import PrideDopplerCharacterization
import os

# %%
# Initialize classes
pride = PrideDopplerCharacterization()
process_fdets = pride.ProcessFdets()
utilities = pride.Utilities()
analysis = pride.Analysis(process_fdets, utilities)

# %%
# Define experiments to analyze
experiments_to_analyze = {
    'juice': ["ec094b", "ec094a"],
    'mex': ['gr035'],
    'min': ['ed045a','ed045c','ed045d','ed045e','ed045f'],
    'mro': ['ec064'],
    'vex': ['vex_140106','vex_140109','vex_140110','vex_140113','vex_140118','vex_140119','vex_140120','vex_140123','vex_140126','vex_140127', 'vex_140131']
}

# Division by experiments
for mission_name, experiment_names in experiments_to_analyze.items():
    for experiment_name in experiment_names:

        if mission_name == 'vex':
            fdets_folder_path = f'/Users/lgisolfi/Desktop/PRIDE_DATA/analysed_pride_data/{mission_name}/usable/converted_old_format_files/vex_1401/{experiment_name}/input/complete' #or insert your path
            output_dir = f'/Users/lgisolfi/Desktop/PRIDE_DATA/analysed_pride_data/{mission_name}/usable/converted_old_format_files/vex_1401/{experiment_name}/output/' #or insert your path

        else:
            fdets_folder_path = f'/Users/lgisolfi/Desktop/PRIDE_DATA/analysed_pride_data/{mission_name}/{experiment_name}/input/complete' #or insert your path
            output_dir = f'/Users/lgisolfi/Desktop/PRIDE_DATA/analysed_pride_data/{mission_name}/{experiment_name}/output/' #or insert your path

        if os.path.exists(output_dir):
            print(f'Deleting directory {output_dir}')
            shutil.rmtree(output_dir)

# Uncomment for division by date.
# Loop over missions and experiments by date
#experiments_to_analyze = {'mex': ['mex_131228', 'mex_131229']}
#for mission_name, dates in experiments_to_analyze.items():
#    for date in dates:
#        corresponding_month_folder = date[:-2]
#        fdets_folder_path = f'/Users/lgisolfi/Desktop/PRIDE_DATA_NEW/analysed_pride_data/{mission_name}/{corresponding_month_folder}/{date}/input/complete' #or insert your path
#        output_dir = f'/Users/lgisolfi/Desktop/PRIDE_DATA_NEW/analysed_pride_data/{mission_name}/{corresponding_month_folder}/{date}/output_test/' #or insert your path
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

# %%
