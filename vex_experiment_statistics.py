
"""
This script performs the same as what is performed in experiment_statistics, but it is tailored to the old structure of vex data.
"""
import sys
sys.path.append('/Users/lgisolfi/ClionProjects/Allan_Features/Analysis_Scripts/') # Adjust this to your actual library location
from pride_characterization_library  import PrideDopplerCharacterization
import os

# Initialize classes
pride = PrideDopplerCharacterization()
process_fdets = pride.ProcessFdets()
utilities = pride.Utilities()
analysis = pride.Analysis(process_fdets, utilities)

# Define experiments to analyze
old_format_files_folder = f"/Users/lgisolfi/Desktop/PRIDE_DATA/analysed_pride_data/vex/usable/converted_old_format_files/" # or change path accordingly
vex_year_month_list  = [name for name in os.listdir(old_format_files_folder) if os.path.isdir(os.path.join(old_format_files_folder, name))]
for vex_year_month in vex_year_month_list:

    if vex_year_month != "vex_1411":
        continue
    parent_folder = f"/Users/lgisolfi/Desktop/data_archiving-2.0/vex/usable/converted_old_format_files/{vex_year_month}"
    # Get the list of folder names inside vex_1401
    experiments_list = [name for name in os.listdir(parent_folder) if os.path.isdir(os.path.join(parent_folder, name))]

    experiments_to_analyze = {
        'vex': experiments_list
    }

    # Loop over missions and experiments
    try:
        for mission_name, experiment_names in experiments_to_analyze.items():
            for experiment_name in experiment_names:
                fdets_folder_path = f'/Users/lgisolfi/Desktop/data_archiving-2.0/{mission_name}/usable/converted_old_format_files/{vex_year_month}/{experiment_name}/input/complete'
                output_dir =  f'/Users/lgisolfi/Desktop/data_archiving-2.0/{mission_name}/usable/converted_old_format_files/{vex_year_month}/{experiment_name}/output'
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
                extracted_data_list =  process_fdets.extract_folder_data(dir_path)

                for extracted_data in extracted_data_list:
                    station_id = extracted_data['receiving_station_name']
                    for file_name in files_list:
                        #if str(extracted_data['utc_date']) not in file_name:
                        #    print(extracted_data['utc_date'])
                        #    continue
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
                                experiment_name =experiment_name ,
                                suppress = True,
                                save_dir = os.path.join(output_dir, 'elevation')
                            )

                            station_ids.append(station_id) # Append station_id to list station_ids

                # Plot SNR and Doppler Noise statistics
                analysis.get_all_stations_statistics(
                    fdets_folder_path = fdets_folder_path,
                    experiment_name = experiment_name,
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
                    experiment_name,
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

    except Exception as e:
        continue