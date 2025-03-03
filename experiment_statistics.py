from pride_characterization_library import PrideDopplerCharacterization
from tudatpy.interface import spice
import os
import glob
spice.load_standard_kernels()
pride = PrideDopplerCharacterization()
process_fdets = pride.ProcessFdets()
utilities = pride.Utilities()
analysis = pride.Analysis(process_fdets, utilities)

#time_interval_minutes= 120
#for fdets_file in os.listdir(f'/Users/lgisolfi/Desktop/data_archiving-1.0/dataset/mex/gr035/complete'):
#    utilities.split_scan_by_time(
#       input_folder=f'/Users/lgisolfi/Desktop/data_archiving-1.0/dataset/mex/gr035/complete',
#        fdets_file = fdets_file,
#        time_interval_minutes= time_interval_minutes,
#        output_folder= f'/Users/lgisolfi/Desktop/data_archiving-1.0/dataset/mex/gr035/complete_{time_interval_minutes}_minutes'
#    )
#exit()

experiments_to_analyze = {'mex': ['gr035']}

for mission_name, experiment_names in experiments_to_analyze.items():
    for experiment_name in experiment_names:
        fdets_folder_path = f'/Users/lgisolfi/Desktop/data_archiving-1.0/dataset/{mission_name}/{experiment_name}/input/complete_badary/'
        output_dir = f'/Users/lgisolfi/Desktop/data_archiving-1.0/dataset/{mission_name}/{experiment_name}/input/complete_badary/'
        horizons_target = utilities.mission_name_to_horizons_target(mission_name)
        print(f'Performing Statistical Analysis for mission: {mission_name} (Horizons Code: {horizons_target})...')
        dir_path = fdets_folder_path
        files_list = []
        station_ids = []
        for file in os.listdir(dir_path):
            if file.startswith('Fdets') and file.endswith('.txt'):
                files_list.append(os.path.join(dir_path, file))

        extracted_data_list =  process_fdets.extract_folder_data(dir_path)

        for extracted_data in extracted_data_list:
            station_id = extracted_data['receiving_station_name']
            for file_name in files_list:
                if str(extracted_data['utc_date']) not in file_name:
                    continue
                if station_id != process_fdets.get_station_name_from_file(file_name):
                    continue
                else:
                    # Plot noise, elevation and combine them
                    analysis.plot_user_defined_parameters(
                        extracted_data,
                        save_dir = os.path.join(output_dir, 'user_defined_parameters'),
                        plot_snr = True,
                        plot_doppler_noise = True,
                        suppress = True
                    )
                    analysis.get_elevation_plot(
                        [file_name],
                        horizons_target,
                        [station_id],
                        experiment_name =experiment_name ,
                        suppress = True,
                        save_dir = os.path.join(output_dir, 'elevation')
                    )

                    station_ids.append(station_id) # Append station_id to list station_ids

        # Plot SNR and Doppler Noise statistics with errorbars
        analysis.get_all_stations_statistics(
            fdets_folder_path = fdets_folder_path,
            experiment_name = experiment_name,
            extracted_parameters_list= extracted_data_list,
            doppler_noise_statistics = True,
            snr_statistics= True,
            save_dir = os.path.join(output_dir, 'statistics')
        )

        # Plot elevation for all available stations
        analysis.get_elevation_plot(
            files_list,
            horizons_target,
            station_ids,
            experiment_name,
            suppress = True,
            save_dir = os.path.join(output_dir, 'statistics')
        )

        #Combine Images
        snr_noise_folder = os.path.join(output_dir, 'user_defined_parameters/snr_noise')
        elevation_folder = os.path.join(output_dir, 'elevation')


        # Get lists of image filenames from both folders
        snr_images = sorted(os.listdir(snr_noise_folder))  # Ensure sorted order
        elevation_images = sorted(os.listdir(elevation_folder))  # Ensure sorted order

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