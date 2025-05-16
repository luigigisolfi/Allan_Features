from pride_characterization_library import PrideDopplerCharacterization
import os
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
import random
from scipy.interpolate import interp1d

"""
This script was created to check the dependence of SNR on elevation, 
but the fact that Horizons does not allow for ephemeris to be retrieved at steps < 1 min, together with the fact that - for some stations -
the elevation plot has a maximum value, limits its applications.
Nevertheless, it might still be useful in some occasions, so we keep it. 
"""

def generate_random_color():
    """Generates a random, well-spaced color in hexadecimal format."""
    r = random.randint(0, 220)  # Avoid extremes (too dark/light)
    g = random.randint(0, 220)
    b = random.randint(0, 220)
    return "#{:02x}{:02x}{:02x}".format(r, g, b)

pride = PrideDopplerCharacterization()
process_fdets = pride.ProcessFdets()
utilities = pride.Utilities()
analysis = pride.Analysis(process_fdets, utilities)

experiments_to_analyze = {'juice': ['ec094b']}
mean_rms_user_defined_parameters = defaultdict(list)
mean_elevations = defaultdict(list)
color_dict = defaultdict(list)
for mission_name, experiment_names in experiments_to_analyze.items():
    for experiment_name in experiment_names:
        color_dict[experiment_name] = generate_random_color()
        fdets_folder_path = f'../small_dataset/{mission_name}/{experiment_name}/input/complete' #or insert your path
        output_dir =  f'../small_dataset/{mission_name}/{experiment_name}/output/' #or insert your path
        if not os.path.exists(output_dir):
            print(f'The folder {output_dir} does not exist. Skipping...')
            continue
        user_defined_parameters_dir = os.path.join(output_dir, 'user_defined_parameters/snr_noise_fdets')
        elevation_dir = os.path.join(output_dir, 'elevation')
        for file in os.listdir(user_defined_parameters_dir):
            station_code = file.split('_')[0]
            if station_code == 'Ib':
                station_code = 'Ir'
            elif station_code == 'O6':
                station_code = 'On'
            elif mission_name == 'mex' and (station_code == 'Mh' or station_code == 'Ht'): # these are filtered out manually cause we have good reasons
                continue
            elif mission_name == 'min' and (station_code == 'Ir'): # these are filtered out manually cause we have good reasons
                continue
            if file.endswith('.txt'):
                parameters_dictionary = analysis.read_user_defined_parameters_file(os.path.join(user_defined_parameters_dir, file))

                snr_array = np.array(parameters_dictionary['SNR'])
                doppler_noise_array = np.array(parameters_dictionary['Doppler_noise'])

                # apply no SNR filters for now (as each mission would require different filters)
                # and we want to keep it as general as possible
                # can always filter later if needed
                # for instance:
                # min and mro - SNR > 30
                # mex and juice - SNR > 100 but some stations require lower (e.g. Ir [EC094A] requires SNR > 50)
                if mission_name in ['juice', 'mex', 'vex']:
                    filter_snr = snr_array > 100
                    if experiment_name == 'ec094a' and station_code in ['Wz']:
                        filter_doppler_noise = np.abs(doppler_noise_array) < 0.1
                        filtered_snr = snr_array[filter_snr & filter_doppler_noise]
                        filtered_doppler_noise = doppler_noise_array[filter_snr & filter_doppler_noise]

                    #elif experiment_name == 'ec094b' and station_code in ['Mh', 'Nt']:
                        #filter_doppler_noise = np.abs(doppler_noise_array) < 0.1
                        #filtered_snr = snr_array[filter_snr & filter_doppler_noise]
                        #filtered_doppler_noise = doppler_noise_array[filter_snr & filter_doppler_noise]

                    else:
                        filtered_snr = snr_array[filter_snr]
                        filtered_doppler_noise = doppler_noise_array[filter_snr]

                elif mission_name in ['mro', 'min']:
                    filter_snr = snr_array > 30
                    if experiment_name == 'ed045a' and station_code in ['Bd', 'T6', 'Hh']:
                        filter_doppler_noise = np.abs(doppler_noise_array) < 0.1
                        filtered_snr = snr_array[filter_snr & filter_doppler_noise]
                        filtered_doppler_noise = doppler_noise_array[filter_snr & filter_doppler_noise]

                    elif experiment_name == 'ed045e' and station_code in ['Wb']:
                        filter_doppler_noise = np.abs(doppler_noise_array) < 0.1
                        filtered_snr = snr_array[filter_snr & filter_doppler_noise]
                        filtered_doppler_noise = doppler_noise_array[filter_snr & filter_doppler_noise]

                    elif experiment_name == 'ed045f' and station_code in ['Wz']:
                        filter_doppler_noise = np.abs(doppler_noise_array) < 0.1
                        filtered_snr = snr_array[filter_snr & filter_doppler_noise]
                        filtered_doppler_noise = doppler_noise_array[filter_snr & filter_doppler_noise]
                    else:
                        filtered_snr = snr_array[filter_snr]
                        filtered_doppler_noise = doppler_noise_array[filter_snr]

                mean_rms_user_defined_parameters[experiment_name].append({station_code:
                                                                              {'snr': snr_array}})

        for file in os.listdir(elevation_dir):
            station_code = file.split('_')[0]
            if station_code == 'Ib':
                station_code = 'Ir'
            elif station_code == 'O6':
                station_code = 'On'
            if file.endswith('.txt'):
                elevations = analysis.get_elevations_from_file(os.path.join(elevation_dir, file))
                # Find and update the correct dictionary entry for this station
                for entry in mean_rms_user_defined_parameters[experiment_name]:
                    if station_code in entry:  # If station exists, update it
                        entry[station_code]['elevation'] = elevations
                        break  # No need to continue searching


                for station, values in entry.items():
                    elevation = np.array(values['elevation'])
                    snr = np.array(values['snr'])

                    # Sort values by elevation (necessary for interpolation)
                    sorted_indices = np.argsort(elevation)
                    elevation = elevation[sorted_indices]
                    snr = snr[sorted_indices]

                    # Create interpolation function (linear or cubic)
                    interp_func = interp1d(elevation, snr, kind='linear', fill_value="extrapolate")

                    plt.plot(elevation, snr, 'o-', label=f'{station} (original)')

                    plt.xlabel("Elevation (deg)")
                    plt.ylabel("SNR (dB)")
                    plt.title("Elevation vs. SNR with Interpolation")
                    plt.legend()
                    plt.grid(True)
                    plt.show()