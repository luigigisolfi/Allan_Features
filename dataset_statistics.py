from pride_characterization_library import PrideDopplerCharacterization
import os

import csv
from collections import defaultdict
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np
import random
import colorsys

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

################################
################################
# DO NOT FORGET TO ADD ed045a!!!
################################
################################

experiments_to_analyze = {'mex': ['gr035'], 'juice': ['ec094a','ec094b'], 'min': ['ed045c','ed045d','ed045e','ed045f']}
mean_rms_user_defined_parameters = defaultdict(list)
color_dict = defaultdict(list)
for mission_name, experiment_names in experiments_to_analyze.items():
    print('new mission')
    for experiment_name in experiment_names:
        color_dict[experiment_name] = generate_random_color()
        fdets_folder_path = f'/Users/lgisolfi/Desktop/data_archiving-1.0/dataset/{mission_name}/{experiment_name}/input/complete/'
        output_dir = f'/Users/lgisolfi/Desktop/data_archiving-1.0/dataset/{mission_name}/{experiment_name}/output'
        user_defined_parameters_dir = os.path.join(output_dir, 'user_defined_parameters/snr_noise_fdets')
        for file in os.listdir(user_defined_parameters_dir):
            station_code = file.split('_')[0]
            if file.endswith('.txt'):
                parameters_dictionary = analysis.read_user_defined_parameters_file(os.path.join(user_defined_parameters_dir, file))
                mean_rms_user_defined_parameters[experiment_name].append({station_code:
                                                 {'mean_snr': np.mean(parameters_dictionary['SNR']),
                                                  'rms_snr': np.std(parameters_dictionary['SNR']),
                                                  'mean_doppler_noise': np.mean(parameters_dictionary['Doppler_noise']),
                                                  'rms_doppler_noise': np.std(parameters_dictionary['Doppler_noise']),
                                                  'mean_fdets': np.mean(parameters_dictionary['Freq_Detection']),
                                                  'rms_fdets': np.std(parameters_dictionary['Freq_Detection'])
                                                  }})


import matplotlib.pyplot as plt

labels_snr = set()
labels_doppler = set()
labels_snr_vs_noise = set()

fig, axes = plt.subplots(3, 1, figsize=(10, 15), sharex=False)  # Three subplots stacked vertically
ax1, ax2, ax3 = axes  # Assign subplots

for experiment_name in mean_rms_user_defined_parameters.keys():
    for station_dict in mean_rms_user_defined_parameters[experiment_name]:
        for station in station_dict.keys():
            mean_snr = [station_dict[station]['mean_snr'] for station in station_dict.keys()]
            rms_snr = [station_dict[station]['rms_snr'] for station in station_dict.keys()]
            mean_doppler = [station_dict[station]['mean_doppler_noise'] for station in station_dict.keys()]
            rms_doppler = [station_dict[station]['rms_doppler_noise'] for station in station_dict.keys()]

            # Plot SNR on the first subplot
            ax1.errorbar(station, mean_snr, yerr=rms_snr, linewidth=2, fmt='o', markersize=5, alpha=0.5,
                         color=color_dict[experiment_name], label=experiment_name if experiment_name not in labels_snr else None)
            labels_snr.add(experiment_name)

            # Plot Doppler Noise on the second subplot
            print(station, experiment_name, mean_doppler)
            ax2.errorbar(station, mean_doppler, yerr=rms_doppler/np.max(rms_doppler), linewidth=2, fmt='s', markersize=5, alpha=0.5,
                         color=color_dict[experiment_name], label=experiment_name if experiment_name not in labels_doppler else None)
            labels_doppler.add(experiment_name)

            # Plot SNR vs. Doppler Noise on the third subplot
            ax3.errorbar(mean_snr, mean_doppler, fmt='o', markersize=5, alpha=0.6,
                         color=color_dict[experiment_name], label=experiment_name if experiment_name not in labels_snr_vs_noise else None)
            labels_snr_vs_noise.add(experiment_name)

# Formatting SNR subplot
ax1.set_ylabel('SNR')
ax1.set_yscale('log')
ax1.grid()
ax1.legend()
ax1.set_title('SNR per Station (No Normalization)')

# Formatting Doppler Noise subplot
ax2.set_ylabel(f'Normalized Doppler Noise')
ax2.grid()
ax2.set_title('Doppler Noise per Station (Normalized)')

# Formatting SNR vs. Doppler Noise subplot
ax3.set_xlabel('SNR')
ax3.set_ylabel('Doppler Noise (Hz)')
ax3.set_xscale('log')
ax3.set_yscale('log')
ax3.grid()
ax3.set_title('SNR vs. Doppler Noise (No Errorbars)')

plt.tight_layout()  # Adjust layout to prevent overlapping
plt.show()
