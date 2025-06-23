# %%
# %matplotlib inline

# %% [markdown]
"""
This script analyzes radio science experiments for JUICE
computing and plotting Key Performance Indicators (foms) such as mean Signal-to-Noise Ratio (SNR),
rms SNR, mean Doppler noise, and mean elevation angle across multiple ground stations.

The workflow includes:
- Loading user-defined parameters (SNR, Doppler noise) and elevation data
- Filtering based on SNR and Doppler noise thresholds
- Computing weighted and unweighted averages
- Visualizing results in a set of subplots

Requirements:
    - The data files where the foms are read from are must be present in the output_dir folder, and they are created via experiment_statistics.py .
"""
# %%
from pride_characterization_library import PrideDopplerCharacterization
import os
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
import random
import re
from scipy.stats import norm

# %%
def generate_random_color():
    """Generates a random, well-spaced color in hexadecimal format."""
    r = random.randint(0, 220)  # Avoid extremes (too dark/light)
    g = random.randint(0, 220)
    b = random.randint(0, 220)
    return "#{:02x}{:02x}{:02x}".format(r, g, b)
# %%
pride = PrideDopplerCharacterization() # create PRIDE Object
process_fdets = pride.ProcessFdets() # create Process Fdets Object
utilities = pride.Utilities() # create Utilities Object
analysis = pride.Analysis(process_fdets, utilities) # create Analysis Object

# Select the experiment(s) for which data analysis will be performed
missions_to_analyse = ['juice', 'mex', 'min', 'mro', 'vex']

allowed_mean_doppler_filter = 0.01 #Hz = 100 mHz
bad_observations_mean_doppler_filter = 0.005 #Hz = 5 mHz
bad_obs_flag = False # if set to true, it 1) plots the observations as flagged and 2) removes them from the final statistics for mean FoM computation

# Create empty dictionaries to be filled with meaningful values
mean_rms_user_defined_parameters = defaultdict(list)
mean_elevations = defaultdict(list)
color_dict = defaultdict(list)

# Loop through missions and experiments
for mission_name in missions_to_analyse:
    count = 0
    if mission_name == 'vex':
        color_dict[mission_name] = 'red'
    else:
        color_dict[mission_name] = generate_random_color()

    if mission_name == 'mro':
        bad_observations_mean_doppler_filter = 0.03 #Hz = 30 mHz

    root_dir = f'/Users/lgisolfi/Desktop/PRIDE_DATA_NEW/analysed_pride_data/'
    mission_root = os.path.join(root_dir, mission_name)

    if not os.path.exists(mission_root):
        print(f"Mission root path '{mission_root}' does not exist.")
        exit()

    for yymm_folder in sorted(os.listdir(mission_root)):
        if mission_name == 'vex' and yymm_folder != 'vex_1401':
            continue
        yymm_path = os.path.join(mission_root, yymm_folder)


        if not os.path.isdir(yymm_path) or not yymm_folder.startswith(mission_name + '_'):
            continue

        for yymmdd_folder in sorted(os.listdir(yymm_path)):

            experiment_name = yymmdd_folder
            full_folder_path = os.path.join(yymm_path, yymmdd_folder)

            if not os.path.isdir(full_folder_path) or not yymmdd_folder.startswith(mission_name):
                continue
            # Define subfolder paths
            fdets_folder_path = os.path.join(full_folder_path, 'input', 'complete')
            output_dir = os.path.join(full_folder_path, 'output')
            if not os.path.exists(output_dir):
                print(f'The folder {output_dir} does not exist. Skipping...')
                continue

            # Assign color based on the month (yymm), not the experiment date (yymmdd)
            if yymm_folder not in color_dict:
                yymm = re.sub(r'^[^_]+_', '', yymm_folder)
                color_dict[yymm] = generate_random_color()

            current_color = color_dict[yymm]
            print(f'✅ Processing: {fdets_folder_path} | Color: {current_color}')

            # Paths for user-defined parameters and elevation data
            user_defined_parameters_dir = os.path.join(output_dir, 'user_defined_parameters/snr_noise_fdets')
            elevation_dir = os.path.join(output_dir, 'elevation')

            # Analyze user-defined parameter files
            for file in os.listdir(user_defined_parameters_dir):
                station_code = file.split('_')[0]

                # Some station codes in the fdets name are wrong.
                if station_code == 'Ib':
                    station_code = 'Ir'
                elif station_code == 'O6':
                    station_code = 'On'

                # Process only TXT files
                if file.endswith('.txt'):
                    parameters_dictionary = analysis.read_user_defined_parameters_file(os.path.join(user_defined_parameters_dir, file))

                    dates_and_snr = np.array(parameters_dictionary['SNR'])
                    dates_and_doppler_noise = np.array(parameters_dictionary['Doppler_noise'])
                    dates, doppler_noise_list = zip(*dates_and_doppler_noise)
                    _, snr_list = zip(*dates_and_snr)

                    snr_array = np.array(snr_list)
                    doppler_noise_array = np.array(doppler_noise_list)

                    median_doppler_noise = np.median(doppler_noise_array)
                    mean_doppler_noise = np.mean(doppler_noise_array)

                    median = np.median(doppler_noise_array)
                    mad = np.median(np.abs(doppler_noise_array - median))
                    modified_z = 0.6745 * (doppler_noise_array - median) / mad

                    threshold = 3.5  # Standard threshold for outlier detection
                    z_score_filter = np.abs(modified_z) < threshold

                    average_filter_doppler_noise = np.abs(doppler_noise_array) < allowed_mean_doppler_filter

                    filter_indexes = np.where(z_score_filter)[0]
                    average_filter_indexes = np.where(average_filter_doppler_noise)[0]

                    filtered_dates = np.array(dates)[filter_indexes]
                    averaged_filtered_dates = np.array(dates)[average_filter_indexes]
                    # Fit Gaussian
                    mu, std = norm.fit(doppler_noise_array)

                    median_retained = len(filtered_dates)*100/len(dates)
                    average_retained = len(averaged_filtered_dates)*100/len(dates)

                    print(f'{station_code},Median Retained %: {median_retained}')
                    print(f'{station_code},Average Retained %: {average_retained}')

                    # Plot Gaussian
                    #fig, ax = plt.subplots(1, 1, figsize=(10, 5), sharex=True)
                    #ax.hist(doppler_noise_array, bins=30, density=True, alpha=0.6)
                    #xmin, xmax = ax.get_xlim()  # Get the range from the histogram plot
                    #x = np.linspace(xmin, xmax, 100)  # Create 100 evenly spaced points between
                    #p = norm.pdf(x, mu, std)
                    #ax.plot(x, p, 'k', linewidth=2, label=f'Gaussian fit: μ={mu:.3f}, σ={std:.3f}')
                    #ax.legend()
                    #plt.show()

                    # Apply median filtering
                    filtered_snr = snr_array[z_score_filter]
                    filtered_doppler_noise = doppler_noise_array[z_score_filter]
                    average_filtered_doppler_noise = doppler_noise_array[average_filter_doppler_noise]


                    # Plot difference between mean filter and z-score filter
                    #plt.plot(dates, doppler_noise_array, label='Original', marker='+', linestyle='-',
                    #            color='blue', markersize=5, linewidth=0.5)
                    #plt.plot(averaged_filtered_dates, average_filtered_doppler_noise, label=f'Window-Filtered, Retained: {average_retained}', marker='+', linestyle='-',
                    #         color='red', markersize =5, linewidth=0.5)
                    #plt.plot(filtered_dates, filtered_doppler_noise, label=f'Modified Z-Score, Retained: {median_retained}', marker='+', linestyle='-',
                    #            color='orange', markersize =5, linewidth=0.5)

                    #plt.title(f'Doppler Noise for {mission_name}, Station: {station_code}')
                    #plt.legend()
                    #plt.show()
                    # Save mean and RMS values
                    mean_rms_user_defined_parameters[experiment_name].append({station_code:
                                                                                  {'mean_snr': np.mean(filtered_snr),
                                                                                   'rms_snr': np.std(filtered_snr),
                                                                                   'mean_doppler_noise': np.mean(filtered_doppler_noise),
                                                                                   'rms_doppler_noise': np.std(filtered_doppler_noise),
                                                                                   }})

            # Analyze elevation files
            for file in os.listdir(elevation_dir):
                station_code = file.split('_')[0]
                if station_code == 'Ib':
                    station_code = 'Ir'
                elif station_code == 'O6':
                    station_code = 'On'
                if file.endswith('.txt'):
                    mean_elevation = analysis.get_mean_elevation_from_file(os.path.join(elevation_dir, file))


                    # Update mean elevation in the main dictionary
                    for entry in mean_rms_user_defined_parameters[experiment_name]:
                        if station_code in entry:
                            entry[station_code]['mean_elevation'] = mean_elevation

# Initialize labels and plot
labels_snr = set()
labels_doppler = set()
labels_snr_vs_noise = set()
fig, axes = plt.subplots(4, 1, figsize=(10, 15), sharex=False)
ax1, ax2, ax3, ax4 = axes  # Assign subplots

# Plotting results
count = 0
count_bad = 0
for experiment_name in mean_rms_user_defined_parameters.keys():
    mission_name = experiment_name.split('_')[0]
    for station_dict in mean_rms_user_defined_parameters[experiment_name]:
        for station in station_dict.keys():
            try:
                antenna_diameter = utilities.antenna_diameters[station]
            except:
                continue
            mean_snr = 10 * np.log10(station_dict[station]['mean_snr'])  # dB
            rms_snr = 10 * np.log10(station_dict[station]['rms_snr'])    # dB
            mean_doppler = station_dict[station]['mean_doppler_noise'] * 1000  # mHz
            rms_doppler = station_dict[station]['rms_doppler_noise'] * 1000    # mHz
            mean_elevation = station_dict[station]['mean_elevation']
            # === Prepare label and color ===
            yymmdd = re.sub(r'^[^_]+_', '', experiment_name)
            yymm = yymmdd[:4]
            year = yymm[:2]
            month = yymm[2:]
            label_str = f'{mission_name.upper()} 20{year}/{month}'
            color = color_dict[yymm]
            # Determine label for legend (once per plot)
            label = label_str if label_str not in labels_snr else None
            count+=1
            # Only plot valid observations
            if np.abs(mean_doppler) > bad_observations_mean_doppler_filter*1000: # in mHz:
                if bad_obs_flag:
                    count_bad +=1
                    print(f'Bad Observation: {station}, {mission_name}, {yymmdd}')
                    ax1.errorbar(station, mean_snr, fmt='x', markersize=6, alpha=0.7, color='red')
                    ax1.errorbar(station, mean_snr, fmt='o', linewidth=2, markersize=6, alpha=0.5,
                                 color=color, label=label)
                    ax2.errorbar(station, mean_doppler, fmt='x', markersize=7, alpha=0.7, color='red')
                    ax2.errorbar(station, mean_doppler, fmt='o', linewidth=2, markersize=6, alpha=0.5,
                                 color=color, label=label)
                    ax3.errorbar(mean_snr, rms_doppler, fmt='x', markersize=6, alpha=0.7, color='red')
                    ax3.errorbar(mean_snr, rms_doppler, fmt='o', markersize=6, alpha=0.6,
                                 color=color, label=label)
                    ax3.annotate(station, (mean_snr, rms_doppler), fontsize=7, alpha=0.7)
                    ax4.errorbar(mean_elevation, mean_snr, fmt='x', markersize=3 * antenna_diameter / 10, alpha=0.7, color='red')
                    ax4.errorbar(mean_elevation, mean_snr, fmt='o', markersize=3 * antenna_diameter / 10, alpha=0.6,
                                 color=color, label=label)
                    ax4.annotate(station, (mean_elevation, mean_snr), fontsize=6, alpha=0.7)
                continue

            # Add label string to memory (after using it)
            labels_snr.add(label_str)
            labels_snr_vs_noise.add(label_str)

            # === Plotting ===
            ax1.errorbar(station, mean_snr, fmt='o', linewidth=2, markersize=6, alpha=0.5,
                         color=color, label=label)
            ax2.errorbar(station, mean_doppler, fmt='o', linewidth=2, markersize=6, alpha=0.5,
                         color=color, label=label)
            ax3.errorbar(mean_snr, rms_doppler, fmt='o', markersize=6, alpha=0.6,
                         color=color, label=label)
            ax3.annotate(station, (mean_snr, rms_doppler), fontsize=7, alpha=0.7)
            ax4.errorbar(mean_elevation, mean_snr, fmt='o', markersize=3 * antenna_diameter / 10, alpha=0.6,
                         color=color, label=label)
            ax4.annotate(station, (mean_elevation, mean_snr), fontsize=7, alpha=0.7)

bad_observations_percentage = (count_bad/count)*100
print(f"Percentage of Flagged Scans: {bad_observations_percentage} %")

# Final plot adjustments
ax1.set_ylabel('SNR [dB]')
ax1.set_xlabel(f'Station Code')
ax1.grid()
ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

# Formatting Doppler Noise subplot
ax2.set_ylabel(f'Doppler Noise [mHz]')
ax2.set_xlabel(f'Station Code')
ax2.grid()

# Formatting SNR vs. Doppler Noise subplot
ax3.set_xlabel('SNR [dB]')
ax3.set_ylabel('RMS Doppler Noise [mHz]')
ax3.grid()

# Formatting SNR vs. Doppler Noise subplot
ax4.set_xlabel('Elevation [deg]')
ax4.set_ylabel('SNR [dB]')
ax4.grid()

plt.savefig(os.path.join(output_dir, experiment_name) + '_dataset_statistics.png')

plt.tight_layout(pad=2)
plt.show()

### This part of the code computes the average values of the Figures of Merit (FoMs)
### among all stations belonging to a given experiment,
### and it does so by weighting the noise by its SNR.

mission_aggregates = defaultdict(lambda: {
    'mean_snr': [],
    'mean_doppler_noise': [],
    'rms_doppler_noise': []
})

for experiment_name in mean_rms_user_defined_parameters.keys():

    mission_name = experiment_name.split('_')[0]

    # First step: Unweighted means
    aggregate = defaultdict(lambda: {'sum': 0, 'count': 0})
    for entry in mean_rms_user_defined_parameters[experiment_name]:
        for station, values in entry.items():
            for key, value in values.items():
                if key not in ['mean_doppler_noise', 'rms_doppler_noise']:
                    aggregate[key]['sum'] += value
                    aggregate[key]['count'] += 1

    mean_values_old = {key: aggregate[key]['sum'] / aggregate[key]['count'] for key in aggregate}

    # Second step: Weighted means using SNR
    aggregate_weighted = defaultdict(lambda: {'sum': 0, 'count': 0})
    for entry in mean_rms_user_defined_parameters[experiment_name]:
        for station, values in entry.items():
            snr_weight = values.get('mean_snr', 0)
            doppler = values.get('mean_doppler_noise', 0)
            if np.abs(doppler) > bad_observations_mean_doppler_filter:
                if bad_obs_flag and mission_name != 'mro': # exclude mro because of one-way doppler (higher RMS expected)
                    print(f'Skipping station {station} for mission {mission_name}')
                    continue
            for key, value in values.items():
                if key in ['mean_doppler_noise', 'rms_doppler_noise']:
                    aggregate_weighted[key]['sum'] += value * snr_weight
                    aggregate_weighted[key]['count'] += snr_weight

    mean_values_new = {
        key: (aggregate_weighted[key]['sum'] / aggregate_weighted[key]['count']) * 1000
        for key in aggregate_weighted
    }

    mission_aggregates[mission_name]['mean_snr'].append(mean_values_old['mean_snr'])
    for key in ['mean_doppler_noise', 'rms_doppler_noise']:
        if key in mean_values_new:
            mission_aggregates[mission_name][key].append(mean_values_new[key])


final_means_per_mission = {
    mission: {
        key: np.mean(values) if values else None
        for key, values in foms.items()
    }
    for mission, foms in mission_aggregates.items()
}

# Optional: Print results
for mission, foms in final_means_per_mission.items():
    print(f"\nMission: {mission}")
    print(f"Mean SNR (dB): {10 * np.log10(foms['mean_snr']) if foms['mean_snr'] else 'N/A'}")
    print(f"Mean Doppler Noise (mHz): {foms['mean_doppler_noise']}")
    print(f"RMS Doppler Noise (mHz): {foms['rms_doppler_noise']}")