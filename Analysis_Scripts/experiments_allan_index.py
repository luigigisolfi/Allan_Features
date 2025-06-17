# %%
"""
This script allows for plotting of the Modified Allan Deviation per mission_name, given a date_folder, for all stations.
"""
from pride_characterization_library import PrideDopplerCharacterization
pride = PrideDopplerCharacterization() # Create PrideDopplerCharacterization Object
process_fdets = pride.ProcessFdets() # Create ProcessFdets Object
utilities = pride.Utilities() # Create Utilities Object
analysis = pride.Analysis(process_fdets, utilities)  # Create Analysis Object

# %%
# Experiments to analyze
experiments_to_analyze = {'juice': ['juice_231019']}

# %%
for mission_name, dates in experiments_to_analyze.items():
    for date in dates:
        corresponding_month_folder = date[:-2]
        fdets_folder_path = f'/Users/lgisolfi/Desktop/PRIDE_DATA_NEW/analysed_pride_data/{mission_name}/{corresponding_month_folder}/{date}/input/complete' #or insert your path
        output_dir = f'/Users/lgisolfi/Desktop/PRIDE_DATA_NEW/analysed_pride_data/{mission_name}/{corresponding_month_folder}/{date}/output/' #or insert your path
        tau_min = 0
        tau_max = 100
        save_dir = output_dir
        suppress = False

        analysis.get_all_stations_madev_plot(fdets_folder_path, mission_name, tau_min = tau_min, tau_max = tau_max, suppress = False, plot_madev_only = True, save_dir = output_dir)

