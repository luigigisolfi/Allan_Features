# %%
"""
This script allows for plotting of the Modified Allan Deviation per experiment, for all stations
(equivalent to the MAD plots in Vidhya's paper)
"""
from pride_characterization_library import PrideDopplerCharacterization
pride = PrideDopplerCharacterization() # Create PrideDopplerCharacterization Object
process_fdets = pride.ProcessFdets() # Create ProcessFdets Object
utilities = pride.Utilities() # Create Utilities Object
analysis = pride.Analysis(process_fdets, utilities)  # Create Analysis Object

# %%
# Experiments to analyze
experiments_to_analyze = {'mex': ['gr035']}

# %%
for mission_name, experiment_names in experiments_to_analyze.items():
    for experiment_name in experiment_names:
        fdets_folder_path = '/Users/lgisolfi/Desktop/data_archiving-1.0/dataset/mex/gr035/input/complete'
        output_dir =  '/Users/lgisolfi/Desktop/data_archiving-1.0/dataset/mex/gr035/output'
        tau_min = 0
        tau_max = 100
        save_dir = output_dir
        suppress = False

        analysis.get_all_stations_madev_plot(fdets_folder_path, experiment_name, tau_min = tau_min, tau_max = tau_max, suppress = False, plot_madev_only = True, save_dir = output_dir)

# %%
