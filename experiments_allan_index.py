from pride_characterization_library import PrideDopplerCharacterization
from tudatpy.interface import spice
import os
import glob
spice.load_standard_kernels()
pride = PrideDopplerCharacterization()
process_fdets = pride.ProcessFdets()
utilities = pride.Utilities()
analysis = pride.Analysis(process_fdets, utilities)

experiments_to_analyze = {'mex': ['gr035']}

for mission_name, experiment_names in experiments_to_analyze.items():
    for experiment_name in experiment_names:
        fdets_folder_path = f'/Users/lgisolfi/Desktop/data_archiving-1.0/dataset/{mission_name}/{experiment_name}/input/complete/'
        output_dir = f'/Users/lgisolfi/Desktop/data_archiving-1.0/dataset/{mission_name}/{experiment_name}/output'
        tau_min = 0
        tau_max = 20
        save_dir = output_dir
        suppress = False

        analysis.get_all_stations_madev_plot(fdets_folder_path, experiment_name, tau_min = tau_min, tau_max = tau_max, suppress = False, plot_madev_only = True, save_dir = output_dir)
