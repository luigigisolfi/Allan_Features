#This script allows for plotting of the Modified Allan Deviation per mission_name, given a date_folder, for all stations.
from pride_characterization_library import PrideDopplerCharacterization
pride = PrideDopplerCharacterization() # Create PrideDopplerCharacterization Object
process_fdets = pride.ProcessFdets() # Create ProcessFdets Object
utilities = pride.Utilities() # Create Utilities Object
analysis = pride.Analysis(process_fdets, utilities)  # Create Analysis Object

# %%
# Experiments to analyze
experiments_to_analyze = {
    'juice': ["ec094a", "ec094b"],
    'mex': ['gr035'],
    'min': ['ed045a','ed045c','ed045d','ed045e','ed045f'],
    'mro': ['ec064'],
    'vex': ['vex_140106','vex_140109','vex_140110','vex_140113','vex_140118','vex_140119','vex_140120','vex_140123','vex_140126','vex_140127', 'vex_140131']
}

# %%
for mission_name, experiment_names in experiments_to_analyze.items():
        for experiment_name in experiment_names:
            fdets_folder_path = f'/Users/lgisolfi/Desktop/PRIDE_DATA/analysed_pride_data/{mission_name}/{experiment_name}/input/complete' #or insert your path

            if mission_name == 'vex':
                output_dir = f'/Users/lgisolfi/Desktop/PRIDE_DATA/analysed_pride_data/{mission_name}/usable/converted_old_format_files/vex_1401/{experiment_name}/output/' #or insert your path
            else:
                output_dir = f'/Users/lgisolfi/Desktop/PRIDE_DATA/analysed_pride_data/{mission_name}/{experiment_name}/output/' #or insert your path

            tau_min = 0
            tau_max = 100
            save_dir = output_dir
            suppress = False

            analysis.get_all_stations_oadev_plot(fdets_folder_path, mission_name, experiment_name, tau_min = tau_min, tau_max = tau_max, z_score_filter = True,save_dir = output_dir)


#experiments_to_analyze = {'juice': ['juice_231019']}
#for mission_name, dates in experiments_to_analyze.items():
#    for date in dates:
#        corresponding_month_folder = date[:-2]
#        fdets_folder_path = f'/Users/lgisolfi/Desktop/PRIDE_DATA_NEW/analysed_pride_data/{mission_name}/{corresponding_month_folder}/{date}/input/complete' #or insert your path
#        output_dir = f'/Users/lgisolfi/Desktop/PRIDE_DATA_NEW/analysed_pride_data/{mission_name}/{corresponding_month_folder}/{date}/output_test/' #or insert your path
#        tau_min = 0
#        tau_max = 100
#        save_dir = output_dir
#        suppress = False
#
#        analysis.get_all_stations_oadev_plot(fdets_folder_path, mission_name, tau_min = tau_min, tau_max = tau_max, suppress = False, plot_oadev_only = True, save_dir = output_dir)
