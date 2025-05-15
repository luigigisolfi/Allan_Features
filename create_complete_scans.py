from tudatpy.interface import spice
import os
import sys
sys.path.append('/Users/lgisolfi/ClionProjects/Allan_Features/Analysis_Scripts/') # Adjust this to your actual library location
from pride_characterization_library import PrideDopplerCharacterization
import glob

spice.load_standard_kernels()
pride = PrideDopplerCharacterization()
process_fdets = pride.ProcessFdets()
utilities = pride.Utilities()
analysis = pride.Analysis(process_fdets, utilities)

# Define the path to the parent folder
old_format_files_folder = f"/Users/lgisolfi/Desktop/data_archiving-2.0/vex/usable/converted_old_format_files/"
vex_year_month_list  = [name for name in os.listdir(old_format_files_folder) if os.path.isdir(os.path.join(old_format_files_folder, name))]

for vex_year_month in vex_year_month_list:
    if vex_year_month not in ["vex_1406"] :
        continue
    parent_folder = f"/Users/lgisolfi/Desktop/data_archiving-2.0/vex/usable/converted_old_format_files/{vex_year_month}"

    # Get the list of folder names inside vex_1401
    folder_names = [name for name in os.listdir(parent_folder) if os.path.isdir(os.path.join(parent_folder, name))]

    for folder_name in folder_names:
        print(folder_name)
        files = glob.glob(f'/Users/lgisolfi/Desktop/data_archiving-2.0/vex/usable/converted_old_format_files/{vex_year_month}/{folder_name}/input/single/Fdets.*.r2i.txt')
        output_folder = f'/Users/lgisolfi/Desktop/data_archiving-2.0/vex/usable/converted_old_format_files/{vex_year_month}/{folder_name}/input/complete'
        utilities.create_complete_scan_from_single_scans(files, output_folder)
