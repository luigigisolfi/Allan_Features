import os
import subprocess

# Define the path to the parent folder
old_format_files_folder = f"/Users/lgisolfi/Desktop/data_archiving-2.0/vex/usable/converted_old_format_files/"
vex_year_month_list  = [name for name in os.listdir(old_format_files_folder) if os.path.isdir(os.path.join(old_format_files_folder, name))]

for vex_year_month in vex_year_month_list:
    if vex_year_month in ["vex_1401","vex_1402"] :
        continue
    parent_folder = f"/Users/lgisolfi/Desktop/data_archiving-2.0/vex/usable/converted_old_format_files/{vex_year_month}"

    # Get the list of folder names inside vex_1401
    folder_names = [name for name in os.listdir(parent_folder) if os.path.isdir(os.path.join(parent_folder, name))]

    # Iterate over each folder name and call the function
    for name_folder in folder_names:
        input_path = f"{parent_folder}/{name_folder}/input/single/*.txt"
        command = f"python old2new_fdets.py {input_path} --experiment_name {name_folder}"
        subprocess.run(command, shell=True)
