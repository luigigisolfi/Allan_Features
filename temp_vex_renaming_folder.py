import os
import re

"""
This code was used on the old vex data folder structure (old structure can be found in PRIDE DATA folder). 
The script temp_vex_folder_structure is also applied to reach the final version of the vex folder structure. 
The new final structure for vex folders, which makes vex structure homogenous with all other missions, can be found in PRIDE DATA NEW.
"""
def rename_vex_folders(root_directory: str) -> None:
    """
    Rename folders inside the given root directory that match the pattern 'vYYMMDD'
    to the format 'vex_YYMMDD'.

    Parameters
    ----------
    root_directory : str
        The root directory path where the folders to rename are located.

    Behavior
    --------
    - Recursively walks through all subdirectories inside root_directory, starting from the bottom.
    - Identifies folders whose names start with 'v' followed by exactly 6 digits (YYMMDD).
    - Renames each matching folder by prefixing it with 'vex_' instead of just 'v'.
    - Prints the old and new folder paths during renaming.

    Notes
    -----
    - Uses a regex pattern to strictly match folder names with the specified format.
    - The renaming occurs in-place in the same directory as the original folder.
    """
    pattern = re.compile(r'^v(\d{6})$')  # matches 'v' followed by exactly 6 digits

    for dirpath, dirnames, filenames in os.walk(root_directory, topdown=False):
        for dirname in dirnames:
            match = pattern.match(dirname)
            if match:
                yymmdd = match.group(1)
                new_name = f"vex_{yymmdd}"
                old_path = os.path.join(dirpath, dirname)
                new_path = os.path.join(dirpath, new_name)
                print(f"Renaming {old_path} -> {new_path}")
                os.rename(old_path, new_path)

# Usage example:
root_dir = f'/Users/lgisolfi/Desktop/PRIDE_DATA_NEW/analysed_pride_data/vex/usable/converted_old_format_files'
rename_vex_folders(root_dir)
