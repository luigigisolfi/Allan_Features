import os
import subprocess

"""
This code was used to iteratively convert all old_format vex data to the new fdets format. 
"""
def run_old2new_conversion(base_folder: str, skip_folders=None) -> None:
    """
    Runs the 'old2new_fdets.py' script on all subfolders inside the specified base folder,
    excluding specified folders, passing all .txt files in a subpath as input.

    Parameters
    ----------
    base_folder : str
        The root directory containing 'vex_YYMM' folders with data to process.
    skip_folders : list[str], optional
        List of folder names to skip during processing (default is None).

    Behavior
    --------
    - Lists all subdirectories inside base_folder.
    - Skips folders listed in skip_folders.
    - For each remaining folder, lists its subdirectories.
    - For each subdirectory, constructs a command line to run
      'old2new_fdets.py' with all .txt files inside 'input/single/' subfolder.
    - Runs the command using subprocess.run with shell=True.
    - Prints nothing explicitly, but subprocess.run output will be visible on console.

    Notes
    -----
    - Assumes that the script 'old2new_fdets.py' is in the PATH or current directory.
    - Uses glob wildcard '*.txt' in the command to process all .txt files in that folder.
    - Uses shell=True, which may have security implications if inputs are not trusted.
    """
    if skip_folders is None:
        skip_folders = []

    vex_year_month_list = [
        name for name in os.listdir(base_folder)
        if os.path.isdir(os.path.join(base_folder, name))
    ]

    for vex_year_month in vex_year_month_list:
        if vex_year_month in skip_folders:
            continue

        parent_folder = os.path.join(base_folder, vex_year_month)

        folder_names = [
            name for name in os.listdir(parent_folder)
            if os.path.isdir(os.path.join(parent_folder, name))
        ]

        for name_folder in folder_names:
            input_path = os.path.join(parent_folder, name_folder, "input", "single", "*.txt")
            command = f"python old2new_fdets.py {input_path} --experiment_name {name_folder}"
            subprocess.run(command, shell=True)

# Example usage
old_format_files_folder = "/Users/lgisolfi/Desktop/data_archiving-2.0/vex/usable/converted_old_format_files/"
skip_list = ["vex_1401", "vex_1402"]
run_old2new_conversion(old_format_files_folder, skip_list)
