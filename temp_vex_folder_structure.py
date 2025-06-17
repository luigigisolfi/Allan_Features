import os
import shutil

# Paths
old_base = f'/Users/lgisolfi/Desktop/PRIDE_DATA_NEW/analysed_pride_data/vex/usable/converted_old_format_files'
new_base = f'/Users/lgisolfi/Desktop/PRIDE_DATA_NEW/analysed_pride_data/vex'

def move_vex_yymm_folders(old_base_path: str, new_base_path: str) -> None:
    """
    Move all subdirectories (named with 'vex_yymm' format) from the old base directory
    to a new base directory.

    Parameters
    ----------
    old_base_path : str
        The source directory containing the 'vex_yymm' folders to move.
    new_base_path : str
        The destination directory where the 'vex_yymm' folders will be moved.

    Behavior
    --------
    - Iterates over all entries in the old_base_path.
    - Checks if the entry is a directory; if not, skips it.
    - Moves the entire directory to the new_base_path, preserving the folder name.
    - Prints progress of each folder being moved.
    - After moving all folders, attempts to remove the old base directory if it is empty.
    - Prints whether the removal was successful or not.

    Raises
    ------
    OSError
        If the old_base_path cannot be removed due to it not being empty.
    """
    for yymm_folder in os.listdir(old_base_path):
        old_yymm_path = os.path.join(old_base_path, yymm_folder)

        if not os.path.isdir(old_yymm_path):
            continue

        new_yymm_path = os.path.join(new_base_path, yymm_folder)

        print(f"Moving {old_yymm_path} -> {new_yymm_path}")
        shutil.move(old_yymm_path, new_yymm_path)

    try:
        os.removedirs(old_base_path)  # Removes only if empty
        print(f"Removed empty directory {old_base_path}")
    except OSError:
        print(f"Could not remove {old_base_path} (not empty)")

# Call the function with the specified paths
move_vex_yymm_folders(old_base, new_base)
