import argparse
import shutil
import numpy as np
import pandas as pd
from astropy.time import Time, TimeDelta
import glob
import os
import re
from pride_characterization_library  import PrideDopplerCharacterization

pride = PrideDopplerCharacterization()
process_fdets = pride.ProcessFdets()
utilities = pride.Utilities()
analysis = pride.Analysis(process_fdets, utilities)
process_vex_files = pride.ProcessVexFiles(process_fdets, utilities)
format_fdets = pride.FormatFdets(process_fdets, utilities, process_vex_files)


"""
Script for processing spectra input files and optionally correcting baseband frequencies.

This script takes one or more spectra input files (or wildcard patterns to match files) and processes them to 
convert their formats or correct bad baseband frequencies based on the specified arguments.

Command-line Arguments:
-----------------------
filename : str
    Required. Input file(s) or wildcard pattern specifying spectra files to process. 
    Multiple files can be specified using space-separated filenames or a wildcard pattern (e.g., '*.txt').

--correct_baseband : flag (optional)
    If provided, enables baseband frequency correction for the specified input files. 
    Requires the --experiment_name argument.

--experiment_name : str (optional)
    Name of the experiment to be used when correcting baseband frequencies. 
    Must be specified when --correct_baseband is used.

Features:
---------
1. Handles both single files and wildcard patterns for batch processing.
2. Creates backup copies of the original files with the `.old` extension.
3. Determines the file format (old or new) based on the number of columns in the data:
    - Old format: 6 columns.
    - New format: 5 columns.
4. Converts files in the old format to a new timestamped format, preserving key data fields:
    - Signal-to-Noise Ratio (SNR)
    - Spectral Max
    - Doppler
    - Residual Fit
5. Performs baseband frequency correction using an external script (`correct_baseband_fdets.py`) 
   if the `--correct_baseband` flag is provided.

Example Usage:
--------------
1. Process files without baseband correction:
    ```
    python script.py file1.txt file2.txt
    ```

2. Process all files matching a wildcard pattern:
    ```
    python script.py "*.txt"
    ```

3. Perform baseband frequency correction for a set of files:
    ```
    python old2new_fdets.py file1.txt file2.txt --correct_baseband --experiment_name EXP001
    ```

Requirements:
-------------
- Python packages: argparse, pandas, numpy, astropy
- External script: `correct_baseband_fdets.py` must be accessible in the current working directory.
- Input files must have a valid spectra format (6 or 5 columns).

Output:
-------
- Converted files are saved in the same location as the input files, overwriting the originals.
- Original files are renamed with the `.old` extension for backup.
- Baseband correction results are handled by the external script if enabled.

"""

parser = argparse.ArgumentParser()
parser.add_argument("filename", nargs='+',
                    help="spectra input file or wildcard pattern", default="check_string_for_empty")
parser.add_argument("--correct_baseband", action="store_true",  # Add flag for baseband correction
                    help="Correct bad baseband frequencies and replace them with '8412.00 MHz'")
parser.add_argument("--experiment_name", required= False,  # Add flag for baseband correction
                    help="Name of the experiment file to correct bad baseband frequencies'")

args = parser.parse_args()


# Validate that experiment_name is provided when correct_baseband is set
if args.correct_baseband and not args.experiment_name:
    parser.error("--experiment_name is required when --correct_baseband is specified.")

# If the argument is a wildcard pattern, use glob to get the list of files
if "*" in args.filename:
    files = glob.glob(args.filename)  # Get list of files matching the pattern
else:
    files = args.filename  # Single file

if not args.correct_baseband:
    old_files = [shutil.copy(file, f'{file}'.replace("txt", "old")) for file in files]

for filename in files:
    print(filename)
    mission_name = re.sub(r'\d', '', filename.split('/')[-1].split('.')[1])  # Remove all digits (if present) from mission name (ex. in mex2013, remove 2013)

if args.correct_baseband:
    print('Correcting Baseband Fdets...')
    file_list = ' '.join(files) # convert the list into a long string of files to correct baseband frequencies

    if args.experiment_name:
        experiment_name = args.experiment_name
        if mission_name.lower() != utilities.get_mission_from_experiment(experiment_name):
            # Filter experiments where the mission name matches
            experiments_list_for_mission = [
                exp_key
                for exp_key, exp_details in utilities.experiments.items()
                if exp_details["mission_name"].lower() == mission_name.lower()
            ]
            print(f'Wrong experiment name: {experiment_name} for mission: {mission_name}.\nList of available {mission_name} experiments: {experiments_list_for_mission}.\nAborting...')
            exit()
        else:
            os.system(f'python correct_baseband_fdets.py {file_list} --experiment_name {experiment_name}')
    else:
        print(f'Argument --experiment_name not specified. Always specify it when using --correct_baseband argument.\n Aborting...')
        exit()

for file in files:
    delete_flag = False
    try:
        print(f"Processing file: {file}")

        # Safely open and process file
        with open(file) as fd:
            lines = [fd.readline() for _ in range(5)]

        num_columns = len(lines[-1].split())
        new_format = num_columns == 5

        if not new_format:
            fdets = pd.read_csv(
                file,
                skiprows=4,
                sep='\\s+',
                names=['mjd', 'time', 'snr', 'spectral_max', 'doppler', 'residuals'],
                on_bad_lines='warn'
            )

            with open(file, "r") as fd:
                fdets_header = [fd.readline() for _ in range(2)]


            fdets_header[0] = re.sub(r'^[^\w]*', '', fdets_header[0])
            fdets_header[1] = re.sub(r'^[^\w]*', '', fdets_header[1])

            fheader = (
                    fdets_header[0] +
                    fdets_header[1] +
                    'Format: UTC Time     |    Signal-to-Noise     |      Spectral max      | '
                    'Freq detection [Hz] |   Doppler noise [Hz]  |\n'
            )

            epoch = Time(fdets.mjd, format='mjd', scale='utc')
            seconds = TimeDelta(fdets.time, format='sec')
            time = epoch + seconds

            data = {
                'timestamp': time.isot,
                'SNR': fdets.snr,
                'Spectral Max': fdets.spectral_max,
                'Doppler': fdets.doppler,
                'Residual fit': fdets.residuals
            }

            np.savetxt(file, pd.DataFrame(data=data),
                       fmt='%15s %.18e %.18e %21.12f %+.16e', header=fheader)

            print(f'Converted File Saved to: {file}\n')

        print(f"Processing completed for: {file}\n")

    except Exception as e:
        delete_flag = True
        print(f"Error processing {file}: {e}")

    if delete_flag:
        try:
            print(f"Deleting problematic file: {file}")
            os.remove(file)
            print(f"Removed: {file}")
        except Exception as del_error:
            print(f"Failed to delete {file}: {del_error}")
