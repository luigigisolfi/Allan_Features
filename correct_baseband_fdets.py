import argparse
import shutil
import pandas as pd
import glob
import re
from pride_characterization_library  import PrideDopplerCharacterization

pride = PrideDopplerCharacterization()
process_fdets = pride.ProcessFdets()
utilities = pride.Utilities()
analysis = pride.Analysis(process_fdets, utilities)
process_vex_files = pride.ProcessVexFiles(process_fdets, utilities)
format_fdets = pride.FormatFdets(process_fdets, utilities, process_vex_files)

"""
This script processes spectral data files, checks for proper baseband frequencies in headers,
and applies corrections when necessary. It supports handling single files or wildcard patterns
for batch processing.

Main functionalities:
- Parses command-line arguments to accept input files and experiment name.
- Copies original files to backup with a ".old" extension.
- Reads and validates the number of columns in data files.
- Corrects baseband frequencies in file headers if missing or incorrect.
- Saves updated headers back to the original files.

Functions used:
- `get_baseband_frequency`: Placeholder for a function that computes or retrieves the correct
  baseband frequency based on mission and experiment details.

Dependencies:
- `argparse`, `glob`, `shutil`, `re`, `pandas` (as `pd`)

Arguments:
- `filename`: One or more input files or wildcard patterns specifying spectral data files.
- `--experiment_name`: Mandatory argument specifying the experiment name for baseband corrections.

Outputs:
- Updates the input files in place with corrected headers, preserving the data format.

########################################################################################################################
 NOTE: please note some stations in the fdets are not specified in the $FREQ block of the vex files, 
       but they are found in the $IF block. The current implementation does not consider stations in the $IF block.
########################################################################################################################

"""

parser = argparse.ArgumentParser()
parser.add_argument("filename", nargs='+',
                    help="spectra input file or wildcard pattern", default="check_string_for_empty")
parser.add_argument("--experiment_name", required="True",  # Add flag for baseband correction
                    help="Name of the experiment file to correct bad baseband frequencies'")
args = parser.parse_args()

if args.experiment_name:
    experiment_name =  args.experiment_name
else:
    print(f'Argument: --experiment_name needs to be specified.\n Aborting...')

# If the argument is a wildcard pattern, use glob to get the list of files
if "*" in args.filename:
    files = glob.glob(args.filename)  # Get list of files matching the pattern
else:
    files = args.filename  # Single file

old_files = [shutil.copy(file, f'{file}'.replace("txt", "old")) for file in files]

# Loop through each file and process it
for file in files:
    try:
        mission_name = re.sub(r'\d', '', file.split('/')[-1].split('.')[1])

        print(experiment_name)

        if mission_name.lower() != utilities.get_mission_from_experiment(experiment_name):
            print(f'Wrong experiment name: {experiment_name} for mission: {mission_name}. Aborting...')
            exit()

        process_vex_files.get_baseband_frequencies_file(mission_name, experiment_name) # get frequency file
        format_fdets.assign_missing_baseband_frequencies(mission_name, experiment_name)

        print(f"Processing file: {file}")
        fd = open(file)
        line = fd.readlines()[4]

        num_columns = len(line.split())
        new_format = False
        if num_columns == 6:
            new_format = False
        elif num_columns == 5:
            new_format = True

        if new_format:
            fdets = pd.read_csv(file, skiprows=4, sep = '\\s+',
                                names=['time', 'snr', 'spectral_max', 'doppler', 'residuals'], on_bad_lines='warn'
                                )

            fd = open(file, 'r')
            fdets_header = []
            for ip in range(2):
                fdets_header.append(fd.readline())

            # Removing non-alphanumeric characters from header
            fdets_header[0] = fdets_header[0] = re.sub(r'^[^\w]*', '', fdets_header[0])
            fdets_header[1] = fdets_header[1] = re.sub(r'^[^\w]*', '', fdets_header[1])

            # Check for baseband frequencies using regex and replace if necessary
            header_str = fdets_header[1]
            #baseband_patterns = [r"Base frequency: \d+(.\d+)? MHz", r"Base frequency: \d+(.\d+)? MHz dF: \d+(.\d+) Hz dT: \d+(.\d+) s"]
            baseband_patternssss = [r"Base frequency: (?!0\.00)\d+(\.\d+)? MHz(?: dF: \d+(\.\d+)? Hz dT: \d+(\.\d+)? s)?"]
            match_flag = False

            try:
                match = re.match(baseband_patternssss, fdets_header[1])
                match_flag = True
            except:
                print('not match, correct_baseband_frequencies')
                correct_baseband_value = process_vex_files.get_baseband_frequency(mission_name, experiment_name, file)
                print('correct baseband', correct_baseband_value)
                try:
                    corrected_header = f"Base frequency: {correct_baseband_value} MHz\n"
                    fdets_header[1] = corrected_header  # Update the second header line with corrected frequency
                except:
                    continue

            fheader = '# '  + fdets_header[0] + '# ' + fdets_header[1] + \
                      '# ' 'Format: UTC Time    |    Signal-to-Noise     |      Spectral max      | Freq detection [Hz] |   Doppler noise [Hz]  |\n#\n'
            fd.close()

        else:
            fdets = pd.read_csv(file, skiprows=4,
                                sep = '\\s+',
                                names=['mjd', 'time', 'snr', 'spectral_max', 'doppler', 'residuals'], on_bad_lines='warn')

            fd = open(file, "r")
            fdets_header = []
            for ip in range(2):
                fdets_header.append(fd.readline())

            # Removing non-alphanumeric characters from header
            fdets_header[0] = fdets_header[0] = re.sub(r'^[^\w]*', '', fdets_header[0])
            fdets_header[1] = fdets_header[1] = re.sub(r'^[^\w]*', '', fdets_header[1])
            baseband_patterns = [r"Base frequency: (?!0\.00)\d+(\.\d+)? MHz(?: dF: \d+(\.\d+)? Hz dT: \d+(\.\d+)? s)?"]

            match_flag = False
            for baseband_pattern in baseband_patterns:
                try:
                    match = re.match(baseband_patterns, fdets_header[1])
                    match_flag = True
                    print('having a match', match)
                except:
                    continue

            if not match_flag:
                correct_baseband_value = process_vex_files.get_baseband_frequency(mission_name, experiment_name, file)
                print('correct baseband', correct_baseband_value)
                try:
                    corrected_header = f"Base frequency: {correct_baseband_value} MHz\n"
                    fdets_header[1] = corrected_header  # Update the second header line with corrected frequency
                except:
                    continue

            fheader = '# ' + fdets_header[0] + '# ' + fdets_header[1] + \
                      '# ' + 'Format: Modified JD   |       Time(UTC) [s]    |    Signal-to-Noise     |      Spectral max      | Freq detection [Hz] |   Doppler noise [Hz]  |\n#\n'
            fd.close()

        with open(file, 'r') as fd:
            content = fd.readlines()[4:]  # Skip the first 4 lines (header) that you've already processed

        # Save the updated file with the new header
        with open(file, 'w') as fd:
            fd.write(fheader)  # Write the new header
            fd.writelines(content)  # Write the rest of the data unchanged

    except Exception as e:
        print(f'Could not convert file: {e}')
