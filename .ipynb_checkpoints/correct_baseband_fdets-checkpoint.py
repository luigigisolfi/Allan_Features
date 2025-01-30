import argparse
import shutil
import numpy as np
import pandas as pd
from astropy.time import Time, TimeDelta
import glob
import re

# Script to convert the baseband for new-type Fdets (JUICE-like)

parser = argparse.ArgumentParser()
parser.add_argument("filename", nargs='+',
                    help="spectra input file or wildcard pattern", default="check_string_for_empty")
parser.add_argument("--correct_baseband", action="store_true",  # Add flag for baseband correction
                    help="Correct bad baseband frequencies and replace them with '8412.00 MHz'")

args = parser.parse_args()

# If the argument is a wildcard pattern, use glob to get the list of files
if "*" in args.filename:
    files = glob.glob(args.filename)  # Get list of files matching the pattern
else:
    files = args.filename  # Single file

# Loop through each file and process it
for file in files:
    try:
        print(f"Processing file: {file}")
        fd = open(file)
        for ip in range(5):
            line = fd.read()

        num_columns = len(line.split(' '))
        new_format = False
        if num_columns == 6:
            new_format = False
        elif num_columns == 5:
            new_format = True

        if new_format:
            fdets = pd.read_csv(file, skiprows=4, sep = '\\s+',
                                names=['time', 'snr', 'spectral_max', 'doppler', 'residuals'], on_bad_lines='warn'
                                )

            print(fdets)
            fd = open(file, 'r')
            fdets_header = []
            for ip in range(2):
                fdets_header.append(fd.readline())
                fdets_header[0] = fdets_header[0].replace('# ','').replace('// ','')
                fdets_header[1] = fdets_header[1].replace('# ','').replace('//','')
                print(fdets_header[0])
                print(fdets_header[1])
            # If the flag is set, correct the baseband frequency in the header
            if args.correct_baseband:
                print('correcting')
                # Check for baseband frequencies using regex and replace if necessary
                header_str = fdets_header[1]
                baseband_pattern = r"Base frequency: \d{4}.\d{2} MHz"
                match = re.match(baseband_pattern, fdets_header[1])

                if not match:
                    correct_baseband_value = 8412.00 #to be changed with a function that computes it based on the mission and station
                    try:
                        corrected_header = re.sub(baseband_pattern, f"Base frequency: {correct_baseband_value} MHz", header_str) #to be changed with a function
                        fdets_header[1] = corrected_header  # Update the second header line with corrected frequency
                        print(fdets_header[1])
                    except:
                        print('nope')

            fheader = fdets_header[0] + fdets_header[1] + \
                      'Format: Modified JD   |       Time(UTC) [s]    |    Signal-to-Noise     |      Spectral max      | Freq detection [Hz] |   Doppler noise [Hz]  |\n'
            fd.close()

        else:
            fdets = pd.read_csv(file, skiprows=4,
                                sep = '\\s+',
                                names=['mjd', 'time', 'snr', 'spectral_max', 'doppler', 'residuals'], on_bad_lines='warn')

            fd = open(file, "r")
            fdets_header = []
            for ip in range(2):
                fdets_header.append(fd.readline())
            if args.correct_baseband:
                print('correct')
                # Check for baseband frequencies using regex and replace if necessary
                fdets_header[0] = fdets_header[0].replace('# ','').replace('// ','')
                fdets_header[1] = fdets_header[1].replace('# ','').replace('//','')
                print(fdets_header[0])
                print(fdets_header[1])
                baseband_pattern = r"Base frequency: \d{4}.\d{2} MHz"
                match = re.match(baseband_pattern, fdets_header[1])

                if not match:
                    print('no match')
                    correct_baseband_value = 8412.00 #to be changed with a function that computes it based on the mission and station
                    try:
                        corrected_header = f"# Base frequency: {correct_baseband_value} MHz\n" #to be changed with a function
                        fdets_header[1] = corrected_header  # Update the second header line with corrected frequency
                        print(fdets_header[1])
                    except:
                        print('nope')

            fheader = '# ' + fdets_header[0] + fdets_header[1] + \
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
