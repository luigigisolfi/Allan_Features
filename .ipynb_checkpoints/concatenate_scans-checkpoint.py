# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.16.4
#   kernelspec:
#     display_name: Python (tudat-bundle)
#     language: python
#     name: tudat-bundle_updated
# ---

# +
import glob
import re
from collections import defaultdict


# Define the base pattern for your files (without the station part)
file_pattern = r"Fdets\.\w+(\d{4})\.(\d{2})\.(\d{2})\.(\w{2})\.(\d{4})\.r2i\.txt"

# Find all matching files (using the glob pattern)
files = glob.glob('/Users/lgisolfi/Desktop/data_archiving-1.0/dataset/juice/juice_cruise/EC094A/Fdets.*.r2i.txt')  # Get all files first, or you can limit this further

print(files)
# Create a dictionary to store file handles for each station
station_files = defaultdict(list)

# Loop over each file to group them by station
for file in files:
    # Use regex to extract station name from the filename
    match = re.search(file_pattern, file)
    if match:
        print('match')
        date =f'{match.group(1)}.{match.group(2)}.{match.group(3)}'
        station_name = match.group(4)  # Extract station name (2 letters, e.g., Ef)
        scan_number = match.group(5)
        
        # Add the file and its scan number to the corresponding station's list
        station_files[station_name].append((file, int(scan_number)))  # Store as tuple (file, scan_number)

# Now, create a Complete file for each station
for station, files in station_files.items():
    # Sort the files based on the scan number (second element of the tuple)
    files.sort(key=lambda x: x[1])  # Sorting by scan_number (second element of tuple)
    # Construct the output filename based on the station name
    output_filename = f"/Users/lgisolfi/Desktop/data_archiving-1.0/dataset/juice/juice_cruise/EC094A/Complete/Fdets.jc{date}.{station}.Complete.r2i.txt"
    # Open the output file in write mode
    with open(output_filename, 'w') as output_file:
        # Flag to check if the header has already been written
        header_written = False

        # Loop over each file for the current station, sorted by scan number
        for file, _ in files:
            # Open the file in read mode
            with open(file, 'r') as f:
                # Read the lines of the file
                lines = f.readlines()

                # Skip the header in all files except the first one
                if not header_written:
                    # Write the first file's header
                    output_file.write(''.join(lines[:5]))  # Assuming the first 3 lines are the header
                    header_written = True

                # Write the rest of the content, skipping the header lines
                output_file.write(''.join(lines[5:]))  # Skip the first 3 lines (the header)

    print(f"Created {output_filename}")

# -


