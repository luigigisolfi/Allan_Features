"""
Experience with the dataset showed that some stations might have suffered from a constant time offset in the time tag, while still keeping the correct header.
This script
1) checks for the observation time in the header
2) checks that the day in the time tag is the same as the header
3) corrects if needed.

Note:
Venus express data (especially Wetzell data) heavily suffered from this.
"""
import glob
import os
from datetime import datetime
import re

def correct_time_tags(input_file):
    corrected_lines = []
    observed_date = None

    with open(input_file, 'r') as f:
        lines = f.readlines()

    for line in lines:
        if line.startswith('#'):
            corrected_lines.append(line)
            if 'Observation conducted on' in line:
                match = re.search(r'(\d{4})\.(\d{2})\.(\d{2})', line)
                if match:
                    observed_date = datetime.strptime(match.group(0), "%Y.%m.%d").date()
        else:
            try:
                parts = line.strip().split()
                timestamp_str = parts[0]
                timestamp = datetime.strptime(timestamp_str, "%Y-%m-%dT%H:%M:%S.%f")
                if observed_date and timestamp.date() != observed_date:
                    corrected_timestamp = timestamp.replace(
                        year=observed_date.year,
                        month=observed_date.month,
                        day=observed_date.day
                    )
                    parts[0] = corrected_timestamp.isoformat(timespec='milliseconds')
                    corrected_line = " ".join(parts) + "\n"
                    corrected_lines.append(corrected_line)
                else:
                    corrected_lines.append(line)
            except Exception:
                corrected_lines.append(line)

    with open(input_file, 'w') as f:
        f.writelines(corrected_lines)

def apply_correction_to_files(wildcard_pattern):
    matched_files = glob.glob(wildcard_pattern, recursive=True)
    print(f"Found {len(matched_files)} matching files.")
    for file in matched_files:
        print(f"Correcting: {file}")
        correct_time_tags(file)

# Example usage:
# Correct all files under the 'data/' directory containing 'Wz' in the filename
apply_correction_to_files(f"/Users/lgisolfi/Desktop/data_archiving-2.0/vex/usable/converted_old_format_files/**/single/*Wz*.txt")
