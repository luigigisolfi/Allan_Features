import os
import csv
import sys
from datetime import datetime

def extract_data_from_txt(txt_file_path):
    with open(txt_file_path, 'r') as file:
        lines = file.readlines()

    # Extracting the required data
    file_name = os.path.basename(txt_file_path)

    # Modifying creation_date_time format
    creation_date_time = datetime.now().strftime("%Y-%m-%dT%H:%M:%S")
    
    # Modifying start_date_time and stop_date_time format
    start_date_time = lines[4][:23] + 'Z'  # Row 5, columns 1-23 with 'Z' added at the end
    stop_date_time = lines[-1][:23] + 'Z'  # Last row, columns 1-23 with 'Z' added at the end

    return {
        "file_name": file_name,
        "creation_date_time": creation_date_time,
        "start_date_time": start_date_time,
        "stop_date_time": stop_date_time
    }

def get_file_size(file_path):
    """
    Returns the size of the file in bytes.
    
    :param file_path: Path to the file
    :return: Size of the file in bytes
    """
    return os.path.getsize(file_path)

def find_row_offset(filename, target_row):
    offset = 0
    with open(filename, 'rb') as file:
        for current_row, line in enumerate(file, start=1):
            if current_row == target_row:
                return offset
            offset += len(line)
    raise ValueError(f"Row {target_row} does not exist in the file")

def calculate_record_rows(file_path, header_rows):
    # Read the CSV file and count the rows
    with open(file_path, 'r') as csvfile:
        reader = csv.reader(csvfile)
        rows = list(reader)
    
    total_rows = len(rows)
    
    # Identify empty rows at the end
    empty_rows_at_end = 0
    for row in reversed(rows):
        if not any(row):  # Check if the row is empty
            empty_rows_at_end += 1
        else:
            break
    
    # Calculate the record length
    record_rows = total_rows - empty_rows_at_end - header_rows
    return record_rows
def get_row_details(filename, target_row):
    with open(filename, 'rb') as file:
        for current_row, line in enumerate(file, start=1):
            if current_row == target_row:
                # Keep the original line length
                record_length = len(line)
                
                # Split the line into columns while ignoring multiple spaces
                columns = line.split()
                
                # Calculate starting position and length of each column
                start_positions = []
                lengths = []
                index = 0
                for col in columns:
                    start_pos = line.index(col, index)
                    start_positions.append(start_pos)
                    col_length = len(col)
                    lengths.append(col_length)
                    index = start_pos + col_length  # Update index for next search

                return record_length, list(zip(start_positions, lengths))
    
    raise ValueError(f"Row {target_row} does not exist in the file")

def write_csv(csv_file_path, data):
    with open(csv_file_path, 'w', newline='') as csvfile:
        fieldnames = [
            "file_name",
            "creation_date_time",
            "start_date_time",
            "stop_date_time",
            "file_size",
            "record_rows",
            "record_length",
            "row_offset",
            "field1_location",
            "field2_location",
            "field3_location",
            "field4_location",
            "field5_location",
            "field1_length",
            "field2_length",
            "field3_length",
            "field4_length",
            "field5_length",
            ]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        writer.writerow(data)

def process_directory(directory):
    for root, _, files in os.walk(directory):
        for file_name in files:
            if file_name.endswith('.txt'):
                txt_file_path = os.path.join(root, file_name)
                
                try:
                    # Extract data from the text file
                    data = extract_data_from_txt(txt_file_path)
                    print(data)
                    # Check if start_date_time contains a space
                    if ' ' in data["start_date_time"]:
                        print(f"Skipping file {file_name} due to invalid start_date_time format.")
                        continue
                    
                    # Add file size to the data
                    data["file_size"] = get_file_size(txt_file_path)
                    print(data["file_size"])
                    
                    # Given number of header rows
                    header_rows = 4

                    # Calculate the record rows
                    data["record_rows"] = calculate_record_rows(txt_file_path, header_rows)
                    print(data["record_rows"])
                    target_row = 5
                    data["row_offset"] = find_row_offset(txt_file_path, target_row)
                    print(data["row_offset"])
                    # Call get_row_details and store the results in the data dictionary
                    record_length, field_info = get_row_details(txt_file_path, target_row)
                    data["record_length"] = record_length

                    # Assign start positions and lengths to specific keys in the data dictionary
                    for i, (start, length) in enumerate(field_info[:5], start=1):
                        data[f"field{i}_location"] = start + 1
                        data[f"field{i}_length"] = length

                    # Generate the output CSV file path
                    csv_file_path = os.path.splitext(txt_file_path)[0] + '.csv'
                    
                    # Write the extracted data to the CSV file
                    write_csv(csv_file_path, data)
                    # print(f"CSV file generated successfully: {csv_file_path}")
                
                except IndexError:
                    print(f"Skipping file {file_name} due to insufficient lines in the file.")
                except Exception as e:
                    print(f"Error processing file {file_name}: {e}")

def clean_csv_files(directory):
    for root, _, files in os.walk(directory):
        for file_name in files:
            if file_name.endswith('.csv'):
                csv_file_path = os.path.join(root, file_name)
                try:
                    os.remove(csv_file_path)
                    # print(f"Deleted file: {csv_file_path}")
                except Exception as e:
                    print(f"Error deleting file {csv_file_path}: {e}")

if __name__ == "__main__":
    if len(sys.argv) < 2 or len(sys.argv) > 3:
        print("Usage: python generate_labels.py <directory_path> [clean]")
    else:
        directory_path = sys.argv[1]
        
        if len(sys.argv) == 3 and sys.argv[2] == 'clean':
            clean_csv_files(directory_path)
            print("Cleanup complete.")
        else:
            process_directory(directory_path)

# python3 2_generate_metadata.py ./dataset/
# python3 2_generate_metadata.py ./dataset/ clean