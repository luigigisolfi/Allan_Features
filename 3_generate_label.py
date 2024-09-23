import csv
import os
import json
import sys

def load_config(config_file):
    with open(config_file, 'r') as f:
        config = json.load(f)
    return config

def generate_label_from_csv(csv_file, template_content, mappings):
    # Read metadata from CSV file
    with open(csv_file, 'r') as metadata_file:
        reader = csv.DictReader(metadata_file)
        row = next(reader)  # Assuming each CSV file has one row of metadata
        
        # Fill in the template with metadata values
        label_content = template_content.replace('${file_name}', row[mappings['file_name']])
        label_content = label_content.replace('${file_size}', row[mappings['file_size']])
        label_content = label_content.replace('${row_offset}', row[mappings['row_offset']])
        label_content = label_content.replace('${record_rows}', row[mappings['record_rows']])
        label_content = label_content.replace('${record_length}', row[mappings['record_length']])
        label_content = label_content.replace('${creation_date_time}', row[mappings['creation_date_time']])
        label_content = label_content.replace('${start_date_time}', row[mappings['start_date_time']])
        label_content = label_content.replace('${stop_date_time}', row[mappings['stop_date_time']])
        label_content = label_content.replace('${field1_location}', row[mappings['field1_location']])
        label_content = label_content.replace('${field2_location}', row[mappings['field2_location']])
        label_content = label_content.replace('${field3_location}', row[mappings['field3_location']])
        label_content = label_content.replace('${field4_location}', row[mappings['field4_location']])
        label_content = label_content.replace('${field5_location}', row[mappings['field5_location']])
        label_content = label_content.replace('${field1_length}', row[mappings['field1_length']])
        label_content = label_content.replace('${field2_length}', row[mappings['field2_length']])
        label_content = label_content.replace('${field3_length}', row[mappings['field3_length']])
        label_content = label_content.replace('${field4_length}', row[mappings['field4_length']])
        label_content = label_content.replace('${field5_length}', row[mappings['field5_length']])
        
        # Generate output filename with .lblx extension in the same directory as csv_file
        csv_file_dir = os.path.dirname(csv_file)
        csv_file_name = os.path.basename(csv_file)
        filename_base, _ = os.path.splitext(csv_file_name)
        label_filename = f"{filename_base}.lblx"
        label_file_path = os.path.join(csv_file_dir, label_filename)
        
        # Save the filled template as a new XML file in the same directory as the CSV file
        with open(label_file_path, 'w') as label_file:
            label_file.write(label_content)
        
        # print(f"Generated label file: {label_file_path}")

def clean_generated_lblx_files(directory):
    # Iterate over each file in the directory and its subdirectories
    for root, _, files in os.walk(directory):
        for file_name in files:
            if file_name.endswith('.lblx'):
                file_path = os.path.join(root, file_name)
                os.remove(file_path)
                # print(f"Deleted: {file_path}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python generate_labels.py <config_file> [clean]")
        sys.exit(1)
    
    config_file = sys.argv[1]
    config = load_config(config_file)
    
    # Check if the second argument is 'clean' to delete generated .lblx files
    if len(sys.argv) == 3 and sys.argv[2] == 'clean':
        clean_generated_lblx_files(config['output_dir'])
        print("Cleanup complete.")
        sys.exit(0)
    
    # Read template file
    with open(config['template'], 'r') as template_file:
        template_content = template_file.read()
    
    dataset_folder = config['output_dir']  # Use output_dir from config
    mappings = config['mappings']
    
    # Iterate over each CSV file in the dataset folder and its subdirectories
    for root, _, files in os.walk(dataset_folder):
        for file_name in files:
            if file_name.endswith('.csv'):
                csv_file = os.path.join(root, file_name)
                
                # Generate label file from current CSV file
                generate_label_from_csv(csv_file, template_content, mappings)
    
    print("Label files generated successfully.")

# python3 3_generate_label.py config.json
# python3 3_generate_label.py config.json clean