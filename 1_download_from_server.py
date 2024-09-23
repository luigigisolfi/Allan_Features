import subprocess
import os
import shutil
import sys

def download_txt_files(url='https://space.phys.utas.edu.au/docs/juice/', destination='./dataset'):
    """
    Downloads .txt files recursively from a given URL and saves them to the specified destination directory.

    Parameters:
    url (str): The URL of the directory to download. Default is 'https://space.phys.utas.edu.au/docs/juice/'.
    destination (str): The destination directory to save the downloaded files. Default is './dataset'.
    
    Returns:
    None
    """
    # Create the destination directory if it does not exist
    if not os.path.exists(destination):
        os.makedirs(destination)

    # Define the command with recursive and file type inclusion options
    command = [
        'wget', 
        '-r',        # Recursive download
        '-np',       # No parent
        '-nH',       # No host directories
        '--cut-dirs=1',  # Skip the first directory component
        '-A', 'txt', # Accept only .txt files
        '-P', destination,  # Destination directory
        url
    ]

    # Execute the command
    result = subprocess.run(command, capture_output=True, text=True)

    # Check the result
    if result.returncode == 0:
        print("Directory and subdirectories with .txt files downloaded successfully")
        print("Output:", result.stdout)
    else:
        print("Command failed with return code", result.returncode)
        print("Error:", result.stderr)

def delete_directory(destination='./dataset'):
    """
    Deletes the specified destination directory.

    Parameters:
    destination (str): The destination directory to delete. Default is './dataset'.
    
    Returns:
    None
    """
    if os.path.exists(destination):
        shutil.rmtree(destination)
        print(f"Deleted the directory: {destination}")
    else:
        print(f"Directory does not exist: {destination}")

if __name__ == "__main__":
    url = 'https://space.phys.utas.edu.au/docs/juice/'
    destination_path = './dataset'
    delete = False

    if len(sys.argv) > 1:
        url = sys.argv[1]
    if len(sys.argv) > 2:
        destination_path = sys.argv[2]
    if len(sys.argv) > 3 and sys.argv[3] == 'delete':
        delete = True

    if delete:
        delete_directory(destination_path)
    else:
        download_txt_files(url, destination_path)