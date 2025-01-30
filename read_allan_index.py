def read_allan_index(file_path):
    """
    Reads a file and extracts the value of 'allan_index'.

    Args:
        file_path (str): Path to the text file.

    Returns:
        float: The value of 'allan_index', or None if not found.
    """
    try:
        with open(file_path, 'r') as file:
            for line in file:
                if line.startswith('allan_index'):
                    # Split the line and extract the value
                    key, value = line.split('=')
                    return float(value.strip())
    except FileNotFoundError:
        print(f"File not found: {file_path}")
    except ValueError:
        print(f"Error parsing 'allan_index' value in file: {file_path}")
    return None