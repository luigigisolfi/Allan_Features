import re

def extract_ground_station_dict(file_path):
    """
    Extracts the ground station configuration from a file and parses it into a dictionary.

    Args:
        file_path (str): Path to the input file containing ground station data.

    Returns:
        dict: A nested dictionary mapping station names to their key-value attributes,
              or None if no ground station section is found.
    """
    dictionary = {}
    ground_station_dict = {}
    inside_section = False
    current_station = None

    try:
        with open(file_path, 'r') as file:
            for line in file:
                # Start of a section
                if line.strip().startswith('$SITE'):
                        inside_section = True
                # Parse key-value pairs if inside the ground station section
                if inside_section:
                    # Skip the def line itself
                    if line.strip().startswith('def'):
                        match = re.search(r'def\s+(\w+);', line.strip())
                        if match:
                            current_station = match.group(1)
                            ground_station_dict[current_station] = {}

                    # End of the section
                    if line.strip().startswith('$ANTENNA'):
                        print('end')
                        current_station = None
                        break

                    # Parse key-value pairs (assumes "key = value" format)
                    if current_station and '=' in line:
                        # Use a regex to find all key-value pairs in the line
                        pairs = re.findall(r'(\w+)\s*=\s*([^\s;]+)', line)
                        for key, value in pairs:
                            if key != 'site_name' and key != 'site_position' and key != 'horizon_map_az' and  key != 'horizon_map_el':
                                ground_station_dict[current_station][key] = value
# Handle line with elev, long, and lat

        # Return the dictionary or None if the section was not found
        return ground_station_dict if ground_station_dict else None

    except FileNotFoundError:
        print(f"File not found: {file_path}")
        return None

def site_to_ID(site):
    """
    Maps a ground station name to its corresponding site ID.

    Args:
        site (str): The name of the ground station (case-insensitive).

    Returns:
        str: The site ID if the site name is found, or None otherwise.
    """
    site_to_id_mapping = {
        'CEDUNA': 'Cd',
        'HOBART12': 'Hb',
        'YARRA12M': 'Yg',
        'KATH12M': 'Ke',
        'WARK': 'Ww',
        'YAMAGU32': 'Ym',
        'TIANMA65': 'T6',
        'KUNMING': 'Km',
        'KVNUS': 'Ku',
        'BADARY': 'Bd',
        'URUMQI': 'Ur',
        'ZELENCHK': 'Zc',
        'HARTRAO': 'Hh',
        'WETTZELL': 'Wz',
        'SVETLOE': 'Sv',
        'MEDICINA': 'Mc',
        'WSTRBORK': 'Wb',
        'ONSALA60': 'On',
        'YEBES40M': 'Ys',
        'VLBA_SC': 'Sc',
        'VLBA_HN': 'Hn',
        'VLBA_NL': 'Nl',
        'VLBA_FD': 'Fd',
        'VLBA_LA': 'La',
        'VLBA_KP': 'Kp',
        'PIETOWN': 'Pt',
        'VLBA_BR': 'Br',
        'VLBA_OV': 'Ov',
        'VLBA_MK': 'Mk',
        'HART15M': 'Ht',
        'METSAHOV': 'Mh'
    }

    # Convert input to uppercase for case-insensitive matching
    site = site.upper()
    return site_to_id_mapping.get(site, None)    


def ID_to_site(site_ID):
    """
    Maps a site ID to its corresponding ground station name.

    Args:
        site_ID (str): The site ID.

    Returns:
        str: The ground station name if the site ID is found, or None otherwise.
    """
    id_to_site_mapping = {
        'Cd': 'CEDUNA',
        'Hb': 'HOBART12',
        'Yg': 'YARRA12M',
        'Ke': 'KATH12M',
        'Ww': 'WARK',
        'Ym': 'YAMAGU32',
        'T6': 'TIANMA65',
        'Km': 'KUNMING',
        'Ku': 'KVNUS',
        'Bd': 'BADARY',
        'Ur': 'URUMQI',
        'Zc': 'ZELENCHK',
        'Hh': 'HARTRAO',
        'Ht': 'HART15M',
        'Mh': 'METSAHOV',
        'Wz': 'WETTZELL',
        'Sv': 'SVETLOE',
        'Mc': 'MEDICINA',
        'Wb': 'WSTRBORK',
        'On': 'ONSALA60',
        'Ys': 'YEBES40M',
        'Sc': 'SC-VLBA',
        'Hn': 'HN-VLBA',
        'Nl': 'NL-VLBA',
        'Fd': 'FD-VLBA',
        'La': 'LA-VLBA',
        'Kp': 'KP-VLBA',
        'Br': 'BR-VLBA',
        'Ov': 'OV-VLBA',
        'Mk': 'MK-VLBA',
        'Pt': 'PIETOWN',
    }

    # Return the corresponding site name or None if the site_ID is not found
    return id_to_site_mapping.get(site_ID, None)

    
