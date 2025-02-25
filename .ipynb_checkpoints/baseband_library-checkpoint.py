import re
from collections import defaultdict
import os
import numpy as np

"""
Define dictionary for 
1) a bunch of radio-tracking experiments experiments 
2) all missions present in the Mas Said table file.
########################################################################################################################
 NOTE: please note there is many more vex files for MEX tracking tests in the VEX FILES folder provided by G. Cimo'.
 NOTE: please note there is various vex files associated with the gr035 experiment in the VEX FILES folder provided by G. Cimo'.
 NOTE: please note some stations in the fdets are not specified in the $FREQ block of the vex files, 
       but they are found in the $IF block. The current implementation does not consider stations in the $IF block.
########################################################################################################################
"""

experiments = {
    "gr035": {
        "mission_name": "mex",
        "vex_file_name": "gr035.vix",
        "exper_description": "mars_express tracking test",
        "exper_nominal_start": "2020y053d01h30m00s",
        "exper_nominal_stop": "2020y053d03h00m00s"
    },

    "m0303": {
        "mission_name": "mex",
        "vex_file_name": "m0303.vex",
        "exper_description": "mars express tracking test",
        "exper_nominal_start": "2010y062d20h00m00s",
        "exper_nominal_stop": "2010y062d21h59m00s"
    },

    "m0325": {
        "mission_name": "mex",
        "vex_file_name": "m0325.vex",
        "exper_description": "mars express tracking test",
        "exper_nominal_start": "2012y085d13h00m00s",
        "exper_nominal_stop": "2012y085d13h59m00s"
    },

    "m0327": {
        "mission_name": "mex",
        "vex_file_name": "m0327.vex",
        "exper_description": "mars express tracking test",
        "exper_nominal_start": "2012y087d01h30m00s",
        "exper_nominal_stop": "2012y087d02h49m00s"
    },

    "m0403": {
        "mission_name": "mex",
        "vex_file_name": "m0403.vex",
        "exper_description": "mars express tracking test",
        "exper_nominal_start": "2012y087d01h30m00s",
        "exper_nominal_stop": "2012y087d02h49m00s"
    },

    "ed045a": {
        "mission_name": "min",
        "vex_file_name": "ed045a.vix",
        "exper_description": "min tracking",
        "exper_nominal_start": "2020y053d01h30m00s",
        "exper_nominal_stop": "2020y053d03h00m00s"
    },
    "ed045c": {
        "mission_name": "min",
        "vex_file_name": "ed045c.vix",
        "exper_description": "min tracking",
        "exper_nominal_start": "2020y150d08h00m00s",
        "exper_nominal_stop": "2020y150d09h30m00s"
    },
    "ed045d": {
        "mission_name": "min",
        "vex_file_name": "ed045d.vix",
        "exper_description": "min tracking",
        "exper_nominal_start": "2020y151d08h30m00s",
        "exper_nominal_stop": "2020y151d10h00m00s"
    },
    "ed045e": {
        "mission_name": "min",
        "vex_file_name": "ed045e.vix",
         "exper_description": "min tracking",
        "exper_nominal_start": "2020y295d02h45m00s",
        "exper_nominal_stop": "2020y295d04h15m00s"
    },
    "ed045f": {
        "mission_name": "min",
        "vex_file_name": "ed045f.vix",
        "exper_description": "min tracking",
        "exper_nominal_start": "2020y296d02h45m00s",
        "exper_nominal_stop": "2020y296d04h15m00s"
    },

    "ec094a": {
        "mission_name": "juice",
        "vex_file_name": "ec094a.vex",
        "exper_description": "JUICE tracking",
        "exper_nominal_start": "2023y292d14h00m00s",
        "exper_nominal_stop": "2023y292d16h00m00s"
    },

    "ec094b": {
        "mission_name": "juice",
        "vex_file_name": "ec094b.vex",
        "exper_description": "JUICE tracking",
        "exper_nominal_start": "2024y066d05h30m00s",
        "exper_nominal_stop": "2024y066d07h30m00s"
    },

    "ec064": {
        "mission_name": "mro_tgo_mex",
        "vex_file_name": "ec064.vex",
        "experiment_description": "MRO-TGO-MEX tracking",
        "exper_nominal_start": "2020y053d01h30m00s",
        "exper_nominal_stop": "2020y053d03h00m00s"
    },

    "v0314": {
        "mission_name": "vex",
        "vex_file_name": "v0314.vex",
        "experiment_description": "VEX tracking",
        "exper_nominal_start": "2014y073d08h30m00s",
        "exper_nominal_stop": "2014y073d11h29m00s"
    }

}

spacecraft_data = {
    "aka": {
        "mission_name": "akatsuki",
        "frequency_MHz": 8410.926,
        "antenna": "Ww",
        "snr": 0,
        "stochastic_noise": "0",
        "updated": "00.00.2010"
    },
    "bco": {
        "mission_name": "bepi colombo",
        "frequency_MHz": 8420.293,
        "antenna": "Cd",
        "snr": 6317,
        "stochastic_noise": "84 mHz",
        "updated": "17.02.2021"
    },
    "tgo": {
        "mission_name": "exomars - trace gas orbiter",
        "frequency_MHz": 8410.710,
        "antenna": "Ht",
        "snr": 11500,
        "stochastic_noise": "90 mHz",
        "updated": "29.07.2017"
    },
    "exo": {
        "mission_name": "exomars",
        "frequency_MHz": 8424.592,
        "antenna": "Xx",
        "snr": 0,
        "stochastic_noise": "0",
        "updated": "00.00.0000"
    },
    "gai": {
        "mission_name": "gaia",
        "frequency_MHz": "xxxx.xxx",
        "antenna": "Xx",
        "snr": 0,
        "stochastic_noise": "0",
        "updated": "00.00.0000"
    },

    "her": {
        "mission_name": "herschell",
        "frequency_MHz": "xxxx.xxx",
        "antenna": "Xx",
        "snr": 0,
        "stochastic_noise": "0",
        "updated": "00.00.0000"
    },
    "huy": {
        "mission_name": "huygens",
        "frequency_MHz": "xxxx.xxx",
        "antenna": "Xx",
        "snr": 0,
        "stochastic_noise": "0",
        "updated": "00.00.0000"
    },
    "ika": {
        "mission_name": "ikaros",
        "frequency_MHz": 8431.296,
        "antenna": "Ww",
        "snr": 0,
        "stochastic_noise": "0",
        "updated": "00.00.2010"
    },
    "jno": {
        "mission_name": "juno",
        "frequency_MHz": 8403.523,
        "antenna": "Ho",
        "snr": 450,
        "stochastic_noise": "150 mHz",
        "updated": "07.09.2020"
    },
    "mex": {
        "mission_name": "mars express",
        "frequency_MHz": 8420.750,
        "antenna": "Cd",
        "snr": 10000,
        "stochastic_noise": "30 mHz",
        "updated": "22.02.2020"
    },
    "min": {
        "mission_name": "mars insight",
        "frequency_MHz": 8404.502,
        "antenna": "T6",
        "snr": 120,
        "stochastic_noise": "300 mHz",
        "updated": "25.07.2020"
    },
    "mod": {
        "mission_name": "mars odyssey",
        "frequency_MHz": 8407.250,
        "antenna": "Cd",
        "snr": 360,
        "stochastic_noise": "830 mHz",
        "updated": "22.02.2020"
    },
    "mom": {
        "mission_name": "mars orbiter mission",
        "frequency_MHz": "xxxx.xxx",
        "antenna": "Xx",
        "snr": 0,
        "stochastic_noise": "0",
        "updated": "00.00.0000"
    },
    "perseverance": {
        "mission_name": "perseverance",
        "frequency_MHz": 8435.550,
        "antenna": "Cd",
        "snr": 115000,
        "stochastic_noise": "88 mHz",
        "updated": "18.02.2021"
    },
    "m20": {
        "mission_name": "mars2020",
        "frequency_MHz": 8415.000,
        "antenna": "??",
        "snr": 0,
        "stochastic_noise": "0",
        "updated": "00.00.0000"
    },
    "mro": {
        "mission_name": "mars reconnaissance orbiter",
        "frequency_MHz": 8439.750,
        "antenna": "Cd",
        "snr": 30,
        "stochastic_noise": "23.70 Hz",
        "updated": "22.02.2020"
    },
    "mvn": {
        "mission_name": "maven",
        "frequency_MHz": 8446.235,
        "antenna": "Xx",
        "snr": 0,
        "stochastic_noise": "0",
        "updated": "00.00.0000"
    },
    "curiosity": {
        "mission_name": "curiosity",
        "frequency_MHz": 8402.777,
        "antenna": None,
        "snr": None,
        "stochastic_noise": None,
        "updated": None
    },
    "hope": {
        "mission_name": "hope",
        "frequency_MHz": 8401.419,
        "antenna": None,
        "snr": None,
        "stochastic_noise": None,
        "updated": None
    },
    "ras": {
        "mission_name": "radio astron",
        "frequency_MHz": 8399.700,
        "antenna": "Wz",
        "snr": 32000,
        "stochastic_noise": "3 mHz",
        "updated": "29.09.2012"
    },
    "ros": {
        "mission_name": "rosetta",
        "frequency_MHz": 8421.875,
        "antenna": "Mh",
        "snr": 898,
        "stochastic_noise": "40 mHz",
        "updated": "01.04.2016"
    },
    "sta": {
        "mission_name": "stereo ab",
        "frequency_MHz": 8446.230,
        "antenna": "Wz",
        "snr": 3000,
        "stochastic_noise": "200 mHz",
        "updated": "06.02.2012"
    },
    "tiw": {
        "mission_name": "tianwen-1",
        "frequency_MHz": 8429.938,
        "antenna": "Cd",
        "snr": None,
        "stochastic_noise": "1.0",
        "updated": "01.04.2021"
    },
    "uly": {
        "mission_name": "ulysses",
        "frequency_MHz": "xxxx.xxx",
        "antenna": "Xx",
        "snr": 0,
        "stochastic_noise": "0",
        "updated": "00.00.0000"
    },
    "vex": {
        "mission_name": "venus express",
        "frequency_MHz": 8418.100,
        "antenna": "Ht",
        "snr": 9500,
        "stochastic_noise": "26 mHz",
        "updated": "21.04.2011"
    },
    "juc": {
        "mission_name": "juice",
        "frequency_MHz": 8435.5, #mas said file says: 8435-8436, so i put 8435.5
        "antenna": None,
        "snr": None,
        "stochastic_noise": None,
        "updated": None
    }
}

def parse_vex_freq_block(file_path):
    """
    Parses the $FREQ block from a vex file and returns a dictionary of stations with their channels.

    Parameters:
        file_path (str): Path to the vex file.

    Returns:
        dict: A dictionary where each station has its sub-dictionary of channels.
    """
    # Dictionary to store the parsed data
    stations_dict = defaultdict(dict)

    # Read the file
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Flag to track the $FREQ block
    in_freq_block = False
    current_stations = []

    # Pattern to match chan_def lines. Two patterns are defined, since the vex file format slightly changed over the years.

    chan_def_pattern = re.compile(r"chan_def\s*=\s*:\s*(\d+(?:\.\d+)? MHz)\s*:\s*(\w+)\s*:\s*(\d+\.\d+ MHz)\s*:\s*(&CH\d+)\s*:\s*(&BBC\d+)\s*:\s*(&\w+);")
    chan_def_pattern_new = re.compile(r"chan_def\s*=\s*:\s*(\d+(?:\.\d+)? MHz)\s*:\s*(\w+)\s*:\s*(\d+\.\d+ MHz)\s*:\s*(&CH\d+)\s*:\s*(&BBC\d+)\s*:\s*(&\w+);\s*\*\s*(\w+)")
    for line in lines:
        line = line.strip()

        # Check for the start of the $FREQ block
        if line.startswith("$FREQ;"):
            in_freq_block = True
            continue

        # Check for the end of a definition
        if line.startswith("enddef;") and in_freq_block:
            current_stations = []
            continue

        # Skip lines outside the $FREQ block
        if not in_freq_block:
            continue

        # Extract stations from the "stations = ..." line
        if "stations =" in line:
            stations_part = line.split("stations =")[1].strip()
            current_stations = [station.strip() for station in stations_part.split(":")]
            current_stations[-1] = current_stations[-1].rstrip('\\') # to remove the slash for last station in the vex file
            continue

        else:
            # Handle a new frequency definition
            if line.startswith("def"):
                current_def = {"name": line.split()[1], "stations": [], "chan_defs": [], "sample_rate": None}
                continue

            # Extract station information from the line
            if "evn+global" in line:
                if ":" in line:
                    stations_part = line.split(":")[1].strip()
                    current_stations = [station.strip() for station in stations_part.split(",")]
                    continue

        # Parse chan_def lines
        try:
            match = chan_def_pattern.match(line)
        except:
            try:
                match = chan_def_pattern_new.match(line)
            except:
                print('It was not possible to retrieve any baseband frequency.')

        if match:
            frequency, polarization, bandwidth, channel, bbc, cal = match.groups()

            # Add the channel details to each current station
            for station in current_stations:
                stations_dict[station][channel] = {
                    "frequency": frequency,
                    "polarization": polarization,
                    "bandwidth": bandwidth,
                    "bbc": bbc,
                    "cal": cal,
                }

    return dict(stations_dict)

def get_vex_file_path(experiment_name):

    """
    Constructs the file path for a VEX file associated with a specific experiment.

    This function generates the complete path to the VEX file by using the provided
    experiment name. It looks up the experiment details from a predefined dictionary
    and combines the base folder, mission folder, and the VEX file name to create
    the full file path.

    Parameters:
    experiment_name (str): The name of the experiment for which the VEX file path is to be generated.

    Returns:
    str: The full file path to the corresponding VEX file.

    Example:
    vex_file_path = get_vex_file_path('experiment_1')
    """
    main_vex_folder = "vex_files"
    vex_file_mission_folder = experiments[experiment_name]['mission_name']
    vex_file_name = experiments[experiment_name]['vex_file_name']

    vex_file_path = os.path.join(main_vex_folder, vex_file_mission_folder, vex_file_name)
    return vex_file_path


def get_baseband_frequencies_file(mission, experiment_name):

    """
    Retrieves the baseband frequencies for all stations from a VEX file.

    This function parses the frequency and bandwidth information from a VEX file
    associated with a given experiment and extracts the baseband frequency for
    each station that falls within the specified frequency range based on the
    spacecraft's X-band frequency.

    Parameters:
    mission (str): The name of the mission (e.g., 'JUICE', 'MRO') used to extract spacecraft data.
    experiment_name (str): The name of the experiment to fetch the corresponding VEX file.

    Returns:
    tuple: A tuple containing:
        - baseband_frequency (dict): A dictionary mapping each station to its baseband frequency
          (in MHz), rounded to 1 decimal place.
        - mas_x_band (float): The X-band frequency (in MHz) for the spacecraft.

    Example:
    baseband_frequencies, x_band = get_baseband_frequencies_file('JUICE', 'experiment_1')
    """

    file_path = get_vex_file_path(experiment_name)
    result = parse_vex_freq_block(file_path)
    mas_x_band = spacecraft_data[mission]['frequency_MHz']

    baseband_frequency = {mission: {}}
    # Print parsed dictionary
    for station, channels in result.items():
        for channel in channels:
            freq_string = channels[channel]['frequency']
            bandwidth_string = channels[channel]['bandwidth']

            # Remove " MHz" and convert to float
            float_frequency = float(freq_string.replace(" MHz", ""))
            float_bandwidth = float(bandwidth_string.replace(" MHz", ""))
            if  float_frequency < mas_x_band < float_frequency + float_bandwidth :
                # create baseband_frequency mission dictionary with baseband frequency for each station
                baseband_frequency[mission][station] = np.round(float_frequency,1)

    # Write to file
    output_file = f"baseband_frequencies/{mission}/{experiment_name}_baseband_frequencies.txt"
    with open(output_file, 'w') as f:
        # Write the header
        f.write(f"# Mission: {mission}\n")
        f.write(f"# X-band Observable: {mas_x_band}\n")
        f.write(f"# VEX file name: {file_path}\n")
        f.write("# Station | Baseband Frequency (MHz)\n\n")

        # Write the station and baseband frequencies
        for station, freq in baseband_frequency[mission].items():
            f.write(f"{station}: {freq}\n")


def get_baseband_frequency(mission, experiment_name, fdets_file):

    """
    Extracts the baseband frequency for a given station from a VEX file.

    This function retrieves the frequency and bandwidth information for each
    station from the VEX file, parses the data, and checks if the frequency
    of a given station falls within the specified range. If the condition is met,
    the baseband frequency is returned for that station.

    Parameters:
    mission (str): The name of the mission (e.g., 'JUICE', 'MRO') used to extract spacecraft data.
    experiment_name (str): The name of the experiment to fetch the corresponding VEX file.
    fdets_file (str): The file path of the FDets file, which contains the station information.

    Returns:
    float: The baseband frequency in MHz for the given station if found and within the correct range.
           If not found, the function will print a message and return None.

    Example:
    baseband_frequency = get_baseband_frequency('JUICE', 'experiment_1', 'fdets_data_file.txt')
    """
    station_fdets = fdets_file.split('/')[-1].split('.')[4]
    file_path = get_vex_file_path(experiment_name)
    result = parse_vex_freq_block(file_path)
    mas_x_band = spacecraft_data[mission]['frequency_MHz']

    if station_fdets not in result.keys():
        print(f'Fdets station: {station_fdets} not found in the vex file $FREQ block stations:\n'
              f'{result.keys()}.\n')
        frequencies_file_name = f'baseband_frequencies/{mission}/{experiment_name}_baseband_frequencies.txt'
        print(f'Trying retrieval from baseband frequencies file: {frequencies_file_name}...')
        if os.path.exists(frequencies_file_name):
            with open(frequencies_file_name, 'r') as f:
                lines = f.readlines()
                for line in lines:

                    if station_fdets in line:
                        baseband_frequency = float(line.split()[1])

        else:
            print(f'Error in get_baseband_frequency: frequencies_file_name: {frequencies_file_name} does not exist.')

    else:
        for station, channels in result.items():
            for channel in channels:
                freq_string = channels[channel]['frequency']
                bandwidth_string = channels[channel]['bandwidth']

                # Remove " MHz" and convert to float
                float_frequency = float(freq_string.replace(" MHz", ""))
                float_bandwidth = float(bandwidth_string.replace(" MHz", ""))
                if  float_frequency < mas_x_band < float_frequency + float_bandwidth :
                    if station == station_fdets:
                        # set baseband_frequency for station
                        baseband_frequency = float_frequency

    return(baseband_frequency)


def extract_vex_block(vex_file_path, block_name):
    """
    Reads the VEX file, extracts the IF block and returns it as a string.
    """
    with open(vex_file_path, 'r') as file:
        vex_content = file.read()

    # Match the IF block (assuming it starts with "$IF" and ends before the next block)
    if_block_match = re.search(f"\${block_name}\s*;.*?(\$|$)", vex_content, re.DOTALL)

    if if_block_match:
        return if_block_match.group(0)
    else:
        return None

def extract_station_if_mapping(vex_block):
    """
    Extracts a dictionary mapping station names to IF names from the IF block.
    """
    # Initialize dictionary to hold station to IF mapping
    station_if_mapping = defaultdict(list)

    # Match each `def` block
    def_blocks = re.findall(r"def\s+([^\s;]+);\n\*\s*.*?stations\s*=\s*([\w:]+)(.*?enddef;)", if_block, re.DOTALL)

    for block in def_blocks:
        lo_name, stations, content = block
        # Extract IF definitions
        if_defs = re.findall(r"if_def\s*=\s*&([^:]+)", content)

        # Split station names and map to IF names
        station_list = stations.split(":")
        for station in station_list:
            station_if_mapping[station].extend(if_defs)

    # Convert defaultdict to a regular dictionary
    return {station: if_list for station, if_list in station_if_mapping.items()}

def extract_station_mapping(block, block_type):
    """
    Extracts a dictionary mapping station names to corresponding names from the given block.
    Supports both 'IF' and 'FREQ' blocks.
    """
    # Initialize dictionary to hold station to mapping
    station_mapping = defaultdict(list)

    if block_type == 'IF':
        stations_dict = dict()

        # Flag to track the IF block
        in_if_block = False
        current_stations = []

        # Regex pattern to match the chan_def lines
        if_def_pattern = re.compile(r"if_def\s*=\s*(&IF_\w+)\s*:\s*\w+\s*:\s*[A-Za-z]+\s*:\s*[\d\.]+ MHz\s*:\s*[A-Za-z]+\s*;\s*\* PCall off!.*\s*[\d\.]+\s*[\d\.]+\s*\d+cm\s*\d+\s*")
        for line in block.splitlines():
            line = line.strip()

            # Check for the start of the $FREQ block
            if line.startswith("$IF;"):
                in_if_block = True
                continue

            # Check for the end of a definition
            if line.startswith("enddef;") and in_if_block:
                current_stations = []
                continue

            # Skip lines outside the $FREQ block
            if not in_if_block:
                continue

            # Extract stations from the "stations = ..." line
            if "stations =" in line:
                stations_part = line.split("stations =")[1].strip()
                current_stations = [station.strip() for station in stations_part.split(":")]
                current_stations[-1] = current_stations[-1].rstrip('\\') # to remove the slash for last station in the vex file
                continue
            else:
                # Handle a new frequency definition
                if line.startswith("def"):
                    current_def = {"name": line.split()[1], "stations": [], "chan_defs": [], "sample_rate": None}
                    continue

                # Extract station information from the line
                if "evn+global" in line:
                    if ":" in line:
                        stations_part = line.split(":")[1].strip()
                        current_stations = [station.strip() for station in stations_part.split(",")]
                        continue

            # Parse chan_def lines
            try:
                match = if_def_pattern.match(line)
            except:
                try:
                    match = if_def_pattern.match(line)
                except:
                    print('It was not possible to retrieve any baseband frequency.')

            if match:
                IF_name = match.groups()[0]

                # Add the channel details to each current station
                for station in current_stations:
                    stations_dict[station] = {
                        "IF_name": IF_name,
                    }

        return stations_dict


    elif block_type == 'FREQ':
        # Dictionary to store the parsed data
        stations_dict = dict()

        # Flag to track the $FREQ block
        in_freq_block = False
        current_stations = []

        # Regex pattern to match the chan_def lines
        chan_def_pattern = re.compile(r"chan_def\s*=\s*:\s*(\d+(?:\.\d+)? MHz)\s*:\s*(\w+)\s*:\s*(\d+\.\d+ MHz)\s*:\s*(&CH\d+)\s*:\s*(&BBC\d+)\s*:\s*(&\w+);")
        chan_def_pattern_new = re.compile(r"chan_def\s*=\s*:\s*(\d+(?:\.\d+)? MHz)\s*:\s*(\w+)\s*:\s*(\d+\.\d+ MHz)\s*:\s*(&CH\d+)\s*:\s*(&BBC\d+)\s*:\s*(&\w+);\s*\*\s*(\w+)")
        for line in block.splitlines():
            line = line.strip()

            # Check for the start of the $FREQ block
            if line.startswith("$FREQ;"):
                in_freq_block = True
                continue

            # Check for the end of a definition
            if line.startswith("enddef;") and in_freq_block:
                current_stations = []
                continue

            # Skip lines outside the $FREQ block
            if not in_freq_block:
                continue

            # Extract stations from the "stations = ..." line
            if "stations =" in line:
                stations_part = line.split("stations =")[1].strip()
                current_stations = [station.strip() for station in stations_part.split(":")]
                current_stations[-1] = current_stations[-1].rstrip('\\') # to remove the slash for last station in the vex file
                continue

            else:
                # Handle a new frequency definition
                if line.startswith("def"):
                    current_def = {"name": line.split()[1], "stations": [], "chan_defs": [], "sample_rate": None}
                    continue

                # Extract station information from the line
                if "evn+global" in line:
                    if ":" in line:
                        stations_part = line.split(":")[1].strip()
                        current_stations = [station.strip() for station in stations_part.split(",")]
                        continue

            # Parse chan_def lines
            try:
                match = chan_def_pattern.match(line)
            except:
                try:
                    match = chan_def_pattern_new.match(line)
                except:
                    print('It was not possible to retrieve any baseband frequency.')

            if match:
                frequency, polarization, bandwidth, channel, bbc, cal = match.groups()

                # Add the channel details to each current station
                for station in current_stations:
                    stations_dict.setdefault(station, {})[channel] = {
                        "frequency": frequency,
                        "polarization": polarization,
                        "bandwidth": bandwidth,
                        "bbc": bbc,
                        "cal": cal,
                    }

        return stations_dict

    elif block_type == 'BBC':
        # Dictionary to store the parsed data
        stations_dict = dict()

        # Flag to track the $FREQ block
        in_freq_block = False
        current_stations = []

        # Regex pattern to match the chan_def lines
        bbc_assign_pattern = re.compile(r"\*?\s*BBC_assign\s*=\s*(&BBC\d{2})\s*:\s*\d+\s*:\s*(&IF_\w+);")
        #bbc_assing_pattern_new = re.compile(r"chan_def\s*=\s*:\s*(\d+(?:\.\d+)? MHz)\s*:\s*(\w+)\s*:\s*(\d+\.\d+ MHz)\s*:\s*(&CH\d+)\s*:\s*(&BBC\d+)\s*:\s*(&\w+);\s*\*\s*(\w+)")
        for line in block.splitlines():
            line = line.strip()

            # Check for the start of the $FREQ block
            if line.startswith("$BBC;"):
                in_bbc_block = True
                continue

            # Check for the end of a definition
            if line.startswith("enddef;") and in_freq_block:
                current_stations = []
                continue

            # Skip lines outside the $FREQ block
            if not in_bbc_block:
                continue

            # Extract stations from the "stations = ..." line
            if "stations =" in line:
                stations_part = line.split("stations =")[1].strip()
                current_stations = [station.strip() for station in stations_part.split(":")]
                current_stations[-1] = current_stations[-1].rstrip('\\') # to remove the slash for last station in the vex file
                continue

            else:
                # Handle a new frequency definition
                if line.startswith("def"):
                    current_def = {"name": line.split()[1], "stations": [], "chan_defs": [], "sample_rate": None}
                    continue

                # Extract station information from the line
                if "evn+global" in line:
                    if ":" in line:
                        stations_part = line.split(":")[1].strip()
                        current_stations = [station.strip() for station in stations_part.split(",")]
                        continue

            # Parse chan_def lines
            try:
                match = bbc_assign_pattern.match(line)
            except:
                print('It was not possible to retrieve any baseband frequency.')

            if match:
                bbc_name, if_name = match.groups()[0], match.groups()[1]
                # Add the channel details to each current station
                for station in current_stations:
                    stations_dict[station] = {
                        "BBC_name": bbc_name,
                        "IF_name": if_name,
                    }

        return stations_dict

def write_missing_baseband_frequency(mission, experiment, station, missing_station):
    find_station_flag = False
    mission = mission.lower()
    experiment = experiment.lower()
    frequencies_file_name = f'baseband_frequencies/{mission}/{experiment}_baseband_frequencies.txt'
    if os.path.exists(frequencies_file_name):
        with open(frequencies_file_name, 'r') as f:
            lines = f.readlines()
            for line in lines:
                if station in line:
                    find_station_flag = True
                    missing_station_frequency = line.split()[1]
                    print('missing_station frequency', missing_station_frequency)
                    with open(frequencies_file_name, 'a') as a:
                        print('appending missing station')
                        a.write(f'{missing_station}: {missing_station_frequency} # assigned, missing from $FREQ block\n')
                        f.close()
                        break
                else:
                    continue

            if not find_station_flag:
                print(f'Could not assign any baseband frequency to missing station: {missing_station}')
        f.close()


    else:
        print(f'The frequencies file: {frequencies_file_name} does not exist.')


#print(parse_vex_freq_block('./vex_files/mex/gr035.vix').keys())
#print(get_baseband_frequencies_file('mex', 'gr035'))
#print(get_baseband_frequency('mex', 'gr035', './mex_dataset/Phobos/GR035/Fdets.mex2013.12.28.Nl.r2i.old'))
#print(extract_vex_block('./vex_files/mex/gr035.vix', 'BBC'))

def assign_missing_baseband_frequencies(mission_name, experiment_name):

    mission_name = mission_name.lower()
    experiment_name = experiment_name.lower()

    BBC_dict = extract_station_mapping(extract_vex_block('./vex_files/mex/gr035.vix', 'BBC'), 'BBC')
    IF_dict =  extract_station_mapping(extract_vex_block('./vex_files/mex/gr035.vix', 'IF'), 'IF')

    IF_keys, IF_values = IF_dict.keys(), IF_dict.values()
    FREQ_dict =  extract_station_mapping(extract_vex_block('./vex_files/mex/gr035.vix', 'FREQ'), 'FREQ')
    FREQ_keys, FREQ_values = FREQ_dict.keys(), FREQ_dict.values()

    IF_keys_not_in_FREQ_keys = [key_IF for key_IF in IF_keys if key_IF not in FREQ_keys]

    if len(IF_keys_not_in_FREQ_keys) >= 1:
        for IF_station, IF_inner_dict in IF_dict.items():  # Iterate over IF_dict to get station
            if IF_station in IF_keys_not_in_FREQ_keys:
                for IF_name in IF_inner_dict.values():
                    for BBC_station, BBC_inner_dict in BBC_dict.items():
                        if BBC_station in FREQ_keys:
                            if IF_name not in BBC_inner_dict.values():
                                missing_station_frequency_flag = False
                                continue
                            else:
                                # Printing the station name by referencing IF_station
                                print(f'IF_name: {IF_name} for station: {IF_station} was found in BBC dictionary.')
                                missing_station_frequency_flag = True
                                write_missing_baseband_frequency(mission_name, experiment_name, BBC_station, IF_station)
                                break

            else:
                continue

            if not missing_station_frequency_flag:
                print(f'Could not assign a baseband frequency to station: {IF_station} ')

def get_mission_from_experiment(experiment_name):
    experiment_name = experiment_name.lower()
    for experiment, values in experiments.items():
        if experiment == experiment_name:
            mission_name = values['mission_name']
            return mission_name
