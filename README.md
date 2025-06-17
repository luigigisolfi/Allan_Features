# Content
This repository contains:
- Conversion Scripts (like old2newfdets and correct_baseband_fdets) allowing for conversion between old format PRIDE data to the new (JUICE-like) format. The new format is: UTC TIME, SNR, SPECTRAL MAX, FREQUENCY DETECTION, DOPPLER NOISE
- Analysis Scripts (found in the Analysis_Script folder). These allow to perform a thorough analysis of the PRIDE dataset, and to facilitate tasks. These scripts heavily rely on the following library: pride_characterization_library.py . 
- Folders containing sample data in order to test the above-mentioned scripts.
- roaming scripts in the upper level folder, constiting of scripts that were used to convert the original PRIDE data into the converted old format files (new format). These scripts are only to be used in the case one wishes to reproduce the conversion from original PRIDE data to converted old format files (new format).*

Each script is commented so as to facilitate reproducibility. 
For any questions, please contact Luigi Gisolfi: l.gisolfi@tudelft.nl

## Useful (Provisional) Definitions
1) **riginal PRIDE data**
We refer to the PRIDE data retrieved from either JIVE's Marcopolo server, or files retrieved from University of Tasmania, or files from another source within the PRIDE team domain, as "original PRIDE data".
2) **new fdets format**
The new fdets format is (header): UTC TIME, SNR, SPECTRAL MAX, FREQUENCY DETECTION, DOPPLER NOISE
3) converted old format files
We refer to the converted old format files as the converted (into the new fdets format) original PRIDE data. The conversion operation was needed, since each mission had somehow different formats, and this would make parsing and processing non-sustainable. Conversion from original PRIDE data to converted old format files is achieved through both old2new_fdets.py and correct_baseband_fdets.py
4) **mission name**
This is actually the mission code. It sometimes correspond to the actual mission name, like for juice. Most times, though, it does not. For instance, min = mars insight and mex = mars express.
5) **single scan**
Single scans represent single passes recorded by a given PRIDE station.
6) **complete scan**
It is an union of single scans over a given coherent session. For instance, the union of all single scans recorded by Onsala during the GR035 experiment constitute a complete scan for Onsala.

*This is anyway not recommended, since the original PRIDE data contained various issues that were solved in itinere (for instance, some Wetzell Time Tags for Venus Express suffered from a time-bias of around 1000 days. This could only be spotted through post-processing visualization and we used the correct_timetags.py script to correct those badly formatted dates). 

## List of Scripts and brief explanation

We hereby provide a brief description of each script present in the repo, although more accurate docstrings and comments are present in each .py script.  

The following scripts are ready to be used for PRIDE Data Analysis. 

- **./Analysis_Scripts/experiment_statistics.py**
Script to perform statistical analysis on PRIDE Doppler FDETS files.

This script:
- Loads standard SPICE kernels.
- Initializes PRIDE characterization tools.
- Extracts Doppler noise, SNR, and elevation data from FDETS input files.
- Generates user-defined parameter plots (SNR, Doppler noise, raw FDETS).
- Generates elevation plots for each station and the whole experiment.
- Computes statistical summaries (mean, std) of SNR and Doppler noise.
- Combines SNR and elevation plots into single images for each station/date.

This specific demo is configured to analyze:
- Mission: 'juice'
- Experiment: 'ec094b'

**Note**:
This code has to be run before being able to run dataset_statistics.py .

- **./Analysis_Scripts/dataset_statistics.py**
This script can be used on the PRIDE_DATA folder to reproduce Vidhya's paper plots.
This script has the division by EXPERIMENT NAMES. If you want the division by times, check the script: dataset_statistics_by_time.
This script analyzes radio science experiments for: JUICE, MEX, MRO, MARS INSIGHT.
It computes and plots Key Performance Indicators (KPIs) such as mean Signal-to-Noise Ratio (SNR),
rms SNR, mean Doppler noise, and mean elevation angle across multiple PRIDE ground stations.

The workflow includes:
- Loading user-defined parameters (SNR, Doppler noise) and elevation data
- Filtering based on SNR and Doppler noise thresholds
- Computing weighted and unweighted averages
- Visualizing results in a set of subplots

**Requirements**:
    - The data files where the KPIs are read from are must be present in the output_dir folder, and they are created via experiment_statistics.py 

- **./Analysis_Scripts/pride_characterization_library**
This is the main library, used throughout the analysis scripts mentioned below.
It contains the following classes:

1) Pride_Doppler_Characterization.Process_Fdets
2) Pride_Doppler_Characterization.Utilities
3) Pride_Doppler_Characterization.Analysis
4) Pride_Doppler_Characterization.ProcessVexFiles
5) Pride_Doppler_Characterization.FormatFdets

- **./Analysis_Scripts/dataset_statistics_by_time.py**
This script analyzes radio science experiments computing and plotting Key Performance Indicators (KPIs) such as mean Signal-to-Noise Ratio (SNR),
rms SNR, mean Doppler noise, and mean elevation angle across multiple ground stations.

The workflow includes:
- Loading user-defined parameters (SNR, Doppler noise) and elevation data
- Filtering based on SNR and Doppler noise thresholds
- Computing weighted and unweighted averages
- Visualizing results in a set of subplots

**Requirements**:
    - The data files where the KPIs are read from are must be present in the output_dir folder, and they are created via experiment_statistics.py .

- **./Analysis_Scripts/experiments_allan_index.py**
This script allows for plotting of the Modified Allan Deviation per mission_name, given a date_folder, for all stations.


- **./Analysis_Scripts/split_scan_by_time**
This script might be useful for splitting long scans for quicker and targeted analysis (for examples, mex gr035 are very long and dense).
To use it with other missions, just change the paths accordingly.


**NOTE**
A function performing the opposite operation (i.e. attaching multiple single scans, creating a complete scan) can also be found in
the pride_characterization_library.py, under the utilities class.


The following scripts shouldnt be used, unless necessary (i.e. unless users only own old format, original PRIDE data):

- **./old2newfdets.py**
Script for processing spectra input files and optionally correcting baseband frequencies.
This script takes one or more spectra input files (or wildcard patterns to match files) and processes them to 
convert their formats or correct bad baseband frequencies based on the specified arguments.

- **./correct_baseband_fdets.py**
This script processes spectral data files, checks for proper baseband frequencies in headers, and applies corrections when necessary. It supports handling single files or wildcard patterns
for batch processing.

**Main functionalities**:
- Parses command-line arguments to accept input files and experiment name.
- Copies original files to backup with a ".old" extension.
- Reads and validates the number of columns in data files.
- Corrects baseband frequencies in file headers if missing or incorrect.
- Saves updated headers back to the original files.

- **./iterative_conversion_vex.py**
This code was used to iteratively convert all original PRIDE vex data to the new fdets format.

- **./temp_vex_folder_structure.py**
Move all subdirectories (named with 'vex_yymm' format) from the old base directory to a new base directory.

- **./temp_vex_renaming_folder.py**
This code was used on the old vex data folder structure (old structure can be found in PRIDE DATA folder). 
The script temp_vex_folder_structure is also applied to reach the final version of the vex folder structure. 
The new final structure for vex folders, which makes vex structure homogenous with all other missions, can be found in PRIDE DATA NEW.

- **./vex_experiment_statistics.py**
This script performs the same as what is performed in experiment_statistics, but it is tailored to the old structure of vex data.

- **./correct_timetags**
Experience with the original dataset showed that some stations might have suffered from a constant time offset in the time tag, while still keeping the correct header.
This script
1) checks for the observation time in the header
2) checks that the day in the time tag is the same as the header
3) corrects if needed.
