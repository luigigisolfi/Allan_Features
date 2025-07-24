# Content

This repository contains:

- **Conversion Scripts** (e.g., `old2newfdets.py` and `correct_baseband_fdets.py`)  
  These scripts convert original PRIDE data into the new standardized (JUICE-like) format:  
  `UTC TIME, SNR, SPECTRAL MAX, FREQUENCY DETECTION, DOPPLER NOISE`.

- **Analysis Scripts** (found in the `Analysis_Scripts/` folder)  
  These scripts enable detailed analysis of PRIDE datasets and facilitate various tasks. They heavily rely on the library `pride_characterization_library.py`.

- **Sample Data Folders**  
  For instance, `analysed_pride_data`, provided for testing the aforementioned scripts.

- **Vex and baseband frequencies folders**
  These folders contain vex files (.vix) for the following experiments:
  juice: ec094a,ec094b,
  mro: ec064,
  mex: gr035, m0325,m0327, m0403
  min: ed045a, ed045c, ed045d, ed045e, ed045f
  vex: v0314

  and the corresponding base frequencies (per station). 
  
  NOTE
  The `baseband_frequencies folder` is populated upon calling the `correct_baseband_frequencies.py` function.
  The files populating it are used to correct the original PRIDE data fdets for which a given baseband frequency was missing in the header, 
  and are never used within the Analysis Scripts (the user will virtually never interact with them, 
  unless they want to convert original PRIDE data into the new format).
  For more information, refer to `correct_baseband_frequencies.py` or to `old2new_fdets.py`

- **Roaming Scripts** (found in the root directory)  
  These were used to convert the original PRIDE data into the "converted old format" (i.e., the new standardized format).
  **Note:** These scripts are **not deprecated**, but users should exercise caution when using them. In many cases, **manual steps or corrections** may be necessary 
  during the conversion process (e.g., Wetzell time tag errors for Venus Express were corrected after analysis using `correct_timetags.py`). 
  They are primarily intended for **reproducibility** and **careful data handling**.

Each script is commented to facilitate reproducibility.  
For any questions, please contact Luigi Gisolfi: l.gisolfi@tudelft.nl

---
## Useful (Provisional) Definitions

1. **Original PRIDE data**  
   Refers to raw PRIDE data obtained from sources such as JIVE's Marcopolo server, the University of Tasmania, or any source internal to the PRIDE team.

2. **New fdets format**  
   A standardized header format used throughout this repo:  
   `UTC TIME, SNR, SPECTRAL MAX, FREQUENCY DETECTION, DOPPLER NOISE`

3. **Converted old format files**  
   Original PRIDE data converted into the new fdets format. This conversion was required due to inconsistencies in file formats across missions.  
   The conversion process uses `old2new_fdets.py` and `correct_baseband_fdets.py`.

4. **Mission name**  
   Refers to a mission code, which may or may not correspond to the actual mission name.  
   Example:
    - `juice` → JUICE
    - `min` → Mars InSight
    - `mex` → Mars Express

5. **Single scan**  
   A pass recorded by a specific PRIDE ground station.

6. **Complete scan**  
   The union of single scans from a coherent observation session.  
   Example: All Onsala station scans during the GR035 experiment make up a complete scan for Onsala.

---
## List of Scripts and Brief Explanation

Each script in this repository is well-documented via comments and docstrings.  
Below is a high-level summary of the available scripts.

### Analysis Scripts

- **`./Analysis_Scripts/pride_data_analysis.py`**  
  Performs statistical analysis on PRIDE Doppler FDETS files.  
  Features:
    - Loads standard SPICE kernels.
    - Initializes PRIDE characterization tools.
    - Extracts Doppler noise, SNR, and elevation.
    - Generates plots for user-defined parameters and station elevations.
    - Computes statistical summaries (mean, std) for each metric.
    - Combines plots per station and date
  Computes and plots:
    - Mean and RMS of SNR (FoMs as in V. Pallichadath's paper)
    - Mean Doppler noise
    - Mean elevation angle
    - Allan Deviation for all stations within a given experiment
f
    Supported missions: JUICE, MEX, MRO, MARS INSIGHT, VENUS EXPRESS
      **Requirements:**
  1) Fdets input files must be present in the root_dir, hereby assumed to be: `analysed_pride_data/`.
  2) The folder structure should be: `analysed_pride_data/<mission_code>/<mission_code>_YYMM/<mission_code>_YYMMDD`. 


- **`./Analysis_Scripts/pride_characterization_library.py`**  
  Core analysis library used across all scripts.  
  Contains:
    - `Pride_Doppler_Characterization.Process_Fdets`
    - `Pride_Doppler_Characterization.Utilities`
    - `Pride_Doppler_Characterization.Analysis`
    - `Pride_Doppler_Characterization.ProcessVexFiles`
    - `Pride_Doppler_Characterization.FormatFdets`

- **`./Analysis_Scripts/split_scan_by_time.py`**  
  Splits long scans for more focused analysis (e.g., GR035 from MEX).  
  To use it with other missions, update the file paths accordingly.

  **Note:**  
  The inverse operation (combining multiple scans into one) is available in the `Utilities` class of `pride_characterization_library.py`.

---

### Conversion & Legacy Scripts

These scripts should only be used when working with **original PRIDE data** in old formats:

- **`./old2newfdets.py`**  
  Converts raw spectral files to the new fdets format and optionally corrects baseband frequencies.  
  Features:
    - Accepts multiple files or wildcard patterns.
    - Applies format conversion and/or frequency corrections.

- **`./correct_baseband_fdets.py`**  
  Similar to above, this script:
    - Parses input files.
    - Validates data column consistency.
    - Backs up original files.
    - Fixes missing/incorrect baseband frequencies in headers.

- **`./iterative_conversion_vex.py`**  
  Batch converts original VEX-format data to the new fdets format.

- **`./temp_vex_folder_structure.py`**  
  Moves subdirectories (named like `vex_yymm`) from an old base path to a new one.

- **`./temp_vex_renaming_folder.py`**  
  Used to rename and restructure old VEX data folders to match the final format in `PRIDE DATA NEW`.

- **`./vex_experiment_statistics.py`**  
  Tailored variant of `experiment_statistics.py` for legacy VEX-format data.

- **`./correct_timetags.py`**  
  Corrects constant time offsets in the observation time tags.  
  Process:
    1. Reads observation time from headers.
    2. Verifies consistency with time tags.
    3. Applies correction if a mismatch is detected.
