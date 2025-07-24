"""
This script might be useful for splitting long scans for quicker and targeted analysis (for examples, mex gr035 are very long and dense).
To use it with other missions, just change the paths accordingly.

NOTE
A function performing the opposite operation (i.e. attaching multiple single scans, creating a complete scan) can also be found in
the pride_characterization_library.py, under the utilities class.
"""
from tudatpy.interface import spice
import os
import sys
sys.path.append('/Users/lgisolfi/ClionProjects/Allan_Features/Analysis_Scripts/') # Adjust this to your actual library location
from pride_characterization_library import PrideDopplerCharacterization

spice.load_standard_kernels()
pride = PrideDopplerCharacterization()
process_fdets = pride.ProcessFdets()
utilities = pride.Utilities()
analysis = pride.Analysis(process_fdets, utilities)

time_interval_minutes= 120
for fdets_file in os.listdir(f'/Users/lgisolfi/Desktop/data_archiving-1.0/dataset/mex/gr035/complete'):
    utilities.split_scan_by_time(
       input_folder=f'/Users/lgisolfi/Desktop/data_archiving-1.0/dataset/mex/gr035/complete',
        fdets_file = fdets_file,
        time_interval_minutes= time_interval_minutes,
        output_folder= f'/Users/lgisolfi/Desktop/data_archiving-1.0/dataset/mex/gr035/complete_{time_interval_minutes}_minutes'
    )
exit()
