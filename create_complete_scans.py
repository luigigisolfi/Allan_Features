from tudatpy.interface import spice
import os
import sys
sys.path.append('/Users/lgisolfi/ClionProjects/Allan_Features/Analysis_Scripts/') # Adjust this to your actual library location
from pride_characterization_library import PrideDopplerCharacterization
import glob

spice.load_standard_kernels()
pride = PrideDopplerCharacterization()
process_fdets = pride.ProcessFdets()
utilities = pride.Utilities()
analysis = pride.Analysis(process_fdets, utilities)

files = glob.glob('/Users/lgisolfi/Desktop/data_archiving-2.0/min/ed045c/input/single/Fdets.*.r2i.txt')
output_folder = '/Users/lgisolfi/Desktop/data_archiving-2.0/min/ed045c/input/complete/'
utilities.create_complete_scan_from_single_scans(files, output_folder)
