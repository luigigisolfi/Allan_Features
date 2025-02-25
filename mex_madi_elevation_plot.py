"""
Based on T. Bocanegra's paper for Mars Express, I tried to test the hypothesis of Doppler Frequency Noise Related to Elevation.
According to the plots in the paper, Ht station has the best SNR, while Ur has the worst SNR, and the reason lies in each antenna's elevation.
If this relation holds for the doppler noise as well, we would expect Ur and On to have a lower allan index than Ht, but this is not the case...

Maybe we can relax the condition on the slope and see what happens...
"""

from Allan_Slope_Library import Allan_Utility_Functions
import os
object = Allan_Utility_Functions()
object_fdets = object.ProcessFdets()


allan_folder = './mex_dataset/gr035/fdets/complete/allan_indexes'
os.makedirs(allan_folder, exist_ok=True)
extracted_data_ur = object_fdets.extract_parameters('./mex_dataset/gr035/fdets/complete/Fdets.mex2013.12.28.Ur.complete.r2i.txt_TATI', remove_outliers = True)
object_fdets.plot_parameters_error_bounds(extracted_data_ur, tau_min = 0, tau_max = 10, suppress = False)
ur_madi = object_fdets.get_allan_index(extracted_data_ur, tau_min = 0, tau_max = 10, save_dir = os.path.join(allan_folder, 'Ur'),suppress = False)
print(f'Ur Allan Index: {ur_madi}\n')
print('finished extracting Ur')
extracted_data_on = object_fdets.extract_parameters('./mex_dataset/gr035/fdets/complete/Fdets.mex2013.12.28.On.complete.r2i.txt_TATI', remove_outliers = True)
object_fdets.plot_parameters_error_bounds(extracted_data_on, suppress = False)
print('finished extracting onsala')
extracted_data_ht = object_fdets.extract_parameters('./mex_dataset/gr035/fdets/complete/Fdets.mex2013.12.28.Ht.complete.r2i.txt_TATI', remove_outliers = True)
object_fdets.plot_parameters_error_bounds(extracted_data_ht, suppress = False)
print('finished extracting ht')
on_madi = object_fdets.get_allan_index(extracted_data_on, tau_min = 0, tau_max = 10, save_dir = os.path.join(allan_folder, 'On'),suppress = False)
print(f'Onsala Allan Index: {on_madi}\n')
ht_madi = object_fdets.get_allan_index(extracted_data_ht, tau_min = 0, tau_max = 10, save_dir = os.path.join(allan_folder, 'ht'),suppress = False)
print(f'Ht Allan Index: {ht_madi}\n')
