#!/usr/bin/env python
# coding: utf-8

import matplotlib.pyplot as plt
from datetime import datetime
import argparse
import os
import numpy as np
import allantools  # Importing AllanTools for modified Allan deviation
from astropy.time import Time
from collections import Counter
from scipy.interpolate import CubicSpline
import sys
import regex as re
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
from collections import Counter
from astropy.time import Time
import allantools
from matplotlib import cm 
from datetime import datetime, timedelta
plt.rcParams.update({'figure.max_open_warning': 0})

########################################################################################################################################
##################################################### Main Class Definition ############################################################
########################################################################################################################################
class Allan_Utility_Functions():

    """
    This class allows fdets data characterization and processing for the following ESA missions:

    JUICE:   [experiments: EC094A, EC094B] + [data between Apr - May 2023]
    MEX:     [experiment: GR035 (Phobos Flyby)])
    INSIGHT: [experiments: ED045A,ED045C,ED045D,ED045E,ED045F]
    
    To be included: 
    VEX, 
    MRO, 
    what else...?
    
    """
    
    def __init__(self):
        self.result = 0
        
########################################################################################################################################
##################################################### ProcessFdets class ###############################################################
########################################################################################################################################
    
    # Function to extract data from the file
    class ProcessFdets:
        
        def __init__(self):
            self.result = 0

########################################################################################################################################
########################################################################################################################################
        
        def get_observation_date(self, filename):
            with open(filename, 'r') as file:        
                lines = file.readlines()
                header_match = re.search(r'Observation conducted on',lines[0])
                if not header_match:
                    print(f'Invalid File: {filename}. Reason: Invalid Observation Header: {header_match}. Skipping...\n')

                else:
            
                    date_match = re.search(r'Observation conducted on (\d{4}\.\d{2}\.\d{2})', lines[0])
                    if date_match:
                        # Return the date as a string in YYYY.MM.DD format
                        self.observation_date = date_match.group(1)
                        return (self.observation_date)
                    else:
                        print(f'Invalid File: {filename}. Reason: Invalid Observation Date Header. First Header Line is:\n{lines[0]}\n Trying Another Way...\n')

                    if not self.observation_date:
                        self.first_col_name = self.get_columns_names(filename)['first_col_name']
                        if self.first_col_name.strip() == 'UTC Time':
                                lines = file.readlines()
                                parts = lines[5].strip().split()
                                utc_time = parts[0]
                                utc_datetime = datetime.strptime(utc_time, "%Y-%m-%dT%H:%M:%S.%f")
                                self.observation_date = utc_datetime[0].strftime("%Y-%m-%d")
                                return(self.observation_date)
                            
                            
                        elif self.first_col_name.strip() == 'Modified Julian Date' or self.first_col_name.strip() == 'Modified JD':
                                lines = file.readlines()
                                parts = lines[5].strip().split()
                                mjd = float(parts[0])
                                self.observation_date = self.mjd_to_utc(mjd)
                                return(self.observation_date)
    
                        else:
                            print(f'Could Not Retrieve Observation Date from File: {filename}. Skipping...\n')

########################################################################################################################################
########################################################################################################################################
                            
        def get_base_frequency(self, filename):   
            observation_date_flag = self.get_observation_date(filename)
            if observation_date_flag:
                # Open the file and read lines
                with open(filename, 'r') as file:
                    lines = file.readlines()
                    
                    try:
                        if float(lines[1].split(' ')[3])*1e6 != 0.0: #the base frequency in the fdets is expressed as MHz
                            self.base_frequency = float(lines[1].split(' ')[3])*1e6 #the base frequency in the fdets is expressed as MHz
    
                        #BEWARE THIS ELIF!!!!!!!! We do not know exactly the InSight frequencies
                        elif float(lines[1].split(' ')[3]) == 0.0:
                            self.base_frequency =  8416.49*1e6 ##Handle InSight data, which has bad headers, 
                                                               ##but in the VEX file you can see the channel frequencies .  
                    except:
                        if lines[1].split(' ')[3] == '2xxx.xx':
                            self.base_frequency = 8412*1e6 #Handle MEX data, which has some bad headers. 
                                                            #Assuming all observations where at 8412 MHz (based on the good headers)
                return self.base_frequency                
            else: 
                print(f'Pointless to retrieve base frequency, as there is no valid observation date in the header.')    

########################################################################################################################################
########################################################################################################################################
                    
        def get_columns_names(self, filename):
            with open(filename, 'r') as file:   
                lines = file.readlines()

                # Extract column names from row 3 (index 2)
                columns_header = lines[2].strip()
     
                # Split the header line by '|'
                columns_header_parts = columns_header.split('|')

                if columns_header_parts[0] == None:
                    self.n_columns = 0

                if columns_header_parts[-1] == '':
                    self.n_columns = len(columns_header_parts) - 1 #some column headers end with |
                else:
                    self.n_columns = len(columns_header_parts)  # others dont ...           
                
                if self.n_columns == 5:  # Assign column names from the extracted parts
                    
                    try: #sometimes there is Format:, others there is not...
                        self.first_col_name = columns_header_parts[0].split(':',1)[1] #Typically UTC time (in YY-MM-DD for JUICE) 
                                                                                    # or Time(UTC) [s] for VEX
                    except: 
                        self.first_col_name = columns_header_parts[0]   # Typically Modified Julian Date or Modified JD 
                    self.second_col_name = columns_header_parts[1]  # Signal-to-Noise
                    self.third_col_name = columns_header_parts[2]  # Spectral Max
                    self.fourth_col_name = columns_header_parts[3]  # Freq. Detection
                    self.fifth_col_name = columns_header_parts[4]  # Doppler Noise
                    
                    return {
                        'number_of_columns': self.n_columns,
                        'first_col_name': self.first_col_name,
                        'second_col_name': self.second_col_name,
                        'third_col_name': self.third_col_name,
                        'fourth_col_name': self.fourth_col_name,
                        'fifth_col_name': self.fifth_col_name,
                    }
                
                elif self.n_columns == 6:
                    
                    try:  #sometimes there is Format:, others there is not...
                        self.first_col_name = columns_header_parts[0].split(':',1)[1]   # Typically Modified Julian Date or Modified JD 
                    except: 
                        self.first_col_name = columns_header_parts[0]   # Typically Modified Julian Date or Modified JD 

                    try:
                        self.second_col_name = columns_header_parts[1].split(':',1)[1]  # Typically Time(UTC) [s]
                    except: 
                        self.second_col_name = columns_header_parts[1] # Typically Time(UTC) [s]

                    self.third_col_name = columns_header_parts[2]  # Signal-to-Noise
                    self.fourth_col_name = columns_header_parts[3]  # Spectral Max
                    self.fifth_col_name = columns_header_parts[4]  # Freq. Detection 
                    self.sixth_col_name = columns_header_parts[5]  # Doppler Noise  

                    return {
                        'number_of_columns': self.n_columns,
                        'first_col_name': self.first_col_name,
                        'second_col_name': self.second_col_name,
                        'third_col_name': self.third_col_name,
                        'fourth_col_name': self.fourth_col_name,
                        'fifth_col_name': self.fifth_col_name,
                        'sixth_col_name': self.sixth_col_name,
                    }

                else:
                    print(f'Invalid number of columns: {self.n_columns}. Skipping File: {filename}...\n')

########################################################################################################################################
########################################################################################################################################

        def mjd_to_utc(self,mjd):
            # Define the starting reference date for MJD, which is 17 November 1858
            mjd_reference = datetime(1858, 11, 17, 0, 0, 0)
            
            # Calculate the UTC date by adding the number of days in the MJD to the reference date
            utc_date = mjd_reference + timedelta(days=mjd)
            return (utc_date.date)

########################################################################################################################################
########################################################################################################################################

        def utc_to_mjd(self, utc_date_str):
            # Convert the UTC date string (in YYYY.MM.DD format) to a datetime object
            utc_date = datetime.strptime(utc_date_str, "%Y.%m.%d")
            
            # MJD starts at midnight on November 17, 1858 (Julian Date 2400000.5)
            mjd_start = datetime(1858, 11, 17)
            
            # Calculate the difference in days between the UTC date and the MJD start
            delta = utc_date - mjd_start
            
            # Return the MJD date
            return delta.days + (delta.seconds / 86400.0)

########################################################################################################################################
########################################################################################################################################
        
        def mjd_utc_seconds_to_utc(self,mjd, utc_seconds):

            """
            
            Utility function to convert MEX/VEX type fdets data in UTC time
            
            """
            
            # Convert MJD to JD
            jd = mjd + 2400000.5
            
            # JD to Gregorian date (subtract 2400000.5 to get the fractional day)
            jd_days = int(jd)  # integer part is the day
            jd_fraction = jd - jd_days  # fractional part is the time of the day

            # Reference date for Julian Day 0 in proleptic Gregorian calendar
            jd_start = datetime(1858, 11, 17, 0, 0, 0) 
            
            # Convert JD days to datetime
            utc_date = jd_start + timedelta(days=jd_days - 2400000.5)
            
            # Handle fractional part of the day (converts to time in hours, minutes, seconds)
            utc_time = utc_date + timedelta(days=jd_fraction)
            
            # Add the UTC seconds from the second column
            final_utc = utc_time + timedelta(seconds=utc_seconds)
            
            # Return the full UTC date and time
            return final_utc

  ################################################################################################################################################################################################################################################################################

        def format_observation_time(self, observation_date, time_in_seconds):
            # Convert self.observation_date from "%Y.%m.%d" to a datetime object
            observation_date = datetime.strptime(observation_date, "%Y.%m.%d")
            
            # Split seconds into integer seconds and fractional microseconds
            int_seconds, frac_seconds = divmod(time_in_seconds, 1)
            
            # Convert to hours, minutes, seconds
            hours, remainder = divmod(int(int_seconds), 3600)
            minutes, seconds = divmod(remainder, 60)
            
            # Convert fractional part to microseconds
            microseconds = round(frac_seconds * 1_000_000)
            
            # Create a timedelta object for the time of day
            time_delta = timedelta(hours=hours, minutes=minutes, seconds=seconds, microseconds=microseconds)
            
            # Combine the observation date and the time of day
            full_datetime = observation_date + time_delta
            # Return the formatted date-time in the desired format
            return datetime.strptime(str(full_datetime), "%Y-%m-%d %H:%M:%S")
            

########################################################################################################################################
########################################################################################################################################
            
        def extract_parameters(self, filename):

            """
            Function to extract relevant parameters from the fdets, ready to be used.

            Inputs: filename [required]
                    JUICE [optional] (if JUICE == True, JUICE format is assumed. if JUICE == False, other missions format is assumed)

            Output: Dictionary with relevant information
            
                    dict{'utc_datetime': self.utc_datetime,
                    'signal_to_noise': self.signal_to_noise,
                    'doppler_noise_hz': self.doppler_noise_hz,
                    'base_frequency': self.base_frequency,
                    'frequency_detection': self.frequency_detection,
                    'first_col_name': self.first_col_name,
                    'second_col_name': self.second_col_name,
                    'fifth_col_name': self.fifth_col_name,
                    'utc_date': self.utc_date}
                    
            """

            self.observation_date = self.get_observation_date(filename)

            if self.observation_date != None:
                utc_time = []
                self.signal_to_noise = []
                self.doppler_noise_hz = []
                self.base_frequency = self.get_base_frequency(filename)
                self.frequency_detection = []
                   
                # Assign column names from the extracted parts
                columns_dict =self.get_columns_names(filename)
                self.n_columns = columns_dict['number_of_columns']

                if self.n_columns == 5:
                    self.first_col_name = columns_dict['first_col_name'] # Should be UTC time (in YY-MM-DD for JUICE) or mjd
                    self.second_col_name = columns_dict['second_col_name']  # Signal-to-Noise
                    self.third_col_name = columns_dict['third_col_name']   # Doppler noise [Hz]
                    self.fourth_col_name = columns_dict['fourth_col_name']
                    self.fifth_col_name = columns_dict['fifth_col_name'] 


                    with open(filename, 'r') as file:   
                        lines = file.readlines()                
                        # Skip the first 4 lines (index 0 to 3)
                        for line in lines[4:]:
                            parts = line.strip().split()
                            # Extract the necessary columns

                            utc_time.append(parts[0])
                            self.signal_to_noise.append(float(parts[1]))
                            self.doppler_noise_hz.append(float(parts[4]))
                            self.frequency_detection.append(float(parts[3]))
                    
                            # Convert UTC time to datetime objects for plotting
                            if self.first_col_name.strip() == 'UTC Time':
                                self.utc_datetime = [datetime.strptime(t, "%Y-%m-%dT%H:%M:%S.%f") for t in utc_time]  
                            else: #else, it is Time(UTC) [s]
                                try:
                                    self.utc_datetime = [datetime.strptime(t, "%Y-%m-%dT%H:%M:%S.%f") for t in utc_time]
                                except:
                                    self.utc_datetime = [self.format_observation_time(self.observation_date, float(t)) for t in utc_time]  
                                
                    return {
                        'utc_datetime': self.utc_datetime,
                        'signal_to_noise': self.signal_to_noise,
                        'doppler_noise_hz': self.doppler_noise_hz,
                        'base_frequency': self.base_frequency,
                        'frequency_detection': self.frequency_detection,
                        'first_col_name': self.first_col_name,
                        'second_col_name': self.second_col_name,
                        'fifth_col_name': self.fifth_col_name,
                        'utc_date': self.observation_date
                    }

                else:
                    self.first_col_name = columns_dict['first_col_name'] # Should be UTC time (in YY-MM-DD for JUICE)
                    self.second_col_name = columns_dict['second_col_name']  # Signal-to-Noise
                    self.third_col_name = columns_dict['third_col_name']   # Doppler noise [Hz]
                    self.fourth_col_name = columns_dict['fourth_col_name']
                    self.fifth_col_name = columns_dict['fifth_col_name'] 
                    self.sixth_col_name = columns_dict['sixth_col_name']    
                    
                    mjd_day = self.utc_to_mjd(self.observation_date)
                    utc_seconds = []
                    with open(filename, 'r') as file:   
                        lines = file.readlines()  
                        # Skip the first 4 lines (index 0 to 3)
                        for line in lines[4:]:
                            parts = line.strip().split()
                            # Extract the necessary columns
                            utc_seconds.append(float(parts[1])) 
                            self.signal_to_noise.append(float(parts[2]))
                            self.doppler_noise_hz.append(float(parts[5]))
                            self.frequency_detection.append(float(parts[4]))
        
                        utc_dates = [self.mjd_utc_seconds_to_utc(mjd_day, utc_second) for utc_second in  utc_seconds]
                       
                        self.utc_datetime = [self.parse_datetime(t) for t in utc_dates]
                        
                    return {
                        'utc_datetime': self.utc_datetime,
                        'signal_to_noise': self.signal_to_noise,
                        'doppler_noise_hz': self.doppler_noise_hz,
                        'base_frequency': self.base_frequency,
                        'frequency_detection': self.frequency_detection,
                        'first_col_name': self.first_col_name,
                        'second_col_name': self.second_col_name,
                        'fifth_col_name': self.fifth_col_name,
                        'utc_date': self.observation_date
                    }                

        ########################################################################################################################################
########################################################################################################################################

        def parse_datetime(self,t):

            """
            
            Utility function trying to parse different time formats, 
            as different missions/experiments might use different formats
            
            """
            try:
                return datetime.strptime(str(t), "%Y-%m-%d %H:%M:%S.%f")
            except ValueError:
                return datetime.strptime(str(t), "%Y-%m-%d %H:%M:%S")
            except ValueError:
                return datetime.strptime(str(t), "%Y-%m-%dT%H:%M:%S.%f")
            except ValueError:
                return datetime.strptime(str(t), "%Y-%m-%dT%H:%M:%S")

        ########################################################################################################################################
########################################################################################################################################

        def plot_parameters(self, extracted_data, save_dir=None, suppress=False): 

            """ 
            This function plots thos regions in the mADEV of our data for which the slope is closest
            to the white noise (target_slope = -1).

            orange regions: abs(slope - target_slope) <= 0.01
            green regions: abs(slope - target_slope) <= 0.01

            A region is created only if two or more consecutive data points satisfy the orange/green slope condition. 

            Input: object.ProcessFdets().extract_parameters [required]
                   save_dir (to save the plot) [optional, default is None]
                   suppress (to suppress the plot show) [optional, default is False]
                   
            Output: mADEV & SNR plots
            
            NOTES: Please note that, within this function, the slope is computed between data points only 
                   (differently from the Get_Slope_At_Tau function, where a cubic spline interpolation is used,
                   and hence one can retrieve slopes at any time). 
        
            """
            target_slope = -1
            if extracted_data != None:
                
                # Extract data
                utc_datetime = extracted_data['utc_datetime']
                signal_to_noise = extracted_data['signal_to_noise']
                doppler_noise_hz = extracted_data['doppler_noise_hz']
                first_col_name = extracted_data['first_col_name']
                second_col_name = extracted_data['second_col_name']
                fifth_col_name = extracted_data['fifth_col_name']
                utc_date = extracted_data['utc_date']
                base_frequency = extracted_data['base_frequency']
                frequency_detection = extracted_data['frequency_detection']
                
                # Prepare the x-ticks to change on minute changes
                xticks = []
                xtick_labels = []
                last_minute = None
                for i, dt in enumerate(utc_datetime):
                    minute = dt.strftime("%H:%M")
                    if minute != last_minute:
                        xticks.append(i)
                        xtick_labels.append(dt.strftime("%H:%M:%S"))
                        last_minute = minute
            
                # Set up the 3x1 subplot figure
                fig, axs = plt.subplots(3, 1, figsize=(10, 15))
                fig.subplots_adjust(hspace=0.5)
    
                # Plot Signal-to-Noise vs UTC time
                axs[1].plot(range(len(utc_datetime)), signal_to_noise, label=f'{second_col_name}', marker='o', linestyle='-', color='blue')
                axs[1].set_xlabel(f'UTC Time (HH:MM:SS) on {utc_date}')
                axs[1].set_ylabel(second_col_name)
                axs[1].set_title(f'{second_col_name} vs {first_col_name}')
                axs[1].set_xticks(xticks)
                axs[1].set_xticklabels(xtick_labels, rotation=45)
                axs[1].legend()
                axs[1].grid(True)
            
                # Plot Doppler noise [Hz] vs UTC time
                axs[2].plot(range(len(utc_datetime)), doppler_noise_hz, label=f'{fifth_col_name}', marker='o', linestyle='-', color='orange')
                axs[2].set_xlabel(f'UTC Time (HH:MM:SS) on {utc_date}')
                axs[2].set_ylabel(fifth_col_name)
                axs[2].set_title(f'{fifth_col_name} vs {first_col_name}')
                axs[2].set_xticks(xticks)
                axs[2].set_xticklabels(xtick_labels, rotation=45)
                axs[2].legend()
                axs[2].grid(True)
                    
                # Calculate sampling rate in Hz
                t_jd = [Time(time).jd for time in utc_datetime]
                # Calculate the differences
                diffs = np.diff(t_jd)
                
                # Get the most common difference
                most_common_diff = Counter(diffs).most_common(1)
                if most_common_diff and most_common_diff[0][0] != 0:
                    rate_fdets = 1 / (most_common_diff[0][0] * 86400)
                else:
                    print("Most common difference is zero or not found; cannot calculate rate_fdets. No plot available.")
                    return(None)
                
                # Proceed to calculate modified Allan deviation, ensuring no invalid values
                try:
                    taus_doppler, mdev_doppler, errors, ns = allantools.mdev(
                        np.array(doppler_noise_hz) / (np.array(frequency_detection) + base_frequency),
                        rate=rate_fdets,
                        data_type='freq',
                        taus='all'
                    )
                except Exception as e:
                    print(f"An error occurred: {e}")
                # Calculate Modified Allan Deviation for Doppler noise
                taus_doppler, mdev_doppler, errors, ns = allantools.mdev(
                    np.array(doppler_noise_hz) / (np.array(frequency_detection) + base_frequency),
                    rate=rate_fdets,
                    data_type='freq',
                    taus='all'
                )
    
                axs[0].errorbar(taus_doppler, mdev_doppler, yerr=errors)
                # Assuming you have your RMS values in a list or array
                rms_values = errors # Replace this with your actual RMS values
                weights = 1 / (np.array(rms_values) ** 2)
                #weights = np.log(weights + 1e-6)  # Adding a small value to avoid log(0)
                
                # Normalize the weights to a range between 0 and 1
                norm_weights = (weights - np.max(weights)) / (np.min(weights) - np.max(weights))
                cmap = cm.get_cmap('plasma')  # You can choose any colormap
    
                # Generate the white noise reference line
                mdev_white = [mdev_doppler[0]]  # Initialize with the first value
                for i in range(1, len(taus_doppler)):
                    mdev_white.append(mdev_doppler[0] * (taus_doppler[i] / taus_doppler[0])**(target_slope))
            
                # Plot the Modified Allan Deviation and the white noise line
                axs[0].loglog(taus_doppler, mdev_white, linestyle='--', color='black', label=r'White Noise $\propto{\tau^{-0.5}}$')
                axs[0].loglog(taus_doppler, mdev_doppler, marker='o', linestyle='dashed', color='b', label='Doppler Noise')
                axs[0].axhline(mdev_doppler[0], color='fuchsia', linestyle='dashdot', label=r'Pink Noise $\propto{\tau^0}$')
            
                # Set labels and grid
                axs[0].set_xlabel('Averaging Time (s)')
                axs[0].set_ylabel('Modified Allan Deviation')
                axs[0].set_title('Modified Allan Deviation Plot')
                axs[0].grid(True, which="both", ls="--")
            
                ### Slope Calculation and Region Coloring Logic ###
                
                # Calculate slopes for given intervals
                slopes = []
                slope_errors = []
                target_slope = -1
                
                # Detect consecutive regions where the slope is close to -0.5 (white noise)
                threshold_10 = 0.1 # 10% away from white noise
                threshold_1 = 0.01   # 1% away from white noise
                regions_to_color_10 = []  # Regions for 10% threshold
                regions_to_color_1 = []  # Regions for 1% threshold
                start_idx_10 = None
                start_idx_1 = None
    
                for i in range(len(taus_doppler) - 1):
                    slope = (np.log10(mdev_doppler[i + 1]) - np.log10(mdev_doppler[i])) / (np.log10(taus_doppler[i + 1]) - np.log10(taus_doppler[i]))
                    slopes.append(slope)           
                for i, slope in enumerate(slopes):
                    #print(f"Interval {taus_doppler[i]:.3f} - {taus_doppler[i+1]:.3f}: abs(Slope - target_slope)= {np.abs(slope - target_slope)}")
    
                    if np.abs(slope - target_slope) < threshold_10:
                        if np.abs(slope - target_slope) < threshold_1:
                            if start_idx_1 is None:
                                start_idx_1 = i  # Start of a region within 1%
                        else:
                            if start_idx_1 is not None:
                                regions_to_color_1.append((start_idx_1, i))  # End of 1% region
                                start_idx_1 = None
            
                        if start_idx_10 is None:
                            start_idx_10 = i  # Start of a region within 10%
                    else:
                        if start_idx_10 is not None:
                            regions_to_color_10.append((start_idx_10, i))  # End of 10% region
                            start_idx_10 = None
                        if start_idx_1 is not None:
                            regions_to_color_1.append((start_idx_1, i))  # End of 1% region
                            start_idx_1 = None
                
                if start_idx_10 is not None:
                    regions_to_color_10.append((start_idx_10, len(slopes)))
                if start_idx_1 is not None:
                    regions_to_color_1.append((start_idx_1, len(slopes)))
            
                # Initialize flags to track if labels have been added
                added_legend_1_percent = False
                added_legend_10_percent = False
                
                # Keep track of the regions where 1% has been satisfied
                regions_covered_by_1_percent = []
                
                # Plot the colored regions for 1% threshold
                for start, end in regions_to_color_1:
                    if start < len(taus_doppler) and end <= len(taus_doppler):  # Check bounds
                        axs[0].axvspan(taus_doppler[start], taus_doppler[end], color='orange', alpha=0.5, label='1% Threshold' if not added_legend_1_percent else "")
                        added_legend_1_percent = True  # Set flag to True after first label
                        regions_covered_by_1_percent.append((start, end))  # Add this region to the covered list
                
                # Plot the colored regions for 10% threshold, only if they are not covered by the 1% regions
                for start, end in regions_to_color_10:
                    # Check if this region overlaps with any 1% region
                    overlap_with_1_percent = False
                    for start_1, end_1 in regions_covered_by_1_percent:
                        if not (end < start_1 or start > end_1):  # Check for any overlap
                            overlap_with_1_percent = True
                            break
                        
                    # If no overlap with 1% regions, plot the 10% region
                    if not overlap_with_1_percent and start < len(taus_doppler) and end <= len(taus_doppler):
                        axs[0].axvspan(taus_doppler[start], taus_doppler[end], color='green', alpha=0.6, label='10% Threshold' if not added_legend_10_percent else "")
                        added_legend_10_percent = True  # Set flag to True after first label
            
                # Tick the x-axis values at the center of each region
                for start, end in regions_to_color_1 + regions_to_color_10:
                    center = (taus_doppler[start] + taus_doppler[end]) / 2
                    axs[0].axvline(center, color='red', linestyle=':', linewidth=1)
        
                    # Prepare valid xticks based on the length of taus_doppler
                    log_ticks = np.logspace(np.log10(taus_doppler[0]), np.log10(taus_doppler[-1]), num=10)
                    # Find the corresponding indices for these ticks in taus_doppler
                    valid_xticks = np.searchsorted(taus_doppler, log_ticks)
                    
                    # Filter valid_xticks to ensure they are within bounds
                    valid_xticks = valid_xticks[valid_xticks < len(taus_doppler)]
                    # Add ticks at the center of each region on the log-log scale
                    center_ticks = []
                    for start, end in regions_to_color_1 + regions_to_color_10:
                        center = (taus_doppler[start] + taus_doppler[end]) / 2
                        center_ticks.append(center)
                    
                    all_ticks = np.unique(np.concatenate([taus_doppler[valid_xticks], center_ticks]))
                    
                    # Set x-ticks on log scale
                    axs[0].set_xticks(all_ticks)
                    axs[0].set_xticklabels([f"{int(tick)}" for tick in all_ticks])
                    axs[0].legend()
                    axs[0].tick_params(axis='x', which='major', labelsize=8, rotation = 45) 
    
                # Save the figure if a directory is specified
                if save_dir:
                    if os.path.exists(os.path.join(save_dir, 'parameters_plot.png')):
                        print(f'Removing {os.path.join(save_dir, "parameters_plot.png")}')
                        os.remove(os.path.join(save_dir, 'parameters_plot.png'))
                    plt.savefig(os.path.join(save_dir, 'parameters_plot.png'))
                    print(f'Saving {os.path.join(save_dir, "parameters_plot.png")}')
                    plt.close(fig)
                
                if suppress == True:
                    plt.close(fig)  # Close the figure to free memory
                else:
                    plt.show()
                    plt.close(fig)

        ########################################################################################################################################
########################################################################################################################################
                
        def plot_parameters_error_bounds(self, extracted_data, save_dir=None, suppress=False): 

            """ 
            This function is supposed to be an improvement of plot_parameters, as 
            it accounts for the errorbars in the slope computation. The "acceptable regions"
            are either:
            
            1) those regions for which, taking the error bars into account, satisfy 

            -0.5 belongs to [slope_min, slope_max]
            
            with 
            
            slope_min = [(slope[i+1] - error_plus[i+1]) - (slope[i] + error_minus[i])]/(tau[i+1] -tau[i]) 
            slope_max = [(slope[i+1] + error_plus[i+1]) - (slope[i] - error_minus[i])]/(tau[i+1] -tau[i]) 

            
            2) those regions for which the slope at the data point satisfies: abs(slope - target_slope) <= 0.1

            I think it is good to have this second condition, as, for small taus, the error bars are small. 

            Input: object.ProcessFdets().extract_parameters [required]
                   save_dir (to save the plot) [optional, default is None]
                   suppress (to suppress the plot show) [optional, default is False]
                   
            Output: mADEV & SNR plots

            NOTES: Please note that, within this function, the slope is computed between data points only 
                   (differently from the Get_Slope_At_Tau function, where a cubic spline interpolation is used,
                   and hence one can retrieve slopes at any time). 
                   The reason for this is that the error bars (which are used to compute 
                   the minimum and maximum slope (two-sided triangle) are only available for the data points, 
                   so it would be meaningless to interpolate the data, 
                   as we do not actually know how to interpolate the error bars.
        
            """
            target_slope = -1
            if extracted_data != None:
                # Extract data
                utc_datetime = extracted_data['utc_datetime']
                signal_to_noise = extracted_data['signal_to_noise']
                doppler_noise_hz = extracted_data['doppler_noise_hz']
                first_col_name = extracted_data['first_col_name']
                second_col_name = extracted_data['second_col_name']
                fifth_col_name = extracted_data['fifth_col_name']
                utc_date = extracted_data['utc_date']
                base_frequency = extracted_data['base_frequency']
                frequency_detection = extracted_data['frequency_detection']
                
                # Prepare the x-ticks to change on minute changes
                xticks = []
                xtick_labels = []
                last_minute = None
                for i, dt in enumerate(utc_datetime):
                    minute = dt.strftime("%H:%M")
                    if minute != last_minute:
                        xticks.append(i)
                        xtick_labels.append(dt.strftime("%H:%M:%S"))
                        last_minute = minute
            
                # Set up the 3x1 subplot figure
                fig, axs = plt.subplots(3, 1, figsize=(10, 15))
                fig.subplots_adjust(hspace=0.5)
            
                # Plot Signal-to-Noise vs UTC time
                axs[1].plot(range(len(utc_datetime)), signal_to_noise, label=f'{second_col_name}', marker='o', linestyle='-', color='blue')
                axs[1].set_xlabel(f'UTC Time (HH:MM:SS) on {utc_date}')
                axs[1].set_ylabel(second_col_name)
                axs[1].set_title(f'{second_col_name} vs {first_col_name}')
                axs[1].set_xticks(xticks)
                axs[1].set_xticklabels(xtick_labels, rotation=45)
                axs[1].legend()
                axs[1].grid(True)
            
                # Plot Doppler noise [Hz] vs UTC time
                axs[2].plot(range(len(utc_datetime)), doppler_noise_hz, label=f'{fifth_col_name}', marker='o', linestyle='-', color='orange')
                axs[2].set_xlabel(f'UTC Time (HH:MM:SS) on {utc_date}')
                axs[2].set_ylabel(fifth_col_name)
                axs[2].set_title(f'{fifth_col_name} vs {first_col_name}')
                axs[2].set_xticks(xticks)
                axs[2].set_xticklabels(xtick_labels, rotation=45)
                axs[2].legend()
                axs[2].grid(True)
                    
                # Calculate sampling rate in Hz
                t_jd = [Time(time).jd for time in utc_datetime]
                # Calculate the differences
                diffs = np.diff(t_jd)
                
                # Get the most common difference
                most_common_diff = Counter(diffs).most_common(1)
                if most_common_diff and most_common_diff[0][0] != 0:
                    rate_fdets = 1 / (most_common_diff[0][0] * 86400)
                else:
                    print("Most common difference is zero or not found; cannot calculate rate_fdets. Hence, no plot is available.")
                    return(None)
                
                # Proceed to calculate modified Allan deviation, ensuring no invalid values
                try:
                    taus_doppler, mdev_doppler, errors, ns = allantools.mdev(
                        np.array(doppler_noise_hz) / (np.array(frequency_detection) + base_frequency),
                        rate=rate_fdets,
                        data_type='freq',
                        taus='all'
                    )
                except Exception as e:
                    print(f"An error occurred: {e}")
    
                # Defining Weights
                rms_values = errors 
                weights = 1 / (np.array(rms_values) ** 2)
    
                # Normalize the weights to a range between 0 and 1
                norm_weights = (weights - np.max(weights)) / (np.min(weights) - np.max(weights))
                cmap = cm.get_cmap('plasma')  # You can choose any colormap
                axs[0].errorbar(taus_doppler, mdev_doppler, yerr=errors)
    
                # Generate the white noise reference line
                mdev_white = [mdev_doppler[0]]  # Initialize with the first value
                for i in range(1, len(taus_doppler)):
                    mdev_white.append(mdev_doppler[0] * (taus_doppler[i] / taus_doppler[0])**(target_slope))
                
                # Plot the Modified Allan Deviation and the white noise line
                axs[0].loglog(taus_doppler, mdev_white, linestyle='--', color='black', label=r'White Noise $\propto{\tau^{-0.5}}$')
                axs[0].loglog(taus_doppler, mdev_doppler, marker='o', linestyle='dashed', color='b', label='Doppler Noise')
                axs[0].axhline(mdev_doppler[0], color='fuchsia', linestyle='dashdot', label=r'Pink Noise $\propto{\tau^0}$')
                
                # Set labels and grid
                axs[0].set_xlabel('Averaging Time (s)')
                axs[0].set_ylabel('Modified Allan Deviation')
                axs[0].set_title('Modified Allan Deviation Plot')
                axs[0].grid(True, which="both", ls="--")
                axs[0].tick_params(axis='x', which='major', labelsize=8, rotation = 45)            
                ### Slope Calculation and Region Coloring Logic ###
                
                # Calculate slopes and determine is_close
                slopes = []
                is_close = []
                mean_weights = []
                target_slope = -1
                threshold = 0.01  # Define a threshold for closeness
                
                for i in range(len(taus_doppler) - 1):
                    slope = (np.log10(mdev_doppler[i + 1]) - np.log10(mdev_doppler[i])) / (np.log10(taus_doppler[i + 1]) - np.log10(taus_doppler[i]))
                    slopes.append(slope)
                
                    # here, for slope_error_minus, we are creating a "double triangle" connecting 
                    # the lowest value at i with the highest at i+1 and checking the slope.
                    # same for slope_error_plus, but with interchanged endpoints
                    
                    slope_error_plus = ((np.log10(mdev_doppler[i + 1] - errors[i + 1]) - np.log10(mdev_doppler[i] + errors[i])) /
                                       (np.log10(taus_doppler[i + 1]) - np.log10(taus_doppler[i])))
                    slope_error_minus = ((np.log10(mdev_doppler[i + 1] + errors[i + 1]) - np.log10(mdev_doppler[i] - errors[i])) /
                                       (np.log10(taus_doppler[i + 1]) - np.log10(taus_doppler[i])))
                    
                    # Check if the slope is close to the target slope
                    is_close.append((slope_error_plus <= target_slope and target_slope <= slope_error_minus) or (slope_error_minus <= target_slope and target_slope <= slope_error_plus) or abs(slope - target_slope) <= 0.1)
                    #print(f"Interval {taus_doppler[i]:.3f} - {taus_doppler[i+1]:.3f}: Slope Bounds= {slope_error_plus:.3f},{slope_error_minus:.3f} Target Slope Falls Within Slope Bounds, or the slope is closer than 0.1: {'Yes' if is_close[i] else 'No'}")
                    mean_weights.append((norm_weights[i] + norm_weights[i+1])/2) #the mean of the two weights at endpoints 
                                                                                # is considered for the interval.
                    
                # Plot the regions where is_close is True
                start_idx = None
                for i in range(len(is_close)):
                    if is_close[i]:
                        if start_idx is None:
                            start_idx = i  # Start of a region
    
                        
                        # For each interval in the region, use the specific weight of the interval
                        weight = mean_weights[i]  # Use weight for the current interval
                        color = cmap(weight)  # Get color based on current interval's weight
                        
                        # Print the weight and corresponding color for debugging
                        #print(f"Interval {taus_doppler[i]:.3f} - {taus_doppler[i+1]:.3f}, Weight: {weight:.3f}, Color: {color}")
                        
                        # Plot the interval using axvspan
                        axs[0].axvspan(taus_doppler[i], taus_doppler[i + 1], color=color, alpha=0.1,
                                       label='Acceptable' if 'Acceptable' not in axs[0].get_legend_handles_labels()[1] else "")
                    
                    else:
                        if start_idx is not None:
                            start_idx = None  # Reset start_idx when region ends
                
                # Check for any remaining open region at the end
                if start_idx is not None:
                    for i in range(start_idx, len(is_close) - 1):
                        if is_close[i]:
                            weight = mean_weights[i]
                            color = cmap(weight)
                            axs[0].axvspan(taus_doppler[i], taus_doppler[i + 1], color=color, alpha=0.1,
                                           label='Acceptable' if 'Acceptable' not in axs[0].get_legend_handles_labels()[1] else "")
                                  
                # Prepare valid xticks based on the length of taus_doppler
                log_ticks = np.logspace(np.log10(taus_doppler[0]), np.log10(taus_doppler[-1]), num=10)
                valid_xticks = np.searchsorted(taus_doppler, log_ticks)
                valid_xticks = valid_xticks[valid_xticks < len(taus_doppler)]
    
                # Set x-ticks on log scale
                axs[0].set_xticks(taus_doppler[valid_xticks])
                axs[0].set_xticklabels([f"{int(tick)}" for tick in taus_doppler[valid_xticks]])
                axs[0].legend()
                
                sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=np.min(norm_weights), vmax=np.max(norm_weights)))
                sm.set_array([])  # Required for ScalarMappable
                cbar = fig.colorbar(sm, ax=axs[0], alpha = 0.4)
                cbar.set_label('Weight')
    
                # Save the figure if a directory is specified
                if save_dir:
                    if os.path.exists(os.path.join(save_dir, 'parameters_error_bounds_plot.png')):
                        print(f'Removing {os.path.join(save_dir, "parameters_error_bounds_plot.png")}')
                        os.remove(os.path.join(save_dir, 'parameters_error_bounds_plot.png'))
                    plt.savefig(os.path.join(save_dir, 'parameters_error_bounds_plot.png'))
                    print(f'Saving {os.path.join(save_dir, "parameters_error_bounds_plot.png")}\n')
                    plt.close(fig)
    
                if suppress == True:
                    plt.close(fig)  # Close the figure 
                else:
                    plt.show()
                    plt.close(fig)
            else:
                print(f'Cannot plot the data for File due to the reasons explained above. Skipping...\n')

        ########################################################################################################################################
########################################################################################################################################
                
        def get_allan_index(self, extracted_data, save_dir=None, suppress=False):

            """ 
            This function computes the Allan Index, based on the two arrays:

            1) is_close = array made of as many boolean values (true or false) as the number of the fdets data points
            2) average_weights = average weight between two consecutive data points

            Input: object.ProcessFdets().extract_parameters [required]
                   save_dir (to create and save a .txt file report) [optional, default is None]
                   suppress (to suppress the function output, useful when used within Get_All_Plots) [optional, default is False]
                   
            Output: allan_index.txt and/or allan_index value

            NOTES: Please note that, within this function, the slope is computed between data points only 
                   (differently from the Get_Slope_At_Tau function, where a cubic spline interpolation is used,
                   and hence one can retrieve slopes at any time). 
                   The reason for this is that the error bars (which are used to compute 
                   the minimum and maximum slope (two-sided triangle) are only available for the data points, 
                   so it would be meaningless to interpolate the data, 
                   as we do not actually know how to interpolate the error bars.
        
            """

            if extracted_data != None:

                target_slope = -1
            
                # Extract data
                utc_datetime = extracted_data['utc_datetime']
                signal_to_noise = extracted_data['signal_to_noise']
                doppler_noise_hz = extracted_data['doppler_noise_hz']
                first_col_name = extracted_data['first_col_name']
                second_col_name = extracted_data['second_col_name']
                fifth_col_name = extracted_data['fifth_col_name']
                utc_date = extracted_data['utc_date']
                base_frequency = extracted_data['base_frequency']
                frequency_detection = extracted_data['frequency_detection']
                            
                # Calculate sampling rate in Hz
                t_jd = [Time(time).jd for time in utc_datetime]

        
                # Calculate the differences
                diffs = np.diff(t_jd)
                
                # Get the most common difference
                most_common_diff = Counter(diffs).most_common(1)
                if most_common_diff and most_common_diff[0][0] != 0:
                    rate_fdets = 1 / (most_common_diff[0][0] * 86400)
                else:
                    print("Most common difference is zero or not found; cannot calculate rate_fdets. Hence, no plot is available.")
                    return(None)
                # Proceed to calculate modified Allan deviation, ensuring no invalid values
                try:
                    taus_doppler, mdev_doppler, errors, ns = allantools.mdev(
                        np.array(doppler_noise_hz) / (np.array(frequency_detection) + base_frequency),
                        rate=rate_fdets,
                        data_type='freq',
                        taus='all'
                    )
                except Exception as e:
                    print(f"An error occurred: {e}")

                # Defining Weights
                rms_values = errors 
                weights = 1 / (np.array(rms_values) ** 2)
    
                # Normalize the weights to a range between 0 and 1
                norm_weights = (weights - np.max(weights)) / (np.min(weights) - np.max(weights))
                cmap = cm.get_cmap('plasma')  # You can choose any colormap
    
                # Generate the white noise reference line
                mdev_white = [mdev_doppler[0]]  # Initialize with the first value
                for i in range(1, len(taus_doppler)):
                    mdev_white.append(mdev_doppler[0] * (taus_doppler[i] / taus_doppler[0])**(target_slope))
    
                # Calculate slopes and determine is_close
                slopes = []
                is_close = []
                mean_weights = []
                target_slope = -1
                threshold = 0.01  # Define a threshold for closeness
    
                num_points = len(taus_doppler) - 1
                
                for i in range(num_points):
                    slope = (np.log10(mdev_doppler[i + 1]) - np.log10(mdev_doppler[i])) / (np.log10(taus_doppler[i + 1]) - np.log10(taus_doppler[i]))
                    slopes.append(slope)
                
                    # here, for slope_error_minus, we are creating a "double triangle" connecting 
                    # the lowest value at i with the highest at i+1 and checking the slope.
                    # same for slope_error_plus, but with interchanged endpoints
                    
                    slope_error_plus = ((np.log10(mdev_doppler[i + 1] - errors[i + 1]) - np.log10(mdev_doppler[i] + errors[i])) /
                                       (np.log10(taus_doppler[i + 1]) - np.log10(taus_doppler[i])))
                    slope_error_minus = ((np.log10(mdev_doppler[i + 1] + errors[i + 1]) - np.log10(mdev_doppler[i] - errors[i])) /
                                       (np.log10(taus_doppler[i + 1]) - np.log10(taus_doppler[i])))
                    
                    # Check if the slope is close to the target slope
                    is_close.append((slope_error_plus <= target_slope and target_slope <= slope_error_minus) or (slope_error_minus <= target_slope and target_slope <= slope_error_plus) or abs(slope - target_slope) <= 0.1)
                    
                    mean_weights.append((norm_weights[i] + norm_weights[i+1])/2) #the mean of the two weights at endpoints 
                                                                                # is considered for the interval.
                    
                    # Convert lists to NumPy arrays
                    is_close_array = np.array(is_close)
                    mean_weights_array = np.array(mean_weights)
    
                
                # Convert is_close to float (1 for True, 0 for False) for multiplication
                is_close_float = is_close_array.astype(float)
                
                # Element-wise multiplication
                product_array = is_close_float * mean_weights_array
                self.allan_index = np.sum(product_array)/(num_points -1)
    
    
                # Prepare data to save in the required format
                #is_close_str = f"is_close: [{', '.join(map(str, is_close_array))}]"
                #mean_weights_str = f"average_weights:  [{', '.join(map(str, mean_weights_array))}]"
                allan_index_str = f"allan_index = {self.allan_index}"
    
            
                # Save to file if save_dir is specified
                if save_dir:
                    # Ensure the directory exists
                    os.makedirs(save_dir, exist_ok=True)
                    file_path = os.path.join(save_dir, 'allan_index.txt')
                    if os.path.exists(file_path):
                        print(f'Removing {file_path}')
            
                    # Save data to file
                    with open(file_path, 'w') as f:
                        #f.write(f"{is_close_str}\n")
                        #f.write(f"{mean_weights_str}\n")
                        f.write(f"{allan_index_str}\n")
                        
                    print(f'Saved Allan index data to: {file_path}\n')
  
                if suppress == False:
                    return (self.allan_index)
            else:
                print(f'Cannot compute Allan Index due to the reasons explained above. Skipping ...')


########################################################################################################################################
########################################################################################################################################

        def Get_All_Outputs(self, root_folder, save_index = False, save_plots = False): # Function to process and save plots for each TXT file

            """
            
            This function iterates over all folders for a given mission dataset and retrieves the 
            plot_parameters and plot_parameters_error_bounds plots for (almost all) fdets. 

            Some exceptions might be:
            - fdets with a bad formatting
            - fdets containing multiple scans
            - fdets with a bad header

            Inputs: root_folder (namely, dataset) [required]
                    JUICE [optional] (if JUICE == True, JUICE format is assumed. if JUICE == False, other missions format is assumed)

            Please note: For mex Phobos flyby, only the file complete.r*i.txt are processed
                    
            """

            print(f'Getting Outputs From Folder: {root_folder} ...\n')
            # Compile the regex pattern to match filenames
            pattern = re.compile(r'00\d.+r.*i.txt$') 
            pattern_mex = re.compile(r'complete.+r.*i.txt$') 

            # Iterate through all directories and subdirectories
            for dirpath, _, filenames in os.walk(root_folder):
                for filename in filenames:
                    # Skip files based on certain conditions
                    if 'Phases' in filename or 'sum' in filename:
                        continue
                    # Check if the filename matches the pattern and ends with .txt
                    if pattern.search(filename):
                        print(f'Processing Directory:{dirpath}')
                        if filename.endswith('.txt'):
                        # Full path to the TXT file
                            txt_file_path = os.path.join(dirpath, filename)
    
                            
                            try:
                                # Extract data from the text file
                                extracted_data = self.extract_parameters(txt_file_path)
                            
                                # Check if start_date_time contains a space
                           #     if ' ' in extracted_data["utc_datetime"]:
                           #         print(f"Skipping file {filename} due to invalid start_date_time format.\n")
                           #         continue
                           #         
                           # except IndexError:
                           #     print(f"Skipping file: {filename} due to insufficient lines in the file.\n")
                           #     continue
                            except Exception as e:
                                print(f"Error processing file {filename}: {e}")
                                continue
                            # Create a directory for saving plots
                            plot_dir = os.path.join(dirpath, os.path.splitext(filename)[0])  # Use the TXT file name without extension
                            os.makedirs(plot_dir, exist_ok=True)
            
                            if save_plots == True:  
                                # Generate and save plots, suppresses output so as not to show all the plots
                                print(f'Generating Plot File...\n')
                                self.plot_parameters(extracted_data, plot_dir, suppress = True)
                                self.plot_parameters_error_bounds(extracted_data, plot_dir, suppress = True)

                            if save_index == True:
                                # Generate and save allan index report
                                print(f'Generating Allan Index File...')
                                self.get_allan_index(extracted_data, plot_dir, suppress = True)
                            
                    if pattern_mex.search(filename): #mex phobos flyby has slightly different names and problematic files
                        print(f'Processing Directory:{dirpath}')
                        if filename.endswith('.txt'):
                        # Full path to the TXT file
                            txt_file_path = os.path.join(dirpath, filename)
                            
                            try:
                                # Extract data from the text file
                                extracted_data = self.extract_parameters(txt_file_path)
                            
                           #     # Check if start_date_time contains a space
                           #     if ' ' in extracted_data["utc_datetime"]:
                           #         print(f"Skipping file {filename} due to invalid start_date_time format.")
                           #         continue
                           #         
                           # except IndexError:
                           #     print(f"Skipping file {filename} due to insufficient lines in the file.")
                           #     continue
                            except Exception as e:
                                print(f"Error processing file {filename}: {e}")
                                continue
                            # Create a directory for saving plots
                            plot_dir = os.path.join(dirpath, os.path.splitext(filename)[0])  # Use the TXT file name without extension
                            os.makedirs(plot_dir, exist_ok=True)

                            if save_plots == True:  
                                # Generate and save plots, suppresses output so as not to show all the plots
                                print(f'Generating Plot File...\n')
                                self.plot_parameters(extracted_data, plot_dir, suppress = True)
                                self.plot_parameters_error_bounds(extracted_data, plot_dir, suppress = True)

                            if save_index == True:
                                # Generate and save allan index report
                                print(f'Generating Allan Index File...')
                                self.get_allan_index(extracted_data, plot_dir, suppress = True)
 
                        
            print(f'...Done.\n') 

########################################################################################################################################
########################################################################################################################################
            
        #Function to compute the averaged Modified Allan Deviation Slope for a given tau (in seconds)
        def Get_Slope_At_Tau(self, extracted_data, tau, delta_tau):

            """

            This function computes the mADEV slope at a given tau, using numerical derivatives of the type:

            slope = [mADEV(tau) + mADEV(tau+dt)]/dt

            The choice of tau and dt is somehow arbitrary and is given as input

            Inputs: object.ProcessFdets().extract_parameters [required]
                    tau [required] (in seconds)
                    dt [required] (in seconds)
            
            """
        
            utc_datetime = extracted_data['utc_datetime']
            base_frequency = extracted_data['base_frequency']
            frequency_detection = extracted_data['frequency_detection']
            doppler_noise_hz = extracted_data['doppler_noise_hz']
        
          # Convert UTC datetime to Julian Date
            t_jd = [Time(time).jd for time in utc_datetime]
            # Calculate the differences
            diffs = np.diff(t_jd)
            
            # Get the most common difference
            most_common_diff = Counter(diffs).most_common(1)
            if most_common_diff and most_common_diff[0][0] != 0:
                rate_fdets = 1 / (most_common_diff[0][0] * 86400)
            else:
                print("Most common difference is zero or not found; cannot calculate rate_fdets. Hence, no plot is available.")
                return(None)
            
            
            # Proceed to calculate modified Allan deviation, ensuring no invalid values
            try:
                taus_doppler, mdev_doppler, errors, ns = allantools.mdev(
                    np.array(doppler_noise_hz) / (np.array(frequency_detection) + base_frequency),
                    rate=rate_fdets,
                    data_type='freq',
                    taus='all'
                )
            except Exception as e:
                print(f"An error occurred: {e}")
                
            # Calculate Modified Allan Deviation for Doppler noise
            taus_doppler, mdev_doppler, errors, ns = allantools.mdev(
                np.array(doppler_noise_hz) / (np.array(frequency_detection) + base_frequency),
                rate=rate_fdets,
                data_type='freq',
                taus='all'
            )
            # Ensure tau and interpol_time are within the range of taus_doppler
            if tau >= max(taus_doppler) or (tau + delta_tau) >= max(taus_doppler):
                print("Most common difference is zero or not found; cannot calculate rate_fdets. Hence, no plot is available.")
                return(None)
                
            
            # Interpolation
            interpolated_data = CubicSpline(taus_doppler, mdev_doppler, extrapolate=True)
            
            # Compute values at tau and tau + delta_tau
            mdev_tau = interpolated_data(tau)
            mdev_tau_delta = interpolated_data(tau + delta_tau)
            
            # Compute the slope in the log-log plot
            if mdev_tau == 0 or mdev_tau_delta == 0:
                print("MDEV values at tau or tau + delta_tau cannot be zero.")
                return(None)
            
            self.slope = (np.log10(mdev_tau_delta) - np.log10(mdev_tau)) / (np.log10(tau + delta_tau) - np.log10(tau))
        
            return(self.slope)
