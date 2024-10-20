# this program has been implemented to run under python 3.12. If not working, first check your python version.
# make sure all the python modules used are installed on your computer.
#
# PURPOSE :
#       using the bias corrected and delta changed data files computed by BiasCorrectPrecip.R and BiasCorrectTemperature.R,
#       this program will automatically create the text setup files to run simulations under PINEHBV. This is the same purpose
#       than csv_to_txt_Sildre_SeNorge.py but with a different input.
#
# INPUT :
#   variables :
#       catchment : name of the catchment. Be careful with upper case letters.
#       altitude : altitude of the measurement point of the catchment. Can be found by csv_to_txt_Sildre_SeNorge.py or at
#                  https://sildre.nve.no/ section 'about the station'
#
#       biascorrect_output : the output of the BiasCorrect functions.
#       pine_input_corr : name of the repository (under PINE/PineProj/Catchment) where will be created the pine text input files of
#                   all climate models, scenario, periods for the corrected data sets. If pine_input_corr doesn't exist already,
#                  it will be created by this program.
#       pine_input_delta : same as pine_input_corr but for the delta changed data sets.
#       biascorr_idd : file name identifier for data sets corrected by BiasCorrect (have to be the same for temp and prec)
#       deltachanged_idd : file name identifier for data sets deltachanged by BiasCorrect (have to be the same for temp and prec)
#
#   files :
#       temp and prec files coming from biascorrect.R, containing the data of the studied climate models :
#       CorrectedP_histo_histperiod.csv
#       CorrectedT_histo_histperiod.csv (histperiod has the format YYYY-YYYY)
#       CorrectedP_scenario_period.csv
#       CorrectedT_scenario_period.csv (scenario has lenght 5 and period has the format YYYY-YYYY)
#       DeltaChangedP_scenario_period.csv
#       DeltaChangedT_scenario_period.csv
#
# OUTPUT :
#   files : under pine_input_path_corr and pine_input_path_delta
#       txt files with the good format to setup a PINE model. The discharge is unknown in the future so it is set to zero.
#       do it for near future (starY = 2041, endY = 2070) and far future (startY = 2071, endY = 2100),
#       just as the files coming from biascorrect are named.
#
# CAUTION :
#       identifiers for temp and prec files have to be at the beginning of the filenames
#
# NB :
#       The baseline (= historical values) will have the following name : Catchment_Observed_histo_histperiod.txt
#       histperiod has the following format : YYYY-YYYY
#       the discharge is set top zero as we do not have climate values for discharge. Simulated discharge will be computed by PINEHBV.
#       minimum and maximum temperature are both equal to daily average temperature.
#
# ------------------------------------ SETUP -----------------------------------
catchment = "Myglevatn"
altitude = 256

biascorrect_output = 'Output_BiasCorrect'
pine_input_corr = "input_corr"
pine_input_delta = "input_delta"
biascorr_idd = "corrected"
deltachanged_idd = "deltachanged"
# ------------------------------------ SETUP -----------------------------------

import os
import inspect

# Get the current frame
current_frame = inspect.currentframe()
# Get the file name of the current frame
script_path = inspect.getfile(current_frame)
# Get the absolute path of the script
absolute_script_path = os.path.abspath(script_path)
# Get the directory name of the script
script_dir = os.path.dirname(absolute_script_path)
# Replace backslashes with forward slashes for uniformity
script_dir = script_dir.replace('\\', '/')
# Delete all characters after the last '/'
klima = script_dir[:script_dir.rfind('/') + 1]

input_path = f'{klima}{biascorrect_output}/{catchment}' # Base directory path with a placeholder for the capital letter
# pattern for the name of the folder where the setups will be stored
pine_input_path_corr = f'{klima}PINE/PineProj/{catchment}/{pine_input_corr}'
pine_input_path_delta = f'{klima}PINE/PineProj/{catchment}/{pine_input_delta}'

def create_directory(directory_path):
    if not os.path.exists(directory_path):
        os.makedirs(directory_path)
        print(f"Directory '{directory_path}' created.")
    else:
        print(f"Directory '{directory_path}' already exists.")
    print("")

# print basic info
print("Processing the data of the following catchment : ", catchment)
print('data coming from : ', input_path)
print(f'directory for biascorrected setup creation : ', pine_input_path_corr)
print(f'directory for deltachanged setup creation : ', pine_input_path_delta)
print("")

# Directories for PINE text input path
create_directory(pine_input_path_corr)
create_directory(pine_input_path_delta)

import re
import numpy as np
import warnings

def find_data_files(input_path_, idd):
    '''input  : the folder path containing the output data files of BiasCorrect for both temp and precip. All the simulations should have been done (if the temperature file for a given scenario and period, the associated precipitation file should be created as well).
       output : list of lists. Each sulist contains the temp file name and the prec file name for a given scenario and period (this includes historical values)
       NB : the temp and prec files should have formated names like these :
            CorrectedT_scenario_period.csv (scenario of lenght 5 and period to the format YYYY-YYYY)
            CorrectedP_histo_period.csv
    '''
    # list all the files coming from the directory
    all_files = [f for f in os.listdir(input_path) if os.path.isfile(os.path.join(input_path, f))]

    # filter the temp and prec files according to the names of the output files of BiasCorrect (if BiasCorrect is changed, this function needs to be changed as well)
    temp_files = [element for element in all_files if (idd) in element]
    temp_files = [element for element in temp_files if ('temp') in element]
    prec_files = [element for element in all_files if (idd) in element]
    prec_files = [element for element in prec_files if ('prec') in element]

    # matching temp and prec files by scenario and period, warning if missing files
    associated_temp_prec = []
    nb_temp_f = len(temp_files)
    nb_prec_f = len(prec_files)
    if nb_temp_f < nb_prec_f:
        for i in range(nb_prec_f):
            prec_file = prec_files[i]
            scenario = prec_file[-19:-14]
            period = prec_file[-13:-4]
            file_extension = prec_file[-4:]

            associated_temp = [temp_file for temp_file in temp_files if temp_file.endswith(scenario + '_' + period + file_extension)]
            # only one temperature file should exist for a given scenario and period
            if len(associated_temp) == 0 :
                warnings.warn("No associated temperature file for " + scenario + '_' + period + "The setup will be empty for temperature")
            elif len(associated_temp) > 1 :
                warnings.warn("To much associated temperature files for " + scenario + '_' + period)
            else :
                associated_temp_prec.append([associated_temp[0], prec_file])
    else :
        for i in range(nb_temp_f):
            temp_file = temp_files[i]
            scenario = temp_file[-19:-14]
            period = temp_file[-13:-4]
            file_extension = temp_file[-4:]

            associated_prec = [prec_file for prec_file in prec_files if prec_file.endswith(scenario + '_' + period + file_extension)]
            # only one precipitation file should exist for a given scenario and period
            if len(associated_prec) == 0 :
                warnings.warn("No associated precipitation file for " + scenario + '_' + period + "The setup will be empty for precipitation")
            elif len(associated_prec) > 1 :
                warnings.warn("To much associated precipitation files for " + scenario + '_' + period)
            else :
                associated_temp_prec.append([temp_file, associated_prec[0]])
    # only complete datasets are returned
    return associated_temp_prec

import csv
from datetime import datetime, timedelta

def csv_to_array(file_path):
    '''input : path to a csv file containing the corrected data from BiasCorrectPrecip.R or BiasCorrectTemperature.R
       output : array of (in column) : the date, the data of the climate models sorted alphabetically. The data is being converted as floats and dates as datetime objects'''
    with open(file_path, newline='', encoding='utf-8') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        data_list = list(reader)

    # Extract headers
    headers = data_list[0][1:]

    # Define dtype with correct types (the type of data in each column)
    dtype = [('date', 'U10')] + [(header, float) for header in headers[1:]]

    # Convert data to appropriate types : change the format of the date and convert string to float
    date_list = [datetime.strptime(row[1], '%Y-%m-%d') for row in data_list[1:]]
    numeric_data = [[float(val) for val in row[2:]] for row in data_list[1:]]

    # Create a structured array with the specified dtype
    structured_data = np.array([(date,) + tuple(data) for date, data in zip(date_list, numeric_data)], dtype=dtype)

    # Sort the data columns based on the headers
    sorted_indices = np.argsort(headers[1:])
    sorted_headers = np.array(headers[1:])[sorted_indices]
    sorted_data = structured_data[list(sorted_headers)]
    final_headers = np.insert(sorted_headers, 0, 'date')

    # Combine the date_list and sorted_array
    concatenated_array = np.column_stack((date_list, np.column_stack([sorted_data[header] for header in sorted_headers])))
    concatenated_array = np.vstack((final_headers[np.newaxis, :], concatenated_array))
    start, end = concatenated_array[1,0], concatenated_array[-1,0]

    return concatenated_array, start, end

def find_info(filename, scenarios):
    '''
    input :
      filename : string
      scenarios : list of studied scenarios to match filename's scenario
    output :
      the scenario, the period, the data type, and the modification type of filename
    '''

    # Use a regular expression to find the period pattern
    period_pattern = r"\d{4}-\d{4}"

    # Extract the period using the regular expression
    period_match = re.search(period_pattern, filename)
    period = period_match.group() if period_match else ""

    # Replace the underscore with a hyphen
    period_formatted = period.replace("_", "-")

    # Initialize variable for the scenario
    found_scenario = None

    # Loop through each scenario to find a match in the input string
    for scenario in scenarios:
        if scenario in filename:
            found_scenario = scenario
            break

    corr_delta = ['Corrected', 'DeltaChanged', 'corrected', 'deltachanged']
    # Initialize variable for the modification type
    found_modif = None

    # Loop through each modification type to find a match in the input string
    for modif in corr_delta:
        if modif in filename:
            found_modif = modif
            break

    # Return scenario, formatted period, data type and 'corrected'/'deltachanged'
    return found_scenario, period_formatted, found_modif


def array_to_txtsetups(pine_input_path, temp_file_path, prec_file_path, catchment_, idd, altitude_=''):
    ''' input :
            temp_file_path : the path of the temperature daily data for every climate model coming from BiasCorrectTemperature.R
            prec_file_path : the path of the precipitation daily data for every climate model coming from BiasCorrectPrecip.R
            catchment : name or iddentifier of the catchment
            idd : biascorrected or deltachanged
        output :
            pine_input_path : the path to the file that will be created for pine setup input data
        '''
    if (altitude_ == ''):
        warnings.warn("The altitude of the catchment is missing")
    temp_array, start_T,end_T = csv_to_array(temp_file_path)
    prec_array, start_P, end_P = csv_to_array(prec_file_path)
    if start_P != start_T or end_P != end_T:
        warnings.warn("The period for Temperature data and Precipitation data is not the same. Please check the provided output of BiasCorrect :" + temp_file_path + " and " + prec_file_path)
    scenario, period, modif = find_info(temp_file_path, ['rcp26','rcp45', 'rcp85', 'histo'])
    # Get the name of the models and delete the '_', '.' and '-' to get shorter names
    models = [model.replace('_','').replace('.','').replace('-','') for model in temp_array[0]]

    # Assign data from temp_array and prec_array at setup_path
    for s in range(1, len(models)):
        # Take the name of the model but cut it if it's too long (>15 char)
        current_model = models[s][:min(15,len(models[s]))]
        setup_name = f"{catchment_[:5]}_{idd[:5]}_{current_model}_{scenario}_{period}.txt".lower()
        setup_path = f"{pine_input_path}/{setup_name}"
        if os.path.isfile(setup_path):
            print("Pine input text file already exists       : ", setup_name)
        else :
            with open(setup_path, 'w') as file:
                pass
            with open(setup_path, 'w') as file:
                # Write the headers and the second line
                file.write(f"Dato\tTime\tP_mm_{altitude}moh\tT_minC_{altitude}moh\tT_maxC_{altitude}moh\tQ_{catchment}\n")
                file.write("dd.mm.yyyy\thh.mm.ss\tmm\tgrC\tgrC\tm3/s\n")
                for i in range(1,len(temp_array[:,0])):
                    # Convert each field to a string and join with tabs. Discharge is set to zero and time to '00.00.00'.
                    row = [temp_array[i,0].strftime('%d.%m.%Y'), '00.00.00', round(prec_array[i,s], 6), round(temp_array[i,s],6), round(temp_array[i,s],6), 0.]
                    file.write('\t'.join(map(str, row)) + '\n')
            print(f'Pine input text file created successfully : ', setup_name)

    print('number of setups for this scenario and period : ', len(models)-1)
    return len(models)-1

# Corrected setup creation
# look for the data files with the good format of name using the identifiers
corrected_temp_prec_filenames = find_data_files(input_path, biascorr_idd) # list of list of filenames [[T,P],[T,P],...]
nb_setups_corr = 0
# loop over all the files found at the given address
for i in range(len(corrected_temp_prec_filenames)):
    print(corrected_temp_prec_filenames[i][0][-19:-4])
    temp_path = input_path + '/' + corrected_temp_prec_filenames[i][0]
    prec_path = input_path + '/' + corrected_temp_prec_filenames[i][1]
    count_setups = array_to_txtsetups(pine_input_path_corr, temp_path, prec_path, catchment, biascorr_idd, altitude)
    nb_setups_corr += count_setups
    print('')
print('total of biascorrected setups found or created : ', nb_setups_corr)

# DeltaChanged setup creation
# look for the data files with the good format of name using the identifiers
deltachanged_temp_prec_filenames = find_data_files(input_path, deltachanged_idd) # list of list of filenames [[T,P],[T,P],...]
nb_setups_delta = 0
# loop over all the files found at the given address
for i in range(len(deltachanged_temp_prec_filenames)):
    print(deltachanged_temp_prec_filenames[i][0][-19:-4])
    temp_path = input_path + '/' + deltachanged_temp_prec_filenames[i][0]
    prec_path = input_path + '/' + deltachanged_temp_prec_filenames[i][1]
    count_setups = array_to_txtsetups(pine_input_path_delta, temp_path, prec_path, catchment, deltachanged_idd, altitude)
    nb_setups_delta += count_setups
    print('')
print('total of deltachanged setups found or created : ', nb_setups_delta)


### test
temp_array = csv_to_array(temp_path)
prec_array = csv_to_array(prec_path)














