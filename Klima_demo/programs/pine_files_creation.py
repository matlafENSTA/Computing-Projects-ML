# this program has been implemented to run under python 3.12. If not working, first check your python version.
# make sure all the python modules used are installed on your computer.
#
# PURPOSE :
#       automatically create the parameter files, pine setup files and pine output files for all scenarios, periods and climate
#       models and fill the evapotranspiration values in the parameter files. It saves time for the PINE setup creation routine.
#       all the values should be computed by TemperatureBiasCorrect.R
#
# INPUT :
#   variables :
#       catchment : name of the catchment. Be careful with the uppercase letters.
#       setuphps_corr_dir : name of the directory for pine setup files (biascorrected).
#       setuphps_delta_dir : name of the directory for pine setup files (deltachanged).
#       param_corr_dir : directory for PINE parameter files (can be the same than setuphps_corr_dir)
#       param_delta_dir : directory for PINE parameter files (can be the same than setuphps_delta_dir)
#       pineinput_corr_dir : directory containing the input text setup files to run simulations. Data comes from BiasCorrect.
#       setup_file_baseline : name of the baseline .hps setup file. Put it under pine_catchment_path.
#       param_file_baseline : name of the baseline parameter file. Have to be completed during the PINE routine for the baseline.
#               Put it under pine_catchment_path.
#
#       pineoutput_corr_dir : directory for pine output files using bias corrected setups
#       pineoutput_delta_dir : directory for pine output files using delta changed setups
#                              (will contain the final values for snowpacks and discharge)
#       input_folder : folder with the monthly evapo files
#       recapfile : recap of all the setup files handled in
#
#   files :
#       the pine parameter file file for the baseline (under pine_catchment_path).
#       all the pine text input files created by corrected_txt_pine_input.py (under pine_catchment_path/pineinput_corr_dir for
#       bias corrected values and pine_catchment_path/pineinput_delta_dir for delta changed values)
#
# OUTPUT :
#       recapfile : recap of all the setup files handled (located at pine_catchment_path)
#       in setup_path_corr, the .hps setup files corresponding to the bias corrected .txt setups created by corrected_txt_pine_input.py
#       in setup_path_delta, the .hps setup files corresponding to the delta changed .txt setups created by corrected_txt_pine_input.py
#       in param_path_corr, the .top parameter files associated to the setup of the same name in setup_path_corr
#       in param_path_delta, the .top parameter files associated to the setup of the same name in setup_path_delta
#       in pineoutput_path_corr, the pine output files (.txt and .dat) associated to the setup of the same name in setup_path_corr
#       in pineoutput_path_delta, the pine output files (.txt and .dat) associated to the setup of the same name in setup_path_delta
#
# CAUTION :
#       to run, this script needs :
#           - the output of TemperatureBiasCorrect.R (DailyPET files for the concerned periods and scenarii)
#           - the txt setups computed by corrected_txt_pine_input.py
#           - an example of parameter file, with the catchment values filled, to be found at pine_catchment_path
#
# NB :
#       parameter files are the same for biascorrected setups than for deltachanged setups,
#       but pine requires one parameter file for one setup ; that is why we copy them.
#       pine output files are empty, you need to run PINE to fill them with the simulated data.
#       creating these empty files just accelerates the process of pine setup creation.
#
# ------------------------------------ SETUP -----------------------------------
catchment = "Myglevatn"

pineinput_corr_dir = "input_corr"
pineinput_delta_dir = "input_delta"
setuphps_corr_dir = "setups_hps_corr"
setuphps_delta_dir = 'setups_hps_delta'
param_corr_dir = "setups_hps_corr"
param_delta_dir = "setups_hps_delta"
setup_file_baseline = "Myglevatn_obs_setup.hps"
param_file_baseline = "Myglevatn_obs_par.top"

pineoutput_corr_dir = "output_corr"
pineoutput_delta_dir = 'output_delta'

input_folder = 'Output_BiasCorrect'
recapfile = "list of setups.txt"

# in case of bug, try to look at these variables :
# month identifier for evapotranspiration values in the parameter files :
months_identifier = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OKT', 'NOV', 'DES']
pattern_setup = r'.*0\.txt$' # to modify in consequence if you choose to name the parameter files differently
pattern_param_files = fr'^{catchment[:5].lower()}.*par\.top$' # to modify in consequence if you choose to name the PINE parameter files differently
pattern_evapo = r'^DailyPET.*0\.csv$' # to modify in consequence if you choose to name the evapotranspiration files differently,
# or compute a period that doesn't end with 0 (like 2041-2070 for example)
# ------------------------------------ SETUP -----------------------------------

import inspect
import os

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

# This part puts back upper letters in the path klima. If not working, just delete it.
# Convert the path to match desired casing
klima_parts = klima.split('/')
klima_parts = [part.capitalize() if part.lower() == 'klima' else part for part in klima_parts]
klima = '/'.join(klima_parts)
# Handle drive letter casing separately
if klima[1] == ':':
    klima = klima[0].upper() + klima[1:]

input_path = f'{klima}{input_folder}/{catchment}/' # Base directory path with a placeholder for the capital letter
# pattern for the name of the folder where the setups will be stored
pine_catchment_path = f'{klima}PINE/PineProj/{catchment}'
setup_path_corr = f'{pine_catchment_path}/{pineinput_corr_dir}/'
param_path_corr = f'{pine_catchment_path}/{param_corr_dir}/'
setuphps_corr_path = f'{pine_catchment_path}/{setuphps_corr_dir}/'
pineoutput_path_corr = f'{pine_catchment_path}/{pineoutput_corr_dir}/'
setup_path_delta = f'{pine_catchment_path}/{pineinput_delta_dir}/'
param_path_delta = f'{pine_catchment_path}/{param_delta_dir}/'
setuphps_delta_path = f'{pine_catchment_path}/{setuphps_delta_dir}/'
pineoutput_path_delta = f'{pine_catchment_path}/{pineoutput_delta_dir}/'
# if the working directory containing all the setup files has changed, just rename it below. Same for all other folders
# setup_path =
# param_path =
# setuphps_path =
# pineoutput_path =

# print basic info
print("Processing the data of the following catchment : ", catchment)
print('data coming from : ', input_path)
print("")

def create_directory(directory_path):
    if not os.path.exists(directory_path):
        os.makedirs(directory_path)
        print(f"Directory '{directory_path}' created.")
    else:
        # print(f"Directory '{directory_path}' already exists.")
        pass

# directory for setup parameter files, output
create_directory(param_path_corr)
create_directory(setuphps_corr_path)
create_directory(pineoutput_path_corr)
print("")
create_directory(param_path_delta)
create_directory(setuphps_delta_path)
create_directory(pineoutput_path_delta)
print("")

import re

with open(pine_catchment_path + '/' + recapfile , "w") as setup_recap:
    setup_recap.write(f"{catchment} recap file \nList of setups successfully created for this catchment :")

def find_files(input_path_, pattern = ''):
    '''input_path_  : the folder path with files you want to work with
       pattern : a string to recognize the files (only the files matching the pattern will be selected)
       output : list of concerned files in the folder
    '''
    pattern = re.compile(pattern)
    #list all the files in the directory
    file_list = [f for f in os.listdir(input_path_) if os.path.isfile(os.path.join(input_path_, f))]
    if pattern == None:
        return sorted(file_list)
    else :
        filtered_list = [filename for filename in file_list if pattern.search(filename)]
        return sorted(filtered_list)

# list of txt setups created by corrected_txt_pine_input.py using TemperatureBiasCorrect.R output
txt_setups_corr = find_files(setup_path_corr, pattern_setup)
txt_setups_corr = [setup.lower() for setup in txt_setups_corr] # lower cases only to fit PINEHBV's typo
txt_setups_delta = find_files(setup_path_delta, pattern_setup)
txt_setups_delta = [setup.lower() for setup in txt_setups_delta]

evapo_files_list = find_files(input_path, pattern_evapo) # list of files containing the evapotranspiration coming from TemperatureBiasCorrect.R
print("number of text files for pine setup found (biascorrected): ", len(txt_setups_corr))
print("number of text files for pine setup found (deltachanged): ", len(txt_setups_delta))
print("number of evapotranspiration data files : ", len(evapo_files_list))


import shutil
import warnings

def update_setuphps(setuphps_dir, setuphps_name, param_filepath, dat_filepath, outtxt_filepath, outdat_filepath, end_filepath):
    '''INPUT :
            setuphps_dir : path to the directory where the setup you want to modify is in
            setuphps_name : name of the hps setup you want to complete (has to be already created by copying an existing setup)
            param_filepath : path of the parameter file for this setup
            dat_filepath : path of the binary pine input files for this setup
            outdat_filepath : path of the binary pine output files for this setup
            outtxt_filepath : path of the text pine output files for this setup
            end_filepath : path of the .top pine end files for this setup
        OUTPUT :
            the setuphps file has been uploaded to be ready for a PINE simulation'''
    # Read the file content into memory
    setuphps_path = setuphps_dir + '/' + setuphps_name
    with open(setuphps_path, 'r') as file:
        lines = file.readlines()

    # Update the name of the setup and the catchment
    # Find concerned lines
    setupname_index = next((i for i, s in enumerate(lines) if 'Catchment' in s), None)
    catchname_index = next((i for i, s in enumerate(lines) if 'Catchment: ' in s), None)
    # Change their values
    lines[setupname_index] = 'Catchment and data files for ' + setuphps_name[:-4].upper() + '-forecast setup.\n'
    lines[catchname_index] = 'Catchment: ' + setuphps_name[:-10] + '\n'

    # Function to find the index and leading spaces
    def find_index_and_spaces(keyword):
        for i, s in enumerate(lines):
            if keyword in s:
                spaces = ''.join(char for char in s if char.isspace() and char == ' ')
                return i, spaces
        return None, None

    # List of file keywords and new paths
    file_updates = {
        'par.top': param_filepath,
        'out.dat': outdat_filepath,
        'out.txt': outtxt_filepath,
        'end.top': end_filepath,
        '.dat': dat_filepath  # Ensure .dat is updated correctly
    }

    for keyword, new_path in file_updates.items():
        index, spaces = find_index_and_spaces(keyword)
        if index is not None:
            lines[index] = spaces + new_path.replace('/', '\\') + '\n'

    # Write the updated content back to the file
    with open(setuphps_path, 'w') as file:
        file.writelines(lines)

    return lines, catchname_index, setupname_index


def create_setup_files(txt_setups, pine_catchment_path_, param_path_, setuphps_path, setuptxt_path, pineoutput_path):
    '''copy the parameter file of the baseline into a different parameter file for all climate models, scenario, period (i.e. for all the txt setup files in the directory for setups). At the same time, create empty files for the outputs in order to simplify the tasks with PINE
'''
    count_setupshps = 0
    count_pinetextoutput = 0
    count_pinebinoutput = 0
    count_pineendoutput = 0
    # Split the setup name by '_'
    parts = txt_setups[0].split('_')
    # Join the parts before the second underscore to get 'corrected' or 'deltachanged'
    if len(parts) > 2:
        corr_delta = ''.join(parts[1])
    else:
        corr_delta = txt_setups[0]

    with open(pine_catchment_path_ + '/' + recapfile , "a") as setup_recap:
        setup_recap.write("\n" + corr_delta + " setup handled\n")
        for current_setup_name in txt_setups:
            current_idd = current_setup_name[:-4]
            setup_recap.write(current_idd + "\n")

            current_param_name = current_idd + "_par.top"
            current_setuphps_name = current_idd + "_setup.hps"
            current_pinetextoutput_name = current_idd + "_out.txt"
            current_pinebinoutput_name = current_idd + "_out.dat"
            current_pineendoutput_name = current_idd + "_end.top"

            # do a copy of the parameter file of the baseline for every model
            # evapo values will then be uploaded using update_evapo_values
            try :
                shutil.copy(pine_catchment_path_ + '/' + param_file_baseline, param_path_ + current_param_name)
            except FileExistsError:
                pass
            try :
                shutil.copy(pine_catchment_path_ + '/' + setup_file_baseline, setuphps_path + current_setuphps_name)
                # make the hps setup
                update_setuphps(setuphps_path, current_setuphps_name, param_path_ + current_param_name, setuptxt_path + current_setup_name.replace('.txt','.dat'), pineoutput_path + current_pinetextoutput_name, pineoutput_path + current_pinebinoutput_name, pineoutput_path + current_pineendoutput_name)
                count_setupshps += 1
            except FileExistsError:
                pass
            try:
                # Open the file in exclusive creation mode
                with open(pineoutput_path + '/' + current_pinetextoutput_name, "x"):
                    pass  # File created successfully, do nothing
                count_pinetextoutput += 1
            except FileExistsError:
                pass
            try:
                # Open the file in exclusive creation mode
                with open(pineoutput_path + '/' + current_pinebinoutput_name, "x"):
                    pass  # File created successfully, do nothing
                count_pinebinoutput += 1
            except FileExistsError:
                pass
            try:
                # Open the file in exclusive creation mode
                with open(pineoutput_path + '/' + current_pineendoutput_name, "x"):
                    pass  # File created successfully, do nothing
                count_pineendoutput += 1
            except FileExistsError:
                pass

    param_files_list = find_files(param_path_, pattern_param_files) # list of parameter files created just above

    print("number of parameter files in the folder : ", len(param_files_list), "(includes the baseline parameter file)")
    print("number of empty setup files just added : ", count_setupshps)
    print("number of pine text output files just added : ", count_pinetextoutput)
    print("number of pine binary output files just added : ", count_pinebinoutput)
    print("number of pine end files just added : ", count_pineendoutput)
    # If no file was created, it means that the files in question already existed.
    print("")
    if len(txt_setups) != len([param_file for param_file in param_files_list if param_file_baseline not in param_file]): # the parameter folder already contains an example parameter file (usually baseline)
        warnings.warn("There shoud be one parameter file for each setup files, no less, no more")
    count_list = [len(param_files_list), count_setupshps, count_pinetextoutput, count_pinebinoutput, count_pineendoutput]
    return(param_files_list, count_list)

os.remove(pine_catchment_path + '/' + recapfile)
param_files_list_corr, count_list_corr = create_setup_files(txt_setups_corr, pine_catchment_path, param_path_corr, setuphps_corr_path, setup_path_corr, pineoutput_path_corr)
param_files_list_delta, count_list_delta = create_setup_files(txt_setups_delta, pine_catchment_path, param_path_delta, setuphps_delta_path, setup_path_delta, pineoutput_path_delta)

import csv
import numpy as np

def csv_to_array(evapo_file_path):
    '''input : path of a csv file containing daily average evapotranspiration values per month
       output : tabular of these values
    '''
    with open(evapo_file_path, newline='', encoding='utf-8') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        data_list = list(reader)
    return np.array(data_list)

# write the good evapo for each month parameter file by parameter file
def update_evapo_values(param_path, param_file_name, months_, new_evapo_values):
    ''' INPUT :
        param_file_name looks like Myglevatn_CNRM_CCLM_histo_1971-2000_par.top
        months_ : the identifier of the month in the lines of evapotranspiration values (list of 12 strings)
        new_evapo_values : the daily value (monthly mean) of evapotranspiration in the same order than months_ (list of 12 integers)
        OUTPUT :
        None. The parameter file given should have been modified with new evapotranspiration values.
    '''
    assert len(months_) == len(new_evapo_values)

    # Format the new evapotranspiration values
    new_evapo_values = [f"{value:.3f}" for value in new_evapo_values]

    # Read the file content into memory
    with open(param_path + '/' + param_file_name, 'r') as file:
        lines = file.readlines()

    count_new_values = 0
    for j, pattern in enumerate(months_):
        # Create the current pattern
        current_pattern = r'\s*EP' + pattern
        #print(f"Searching for pattern: {current_pattern}")

        # Replace the specific line matching the pattern
        for i, line in enumerate(lines):
            if re.search(current_pattern, line):
                #print(f"Match found in line: {line.strip()}")
                # Build the new line with the new value of evapo
                # The position to insert the new value should be determined accurately
                # Assuming the value to replace is right-aligned in a fixed-width column
                lines[i] = re.sub(r'\d{1,2}\.\d{3}', new_evapo_values[j], line.rstrip()) + '\n'
                count_new_values += 1
                break
        #else:
            # If no break occurs (no match found)
            #print(f"No match found for pattern: {current_pattern}")

    # Write the updated content back to the file
    with open(param_path + '/' + param_file_name, 'w') as file:
        file.writelines(lines)
    if count_new_values == 12:
        print("evapotranspiration values updated successfully for ", param_file_name)
    else :
        print("only ", count_new_values, "values for evapotranspiration have been changed")
    return None

def param_files_modif_loop(input_path_, param_path, param_files_list, evapo_files_list):
    ''' use update_evapo_values to fill the evapotranspiration values for every parameter file in param_files_list'''
    count_modified_param_files = 0
    # browse all outputs of TemperatureBiasCorrect.R
    for evapo_file in evapo_files_list:
        scenario = evapo_file[-19:-14]
        print(scenario)
        period = evapo_file[-13:-4]

        #extract data into array
        daily_evapo = csv_to_array(input_path_ + evapo_file)
        print(evapo_file, np.shape(daily_evapo))
        months = daily_evapo[1:,0]
        nb_models = np.shape(daily_evapo)[1]
        # all the climate models for current file (containing the data for given scenario and period)
        for i in range(1,nb_models):
            current_evapo_values = [float(daily_value) for daily_value in daily_evapo[1:,i]]
            model = daily_evapo[0,i].lower()
            #find the parameter file associated to the current scenario, period and climate model
            filtered_list_by_model = list(filter(lambda s: model in s, param_files_list))
            filtered_list_by_scenario = list(filter(lambda s: scenario in s, filtered_list_by_model))
            filtered_list_by_period = list(filter(lambda s: period in s, filtered_list_by_scenario))
            if len(filtered_list_by_period) > 1:
                print("model : ", filtered_list_by_model)
                print("scenario : ", filtered_list_by_scenario)
                print("period : ", filtered_list_by_period)
                warnings.warn(f"Too many parameter files for {param_path}, {model}, {scenario}, {period}, you should delete the unnecessary ones")
            elif len(filtered_list_by_period) == 0:
                warnings.warn(f"No parameter file found for {model}, {scenario}, {period}")
                continue
            else :
                param_file = filtered_list_by_period[0]
                update_evapo_values(param_path, param_file, months_identifier, current_evapo_values)
                count_modified_param_files += 1
    print("amount of parameter files successfully modified : ", count_modified_param_files)
    return count_modified_param_files

print("BIAS CORRECTED SETUPS : ")
c_corr = param_files_modif_loop(input_path, param_path_corr, param_files_list_corr, evapo_files_list)
print("DELTA CHANGED SETUPS : ")
c_delta = param_files_modif_loop(input_path, param_path_delta, param_files_list_delta, evapo_files_list)

###test
evapo = csv_to_array(input_path + evapo_files_list[0])

new_evapo_values_test = [0.100000,0.2000,1.5,1.7000,0.100000,0.2000,1.5,1.7000,0.100000,0.2000,1.5,1.7000]
param_lines = update_evapo_values("D:/Klima/PINE/PineProj/Myglevatn/parameter_files", "Myglevatn_CNRM_CCLM_histo_1971-2000_par.top", months_identifier, new_evapo_values_test)

# Test the function update_setuphps
test_setuphps_dir = "D:/Klima/others/test"
test_setuphps_name = "test_setup.hps"
test_paramfile_par = "D:/Klima/others/test/test_paramfile_par.top"
test_input_dat = "D:/Klima/others/test/test_input_dat.dat"
test_output_out_txt = "D:/Klima/others/test/test_output_out.txt"
test_output_out_dat = "D:/Klima/others/test/test_output_out.dat"
test_endfile_end = "D:/Klima/others/test/test_endfile_end.top"

test_lines, catchtest, setuptest = update_setuphps(test_setuphps_dir, test_setuphps_name, test_paramfile_par, test_input_dat, test_output_out_txt, test_output_out_dat, test_endfile_end)
print(test_lines, catchtest, setuptest)








