# this program has been implemented to run under python 3.12. If not working, first check your python version.
# make sure all the python modules used are installed on your computer.
#
# PURPOSE :
#       link the output of climate_models_RCP_sfVersion.R to the input of BiasCorrect R functions
#       the data coming from the climate models and extracted with climate_models_RCP_sfVersion.R have to be classified
#       before they can be put into BiasCorrectPrecip and BiasCorrectTemperature. We suppose here that the data haven't
#       any errors (missing values or dates, different period for each file).
#
# INPUT :
#   variables :
#       catchment : name of the catchment. Be careful with upper case letters.
#       geopoint_index : the index of the geographical point where the data has been measured, it can be found using Norgeskart.no
#       scenarii : list of strings representing the names of scenarii handled by climate_models_RCP_sfVersion.R
#                  It usually is the name of the output subfolders of climate_models_RCP_sfVersion.R : ["hist","rcp45","rcp85"]
#                  !! the lenght of the string representing a scenario has to be superior or equal to 4 !!
#       startY_models, endY_models : starting year and ending year of the period for climate models data (future)
#       startY_hist, endY_hist : period for the historical data. has to be the same duration as the climate models period
#
#       input_folder : name of the output folder of climate_models_RCP_sfVersion.R
#       output_folder : name of the input folder of BiasCorrectPrecip.R and BiasCorrectTemperature.R
#
#   files :
#       output files of climate_models_RCP_sfVersion.R (under input_path):
#       temperature : scenario_climatemodel_daily_date_v4-NA_T.txt (scenario is one of the scenarii, length(date)<5)
#       precipitation : scenario_climatemodel_daily_date_v4-NA_RR.txt (see CAUTION for more details)
#
# OUTPUT :
#   for scenarii = ["hist","rcp45"], the output will be the following :
#       hist_prec.csv : precipitation historical values for the period startY_hist-endY_hist
#       hist_temp.csv : temperature historical values for the period startY_hist-endY_hist
#       rcp45_prec.csv : precipitation future values for the period startY_models-endY_models
#       rcp45_temp.csv : temperature future values for the period startY_models-endY_models
#       the first line is the header ('Date','scenario1','scenario2',...)
#       the first column is the date, every other column is a climate model.
#
# CAUTION :
#       will work if the files from climate_models_RCP_sfVersion.R have names like :
#       hist_climatemodelname_RR_daily_****.txt
#       rcp45_climatemodelname_TM_daily_****.txt
#       rcp85_climatemodelname_RR_daily_****.txt
#       **** usually looks like '2006_v-NA_RR' or '20-NA_T'
#
# ------------------------------------ SETUP -----------------------------------
catchment = 'Myglevatn'
geopoint_index = 184
scenarii = ["hist","rcp45","rcp85"]
startY_models, endY_models = 2041, 2070
startY_hist, endY_hist = 1971, 2000

input_folder = 'Output_RCP'
output_folder ='Input_BiasCorrect' # the input of BiasCorrect
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

import csv
from collections import defaultdict
from datetime import datetime, timedelta

#create a new repository for the catchment if not created yet
input_path = f'{klima}{input_folder}/{catchment}/'
output_path = f'{klima}{output_folder}/{catchment}/'

if not(os.path.exists(output_path)):
    os.mkdir(output_path)

# if the names of the input files have changed, you can modify it in consequence below
#input_path = 'D:/Klima/Output_RCP' #should end by '/Klima/Output_RCP'
#output_path = 'D:/Klima/Input_BiasCorrect'

# print basic info
print("Processing the data of the following catchment : ", catchment)
print('data coming from : ', input_path)
print(f'directory for setup creation : ', output_path)
print("")

# The historical and models data are being put together so they need to have the same duration
assert(endY_models - startY_models == endY_hist - startY_hist)

# basic function to take the date and the data from the interesting geographical point
def extract_dates_data_from_txt(file_path, geopoint_index_, scenarii_):
    '''input : .txt file computed by climate_models_RCP_sfVersion.R
       output : 2 lists, one with the dates the other with the corresponding data'''
    # Open the text file
    with open(file_path, 'r') as file:
        # Read the file line by line
        lines = file.readlines()
        # Separate each line into a list of strings using the space separator
        separated_lines = [line.split() for line in lines]

    tabular = np.array(separated_lines)
    # Pattern to match the name of the climate data
    patterns = [rf'{scenario}_(.*?)_daily_' for scenario in scenarii_]

    # Extract the model if its pattern is found
    for pattern in patterns :
        # Search for the pattern in the file_path
        match = re.search(pattern, file_path)
        if match:
            model = match.group(1)[:-3]

    dates = []
    dates.append('Date')
    data_list = []
    # Write the name of the models avoiding '_', '.' and '-' to get shorter names
    # (if the path for the path of the PINE parameter file is too long, PINE will crash)
    data_list.append(model.replace('_','').replace('.','').replace('-','')[:min(15,len(model))])

    for i in range(len(tabular[3:,0])):
        date = tabular[i+3,0]
        data = float(tabular[i+3,geopoint_index_])
        dates.append(datetime.strptime(date.replace("\"", ""), '%Y-%m-%d').date()) # the dates have quotation marks that we need to remove
        data_list.append(data)
        #print(i, data, date, len(dates))
    assert(len(dates) == len(data_list))
    return dates, data_list

import re
import numpy as np

def find_data_files(input_path_, P_T, scenario, catchment_):
    '''input  : the folder path with folders named catchment_scenario containing .txt data files generated by climate_models_RCP_sfVersion.R
                P_T = 'P' or 'T' in function of the desired data : prec or temp
                scenario = 'scenario_name'
       output : list of file names containing the temp or prec values for the grided catchment and its lenght
    '''
    input_path_s = input_path_ + catchment_ + '_' + scenario
    file_list = [f for f in os.listdir(input_path_s) if os.path.isfile(os.path.join(input_path_s, f))]
    if P_T == "P":
        pattern_P = re.compile(r'.*-NA_RR\.txt$')
        filtered_list_ = [filename for filename in file_list if pattern_P.search(filename)]
    elif P_T == "T":
        pattern_T = re.compile(r'.*-NA_T\.txt$')
        filtered_list_ = [filename for filename in file_list if pattern_T.search(filename)]
    nb_files = len(filtered_list_)
    #sort it because BiasCorrect needs the exact same order
    sorted_list = sorted(filtered_list_)
    return sorted_list, nb_files

def build_csv_loop(filtered_list_, input_path_, P_T, scenario, geopoint_index_):
    ''' create the temp or prec input file for BiasCorrect'''
    output_tab = [0 for i in range(len(filtered_list_) + 1)]

    for i in range(len(filtered_list_)):
        if scenario in filtered_list_[i]:
            file_path = input_path_ + catchment + '_' + scenario + '/' + filtered_list_[i]
            dates, data = extract_dates_data_from_txt(file_path, geopoint_index_, scenarii)
            output_tab[0] = dates
            output_tab[i+1] = data
    output_tab = np.array(output_tab)
    output_tab = np.transpose(output_tab)
    if filtered_list_[0][:4] in scenario and P_T == "P":
        output_file_path = output_path + scenario + '_prec.csv'
        np.savetxt(output_file_path, output_tab, delimiter=',', fmt='%s')
        print('csv file created successfully : ', output_file_path)
    elif filtered_list_[0][:4] in scenario and P_T == "T":
        output_file_path = output_path + scenario + '_temp.csv'
        np.savetxt(output_file_path, output_tab, delimiter=',', fmt='%s')
        print('csv file created successfully : ', output_file_path)
    else :
        print("file name not matched")
        output_file_path = None
    return output_tab

# Use build_csv_loop() to create all the input of BiasCorrect for the climate models
for scenario in scenarii:
    print('scenario : ', scenario)
    # filter the data files for both prec and temp and for a given scenario
    filtered_list_P, nb_files_P = find_data_files(input_path,"P",scenario, catchment)
    filtered_list_T, nb_files_T = find_data_files(input_path,"T", scenario, catchment)
    assert(nb_files_P == nb_files_T)
    print('number of input files : ', nb_files_P, 'for both temperature and precipitation')
    current_output_tabP = build_csv_loop(filtered_list_P, input_path, 'P', scenario, geopoint_index)
    current_output_tabT = build_csv_loop(filtered_list_T, input_path, 'T', scenario, geopoint_index)
    print("")














