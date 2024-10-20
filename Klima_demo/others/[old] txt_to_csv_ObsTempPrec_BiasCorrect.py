# USELESS

# this program has been implemented to run under python 3.12. If not working, check first your python version.
# make sure all the python modules used are installed on your computer.
# PURPOSE : format the data computed by csv_to_txt_Sildre_Norge.py for the temperature and precipitation before running BiasCorrectTemperature.R and BiasCorrectPrecip.R
# The main work here is to split the data into two .csv files (one for prec, the other for temp) and change the format of the date from %d.%m.%Y to %Y-%m-%d
#why not to just take the files downloaded at Sildre or SeNorge ? Because data can be missing and this problem is fixed by csv_to_txt_Sildre_Norge.py
### How to use
# modify catchment name, input_path and output_path according to your device

# ------------------------------------ SETUP -----------------------------------
catchment = 'Myglevatn'
#input_path = "C:/PINE/PineProj"
#output_path = 'D:/Klima/Input_BiasCorrect'
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
script_dir = os.path.dirname(absolute_script_path) + '/'
script_dir = script_dir.replace('\\','/')

#create a new repository for the catchment if not created yet

input_path = script_dir + 'PINE/PineProj/'
output_path = script_dir + 'Input_BiasCorrect/' + catchment + '/'
if not(os.path.exists(output_path)):
    os.mkdir(output_path)

# print basic info
print("Processing the data of the following catchment : ", catchment)
print('data coming from : ', input_path)
print(f'directory for setup creation : ', output_path)
print("")

import csv
from collections import defaultdict
from datetime import datetime, timedelta
import numpy as np

def extract_dates_data_from_txt(file_path):
    '''input : .txt file computed by csv_to_txt_Sildre_SseNorge.py
       output : 3 lists, one with the dates the others with the corresponding precipitation and temperature values'''
    # Open the text file
    with open(file_path, 'r') as file:
        # Read the file line by line
        lines = file.readlines()
        # Separate each line into a list of strings using the space separator
        separated_lines = [line.split() for line in lines]

    tabular = np.array(separated_lines)
    # Pattern to match the name of the climate data

    dates = []
    prec_list = []
    temp_list = []
    units = tabular[1,:]#to be printed

    for i in range(len(tabular[2:,0])):
        date = tabular[i+2,0]
        prec = float(tabular[i+2,2])
        temp = (float(tabular[i+2,3])+float(tabular[i+2,4]))*0.5
        dates.append(datetime.strptime(date, '%d.%m.%Y').date())
        prec_list.append(prec)
        temp_list.append(temp)
        #print(i, data, date, len(dates))
    assert(len(dates) == len(prec_list) and len(dates) == len(temp_list))
    return dates, prec_list, temp_list

dates, prec_list, temp_list = extract_dates_data_from_txt(input_path + catchment + '/' + catchment + '.txt')

def build_csv(dates_, data, output_path_, P_T):
    output_tab = [0,0]
    output_tab[0] = dates
    output_tab[1] = data
    output_tab[0].insert(0,'Date')
    output_tab[1].insert(0,'Obs')
    # Transpose the table
    output_tab = list(map(list, zip(*output_tab)))
    output_tab = np.array(output_tab)
    if P_T == "P":
        output_file_path = output_path_ + 'observed_prec.csv'
        np.savetxt(output_file_path, output_tab, delimiter=',', fmt='%s')
    elif P_T == "T":
        output_file_path = output_path_ + 'observed_temp.csv'
        np.savetxt(output_file_path, output_tab, delimiter=',', fmt='%s')
    else :
        print("file name not matched")
        output_file_path = None
    print('csv file created successfully : ', output_file_path)
    return None

build_csv(dates, temp_list, output_path, 'T')
build_csv(dates, prec_list, output_path, 'P')




