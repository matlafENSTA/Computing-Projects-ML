# this program has been implemented to run under python 3.12. If not working, first check your python version.
# make sure all the python modules used are installed on your computer.
#
# PURPOSE :
#       automatise the creation of PINE input data file to create the setup for observed values and calibrate the model.
#       you will find more information in PINE_tutorial.pdf
#
# INPUT :
#   variables :
#       catchment : name of the studied catchment (and PINE/PineProj repository). Be careful with upper case letters.
#       input_dirname : name of the subfolder you put the downloaded files in. Has to be under Klima/PINE/PineProj/catchment
#       input_bias_correct : name of the input of biascorrect, it will create the observed input files in it
#                            (in a subfolder of the name of the catchment).
#
#       extreme_value : if the absolute value of some data is superior to extreme_value, it will be considerated it as
#                       an error and will be corrected linearly. 500 is a good value.
#       overview_len : (optional) the overview_len first values extracted from the files will be displayed on the shell for checking.
#
#   files : have to be under input_path (run the beginning of this script to display it)
#       filename_discharge : name of the discharge data file coming from Sildre
#       filename_temperature : name of the temperature data file coming from SeNorge
#       filename_precipitation : name of the precipitation data file coming from SeNorge
#
# OUTPUT :
#       txt_tabular : text file containing the data to create a setup with PINE.
#       corrected_discharge_path : the data corrected to fit the data type taken by PINE : daily data without missing days or
#                                  nonsensical values.
#       corrected_temperature_path : same for temperature
#       corrected_precipitation_path : same for precipitation
#       input_bias_correct + "observed_temp.csv" : copy of corrected_temperature_path for the input of BiasCorrectTemperature.R
#       input_bias_correct + "observed_prec.csv" : copy of corrected_precipitation_path for the input of BiasCorrectPrecip.R
#
# CAUTION :
#       if there is a problem of delimiter, open each file with Excel > register as > same name, but with the delimiter normally
#       used by your computer (usually a comma).
#
# NB :
#       extreme_value is the same for temperature, precipitation and discharge. feel free to implement different extreme
#       values for every data type (precipitation and discharge cannot be negative for example).
#
# ------------------------------------ SETUP -----------------------------------
catchment = "Myglevatn"
input_dirname = 'input_obs'
input_bias_correct_dirname = 'Input_BiasCorrect'

filename_discharge = catchment + "_downloaded_Discharge.csv"
filename_temperature = catchment + "_downloaded_Temperature.csv"
filename_precipitation = catchment + "_downloaded_Precipitation.csv"

extreme_value = 500
overview_len = 20
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

input_path = f'{klima}PINE/PineProj/{catchment}/{input_dirname}/'
file_path_discharge = input_path + filename_discharge
file_path_temperature = input_path + filename_temperature
file_path_precipitation = input_path + filename_precipitation
# path to the pine input file that will be created by this script :
txt_tabular = f'{klima}PINE/PineProj/{catchment}/{catchment}_obs.txt'

#creating a subfolder of the name of the catchment in the input folder of BiasCorrect
input_bias_correct = f'{klima}/{input_bias_correct_dirname}/{catchment}/'

def create_directory(directory_path):
    if not os.path.exists(directory_path):
        os.makedirs(directory_path)
        print(f"Directory '{directory_path}' created.")
    else:
        print(f"Directory '{directory_path}' already exists.")
    print("")

create_directory(input_bias_correct)

# print basic info
print("Processing the data of the following catchment : ", catchment)
print("data coming from : ", input_path)
print("")

# Extraction of data from the downloaded files
import csv
from datetime import datetime, timedelta
from colorama import Fore, Style # just to print green success comments
import warnings

def format_date(date_string):
    ''' find the format of a date string trying a few different date format patterns'''
    try :
        date_YYYYmmdd = datetime.strptime(date_string, '%Y-%m-%d %H:%M')
    except ValueError:
        try :
            date_YYYYmmdd = datetime.strptime(date_string, '%Y.%m.%d %H:%M')
        except ValueError:
            try :
                date_YYYYmmdd = datetime.strptime(date_string, '%d-%m-%Y %H:%M')
            except ValueError :
                date_YYYYmmdd = datetime.strptime(date_string, '%d.%m.%Y %H:%M')
    return date_YYYYmmdd

def find_altitude(row_list):
    for elem in row_list :
        if elem.find("Altitude: ") == 0:
            altitude = int(elem[10:])
            return altitude
    return "_"

def extract_data_from_csv(file_path_, extreme_value_, catchment_):
    '''just extract the data as it is AND add the missing dates'''
    with open(file_path_, mode='r') as csvfile:
        # Detect the delimiter used in the file
        sample = csvfile.read(1024)
        csvfile.seek(0)
        dialect = csv.Sniffer().sniff(sample)

        # Create a CSV reader object with the detected dialect
        csvreader = csv.reader(csvfile, dialect)

        # Read and ignore the title row
        title_row = next(csvreader)
        altitude = find_altitude(title_row)
        # Read and skip the header row
        header = next(csvreader)

        # Read the first data row to initialize `prev_time_value`
        first_row = next(csvreader)
        prev_time_value = format_date(first_row[0][:16])  # Get only the date+hour part (YYYY-MM-DD HH-MM)
        try:
            prev_data_value = float(first_row[1])
        except ValueError:
            prev_data_value = 1.0

        # Initialize variables for daily data aggregation
        time = []
        time.append(prev_time_value)
        data = []
        data.append(prev_data_value)

        # Loop through each row in the CSV file
        count_measures = 0 # all the measures should have = time interval * acquisition frequence
        count_missing_data = 0 # missing measurement
        count_missing_days = 0 # missing days in the document (no date, no data)
        r = 1
        while True:
            try:
                row = next(csvreader)
            except StopIteration:
                break

            # Extract the date+hour part (first 16 characters) and convert it to date format
            current_date = format_date(row[0][:16])

            if current_date - time[-1] > timedelta(1):
                # if the date is missing, we increment it anyway, set the data value to "" (fix that later using linear)
                missing_day = time[-1]
                while missing_day < current_date:
                    missing_day = missing_day + timedelta(1)
                    time.append(missing_day)
                    data.append("")
                    count_missing_days += 1
                    count_missing_data += 1
                    count_measures += 1
            else : # we have the data of the following day
                time.append(current_date)
                current_data = row[1]
                if row[1] == "":
                    data.append("")
                    count_missing_data += 1
                else :
                    float_data = float(current_data)
                    # for discharge, temperature and precipitation cannot exceed 500 : when there are such unexpected values, we treat it as a missing data.
                    if abs(float_data) > extreme_value_:
                        data.append("")
                        count_missing_data += 1
                    else :
                        data.append(float_data)

                count_measures += 1

            #Update indents
            prev_time_value = current_date
            r +=1

        first_day = min(time)
        last_day = max(time)
        assert(first_day == time[0] and last_day == time[-1])
        assert(len(time) == len(data))

        # Return the first day, last day, and the list of data values, sorted chronologically. Also print information about the missing values
        print("  - ", catchment, file_path_[file_path_.rindex('_') + 1:-4])
        print("time interval : ", first_day.date(), " - ", last_day.date())
        print("missing days : ", count_missing_days, "/", (last_day - first_day).days)
        percent_missing = 100*count_missing_data/count_measures
        print("missing data : ", count_missing_data, "/", count_measures, "being ", percent_missing, "% of empty values")
        if 5 < percent_missing:
            warnings.warn("!! A lot of data values are missing !!")
        print("first data values : ", data[:overview_len])
        return data, time, altitude

def linear(d1,y1,d2,y2):
    '''input : d1 < d2 two measurement time (datetime type) and y1, y2 their associated data
    output : a,b the coefficents of the linear function f(x) = ax + b that meets these two points'''
    #convert the times to float, the zero is set as reference_day
    x1 = 0
    x2 = (d2 - d1).total_seconds()/60
    a = (y2 - y1)/(x2 - x1)
    b = y1 - a*x1
    return a, b

def fill_data(data, time):
    '''fill the empty values linearly from the last known value to the next known value'''
    numeric_data = [value for value in data if isinstance(value, float)]
    maxi = max(numeric_data)
    mini = min(numeric_data)
    time_indent = 0
    end = len(time)
    count_missing_values = 0
    while time_indent < end:
        if data[time_indent] == "":
            indent_before_missing = time_indent - 1 #problem if the very first day has no value
            time_before_missing = time[indent_before_missing]
            data_before_missing = data[indent_before_missing]
            while data[time_indent] == "":
                time_indent +=1
                count_missing_values += 1
            a,b = linear(time_before_missing, data_before_missing, time[time_indent], data[time_indent])
            for m in range(indent_before_missing + 1, indent_before_missing + count_missing_values + 1):
                data[m] = a*(time[m] - time_before_missing).total_seconds()/60 + b
            count_missing_values = 0
        time_indent += 1
    # verify that there is no bad values (the maximum cannot be superior with the linear filling, mathematical property)
    assert(maxi <= max(data))
    assert(mini <= min(data))
    return data, mini, maxi

def daily_average_data(data, times):
    '''calculates the average of the data by day (useful for hourly data)'''
    if len(times) != len(data):
        raise ValueError("Both lists must have the same length")

    # Dictionary to store the sum and count of data per day
    data_per_day = {}

    # Iterate over the time and data lists
    for time, value in zip(times, data):
        date = time.date()
        if date not in data_per_day:
            data_per_day[date] = {'sum': 0, 'count': 0}
        data_per_day[date]['sum'] += value
        data_per_day[date]['count'] += 1

    # Prepare the results lists
    processed_times = []
    processed_data = []

    # Calculate the average for each day
    for date, data in sorted(data_per_day.items()):
        processed_times.append(date)
        processed_data.append(round(data['sum'] / data['count'], 6))

    return processed_data, processed_times

# for precipitation, we have to take the daily sum not the daily average
def daily_sum_data(data, times):
    '''calculates the sum of the data by day (useful for hourly data)'''
    if len(times) != len(data):
        raise ValueError("Both lists must have the same length")

    # Dictionary to store the sum and count of data per day
    data_per_day = {}

    # Iterate over the time and data lists
    for time, value in zip(times, data):
        date = time.date()
        if date not in data_per_day:
            data_per_day[date] = {'sum': 0}
        data_per_day[date]['sum'] += value

    # Prepare the results lists
    processed_times = []
    processed_data = []

    # Calculate the average for each day
    for date, data in sorted(data_per_day.items()):
        processed_times.append(date)
        processed_data.append(round(data['sum'], 6))

    return processed_data, processed_times

def process_data(file_path, sum_or_average):
    '''gather the work of the last functions : extract the data of the downloaded files, fill missing values linearly and remove unrealisitc values
    input : sum_or_average = 'sum' or 'avg' to precise if you want daily sum or daily average for data (precip : sum, discharge : avg)'''
    data, time, altitude_ = extract_data_from_csv(file_path, extreme_value, catchment)
    filled_data, mini, maxi = fill_data(data, time)
    print('extreme values before average or sum      : ', mini, maxi)
    if sum_or_average == 'sum':
        daily_data, timeline = daily_sum_data(filled_data,time)
    else :
        daily_data, timeline = daily_average_data(filled_data,time)
    print('extreme values after daily average or sum : ', min(daily_data), max(daily_data))
    print()
    return timeline, daily_data, altitude_

# Discharge
timeline_D, discharge_all, altitude_D = process_data(file_path_discharge, 'avg') #the altitude is not written in the files from Sildre so altitude_D will always be "_". The good data comes from SeNorge (temp and prec files). So do not worry if altitude_D = '_'.

# Temperature
timeline_T, temperature_all, altitude_T = process_data(file_path_temperature, 'avg')

# Precipitation
timeline_P, precipitation_all, altitude_P = process_data(file_path_precipitation, 'sum')

assert(altitude_P == altitude_T)
if altitude_P == "_":
    print("Altitude not found for this catchment.")
else :
    print("Altitude : ", altitude_P, "meters above sea level")

def save_corrected_dailydata(timeline, data, data_name, catchment_):
    ''' save the full observed daily value for a variable (for the setup creation it is being cut to match the period of each data set)
    input : the full data set for a variable and the associated timeline (datetime list). data_name and catchment_ are strings.
    what for ? to run BiasCorrect R function'''
    # Convert dates to string format "YYYY-MM-DD"
    timeline_str = [date.strftime("%Y-%m-%d") for date in timeline]

    # Define the filename (and path)
    filename = f'{input_path}{catchment_}_corrected_{data_name}.csv'

    # Write to CSV
    with open(filename, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["Date", "Obs"])
        writer.writerows(zip(timeline_str, data))
    return filename

corrected_discharge_path = save_corrected_dailydata(timeline_D, discharge_all, 'Discharge', catchment)
corrected_temperature_path = save_corrected_dailydata(timeline_T, temperature_all, 'Temperature', catchment)
corrected_precipitation_path = save_corrected_dailydata(timeline_P, precipitation_all, 'Precipitation', catchment)

import shutil

# copy the corrected files for temprerature and precipitation into the input folder of BiasCorrect R functions.
shutil.copy(corrected_temperature_path, input_bias_correct + "observed_temp.csv")
shutil.copy(corrected_precipitation_path, input_bias_correct + "observed_prec.csv")

print(f"File copied from {corrected_temperature_path} to {input_bias_correct}")

def create_pine_input_file(filename, timeline_D_, timeline_T_, timeline_P_, Dis, Temp, Prec, catchment_, altitude = ''):
    firstday_D, lastday_D = timeline_D_[0], timeline_D_[-1]
    firstday_T, lastday_T = timeline_T_[0], timeline_T_[-1]
    firstday_P, lastday_P = timeline_P_[0], timeline_P_[-1]
    first_day, last_day = max(firstday_D, firstday_T, firstday_P), min(lastday_D, lastday_T, lastday_P)
    if first_day > last_day:
        raise ValueError('There is no period in common for the data provided')
    duration = abs(first_day - last_day)  # Calculate duration based on first_day and last_day
    Dis = Dis[abs(first_day - firstday_D).days: abs(first_day - firstday_D).days + duration.days]
    Temp = Temp[abs(first_day - firstday_T).days: abs(first_day - firstday_T).days + duration.days]
    Prec = Prec[abs(first_day - firstday_P).days: abs(first_day - firstday_P).days + duration.days]
    #print(len(Dis), len(Temp), len(Prec), first_day.date(), last_day.date(), duration

    # Open a text file for writing
    with open(filename, 'w') as file:
        # Write headers
        file.write("Dato\tTime\tP_mm_{}moh\tT_minC_{}moh\tT_maxC_{}moh\tQ_{}\n".format(altitude, altitude, altitude, catchment_))
        file.write("dd.mm.yyyy\thh.mm.ss\tmm\tgrC\tgrC\tm3/s\n")

        # Generate dates between first_day and last_day
        current_day = first_day
        current_row = 0  # Start from the first row of data

        time_str = "00.00.00"
        while current_day < last_day:
            # Append date and an arbitrary hour to the corresponding columns
            date_str = current_day.strftime("%d.%m.%Y")
            # Retrieve data from lists
            prec_value = Prec[current_row]
            temp_min_value = Temp[current_row]
            temp_max_value = Temp[current_row]
            dis_value = Dis[current_row]

            # Write data to the file
            file.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(date_str, time_str, prec_value, temp_min_value, temp_max_value, dis_value))

            # Move to the next day
            current_day += timedelta(days=1)
            current_row += 1
    return None

create_pine_input_file(txt_tabular, timeline_D, timeline_T, timeline_P, discharge_all, temperature_all, precipitation_all, catchment, altitude_P)

print(Fore.GREEN + "data corrected and written successfully" + Style.RESET_ALL)