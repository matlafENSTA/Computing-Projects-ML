# FINAL PLOTS OF AVERAGE YEARS FOR FUTURE PERIODS AND COMPARISON WITH HISTORICAL VALUES
#
# PURPOSE :
#     Extract the average snowpack and simulated discharge from PINE output files and regroup these data sets 
#     by scenario and period to create the input files for Plotting.R and AnnualMinFlows.R
#
# INPUT  :
#   variables :
#     catchment : name of the studied catchment. Be careful with capital letters or such
#     input_dir : name of the folder (not the path !) containing the output text files computed by PINE for all 
#                 climate models. It should be under catchment_path and named 'output_corr' or 'output_delta'.
#     output_dir : name of the directory where will be stored the concatenated results from PINE simulations for every climate model.
#     baseline : name of the corrected historical text output file of PINEHBV (has to be located in catchment_path)
#     end_identifier : the end of the output file names that all the results of the PINEHBV simulations have in 
#                 common (should be "_out.txt$")
#     corr_or_delta : string to make the difference between bias corrected and delta changed data sets.
#                     Should be consistent with input_dir directory. Usually "corrected" or "deltachanged".
#     snowpack_idd : name of the column with the daily sum of snowpacks for the catchment in the .txt PINE output file.
#     simdschrg_idd : name of the column with the daily data of simulated discharge for the catchment in the .txt PINE output file
#     date_format : string. Usually '%d.%m.%Y'.
#     len_catch_idd : integer. length of the string that identifies the catchment in the PINE output file names. 
#     Usually the 5 first letters. Example : the identifier of 'catchmentexample' is 'catch'.
#
#   files :
#     the pine text output files that are generated after you've done a simulation for every climate model and scenario
#     (have to be under catchment_path/input_dir), and for the baseline (under catchment_path).
#
# OUTPUT :
#   files : they will be created under pine_output_path/catchment.
#     catchment corr_or_delta discharge averageyear scenario period.png : daily average for discharge over a year
#         for a given future period and scenario. Comparison with historical values (Baseline).
#     catchment corr_or_delta discharge averageyear scenario period.csv : associated data file.
#     catchment corr_or_delta snowpack averageyear scenario period.png :  daily average for snowpacks over a year
#         for a given future period and scenario. Comparison with historical values (Baseline).
#     catchment corr_or_delta snowpack averageyear scenario period.csv : associated data file.
#
# WARNING : 
#     if there are already files in the output directory that have the same names as the output files of this script, 
#     they might be overwritten.
#     if the PINE output files are not named as following, this program may not work :
#         catchmentidd_corrordelta_model_scenario_period_out.txt
#         catchmentidd, corrordelta and scenario have 5 characters, period has 9 characters. if not, output files won't have the good names.
#     if you change the name of the dataframes ('snowpack_deltachanged_rcp45_2041_2070', simdschrg_deltachanged_rcp85_2071_2100'...),
#     created in DATA EXTRACTION,you have to adapt average_func() (reading the date and the scenario).
#
# NB :
#     if you are studying both bias corrected and delta changed data sets, you'll have to run this code twice.
#
# Load necessary library
library(ggplot2)
library(lubridate)
library(readxl)
library(scales)
library(showtext)
library(thematic)
library(ggokabeito)
library(janitor)
library(tidyr) # for the use of pivot_longer
library(sfsmisc) # to integrate the volume of water
library(scales) # for scientific notation formatting
library(dplyr)

# --------------------------------------------- SETUP ---------------------------------
catchment <- 'catchmentexample'
input_dir <- 'output_corr'
output_dir <- 'Input_plots'
baseline <- 'catchmentexample_obs_out.txt'
end_identifier <- "_out.txt$" # to identify the output text files of PINE in input_dir.

corr_or_delta = "corrected" # match with input_dir !
snowpack_idd <- 'FLDSNWPCK'
simdschrg_idd <- 'SIMDSCHRG'
date_format <- '%d.%m.%Y'
len_catch_idd <- 5
# --------------------------------------------- SETUP ---------------------------------

# Get the directory where the script is, then go back up to klima.
script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
# Delete all characters after the last '/'
klima <- sub("/[^/]*$", "", script_dir)
setwd(klima)
# if the automatic directory setting doesn't work, just precise the working directory below :
# setwd()
print(paste("wd: ", getwd()))

# Define the folder path
catchment_path <- paste(klima, 'PINE/PineProj', catchment, sep='/')
pine_output_path <- paste(catchment_path, input_dir, sep='/')
output_path <- paste(klima, output_dir, catchment, sep='/')
#print(paste('files coming from ', pine_output_path))
#print(paste('plots created in ', output_path))

if (!file.exists(output_path)) {
  dir.create(output_path)
  print(paste("Output directory created successfully for",catchment))
}

# --------------------------------------------- DATA EXTRACTION -----------------------

# List all files in the input folder
pine_output_all <- list.files(pine_output_path)

# Filter files ending with end_identifier and extract the different periods and scenarii
out_files <- pine_output_all[grep(end_identifier, pine_output_all)]

# Find automatically the computed scenarii and periods. only works if the name of the files have the good format
# (see WARNING above)
scenarii <- unique(substr(out_files, nchar(out_files) - 22, nchar(out_files) - 18))
periods <- unique(substr(out_files, nchar(out_files) - 16, nchar(out_files) - 8))
# Remove periods before year 2000
periods <- periods[as.integer(substr(periods, 1, 4)) >= 2020]
# if the lines above don't work, enter the lists of periods and scenarii manually :
# periods <- c("2041-2070","2071-2100")
# scenarii <- c("histo", "rcp45", "rcp85")

# Check if there are any files matching the pattern
if(length(out_files) == 0) {
  stop(paste("No files ending with",end_identifier,"found in the folder. Try another pattern"))
}

#' INPUT :
#'    strings : list of strings
#'    patterns : list of strings (patterns that can be found in the elements of strings)
#' OUTPUT :
#'    for a given list of strings, return a list of lists where the nth sublist groups all the 
#' elements of strings matching patterns[i]
group_strings_by_patterns <- function(strings, patterns) {
  # Initialize a list to hold sublists for each pattern
  result <- vector("list", length(patterns))
  names(result) <- patterns
  
  # Loop through each pattern
  for (i in seq_along(patterns)) {
    pattern <- patterns[i]
    
    # Find strings that match the current pattern
    matching_strings <- strings[grepl(pattern, strings, ignore.case = TRUE)]
    
    # Add matching strings to the corresponding sublist
    result[[pattern]] <- matching_strings
  }
  return(result)
}

# make sublists of out_files for each period that has been studied (near and far future for example)
out_by_period <- group_strings_by_patterns(out_files, periods)
out_by_scenario <- vector("list", length(periods))
for (i in seq_along(out_by_period)){
  out_by_scenario[[i]] <- group_strings_by_patterns(out_by_period[[i]], scenarii)
}

#' extract_model(input_string, catchment_name, end_identifier_str, period_str, len_catch_idd_)
#' knowing the name of the catchment, the period format and the end of the filename (=input_string), 
#' finds the name of the climate model of the input_string
extract_model <- function(input_string, catchment_name, end_identifier_str, period_str, scenario_str, len_catch_idd_) {
  # Append an underscore to the catchment_name and only keep the len_catch_idd_ first characters
  catchment_name <- paste0(substr(catchment_name,1,len_catch_idd_), "_")
  
  # Find the start index of catchment_name in the input_string
  start_index <- regexpr(catchment_name, input_string, ignore.case = FALSE)
  
  if (start_index == -1) {
    # If not found, convert input_string to lowercase and try again
    catchment_name <- tolower(catchment_name)
    start_index <- regexpr(catchment_name, input_string, ignore.case = TRUE)
    if (start_index == -1) {
      print(paste0("Substring ", catchment_name, " not found"))
      return("UNKNOWN_MODEL")
    }
  }
  
  # Calculate start and end indices for the substring to extract (period and scenario are written at the end of the file name)
  start_index <- start_index + nchar(catchment_name)
  end_index <- nchar(input_string) - nchar(end_identifier_str) - nchar(period_str) - nchar(scenario_str) - 1
  
  # Extract the substring
  extracted_str <- substr(input_string, start_index, end_index)
  
  return(extracted_str)
}

#' process_pine_output(input_path, file_name_str, catchment_str, end_identifier_str, period_str, date_format_str)
#' with the information provided, this program extracts the snowpack and the discharge values simulated by PINE.
process_pine_output <- function(input_path, file_name_str, catchment_str, end_identifier_str, period_str, scenario_str, date_format_str) {
  # Construct the full file path
  file_to_read <- file.path(input_path, file_name_str)
  
  # Check if the file exists before attempting to read it
  if (!file.exists(file_to_read)) {
    warning(paste("File does not exist:", file_to_read))
    return(NULL)
  }
  
  # Read the header line to get column names
  header_line <- readLines(file_to_read, n = 1)
  column_names <- unlist(strsplit(header_line, "\\s+"))
  # Read the data while skipping the first two header lines
  data <- read.table(file_to_read, header = FALSE, skip = 2, col.names = column_names, fill = TRUE)
  # Remove 'X' from column names
  colnames(data) <- gsub("^X", "", colnames(data))
  
  if (inherits(data, "try-error")) {
    warning(paste("Failed to read file:", file_to_read))
    return(NULL)
  }
  
  # Remove the 'X' and the '1' prefixes from column names
  names(data) <- gsub("^X", "", names(data))
  names(data) <- gsub("^1", "", names(data))
  
  # Check if the required columns exist
  # should match the names of the columns of pine's output files
  if (!all(c("date", snowpack_idd, simdschrg_idd) %in% names(data))) {
    warning(paste("Required columns missing in file :", file_to_read))
    return(NULL)
  }
  
  # Extract the date, the target columns and the units
  date_column <- as.Date(data$date, date_format_str)
  snwpcks_column <- data[[snowpack_idd]]
  simdschrg_column <- data[[simdschrg_idd]]
  snwpck_unit <- snwpcks_column[1]
  simdschrg_unit <- simdschrg_column[1]
  
  # Extract the model name from the file name
  model_name <- extract_model(file_name_str, catchment_str, end_identifier_str, period_str, scenario_str, len_catch_idd)
  
  # Create data frames with date and the target columns, setting dynamic column names
  snwpcks_df <- data.frame(Date = date_column, snwpcks_column)
  colnames(snwpcks_df)[2] <- model_name
  
  simdschrg_df <- data.frame(Date = date_column, simdschrg_column)
  colnames(simdschrg_df)[2] <- model_name
  
  return(list(snowpack = snwpcks_df, simdschrg = simdschrg_df, model_name = model_name, snowpack_unit = snwpck_unit, discharge_unit = simdschrg_unit))
}

# Baseline processing
baseline_data <- process_pine_output(catchment_path, baseline, catchment, end_identifier, periods[1], scenarii[1], date_format)
snowpack_baseline <- baseline_data$snowpack
colnames(snowpack_baseline)[2] <- 'Baseline'
simdschrg_baseline <- baseline_data$simdschrg
colnames(simdschrg_baseline)[2] <- 'Baseline'

# Initialize lists to store the names of the data frames
snowpack_all <- list()
simdschrg_all <- list()

# Iterate through each period and scenario to extract simdschrg and snwpck from pine output files
for (p in seq_along(out_by_period)) {
  current_period <- periods[p]
  for (s in seq_along(out_by_scenario[[p]])){
    current_files <- out_by_scenario[[p]][[s]]
    current_scenario <- scenarii[[s]]
    
    # Filter the non-empty files into current_filled_files
    current_filled_files <- c()
    # Iterate through each file
    for (file_name in current_files) {
      # Check if the file is empty
      if (file.size(paste(pine_output_path,file_name,sep='/')) == 0) {
        cat(sprintf("File %s is empty. Skipping...\n", file_name))
        next
      }
      # Append non-empty file path to the list
      current_filled_files <- c(current_filled_files, file_name)
    }
    if (length(current_filled_files) == 0){ # all pine simulations concerning this scenario and this period are empty
      next
    }
    
    # Initialize empty lists to store the data frames and dates
    snwpcks_data <- list()
    simdschrg_data <- list()
    
    # Processing of the first file to set the correct dates
    first_file <- process_pine_output(pine_output_path, current_filled_files[1], catchment, end_identifier, current_period, current_scenario, date_format)
    
    # Convert date list to character vector
    snwpcks_date_char <- as.character(first_file[['snowpack']][['Date']])
    
    # Replace NA with format
    snwpcks_date_char[is.na(snwpcks_date_char)] <- "%Y-%m-%d" #useless
    snwpcks_data[['Date']] <- as.character(snwpcks_date_char)
    simdschrg_data[['Date']] <- as.character(first_file[['simdschrg']][['Date']])
    snwpcks_data[['Baseline']] <- snowpack_baseline[['Baseline']]
    simdschrg_data[['Baseline']] <- simdschrg_baseline[['Baseline']]
    
    # Iterate through each climate model and scenario for the given period
    for (file_name in current_filled_files) {
      results <- process_pine_output(pine_output_path, file_name, catchment, end_identifier, current_period, current_scenario, date_format)
      current_model <- results$model_name
      snowpack_unit <- results$snowpack_unit
      discharge_unit <- results$discharge_unit
      
      # Add the data frames to the respective lists (excluding date columns)
      snwpcks_data[[current_model]] <- results$snowpack[, -1]
      simdschrg_data[[current_model]] <- results$simdschrg[, -1]
    }
    
    # Combine all data frames for snowpack and simulated discharge
    combined_snwpcks <- paste("snowpack",corr_or_delta, current_scenario, gsub("-", "_", current_period), sep = "_")
    combined_simdschrg <- paste("simdschrg",corr_or_delta, current_scenario, gsub("-", "_", current_period), sep = "_")
    print(combined_snwpcks)
    
    # Create combined data frames
    snwpcks_combined <- do.call(cbind, snwpcks_data)[-1,]
    simdschrg_combined <- do.call(cbind, simdschrg_data)[-1,]
    
    # Assign final data frames with correct date columns and variable names
    assign(combined_snwpcks, snwpcks_combined)
    assign(combined_simdschrg, simdschrg_combined)
    
    # Register the names of combined data frames in lists
    snowpack_all <- c(snowpack_all, combined_snwpcks)
    simdschrg_all <- c(simdschrg_all, combined_simdschrg)
    
    # Save the dataframes to text files
    snwpck_filepath = paste0(output_path,'/',catchment,'_daily_',combined_snwpcks,".txt")
    simdschrg_filepath = paste0(output_path,'/',catchment,'_daily_',combined_simdschrg,".txt")
    write.table(snwpcks_combined, file=snwpck_filepath, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
    write.table(simdschrg_combined, file=simdschrg_filepath, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
  }
}

print("data extracted successfully")
