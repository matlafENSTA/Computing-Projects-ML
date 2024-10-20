# CHECK THE FORMAT OF AN INPUT TEXT FILE FOR PINE
#
# PURPOSE :
#     PINEHBV requires some data to work but the input data files are really specific. If the conversion of your
#     .txt text input file into a .dat binary file doesn't work (PINEHBV>Import data>TXT2PINE), this program could help you find the problem.
#
# INPUT :
#     file_path : string. The path to the text file you want to check. This should be a valid path to a .txt file.
#                 Example: 'C:/Klima/PINE/PineProj/catchmentexample/catchmentexample_obs.txt'
#     input_path : if you want to check a lot of pine input files at the same time, put them at input_path, ensure there are the only text files 
#                 in this repository and run the script. If you just want to check one file, enter input_path = NULL.
#     max_value : numeric. The maximum allowed value for numeric columns in the data file. Values exceeding this threshold will trigger a warning.
#                 Example: 100
#     min_value : numeric. The minimum allowed value for numeric columns in the data file. Values below this threshold will trigger a warning.
#                 Example: -30
#     date_format : string. The expected format for date values in the data file. Should be in the format 'dd.mm.yyyy'.
#     time_format : string. The expected format for time values in the data file. Should be in the format 'hh.mm.ss'.
#     date_pattern : string. A regular expression pattern to validate the date format, ensuring day, month, and year are within valid ranges.
#                 Usually : "^([0-2][0-9]|3[0-1])\\.([0][1-9]|1[0-2])\\.([1][9][0-9][0-9]|[2][0-9][0-9][0-9])$"
#     time_pattern : string. A regular expression pattern to validate the time format, ensuring hours, minutes, and seconds are within valid ranges.
#                 Usually : "^(0[0-9]|1[0-9]|2[0-4])\\.(0[0-9]|[1-5][0-9]|60)\\.(0[0-9]|[1-5][0-9]|60)$"
#
# OUTPUT :
#     validity : 'VALID' if the file can be used as a PINE input and be converted into a .dat data file. The program will stop and raise an error
#     if any issues are detected.
#
# NB : 
#     Do not forget to replace the '/' by '/' in the path when specifying the file path.
#     Example of a valid PINE input file: 'C:/Klima/PINE/PineProj/catchmentexample/catchmentexample.txt'
#     This script is designed for PINEHBV Version 1.0

library(lubridate)

# --------------------------------------- SETUP ---------------------------------------
file_path <- "D:/Klima/PINE/PineProj/Myglevatn/Myglevatn_obs.txt"
input_path <- "D:/Klima/PINE/PineProj/Myglevatn/input_delta"

max_value <- 200
min_value <- -30

# do not modify this if you are running the Version 1.0 of PINEHBV :
date_format <- "dd.mm.yyyy"
time_format <- "hh.mm.ss"
# Date pattern: dd.mm.yyyy (where day is 01-31, month is 01-12, and year is > 1900)
date_pattern <- "^(0[1-9]|[12][0-9]|3[01])\\.(0[1-9]|1[0-2])\\.(19[0-9][0-9]|20[0-9][0-9])$"
# Time pattern: hh.mm.ss (where hh is between 00 and 24, mm and ss are between 00 and 60)
time_pattern <- "^(0[0-9]|1[0-9]|2[0-4])\\.(0[0-9]|[1-5][0-9]|60)\\.(0[0-9]|[1-5][0-9]|60)$"
# --------------------------------------- SETUP ---------------------------------------

checking <- function(file_path){
  print(paste("Checking", file_path))
  # Ensure the file is a .txt file
  len <- nchar(file_path)
  if (substr(file_path, len - 3, len) != '.txt') {
    stop('This is not a text file. If it is a .csv file, convert it just by changing the extension to .txt')
  }
  
  # Read the tabular data
  raw_tabular <- read.table(file_path, stringsAsFactors = FALSE)
  
  # ---------- STEP1: Check column names
  acceptable_date_names <- c("dato", "date")
  acceptable_time_names <- c("time", "tid")
  
  header_names <- as.character(raw_tabular[1, ])
  if (!tolower(header_names[1]) %in% acceptable_date_names) {
    warning("The first column should have a name like 'date' or 'Dato'")
  }
  if (!tolower(header_names[2]) %in% acceptable_time_names) {
    warning("The second column should have a name like 'Time' or 'tid'.")
  }
  
  # Extract the first row and store it in 'header', then remove it
  header <- raw_tabular[1, ]
  without_header <- raw_tabular[-1, ]
  
  # ---------- STEP2: Check the units
  # Extract the second row and store it in 'units', then remove it
  units <- without_header[1, ]
  without_units <- without_header[-1, ]
  if (units$V1 != date_format) {
    warning(paste("The format of the date unit is wrong, it should be", date_format, "not", units$V1))
  }
  if (units$V2 != time_format) {
    warning(paste("The format of the time unit is wrong, it should be", time_format, "not", units$V2))
  }
  
  # ---------- STEP3: Convert and check date format for every date
  valid_dates <- grepl(date_pattern, without_units$V1)
  
  if (!all(valid_dates)) {
    invalid_dates <- which(!valid_dates)
    stop(paste("Invalid date format found in date column at row(s):", paste(invalid_dates, collapse = ", ")))
  }
  
  # ---------- STEP4: Convert and check time format for every time
  valid_times <- grepl(time_pattern, without_units$V2)
  
  if (!all(valid_times)) {
    invalid_times <- which(!valid_times)
    stop(paste("Invalid time format found in time column at row(s):", paste(invalid_times, collapse = ", ")))
  }
  
  # ---------- STEP5: Check if data columns contain only numeric values
  only_data <- without_units[, !(names(without_units) %in% c("V1", "V2"))]
  
  for (col in names(only_data)) {
    # Convert the column to numeric, handling character representations of numbers
    numeric_col <- as.numeric(as.character(only_data[[col]]))
    
    # Check for NA values that indicate non-numeric entries
    if (any(is.na(numeric_col))) {
      invalid_rows <- which(is.na(numeric_col))
      stop(paste("Column", col, "contains non-numeric values at row(s):", paste(invalid_rows, collapse = ", ")))
    }
  }
  
  # ---------- STEP6: Check for values greater than max_value and less than min_value
  for (i in seq_along(names(only_data))) {
    col <- names(only_data)[i]
    numeric_col <- as.numeric(as.character(only_data[[col]]))
    
    # Check for values greater than max_value
    if (any(numeric_col > max_value)) {
      invalid_rows_max <- which(numeric_col > max_value)
      warning(paste("Column", header_names[i + 2], "contains values greater than", max_value, "at row(s):", paste(invalid_rows_max, collapse = ", ")))
    }
    
    # Check for values less than min_value
    if (any(numeric_col < min_value)) {
      invalid_rows_min <- which(numeric_col < min_value)
      warning(paste("Column", header_names[i + 2], "contains values less than", min_value, "at row(s):", paste(invalid_rows_min, collapse = ", ")))
    }
    
    # Report max and min for each column
    max_val <- max(numeric_col, na.rm = TRUE)
    min_val <- min(numeric_col, na.rm = TRUE)
    
    cat(sprintf("Column '%s': Max = %.2f, Min = %.2f\n", header_names[i + 2], max_val, min_val))
  }
  
  # If all checks pass
  print(paste("All checks passed. This PINE input file is valid :", file_path))
  return(raw_tabular)
}

if (!is.null(file_to_check)){
  raw_data <- checking(file_to_check)
}

if (!is.null(input_path)){
  # List all the files stored at input_path
  input_all <- list.files(input_path)
  # Filter the list to include only text files
  filtered_input <- grep(".txt", input_all, value = TRUE)
  for (current_pine_input in filtered_input){
    raw_data <- checking(paste(input_path,current_pine_input, sep='/'))
  }
  print(paste("every pine input file found at", input_path, "is valid"))
}
