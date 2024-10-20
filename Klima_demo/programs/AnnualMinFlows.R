# ANNUAL MIN FLOWS
#
# PURPOSE :
#     Calculate the 7-days average flows for the input periods (every daily data is replaced by the average of the 7 next days)
#     Identify, for each year, the lowest flows and plot them for every climate model, the baseline and the average of the climate models
#
# INPUT : 
#   variables :
#     catchment : name of the studied catchment. Be careful with capital letters or such
#     input_dir : name of the folder (not the path !) containing the output text files computed by PINE for all 
#                 climate models. It should be under catchment_path and named 'output_corr' or 'output_delta'.
#     output_dir : name of the folder you want the output to be registered in.
#     output_subdir : name of the subdirectory of output_dir where will be stored the output plots and datafiles for low flows.
#     discharge_idd : string. identifier for the discharge file names (usually 'simdschrg').
#     scenarii : list of strings. Scenarii data computed by climate_models_RCP_sfVersion.R (usually rcp45 and rcp85)
#     date_format : string. Usually '%d.%m.%Y'.
#
#   files :
#     daily flow series for Baseline + 10 models
#
# OUTPUT : 
#   files :
#     catchment_7days_avg_corrordelta_scenario_period.txt : 7-day mean flow series for the whole period
#     catchment_7Qav_min_corrordelta_scenario_period.png : minimal 7-day means per year and model
#     catchment_7Qav_min_corrordelta_scenario_period.txt : associated data frame
#
# WARNING : 
#     if there are already files in the output directory that have the same names as the output files of this script, 
#     they might be overwritten.
#
# Load necessary library
library(data.table)
library(dplyr)
library(ggplot2)
library(lubridate)
library(zoo)
library(tidyr)
library(lmom)

# --------------------------------------------- SETUP ---------------------------------
catchment <- 'Myglevatn'
input_dir <- 'Input_plots'
output_dir <- 'Ouput_plots'
output_subdir <- 'annual min'

discharge_idd <- 'simdschrg'

scenarii <- list('rcp45', 'rcp85')
date_format <- "%Y-%m-%d"
# --------------------------------------------- SETUP ---------------------------------

# Get the directory where the script is, then go back up to klima
script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
# Delete all characters after the last '/'
klima <- sub("/[^/]*$", "", script_dir)
setwd(klima)
# if the automatic directory setting doesn't work, just precise the working directory below :
# setwd()
print(paste("wd: ", getwd()))

# Define the folder path
input_path <- paste(klima, input_dir, catchment, sep='/')
output_path <- paste(klima, output_dir, catchment, output_subdir, sep='/')
#print(paste('files coming from ', pine_output_path))
#print(paste('plots created in ', plot_output_path))

if (!file.exists(paste(klima, output_dir, catchment, sep='/'))) {
  dir.create(paste(klima, output_dir, catchment, output_subdir, sep='/'))
  print(paste("Output directory created successfully for",catchment))
}
if (!file.exists(output_path)) {
  # If not, create directory
  dir.create(output_path)
  print(paste("Output directory created successfully for",catchment))
}

# Read the input data files
input_all <- list.files(input_path)
# Filter the list to include only elements that contain the discharge  identifier in their name
filtered_input <- grep(discharge_idd, input_all, value = TRUE)

# --------------------------------- MINFLOWS COMPUTATION ------------------------------

#' find_data_info(data_frame_name, scenarios)
#' input :
#'  data_frame_name : string
#'  scenarios : list of studied scenarios to match data_frame_name's scenario
#' output :
#'  the scenario and the period of data_frame_name
find_data_info <- function(data_frame_name, scenarios) {
  # Use a regular expression to find the period pattern
  period_pattern <- "\\d{4}_\\d{4}"
  
  # Extract the period using the regular expression
  period_match <- regexpr(period_pattern, data_frame_name)
  period <- regmatches(data_frame_name, period_match)
  
  # Replace the underscore with a hyphen
  period_formatted <- gsub("_", "-", period)
  
  # Initialize variable for the scenario
  found_scenario <- NULL
  
  # Loop through each scenario to find a match in the input string
  for (scenario in scenarios) {
    if (grepl(scenario, data_frame_name)) {
      found_scenario <- scenario
      break
    }
  }
  
  # Assign the appropriate datatype based on the presence of 'simdschrg' or 'snowpack' in data_frame_name
  data_type <- ifelse(grepl("simdschrg", data_frame_name), "discharge", 
                      ifelse(grepl("snowpack", data_frame_name), "snowpack", "unknown"))
  
  corr_delta <- list('corrected', 'deltachanged')
  # Loop through each scenario to find a match in the input string
  for (modif in corr_delta) {
    if (grepl(modif, data_frame_name)) {
      found_modif <- modif
      break
    }
  }
  
  # Return scenario, formatted period, data type and 'corrected'/'deltachanged'
  return(list(scenario = found_scenario, period = period_formatted, datatype = data_type, modiftype = found_modif))
}

# Function to compute 7-day moving averages for each model within each year
compute_7day_moving_averages <- function(data) {
  data %>%
    arrange(Date) %>%
    group_by(Model, Year = year(Date)) %>%
    mutate(Value_7day = rollapply(Value, width = 7, by = 1, FUN = function(x) {
      if (length(x) < 7) return(NA) # If less than 7 days, return NA
      mean(x, na.rm = TRUE)
    }, align = "left", fill = NA, partial = TRUE)) %>%
    ungroup()
}

# Function to add sequence numbers within each year
add_sequence_numbers <- function(data) {
  data %>%
    group_by(Model, Year) %>%
    mutate(Seq = seq_along(Value_7day)) %>%
    ungroup()
}

# Function to compute minimum 7-day average values per year and model
compute_7day_lowest <- function(data) {
  data %>%
    group_by(Year, Model) %>%
    summarize(Min_7day = min(Value_7day, na.rm = TRUE), .groups = 'drop')
}

# Function to remove rows with NA or Inf values
remove_na_inf_rows <- function(dataframe) {
  # Check for NA or Inf in each row
  valid_rows <- apply(dataframe, 1, function(row) all(is.finite(row)) && all(!is.na(row)))
  # Subset the dataframe to include only valid rows
  cleaned_dataframe <- dataframe[valid_rows, ]
  return(cleaned_dataframe)
}

for (file_name in filtered_input){
  info <- find_data_info(file_name, scenarii)
  scenario <- info$scenario
  period <- info$period
  corr_or_delta <- info$modiftype
  print(paste('data processing :', corr_or_delta, scenario, period, file_name))
  
  # Extract some information from the name of the input file
  file_path <- paste(input_path, file_name, sep='/')
  # Read the data
  data <- read.table(file_path, header = TRUE)
  data$Date <- as.Date(data$Date, format = date_format)
  setDT(data)
  
  long_data <- melt(data, id.vars = "Date", variable.name = "Model", value.name = "Value")

  model_names <- colnames(data)[!(colnames(data) %in% c("Date", "Baseline"))]
  
  # Compute 7-day moving averages
  moving_averages <- compute_7day_moving_averages(long_data)
  
  # Add sequence numbers within each year
  moving_averages_seq <- add_sequence_numbers(moving_averages)
  
  moving_avg_series <- moving_averages_seq %>% dplyr::select(Year, Seq, Model, Value_7day) # reduce table
  
  avg_wide <- moving_avg_series %>%
    spread(key = Model, value = Value_7day)
  
  # remove rows with NA or Inf values
  avg_wide <- remove_na_inf_rows(avg_wide)
  
  # Save the moving averages with sequence numbers
  title1 <- paste0(catchment, "_7days_avg_", corr_or_delta, '_', scenario, '_',  period, ".txt")
  write.table(avg_wide, file = paste(output_path, title1, sep='/'), row.names = FALSE, sep = "\t", quote = FALSE)
  
  # Select the lowest 7-day mean flows
  # Compute the minimum 7-day average flow  per year and model
  
  nmq7 <- compute_7day_lowest(moving_averages_seq)
  seven_day_lowest <- nmq7 %>%
    rename(Min7day = Min_7day)
  
  q7 <- seven_day_lowest %>%  ## reshape
    spread(key = Model, value = Min7day)
  
  # remove rows with NA or Inf values
  q7 <- remove_na_inf_rows(q7)
  
  # Save the annual minimal 7-day average values
  title2 <- paste0(catchment, "_7Qav_min_", corr_or_delta, '_', scenario, '_', period, ".txt")
  write.table(q7, file = paste(output_path, title2, sep='/'), row.names = FALSE, sep = "\t", quote = FALSE)
  
  # Data Preparation
  q7_summary <- q7 %>%
    select(-Baseline) %>%
    pivot_longer(cols = -Year, names_to = "variable", values_to = "value") %>%
    group_by(Year) %>%
    summarise(min_value = min(value), max_value = max(value), mean_value = mean(value))
  
  baseline_data <- q7 %>%
    select(Year, Baseline)
  
  # Plotting
  p <- ggplot() +
    geom_ribbon(data = q7_summary, aes(x = Year, ymin = min_value, ymax = max_value, fill = "all models"), alpha = 0.5) +
    geom_line(data = q7_summary, aes(x = Year, y = mean_value, color = "Mean"), size = 0.5) +
    geom_line(data = baseline_data, aes(x = Year, y = Baseline, color = "Baseline"), size = 0.5) +
    theme_minimal() +
    labs(x = "Year", y = "lowest flow (m3/s)", title = "Yearly minimum 7-days-average flow") +
    scale_color_manual(name = NULL, values = c("Mean" = "blue", "Baseline" = "red")) +
    scale_fill_manual(name = NULL, values = c("all models" = "grey"))
  
  print(p)
  
  # Save plot under output_path directory, the name includes the catchment and the period
  plot_path <- paste0(output_path, '/', substr(title2, 1, nchar(title2)-4), '.png')
  ggsave(plot_path, plot=p, device="png", width=14,height=10,unit="in",dpi=200)
}


