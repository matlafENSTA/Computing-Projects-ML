# FINAL PLOTS OF AVERAGE YEARS FOR FUTURE PERIODS AND COMPARISON WITH HISTORICAL VALUES
#
# PURPOSE :
#     Extract the average snowpack and simulated discharge from PINE output files (DATA EXTRACTION)
#     visualise the simulated discharge for all given climate models, the baseline and the average of all climate models
#     on a plot. (PLOTS)
#     Visualise the average of the all the snowpacks (there is a snowpack value for each of the 10 elevation levels). 
#     Plots include the average of the climate models and the baseline (historical values).
#     Includes the computation of the drained water volume per year.
#
# INPUT  :
#   variables :
#     catchment : name of the studied catchment. Be careful with capital letters or such
#     input_dir : name of the folder (not the path !) containing the output text files computed by PINE for all 
#                 climate models. It should be under catchment_path and named 'output_corr' or 'output_delta'.
#     output_dir : name of the folder you want the output to be registered in.
#   subsidiary :
#     snowpack_idd : string. identifier for the snowpack file names (usually 'snowpack').
#     discharge_idd : string. identifier for the discharge file names (usually 'simdschrg').
#     snowpack_unit : string. unit of simulated snowpack values.
#     discharge_unit : string. unit of simulated discharge values.
#     Usually the 5 first letters. Example : the identifier of 'catchmentexample' is 'catch'.
#     output_plot : directory where will be stored the averageyear plots for simdschrg and snowpacks.
#     output_avgdata : the exact associated data frames (.csv) to the plots under output_plot (file names are the same).
#     scenarii : list of strings. Scenarii data computed by climate_models_RCP_sfVersion.R (usually rcp45 and rcp85)
#
#   files :
#     daily flow series for Baseline + 10 models
#
# OUTPUT :
#   files : they will be created under pine_output_path/catchment/output_plot and pine_output_path/catchment/output_avgdata.
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
input_dir <- 'Input_plots'
output_dir <- 'Ouput_plots'

# subsidiary
snowpack_idd <- 'snowpack'
discharge_idd <- 'simdschrg'
snowpack_unit <- 'mm'
discharge_unit <- 'm3/s'
output_plot <- 'average years plots'
output_avgdata <- 'average years data'
scenarii <- list('rcp45', 'rcp85')
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
input_path <- paste(input_dir, catchment, sep='/')
output_path <- paste(klima, output_dir, catchment, sep='/')
#print(paste('files coming from ', pine_output_path))
#print(paste('plots created in ', output_path))

if (!file.exists(output_path)) {
  dir.create(output_path)
  print(paste("Output directory created successfully for",catchment))
}
if (!file.exists(paste(output_path, output_plot, sep='/'))) {
  dir.create(paste(output_path, output_plot, sep='/'))
  print(paste("Subfolder created successfully : ",paste(output_path, output_plot, sep='/')))
}
if (!file.exists(paste(output_path, output_avgdata, sep='/'))) {
  dir.create(paste(output_path, output_avgdata, sep='/'))
  print(paste("Subfolder created successfully : ",paste(output_path, output_avgdata, sep='/')))
}

# Read the input data files
input_all <- list.files(input_path)
input_txt <- grep('.txt', input_all, value = TRUE)
# Filter the list to include only elements that contain the discharge  identifier in their name
discharge_all <- grep(discharge_idd, input_txt, value = TRUE)
snowpack_all <- grep(snowpack_idd, input_txt, value = TRUE)

# --------------------------------------------- PLOTS ---------------------------------

# Read the data from the text files
# data <- read.table("D:/Klima/Input_plotKlima/Runoff_long_term_Rinda (plotKlima example file).txt", header = TRUE)

#' find_data_info(filename, scenarios)
#' input :
#'  filename : string.
#'  scenarios : list of studied scenarios to match filename's scenario
#' output :
#'  the scenario and the period of filename
find_data_info <- function(filename, scenarios) {
  # Use a regular expression to find the period pattern
  period_pattern <- "\\d{4}_\\d{4}"
  
  # Extract the period using the regular expression
  period_match <- regexpr(period_pattern, filename)
  period <- regmatches(filename, period_match)
  
  # Replace the underscore with a hyphen
  period_formatted <- gsub("_", "-", period)
  
  # Initialize variable for the scenario
  found_scenario <- NULL
  
  # Loop through each scenario to find a match in the input string
  for (scenario in scenarios) {
    if (grepl(scenario, filename)) {
      found_scenario <- scenario
      break
    }
  }
  
  # Assign the appropriate datatype based on the presence of 'simdschrg' or 'snowpack' in filename
  data_type <- ifelse(grepl("simdschrg", tolower(filename)), "discharge", 
                      ifelse(grepl("snowpack", tolower(filename)), "snowpack", "unknown"))
  
  corr_delta <- list('corrected', 'deltachanged')
  # Loop through each scenario to find a match in the input string
  for (modif in corr_delta) {
    if (grepl(modif, filename)) {
      found_modif <- modif
      break
    }
  }
  
  # Return scenario, formatted period, data type and 'corrected'/'deltachanged'
  return(list(scenario = found_scenario, period = period_formatted, datatype = data_type, modiftype = found_modif))
}

#' average_func(data_frame)
#' INPUT : string. Name of a frame (large matrix of characters to be precise) that has just been created
#'         in the DATA EXTRACTION part above.
#' OUPUT : data ready to be plotted 
average_func <- function(data_frame) {
  # Read the data frame
  dataframe <- as.data.frame(data_frame)
  # Convert Date to proper Date format
  dataframe$Date <- as.Date(dataframe$Date, "%Y-%m-%d")
  
  # Convert character columns to numeric
  dataframe <- dataframe %>% mutate_if(is.character, as.numeric)
  
  # Extract day and month from Date
  dataframe$DayMonth <- format(dataframe$Date, "%m-%d")
  
  # Add a column for month names
  dataframe$Month <- month(dataframe$Date)
  
  # Calculate the mean discharge for each day and month combination, excluding non-numeric columns
  avg <- dataframe %>% group_by(DayMonth) %>% summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))
  
  # Initialize baseline  with the number of days column
  base <- data.frame(DayMonth = avg$DayMonth,
                     Baseline = avg$Baseline)
  base <- base %>% mutate(DayNumber = row_number())
  
  # Remove the last row
  base <- base[-nrow(base), ]
  # Remove the 29th of February which is badly handled
  base <- subset(base, DayMonth != '02-29')
  
  # Add Month column with month names
  avg <- avg %>% 
    mutate(Month = month(as.Date(paste("2000", DayMonth, sep = "-"), format = "%Y-%m-%d")))
  
  # Add a column for month names
  ##dataframe$Month <- month(dataframe$Date, label = TRUE)
  
  # Exclude Baseline values from the main dataframe (GCM models)
  avg <- avg[, !colnames(avg) %in% "Baseline"]
  
  # Assign a sequential number to each date
  avg <- avg %>% mutate(DayNumber = row_number())
  
  # Remove the last row
  avg <- avg[-nrow(avg), ]
  # Remove the 29th of February which is badly handled
  avg <- subset(avg, DayMonth != '02-29')
  
  # Create a new dataframe for max curve
  max_curve_AVG <- data.frame(Month = avg$Month, DayNumber = avg$DayNumber, MaxDischarge = apply(avg[, !(names(avg) %in% c("DayMonth", "Month", "DayNumber"))], 1, max, na.rm = TRUE))
  # Create a new dataframe for min curve
  min_curve_AVG <- data.frame(Month = avg$Month,
                              DayNumber = avg$DayNumber, 
                              MinDischarge = apply(avg[, !(names(avg) %in% c("DayMonth", "Month", "DayNumber"))], 1, min, na.rm = TRUE))
  # Create a new dataframe for mean curve
  mean_curve_AVG <- data.frame(Month = avg$Month,
                           DayNumber = avg$DayNumber, 
                           MeanDischarge = rowMeans(avg[, !(names(avg) %in% c("DayMonth", "DayNumber", "Month"))], na.rm = TRUE))
  
  return(list(baseline = base, average = avg, mean_curve = mean_curve_AVG, min_curve = min_curve_AVG, max_curve = max_curve_AVG))
}

#' final_plot_and_datafile(base, avg, mean_curve_AVG, min_curve_AVG, max_curve_AVG, data_type, scenario_, period_, corr_or_delta)
#' OUTPUT : plots under the format .png and the associated data frames in .csv stored under output_plot and output_avgdata
final_plot_and_datafile <- function(base, avg, mean_curve, min_curve, max_curve, data_type, scenario_, period_, corr_or_delta){
  # find the appropriate unit
  data_unit <- ifelse(data_type == "discharge", discharge_unit, ifelse(data_type == "snowpack", snowpack_unit, "unknown unit"))
  
  # Integrate through the all year to get the volume of drained water. Do it for mean curve and baseline.
  historical_volume <- integrate.xy(base$DayNumber, base$Baseline) * 24 * 3600 # discharge is in m3/s
  future_volume <- integrate.xy(mean_curve$DayNumber, mean_curve$MeanDischarge) * 24 * 3600
  
  # Create the annotation text (adapt for discharge and snowpacks cases)
  annot_unit <- ifelse(grepl("discharge", data_type), "drained", 
                       ifelse(grepl("snowpack", data_type), "snowpack", "unknown"))
  annotation_text <- paste("baseline",annot_unit, "volume :", scientific(historical_volume, digits = 3), "m3/year\nfuture",annot_unit, "volume:", scientific(future_volume, digits = 3), "m3/year")
  
  # Plot the curves with months on the x-axis (just a year)
  p1 <- ggplot() +
    geom_ribbon(data = max_curve, aes(x = DayNumber, ymin = MaxDischarge, ymax = min_curve$MinDischarge, fill = "all models")) +
    geom_line(data = mean_curve, aes(x = DayNumber, y = MeanDischarge, color = "Mean"), size = 0.5) +
    geom_line(data = base, aes(x = DayNumber, y = Baseline, color = "Baseline"), size = 0.5) +
    scale_x_continuous(breaks = seq(1, 366, length.out = 6), 
                       labels = unique(avg$Month)[seq(1, length(unique(avg$Month)), length.out = 6)]) +
    labs(x = " ", y = paste0(data_type, ' (', data_unit,')'), 
         title = paste0(corr_or_delta, ' ', data_type, " comparison\n", scenario_, ' ', period_, " projections and Baseline - ", catchment)) +
    theme_minimal() +
    scale_color_manual(name = NULL, values = c("Mean" = "blue", "Baseline" = "red")) +
    scale_fill_manual(name = NULL, values = c("all models" = "lightgray"), labels = c("all models" = "all models")) +
    annotate("text", x = Inf, y = Inf, label = annotation_text, hjust = 1.1, vjust = 1.1, size = 3.5, color = "black")
  
  print(p1)
  
  
  # Save plot under output_path directory, the name includes the catchment and the period
  plot_name <- paste0(catchment, ' ', data_type, ' ', corr_or_delta, ' averageyear ', scenario_, ' ', period_, '.png')
  ggsave(paste(output_path, output_plot, plot_name, sep='/'), plot=p1, device="png", width=14,height=10,unit="in",dpi=200)
  
  ggplot() +
    geom_line(data = , aes(x = Date, 
                           y = CNRM_CCLM), size = 0.5)
  
  # csv file associated to the plot ; includes the baseline data at the end
  # Select columns from avg excluding 'Month' and 'DayNumber'
  avg_filtered <- avg[, !(colnames(avg) %in% c('Month', 'DayNumber'))]
  # Select the 'Baseline' column from base
  base_filtered <- base[, 'Baseline', drop = FALSE]
  # Concatenate the filtered data frames using cbind
  average_year <- cbind(avg_filtered, base_filtered)
  write.csv(average_year,paste0(output_path,'/', output_avgdata, '/', substr(plot_name, 1, nchar(plot_name) - 4), '.csv'))
  
  # # Plot average of all curves projections
  # """
  # ggplot() +
  #   geom_ribbon(data = max_curve, aes(x = DayNumber, 
  #                                     ymin = MaxDischarge, 
  #                                     ymax = min_curve$MinDischarge), fill = "lightgray") +
  #   geom_line(data = mean_curve, aes(x = DayNumber, y = MeanDischarge, color = "Mean"), size = 0.5) +
  #   geom_line(data = base, aes(x = DayNumber, y = Baseline, color = "Baseline"), size = 0.5) +
  #   scale_x_continuous(breaks = seq(1, max($DayNumber), length.out = 6), 
  #                      labels = unique($Month)[seq(1, length(unique($Month)), length.out = 6)]) +
  #   labs(x = " ", y = "Discharge m3/s", title = "Daily Discharge from 01-09-2041 to 31-08-2042 - Surna", color = "Legend") + theme_minimal() +
  #   scale_color_manual(values = c("Mean" = "blue", "Baseline" = "red"), labels = c("Baseline", "Mean"))
  
}

# make the averageyear plots and the associated .csv for all the data sets (snowpacks and simulated discharge)
# DISCHARGE
for (current_discharge_filename in discharge_all){
  current_discharge_dataframe <- read.table(paste(input_path, current_discharge_filename, sep='/'), header = TRUE)
  data_info <- find_data_info(current_discharge_filename, scenarii)
  datatype_discharge <- data_info$datatype
  period_discharge <- data_info$period
  scenario_discharge <- data_info$scenario
  corr_or_delta <- data_info$modiftype
  print(paste('data processing :', datatype_discharge, corr_or_delta, scenario_discharge, period_discharge))
  result_discharge <- average_func(current_discharge_dataframe)
  base_discharge <- result_discharge$baseline
  avg_discharge <- result_discharge$average
  mean_curve_discharge <- result_discharge$mean_curve
  min_curve_discharge <- result_discharge$min_curve
  max_curve_discharge <- result_discharge$max_curve
  final_plot_and_datafile(base_discharge, avg_discharge, mean_curve_discharge, min_curve_discharge, max_curve_discharge, datatype_discharge, scenario_discharge, period_discharge, corr_or_delta)
}
# SNOWPACKS
for (current_snowpack_filename in snowpack_all){
  current_snowpack_dataframe <- read.table(paste(input_path, current_snowpack_filename, sep='/'), header = TRUE)
  data_info <- find_data_info(current_snowpack_filename, scenarii)
  datatype_snowpack <- data_info$datatype
  period_snowpack <- data_info$period
  scenario_snowpack <- data_info$scenario
  corr_or_delta <- data_info$modiftype
  print(paste('data processing :', datatype_snowpack, corr_or_delta, scenario_snowpack, period_snowpack))
  result_snowpack <- average_func(current_snowpack_dataframe)
  base_snowpack <- result_snowpack$baseline
  avg_snowpack <- result_snowpack$average
  mean_curve_snowpack <- result_snowpack$mean_curve
  min_curve_snowpack <- result_snowpack$min_curve
  max_curve_snowpack <- result_snowpack$max_curve
  final_plot_and_datafile(base_snowpack, avg_snowpack, mean_curve_snowpack, min_curve_snowpack, max_curve_snowpack, datatype_snowpack, scenario_snowpack, period_snowpack, corr_or_delta)
}
