########################################################
#######            PET calculator               ########
########################################################


##   Input: TXT file with monthly average temperature values per model in columns

##   Output: PET values in cm/month and in mm/day of all models

###   NOTE: VALID FOR LATITUDE 60. If different, K values should be changed 
####      accordingly --->  (see https://ponce.sdsu.edu/textbookhydrologytablea04.html)

library(dplyr)

######### Function to calculate PET using Thornthwaite method (LATITUDE 60)

PET_cm_month <- function(temperature_vector) {
  
  # Set negative temperature values to 0
  temperature_vector[temperature_vector < 0] <- 0
  
  # Calculate index I for each temperature value
  I <- (temperature_vector / 5)^1.514
  
  # Calculate index J as the sum of I values
  J <- sum(I)
  
  # Calculate c using the provided formula
  c <- 0.000000675 * J^3 - 0.0000771 * J^2 + 0.01792 * J + 0.49239
  
  # Calculate pet0 for each temperature value
  pet0 <- 1.6 * (10 * temperature_vector / J) ^ c
  
  # Define K values corresponding to each month ->>> VALID FOR LAT 60 !!!!
  K_values <- c(0.54, 0.67, 0.97, 1.19, 1.33, 1.56, 1.55, 1.33, 1.07, 0.84, 0.58, 0.48)
  
  # Get the corresponding K value for each temperature value based on the position in the temperature vector
  K <- K_values[seq_along(temperature_vector)]
  
  # Calculate PET per temperature value
  PET <- K * pet0
  
  # Return the set of resulting PET values
  return(PET)
}

###################################### Read input file ####################################################

tm <- read.table("TM_sju_LT.txt", header = TRUE)

# Extract temperature values (excluding month column)
tmp <- tm[, -1]

# Initialize list to store PET results for each model
PET_results <- list()

########### Loop through each model
for (model in colnames(tmp)) {
  
  # Calculate PET using Thornthwaite method for the current model
  PET_model <- PET_cm_month(tmp[[model]])
  
  # Store PET results for the current model
  PET_results[[model]] <- PET_model
}

# Combine PET results from all models into a single data frame
combined_PET <- do.call(cbind, PET_results)

# Write combined PET results to a CSV file
write.csv(combined_PET, file = "PET_results_cm_month_LT_Sju.csv", row.names = TRUE)

# Convert PET results to mm/day
PET_mm_day <- combined_PET * (10 / c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31))

# Write PET results converted to mm/day to a separate CSV file
write.csv(PET_mm_day, file = "PET_results_mm_day_LT_Sju.csv", row.names = TRUE)
