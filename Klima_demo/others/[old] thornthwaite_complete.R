#Script fo computation of monthly evapotranspiration

library(readxl)
library(ggplot2)
library(tidyverse)
library(xts)
library(hydroTSM)
library(dplyr)
library(SPEI)

script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
print(script_dir)

# --------------------------------- SETUP --------------------------------------
catchment <- 'Myglevatn'
catchment_lat <- 58.26452 #latitude of the measurement point with gmaps coordinates (spherical)
setwd("/Users/clementbernerd/Desktop/Final Data/PET/Delta 85")
pathh <- "/Users/clementbernerd/Desktop/Final Data/PET/Delta 85/Final PET Delta85.xlsx"
# --------------------------------- SETUP --------------------------------------

data1 <- excel_sheets(pathh)
data1 <- read_excel(pathh,sheet = 1,
                    col_types = c("text", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric"), 
                    skip = 0)
data1$Date <- as.Date(data1$Date,format="%Y-%m-%d")
?read_excel()
#Compute monthly temperature
snf <- as.data.frame(data1)
nc <- ncol(snf)
name <- names(snf)
data1 = as.data.frame(data1)



for (i in 2:nc){
  xtssn <- xts(snf[,i],order.by=snf[,1])
  msnow_check <- xtssn # <- daily2monthly(xtssn,FUN="mean")  #Runoff = mean
  allsnow_check <- monthlyfunction(msnow_check,FUN="mean")
  
  if (i==2){
    print(name[i])
    result <- data.frame(date=index(msnow_check), coredata(msnow_check))
    names(result) <- c(name[1],name[i])
  }
  else {
    print(name[i])
    temp <-  data.frame(coredata(msnow_check))
    names(temp) <- c(name[i])
    result <- cbind(result,temp)
  }
}
#view(result)
write.csv(result,"monthly_temp_precis.csv")




#Compute PET using Thornthwaite at monthly timescale

for (i in 2:nc){
  xtssn <- thornthwaite(result[,i],lat = catchment_lat)
  
  if (i==2){
    print(name[i])
    result2 <- data.frame(date=index(xtssn), coredata(xtssn))
    names(result2) <- c(name[1],name[i])
  }
  else {
    print(name[i])
    temp <-  data.frame(coredata(xtssn))
    names(temp) <- c(name[i])
    result2 <- cbind(result2,temp)
    
  }
}

write.csv(result2,"Valeurs_PET_precis.csv")

##Thornthwaite does not work at daily timescale

cc_1 <- thornthwaite(snf[,2], catchment_lat)
cc_2 <- thornthwaite(snf[,3], catchment_lat)
dai <- cbind.data.frame(cc_1,cc_2)
names(dai) <- c("obs","CNRM_CCLM")
write.csv(dai,"check.csv")

###############################




