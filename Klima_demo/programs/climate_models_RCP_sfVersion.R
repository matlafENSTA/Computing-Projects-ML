# CLIMATE MODEL DATA EXTRACTION
#
# PURPOSE : 
#     for a given catchment, extract the all the local future temperature and precipitation daily data sets (until 2100).
#     reads .nc climate data files and writes on .txt text files. In other words, converts Klima grided daily P & T to text files
#     VERSION USING SF.
#
# INPUT : 
#   variables :
#     catchment : name of the catchment. Be careful with upper case letters.
#     input_folder : the name of the BiasCorrect input folder. Usually 'Input_BiasCorrect'.
#     output_folder : the name of the BiasCorrect output folder. Usually 'Output_BiasCorrect'.
#
#     scenario : string. Name of the scenario you want to study. The data of this scenario is under Klima/RR/... and Klima/TM...
#       Usually 'rcp45' or 'rcp85'.
#     models : list of strings. identifies all the climate models you want to extract the data from. Usually 
#       list("CNRM_CCLM", "CNRM_RCA", "EC-EARTH_CCLM", "EC-EARTH_HIRHAM", "EC-EARTH_RACMO", "EC-EARTH_RCA", "HADGEM_RCA","IPSL_RCA", "MPI_CCLM", "MPI_RCA")
#
#   files : 
#     Klima/input_folder/catchment : catchment shapefiles (named catchment.shp)
#     Klima/RR : all precipitation netcdfs files(and nothing else)
#     Klima/TM : all temperature netcdfs files
#
# OUTPUT :
#   files : Klima/output_folder/catchment
#     hist_climatemodel_RR_daily_1971_v4-NA_RR.txt : precipitation historical values computed by the climate model.
#     hist_climatemodel_TM_daily_1971_v4-NA_T.txt : temperature historical values computed by the climate model.
#     scenario_climatemodel_RR_daily_1971_v4-NA_RR.txt : precipitation future values predicted by the climate model.
#     scenario_climatemodel_TM_daily_1971_v4-NA_T.txt : temperature future values predicted by the climate model.
#
#     catchment_Grid_nodesD.png :  map of the studied catchment.
#     catchment_Grid_nodesD.txt : data file with associated geographical coordinates. CAUTION : the coordinates system is the 
#       Norwegian EUREF89 - UTM33N. To get the spherical coordinates (the ones that are used on google maps for example), 
#       use https://norgeskart.no/ 
#
# CAUTION : 
#     RR and TM directories containing climate data have to be under Klima (see Klima_folder_structure.png)
#
# NB :
#     First run a test for a given scenario and climate model (list of one). Then do it for the list of all your climate models
#     and repeat it for every climate scenario you want to study.
#     Historical values are computed anyway, regardless of the scenario.
#
library(ncdf4)
library(zoo)
library(readr)
library(sf) # for geom_sf() in plotting
library(ggplot2)

# --------------------------------------------- SETUP ---------------------------------
catchment <- 'Myglevatn'
input_folder <- 'Input_RCP'
output_folder <- 'Output_RCP'

scenario <- "hist"
models <- list("test_model")
# --------------------------------------------- SETUP ---------------------------------

# get the directory where the script is. The input and ouput should be in the same directory
script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
# Delete all characters after the last '/'
klima <- sub("/[^/]*$", "/", script_dir)
setwd(klima)
# if the automatic directory setting doesn't work, just precise the working directory below :
# setwd()
print(paste("wd: ", getwd()))
input_folder <- paste(input_folder, catchment, sep='/')
output_folder <- paste(output_folder, catchment, sep='/')

if (!file.exists(output_folder)) {
  # If not, create directory
  dir.create(output_folder)
  print(paste("Output directory created successfully for",catchment))
} else {
  print(paste("Output directory already exists for", catchment))
}

# Record start time
start_time <- Sys.time()
# Read shapefile (the shape file should be in UTM zone 33 because the netCDF file is in UTM zone 33)

# --------------------------------------------- DATA PROCESSING -----------------------

#map <- rgdal::readOGR("Input/Sandungen.shp")
map <- read_sf(paste0(input_folder,'/',catchment,'.shp'))

#Read the first nc file that comes up to get the structure of the points.
test_repository <- paste("RR",models[[1]],"hist", sep="/")
file_list <- list.files(path = test_repository, recursive = TRUE)
mycdf <- nc_open(paste("RR",models[[1]],'hist',file_list[[1]], sep="/"), verbose = T, write = FALSE) #CHANGE NAME

#Grab the longitude (lon) and latitude (lat) and other data (tm, pr) from RCM
lat <- ncvar_get(mycdf,"Yc")
lon <- ncvar_get(mycdf,"Xc")

pr <- ncvar_get(mycdf, varid="precipitation__map_hist_daily", start=NA, count=NA, verbose=FALSE,
                signedbyte=TRUE, collapse_degen=TRUE)

##Construct a data frame from the given data structure, then they will be filled-in with lat, long and .nc other data

latlng1 <- data.frame()
prg_f <- data.frame()

# Define lat and long for each of the points existing in the .nc

latlong2 <- as.matrix(merge(lon,lat))
latlong3 <- rbind(latlng1, latlong2)
pktID <-  rownames(latlong3)

# 1550 should be changed based on the dim of the data #Climate Yc 
for ( i in 1:length(lat) ) {
  
  prg <- pr[, i,1]
  prg <- data.frame(prg)
  prg_f <- rbind(prg_f, prg)
  
}

prg_f <- cbind(prg_f,as.numeric(pktID))
names(prg_f) <-c("Var","pktID") 

#
# Create a shapefile of the gridpoints with data
# New version - rgdal/rgeos removed.
#
n <- cbind(latlong3,prg_f)
nn <- st_as_sf(n,coords = c("x", "y"))
st_crs(nn) <- "+proj=utm +zone=33 +datum=WGS84"

## Write above data as shape file (it goes into a folder node_Shapefile)

st_write(nn, dsn = paste(output_folder, "node_Shapefile", sep="/"),layer="node_Shapefile2", driver="ESRI Shapefile", append=FALSE)

## Read above data for later use (read the data from the created node_Shapefile)

xx <- read_sf(paste(output_folder, "node_Shapefile/node_Shapefile2.shp", sep = "/"))


# Get overlayed coordinates

op1 <- xx[map,]
op2 <- as.data.frame(op1)
idd <- op2$pktID #rownames(op2)
idd <- as.numeric(idd)

rr1 <- cbind(op2,st_coordinates(op1))
rr1 <- rr1[,-c(1,3)]

head(rr1)

colnames(rr1)[2]<-"longitude"
colnames(rr1)[3]<-"latitude"

head(rr1)

p_gridnode <- ggplot() + geom_sf(data = map) + geom_point(data=rr1,aes(x=longitude,y=latitude))
write.table(rr1, file= paste(output_folder, '/', catchment, "_Grid_nodesD.txt", sep = ""),col.names = T, row.names = F)
ggsave(filename = paste(output_folder, '/', catchment, "_Grid_nodesD.png", sep = ""), plot = p_gridnode, width = 8, height = 6, dpi = 300)

#-------------------------------------------------------------------------------------------------------------
# PROCESS ALL CLIMATE FILES
# -------------------------------------------------------------------------------------------------------------

# Precipitation

# change the working directory at each step (to treat all the models in a row)
for (model in models) {
  # Construct the directory path
  dir_path <- paste("RR", model, scenario, sep = "/") # CHANGE NAME
  setwd(dir_path)
  print(paste("input data folder : ", getwd()))
  
  lst <- list.files(getwd(), pattern="\\.nc$")
  
  # lst[5] here 5 is the number of files and it should be changed based on the number of files you have in the folder
  #ORIG nm<- paste(substr(lst[1], 1,25),substr(lst[94], 22,25), sep = '-')
  
  nm<- paste(substr(lst[1], 1,31),substr(lst[132], 29,32), sep = '-') ##CHANGE [5] to last file number
  rr_f1 <- data.frame()
  
  for (k in 1:length(lst)){
    rm(nc)
    print(lst[k])
    nc <- nc_open(lst[k], verbose = F, write = FALSE)
    
    lon <- ncvar_get(nc,"Xc")
    lat <- ncvar_get(nc,"Yc")
    
    latlng1 <- data.frame()
    latlong2 <- as.matrix(merge(lon,lat))   
    latlong3 <- rbind(latlng1, latlong2)  
    
    # Read correct varid
    #
    fn <- lst[k]
    ststring <- substr(fn,1,3)
    
    pr <- ncvar_get(nc, varid=paste("precipitation__map_",scenario,"_daily", sep=""), start=NA, count=NA, verbose=FALSE,
                      signedbyte=TRUE, collapse_degen=TRUE)
    
    
    #  fn <- lst[k]
    #  dt <- substr(fn, (nchar(fn)-19),(nchar(fn)-3))
    #  st <- as.Date(paste(substr(dt,1,4), substr(dt,5,6), substr(dt,7,8) , sep="-"))
    #  ed <- as.Date(paste(substr(dt,10,13), substr(dt,14,15), substr(dt,16,17) , sep="-"))
    #  dtt <- seq(st, ed, by="1 day")	
    aar <- parse_number(substring(fn,10,50))
    st <- as.Date(paste(aar,"01","01",sep="-"))
    ed <- as.Date(paste(aar,"12","31",sep="-"))
    dtt <- seq(st, ed, by="1 day")	
    rr_f <- data.frame()
    
    for (j in 1:dim(pr)[3]) {
      rr <- pr [, , j] [idd]
      rr <- as.data.frame(rr)
      colnames(rr)<- dtt[j]
      rr <- t(rr)
      
      rr_f <-  rbind(rr_f,rr)
    }
    
    rr_f1 <- rbind(rr_f1, rr_f)
    
  }
  
  rr_f1 <- round((rr_f1/10),2)  ##DIVISION BY 10 INSERTED
  
  rr_f11 <- as.data.frame(round(rowMeans(rr_f1, na.rm = TRUE),2))
  colnames(rr_f11) <- "Average"
  
  lat <- round(latlong3$y[idd],0)
  lon <- round(latlong3$x[idd],0)
  rr_f1 <- rbind(lon,lat,rr_f1)
  
  latt <- round(mean(latlong3$y[idd]),0)
  lonn <- round(mean(latlong3$x[idd]),0)
  rr_f11 <- rbind(lonn,latt,rr_f11)
  
  rr_fff <- cbind(rr_f1,rr_f11)
  rownames(rr_fff)[1] <- "Longitude"
  rownames(rr_fff)[2] <- "Latitude"
  
  fn_txt <- paste(nm,"_RR.txt", sep = "")
  
  setwd(klima)
  
  write.table(rr_fff, paste(output_folder, '/',fn_txt,sep=""), row.names = T, col.names = NA)
  
  # Reset working directory to the initial directory to avoid errors in the loop
  setwd(klima) # Assuming you are returning one level up after each iteration for demonstration
}

cat(paste0("\033[32m", "end of precipitation computation", "\033[0m", "\n"))

# -------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------

# Temperature
# see Precipitation for models (l.104)

#setwd(klima)
# change the working directory at each step
for (model in models) {
  dir_path <- paste("TM", model, scenario, sep = "/")
  setwd(dir_path) #direct path, may not be working
  print(paste("input data folder : ", getwd()))
  
  lst <- list.files(getwd(), pattern="\\.nc$")
  
  #ORIG nm<- paste(substr(lst[1], 1,25),substr(lst[94], 22,25), sep = '-')
  nm<- paste(substr(lst[1], 1,31),substr(lst[130], 29,32), sep = '-') ##CHANGE [5] to last file number
  
  rr_f1 <- data.frame()
  for (k in 1:length(lst)){
    rm(nc)
    print(lst[k])
    nc <- nc_open(lst[k], verbose = F, write = FALSE)
    
    lon <- ncvar_get(nc,"Xc")
    lat <- ncvar_get(nc,"Yc")
    
    latlng1 <- data.frame()
    latlong2 <- as.matrix(merge(lon,lat))   
    latlong3 <- rbind(latlng1, latlong2)  
    
    ##
    ## REMEMBER TO CHECK varid for the different files and edit _hist, _rcp45 and _rcp85 accordingly
    ##
    fn <- lst[k]
    ststring <- substr(fn,1,3)
    pr <- ncvar_get(nc, varid=paste("air_temperature__map_",scenario,"_daily", sep=""), start=NA, count=NA, verbose=FALSE,
                      signedbyte=TRUE, collapse_degen=TRUE)
    
    #  fn <- lst[k]
    #  dt <- substr(fn, (nchar(fn)-19),(nchar(fn)-3))
    #  st <- as.Date(paste(substr(dt,1,4), substr(dt,5,6), substr(dt,7,8) , sep="-"))
    #  ed <- as.Date(paste(substr(dt,10,13), substr(dt,14,15), substr(dt,16,17) , sep="-"))
    #  dtt <- seq(st, ed, by="1 day")	
    aar <- parse_number(substring(fn,10,50))
    st <- as.Date(paste(aar,"01","01",sep="-"))
    ed <- as.Date(paste(aar,"12","31",sep="-"))
    dtt <- seq(st, ed, by="1 day")	
    
    rr_f <- data.frame()
    
    for (j in 1:dim(pr)[3]) {
      rr <- pr [, , j] [idd]
      rr <- round(rr/10-273.15, 2) #rounds to two decimals
      rr <- as.data.frame(rr)
      colnames(rr)<- dtt[j]
      rr <- t(rr)
      
      rr_f <-  rbind(rr_f,rr)
    }
    
    rr_f1 <- rbind(rr_f1, rr_f)
    
  }
  
  rr_f1 <- round(rr_f1,2)
  
  rr_f11 <- as.data.frame(round(rowMeans(rr_f1, na.rm = TRUE),2))
  colnames(rr_f11) <- "Average"
  
  lat <- round(latlong3$y[idd],0)
  lon <- round(latlong3$x[idd],0)
  rr_f1 <- rbind(lon,lat,rr_f1)
  
  latt <- round(mean(latlong3$y[idd]),0)
  lonn <- round(mean(latlong3$x[idd]),0)
  rr_f11 <- rbind(lonn,latt,rr_f11)
  
  rr_fff <- cbind(rr_f1,rr_f11)
  rownames(rr_fff)[1] <- "Longitude"
  rownames(rr_fff)[2] <- "Latitude"
  
  fn_txt <- paste(nm,"_T.txt", sep = "")
  
  setwd(klima)
  current_output_folder <- paste0(output_folder, '/', catchment, "_", scenario, '/')  ##NEW DIRECTORY FOR OUTPUT
  if (!dir.exists(current_output_folder)) {
    dir.create(current_output_folder, recursive = TRUE)
  }
  # Define the full file path by adding the separator "/"
  file_path <- paste(current_output_folder, fn_txt, sep="/")
  
  # Write the table to the defined file path
  write.table(rr_fff, file_path, row.names = TRUE, col.names = NA)
  
  # Reset working directory to the initial directory to avoid errors in the loop
  setwd(klima) # Assuming you are returning one level up after each iteration for demonstration
}
cat(paste0("\033[32m", "end of temperature computation", "\033[0m", "\n"))

# ---------------------------------------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------------------------------

# time measurement
end_time <- Sys.time()
execution_time <- end_time - start_time
print(execution_time)
      
