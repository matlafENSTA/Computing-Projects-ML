# BIAS CORRECTION AND DELTACHANGE FOR TEMPERATURE
#
# PURPOSE :
#     Correct the bias and delta change the climate data extracted by climate_models_RCP_sfVersion.R knowing the observed values
#     The program uses QMap for bias correction (Gudmunsson et al). Also computes the daily evapotranspiration (monthly average)
#
# INPUT :
#   variables :
#     catchment : name of the catchment. Be careful with upper case letters.
#     input_folder : the name of the BiasCorrect input folder. Usually 'Input_BiasCorrect'.
#     output_folder : the name of the BiasCorrect output folder. Usually 'Output_BiasCorrect'.
#
#     scenarii : list of strings. Scenarii data computed by climate_models_RCP_sfVersion.R (usually rcp45 and rcp85)
#     startYs : list of integers. The starting years of the future periods you want to compute the data.
#         Usually list(2041,2071) for near future and far future.
#     endYs : list of integers. The ending years of the future periods you want to compute the data.
#         Usually list(2070,2100) for near future and far future. The studied periods will then be 2041-2070 and 2071-2100.
#     latitude : latitude of the measurement station of the catchment. This coordinates is used for evapotranspiration
#         computation (find it on google maps for example). To get this latitude, find on 
#         output_RCP/catchment/catchment_Grid_nodesD.png the nearest point for the measurement station, and convert its coordinates
#         by you can convert them onto https://norgeskart.no/ by pointing the measurement station on the map. The coordinates 
#         system of catchment_Grid_nodesD.png is the Norwegian EUREF89 - UTM33N.
#
#   files : have to be under input_folder/catchment.
#     observed_temp.csv : Observed data from gauge/cell. This file is created by csv_to_txt_Sildre_SeNorge.py
#     hist_temp.csv : Historical data from one or more climate models: Created by climate_models_RCP_to_BiasCorrect.py
#         One column for each climate model, heading model name.
#     scenarioname_temp.csv : Future prediction for one or more climate models. Created by climate_models_RCP_to_BiasCorrect.py
#         One column for each climate model, heading model name.
#     Format: Dates as YYYY-mm-dd and column name "Date". Then, one column for each model with e.g. model name as identifier
#     Missing data should be removed before reading. Stops if NA is found.
#
# OUTPUT :
#   files : will be created under output_folder/catchment.
#     temp_corrected_histo_histperiod.csv -> Bias corrected historical data (histperiod has the format YYYY-YYYY)
#     temp_corrected_scenario_period.csv -> Bias corrected future data (scenario has lenght 5 and period has the format YYYY-YYYY)
#     temp_deltachanged_scenario_period.csv -> observed temperature data adjusted with delta method (for the historical period)
#
#     MonthlyT_histo_histperiod.csv -> The raw historical data, monthly resolution
#     MonthlyT_obser_histperiod.csv -> The raw observed data, monthly resolution
#     MonthlyT_histc_scenario_period.csv -> The corrected historical data, monthly resolution.
#     MonthlyT_scenario_period.csv -> bias corrected scenario data, monthly resolution
#     DailyPET_histo_period.csv -> the average of evapotranspiration by DAY for each MONTH (in mm/day) 
#     (corrected historical values)
#     DailyPET_scenario_period.csv -> the average of evapotranspiration by DAY for each MONTH (in mm/day)
#     why taking the daily average PET for each month ? because PINE takes daily PET in parameter
#
#     Graphs showing monthly data
#     Delta change factors for all periods and scenarii given.
#
# CAUTION :
#     if you change the names of the output files, the python programs to make the setups for PINE may not work !
#     You can adapt them (corrected_txt_setups.py and pine_files_creation.py)
#
# NB :
#     The columns should be in the same order in both files. Do also make sure that the historical and observed 
#     have a period in common big enough (longer than endY - startY).
#     if something goes wrong, the most frequent reason is a wrong folder arrangement, then a wrong format for the input files
#     if you wish to make some changes on this function, you better remove first the biggest 'for' loops (scenarii and period)
#
library(qmap)
library(ggplot2)
library(hydroGOF)
library(xts)
library(hydroTSM)
library(dplyr)
library(tidyr)
library(gt)
library(SPEI) # for thornthwaite function
library(lubridate) # to calculate the daily mean of evapotranspiration for each month

# --------------------------------------------- SETUP ---------------------------------
catchment <- 'Myglevatn'
input_folder <- 'Input_BiasCorrect'
output_folder <- 'Output_BiasCorrect'

scenarii <- list('rcp45', 'rcp85')
startYs <- list(2041,2071)
endYs <- list(2070,2100)
latitude <- 58.26452
# --------------------------------------------- SETUP ---------------------------------

# get the directory where the script is. The input and ouput should be in the same directory
script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
# Delete all characters after the last '/'
klima <- sub("/[^/]*$", "/", script_dir)
setwd(klima)
# if the automatic directory setting doesn't work, just precise the working directory below :
# setwd()
print(paste("wd: ", getwd()))

# automatic subfolder path creation
input_folder <- paste(input_folder, catchment, sep='/')
output_folder <- paste(output_folder, catchment, sep ='/')

if (!file.exists(output_folder)) {
  # If not, create directory
  dir.create(output_folder)
  print(paste("Output directory created successfully for",catchment))
} else {
  print(paste("Output directory already exists for", catchment))
}

# --------------------------------------------- DATA PROCESSING -----------------------

for (scenario in scenarii){
  for (i in 1:length(startYs)){
    startY <- startYs[[i]]
    endY <- endYs[[i]]
    if (endY < startY || startY < 2006 || endY > 2100){
      stop(paste0("endY and startY should verify 2006 < startY < endY < 2100, but you entered ", startY, "-", endY))
    }
    
    observed_temp <- paste(input_folder,'observed_temp.csv', sep='/')
    historical_temp <- paste(input_folder,'hist_temp.csv', sep='/')
    scenario_temp <- paste(input_folder, '/', scenario, '_temp.csv', sep='')
    
    if (!file.exists(output_folder)) {
      # If not, create directory
      dir.create(output_folder)
      print(paste("Output directory created successfully for",catchment))
    } else {
      print(paste("Output directory already exists for", catchment))
    }
    
    #
    # Read observed data, make the date into an R date class. 
    #
    obs <- read.csv(observed_temp, header=TRUE, stringsAsFactors=FALSE)
    obs$Date <- as.Date(obs$Date, format = '%Y-%m-%d')
    
    #
    # Read historical RCM data.
    #
    hist <- read.csv(historical_temp, header=TRUE,sep=",") # make sure the separator is the good one !!
    hist$Date <- as.Date(hist$Date)
    
    #
    # Read scenario RCM data. Set years for subset.
    #
    #make sure startY and endY are well defined according and consistent with the file
    scen <- read.csv(scenario_temp, header=TRUE, sep=",") # make sure the separator is the good one !!
    scen$Date <- as.Date(scen$Date)
    sDate=as.Date(paste(startY-1,"12","31",sep="-"))
    eDate=as.Date(paste(endY,"12","31",sep="-"))
    
    scen <- subset(scen,Date>sDate & Date<=eDate)
    
    period <- paste('_',startY,'-',endY,sep='')
    duration <- as.numeric(eDate - sDate) # duration in days
    
    # Some checks
    # 
    # the historical data and the observed data have not the same lenght, let's fix this :
    
    # matching the periods of historical and observed values
    first_date_obs <- obs$Date[1]
    end_date_obs <- obs$Date[length(obs$Date)]
    first_date_hist <- hist$Date[1]
    end_date_hist <- hist$Date[length(hist$Date)]
    
    first_date <- max(first_date_hist, first_date_obs)
    # the duration of historical data set should match the duration of the climate models data sets. 
    #If not, you need to reduce the duration until it gets right for the historical and observed data set.
    
    end_date <- min(end_date_hist, end_date_obs, first_date + duration)
    # period for historical and observed values (only years)
    period_hist <- paste0("_", substr(as.character(first_date),1,4), "-", substr(as.character(end_date),1,4))
    
    # the historical data and the observed data have not the same lenght, let's fix this :
    if (length(hist$Date) != length(obs$Date)) {
      
      if (end_date - first_date < duration){
        print(paste0("the observed and historical values for this 
                 catchment haven't enough dates in common to simulate on a period as wide as the one you gave : ", duration, "years"))
        stop(paste0("You should reduce the period duration to : ", end_date - first_date, "days"))
      }
      # Find the indices corresponding to the first_date and end_date
      obs_indices <- which(obs$Date >= first_date & obs$Date <= end_date)
      hist_indices <- which(hist$Date >= first_date & hist$Date <= end_date)
      
      if (length(obs_indices) == length(hist_indices)) {
        obs <- subset(obs, Date >= first_date & Date <= end_date)
        hist <- subset(hist, Date >= first_date & Date <= end_date)
      } else {
        stop("The lengths of the filtered observed and historical data do not match. 
         Make sure that the data files have one and only one value per day for example.")
      }
    }
    
    # Find the number of models, loop through them and run Qmap for each.
    
    out <- hist[,1]
    fut <- scen[,1]
    
    for (i in 2:ncol(hist)){
      #Correct and store the adjusted historical series.
      #Then proceed to correct the future series using the same quantile.
      #
      temp <- hist[,c(1,i)]
      fq <- fitQmap(obs[,2],temp[,2],method="QUANT", wet.day=FALSE, qstep=0.01, nboot=3, type="tricube")
      corrP <- doQmap(temp[,2],fq)
      corrP <- as.data.frame(corrP)
      colName <- names(temp)[2]
      names(corrP) <- colName
      print(paste("Historical bias correction for: ", scenario, period, colName))
      
      out <- cbind(out,corrP)
      
      # Correct the future
      #
      futemp <- scen[,c(1,i)]
      fucolName <- names(futemp)[2]
      print(paste("Future bias correction for : ", scenario, period, fucolName))
      if (fucolName != colName) stop("Error: Column names in historical and future does not match")
      
      futureP <- doQmap(futemp[,2],fq)
      futureP <- as.data.frame(futureP)
      names(futureP) <- fucolName
      fut <- cbind(fut,futureP)
    }
    
    # Add observed data to out just for control
    #
    obsP <- as.data.frame(obs[,2])
    names(obsP) <- "Observed"
    out <- cbind(out,obsP)
    #
    # Write the corrected data as csv files
    #
    write.csv(out,paste0(output_folder,"/temp_corrected_histo",period_hist,".csv"))
    write.csv(fut,paste0(output_folder,"/temp_corrected_", scenario,period,".csv"))
    
    #
    # COMPUTE THE MONTHLY MEANS
    # Control and basis for delta changes.
    #
    ############################################################################
    # Observed monthly means
    #
    xtsObs <- xts(obs[,2],order.by=as.Date(obs[,1]))
    mObs <- daily2monthly(xtsObs,FUN="mean")  #Temperarture = average
    allObs <- monthlyfunction(mObs,FUN="mean")
    obsFrame <- data.frame(date=index(allObs), coredata(allObs))
    names(obsFrame) <- c("Date", "Observed")
    
    # observed values for temperature
    write.csv(obsFrame,paste0(output_folder, "/MonthlyT_obser",period_hist,".csv"))
    
    ###########################################################################
    # Monthly means for the raw climate data and the corrected climate data
    # Observed series is added to both for comparison
    #
    rawFrame <- obsFrame
    corrFrame <- obsFrame
    
    for (j in 2:ncol(hist)){
      
      xtsRaw <- xts(hist[,j],order.by=as.Date(hist[,1]))
      mRaw <- daily2monthly(xtsRaw,FUN="mean")  #Temperature = mean
      allRaw <- monthlyfunction(mRaw,FUN="mean")
      tf <- data.frame(date=index(allRaw), coredata(allRaw))
      rawName <- names(hist)[j]
      rawFrame <- mutate(rawFrame, !! rawName := tf[,2])
      
      xtsC <- xts(out[,j],order.by=as.Date(out[,1]))
      mC <- daily2monthly(xtsC,FUN="mean")  #Precip = mean
      allC <- monthlyfunction(mC,FUN="mean")
      tc <- data.frame(date=index(allC), coredata(allC))
      corrName <- names(out)[j]
      corrFrame <- mutate(corrFrame, !! corrName := tc[,2])
    }
    
    write.csv(rawFrame,paste0(output_folder,"/MonthlyT_histo",period_hist,".csv"))
    # bias corrected historical values for temperature
    write.csv(corrFrame,paste0(output_folder,"/MonthlyT_histc",period_hist,".csv"))
    
    # Evapotranspiration (PET) historical (not needed)
    
    ## Keep the names of the climate headers
    #headers_hist <- names(rawFrame)
    ## Initialize a data frame with the Date column for monthly evapotranspiration
    #pet_hist <- data.frame(Date = rawFrame[["Date"]])
    ## apply thornthwaite method to calculate the evapotranspiration for all the climate models
    #for (i in 2:ncol(rawFrame)){
    #pet_hist[[headers_hist[i]]] <- thornthwaite(rawFrame[, i], lat = latitude)
    #}
    
    #evapo_hist <- paste(output_folder,'/',"MonthlyPET_Historical",period_hist,".csv", sep='')
    #write.csv(pet_hist, file = evapo_hist, row.names = FALSE)
    
    # Evapotranspiration (PET) historical corrected
    
    # Mapping months to the number of days
    month_days <- c(Jan = 31, Feb = 28, Mar = 31, Apr = 30, May = 31, Jun = 30,
                    Jul = 31, Aug = 31, Sep = 30, Oct = 31, Nov = 30, Dec = 31)
    
    # Keep the names of the climate headers
    headers_hist_corr <- names(corrFrame)
    
    # Initialize a data frame with the Date column for monthly evapotranspiration
    pet_hist_corr <- data.frame(Date = corrFrame[["Date"]])
    
    # Apply Thornthwaite method to calculate the monthly evapotranspiration for all the climate models (column by column)
    for (i in 2:ncol(corrFrame)) {
      monthly_pet <- thornthwaite(corrFrame[, i], lat = latitude)
      
      # Extract the month from the Date column
      month <- corrFrame[["Date"]]
      
      # Calculate the number of days in each month
      days_in_month <- sapply(month, function(m) month_days[m])
      
      # Convert monthly PET to daily PET
      daily_pet <- monthly_pet / days_in_month
      
      # Add daily PET to the data frame
      pet_hist_corr[[headers_hist_corr[i]]] <- daily_pet
    }
    
    # Define the output file path and write the result to a CSV file
    evapo_hist_corr <- paste0(output_folder, "/DailyPET_histo", period_hist, ".csv")
    write.csv(pet_hist_corr, file = evapo_hist_corr, row.names = FALSE)
    
    ######################################
    # SOME CHECKS AND A PLOT
    #
    colMeans(rawFrame[2:ncol(rawFrame)])
    colMeans(corrFrame[2:ncol(corrFrame)])
    
    of <- rawFrame[,1:2]
    rf <- rawFrame[,-2]
    cf <- corrFrame[,-2]
    lrf <- pivot_longer(rf,cols=2:ncol(rf))
    lrf <- mutate(lrf,Type="Raw")
    lcf <- pivot_longer(cf,cols=2:ncol(cf))
    lcf <- mutate(lcf,Type="Corrected")
    if (ncol(rf) > 2){ #One column of corrected for each climate model.
      tof <- of[rep(2, each = (ncol(rf))-2)]
      of <- cbind(of,tof)
    }
    names(of) <- names(rf)
    lof <- pivot_longer(of,cols=2:ncol(of))
    lof <- mutate(lof,Type="Observed")
    
    plFr <- rbind(lrf,lcf)
    plFr <- rbind(plFr,lof)
    
    p1 <- ggplot(plFr) + geom_bar(aes(x=Date,y=value,fill=Type),stat="identity", position="dodge") +
      facet_wrap(.~name)+
      scale_x_discrete(breaks = function(x){x[c(TRUE, FALSE)]})
    print(p1)
    ggsave(paste0(output_folder,"/TemperatureCompare",period_hist,".png"),plot=p1,device="png",width=14,height=10,unit="in",dpi=200)
    
    ##########################################################################
    # Process the future
    #
    futureFrame <- as.data.frame(obsFrame[,1]) 
    names(futureFrame) <- "Date"
    
    for (j in 2:ncol(fut)){
      
      xtsF <- xts(fut[,j],order.by=as.Date(fut[,1]))
      mF <- daily2monthly(xtsF,FUN="mean")  #Temperature = mean
      allF <- monthlyfunction(mF,FUN="mean")
      tff <- data.frame(date=index(allF), coredata(allF))
      fuName <- names(hist)[j]
      futureFrame <- mutate(futureFrame, !! fuName := tff[,2])
    }
    
    write.csv(futureFrame,paste0(output_folder,"/MonthlyT_", scenario,period,".csv"))
    
    
    # Evapotranspiration (PET) for future
    
    # Keep the names of the climate headers
    headers_future <- names(futureFrame)
    # Initialize a data frame with the Date column for monthly evapotranspiration
    pet_future <- data.frame(Date = futureFrame[["Date"]])
    
    # Apply Thornthwaite method to calculate the monthly evapotranspiration for all the climate models (column by column)
    for (i in 2:ncol(futureFrame)) {
      monthly_pet <- thornthwaite(futureFrame[, i], lat = latitude)
      
      # Extract the month from the Date column
      month <- futureFrame[["Date"]]
      
      # Calculate the number of days in each month
      days_in_month <- sapply(month, function(m) month_days[m])
      
      # Convert monthly PET to daily PET
      daily_pet <- monthly_pet / days_in_month
      
      # Add daily PET to the data frame
      pet_future[[headers_future[i]]] <- daily_pet
    }
    
    # Define the output file path and write the result to a CSV file
    evapo_future <- paste0(output_folder, "/DailyPET_", scenario, period, ".csv")
    write.csv(pet_future, file = evapo_future, row.names = FALSE)
    
    # #monthly evapotranspiration values for the future (not needed)
    # # Keep the names of the climate headers
    # headers_future <- names(futureFrame)
    # # Initialize a data frame with the Date column for monthly evapotranspiration
    # pet_future <- data.frame(Date = futureFrame[["Date"]])
    # # apply thornthwaite method to calculate the evapotranspiration for all the climate models
    # for (i in 2:ncol(futureFrame)){
    #   pet_future[[headers_future[i]]] <- thornthwaite(futureFrame[, i], lat = latitude)
    # }
    # 
    # evapo_future <- paste(output_folder,'/',"MonthlyPET_Future",'_', scenario,period,".csv", sep='')
    # write.csv(pet_future, file = evapo_future, row.names = FALSE)
    
    ########################################################
    # sum of precip on the screen and update to the plot
    colSums(futureFrame[2:ncol(futureFrame)])
    
    ff <- futureFrame
    lff <- pivot_longer(ff,cols=2:ncol(ff))
    lff <- mutate(lff,Type="Scenario")
    plFr <- rbind(plFr,lff)
    
    p1 <- ggplot(plFr) + geom_bar(aes(x=Date,y=value,fill=Type),stat="identity", position="dodge") +
      facet_wrap(.~name)+
      scale_x_discrete(breaks = function(x){x[c(TRUE, FALSE)]})
    print(p1)
    ggsave(paste0(output_folder,'/TemperatureCompareScenario_', scenario,period,".png"),plot=p1,device="png",width=14,height=10,unit="in",dpi=200)
    
    
    #############################################################
    # DELTA CHANGES
    #
    df <- cbind(lcf,lff)
    df <- df[,c(-4,-5,-6,-8)]
    names(df) <- c("Month","Model","Corrected","Scenario")
    
    df <- mutate(df,DeltaCh = round((Scenario - Corrected),3))
    toD <- df[,c(1,2,5)]
    deltaOut <- pivot_wider(toD,names_from=Model,values_from=DeltaCh)
    
    write.csv(deltaOut,paste0(output_folder,'/DeltaChangeFactorsT_', scenario,period,".csv"))
    
    delplot <- pivot_longer(deltaOut,cols=2:11)
    p2 <- ggplot(delplot) + geom_bar(aes(x=Month,y=value),stat="identity", position="dodge") +
      #geom_hline(yintercept=1.0) +
      facet_wrap(.~name)+
      scale_x_discrete(breaks = function(x){x[c(TRUE, FALSE)]})
    print(p2)
    ggsave(paste(output_folder,'/',"DeltaChangeFactorsT",'_', scenario,period,".png", sep=''),plot=p2,device="png",width=14,height=10,unit="in",dpi=200)
    
    tab <- gt(deltaOut)
    #tab <- data_color(tab,1,"red")
    tab <- tab_options(tab,column_labels.font.weight = "bold")
    print(tab)
    # docx format
    gtsave(tab,paste(output_folder,"/DeltaChangeFactorsT_",scenario,period,".docx", sep=''))
    
    ############################################################
    # APPLY DELTA CHANGES TO OBSERVED
    #
    
    rm(Cobs)
    Cobs <- mutate(obs,Month=format(Date,"%m"))

    for (k in 2:ncol(deltaOut)){
      dch <- deltaOut[,c(1,k)]
      modName <- names(dch)[2]
      print(paste("Delta for:", scenario, period, modName))
      Cobs <- mutate(Cobs,!! modName := NA)
      col <- ncol(Cobs)
      for (i in 1:length(Cobs$Obs)){
        mnd <- as.numeric(Cobs[i,"Month"])
        dc <- as.numeric(dch[mnd,2])
        Cobs[i,col] <- Cobs[i,"Obs"] + dc
      }
    }
    # Do not include the observed values column neither the monthly correction column in the final tabular 
    # (no need to create pine input files with these data sets).
    Cobs <- Cobs %>% dplyr::select(-Obs, -Month)
    write.csv(Cobs,paste0(output_folder,'/temp_deltachanged_', scenario,period,".csv"))
  }
}

#
# Plot the results from the qmap scaling
#

sqrtquant <- function(x,qstep=0.01){
  qq <- quantile(x,prob=seq(0,1,by=qstep))
  sqrt(qq)
}
# plot(sqrtquant(obs$Temp),sqrtquant(hist$Temperature), ylab="Obs Temp (mm)", xlab="Hist T (mm)")
# lines(sqrtquant(hist$Temperature),sqrtquant(hist$corrT),col="red")

