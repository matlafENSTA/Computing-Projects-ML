setwd("F:/MyNTNU/NTE_Project/MyStudy/NamsvatnGaugeDaily")
getwd()

#Read the reference flow file
df_flow = read.delim("./SimulatedFlow/NamsvatnDaily_RefFlow_Verify.txt", header = T, sep="\t")

#Read 100 Montecarlo parameters
datfile <- "./ParameterAnalysis/MonteCarlo/NamsvatnGaugeDailyMC100.txt"
pardata <- read.delim(datfile, header = T, sep = "\t" )
N<- nrow(pardata)


df_result<- pardata
row.names(df_result)<- 1:100

newcols<-c("VerifyR2","VerifyACCD")
df_result[,newcols]<- NA


for (i in 1:N){
  #i = 1
  print(i)
  RCORR = pardata[i,"PCOR"]
  SCORR = pardata[i,"SCOR"]
  TX = pardata[i,"TX"]
  CX = pardata[i,"CX"]
  CXN = pardata[i,"CXN"]
  TS = pardata[i,"TS"]
  TSN = pardata[i,"TSN"]
  FC = pardata[i,"FC"]
  BETA = pardata[i,"BETA"]
  FCDEL = pardata[i,"FCDEL"]
  KUZ2 = pardata[i,"KUZ2"]
  KUZ1 = pardata[i,"KUZ1"]
  KUZ = pardata[i,"KUZ"]
  KLZ = pardata[i,"KLZ"]
  PERC = pardata[i,"PERC"]
  UZ2 = pardata[i,"UZ2"]
  UZ1 = pardata[i,"UZ1"]
  
  #initial Soil Moisture = 0.8 * Field Capacity
  SM = FC*0.8
  
  # Type measurement elevation for Pecipitation and Temperature
  EL_PT = 450.0
  
  
  file_catchment= "NamsvatnGaugeDaily"
  
  filename<-".\\ParameterFile\\Parameter.top"
  sink(file=filename)

  cat(sprintf("[Hydrological system description]  : \"%s\"  Tue Dec 04 17:54:23 2018\n", filename))
  cat("\n")
  cat("\n")
  cat(sprintf("[respmodules]\n"))
  cat(sprintf("module: PINE_IDTilsig\n"))
  cat(sprintf("  node: 1   IDTilsig_node\n"))
  cat(sprintf("      Method: IDTILSIG:\n"))
  
  cat(sprintf("              ELEV0           [       masl ]      439.000 : Lower elevation of zone 1\n"))
  cat(sprintf("              ELEV1           [       masl ]      486.000 : Lower elevation of zone 2\n"))
  cat(sprintf("              ELEV2           [       masl ]      559.000 : Lower elevation of zone 3\n"))
  cat(sprintf("              ELEV3           [       masl ]      620.000 : Lower elevation of zone 4\n"))
  cat(sprintf("              ELEV4           [       masl ]      689.000 : Lower elevation of zone 5\n"))
  cat(sprintf("              ELEV5           [       masl ]      743.000 : Lower elevation of zone 6\n"))
  cat(sprintf("              ELEV6           [       masl ]      813.000 : Lower elevation of zone 7\n"))
  cat(sprintf("              ELEV7           [       masl ]      876.000 : Lower elevation of zone 8\n"))
  cat(sprintf("              ELEV8           [       masl ]      949.000 : Lower elevation of zone 9\n"))
  cat(sprintf("              ELEV9           [       masl ]     1036.000 : Lower elevation of zone 10\n"))
  cat(sprintf("              ELEV10          [       masl ]     1675.000 : Upper elevation of zone 10\n"))
  cat(sprintf("              NEDNIV          [          - ]            2 : Number of forested elevation zones\n"))
  cat(sprintf("              FLDAREA         [        km2 ]      701.500 : Catchment area\n"))
  cat(sprintf("              AREA1           [        km2 ]       70.150 : Area of zone 1\n"))
  cat(sprintf("              AREA2           [        km2 ]       70.150 : Area of zone 2\n"))
  cat(sprintf("              AREA3           [        km2 ]       70.150 : Area of zone 3\n"))
  cat(sprintf("              AREA4           [        km2 ]       70.150 : Area of zone 4\n"))
  cat(sprintf("              AREA5           [        km2 ]       70.150 : Area of zone 5\n"))
  cat(sprintf("              AREA6           [        km2 ]       70.150 : Area of zone 6\n"))
  cat(sprintf("              AREA7           [        km2 ]       70.150 : Area of zone 7\n"))
  cat(sprintf("              AREA8           [        km2 ]       70.150 : Area of zone 8\n"))
  cat(sprintf("              AREA9           [        km2 ]       70.150 : Area of zone 9\n"))
  cat(sprintf("              AREA10          [        km2 ]       70.150 : Area of zone 10\n"))
  
  cat(sprintf("              ELEVTMP         [       masl ]    %9.3f : Ref.height temperature\n", EL_PT))
  cat(sprintf("              ELEVPRC         [       masl ]    %9.3f : Ref.height precipitation\n", EL_PT))
  
  cat(sprintf("              GLAC1           [          - ]        0.140 : Glacier portion of zone 1\n"))
  cat(sprintf("              GLAC2           [          - ]        0.140 : Glacier portion of zone 2\n"))
  cat(sprintf("              GLAC3           [          - ]        0.140 : Glacier portion of zone 3\n"))
  cat(sprintf("              GLAC4           [          - ]        0.140 : Glacier portion of zone 4\n"))
  cat(sprintf("              GLAC5           [          - ]        0.140 : Glacier portion of zone 5\n"))
  cat(sprintf("              GLAC6           [          - ]        0.140 : Glacier portion of zone 6\n"))
  cat(sprintf("              GLAC7           [          - ]        0.140 : Glacier portion of zone 7\n"))
  cat(sprintf("              GLAC8           [          - ]        0.140 : Glacier portion of zone 8\n"))
  cat(sprintf("              GLAC9           [          - ]        0.140 : Glacier portion of zone 9\n"))
  cat(sprintf("              GLAC10          [          - ]        0.140 : Glacier portion of zone 10\n"))
  cat(sprintf("              S100            [          - ]        2.000 : Upper limit in snow distribution in upper forest free part\n"))
  cat(sprintf("              S75             [          - ]        1.500 : Upper quartile in snow distribution in upper forest free part\n"))
  cat(sprintf("              S50             [          - ]        1.000 : Median in snow distribution function in upper forest free part\n"))
  cat(sprintf("              S25             [          - ]        0.500 : Lower quartile in snow distribution in upper forest free part\n"))
  cat(sprintf("              S00             [          - ]        0.000 : Lower limit in snow distribution in upper forest free part\n"))
  cat(sprintf("              SL100           [          - ]        1.500 : Upper limit in snow distribution in lower forested part\n"))
  cat(sprintf("              SL75            [          - ]        1.250 : Upper quartile in snow distribution in lower forested part\n"))
  cat(sprintf("              SL50            [          - ]        1.000 : Median in snow distribution function in lower forested part\n"))
  cat(sprintf("              SL25            [          - ]        0.750 : Lower quartile in snow distribution in lower forested part\n"))
  cat(sprintf("              SL00            [          - ]        0.500 : Lower limit in snow distribution in lower forested part\n"))
  cat(sprintf("              MAXUNIFSNW      [         mm ]       20.000 : Threshold for redistribution of snow\n"))
  
  cat(sprintf("              RCORR           [          - ]    %9.3f : Point correction factor for rain\n", RCORR))
  cat(sprintf("              SCORR           [          - ]    %9.3f : Extra point correction factor for snow\n", SCORR))
  cat(sprintf("              TX              [         °C ]    %9.3f : Threshold temperature rain - snow\n", TX))
  
  cat(sprintf("              TCGRAD          [    °C/100m ]       -1.000 : Temperature lapse rate on clear days\n"))
  cat(sprintf("              TPGRAD          [    °C/100m ]       -0.500 : Temperature lapse rate on overcast days\n"))
  cat(sprintf("              PGRAD           [     -/100m ]        0.050 : Precipitation height gradient\n"))
  
  cat(sprintf("              CX              [     mm/d°C ]    %9.3f : Degree-day factor for snow melt in upper forest free part\n", CX))
  cat(sprintf("              CXN             [     mm/d°C ]    %9.3f : Degree-day factor for snow melt in lower forested part\n", CXN))
  cat(sprintf("              TS              [         °C ]    %9.3f : Treshold melt/freeze in upper forest free part\n", TS))
  cat(sprintf("              TSN             [         °C ]    %9.3f : Treshold melt/freeze in lower forested part\n", TSN))
  
  cat(sprintf("              CFR             [     mm/d°C ]        0.010 : Refreeze coefficient\n"))
  cat(sprintf("              LW              [          - ]        0.070 : Max relative portion liquid water in snow\n"))
  cat(sprintf("              NDAG            [          - ]      270.000 : Day no. for snow to ice conversion\n"))
  cat(sprintf("              CBRE            [          - ]        2.000 : Adjustment of CX for glacier melting\n"))
  
  cat(sprintf("              EPJAN           [         mm ]        0.100 : Daily potential evapotranspiration in January\n"))
  cat(sprintf("              EPFEB           [         mm ]        0.200 : Daily potential evapotranspiration in February\n"))
  cat(sprintf("              EPMAR           [         mm ]        0.700 : Daily potential evapotranspiration in March\n"))
  cat(sprintf("              EPAPR           [         mm ]        1.000 : Daily potential evapotranspiration in April\n"))
  cat(sprintf("              EPMAY           [         mm ]        2.300 : Daily potential evapotranspiration in May\n"))
  cat(sprintf("              EPJUN           [         mm ]        3.500 : Daily potential evapotranspiration in June\n"))
  cat(sprintf("              EPJUL           [         mm ]        3.500 : Daily potential evapotranspiration in July\n"))
  cat(sprintf("              EPAUG           [         mm ]        2.300 : Daily potential evapotranspiration in August\n"))
  cat(sprintf("              EPSEP           [         mm ]        1.000 : Daily potential evapotranspiration in September\n"))
  cat(sprintf("              EPOKT           [         mm ]        0.700 : Daily potential evapotranspiration in October\n"))
  cat(sprintf("              EPNOV           [         mm ]        0.200 : Daily potential evapotranspiration in November\n"))
  cat(sprintf("              EPDES           [         mm ]        0.100 : Daily potential evapotranspiration in December\n"))
  
  cat(sprintf("              FC              [         mm ]    %9.3f : Field capacity\n", FC))
  cat(sprintf("              FCDEL           [          - ]    %9.3f : Minimum soil moisture filling for potential evapotranspiration\n", FCDEL))
  cat(sprintf("              BETA            [          - ]    %9.3f : Non-linearity in soil water retention\n", BETA))
  
  cat(sprintf("              INFMAX          [       mm/h ]       50.000 : Infiltration capacity\n"))
  
  cat(sprintf("              KUZ2            [      1/day ]    %9.3f : Outlet coefficient for quickest surface runoff\n", KUZ2))
  cat(sprintf("              KUZ1            [      1/day ]    %9.3f : Outlet coefficient for quick surface runoff\n", KUZ1))
  cat(sprintf("              KUZ             [      1/day ]    %9.3f : Outlet coefficient for slow surface runoff\n", KUZ))
  cat(sprintf("              KLZ             [      1/day ]    %9.3f : Outlet coefficient for groundwater runoff\n", KLZ))
  cat(sprintf("              PERC            [     mm/day ]    %9.3f : Constant percolation rate to groundwater storage\n", PERC))
  cat(sprintf("              UZ2             [         mm ]    %9.3f : Threshold between quickest and quick surface runoff\n", UZ2))
  cat(sprintf("              UZ1             [         mm ]    %9.3f : Threshold between quick and slow surface runoff\n", UZ1))
  
  cat(sprintf("              SJODEL          [          %% ]       12.500 : Lake portion\n"))
  cat(sprintf("              MAGDEL          [          %% ]        0.000 : Reservoir portion\n"))
  cat(sprintf("              OPPTR           [     mm/day ]        1.000 : Draw-up from ground water to root zone\n"))
  cat(sprintf("              EQLAKEAREA      [        km2 ]        0.000 : Equivalent lake area\n"))
  cat(sprintf("              RATINGCNST      [          - ]        0.000 : Rating curve constant\n"))
  cat(sprintf("              RATINGZERO      [          m ]        0.000 : Rating curve zero (m)  (stage for zero flow)\n"))
  cat(sprintf("              RATINGEXP       [          - ]        1.500 : Rating curve exponent\n"))
  cat(sprintf("              LAKECATCHFACT   [          - ]        0.000 : Portion of catchment that drains through lake\n"))
  cat(sprintf("              PARDSP          [          - ]            0 : Display parameteres after model runs [1=yes, 0=no]\n"))
  cat(sprintf("              SNWDRY1         [         mm ]        0.000 : Initial dry snow pack in zone 1\n"))
  cat(sprintf("              SNWDRY2         [         mm ]        0.000 : Initial dry snow pack in zone 1\n"))
  cat(sprintf("              SNWDRY3         [         mm ]        0.000 : Initial dry snow pack in zone 1\n"))
  cat(sprintf("              SNWDRY4         [         mm ]        0.000 : Initial dry snow pack in zone 1\n"))
  cat(sprintf("              SNWDRY5         [         mm ]        0.000 : Initial dry snow pack in zone 1\n"))
  cat(sprintf("              SNWDRY6         [         mm ]        0.000 : Initial dry snow pack in zone 1\n"))
  cat(sprintf("              SNWDRY7         [         mm ]        0.000 : Initial dry snow pack in zone 1\n"))
  cat(sprintf("              SNWDRY8         [         mm ]        0.000 : Initial dry snow pack in zone 1\n"))
  cat(sprintf("              SNWDRY9         [         mm ]        0.000 : Initial dry snow pack in zone 1\n"))
  cat(sprintf("              SNWDRY10        [         mm ]        0.000 : Initial dry snow pack in zone 1\n"))
  cat(sprintf("              SNWWET1         [         mm ]        0.000 : Initial snow water in zone 1\n"))
  cat(sprintf("              SNWWET2         [         mm ]        0.000 : Initial snow water in zone 1\n"))
  
  
  cat(sprintf("              SNWWET3         [         mm ]        0.000 : Initial snow water in zone 1\n"))
  cat(sprintf("              SNWWET4         [         mm ]        0.000 : Initial snow water in zone 1\n"))
  cat(sprintf("              SNWWET5         [         mm ]        0.000 : Initial snow water in zone 1\n"))
  cat(sprintf("              SNWWET6         [         mm ]        0.000 : Initial snow water in zone 1\n"))
  cat(sprintf("              SNWWET7         [         mm ]        0.000 : Initial snow water in zone 1\n"))
  cat(sprintf("              SNWWET8         [         mm ]        0.000 : Initial snow water in zone 1\n"))
  cat(sprintf("              SNWWET9         [         mm ]        0.000 : Initial snow water in zone 1\n"))
  cat(sprintf("              SNWWET10        [         mm ]        0.000 : Initial snow water in zone 1\n"))
  cat(sprintf("              SNWCOV1         [         mm ]        0.000 : Initial snow cover in zone 1\n"))
  cat(sprintf("              SNWCOV2         [         mm ]        0.000 : Initial snow cover in zone 1\n"))
  cat(sprintf("              SNWCOV3         [         mm ]        0.000 : Initial snow cover in zone 1\n"))
  cat(sprintf("              SNWCOV4         [         mm ]        0.000 : Initial snow cover in zone 1\n"))
  cat(sprintf("              SNWCOV5         [         mm ]        0.000 : Initial snow cover in zone 1\n"))
  cat(sprintf("              SNWCOV6         [         mm ]        0.000 : Initial snow cover in zone 1\n"))
  cat(sprintf("              SNWCOV7         [         mm ]        0.000 : Initial snow cover in zone 1\n"))
  cat(sprintf("              SNWCOV8         [         mm ]        0.000 : Initial snow cover in zone 1\n"))
  cat(sprintf("              SNWCOV9         [         mm ]        0.000 : Initial snow cover in zone 1\n"))
  cat(sprintf("              SNWCOV10        [         mm ]        0.000 : Initial snow cover in zone 1\n"))
  cat(sprintf("              ICEBAL1         [         mm ]        0.000 : Initial glacier balance in zone 1\n"))
  cat(sprintf("              ICEBAL2         [         mm ]        0.000 : Initial glacier balance in zone 1\n"))
  cat(sprintf("              ICEBAL3         [         mm ]        0.000 : Initial glacier balance in zone 1\n"))
  cat(sprintf("              ICEBAL4         [         mm ]        0.000 : Initial glacier balance in zone 1\n"))
  cat(sprintf("              ICEBAL5         [         mm ]        0.000 : Initial glacier balance in zone 1\n"))
  cat(sprintf("              ICEBAL6         [         mm ]        0.000 : Initial glacier balance in zone 1\n"))
  cat(sprintf("              ICEBAL7         [         mm ]        0.000 : Initial glacier balance in zone 1\n"))
  cat(sprintf("              ICEBAL8         [         mm ]        0.000 : Initial glacier balance in zone 1\n"))
  cat(sprintf("              ICEBAL9         [         mm ]        0.000 : Initial glacier balance in zone 1\n"))
  cat(sprintf("              ICEBAL10        [         mm ]        0.000 : Initial glacier balance in zone 1\n"))

  cat(sprintf("              SOLMST          [         mm ]    %9.3f : Initial soil moisture content\n", SM))
  cat(sprintf("              SRFWATSTOR      [         mm ]        0.000 : Initial surface water content\n"))
  cat(sprintf("              GRNDWATSTOR     [         mm ]       20.000 : Initial ground water content\n"))
  cat(sprintf("              LAKELEV         [          m ]        0.000 : Initial lake water level\n"))
  
  cat(sprintf("      Method: VOID\n"))
  cat(sprintf("      Transf: 0	        VOID      0\n"))
  cat(sprintf("      Inlets: 0         VOID      0\n"))
  cat(sprintf("      RInlts: 0         VOID      0\n"))
  cat(sprintf("      XtvRef: VOID      0         VOID\n"))
  cat(sprintf("      Direct: SIMRUNOFF\n"))
  cat(sprintf("module: VOID\n"))
  cat(sprintf("[inputvars]\n"))
  cat(sprintf("    node: 1  INPRECIP = PRECIP\n"))
  cat(sprintf("    node: 1  INMAXTMP = TMAX\n"))
  cat(sprintf("    node: 1  INMINTMP = TMIN\n"))
  cat(sprintf("    node: 1  INOBSVAR = OBSVAR\n"))
  cat(sprintf("    node: 0  VOID = VOID\n"))
  cat(sprintf("[outputvars]\n"))
  cat(sprintf("        node: 1     OBSRUNOFF\n"))
  cat(sprintf("        node: 1     SIMRUNOFF\n"))
  cat(sprintf("        node: 1     INOBSVAR\n"))
  cat(sprintf("        node: 1     ACCDIFFERENS\n"))
  cat(sprintf("        node: 1     INPRECIP\n"))
  cat(sprintf("        node: 1     INMAXTMP\n"))
  cat(sprintf("        node: 1     INMINTMP\n"))
  cat(sprintf("        node: 1     FLDPRECIP\n"))
  cat(sprintf("        node: 1     ACCFLDPRECIP\n"))
  cat(sprintf("        node: 1     FLDTEMPER\n"))
  cat(sprintf("        node: 1     FLDSNWPCK\n"))
  cat(sprintf("        node: 1     FLDSNWWET\n"))
  cat(sprintf("        node: 1     FLDSNWCOV\n"))
  cat(sprintf("        node: 1     FLDSNWMLT\n"))
  cat(sprintf("        node: 1     FLDSNWOUT\n"))
  cat(sprintf("        node: 1     FLDICEBAL\n"))
  cat(sprintf("        node: 1     FLDICEOUT\n"))
  cat(sprintf("        node: 1     EPOT\n"))
  cat(sprintf("        node: 1     AET\n"))
  cat(sprintf("        node: 1     FLDEVAPOT\n"))
  cat(sprintf("        node: 1     ACCFLDEVAPOT\n"))
  cat(sprintf("        node: 1     SOLMST\n"))
  cat(sprintf("        node: 1     SRFRUN\n"))
  cat(sprintf("        node: 1     SOLOUT\n"))
  cat(sprintf("        node: 1     AVGRUNOFF\n"))
  cat(sprintf("        node: 1     MOMSIMUPR\n"))
  cat(sprintf("        node: 1     MOMSIMLWR\n"))
  cat(sprintf("        node: 1     MOMSIMRUN\n"))
  cat(sprintf("        node: 1     ACCOBSRUNOFF\n"))
  cat(sprintf("        node: 1     ACCSIMRUNOFF\n"))
  cat(sprintf("        node: 1     SRFWATSTOR\n"))
  cat(sprintf("        node: 1     GRNDWATSTOR\n"))
  cat(sprintf("        node: 1     WATERSTORED\n"))
  cat(sprintf("        node: 1     ACTUALPERCOL\n"))
  cat(sprintf("        node: 1     ACTUALDRAWUP\n"))
  cat(sprintf("        node: 1     LAKELEV\n"))
  cat(sprintf("        node: 1     LAKECONT\n"))
  cat(sprintf("        node: 1     SNWDRY1\n"))
  cat(sprintf("        node: 1     SNWDRY2\n"))
  cat(sprintf("        node: 1     SNWDRY3\n"))
  cat(sprintf("        node: 1     SNWDRY4\n"))
  cat(sprintf("        node: 1     SNWDRY5\n"))
  cat(sprintf("        node: 1     SNWDRY6\n"))
  cat(sprintf("        node: 1     SNWDRY7\n"))
  cat(sprintf("        node: 1     SNWDRY8\n"))
  cat(sprintf("        node: 1     SNWDRY9\n"))
  cat(sprintf("        node: 1     SNWDRY10\n"))
  cat(sprintf("        node: 1     SNWWET1\n"))
  cat(sprintf("        node: 1     SNWWET2\n"))
  cat(sprintf("        node: 1     SNWWET3\n"))
  cat(sprintf("        node: 1     SNWWET4\n"))
  cat(sprintf("        node: 1     SNWWET5\n"))
  cat(sprintf("        node: 1     SNWWET6\n"))
  cat(sprintf("        node: 1     SNWWET7\n"))
  cat(sprintf("        node: 1     SNWWET8\n"))
  cat(sprintf("        node: 1     SNWWET9\n"))
  cat(sprintf("        node: 1     SNWWET10\n"))
  cat(sprintf("        node: 1     SNWPCK1\n"))
  cat(sprintf("        node: 1     SNWPCK2\n"))
  cat(sprintf("        node: 1     SNWPCK3\n"))
  cat(sprintf("        node: 1     SNWPCK4\n"))
  cat(sprintf("        node: 1     SNWPCK5\n"))
  cat(sprintf("        node: 1     SNWPCK6\n"))
  cat(sprintf("        node: 1     SNWPCK7\n"))
  cat(sprintf("        node: 1     SNWPCK8\n"))
  cat(sprintf("        node: 1     SNWPCK9\n"))
  cat(sprintf("        node: 1     SNWPCK10\n"))
  cat(sprintf("        node: 1     SNWCOV1\n"))
  cat(sprintf("        node: 1     SNWCOV2\n"))
  cat(sprintf("        node: 1     SNWCOV3\n"))
  cat(sprintf("        node: 1     SNWCOV4\n"))
  cat(sprintf("        node: 1     SNWCOV5\n"))
  cat(sprintf("        node: 1     SNWCOV6\n"))
  cat(sprintf("        node: 1     SNWCOV7\n"))
  cat(sprintf("        node: 1     SNWCOV8\n"))
  cat(sprintf("        node: 1     SNWCOV9\n"))
  cat(sprintf("        node: 1     SNWCOV10\n"))
  cat(sprintf("        node: 1     ICEBAL1\n"))
  cat(sprintf("        node: 1     ICEBAL2\n"))
  cat(sprintf("        node: 1     ICEBAL3\n"))
  cat(sprintf("        node: 1     ICEBAL4\n"))
  cat(sprintf("        node: 1     ICEBAL5\n"))
  cat(sprintf("        node: 1     ICEBAL6\n"))
  cat(sprintf("        node: 1     ICEBAL7\n"))
  cat(sprintf("        node: 1     ICEBAL8\n"))
  cat(sprintf("        node: 1     ICEBAL9\n"))
  cat(sprintf("        node: 1     ICEBAL10\n"))
  cat(sprintf("        node: 1     INOBSVAR\n"))
  cat(sprintf("        node: 1     SIMDSCHRG\n"))
  cat(sprintf("    node: 0     VOID\n"))
  cat(sprintf("[erroranals]\n"))
  cat(sprintf("          R2        1  SIMRUNOFF              1  OBSRUNOFF\n"))
  cat(sprintf("          ACCDIFF   1  SIMRUNOFF              1  OBSRUNOFF\n"))
  cat(sprintf("          MAE       1  SIMRUNOFF              1  OBSRUNOFF\n"))
  cat(sprintf("          MAAE      1  SIMRUNOFF              1  OBSRUNOFF\n"))
  cat(sprintf("          RWBE      1  SIMRUNOFF              1  OBSRUNOFF\n"))
  cat(sprintf("          RVBE      1  SIMRUNOFF              1  OBSRUNOFF\n"))
  cat(sprintf("          AVERAGE   1  FLDPRECIP\n"))
  cat(sprintf("          AVERAGE   1  FLDEVAPOT\n"))
  cat(sprintf("          AVERAGE   1  SIMRUNOFF\n"))
  cat(sprintf("          AVERAGE   1  OBSRUNOFF\n"))
  cat(sprintf("          CHANGE    1  WATERSTORED\n"))
  cat(sprintf("          CHANGE    1  FLDICEBAL\n"))
  cat(sprintf("      VOID\n"))
  cat(sprintf("[endsystem]\n"))
  
  sink()
  print("Parameterfile generated")

  
  ##Run simulation with Verification period 
  #Change Input file
  #Type Time period for simulation
  system('hbvpine  /t topol="F:/MyNTNU/NTE_Project/MyStudy/NamsvatnGaugeDaily/ParameterFile/Parameter.top" input="F:/MyNTNU/NTE_Project/MyStudy/NamsvatnGaugeDaily/InputData/NamsvatnObsDaily_1416_inp.dat" output="F:/MyNTNU/NTE_Project/MyStudy/NamsvatnGaugeDaily/SimulatedFlow/NamsvatnGaugeDaily_Flow_Verify.txt" start=01.09.2014  00:00:00 stop=31.08.2016  00:00:00')
  print("SimFlowfile generated")
  
  df_sim <- read.delim("./SimulatedFlow/NamsvatnGaugeDaily_Flow_Verify.txt", header = T, sep="")
  name<- paste("SimFlow",as.character(i), sep ="")
  df_flow[,name]<- df_sim[2:732,"X1SIMDSCHRG"]
  
  
  df_NS<-read.delim("./SimulatedFlow/NamsvatnGaugeDaily_Flow_Verify.inf", skip=1, header = F, sep="")
  PresentR2<- round(df_NS[1,3], digits = 2)
  PresentACCD<-round(df_NS[2,3], digits = 2)
  
  I<-as.character(i)
  df_result[I,][newcols] <- c(PresentR2,PresentACCD)
  
}
print("Loop finished")


#ResultFile<- "./SimulatedFlow/NamsvatnGaugeDaily_100SimFlow_Verify.txt"
#write.table(df_flow, ResultFile, sep ="\t", row.names = TRUE, col.names = TRUE)


write.table(df_result,"./SimulatedFlow/NamsvatnGaugeDaily_100NS_Verify.txt" , sep ="\t", row.names = TRUE, col.names = TRUE)
print("ok")









