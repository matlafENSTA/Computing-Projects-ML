ptf #
[Hydrological system description]  : "D:\Klima\PINE\PineProj\Myglevatn\cal1\cal1.top"  Thu Jul 18 13:30:04 2024


[respmodules]
module: PINE_IDTilsig
  node: 1   IDTilsig_node
      Method: IDTILSIG:
              ELEV0           [       masl ] #ELEV0   # : Lower elevation of zone 1
              ELEV1           [       masl ] #ELEV1   # : Lower elevation of zone 2
              ELEV2           [       masl ] #ELEV2   # : Lower elevation of zone 3
              ELEV3           [       masl ] #ELEV3   # : Lower elevation of zone 4
              ELEV4           [       masl ] #ELEV4   # : Lower elevation of zone 5
              ELEV5           [       masl ] #ELEV5   # : Lower elevation of zone 6
              ELEV6           [       masl ] #ELEV6   # : Lower elevation of zone 7
              ELEV7           [       masl ] #ELEV7   # : Lower elevation of zone 8
              ELEV8           [       masl ] #ELEV8   # : Lower elevation of zone 9
              ELEV9           [       masl ] #ELEV9   # : Lower elevation of zone 10
              ELEV10          [       masl ] #ELEV10  # : Upper elevation of zone 10
              NEDNIV          [          - ] #NEDNIV  # : Number of forested elevation zones
              FLDAREA         [        km2 ] #FLDAREA # : Catchment area
              AREA1           [        km2 ] #AREA1   # : Area of zone 1
              AREA2           [        km2 ] #AREA2   # : Area of zone 2
              AREA3           [        km2 ] #AREA3   # : Area of zone 3
              AREA4           [        km2 ] #AREA4   # : Area of zone 4
              AREA5           [        km2 ] #AREA5   # : Area of zone 5
              AREA6           [        km2 ] #AREA6   # : Area of zone 6
              AREA7           [        km2 ] #AREA7   # : Area of zone 7
              AREA8           [        km2 ] #AREA8   # : Area of zone 8
              AREA9           [        km2 ] #AREA9   # : Area of zone 9
              AREA10          [        km2 ] #AREA10  # : Area of zone 10
              ELEVTMP         [       masl ] #ELEVTMP # : Ref.height temperature
              ELEVPRC         [       masl ] #ELEVPRC # : Ref.height precipitation
              GLAC1           [          - ] #GLAC1   # : Glacier portion of zone 1
              GLAC2           [          - ] #GLAC2   # : Glacier portion of zone 2
              GLAC3           [          - ] #GLAC3   # : Glacier portion of zone 3
              GLAC4           [          - ] #GLAC4   # : Glacier portion of zone 4
              GLAC5           [          - ] #GLAC5   # : Glacier portion of zone 5
              GLAC6           [          - ] #GLAC6   # : Glacier portion of zone 6
              GLAC7           [          - ] #GLAC7   # : Glacier portion of zone 7
              GLAC8           [          - ] #GLAC8   # : Glacier portion of zone 8
              GLAC9           [          - ] #GLAC9   # : Glacier portion of zone 9
              GLAC10          [          - ] #GLAC10  # : Glacier portion of zone 10
              S100            [          - ] #S100    # : Upper limit in snow distribution in upper forest free part
              S75             [          - ] #S75     # : Upper quartile in snow distribution in upper forest free part
              S50             [          - ] #S50     # : Median in snow distribution function in upper forest free part
              S25             [          - ] #S25     # : Lower quartile in snow distribution in upper forest free part
              S00             [          - ] #S00     # : Lower limit in snow distribution in upper forest free part
              SL100           [          - ] #SL100   # : Upper limit in snow distribution in lower forested part
              SL75            [          - ] #SL75    # : Upper quartile in snow distribution in lower forested part
              SL50            [          - ] #SL50    # : Median in snow distribution function in lower forested part
              SL25            [          - ] #SL25    # : Lower quartile in snow distribution in lower forested part
              SL00            [          - ] #SL00    # : Lower limit in snow distribution in lower forested part
              MAXUNIFSNW      [         mm ] #MAXUNIF # : Lower limit in snow distribution in lower forested part
              RCORR           [          - ] #RCOR    # : Point correction factor for rain
              SCORR           [          - ] #SCOR    # : Extra point correction factor for snow
              TX              [         °C ] #TX      # : Threshold temperature rain - snow
              TCGRAD          [    °C/100m ] #TCGR    # : Temperature lapse rate on clear days
              TPGRAD          [    °C/100m ] #TPGR    # : Temperature lapse rate on overcast days
              PGRAD           [     -/100m ] #PGRD    # : Precipitation height gradient
              CX              [     mm/d°C ] #CX      # : Degree-day factor for snow melt in upper forest free part
              CXN             [     mm/d°C ] #CXN     # : Degree-day factor for snow melt in lower forested part
              TS              [         °C ] #TS      # : Treshold melt/freeze in upper forest free part
              TSN             [         °C ] #TSN     # : Treshold melt/freeze in lower forested part
              CFR             [     mm/d°C ] #CFR     # : Refreeze coefficient
              LW              [          - ] #LW      # : Max relative portion liquid water in snow
              NDAG            [          - ] #NDAG    # : Day no. for snow to ice conversion
              CBRE            [          - ] #CBRE    # : Adjustment of CX for glacier melting
              EPJAN           [         mm ] #EPJAN   # : Daily potential evapotranspiration in January
              EPFEB           [         mm ] #EPFEB   # : Daily potential evapotranspiration in February
              EPMAR           [         mm ] #EPMAR   # : Daily potential evapotranspiration in March
              EPAPR           [         mm ] #EPAPR   # : Daily potential evapotranspiration in April
              EPMAY           [         mm ] #EPMAY   # : Daily potential evapotranspiration in May
              EPJUN           [         mm ] #EPJUN   # : Daily potential evapotranspiration in June
              EPJUL           [         mm ] #EPJUL   # : Daily potential evapotranspiration in July
              EPAUG           [         mm ] #EPAUG   # : Daily potential evapotranspiration in August
              EPSEP           [         mm ] #EPSEP   # : Daily potential evapotranspiration in September
              EPOKT           [         mm ] #EPOCT   # : Daily potential evapotranspiration in October
              EPNOV           [         mm ] #EPNOV   # : Daily potential evapotranspiration in November
              EPDES           [         mm ] #EPDEC   # : Daily potential evapotranspiration in December
              FC              [         mm ] #FC      # : Field capacity
              FCDEL           [          - ] #FCDEL   # : Minimum soil moisture filling for potential evapotranspiration
              BETA            [          - ] #BETA    # : Non-linearity in soil water retention
              INFMAX          [       mm/h ] #INFMAX  # : Infiltration capacity
              KUZ2            [      1/day ] #KUZ2    # : Outlet coefficient for quickest surface runoff
              KUZ1            [      1/day ] #KUZ1    # : Outlet coefficient for quick surface runoff
              KUZ             [      1/day ] #KUZ     # : Outlet coefficient for slow surface runoff
              KLZ             [      1/day ] #KLZ     # : Outlet coefficient for groundwater runoff
              PERC            [     mm/day ] #PERC    # : Constant percolation rate to groundwater storage
              UZ2             [         mm ] #UZ2     # : Threshold between quickest and quick surface runoff
              UZ1             [         mm ] #UZ1     # : Threshold between quick and slow surface runoff
              SJODEL          [          % ] #SJODEL  # : Lake portion
              MAGDEL          [          % ] #MAGDEL  # : Reservoir portion
              OPPTR           [     mm/day ] #OPPTR   # : Draw-up from ground water to root zone
              EQLAKEAREA      [        km2 ] #EQLAKEA # : Equivalent lake area
              RATINGCNST      [          - ] #RATCNST # : Rating curve constant
              RATINGZERO      [          m ] #RATZERO # : Rating curve zero (m)  (stage for zero flow)
              RATINGEXP       [          - ] #RATEXPN # : Rating curve exponent
              RATINGEXP       [          - ] #CATFACT # : Portion of catchment that drains through lake
              PARDSP          [          - ]      2.000 : Display parameters after model run [1=yes, 0=no]
              SNWDRY1         [         mm ] #ISNW1   # : Initial dry snow pack in zone 1
              SNWDRY2         [         mm ] #ISNW2   # : Initial dry snow pack in zone 1
              SNWDRY3         [         mm ] #ISNW3   # : Initial dry snow pack in zone 1
              SNWDRY4         [         mm ] #ISNW4   # : Initial dry snow pack in zone 1
              SNWDRY5         [         mm ] #ISNW5   # : Initial dry snow pack in zone 1
              SNWDRY6         [         mm ] #ISNW6   # : Initial dry snow pack in zone 1
              SNWDRY7         [         mm ] #ISNW7   # : Initial dry snow pack in zone 1
              SNWDRY8         [         mm ] #ISNW8   # : Initial dry snow pack in zone 1
              SNWDRY9         [         mm ] #ISNW9   # : Initial dry snow pack in zone 1
              SNWDRY10        [         mm ] #ISNW10  # : Initial dry snow pack in zone 1
              SNWWET1         [         mm ] #ISWW1   # : Initial snow water in zone 1
              SNWWET2         [         mm ] #ISWW2   # : Initial snow water in zone 1
              SNWWET3         [         mm ] #ISWW3   # : Initial snow water in zone 1
              SNWWET4         [         mm ] #ISWW4   # : Initial snow water in zone 1
              SNWWET5         [         mm ] #ISWW5   # : Initial snow water in zone 1
              SNWWET6         [         mm ] #ISWW6   # : Initial snow water in zone 1
              SNWWET7         [         mm ] #ISWW7   # : Initial snow water in zone 1
              SNWWET8         [         mm ] #ISWW8   # : Initial snow water in zone 1
              SNWWET9         [         mm ] #ISWW9   # : Initial snow water in zone 1
              SNWWET10        [         mm ] #ISWW10  # : Initial snow water in zone 1
              SNWCOV1         [         mm ] #ISCOV1  # : Initial snow cover in zone 1
              SNWCOV2         [         mm ] #ISCOV2  # : Initial snow cover in zone 1
              SNWCOV3         [         mm ] #ISCOV3  # : Initial snow cover in zone 1
              SNWCOV4         [         mm ] #ISCOV4  # : Initial snow cover in zone 1
              SNWCOV5         [         mm ] #ISCOV5  # : Initial snow cover in zone 1
              SNWCOV6         [         mm ] #ISCOV6  # : Initial snow cover in zone 1
              SNWCOV7         [         mm ] #ISCOV7  # : Initial snow cover in zone 1
              SNWCOV8         [         mm ] #ISCOV8  # : Initial snow cover in zone 1
              SNWCOV9         [         mm ] #ISCOV9  # : Initial snow cover in zone 1
              SNWCOV10        [         mm ] #ISCOV10 # : Initial snow cover in zone 1
              ICEBAL1         [         mm ] #IBAL1   # : Initial glacier balance in zone 1
              ICEBAL2         [         mm ] #IBAL2   # : Initial glacier balance in zone 1
              ICEBAL3         [         mm ] #IBAL3   # : Initial glacier balance in zone 1
              ICEBAL4         [         mm ] #IBAL4   # : Initial glacier balance in zone 1
              ICEBAL5         [         mm ] #IBAL5   # : Initial glacier balance in zone 1
              ICEBAL6         [         mm ] #IBAL6   # : Initial glacier balance in zone 1
              ICEBAL7         [         mm ] #IBAL7   # : Initial glacier balance in zone 1
              ICEBAL8         [         mm ] #IBAL8   # : Initial glacier balance in zone 1
              ICEBAL9         [         mm ] #IBAL9   # : Initial glacier balance in zone 1
              ICEBAL10        [         mm ] #IBAL10  # : Initial glacier balance in zone 1
              SOLMST          [         mm ] #ISM     # : Initial soil moisture content
              SRFWATSTOR      [         mm ] #IUZ     # : Initial surface water content
              GRNDWATSTOR     [         mm ] #ILZ     # : Initial ground water content
              LAKELEV         [         mm ] #LAKLEV  # : Initial lake level
      Method: VOID
      Transf: 0	        VOID      0
      Inlets: 0         VOID      0
      RInlts: 0         VOID      0
      XtvRef: VOID      0         VOID
      Direct: SIMRUNOFF
module: VOID
[inputvars]
    node: 1  INPRECIP = PRECIP
    node: 1  INMAXTMP = TMAX
    node: 1  INMINTMP = TMIN
    node: 1  INOBSVAR = OBSVAR
    node: 0  VOID = VOID
[outputvars]
    node: 1   OBSRUNOFF
    node: 1   SIMRUNOFF
    node: 0     VOID
[erroranals]
          R2        1  SIMRUNOFF              1  OBSRUNOFF
          ACCDIFF   1  SIMRUNOFF              1  OBSRUNOFF
          CHANGE    1  FLDICEBAL
      VOID
[endsystem]
CXN_DCX         [          - ] #CXN_DCX # : Ratio between CXN and CX (help variable used during autocalibration)
KUZ2_D1         [          - ] #KUZ2_D1 # : Ratio between KUZ2 and KUZ1 (help variable used during autocalibration)
KUZ1_D0         [          - ] #KUZ1_D0 # : Ratio between KUZ1 and KUZ (help variable used during autocalibration)
KUZ_DL          [          - ] #KUZ_DL  # : Ratio between KUZ and KLZ (help variable used during autocalibration)
UZ2_D1          [         mm ] #UZ2_D1  # : Ratio between UZ2 and UZ1 (help variable used during autocalibration)
PARDSP          [          - ] #PARDSP  # : Display parameters after model run [1=yes, 0=no]
