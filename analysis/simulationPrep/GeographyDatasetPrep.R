# Build a set of figures showing how MEM works over geography

# Some functions for batch calculating datums and determining lunar nodal cycle
{
  require(lubridate)
  require(tidyverse)
  require(httr)
  require(XML)
  require(TideHarmonics)
  options(tidyverse.quiet = TRUE)
  require(zoo)
  
  # Our transformation function
  scaleFUN <- function(x) sprintf("%.1f", x)
  
  downloadHLData <- function(station_id=8575512, startDate = '2016-01-01 00:00',
                             endDate = '2016-12-31 23:59', datum="NAVD") {
    
    startDate <- ymd_hm(startDate)
    endDate <- ymd_hm(endDate)
    
    queryStart <- startDate
    queryEnd <- as_datetime(ifelse(endDate <= startDate+days(365), 
                                   endDate, startDate+days(365)))
    
    n_iterations <- round(as.numeric(difftime(endDate, startDate, units = "days"))/365+0.5)
    
    for (i in 1:n_iterations) {
      
      queryStartString <- toString(format(queryStart, "%Y%m%d %H:%M"))
      queryEndString <- toString(format(queryEnd, "%Y%m%d %H:%M"))
      
      xml_path <- paste("https://tidesandcurrents.noaa.gov/api/datagetter?product=high_low&application=NOS.COOPS.TAC.WL&begin_date=",
                        queryStartString, 
                        "&end_date=", queryEndString, 
                        "&datum=", datum, 
                        "&station=", toString(station_id), 
                        "&time_zone=GMT&units=metric&format=xml", sep="")
      
      xml_path <- gsub(" ", "%20", xml_path)
      tide_link <- httr::GET(xml_path)
      doc <- xmlParse(tide_link, useInternalNodes = TRUE) ### xmlParse()- is to parse the xml content, the parsed content is stored into doc
      xL <- xmlToList(doc) ###is to convert xml doc into List
      
      if (is.list(xL)) {
        if (exists("observations", where=xL) == T) {
          downloadData <- data.frame(matrix(unlist(xL$observations), ncol=4, byrow=T), 
                                     stringsAsFactors = F)
          names(downloadData) <- c("t", "v", "ty", "f")
          
          if (i == 1) {
            storeData <- downloadData
          } else {
            storeData <- storeData %>%
              bind_rows(downloadData)
          }
        } else {
          stop("Something went wrong with the data availability")
        }
      } else {
        stop("Something went wrong with the query")
      }
      
      queryStart <- queryEnd+minutes(6)
      
      if (i == n_iterations-1) {
        queryEnd <- endDate
      } else if (i < n_iterations-1) {
        queryEnd <- queryStart+days(365*5)
      }
    }
    
    outputData <- storeData %>%
      dplyr::mutate(dateTime = ymd_hm(`t`),
             waterLevel = as.numeric(`v`),
             classifiedTides = ty) %>%
      dplyr::select(dateTime, waterLevel, classifiedTides)
    
    return(outputData)
  }
  
  parseHLdata <- function(hlTable, startDate='2015-01-01 00:00',
                          bufferStart = 15,
                          bufferEnd = 15,
                          endDate='2015-12-31 24:59', graph=F, out_fig_name="temp",
                          gauge_data="biomass/agb_biomass_inundation/NOAA_water_levels") {
    
    classifiedTides <- hlTable %>%
      dplyr::mutate(classifiedTides = str_remove_all(classifiedTides, " "),
             HL = ifelse(classifiedTides %in% c("H", "HH"), "H",
                         ifelse(classifiedTides %in% c("L", "LL"), "L", NA))) %>%
      # filter(! is.na(HL)) %>%
      dplyr::mutate(precedingTideRange = abs(waterLevel-lag(waterLevel)),
             year = year(dateTime),
             month = month(dateTime),
             day = day(dateTime)) %>%
      dplyr::mutate(localMaxTR = suppressWarnings(rollmax(precedingTideRange, 57, "center", na.pad=T)),
             flagMax = ifelse(precedingTideRange == localMaxTR | localMaxTR == lead(precedingTideRange) | localMaxTR == lag(precedingTideRange), "monthlyMax", "not")) %>%
      dplyr::mutate(classifiedTides = ifelse(flagMax=="monthlyMax" & classifiedTides == "HH", "HHS",
                                      ifelse(flagMax=="monthlyMax" & classifiedTides == "LL", "LLS", classifiedTides)),
             timeSpan = difftime(dateTime, lag(dateTime), units="hours")) %>%
      dplyr::mutate(classifiedTides = recode(classifiedTides, "L"="HL", "H"="LH")) %>% 
      dplyr::filter((dateTime >= ymd_hm(startDate, tz="UTC")) & (dateTime <= ymd_hm(endDate, tz="UTC")))
    
    simpleHL <- classifiedTides %>%
      dplyr::group_by(HL) %>%
      dplyr::summarise(wl=mean(waterLevel, na.rm = T),
                n=n()) %>%
      dplyr::mutate(HL = recode(HL, "H"="MHW", "L"="MLW")) %>%
      dplyr::rename(classifiedTides = HL)
    
    datum_summaries <- classifiedTides %>%
      group_by(classifiedTides) %>%
      dplyr::summarise(wl=mean(waterLevel, na.rm = T),
                n = n()) %>%
      bind_rows(simpleHL) %>% 
      filter(complete.cases(.)) %>%
      dplyr::rename(Datum = classifiedTides) %>%
      mutate(Datum = recode(Datum, "LH"="MLHW", "HH"="MHHW", "HHS"="MHHWS", "HL"="MHLW", "LL"="MLLW", "LLS"="MLLWS"))
    
    # wlTableSubset <- allWaterLevels %>% 
    #   filter((dateTime >= ymd_hm(startDate, tz="UTC")) & (dateTime <= ymd_hm(endDate, tz="UTC")))
    # 
    otherDatums1 <- as_tibble(data.frame(Datum = c("MTL"),
                                         wl=c( (simpleHL$wl[simpleHL$classifiedTides == "MHW"] - simpleHL$wl[simpleHL$classifiedTides == "MLW"])/2 + simpleHL$wl[simpleHL$classifiedTides == "MLW"]),
                                         n=c(nrow(classifiedTides)),
                                         stringsAsFactors = F))
    
    datum_summaries <- datum_summaries %>% 
      bind_rows(otherDatums1) %>% 
      arrange(-wl)
    
    otherDatums2 <- as_tibble(data.frame(Datum = c("HOT", "LOT"), 
                                         wl=c(max(hlTable$waterLevel, na.rm=T), min(hlTable$waterLevel, na.rm=T)),
                                         n=c(1, 1),
                                         stringsAsFactors = F))
    
    datum_summaries <- datum_summaries %>% 
      bind_rows(otherDatums2) %>% 
      arrange(-wl)
    
    if (graph == T) {
      print("  ... graphing")
      # Coded datums, observed and predicted WL by fractional month
      tidesPlot <- classifiedTides %>%
        mutate(day_in_month = days_in_month(dateTime),
               f_day_of_month = ((day + (hour(dateTime)/24 + (minute(dateTime)/(60*24))))) /
                 days_in_month(month)
        ) %>% 
        dplyr::rename(wl = waterLevel)
      
      # Add predicted to wl table
      # Change waterLevel name to measured
      # Gather excluding datetime
      
      # wlPlotting <- allWaterLevels %>% 
      #   filter((dateTime >= ymd_hm(startDate, tz="UTC")) & (dateTime <= ymd_hm(endDate, tz="UTC"))) %>%
      #   # left_join(wlTable, by = "dateTime") %>%
      #   # mutate(anomaly = predicted - waterLevel) %>% 
      #   #select(-waterLevel) %>% 
      #   tidyr::gather(key = measuredOrModeled, value = wl, -dateTime) %>%
      #   arrange(dateTime, measuredOrModeled) %>%
      #   mutate(year = year(dateTime), 
      #          month = month(dateTime), 
      #          day_in_month = days_in_month(dateTime),
      #          f_day_of_month = ((day(dateTime) + (hour(dateTime)/24 + (minute(dateTime)/(60*24))))) /
      #            days_in_month(month)
      #   )
      
      wl_plot <- ggplot(data = tidesPlot, aes(x=f_day_of_month, y=wl)) +
        facet_wrap(year~month) +
        geom_line(color = "lightblue") +
        geom_point(data = tidesPlot, aes(shape = classifiedTides)) +
        ylab("Water Level (m)") +
        xlab("Month (fraction)") +
        theme_minimal() +
        scale_x_continuous(labels=scaleFUN) +
        theme(legend.title = element_blank())
      out_fig <- paste(gauge_data, "/", out_fig_name, ".pdf", sep = "")
      ggsave(out_fig, wl_plot, height = 8.5, width = 11)
    }
    
    return(list(classifiedTides, datum_summaries))
    
  }
  
  
  getDatumsForGaugeHl <- function(station_id=9410660, 
                                  startDate = '2015-01-01 00:00',
                                  endDate = '2015-12-31 23:59', 
                                  bufferStart = 15,
                                  bufferEnd = 15,
                                  datum="NAVD",
                                  graph=F,
                                  gauge_data = "biomass/agb_biomass_inundation/NOAA_water_levels") {
    
    print(paste("Analysing gauge ", as.character(station_id), " from ", 
                startDate, " to ", endDate, ".", sep=""))
    
    # We need a buffer start and end because, 1. We need to make sure we have the edges of the analysis period covered for
    # getting MHHWS and stuff.
    # Sometimes the longest period harmonic, the solar annual one, throws an error if you analyse over to small a period.
    bufferStartDate <- lubridate::ymd_hm(startDate) - lubridate::days(bufferStart)
    bufferStartDate <- toString(format(bufferStartDate, "%Y-%m-%d %H:%M"))
    
    bufferEndDate <- lubridate::ymd_hm(endDate) + lubridate::days(bufferEnd)
    bufferEndDate <- toString(format(bufferEndDate, "%Y-%m-%d %H:%M"))
    
    # First check to see if the data was already downloaded
    output_file <- paste(gauge_data, "/", as.character(station_id), ".csv", sep="")
    if (file.exists(output_file)) {
      
      # If it is load it
      station_file <- read_csv(output_file, col_types = "Tnc")
      
      # See if it has the right number of columns
      time_steps_needed <- seq(as.POSIXct(bufferStartDate, tz = "UTC"),
                               as.POSIXct(bufferEndDate, tz = "UTC"), by="6 min")
      
      if (all(time_steps_needed %in% station_file$dateTime)) {
        print("  already downloaded.") 
        wlTable <- station_file %>% 
          filter(dateTime >= ymd_hm(bufferStartDate, tz = "UTC") & dateTime <= ymd_hm(bufferEndDate, tz="UTC"))
      } else {
        # If not download, bind rows, arrange by date time and write to file
        # Run the download function
        print("  downloading and adding to file...")
        wlTable <- downloadHLData(station_id=station_id, 
                                  startDate = bufferStartDate,
                                  endDate = bufferEndDate, 
                                  datum=datum)
        
        # Filter out any obs that are in the file already
        station_file <- station_file %>%
          filter(!(dateTime %in% time_steps_needed))
        
        station_file <- station_file %>% 
          bind_rows(wlTable) %>% 
          arrange(dateTime) %>% 
          distinct()
        
        write_csv(station_file, output_file)
      }
    } else {
      # If not download, bind rows, arrange by date time and write to file
      # Run the download function
      print("  downloading and adding to file...")
      wlTable <- downloadHLData(station_id=station_id, 
                                startDate = bufferStartDate,
                                endDate = bufferEndDate, 
                                datum=datum)
      
      # Filter out any obs that are in the file already
      # station_file <- station_file %>%
      #   filter(!(dateTime %in% time_steps_needed))
      
      station_file <- wlTable %>% 
        arrange(dateTime) %>% 
        distinct()
      
      write_csv(station_file, output_file)
    }
    
    print("  ... fitting custom tidal datum.")
    
    out_fig_name2 <- paste(as.character(station_id), " ", 
                           as.character(startDate),
                           " to ",
                           as.character(endDate),
                           # ".pdf", # DK commented this out on 2020-05-17 
                           # because the file extension is added 
                           # later on in the script
                           sep="")
    
    tidalDatums <- parseHLdata(hlTable = wlTable, startDate = startDate,
                               endDate = endDate,
                               bufferStart = bufferStart,
                               bufferEnd = bufferEnd,
                               graph = graph,
                               out_fig_name = out_fig_name2,
                               gauge_data = gauge_data)
    
    tidalDatums[[2]] <- tidalDatums[[2]] %>% 
      mutate(station_id = station_id,
             startDate = startDate,
             endDate = endDate) %>%
      dplyr::rename(meters = wl) %>% 
      dplyr::select(station_id, startDate, endDate, Datum, meters, n)
    
    # Add to the master file
    # Rewrite master file
    
    return_datums <- tidalDatums[[2]]
    
    wlTable2 <- wlTable %>% 
      filter(dateTime >= ymd_hm(startDate) &
               dateTime <= ymd_hm(endDate))
    
    inundationProfile <- getPctInundationProfile(wlTable2)
    
    return(list(return_datums, inundationProfile))
    
  }
  
  getPctInundationProfile <- function(wlTable) {
    targetElevations <- seq(round(min(wlTable$waterLevel, na.rm=T),2),
                            round(max(wlTable$waterLevel, na.rm=T),2),
                            0.01)
    
    fInundation <- c()
    for (i in 1:length(targetElevations)) {
      
      inundationEvents <- filter(wlTable, waterLevel >= targetElevations[i])
      fInundation <- c(fInundation, 
                       nrow(inundationEvents)/nrow(wlTable))
    }
    
    return(data.frame(elevation = targetElevations,
                      fractionInundation = fInundation))
    
  }
  
  batchCalculateDatumsHl <- function(startYear = 1983,
                                  endYear = 2020,
                                  station_id = 9410660,
                                  temp_file = "temp/gauge_data") {
    j = 1
    all_datums <- data.frame()
    # Lunar Nodal Cycle
    for (i in startYear:endYear) {
      try({
          
        datums <- getDatumsForGaugeHl(station_id = station_id,
                                      startDate = paste(i, "-01-01 00:00", sep=""),
                                      endDate = paste(i, "-12-31 23:59", sep=""),
                                      graph = T,
                                      gauge_data = temp_file)
        
        if (j == 1) {
          all_datums <- datums[[1]]
        } else {
          all_datums <- bind_rows(all_datums,
                                  datums[[1]])
        }
        j = j+1
        
      
      })
    }
    
    if (nrow(all_datums)>=1) {
      return(all_datums)
    } else {
      print("No data found.")
    }
  }
  
  calculateLunarNodalCycle <- function(all_datums, temp_file = "temp/lunar_nodal") {
    if(! all(c("startDate", "Datum", "meters", "station_id") %in% 
             names(all_datums))) {
      break("Input datum table formatted in correctly")
    }
    
    # datumsPlot <- all_datums %>% 
    #   dplyr::select("startDate", "Datum", "meters") %>% 
    #   mutate(startDate = ymd_hm(startDate)) %>% 
    #   filter(Datum %in% c("MTL", "MLHW", "MHHW", "MHHWS"))
    
    # ggplot(data = datumsPlot, aes(x = startDate, y=meters)) +
    #   geom_point(aes(shape = Datum, color = Datum)) +
    #   geom_line(aes(color = Datum))
    
    datumsPlot <- all_datums %>% 
      dplyr::select("startDate", "Datum", "meters") %>% 
      mutate(startDate = ymd_hm(startDate),
             year = year(startDate)+0.5) %>% 
      filter(Datum %in% c("MTL", "MHW")) %>% 
      select(-startDate) %>% 
      spread(key = Datum, value = meters) %>% 
      mutate(amplitude = MHW - MTL)
    
    meanAmp <- mean(datumsPlot$amplitude, na.rm = T)
    
    fit1 <-  nls(amplitude ~ meanAmp + lunarNodalAmp * sin(2 * pi * (year - lunarNodalPhase)/18.61), 
                 data = datumsPlot,
                 start = list(lunarNodalAmp = 1, lunarNodalPhase = 1))
    
    lunarNodalAmp <- coef(fit1)[1]
    lunarNodalPhase <- coef(fit1)[2]
    lonarNodalP <-  summary(fit1)$coefficients[1:2,4]
    
    if(all(lonarNodalP <= 0.05)) {
      print("Significant lunar nodal cycle fit.")
      
      startYear <- 1983
      endYear <- 2020
      
      predict_line <- data.frame(year = (startYear:endYear)+0.5)
      predict_line$amplitude <- meanAmp + lunarNodalAmp * sin(2 * pi * (predict_line$year - lunarNodalPhase)/18.61)
                 
      ggplot() +
        geom_line(data = predict_line, aes(x = year, y=amplitude)) +
        geom_point(data = datumsPlot, color = "lightgrey", aes(x = year, y = amplitude)) +
        ylab("Tidal Amplitude (m)") +
        xlab(NULL) +
        theme_minimal() +
        ggtitle(paste("Station:", station_id, "- 18.6 year cycle", sep = " ")) +
        theme(axis.text.x = element_text(angle = 45, hjust=1))
      
      ggsave(paste0(temp_file, "/", unique(all_datums$station_id), "_lunarNodalCycle.pdf"),
             height = 3.5, width=3.5)
      
        } else {
      print("No lunar nodal activity detected.")
          lunarNodalAmp <- 0
    }
    
    output_file <- data.frame(station_id = rep(first(unique(all_datums$station_id)), 2),
                              variable = c("lunarNodalAmp",
                                            "lunarNodalPhase"),
                              value = c(lunarNodalAmp, 
                                        lunarNodalPhase),
                              p = lonarNodalP,
                              stringsAsFactors = F)
    return(output_file)
    
  }

  
}

gauge_df <- data.frame(Name = c("Seattle, WA", "San Francisco, CA", "Port Isabel, TX",
                                "Pensacola, FL", "Charleston, SC", "Annapolis, MD",
                                "Portland, ME"),
                       station_id = c(9447130, 9414290, 8779770,
                                      8729840, 8665530, 8575512,
                                      8418150),
                       Latitude = c(47+36.1/60, 37+48.4/60, 26+3.7/60,
                                    30+24.3/60, 32+46.8/60, 38+59.0/60,
                                    43+39.5/60),
                       Longitude = c(-(122+2.04/60), -(122+28.0/60), -(97+12.9/60),
                                     -(87+12.7/60), -(79+55.4/60), -(76+28.9/60),
                                     -(70+14.7/60)))

# Notes 
{
  # Seattle, WA - Station ID: 9447130
  # Latitude	47° 36.1 N
  # Longitude	122° 2.04 W
  # MHHW	2.748	Mean Higher-High Water
  # MHW	2.484	Mean High Water
  # MSL	1.309	Mean Sea Level
  # NAVD88	0
  # 2.07 mm/yr
  
  # San Francisco, CA - Station ID: 9414290
  # Latitude	37° 48.4 N
  # Longitude	122° 28.0 W
  # MHHW# 	1.798
  # MHW	1.612
  # MSL	0.969
  # NAVD88	0.000
  # 1.97 mm/yr
  
  # Port Isabel, TX - Station ID: 8779770
  # Latitude	26° 3.7 N
  # Longitude	97° 12.9 W
  # MHHW	0.159	Mean Higher-High Water
  # MHW	0.141	Mean High Water
  # MSL	-0.012	Mean Sea Level
  # NAVD88	0.000	North American Vertical Datum of 1988
  # 4.25 mm/yr
  
  # Dauphin Island, AL - Station ID: 8735180
  # Latitude	30° 15.0 N
  # Longitude	88° 4.5 W
  # MHHW	0.213	Mean Higher-High Water
  # MHW	0.207	Mean High Water
  # MTL	0.026	Mean Tide Level
  # MSL	0.016	Mean Sea Level
  # NAVD88	0.000	North American Vertical Datum of 1988
  # 4.25 mm/yr
  
  # Charleston, Cooper River Entrance, SC - Station ID: 8665530
  # Latitude	32° 46.8 N
  # Longitude	79° 55.4 W
  # MHHW	0.800	Mean Higher-High Water
  # MHW	0.691	Mean High Water
  # MSL	-0.067	Mean Sea Level
  # 3.39 mm/yr
  
  # Annapolis, MD - Station ID: 8575512
  # Latitude	38° 59.0 N
  # Longitude	76° 28.9 W
  # MHHW	0.203	Mean Higher-High Water
  # MHW	0.129	Mean High Water
  # MSL	-0.016	Mean Sea Level
  # 3.73 mm/yr
  
  # Portland, ME - Station ID: 8418150
  # Latitude	43° 39.5 N
  # Longitude	70° 14.7 W
  # MHHW	1.418	Mean Higher-High Water
  # MHW	1.285	Mean High Water
  # MSL	-0.095	Mean Sea Level
  # 1.90 mm/yr
  
}

for (i in 1:nrow(gauge_df)) {
  
  station_id <- gauge_df$station_id[i]
  
  batchDatums <- batchCalculateDatumsHl(startYear = 1983, endYear = 2020,
                                        station_id = station_id)
  
  batchDatumsMTL <- dplyr::filter(batchDatums, Datum == "MTL") %>% 
    dplyr::mutate(startDate = ymd_hm(startDate),
           endDate = ymd_hm(endDate))
  
  batchDatumsMTL_1983to2001 <- batchDatums %>% 
    dplyr::mutate(startDate = ymd_hm(startDate),
           endDate = ymd_hm(endDate)) %>% 
    dplyr::filter(startDate >= ymd("1983-01-01"),
           endDate <= ymd("2001-12-31")) %>% 
    dplyr::mutate(d_year = year(startDate) + 0.5)
  
  batchDatumsMTL_1983to2001_summarized <- batchDatumsMTL_1983to2001 %>% 
    dplyr::group_by(station_id, Datum) %>% 
    dplyr::summarise(meters = mean(meters, na.rm = T)) %>% 
    dplyr::arrange(-meters) %>% 
    dplyr::ungroup() %>% 
    dplyr::rename(variable = Datum,
           value = meters) %>% 
    dplyr::mutate(value = value * 100) # convert from meters to cm
  
  # Filter MTL to 1983 to 2001
  batchDatumsMTL_1983to2001_slr <- batchDatumsMTL_1983to2001 %>% filter(Datum == "MTL")
  
  ggplot(data = batchDatumsMTL_1983to2001_slr, aes(x = startDate, y = meters)) +
    geom_point()
  
  # Calculate sea level rise as a linear model 
  slr <- lm(meters~d_year, data = batchDatumsMTL_1983to2001_slr)
  
  slr_df <- data.frame(station_id = gauge_df$station_id[i],
                       variable = "relSeaLevelRiseInit",
                       value = coefficients(slr)[2]*100,
                       p = summary(slr)$coefficients[2,4])
  
  lunarNodal <- calculateLunarNodalCycle(batchDatums) %>% 
    dplyr::mutate(value = ifelse(variable == "lunarNodalAmp", value * 100, value))
  
  batchDatumsMTL_2000 <- batchDatumsMTL_1983to2001 %>% 
    dplyr::filter(Datum == "MTL" & year(startDate) == 2000) %>% 
    dplyr::mutate(value = meters*100,
           variable = "meanSeaLevel") %>% 
    dplyr::select(station_id, value, variable)
  
  
  summary_df <- dplyr::bind_rows(batchDatumsMTL_1983to2001_summarized,
                          lunarNodal,
                          slr_df,
                          batchDatumsMTL_2000
  )
  
  rownames(summary_df) <- NULL
  
  if (i == 1) {
    all_sites_and_variables <- summary_df
  } else {
    all_sites_and_variables <- dplyr::bind_rows(all_sites_and_variables, summary_df)
  }
  
}

write_csv(all_sites_and_variables, "data/parameterSet/All_sites_and_variables.csv")

all_sites_and_variables <- read_csv("data/parameterSet/All_sites_and_variables.csv")

lunarNodalAmpPlot <- all_sites_and_variables %>% 
  left_join(gauge_df) %>% 
  filter(variable == "lunarNodalAmp")

ggplot(lunarNodalAmpPlot, aes(x = Longitude, y = Latitude)) +
  geom_point(aes(color = value))

lunarNodalAPhasePlot <- all_sites_and_variables %>% 
  left_join(gauge_df) %>% 
  filter(variable == "lunarNodalPhase")

ggplot(lunarNodalAPhasePlot, aes(x = Longitude, y = Latitude)) +
  geom_point(aes(color = value))

kopp_2014 <- read_csv("data/Kopp_et_al_2014/Kopp_AllCompleteNorthAmericanSLR_inMeters.csv") %>% 
  select(gauge_id, rcp45_50_2100) %>% 
  dplyr::rename(station_id = gauge_id, relSeaLevelRiseTotal = rcp45_50_2100) %>% 
  mutate(relSeaLevelRiseTotal = relSeaLevelRiseTotal * 100)

library(rCMEM)

all_sites_and_variables2_wide <- all_sites_and_variables %>% 
  select(-p) %>% 
  spread(key = variable, value = value) %>% # relSeaLevelRiseInit = NA
  left_join(gauge_df) %>% 
  left_join(kopp_2014) %>%   # relSeaLevelRiseTotal = NA,
  dplyr::rename(meanSeaLevelDatum = MTL) %>% 
  dplyr::mutate(meanHighWaterDatum = ifelse(station_id %in% c(8729840, 8779770),
                                     MHW, MLHW),
         meanHighHighWaterDatum = ifelse(station_id %in% c(8729840, 8779770),
                                         NA, MHHW),
         meanHighHighWaterSpringDatum = ifelse(station_id %in% c(8729840, 8779770),
                                               NA, MHHWS)) %>% 
  dplyr::mutate(nFloods = ifelse(station_id %in% c(8729840, 8779770),
                          705.98/2, 705.98),
         captureRate=ifelse(station_id %in% c(8729840, 8779770), 
                            2.8*2, 2.8)) %>% 
  # Let's make initElv == peak biomass for Spartina according to Morris's North Inlet Marsh Organs
  dplyr::mutate(initElv = convertZStarToZ(0.831, meanSeaLevel+(meanHighWaterDatum-meanSeaLevelDatum), meanSeaLevel),
         )  %>%
  select(station_id, Name, relSeaLevelRiseInit, relSeaLevelRiseTotal, initElv,
         meanSeaLevel, meanSeaLevelDatum, meanHighWaterDatum, meanHighHighWaterDatum, meanHighHighWaterSpringDatum,
         lunarNodalAmp, lunarNodalPhase, nFloods, captureRate)
